// Updated version in DSL2

// this consumes fastq files, adapter trims, then aligns reads to the specified reference using bwa-meth
// mode can be set to tile_fastqs or run_fastqs depending on whether the system should map each tile's reads in distinct jobs then combine (tile_fastqs)
// or all reads for a library in a single job (run_fastqs)

nextflow.enable.dsl=2

genome = params.genome
threads = params.threads ?:4

tmp_dir = params.tmp_dir
fastq_glob = params.fastq_glob
flowcell = params.flowcell ?: 'unknown'

outputPath = params.outdir ?:'output'
species = params.species ?:'human'

qc_report = params.qc ?:true

workflow {
    fq_set_channel = Channel.fromFilePairs(fastq_glob, checkIfExists:true).view()
    

    mapping(fq_set_channel)
    sampleID = mapping.out.sampleID
    mergeAndMarkDuplicates(mapping.out.aligned_files)
    
    md_bams = mergeAndMarkDuplicates.out.md_bams
    
    methylDackel_mbias(md_bams)
    methylDackel_extract(md_bams)
    
    combine_mbias_tsv(methylDackel_mbias.out.mbias_output_tsv)
    combine_mbias_svg(methylDackel_mbias.out.mbias_output_svg)
                      
    if(qc_report){
        runFastQC(md_bams)
    
        sum_nonconverted_reads(mapping.out.nonconverted_counts)
        combine_nonconversion(sum_nonconverted_reads.out)

        samtools_flagstats(md_bams)
        samtools_stats(md_bams)

        picard_gc_bias(md_bams)
        picard_stats(md_bams)
        goleft(md_bams)

        if(species == 'human'){
            select_human_reads(md_bams)
            human_gc_bias(select_human_reads.out)
            human_insert_size(select_human_reads.out)
        }else {
            select_mouse_reads(md_bams)
            mouse_insert_size(select_mouse_reads.out)
        }
        
        qc_files = sampleID
                    .join(runFastQC.out.fastqc_results)
                    .join(samtools_flagstats.out.flagstats)
                    .join(samtools_flagstats.out.idxstats)
                    .join(samtools_stats.out.samstats)
                    .join(picard_stats.out.picard_stats)
                    .join(picard_gc_bias.out.picard_gc_stats)
                    .join(goleft.out.goleft_ped)
                    .join(goleft.out.goleft_roc)
                    .join(mergeAndMarkDuplicates.out.samblaster_logs)
                    .join(mapping.out.fastp_log_files)
        multiqc(qc_files)                         
    }     
}





process mapping {
    cpus {threads}
    errorStrategy 'retry'
    tag { [flowcell, library] }
    conda "bwameth=0.2.2 seqtk=1.3 sambamba=0.7.0 fastp=0.20.1 mark-nonconverted-reads=1.1"

    input:
        
        tuple val(library), path(reads)
    
    output:
        val(library), emit: sampleID
        tuple val(library), path("*.aln.bam"), emit: aligned_files
        tuple val(library), path("*.nonconverted.tsv"), emit: nonconverted_counts
        tuple val(library), path("*_fastp.json"), emit: fastp_log_files

    
    shell:
    '''
    inst_name=$(zcat -f '!{reads[0]}' | head -n 1 | cut -f 1 -d ':' | sed 's/^@//')
    fastq_barcode=$(zcat -f '!{reads[0]}' | head -n 1 | sed -r 's/.*://')

    if [[ "${inst_name:0:2}" == 'A0' ]] || [[ "${inst_name:0:2}" == 'NS' ]] || \
       [[ "${inst_name:0:2}" == 'NB' ]] || [[ "${inst_name:0:2}" == 'VH' ]] ; then
       trim_polyg='--trim_poly_g'
       echo '2-color instrument: poly-g trim mode on'
    else
       trim_polyg=''
    fi
    seqtk mergepe !{reads} \
    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis \
            -j "!{library}_fastp.json" 2> fastp.stderr \
    | bwameth.py -p -t !{threads} --read-group "@RG\\tID:${fastq_barcode}\\tSM:!{library}" \
            --reference !{genome} /dev/stdin 2>  "!{library}.log.bwamem" \
    | mark-nonconverted-reads.py 2> "!{library}.nonconverted.tsv" \
    | sambamba view -t !{threads} -S -f bam -o "!{library}.aln.bam" /dev/stdin 2> sambamba.stderr;
    '''

}

process mergeAndMarkDuplicates {
    cpus {threads}
    errorStrategy 'retry'
    tag { library }
    publishDir "${outputPath}", mode: 'copy', pattern: '*.{md.bam}*'
    conda "samtools=1.9 samblaster=0.1.24 sambamba=0.7.0"

    input:
        tuple val(library), path(libraryBam)

    output:
        tuple val(library), path('*.md.bam'), path('*.md.bam.bai'), emit: md_bams
        tuple val(library), path('*.samblaster'), emit: samblaster_logs

    shell:
    '''
    samtools cat  -b <( find . -name '*.aln.bam' ) \
    | samtools view -h /dev/stdin \
    | samblaster 2> !{library}.log.samblaster \
    | sambamba view -t !{threads} -l 0 -S -f bam /dev/stdin \
    | sambamba sort --tmpdir=!{tmp_dir} -t !{threads} -m 20GB -o !{library}.md.bam /dev/stdin

    '''
}


    process methylDackel_mbias {
        cpus {threads}
        errorStrategy 'retry'
        tag {library}
        conda "methyldackel=0.4.0 samtools=1.9"

        input:
            tuple val(library), path(md_file), path(md_bai)

        output:
            tuple val(library), path('*.svg'), emit: mbias_output_svg
            tuple val(library), path('*.tsv'), emit: mbias_output_tsv

        shell:
        '''
        echo -e "chr\tcontext\tstrand\tRead\tPosition\tnMethylated\tnUnmethylated\tnMethylated(+dups)\tnUnmethylated(+dups)" > !{library}_combined_mbias.tsv
        chrs=(`samtools view -H !{md_file} | grep @SQ | cut -f 2 | sed 's/SN://'| grep -v _random | grep -v chrUn | sed 's/|/\\|/'`)

        for chr in ${chrs[*]}; do
            for context in CHH CHG CpG; do
                arg=''
                if [ $context = 'CHH' ]; then
                arg='--CHH --noCpG'
                elif [ $context = 'CHG' ]; then
                arg='--CHG --noCpG'
                fi
                # need two calls to add columns containing the counts without filtering duplicate reads (for rrEM-seq where start/end is constrained)
                # not sure why we need both --keepDupes and -F, probably a bug in mbias
                join -t $'\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \
                <( \
                    MethylDackel mbias --noSVG $arg -@ !{threads} -r $chr !{genome} !{md_file} | \
                    tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
                ) \
                <( \
                    MethylDackel mbias --noSVG --keepDupes -F 2816 $arg -@ !{threads} -r $chr !{genome} !{md_file} | \
                    tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
                ) \
                | sed "s/^/${chr}\t${context}\t/" \
                >> !{library}_combined_mbias.tsv
            done
        done
        # makes the svg files for trimming checks
        MethylDackel mbias -@ !{threads} --noCpG --CHH --CHG -r ${chrs[0]} !{genome} !{md_file} !{library}_chn
        for f in *chn*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CHN <\\/text>/" $f; done;

        MethylDackel mbias -@ !{threads} -r ${chrs[0]} !{genome} !{md_file} !{library}_cpg
        for f in *cpg*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CpG<\\/text>/" $f; done;

        '''

    }

    process methylDackel_extract {
        cpus {threads}
        tag {library}
        publishDir "${outputPath}", mode: 'copy'
        conda "methyldackel=0.4.0 pigz=2.4"

        input:
            tuple val(library), path(md_file), path(md_bai)

        output:
            tuple val(library), path('*.methylKit.gz'), emit: extract_output

        shell:
        '''
        MethylDackel extract --methylKit --nOT 0,0,0,5 --nOB 0,0,5,0 -@ !{threads} --CHH --CHG -o !{library} !{genome} !{md_file}
        pigz -p !{threads} *.methylKit
        '''

    }
    
    
    process runFastQC {
        cpus {threads}
        errorStrategy 'retry'
        tag { library }
        conda "fastqc=0.11.8"

        input:
            tuple val(library), file(md_file), file(md_bai)

        output:
            tuple val(library), path('*_fastqc.zip'), emit: fastqc_results

        shell:
        '''
        fastqc -f bam !{md_file}
        '''

    }

    process sum_nonconverted_reads {	

        input:	
            tuple val(library), path(count_files)

        output:	
            tuple val(library), path('*-nonconverted-counts.tsv')

        shell:	
        '''	
        files=(*.tsv)	
        paste *.tsv | awk -v numFiles=${#files[@]} -v OFS='\t' '	
        {	
        row = sep = ""	
        for(i=1; i < NF/numFiles; ++i) { row = row sep $i; sep = OFS }	
        sum = $(NF/numFiles) # last header col. / (1st) data col. to sum	
        for(i=2; i<=numFiles; ++i) sum += $(NF/numFiles * i) # add other cols.	
        printf "%s%s%s\\n", row, OFS, sum	
        }' > tmp-counts.tsv	
        awk '{print "!{library}\t" $0}' tmp-counts.tsv > !{library}-nonconverted-counts.tsv	
        '''    	
    }	

    process combine_nonconversion {	
        publishDir "${outputPath}", mode: 'copy'	

        input:	
            tuple val(library), path(nonconverted_files)

        output:	
            path("${library}_combined-nonconverted.tsv")	

        shell:	
        '''	
        cat !{nonconverted_files} > !{library}_combined-nonconverted.tsv	
        '''	

    }

    process samtools_flagstats {
        cpus {threads}
        errorStrategy 'retry'
        tag { library }
        conda "samtools=1.9"

        input:
            tuple val(library), path(md_file), path(md_bai)

        output:
            tuple val(library), path('*.flagstat'), emit: flagstats
            tuple val(library), path('*.idxstat'), emit: idxstats
            

        shell:
        '''
        samtools flagstat -@ !{threads} !{md_file} > !{md_file}.flagstat
        samtools idxstats !{md_file} > !{md_file}.idxstat
        '''
    }

    process samtools_stats {
        cpus {threads}
        errorStrategy 'retry'
        tag { library }
        conda "samtools=1.9"

        input:
            tuple val(library), path(md_file),path(md_bai)

        output:

            tuple val(library), path('*.samstat'), emit: samstats

        shell:
        '''
        samtools stats -@ !{threads} !{md_file} > !{md_file}.samstat
        '''
    }

    process picard_gc_bias {
        cpus {threads}
        errorStrategy 'retry'
        tag { library }
        conda "picard=2.20.7"

        input:
            tuple val(library), path(md_file), path(md_bai)

        output:
            tuple val(library), path('*gc_metrics'), emit: picard_gc_stats

        shell:
        '''
        picard -Xmx4g CollectGcBiasMetrics IS_BISULFITE_SEQUENCED=true VALIDATION_STRINGENCY=LENIENT I=!{md_file} O=!{md_file}.gc_metrics S=!{md_file}.gc_summary_metrics CHART=!{md_file}.gc.pdf R=!{genome}
        '''
    }

    process picard_stats {

        cpus {threads}
        errorStrategy 'retry'
        tag { library }
        conda "picard=2.20.7"

        input:
            tuple val(library), path(md_file), path(md_bai)

        output:
            tuple val(library), path('*_metrics'), emit: picard_stats

        shell:
        '''
        picard -Xmx16g CollectInsertSizeMetrics VALIDATION_STRINGENCY=LENIENT  I=!{md_file} O=!{md_file}.insertsize_metrics MINIMUM_PCT=0 HISTOGRAM_FILE=/dev/null
        '''
    }
    
    process goleft {
        cpus {threads}
        conda 'goleft=0.2.0'

        input:
            tuple val(library), path(md_file), path(md_bai)

        output:
            tuple val(library), path("${library}/*-indexcov.ped"), emit: goleft_ped
            tuple val(library), path("${library}/*-indexcov.roc"), emit: goleft_roc

        shell:
        '''
            goleft indexcov --directory !{library} *.bam
        '''
    }
    
    process select_human_reads {
        cpus {threads}
        tag {library}
        conda "sambamba=0.7.0 bedtools=2.29.2"

        input:
            tuple val(library), path(md_file), path(md_bai)

        output:
            tuple val(library), path('*.human.bam')

        shell:
        '''
        sambamba view -t !{threads} -l 0 -f bam !{md_file} chr1 chr2 chr3 chr4 chr5 chr6 \
                                                  chr7 chr8 chr9 chr10 chr11 chr12 \
                                                  chr13 chr14 chr15 chr16 chr17 chr18 \
                                                  chr19 chr20 chr21 chr22 chrX chrY \
        > !{md_file}.human.bam
        '''

    }
    
    process human_gc_bias {
        cpus {threads}
        errorStrategy 'retry'
        tag { library }
        conda "picard=2.20.7"

        input:
            tuple val(library), path(md_file)

        output:
            tuple library, path('*gc_metrics')

        shell:
        '''
        picard -Xmx4g CollectGcBiasMetrics IS_BISULFITE_SEQUENCED=true VALIDATION_STRINGENCY=LENIENT I=!{md_file} O=!{md_file}.gc_metrics S=!{md_file}.gc_summary_metrics CHART=!{md_file}.gc.pdf R=!{genome}
        '''
    }

    process human_insert_size {

        cpus {threads}
        errorStrategy 'retry'
        tag { library }
        conda "picard=2.20.7"

        input:
            tuple val(library), path(md_file)

        output:
            tuple val(library), path('*insertsize_metrics')

        shell:
        '''
        picard -Xmx16g CollectInsertSizeMetrics VALIDATION_STRINGENCY=LENIENT  I=!{md_file} O=!{md_file}.insertsize_metrics MINIMUM_PCT=0.0001 HISTOGRAM_FILE=/dev/null
        '''
    }
    
    process select_mouse_reads {
        cpus {threads}
        tag {library}
        conda "sambamba=0.7.0 bedtools=2.29.2"

        input:
            tuple val(library), path(md_file), path(md_bai)

        output:
            tuple val(library), path('*.mouse.bam')

        shell:
        '''
        sambamba view -t !{threads} -l 0 -f bam !{md_file} chr1 chr2 chr3 chr4 chr5 chr6 \
                                                  chr7 chr8 chr9 chr10 chr11 chr12 \
                                                  chr13 chr14 chr15 chr16 chr17 chr18 \
                                                  chr19 chrX chrY \
        > !{md_file}.mouse.bam
        '''

    }
    
    process mouse_insert_size {

        cpus {threads}
        errorStrategy 'retry'
        tag { library }
        conda "picard=2.20.7"

        input:
            tuple val(library), path(md_file)

        output:
            tuple val(library), path('*insertsize_metrics')

        shell:
        '''
        picard -Xmx16g CollectInsertSizeMetrics VALIDATION_STRINGENCY=LENIENT  I=!{md_file} O=!{md_file}.insertsize_metrics MINIMUM_PCT=0.0001 HISTOGRAM_FILE=/dev/null
        '''
    }
    
    
    process multiqc {
        cpus {threads}
        publishDir "${outputPath}", mode: 'copy'
        tag { library }
        conda "multiqc=1.17"
        
        input:
            tuple val(library), path("*")
        output:
            path("${library}_multiqc_report.html")

        shell:
        '''
        
        cat <<CONFIG > multiqc_config.yaml 
    title: Bwameth Alignment Summary - !{flowcell}
    extra_fn_clean_exts:
        - '.md'
        - '_fastp'
    custom_plot_config:
        picard_insert_size:
            xmax: 1000
    table_columns_placement:
        Samtools Stats:
            raw_total_sequences: 10
            reads_mapped_percent: 20
            reads_properly_paired_percent: 30
            reads_MQ0_percent: 35
        Samblaster:
            pct_dups: 40
        Picard:
            summed_median: 50
    table_columns_visible:
        Picard:
            PCT_PF_READS_ALIGNED: False
            summed_mean: False
        Samtools Stats:
            reads_mapped: False
            mapped_passed: False
            non-primary_alignments: False
            reads_MQ0_percent: True
        Samtools Flagstat:
            mapped_passed: False
        samtools_idxstats_always:
            - plasmid_puc19c
            - phage_lambda
        FastQC:
            percent_duplicates: False
            total_sequences: False
            avg_sequence_length: False
            percent_fails: False
            total_sequences: False
    CONFIG
        multiqc -ip  . -n !{library}_multiqc_report.html
        '''
    }

    process combine_mbias_tsv {
        publishDir "${outputPath}", mode: 'copy', pattern: '*combined*'

        input:
            tuple val(library), path('*_mbias.tsv')

        output:
            path("*_combined-mbias.tsv")

        shell:
        '''
            echo -ne 'flowcell\tlibrary\t' > !{library}_combined-mbias.tsv
            ls *_mbias.tsv | head -n 1 | xargs head -n 1 >> !{library}_combined-mbias.tsv
            for f in *_mbias.tsv; do 
                filebase=`basename "${f}" !{library}_combined_mbias.tsv`
                paste <( yes "!{flowcell}	!{library}" | head -n `nl "$f" | tail -n 1 | cut -f 1` ) "$f" | tail -n +2  >> !{library}_combined-mbias.tsv
            done
        '''
    }

    process combine_mbias_svg {
        publishDir "${outputPath}", mode: 'copy', pattern: '*combined*'
        conda 'cairosvg=2.4.2 ghostscript=9.22'

        input:
            tuple val(library), path('*.svg')

        output:
            path("*_combined-mbias.pdf")

        shell:
        '''
        for f in *.svg; do
            cairosvg <(sed s/-nan/-1/ $f) -o $f.pdf
        done
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=!{library}_combined-mbias.pdf *.pdf
        '''
    }


