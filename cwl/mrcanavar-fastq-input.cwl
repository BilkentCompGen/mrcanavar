class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: mrcanavar_fastq_input
baseCommand:
  - mrcanavar-auto
inputs:
  - id: inputfastq
    type: 'File[]'
    inputBinding:
      position: 0
      prefix: '--input'
      itemSeparator: ','
    label: Input FASTQ file(s)
    doc: The FASTQ files for the sample to be analyzed.
    'sbg:fileTypes': 'FASTQ, FASTQ.GZ'
  - id: maskedref
    type: File
    inputBinding:
      position: 0
      prefix: '--ref'
    label: Repeat masked reference genome
    doc: >-
      This reference genome should be repeat masked, including RepeatMasker and
      TRF masking. The index file for mrsFAST for this reference must also be
      present in the same directory.
    'sbg:fileTypes': FASTA
    secondaryFiles:
      - .fai
      - .index
  - id: config
    type: File
    inputBinding:
      position: 0
      prefix: '--conf'
    label: CNVR config file
    doc: >-
      The CNVR configuration file (generated with mrcanavar --prep). This must
      match the repeat masked reference genome.
    'sbg:fileTypes': CNVR
  - id: genelist
    type: File
    inputBinding:
      position: 0
      prefix: '--gene'
    label: Gene annotations
    doc: >-
      Coordinates of gene models in the genome in BED format. Additional columns
      are permitted.
    'sbg:fileTypes': BED
  - id: mem
    type: int?
    inputBinding:
      position: 0
      prefix: '--mem'
    label: Memory limit for mrsFAST mapping.
    doc: >-
      This value is passed directly to mrsFAST to set highest amount of memory
      used in gigabytes.
  - id: threads
    type: int?
    inputBinding:
      position: 0
      prefix: '--threads'
    label: Number of threads
    doc: >-
      This value is passed to mrsFAST to set the number of threads for read
      mapping.
  - id: cropfactor
    type: int?
    inputBinding:
      position: 0
      prefix: '--kmer'
    label: Read crop length
    doc: >-
      Default is 36. Do not input less than 25. This is the crop length for
      mrsFAST read mapping.
  - id: nocompress
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--no-gz'
    label: Uncompressed SAM file
    doc: >-
      Do not compress mrsFAST output. This option will generate larger files,
      but it will save some run time.
  - id: nosam
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--no-sam'
    label: Do not generate SAM(.gz) files.
    doc: >-
      Pipe the mrsFAST output directly into mrCaNaVaR, substantially reducing
      storage requirements. This option is automatically disabled if there are
      more than one FASTQ(.gz) file.
  - id: dryrun
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--dry-run'
    label: Dry Run
    doc: Do not run. Only print out command lines.
outputs:
  - id: dups
    label: Segmental duplications
    type: File?
    outputBinding:
      glob: '*.dups.bed'
      outputEval: '$(inheritMetadata(self, inputs.inputfastq))'
    'sbg:fileTypes': BED
  - id: dels
    label: Large deletions
    type: File?
    outputBinding:
      glob: '*.dels.bed'
      outputEval: '$(inheritMetadata(self, inputs.inputfastq))'
    'sbg:fileTypes': BED
  - id: genecopy
    label: Gene copy numbers
    type: File?
    outputBinding:
      glob: '*.genes.bed'
      outputEval: '$(inheritMetadata(self, inputs.inputfastq))'
    'sbg:fileTypes': BED
  - id: copynumbers
    label: Genome wide copy numbers
    type: File?
    outputBinding:
      glob: '*.copynumber.bed'
      outputEval: '$(inheritMetadata(self, inputs.inputfastq))'
    'sbg:fileTypes': BED
  - id: normdepth
    label: Normalized read depth counts
    type: 'File[]?'
    outputBinding:
      glob: '*out-*_norm.bed'
      outputEval: '$(inheritMetadata(self, inputs.inputfastq))'
    'sbg:fileTypes': BED
  - id: rawdepthcnt
    label: Raw read depth counts
    type: 'File[]?'
    outputBinding:
      glob: '*.txt'
      outputEval: '$(inheritMetadata(self, inputs.inputfastq))'
    'sbg:fileTypes': TXT
  - id: depthfile
    doc: Generated by mrcanavar --read
    label: DEPTH file
    type: File?
    outputBinding:
      glob: '*.depth'
      outputEval: '$(inheritMetadata(self, inputs.inputfastq))'
    'sbg:fileTypes': DEPTH
  - id: mapfile
    doc: >-
      SAM(.gz) files for mrsFAST mapping. Will not be generated if --no-sam is
      selected.
    label: Read mapping output
    type: 'File[]?'
    outputBinding:
      glob: '*.sam*'
      outputEval: '$(inheritMetadata(self, inputs.inputfastq))'
    'sbg:fileTypes': 'SAM, SAM.GZ'
label: mrCaNaVaR FASTQ input
requirements:
  - class: ResourceRequirement
    ramMin: 31000
    coresMin: 0
  - class: DockerRequirement
    dockerPull: alkanlab/mrcanavar
  - class: InlineJavascriptRequirement
    expressionLib:
      - |-

        var setMetadata = function(file, metadata) {
            if (!('metadata' in file))
                file['metadata'] = metadata;
            else {
                for (var key in metadata) {
                    file['metadata'][key] = metadata[key];
                }
            }
            return file
        };

        var inheritMetadata = function(o1, o2) {
            var commonMetadata = {};
            if (!Array.isArray(o2)) {
                o2 = [o2]
            }
            for (var i = 0; i < o2.length; i++) {
                var example = o2[i]['metadata'];
                for (var key in example) {
                    if (i == 0)
                        commonMetadata[key] = example[key];
                    else {
                        if (!(commonMetadata[key] == example[key])) {
                            delete commonMetadata[key]
                        }
                    }
                }
            }
            if (!Array.isArray(o1)) {
                o1 = setMetadata(o1, commonMetadata)
            } else {
                for (var i = 0; i < o1.length; i++) {
                    o1[i] = setMetadata(o1[i], commonMetadata)
                }
            }
            return o1;
        };
'sbg:categories':
  - Copy Number Polymorphism
  - WGS
  - DNA
  - Variant Calling
'sbg:links':
  - id: 'https://github.com/BilkentCompGen/mrcanavar'
    label: GitHub Repository
  - id: 'https://www.nature.com/articles/ng.437'
    label: Paper
'sbg:toolkitVersion': mrCaNaVaR
