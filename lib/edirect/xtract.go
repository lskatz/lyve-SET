// ===========================================================================
//
//                            PUBLIC DOMAIN NOTICE
//            National Center for Biotechnology Information (NCBI)
//
//  This software/database is a "United States Government Work" under the
//  terms of the United States Copyright Act. It was written as part of
//  the author's official duties as a United States Government employee and
//  thus cannot be copyrighted. This software/database is freely available
//  to the public for use. The National Library of Medicine and the U.S.
//  Government do not place any restriction on its use or reproduction.
//  We would, however, appreciate having the NCBI and the author cited in
//  any work or product based on this material.
//
//  Although all reasonable efforts have been taken to ensure the accuracy
//  and reliability of the software and data, the NLM and the U.S.
//  Government do not and cannot warrant the performance or results that
//  may be obtained by using this software or data. The NLM and the U.S.
//  Government disclaim all warranties, express or implied, including
//  warranties of performance, merchantability or fitness for any particular
//  purpose.
//
// ===========================================================================
//
// File Name:  xtract.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

/*
  test for presence of go compiler, cross-compile xtract executables, and pack into archive, by running:

  if hash go 2>/dev/null
  then
    env GOOS=darwin GOARCH=amd64 go build -o xtract.Darwin -v xtract.go
    env GOOS=linux GOARCH=amd64 go build -o xtract.Linux -v xtract.go
    env GOOS=windows GOARCH=386 go build -o xtract.CYGWIN_NT -v xtract.go
    tar -czf archive.tar.gz xtract.[A-Z]*
    rm xtract.[A-Z]*
  fi
*/

package main

import (
	"bytes"
	"container/heap"
	"fmt"
	"html"
	"io"
	"math"
	"os"
	"runtime"
	"runtime/debug"
	"runtime/pprof"
	"strconv"
	"strings"
	"sync"
	"time"
	"unicode"
)

// VERSION AND HELP MESSAGE TEXT

const xtractVersion = "5.90"

const xtractHelp = `
Overview

  Xtract uses command-line arguments to convert XML data into a tab-delimited table.

  -pattern places the data from individual records into separate rows.

  -element extracts values from specified fields into separate columns.

  -group, -block, and -subset limit element exploration to selected XML subregions.

Processing

  -cleanup         Fix non-ASCII spaces
  -compress        Compress runs of spaces
  -input           Read from file instead of stdin

Exploration Argument Hierarchy

  -pattern         Name of record within set
  -group             Use of different argument
  -block               names allows command-line
  -subset                control of nested looping

Exploration Constructs

  Object           DateCreated
  Parent/Child     Book/AuthorList
  Heterogeneous    "PubmedArticleSet/*"
  Nested           "*/Taxon"
  Recursive        "**/Gene-commentary"

Conditional Execution

  -if              Element [@attribute] required
  -unless          Skip if element matches
  -and             All tests must pass
  -or              Any passing test suffices
  -else            Execute if conditional test failed
  -position        Must be at [first|last] location in list

String Constraints

  -equals          String must match exactly
  -contains        Substring must be present
  -starts-with     Substring must be at beginning
  -ends-with       Substring must be at end
  -is-not          String must not match

Numeric Constraints

  -gt              Greater than
  -ge              Greater than or equal to
  -lt              Less than
  -le              Less than or equal to
  -eq              Equal to
  -ne              Not equal to

Format Customization

  -ret             Override line break between patterns
  -tab             Replace tab character between fields
  -sep             Separator between group members
  -pfx             Prefix to print before group
  -sfx             Suffix to print after group
  -clr             Clear queued tab separator
  -pfc             Preface combines -clr and -pfx
  -rst             Reset -sep, -pfx, and -sfx
  -def             Default placeholder for missing fields
  -lbl             Insert arbitrary text

Element Selection

  -element         Print all items that match tag name
  -first           Only print value of first item
  -last            Only print value of last item
  -encode          URL-encode <, >, &, ", and ' characters
  -upper           Convert text to upper-case
  -lower           Convert text to lower-case
  -NAME            Record value in named variable

-element Constructs

  Tag              Caption
  Group            Initials,LastName
  Parent/Child     MedlineCitation/PMID
  Attribute        DescriptorName@MajorTopicYN
  Recursive        "**/Gene-commentary_accession"
  Object Count     "#Author"
  Item Length      "%Title"
  Element Depth    "^PMID"
  Parent Index     "+"
  XML Subtree      "*"
  Variable         "&NAME"

Numeric Processing

  -num             Count
  -len             Length
  -sum             Sum
  -min             Minimum
  -max             Maximum
  -inc             Increment
  -dec             Decrement
  -sub             Difference
  -avg             Average
  -dev             Deviation

Phrase Processing

  -terms           Partition phrase at spaces
  -words           Split at punctuation marks
  -pairs           Adjacent informative words
  -phrase          Experimental index generation

Sequence Coordinates

  -0-based         Zero-Based
  -1-based         One-Based
  -ucsc            Half-Open

Command Generator

  -insd            Generate INSDSeq extraction commands

-insd Argument Order

  Descriptors      INSDSeq_sequence INSDSeq_definition INSDSeq_division
  Flags            complete or partial [optional]
  Feature(s)       CDS,mRNA
  Qualifiers       INSDFeature_key "#INSDInterval" gene product

Miscellaneous

  -head            Print before everything else
  -tail            Print after everything else
  -hd              Print before each record
  -tl              Print after each record

Reformatting

  -format          [compact|indent|expand]

Modification

  -filter          Object
                     [retain|remove|encode|decode|shrink]
                       [content|cdata|comment|object|attributes]

Validation

  -verify          Report XML data integrity problems

Summary

  -outline         Display outline of XML structure
  -synopsis        Display count of unique XML paths

Documentation

  -help            Print this document
  -examples        Examples of EDirect and xtract usage
  -version         Print version number

Notes

  String constraints use case-insensitive comparisons.

  Numeric constraints and selection arguments use integer values.

  -num and -len selections are synonyms for Object Count (#) and Item Length (%).

  -words, -pairs, and -phrase convert to lower case.

Examples

  -pattern DocumentSummary -element Id -first Name Title

  -pattern "PubmedArticleSet/*" -block Author -sep " " -element Initials,LastName

  -pattern PubmedArticle -block MeshHeading -if "@MajorTopicYN" -equals Y -sep " / " -element DescriptorName,QualifierName

  -pattern GenomicInfoType -element ChrAccVer ChrStart ChrStop

  -pattern Taxon -block "*/Taxon" -unless Rank -equals "no rank" -tab "\n" -element Rank,ScientificName

  -pattern Entrezgene -block "**/Gene-commentary"

  -block INSDReference -position 2

  -if Author -and Title

  -if "#Author" -lt 6 -and "%Title" -le 70

  -if DateCreated/Year -gt 2005

  -if ChrStop -lt ChrStart

  -if CommonName -contains mouse

  -if "&ABST" -starts-with "Transposable elements"

  -if MapLocation -element MapLocation -else -lbl "\-"

  -min ChrStart,ChrStop

  -max ExonCount

  -inc @aaPosition -element @residue

  -1-based ChrStart

  -insd CDS gene product protein_id translation

  -insd complete mat_peptide "%peptide" product peptide

  -insd CDS INSDInterval_iscomp@value INSDInterval_from INSDInterval_to

  -filter ExpXml decode content

  -filter LocationHist remove object
`

const xtractInternal = `
ReadBlocks -> SplitPattern => StreamTokens => ParseXML => ProcessQuery -> MergeResults

Performance Default Overrides

  -proc     Number of CPU processors used
  -cons     Ratio of parsers to processors
  -serv     Concurrent parser instances
  -chan     Communication channel depth
  -heap     Order restoration heap size
  -farm     Node allocation buffer length
  -gogc     Garbage collection tuning knob

Debugging

  -debug    Display run-time parameter summary
  -empty    Flag records with no output
  -index    Print record index numbers
  -stats    Show processing time for each record
  -timer    Report processing duration and rate
  -trial    Optimize -proc value, requires -input

Internal Component Performance

  -chunk    ReadBlocks
  -split    ReadBlocks -> SplitPattern
  -drain    ReadBlocks -> SplitPattern -> ConcurrencyChannel
  -token    ReadBlocks -> StreamTokens

Documentation

  -keys     Keyboard navigation shortcuts
  -unix     Common Unix commands

Sample File Download

  ftp-cp ftp.ncbi.nlm.nih.gov /entrez/entrezdirect/samples carotene.xml.zip
  unzip carotene.xml.zip
  rm carotene.xml.zip

Performance Tuning Script

  XtractTrials() {
    echo -e "<Trials>"
    for tries in {1..5}
    do
      xtract -debug -input "$1" -proc "$2" -pattern PubmedArticle -element LastName
    done
    echo -e "</Trials>"
  }

  for proc in {1..8}
  do
    XtractTrials "carotene.xml" "$proc" |
    xtract -pattern Trials -lbl "$proc" -avg Rate -dev Rate
  done

Processor Titration Results

  1    10455    124
  2    18674    684
  3    23654    1371
  4    30527    521
  5    35132    2349
  6    40289    907
  7    45369    1370
  8    47397    1141

Execution Profiling

  xtract -profile -input carotene.xml -pattern PubmedArticle -element LastName
  go tool pprof --pdf ./xtract ./cpu.pprof > ./callgraph.pdf
`

const xtractExamples = `
Author Frequency

  esearch -db pubmed -query "rattlesnake phospholipase" |
  efetch -format docsum |
  xtract -pattern DocumentSummary -sep "\n" -element Name |
  sort-uniq-count-rank

  39    Marangoni S
  31    Toyama MH
  26    Soares AM
  25    Bon C
  ...

Publications

  efetch -db pubmed -id 6271474,5685784,4882854,6243420 -format xml |
  xtract -pattern PubmedArticle -element MedlineCitation/PMID "#Author" \
    -block Author -position first -sep " " -element Initials,LastName \
    -block Article -element ArticleTitle

  6271474    5    MJ Casadaban     Tn3: transposition and control.
  5685784    2    RK Mortimer      Suppressors and suppressible mutations in yeast.
  4882854    2    ED Garber        Proteins and enzymes as taxonomic tools.
  6243420    1    NR Cozzarelli    DNA gyrase and the supercoiling of DNA.

Formatted Authors

  efetch -db pubmed -id 1413997,6301692,781293 -format xml |
  xtract -pattern PubmedArticle -element MedlineCitation/PMID \
    -block DateCreated -sep "-" -element Year,Month,Day \
    -block Author -sep " " -tab "" \
      -element "&COM" Initials,LastName -COM "(, )"

  1413997    1992-11-25    RK Mortimer, CR Contopoulou, JS King
  6301692    1983-06-17    MA Krasnow, NR Cozzarelli
  781293     1976-10-02    MJ Casadaban

Medical Subject Headings

  efetch -db pubmed -id 6092233,2539356,1937004 -format xml |
  xtract -pattern PubmedArticle -element MedlineCitation/PMID \
    -block MeshHeading \
      -subset DescriptorName -pfc "\n" -sep "|" -element @MajorTopicYN,DescriptorName \
      -subset QualifierName -pfc " / " -sep "|" -element @MajorTopicYN,QualifierName |
  sed -e 's/N|//g' -e 's/Y|/*/g'

  6092233
  Base Sequence
  DNA Restriction Enzymes
  DNA, Fungal / genetics / *isolation & purification
  *Genes, Fungal
  ...

Peptide Sequences

  esearch -db protein -query "conotoxin AND mat_peptide [FKEY]" |
  efetch -format gpc |
  xtract -insd complete mat_peptide "%peptide" product peptide |
  grep -i conotoxin | sort -t $'\t' -u -k 2,2n | head -n 8

  ADB43131.1    15    conotoxin Cal 1b     LCCKRHHGCHPCGRT
  AIC77099.1    16    conotoxin Im1.2      GCCSHPACNVNNPHIC
  AIC77105.1    17    conotoxin Lt1.4      GCCSHPACDVNNPDICG
  AIC77103.1    18    conotoxin Lt1.2      PRCCSNPACNANHAEICG
  AIC77083.1    20    conotoxin Bt14.6     KDCTYCMHSSCSMMYEKCRP
  AIC77085.1    21    conotoxin Bt14.8     NECDNCMRSFCSMIYEKCRLK
  AIC77093.1    22    conotoxin Bt14.16    GDCKPCMHPDCRFNPGRCRPRE
  AIC77154.1    23    conotoxin Bt14.19    VREKDCPPHPVPGMHKCVCLKTC

Chromosome Locations

  esearch -db gene -query "calmodulin [PFN] AND mammalia [ORGN]" |
  efetch -format docsum |
  xtract -pattern DocumentSummary -MAP "(-)" -MAP MapLocation \
    -element Id Name "&MAP" ScientificName

  801       CALM1    14q32.11         Homo sapiens
  808       CALM3    19q13.2-q13.3    Homo sapiens
  805       CALM2    2p21             Homo sapiens
  24242     Calm1    6q31-q32         Rattus norvegicus
  12313     Calm1    12 E             Mus musculus
  326597    CALM     -                Bos taurus
  50663     Calm2    6q11-q12         Rattus norvegicus
  24244     Calm3    1q22             Rattus norvegicus
  12315     Calm3    7 9.15 cM        Mus musculus
  12314     Calm2    17 E4            Mus musculus
  617095    CALM1    -                Bos taurus
  396838    CALM3    6                Sus scrofa
  ...

Gene Regions

  esearch -db gene -query "DDT [GENE] AND mouse [ORGN]" |
  efetch -format docsum |
  xtract -pattern GenomicInfoType -element ChrAccVer ChrStart ChrStop |
  xargs -n 3 sh -c 'efetch -db nuccore -format gb \
    -id "$0" -chr_start "$1" -chr_stop "$2"'

  LOCUS       NC_000076               2142 bp    DNA     linear   CON 09-FEB-2015
  DEFINITION  Mus musculus strain C57BL/6J chromosome 10, GRCm38.p3 C57BL/6J.
  ACCESSION   NC_000076 REGION: complement(75771233..75773374) GPC_000000783
  VERSION     NC_000076.6  GI:372099100
  ...
  FEATURES             Location/Qualifiers
       source          1..2142
                       /organism="Mus musculus"
                       /mol_type="genomic DNA"
                       /strain="C57BL/6J"
                       /db_xref="taxon:10090"
                       /chromosome="10"
       gene            1..2142
                       /gene="Ddt"
       mRNA            join(1..159,462..637,1869..2142)
                       /gene="Ddt"
                       /product="D-dopachrome tautomerase"
                       /transcript_id="NM_010027.1"
       CDS             join(52..159,462..637,1869..1941)
                       /gene="Ddt"
                       /codon_start=1
                       /product="D-dopachrome decarboxylase"
                       /protein_id="NP_034157.1"
                       /translation="MPFVELETNLPASRIPAGLENRLCAATATILDKPEDRVSVTIRP
                       GMTLLMNKSTEPCAHLLVSSIGVVGTAEQNRTHSASFFKFLTEELSLDQDRIVIRFFP
                       ...

Taxonomic Names

  esearch -db taxonomy -query "txid10090 [SBTR] OR camel [COMN]" |
  efetch -format docsum |
  xtract -pattern DocumentSummary -if CommonName \
    -element Id ScientificName CommonName

  57486    Mus musculus molossinus    Japanese wild mouse
  39442    Mus musculus musculus      eastern European house mouse
  35531    Mus musculus bactrianus    southwestern Asian house mouse
  10092    Mus musculus domesticus    western European house mouse
  10091    Mus musculus castaneus     southeastern Asian house mouse
  10090    Mus musculus               house mouse
  9838     Camelus dromedarius        Arabian camel
  9837     Camelus bactrianus         Bactrian camel

Structural Similarity

  esearch -db structure -query "crotalus [ORGN] AND phospholipase A2" |
  elink -related |
  efilter -query "archaea [ORGN]" |
  efetch -format docsum |
  xtract -pattern DocumentSummary \
    -if PdbClass -equals Hydrolase \
      -element PdbAcc PdbDescr

  3VV2    Crystal Structure Of Complex Form Between S324a-subtilisin And Mutant Tkpro
  3VHQ    Crystal Structure Of The Ca6 Site Mutant Of Pro-Sa-Subtilisin
  2ZWP    Crystal Structure Of Ca3 Site Mutant Of Pro-S324a
  2ZWO    Crystal Structure Of Ca2 Site Mutant Of Pro-S324a
  ...

Multiple Links

  esearch -db pubmed -query "conotoxin AND dopamine [MAJR]" |
  elink -target protein -cmd neighbor |
  xtract -pattern LinkSet -if Link/Id -element IdList/Id Link/Id

  23624852    17105332
  14657161    27532980    27532978
  12944511    31542395
  11222635    144922602

Gene Comments

  esearch -db gene -query "rbcL [GENE] AND maize [ORGN]" |
  efetch -format xml |
  xtract -pattern Entrezgene -block "**/Gene-commentary" \
    -if Gene-commentary_type@value -equals genomic \
      -tab "\n" -element Gene-commentary_accession |
  sort | uniq

  NC_001666
  X86563
  Z11973

Vitamin Biosynthesis

  esearch -db pubmed -query "tomato lycopene cyclase" |
  elink -related |
  elink -target protein |
  efilter -organism mammals |
  efetch -format gpc |
  xtract -pattern INSDSeq -if INSDSeq_definition -contains carotene \
    -element INSDSeq_accession-version INSDSeq_definition

  NP_573480.1       beta,beta-carotene 9',10'-oxygenase [Mus musculus]
  NP_001156500.1    beta,beta-carotene 15,15'-dioxygenase isoform 2 [Mus musculus]
  NP_067461.2       beta,beta-carotene 15,15'-dioxygenase isoform 1 [Mus musculus]
  NP_001297121.1    beta-carotene oxygenase 2 [Mustela putorius furo]
  AAS20392.1        carotene-9',10'-monooxygenase [Mustela putorius furo]

Indexed Fields

  einfo -db pubmed |
  xtract -pattern Field \
    -if IsDate -equals Y -and IsHidden -equals N \
      -pfx "[" -sep "]\t" -element Name,FullName |
  sort -t $'\t' -k 2f

  [CDAT]    Date - Completion
  [CRDT]    Date - Create
  [EDAT]    Date - Entrez
  [MHDA]    Date - MeSH
  [MDAT]    Date - Modification
  [PDAT]    Date - Publication

Author Numbers

  esearch -db pubmed -query "conotoxin" |
  efetch -format xml |
  xtract -pattern PubmedArticle -num Author |
  sort-uniq-count -n |
  reorder-columns 2 1 |
  head -n 15 |
  xy-plot auth.png

  0     11
  1     193
  2     854
  3     844
  4     699
  5     588
  6     439
  7     291
  8     187
  9     124
  10    122
  11    58
  12    33
  13    18

  900 +
      |           ********
  800 +           *       **
      |          *          *
  700 +          *          ***
      |          *             **
  600 +         *                *
      |         *                ***
  500 +         *                   **
      |        *                      ***
  400 +       *                          **
      |       *                            *
  300 +       *                            ***
      |      *                                *
  200 +      *                                 ******
      |     *                                        *********
  100 +   **                                                  *
      |  *                                                     **********
    0 + *                                                                ******
        +---------+---------+---------+---------+---------+---------+---------+
        0         2         4         6         8        10        12        14

Record Counts

  echo "diphtheria measles pertussis polio tuberculosis" |
  xargs -n 1 sh -c 'esearch -db pubmed -query "$0 [MESH]" |
  efilter -days 365 -datetype PDAT |
  xtract -pattern ENTREZ_DIRECT -lbl "$0" -element Count'

  diphtheria      18
  measles         166
  pertussis       98
  polio           75
  tuberculosis    1386

Gene Products

  for sym in HBB DMD TTN ATP7B HFE BRCA2 CFTR PAH PRNP RAG1
  do
    esearch -db gene -query "$sym [GENE] AND human [ORGN]" |
    efilter -query "alive [PROP]" | efetch -format docsum |
    xtract -pattern GenomicInfoType \
      -element ChrAccVer ChrStart ChrStop |
    while read acc str stp
    do
      efetch -db nuccore -format gbc \
        -id "$acc" -chr_start "$str" -chr_stop "$stp" |
      xtract -insd CDS,mRNA INSDFeature_key "#INSDInterval" \
        gene "%transcription" "%translation" \
        product transcription translation |
      grep -i $'\t'"$sym"$'\t'
    done
  done

  NC_000011.10    mRNA    3     HBB    626      hemoglobin, beta                     ACATTTGCTT...
  NC_000011.10    CDS     3     HBB    147      hemoglobin subunit beta              MVHLTPEEKS...
  NC_000023.11    mRNA    78    DMD    13805    dystrophin, transcript variant X2    AGGAAGATGA...
  NC_000023.11    mRNA    77    DMD    13794    dystrophin, transcript variant X6    ACTTTCCCCC...
  NC_000023.11    mRNA    77    DMD    13800    dystrophin, transcript variant X5    ACTTTCCCCC...
  NC_000023.11    mRNA    77    DMD    13785    dystrophin, transcript variant X7    ACTTTCCCCC...
  NC_000023.11    mRNA    74    DMD    13593    dystrophin, transcript variant X8    ACTTTCCCCC...
  NC_000023.11    mRNA    75    DMD    13625    dystrophin, transcript variant X9    ACTTTCCCCC...
  ...

Genome Range

  esearch -db gene -query "Homo sapiens [ORGN] AND Y [CHR]" |
  efilter -status alive | efetch -format docsum |
  xtract -pattern DocumentSummary -NAME Name -DESC Description \
    -block GenomicInfoType -if ChrLoc -equals Y \
      -min ChrStart,ChrStop -element "&NAME" "&DESC" |
  sort -k 1,1n | cut -f 2- |
  between-two-genes ASMT IL3RA

  IL3RA           interleukin 3 receptor subunit alpha
  LOC101928032    uncharacterized LOC101928032
  LOC101928055    uncharacterized LOC101928055
  SLC25A6         solute carrier family 25 member 6
  LOC105373102    uncharacterized LOC105373102
  LINC00106       long intergenic non-protein coding RNA 106
  ASMTL-AS1       ASMTL antisense RNA 1
  ASMTL           acetylserotonin O-methyltransferase-like
  P2RY8           purinergic receptor P2Y8
  AKAP17A         A-kinase anchoring protein 17A
  ASMT            acetylserotonin O-methyltransferase

Amino Acid Substitutions

  ApplySNPs() {
    seq=""
    last=""

    while read rsid accn pos res
    do
      if [ "$accn" != "$last" ]
      then
        insd=$(efetch -db protein -id "$accn" -format gbc < /dev/null)
        seq=$(echo $insd | xtract -pattern INSDSeq -element INSDSeq_sequence)
        last=$accn
      fi

      pos=$((pos+1))
      pfx=""
      sfx=""
      echo ">rs$rsid [$accn $res@$pos]"
      if [ $pos -gt 1 ]
      then
        pfx=$(echo ${seq:0:$pos-1})
      fi
      if [ $pos -lt ${#seq} ]
      then
        sfx=$(echo ${seq:$pos})
      fi
      echo "$pfx$res$sfx" | fold -w 50
    done
  }

  esearch -db gene -query "OPN1MW [GENE] AND human [ORGN]" |
  elink -target snp | efetch -format xml |
  xtract -pattern Rs -RSID Rs@rsId \
    -block FxnSet -if @fxnClass -equals missense \
      -sep "." -element "&RSID" @protAcc,@protVer @aaPosition \
      -tab "\n" -element @residue |
  sort -t $'\t' -k 2,2 -k 3,3n -k 4,4 | uniq |
  ApplySNPs

  >rs104894915 [NP_000504.1 K@94]
  maqqwslqrlagrhpqdsyedstqssiftytnsnstrgpfegpnyhiapr
  wvyhltsvwmifvviasvftnglvlaatmkfkklrhplnwilvKlavadl
  aetviastisvvnqvygyfvlghpmcvlegytvslcgitglwslaiiswe
  ...

Amino Acid Composition

  #!/bin/bash -norc

  abbrev=( Ala Asx Cys Asp Glu Phe Gly His Ile \
           Xle Lys Leu Met Asn Pyl Pro Gln Arg \
           Ser Thr Sec Val Trp Xxx Tyr Glx )

  AminoAcidComp() {
    local count
    while read num lttr
    do
      idx=$(printf %i "'$lttr'")
      ofs=$((idx-97))
      count[$ofs]="$num"
    done <<< "$1"
    for i in {0..25}
    do
      echo -e "${abbrev[$i]}\t${count[$i]-0}"
    done |
    sort
  }

  AminoAcidJoin() {
    result=""
    while read acc seq gene
    do
      comp="$(echo "$seq" | tr A-Z a-z | sed 's/[^a-z]//g' | fold -w 1 | sort-uniq-count)"
      current=$(AminoAcidComp "$comp")
      current=$(echo -e "GENE\t$gene\n$current")
      if [ -n "$result" ]
      then
        result=$(join -t $'\t' <(echo "$result") <(echo "$current"))
      else
        result=$current
      fi
    done
    echo "$result" |
    grep -e "GENE" -e "[1-9]"
  }

  ids="NP_001172026,NP_000509,NP_004001,NP_001243779"
  efetch -db protein -id "$ids" -format gpc |
  xtract -insd INSDSeq_sequence CDS gene |
  AminoAcidJoin

  GENE    INS    HBB    DMD    TTN
  Ala     10     15     210    2084
  Arg     5      3      193    1640
  Asn     3      6      153    1111
  Asp     2      7      185    1720
  Cys     6      2      35     513
  Gln     7      3      301    942
  Glu     8      8      379    3193
  Gly     12     13     104    2066
  His     2      9      84     478
  Ile     2      0      165    2062
  Leu     20     18     438    2117
  Lys     2      11     282    2943
  Met     2      2      79     398
  Phe     3      8      77     908
  Pro     6      7      130    2517
  Ser     5      5      239    2463
  Thr     3      7      194    2546
  Trp     2      2      67     466
  Tyr     4      3      61     999
  Val     6      18     186    3184

Word Pair Indexing

  efetch -db pubmed -id 2981625,12857958 -format xml |
  xtract -head "<Set>" -tail "</Set>" -hd "<Rec>" -tl "</Rec>" \
    -pattern PubmedArticle -if ArticleTitle \
      -pfx "<Id>" -sfx "</Id>" -element MedlineCitation/PMID \
      -pfc "<Pairs><Term>" -sfx "</Term></Pairs>" \
      -sep "</Term><Term>" -pairs ArticleTitle |
  xtract -pattern Rec -UID Id -block Term -pfc "\n" -element "&UID",Term

  2981625     recombination site
  2981625     site selection
  2981625     tn3 resolvase
  2981625     resolvase topological
  2981625     topological tests
  2981625     tracking mechanism
  12857958    chirality sensing
  12857958    escherichia coli
  12857958    coli topoisomerase
  12857958    topoisomerase iv
  12857958    type ii
  12857958    ii topoisomerases

Phrase Searching

  entrez-phrase-search -db pubmed -field WORD \
    selective serotonin reuptake inhibitor + monoamine oxidase inhibitor |
  efetch -format xml |
  xtract -pattern PubmedArticle -element MedlineCitation/PMID \
    -block Keyword -pfc "\n  " -element Keyword

  24657329
    Antidepressant
    Organic cation transporter 2
    Piperine
    Uptake 2
  24280122
    5-HIAA
    5-HT
    5-HTP
    5-hydroxyindoleacetic acid
    5-hydroxytryptophan
    ...
`

const pubMedArtSample = `
<PubmedArticle>
<MedlineCitation Owner="NLM" Status="MEDLINE">
<PMID Version="1">6301692</PMID>
<DateCreated>
<Year>1983</Year>
<Month>06</Month>
<Day>17</Day>
</DateCreated>
<DateCompleted>
<Year>1983</Year>
<Month>06</Month>
<Day>17</Day>
</DateCompleted>
<DateRevised>
<Year>2007</Year>
<Month>11</Month>
<Day>14</Day>
</DateRevised>
<Article PubModel="Print">
<Journal>
<ISSN IssnType="Print">0092-8674</ISSN>
<JournalIssue CitedMedium="Print">
<Volume>32</Volume>
<Issue>4</Issue>
<PubDate>
<Year>1983</Year>
<Month>Apr</Month>
</PubDate>
</JournalIssue>
<Title>Cell</Title>
<ISOAbbreviation>Cell</ISOAbbreviation>
</Journal>
<ArticleTitle>Site-specific relaxation and recombination by the Tn3 resolvase: recognition of the DNA path between oriented res sites.</ArticleTitle>
<Pagination>
<MedlinePgn>1313-24</MedlinePgn>
</Pagination>
<Abstract>
<AbstractText Label="RESULTS>We studied the dynamics of site-specific recombination by the resolvase encoded by the Escherichia coli transposon Tn3.
The pure enzyme recombined supercoiled plasmids containing two directly repeated recombination sites, called res sites.
Resolvase is the first strictly site-specific topoisomerase.
It relaxed only plasmids containing directly repeated res sites; substrates with zero, one or two inverted sites were inert.
Even when the proximity of res sites was ensured by catenation of plasmids with a single site, neither relaxation nor recombination occurred.
The two circular products of recombination were catenanes interlinked only once.
These properties of resolvase require that the path of the DNA between res sites be clearly defined and that strand exchange occur with a unique geometry.</AbstractText>
<AbstractText Label="SUMMARY">A model in which one subunit of a dimeric resolvase is bound at one res site,
while the other searches along adjacent DNA until it encounters the second site,
would account for the ability of resolvase to distinguish intramolecular from intermolecular sites,
to sense the relative orientation of sites and to produce singly interlinked catenanes.
Because resolvase is a type 1 topoisomerase, we infer that it makes the required duplex bDNA breaks of recombination one strand at a time.</AbstractText>
</Abstract>
<AuthorList CompleteYN="Y">
<Author ValidYN="Y">
<LastName>Krasnow</LastName>
<ForeName>Mark A</ForeName>
<Initials>MA</Initials>
</Author>
<Author ValidYN="Y">
<LastName>Cozzarelli</LastName>
<ForeName>Nicholas R</ForeName>
<Initials>NR</Initials>
</Author>
</AuthorList>
<Language>eng</Language>
<GrantList CompleteYN="Y">
<Grant>
<GrantID>GM-07281</GrantID>
<Acronym>GM</Acronym>
<Agency>NIGMS NIH HHS</Agency>
<Country>United States</Country>
</Grant>
</GrantList>
<PublicationTypeList>
<PublicationType UI="D016428">Journal Article</PublicationType>
<PublicationType UI="D013487">Research Support, U.S. Gov't, P.H.S.</PublicationType>
</PublicationTypeList>
</Article>
<MedlineJournalInfo>
<Country>UNITED STATES</Country>
<MedlineTA>Cell</MedlineTA>
<NlmUniqueID>0413066</NlmUniqueID>
<ISSNLinking>0092-8674</ISSNLinking>
</MedlineJournalInfo>
<ChemicalList>
<Chemical>
<RegistryNumber>0</RegistryNumber>
<NameOfSubstance UI="D004269">DNA, Bacterial</NameOfSubstance>
</Chemical>
<Chemical>
<RegistryNumber>0</RegistryNumber>
<NameOfSubstance UI="D004278">DNA, Superhelical</NameOfSubstance>
</Chemical>
<Chemical>
<RegistryNumber>0</RegistryNumber>
<NameOfSubstance UI="D004279">DNA, Viral</NameOfSubstance>
</Chemical>
<Chemical>
<RegistryNumber>EC 2.7.7.-</RegistryNumber>
<NameOfSubstance UI="D009713">Nucleotidyltransferases</NameOfSubstance>
</Chemical>
<Chemical>
<RegistryNumber>EC 2.7.7.-</RegistryNumber>
<NameOfSubstance UI="D019895">Transposases</NameOfSubstance>
</Chemical>
<Chemical>
<RegistryNumber>EC 5.99.1.2</RegistryNumber>
<NameOfSubstance UI="D004264">DNA Topoisomerases, Type I</NameOfSubstance>
</Chemical>
</ChemicalList>
<CitationSubset>IM</CitationSubset>
<MeshHeadingList>
<MeshHeading>
<DescriptorName MajorTopicYN="N" UI="D004264">DNA Topoisomerases, Type I</DescriptorName>
<QualifierName MajorTopicYN="N" UI="Q000378">metabolism</QualifierName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N" UI="D004269">DNA, Bacterial</DescriptorName>
<QualifierName MajorTopicYN="Y" UI="Q000378">metabolism</QualifierName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N" UI="D004278">DNA, Superhelical</DescriptorName>
<QualifierName MajorTopicYN="N" UI="Q000378">metabolism</QualifierName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N" UI="D004279">DNA, Viral</DescriptorName>
<QualifierName MajorTopicYN="Y" UI="Q000378">metabolism</QualifierName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N" UI="D008957">Models, Genetic</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="Y" UI="D009690">Nucleic Acid Conformation</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N" UI="D009713">Nucleotidyltransferases</DescriptorName>
<QualifierName MajorTopicYN="N" UI="Q000302">isolation &amp; purification</QualifierName>
<QualifierName MajorTopicYN="Y" UI="Q000378">metabolism</QualifierName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N" UI="D010957">Plasmids</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="Y" UI="D011995">Recombination, Genetic</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N" UI="D012091">Repetitive Sequences, Nucleic Acid</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N" UI="D013539">Simian virus 40</DescriptorName>
</MeshHeading>
<MeshHeading>
<DescriptorName MajorTopicYN="N" UI="D019895">Transposases</DescriptorName>
</MeshHeading>
</MeshHeadingList>
</MedlineCitation>
<PubmedData>
<History>
<PubMedPubDate PubStatus="pubmed">
<Year>1983</Year>
<Month>4</Month>
<Day>1</Day>
</PubMedPubDate>
<PubMedPubDate PubStatus="medline">
<Year>1983</Year>
<Month>4</Month>
<Day>1</Day>
<Hour>0</Hour>
<Minute>1</Minute>
</PubMedPubDate>
<PubMedPubDate PubStatus="entrez">
<Year>1983</Year>
<Month>4</Month>
<Day>1</Day>
<Hour>0</Hour>
<Minute>0</Minute>
</PubMedPubDate>
</History>
<PublicationStatus>ppublish</PublicationStatus>
<ArticleIdList>
<ArticleId IdType="pubmed">6301692</ArticleId>
<ArticleId IdType="pii">0092-8674(83)90312-4</ArticleId>
</ArticleIdList>
</PubmedData>
</PubmedArticle>
`

const insdSeqSample = `
<INSDSeq>
<INSDSeq_locus>AF480315_1</INSDSeq_locus>
<INSDSeq_length>67</INSDSeq_length>
<INSDSeq_moltype>AA</INSDSeq_moltype>
<INSDSeq_topology>linear</INSDSeq_topology>
<INSDSeq_division>INV</INSDSeq_division>
<INSDSeq_update-date>25-JUL-2016</INSDSeq_update-date>
<INSDSeq_create-date>31-DEC-2003</INSDSeq_create-date>
<INSDSeq_definition>four-loop conotoxin preproprotein, partial [Conus purpurascens]</INSDSeq_definition>
<INSDSeq_primary-accession>AAQ05867</INSDSeq_primary-accession>
<INSDSeq_accession-version>AAQ05867.1</INSDSeq_accession-version>
<INSDSeq_other-seqids>
<INSDSeqid>gb|AAQ05867.1|AF480315_1</INSDSeqid>
<INSDSeqid>gi|33320307</INSDSeqid>
</INSDSeq_other-seqids>
<INSDSeq_source>Conus purpurascens</INSDSeq_source>
<INSDSeq_organism>Conus purpurascens</INSDSeq_organism>
<INSDSeq_taxonomy>Eukaryota; Metazoa; Lophotrochozoa; Mollusca; Gastropoda; Caenogastropoda; Hypsogastropoda; Neogastropoda; Conoidea; Conidae; Conus</INSDSeq_taxonomy>
<INSDSeq_references>
<INSDReference>
<INSDReference_reference>1</INSDReference_reference>
<INSDReference_position>1..67</INSDReference_position>
<INSDReference_authors>
<INSDAuthor>Duda,T.F. Jr.</INSDAuthor>
<INSDAuthor>Palumbi,S.R.</INSDAuthor>
</INSDReference_authors>
<INSDReference_title>Convergent evolution of venoms and feeding ecologies among polyphyletic piscivorous Conus species</INSDReference_title>
<INSDReference_journal>Unpublished</INSDReference_journal>
</INSDReference>
<INSDReference>
<INSDReference_reference>2</INSDReference_reference>
<INSDReference_position>1..67</INSDReference_position>
<INSDReference_authors>
<INSDAuthor>Duda,T.F. Jr.</INSDAuthor>
<INSDAuthor>Palumbi,S.R.</INSDAuthor>
</INSDReference_authors>
<INSDReference_title>Direct Submission</INSDReference_title>
<INSDReference_journal>Submitted (04-FEB-2002) Naos Marine Lab, Smithsonian Tropical Research Institute, Apartado 2072, Balboa, Ancon, Panama, Republic of Panama</INSDReference_journal>
</INSDReference>
</INSDSeq_references>
<INSDSeq_comment>Method: conceptual translation supplied by author.</INSDSeq_comment>
<INSDSeq_source-db>accession AF480315.1</INSDSeq_source-db>
<INSDSeq_feature-table>
<INSDFeature>
<INSDFeature_key>source</INSDFeature_key>
<INSDFeature_intervals>
<INSDInterval>
<INSDInterval_from>1</INSDInterval_from>
<INSDInterval_to>67</INSDInterval_to>
<INSDInterval_accession>AAQ05867.1</INSDInterval_accession>
</INSDInterval>
</INSDFeature_intervals>
<INSDFeature_quals>
<INSDQualifier>
<INSDQualifier_name>organism</INSDQualifier_name>
<INSDQualifier_value>Conus purpurascens</INSDQualifier_value>
</INSDQualifier>
<INSDQualifier>
<INSDQualifier_name>isolate</INSDQualifier_name>
<INSDQualifier_value>purpurascens-2c</INSDQualifier_value>
</INSDQualifier>
<INSDQualifier>
<INSDQualifier_name>db_xref</INSDQualifier_name>
<INSDQualifier_value>taxon:41690</INSDQualifier_value>
</INSDQualifier>
<INSDQualifier>
<INSDQualifier_name>clone_lib</INSDQualifier_name>
<INSDQualifier_value>venom duct cDNA library</INSDQualifier_value>
</INSDQualifier>
<INSDQualifier>
<INSDQualifier_name>country</INSDQualifier_name>
<INSDQualifier_value>Panama</INSDQualifier_value>
</INSDQualifier>
<INSDQualifier>
<INSDQualifier_name>note</INSDQualifier_name>
<INSDQualifier_value>isolated from the Bay of Panama</INSDQualifier_value>
</INSDQualifier>
</INSDFeature_quals>
</INSDFeature>
<INSDFeature>
<INSDFeature_key>Protein</INSDFeature_key>
<INSDFeature_location>&lt;1..67</INSDFeature_location>
<INSDFeature_intervals>
<INSDInterval>
<INSDInterval_from>1</INSDInterval_from>
<INSDInterval_to>67</INSDInterval_to>
<INSDInterval_accession>AAQ05867.1</INSDInterval_accession>
</INSDInterval>
</INSDFeature_intervals>
<INSDFeature_partial5 value="true"/>
<INSDFeature_quals>
<INSDQualifier>
<INSDQualifier_name>product</INSDQualifier_name>
<INSDQualifier_value>four-loop conotoxin preproprotein</INSDQualifier_value>
</INSDQualifier>
</INSDFeature_quals>
</INSDFeature>
<INSDFeature>
<INSDFeature_key>mat_peptide</INSDFeature_key>
<INSDFeature_location>41..67</INSDFeature_location>
<INSDFeature_intervals>
<INSDInterval>
<INSDInterval_from>41</INSDInterval_from>
<INSDInterval_to>67</INSDInterval_to>
<INSDInterval_accession>AAQ05867.1</INSDInterval_accession>
</INSDInterval>
</INSDFeature_intervals>
<INSDFeature_quals>
<INSDQualifier>
<INSDQualifier_name>product</INSDQualifier_name>
<INSDQualifier_value>four-loop conotoxin</INSDQualifier_value>
</INSDQualifier>
<INSDQualifier>
<INSDQualifier_name>calculated_mol_wt</INSDQualifier_name>
<INSDQualifier_value>3008</INSDQualifier_value>
</INSDQualifier>
<INSDQualifier>
<INSDQualifier_name>peptide</INSDQualifier_name>
<INSDQualifier_value>PCKKTGRKCFPHQKDCCGRACIITICP</INSDQualifier_value>
</INSDQualifier>
</INSDFeature_quals>
</INSDFeature>
<INSDFeature>
<INSDFeature_key>CDS</INSDFeature_key>
<INSDFeature_location>1..67</INSDFeature_location>
<INSDFeature_intervals>
<INSDInterval>
<INSDInterval_from>1</INSDInterval_from>
<INSDInterval_to>67</INSDInterval_to>
<INSDInterval_accession>AAQ05867.1</INSDInterval_accession>
</INSDInterval>
</INSDFeature_intervals>
<INSDFeature_partial5 value="true"/>
<INSDFeature_quals>
<INSDQualifier>
<INSDQualifier_name>coded_by</INSDQualifier_name>
<INSDQualifier_value>AF480315.1:&lt;1..205</INSDQualifier_value>
</INSDQualifier>
<INSDQualifier>
<INSDQualifier_name>codon_start</INSDQualifier_name>
<INSDQualifier_value>2</INSDQualifier_value>
</INSDQualifier>
</INSDFeature_quals>
</INSDFeature>
</INSDSeq_feature-table>
<INSDSeq_sequence>vvivavlfltacqlitaddsrrtqkhralrsttkratsnrpckktgrkcfphqkdccgraciiticp</INSDSeq_sequence>
</INSDSeq>
`

const geneDocSumSample = `
<DocumentSummary>
<Id>3581</Id>
<Name>IL9R</Name>
<Description>interleukin 9 receptor</Description>
<Status>0</Status>
<CurrentID>0</CurrentID>
<Chromosome>X, Y</Chromosome>
<GeneticSource>genomic</GeneticSource>
<MapLocation>Xq28 and Yq12</MapLocation>
<OtherAliases>CD129, IL-9R</OtherAliases>
<OtherDesignations>interleukin-9 receptor|IL-9 receptor</OtherDesignations>
<NomenclatureSymbol>IL9R</NomenclatureSymbol>
<NomenclatureName>interleukin 9 receptor</NomenclatureName>
<NomenclatureStatus>Official</NomenclatureStatus>
<Mim>
<int>300007</int>
</Mim>
<GenomicInfo>
<GenomicInfoType>
<ChrLoc>X</ChrLoc>
<ChrAccVer>NC_000023.11</ChrAccVer>
<ChrStart>155997580</ChrStart>
<ChrStop>156013016</ChrStop>
<ExonCount>14</ExonCount>
</GenomicInfoType>
<GenomicInfoType>
<ChrLoc>Y</ChrLoc>
<ChrAccVer>NC_000024.10</ChrAccVer>
<ChrStart>57184100</ChrStart>
<ChrStop>57199536</ChrStop>
<ExonCount>14</ExonCount>
</GenomicInfoType>
</GenomicInfo>
<GeneWeight>5425</GeneWeight>
<Summary>The protein encoded by this gene is a cytokine receptor that specifically mediates the biological effects of interleukin 9 (IL9).
The functional IL9 receptor complex requires this protein as well as the interleukin 2 receptor, gamma (IL2RG), a common gamma subunit shared by the receptors of many different cytokines.
The ligand binding of this receptor leads to the activation of various JAK kinases and STAT proteins, which connect to different biologic responses.
This gene is located at the pseudoautosomal regions of X and Y chromosomes.
Genetic studies suggested an association of this gene with the development of asthma.
Multiple pseudogenes on chromosome 9, 10, 16, and 18 have been described.
Alternatively spliced transcript variants have been found for this gene.</Summary>
<ChrSort>X</ChrSort>
<ChrStart>155997580</ChrStart>
<Organism>
<ScientificName>Homo sapiens</ScientificName>
<CommonName>human</CommonName>
<TaxID>9606</TaxID>
</Organism>
</DocumentSummary>
`

const keyboardShortcuts = `
Command History

  Ctrl-n     Next command
  Ctrl-p     Previous command

Move Cursor Forward

  Ctrl-e     To end of line
  Ctrl-f     By one character
  Esc-f      By one argument

Move Cursor Backward

  Ctrl-a     To beginning of line
  Ctrl-b     By one character
  Esc-b      By one argument

Delete

  Del        Previous character
  Ctrl-d     Next character
  Ctrl-k     To end of line
  Ctrl-u     Entire line
  Ctrl-w     Previous word
  Esc-Del    Previous argument
  Esc-d      Next argument

Autocomplete

  Tab        Completes directory or file names

Program Control

  Ctrl-c     Quit running program
  ^x^y       Run last command replacing x with y
`

const unixCommands = `
Process by Contents

 sort      Sorts lines of text

  -f       Ignore case
  -n       Numeric comparison
  -r       Reverse result order

  -k       Field key (start,stop or first)
  -u       Unique lines with identical keys

  -b       Ignore leading blanks
  -s       Stable sort
  -t       Specify field separator

 uniq      Removes repeated lines

  -c       Count occurrences
  -i       Ignore case

  -f       Ignore first n fields
  -s       Ignore first n characters

  -d       Only output repeated lines
  -u       Only output non-repeated lines

 grep      Matches patterns using regular expressions

  -i       Ignore case
  -v       Invert search
  -w       Search expression as a word
  -x       Search expression as whole line

  -e       Specify individual pattern

  -c       Only count number of matches
  -n       Print line numbers

Regular Expressions

 Characters

  .        Any single character (except newline)
  \w       Alphabetic [A-Za-z], numeric [0-9], or underscore (_)
  \s       Whitespace (space or tab)
  \        Escapes special characters
  []       Matches any enclosed characters

 Positions

  ^        Beginning of line
  $        End of line
  \b       Word boundary

 Repeat Matches

  ?        0 or 1
  *        0 or more
  +        1 or more
  {n}      Exactly n

 Escape Sequences

  \n       Line break
  \t       Tab character

Modify Contents

 sed       Replaces text strings

  -e       Specify individual expression

 tr        Translates characters

  -d       Delete character

 rev       Reverses characters on line

Format Contents

 column    Aligns columns by content width

  -s       Specify field separator
  -t       Create table

 expand    Aligns columns to specified positions

  -t       Tab positions

 fold      Wraps lines at a specific width

  -w       Line width

Filter by Position

 cut       Removes parts of lines

  -c       Characters to keep
  -f       Fields to keep
  -d       Specify field separator
  -s       Suppress lines with no delimiters

 head      Prints first lines

  -n       Number of lines

 tail      Prints last lines

  -n       Number of lines

Miscellaneous

 wc        Counts words, lines, or characters

  -c       Characters
  -l       Lines
  -w       Words

 xargs     Constructs arguments

  -n       Number of words per batch

File Compression

 tar       Archive files

  -c       Create archive
  -f       Name of output file
  -z       Compress archive with gzip

 gzip      Compress file

  -k       Keep original file

 unzip     Decompress .zip archive

  -p       Pipe to stdout

 gzcat     Decompress .gz archive and pipe to stdout

Directory and File Navigation

 cd        Changes directory

  /        Root
  ~        Home
  .        Current
  ..       Parent
  -        Previous

 ls        Lists file names

  -1       One entry per line
  -a       Show files beginning with dot (.)
  -l       List in long format
  -R       Recursively explore subdirectories
  -S       Sort files by size
  -t       Sort by most recently modified

 pwd       Prints working directory path
`

// TYPED CONSTANTS

type LevelType int

const (
	_ LevelType = iota
	UNIT
	SUBSET
	SECTION
	BLOCK
	BRANCH
	GROUP
	DIVISION
	PATTERN
)

type IndentType int

const (
	SINGULARITY IndentType = iota
	COMPACT
	FLUSH
	INDENT
	SUBTREE
	WRAPPED
)

type SideType int

const (
	_ SideType = iota
	LEFT
	RIGHT
)

type TagType int

const (
	NOTAG TagType = iota
	STARTTAG
	SELFTAG
	STOPTAG
	ATTRIBTAG
	CONTENTTAG
	CDATATAG
	COMMENTTAG
	OBJECTTAG
	ISCLOSED
)

type OpType int

const (
	UNSET OpType = iota
	ELEMENT
	FIRST
	LAST
	ENCODE
	UPPER
	LOWER
	TERMS
	WORDS
	PAIRS
	PHRASE
	PFX
	SFX
	SEP
	TAB
	RET
	LBL
	CLR
	PFC
	RST
	DEF
	POSITION
	IF
	UNLESS
	MATCH
	AVOID
	AND
	OR
	EQUALS
	CONTAINS
	STARTSWITH
	ENDSWITH
	ISNOT
	GT
	GE
	LT
	LE
	EQ
	NE
	NUM
	LEN
	SUM
	MIN
	MAX
	INC
	DEC
	SUB
	AVG
	DEV
	ZEROBASED
	ONEBASED
	UCSC
	ELSE
	VARIABLE
	VALUE
	STAR
	COUNT
	LENGTH
	DEPTH
	INDEX
	UNRECOGNIZED
)

type ArgumentType int

const (
	_ ArgumentType = iota
	EXPLORATION
	CONDITIONAL
	EXTRACTION
	CUSTOMIZATION
)

type SpecialType int

const (
	NOPROCESS SpecialType = iota
	DOFORMAT
	DOOUTLINE
	DOSYNOPSIS
	DOVERIFY
	DOFILTER
	DOCHUNK
	DOSPLIT
	DODRAIN
	DOTOKEN
)

type SeqEndType int

const (
	_ SeqEndType = iota
	ISSTART
	ISSTOP
	ISPOS
)

type SequenceType struct {
	Based int
	Which SeqEndType
}

// ARGUMENT MAPS

var argTypeIs = map[string]ArgumentType{
	"-unit":        EXPLORATION,
	"-Unit":        EXPLORATION,
	"-subset":      EXPLORATION,
	"-Subset":      EXPLORATION,
	"-section":     EXPLORATION,
	"-Section":     EXPLORATION,
	"-block":       EXPLORATION,
	"-Block":       EXPLORATION,
	"-branch":      EXPLORATION,
	"-Branch":      EXPLORATION,
	"-group":       EXPLORATION,
	"-Group":       EXPLORATION,
	"-division":    EXPLORATION,
	"-Division":    EXPLORATION,
	"-pattern":     EXPLORATION,
	"-Pattern":     EXPLORATION,
	"-position":    CONDITIONAL,
	"-if":          CONDITIONAL,
	"-unless":      CONDITIONAL,
	"-match":       CONDITIONAL,
	"-avoid":       CONDITIONAL,
	"-and":         CONDITIONAL,
	"-or":          CONDITIONAL,
	"-equals":      CONDITIONAL,
	"-contains":    CONDITIONAL,
	"-starts-with": CONDITIONAL,
	"-ends-with":   CONDITIONAL,
	"-is-not":      CONDITIONAL,
	"-gt":          CONDITIONAL,
	"-ge":          CONDITIONAL,
	"-lt":          CONDITIONAL,
	"-le":          CONDITIONAL,
	"-eq":          CONDITIONAL,
	"-ne":          CONDITIONAL,
	"-element":     EXTRACTION,
	"-first":       EXTRACTION,
	"-last":        EXTRACTION,
	"-encode":      EXTRACTION,
	"-upper":       EXTRACTION,
	"-lower":       EXTRACTION,
	"-terms":       EXTRACTION,
	"-words":       EXTRACTION,
	"-pairs":       EXTRACTION,
	"-phrase":      EXTRACTION,
	"-num":         EXTRACTION,
	"-len":         EXTRACTION,
	"-sum":         EXTRACTION,
	"-min":         EXTRACTION,
	"-max":         EXTRACTION,
	"-inc":         EXTRACTION,
	"-dec":         EXTRACTION,
	"-sub":         EXTRACTION,
	"-avg":         EXTRACTION,
	"-dev":         EXTRACTION,
	"-0-based":     EXTRACTION,
	"-zero-based":  EXTRACTION,
	"-1-based":     EXTRACTION,
	"-one-based":   EXTRACTION,
	"-ucsc":        EXTRACTION,
	"-else":        EXTRACTION,
	"-pfx":         CUSTOMIZATION,
	"-sfx":         CUSTOMIZATION,
	"-sep":         CUSTOMIZATION,
	"-tab":         CUSTOMIZATION,
	"-ret":         CUSTOMIZATION,
	"-lbl":         CUSTOMIZATION,
	"-clr":         CUSTOMIZATION,
	"-pfc":         CUSTOMIZATION,
	"-rst":         CUSTOMIZATION,
	"-def":         CUSTOMIZATION,
}

var opTypeIs = map[string]OpType{
	"-element":     ELEMENT,
	"-first":       FIRST,
	"-last":        LAST,
	"-encode":      ENCODE,
	"-upper":       UPPER,
	"-lower":       LOWER,
	"-terms":       TERMS,
	"-words":       WORDS,
	"-pairs":       PAIRS,
	"-phrase":      PHRASE,
	"-pfx":         PFX,
	"-sfx":         SFX,
	"-sep":         SEP,
	"-tab":         TAB,
	"-ret":         RET,
	"-lbl":         LBL,
	"-clr":         CLR,
	"-pfc":         PFC,
	"-rst":         RST,
	"-def":         DEF,
	"-position":    POSITION,
	"-if":          IF,
	"-unless":      UNLESS,
	"-match":       MATCH,
	"-avoid":       AVOID,
	"-and":         AND,
	"-or":          OR,
	"-equals":      EQUALS,
	"-contains":    CONTAINS,
	"-starts-with": STARTSWITH,
	"-ends-with":   ENDSWITH,
	"-is-not":      ISNOT,
	"-gt":          GT,
	"-ge":          GE,
	"-lt":          LT,
	"-le":          LE,
	"-eq":          EQ,
	"-ne":          NE,
	"-num":         NUM,
	"-len":         LEN,
	"-sum":         SUM,
	"-min":         MIN,
	"-max":         MAX,
	"-inc":         INC,
	"-dec":         DEC,
	"-sub":         SUB,
	"-avg":         AVG,
	"-dev":         DEV,
	"-0-based":     ZEROBASED,
	"-zero-based":  ZEROBASED,
	"-1-based":     ONEBASED,
	"-one-based":   ONEBASED,
	"-ucsc":        UCSC,
	"-else":        ELSE,
}

var levelTypeIs = map[string]LevelType{
	"-unit":     UNIT,
	"-Unit":     UNIT,
	"-subset":   SUBSET,
	"-Subset":   SUBSET,
	"-section":  SECTION,
	"-Section":  SECTION,
	"-block":    BLOCK,
	"-Block":    BLOCK,
	"-branch":   BRANCH,
	"-Branch":   BRANCH,
	"-group":    GROUP,
	"-Group":    GROUP,
	"-division": DIVISION,
	"-Division": DIVISION,
	"-pattern":  PATTERN,
	"-Pattern":  PATTERN,
}

var sequenceTypeIs = map[string]SequenceType{
	"INSDSeq:INSDInterval_from":       {1, ISSTART},
	"INSDSeq:INSDInterval_to":         {1, ISSTOP},
	"DocumentSummary:ChrStart":        {0, ISSTART},
	"DocumentSummary:ChrStop":         {0, ISSTOP},
	"DocumentSummary:Chr_start":       {1, ISSTART},
	"DocumentSummary:Chr_end":         {1, ISSTOP},
	"DocumentSummary:Chr_inner_start": {1, ISSTART},
	"DocumentSummary:Chr_inner_end":   {1, ISSTOP},
	"DocumentSummary:Chr_outer_start": {1, ISSTART},
	"DocumentSummary:Chr_outer_end":   {1, ISSTOP},
	"DocumentSummary:start":           {1, ISSTART},
	"DocumentSummary:stop":            {1, ISSTOP},
	"DocumentSummary:display_start":   {1, ISSTART},
	"DocumentSummary:display_stop":    {1, ISSTOP},
	"Entrezgene:Seq-interval_from":    {0, ISSTART},
	"Entrezgene:Seq-interval_to":      {0, ISSTOP},
	"GenomicInfoType:ChrStart":        {0, ISSTART},
	"GenomicInfoType:ChrStop":         {0, ISSTOP},
	"Rs:@aaPosition":                  {0, ISPOS},
	"Rs:@asnFrom":                     {0, ISSTART},
	"Rs:@asnTo":                       {0, ISSTOP},
	"Rs:@end":                         {0, ISSTOP},
	"Rs:@leftContigNeighborPos":       {0, ISSTART},
	"Rs:@physMapInt":                  {0, ISPOS},
	"Rs:@protLoc":                     {0, ISPOS},
	"Rs:@rightContigNeighborPos":      {0, ISSTOP},
	"Rs:@start":                       {0, ISSTART},
	"Rs:@structLoc":                   {0, ISPOS},
}

var isStopWord = map[string]bool{
	"a":             true,
	"about":         true,
	"again":         true,
	"all":           true,
	"almost":        true,
	"also":          true,
	"although":      true,
	"always":        true,
	"among":         true,
	"an":            true,
	"and":           true,
	"another":       true,
	"any":           true,
	"are":           true,
	"as":            true,
	"at":            true,
	"be":            true,
	"because":       true,
	"been":          true,
	"before":        true,
	"being":         true,
	"between":       true,
	"both":          true,
	"but":           true,
	"by":            true,
	"can":           true,
	"could":         true,
	"did":           true,
	"do":            true,
	"does":          true,
	"done":          true,
	"due":           true,
	"during":        true,
	"each":          true,
	"either":        true,
	"enough":        true,
	"especially":    true,
	"etc":           true,
	"for":           true,
	"found":         true,
	"from":          true,
	"further":       true,
	"had":           true,
	"has":           true,
	"have":          true,
	"having":        true,
	"here":          true,
	"how":           true,
	"however":       true,
	"i":             true,
	"if":            true,
	"in":            true,
	"into":          true,
	"is":            true,
	"it":            true,
	"its":           true,
	"itself":        true,
	"just":          true,
	"kg":            true,
	"km":            true,
	"made":          true,
	"mainly":        true,
	"make":          true,
	"may":           true,
	"mg":            true,
	"might":         true,
	"ml":            true,
	"mm":            true,
	"most":          true,
	"mostly":        true,
	"must":          true,
	"nearly":        true,
	"neither":       true,
	"no":            true,
	"nor":           true,
	"obtained":      true,
	"of":            true,
	"often":         true,
	"on":            true,
	"our":           true,
	"overall":       true,
	"perhaps":       true,
	"pmid":          true,
	"quite":         true,
	"rather":        true,
	"really":        true,
	"regarding":     true,
	"seem":          true,
	"seen":          true,
	"several":       true,
	"should":        true,
	"show":          true,
	"showed":        true,
	"shown":         true,
	"shows":         true,
	"significantly": true,
	"since":         true,
	"so":            true,
	"some":          true,
	"such":          true,
	"than":          true,
	"that":          true,
	"the":           true,
	"their":         true,
	"theirs":        true,
	"them":          true,
	"then":          true,
	"there":         true,
	"therefore":     true,
	"these":         true,
	"they":          true,
	"this":          true,
	"those":         true,
	"through":       true,
	"thus":          true,
	"to":            true,
	"upon":          true,
	"use":           true,
	"used":          true,
	"using":         true,
	"various":       true,
	"very":          true,
	"was":           true,
	"we":            true,
	"were":          true,
	"what":          true,
	"when":          true,
	"which":         true,
	"while":         true,
	"with":          true,
	"within":        true,
	"without":       true,
	"would":         true,
}

// DATA OBJECTS

type Tables struct {
	InBlank   [256]bool
	AltBlank  [256]bool
	InFirst   [256]bool
	InElement [256]bool
	ChanDepth int
	FarmSize  int
}

type Node struct {
	Name       string
	Parent     string
	Contents   string
	Attributes string
	Attribs    []string
	Children   *Node
	Next       *Node
}

type Step struct {
	Type   OpType
	Value  string
	Parent string
	Match  string
	Attrib string
	Wild   bool
}

type Operation struct {
	Type   OpType
	Value  string
	Stages []*Step
}

type Block struct {
	Visit      string
	Parent     string
	Match      string
	Working    []string
	Parsed     []string
	Position   string
	Conditions []*Operation
	Commands   []*Operation
	Failure    []*Operation
	Subtasks   []*Block
}

// UTILITIES

func IsNotJustWhitespace(str string) bool {

	for _, ch := range str {
		if ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r' && ch != '\f' {
			return true
		}
	}

	return false
}

func HasAmpOrNotASCII(str string) bool {

	for _, ch := range str {
		if ch == '&' || ch > 127 {
			return true
		}
	}

	return false
}

func IsAllCapsOrDigits(str string) bool {

	for _, rune := range str {
		if !unicode.IsUpper(rune) && !unicode.IsDigit(rune) {
			return false
		}
	}

	return true
}

func CompressRunsOfSpaces(str string) string {

	whiteSpace := false
	var buffer bytes.Buffer

	for _, rune := range str {
		if unicode.IsSpace(rune) {
			if !whiteSpace {
				buffer.WriteRune(' ')
			}
			whiteSpace = true
		} else {
			buffer.WriteRune(rune)
			whiteSpace = false
		}
	}

	return buffer.String()
}

func HasFlankingSpace(str string) bool {

	if str == "" {
		return false
	}

	ch := str[0]
	if ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || ch == '\f' {
		return true
	}

	strlen := len(str)
	ch = str[strlen-1]
	if ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || ch == '\f' {
		return true
	}

	return false
}

func HasBadSpace(str string) bool {

	for _, rune := range str {
		if unicode.IsSpace(rune) && rune != ' ' {
			return true
		}
	}

	return false
}

func CleanupBadSpaces(str string) string {

	var buffer bytes.Buffer

	for _, rune := range str {
		if unicode.IsSpace(rune) {
			buffer.WriteRune(' ')
		} else {
			buffer.WriteRune(rune)
		}
	}

	return buffer.String()
}

func SplitInTwoAt(str, chr string, side SideType) (string, string) {

	slash := strings.SplitN(str, chr, 2)
	if len(slash) > 1 {
		return slash[0], slash[1]
	}

	if side == LEFT {
		return str, ""
	}

	return "", str
}

func ConvertSlash(str string) string {

	if str == "" {
		return str
	}

	length := len(str)
	res := make([]byte, length+1, length+1)

	isSlash := false
	idx := 0
	for _, rune := range str {
		if isSlash {
			switch rune {
			case 'n':
				// line feed
				res[idx] = '\n'
			case 'r':
				// carriage return
				res[idx] = '\r'
			case 't':
				// horizontal tab
				res[idx] = '\t'
			case 'f':
				// form feed
				res[idx] = '\f'
			case 'a':
				// audible bell from terminal (undocumented)
				res[idx] = '\x07'
			default:
				res[idx] = byte(rune)
			}
			idx++
			isSlash = false
		} else if rune == '\\' {
			isSlash = true
		} else {
			res[idx] = byte(rune)
			idx++
		}
	}

	res = res[0:idx]

	return string(res)
}

func ParseFlag(str string) OpType {

	op, ok := opTypeIs[str]
	if ok {
		return op
	}

	if len(str) > 1 && str[0] == '-' && IsAllCapsOrDigits(str[1:]) {
		return VARIABLE
	}

	if len(str) > 0 && str[0] == '-' {
		return UNRECOGNIZED
	}

	return UNSET
}

// CREATE COMMON DRIVER TABLES

// InitTables creates lookup tables to simplify the tokenizer
func InitTables() *Tables {

	tbls := &Tables{}

	for i := range tbls.InBlank {
		tbls.InBlank[i] = false
	}
	tbls.InBlank[' '] = true
	tbls.InBlank['\t'] = true
	tbls.InBlank['\n'] = true
	tbls.InBlank['\r'] = true
	tbls.InBlank['\f'] = true

	// alternative version of InBlank allows newlines to be counted
	for i := range tbls.AltBlank {
		tbls.AltBlank[i] = false
	}
	tbls.AltBlank[' '] = true
	tbls.AltBlank['\t'] = true
	tbls.AltBlank['\r'] = true
	tbls.AltBlank['\f'] = true

	// first character of element cannot be a digit, dash, or period
	for i := range tbls.InFirst {
		tbls.InFirst[i] = false
	}
	for ch := 'A'; ch <= 'Z'; ch++ {
		tbls.InFirst[ch] = true
	}
	for ch := 'a'; ch <= 'z'; ch++ {
		tbls.InFirst[ch] = true
	}
	tbls.InFirst['_'] = true

	// remaining characters also includes colon for namespace
	for i := range tbls.InElement {
		tbls.InElement[i] = false
	}
	for ch := 'A'; ch <= 'Z'; ch++ {
		tbls.InElement[ch] = true
	}
	for ch := 'a'; ch <= 'z'; ch++ {
		tbls.InElement[ch] = true
	}
	for ch := '0'; ch <= '9'; ch++ {
		tbls.InElement[ch] = true
	}
	tbls.InElement['_'] = true
	tbls.InElement['-'] = true
	tbls.InElement['.'] = true
	tbls.InElement[':'] = true

	return tbls
}

// examine structure of parsed arguments (undocumented)
func DebugBlock(blk *Block, depth int) {

	doIndent := func(indt int) {
		for i := 1; i < indt; i++ {
			fmt.Fprintf(os.Stderr, "  ")
		}
	}

	doIndent(depth)

	if blk.Visit != "" {
		doIndent(depth + 1)
		fmt.Fprintf(os.Stderr, "<Visit> %s </Visit>\n", blk.Visit)
	}
	if len(blk.Parsed) > 0 {
		doIndent(depth + 1)
		fmt.Fprintf(os.Stderr, "<Parsed>")
		for _, str := range blk.Parsed {
			fmt.Fprintf(os.Stderr, " %s", str)
		}
		fmt.Fprintf(os.Stderr, " </Parsed>\n")
	}

	if len(blk.Subtasks) > 0 {
		for _, sub := range blk.Subtasks {
			DebugBlock(sub, depth+1)
		}
	}
}

// PARSE COMMAND-LINE ARGUMENTS

// ParseArguments parses nested exploration instruction from command-line arguments
func ParseArguments(args []string, pttrn string) *Block {

	// different names of exploration control arguments allow multiple levels of nested "for" loops in linear command line
	// (capitalized versions for backward-compatibility with original Perl implementation handling of recursive definitions)
	var (
		lcname = []string{
			"",
			"-unit",
			"-subset",
			"-section",
			"-block",
			"-branch",
			"-group",
			"-division",
			"-pattern",
		}

		ucname = []string{
			"",
			"-Unit",
			"-Subset",
			"-Section",
			"-Block",
			"-Branch",
			"-Group",
			"-Division",
			"-Pattern",
		}
	)

	/*
		xtract -pattern PubmedArticle -element MedlineCitation/PMID \
		  -block DateCreated -sep "-" -element Year,Month,Day \
		  -block Author -sep " " -tab "" -element "&COM" Initials,LastName -COM "(, )"

		<Pattern>
		  <Visit> PubmedArticle </Visit>
		  <Parsed> -element MedlineCitation/PMID </Parsed>
		  <Block>
			<Visit> DateCreated </Visit>
			<Parsed> -sep "-" -element Year,Month,Day </Parsed>
		  </Block>
		  <Block>
			<Visit> Author </Visit>
			<Parsed> -sep " " -tab "" -element &COM Initials,LastName -COM "(, )" </Parsed>
		  </Block>
		</Pattern>
	*/

	// parseCommands recursive definition
	var parseCommands func(parent *Block, startLevel LevelType)

	// parseCommands does initial parsing of exploration command structure
	parseCommands = func(parent *Block, startLevel LevelType) {

		// function to find next highest level exploration argument
		findNextLevel := func(args []string, level LevelType) (LevelType, string, string) {

			if len(args) > 1 {

				for {

					if level < UNIT {
						break
					}

					lctag := lcname[level]
					uctag := ucname[level]

					for _, txt := range args {
						if txt == lctag || txt == uctag {
							return level, lctag, uctag
						}
					}

					level--
				}
			}

			return 0, "", ""
		}

		arguments := parent.Working

		level, lctag, uctag := findNextLevel(arguments, startLevel)

		if level < UNIT {

			// break recursion
			return
		}

		// function to group arguments at a given exploration level
		subsetCommands := func(args []string) *Block {

			max := len(args)

			visit := ""

			// extract name of object to visit
			if max > 1 {
				visit = args[1]
				args = args[2:]
				max -= 2
			}

			partition := 0
			for cur, str := range args {

				// record point of next exploration command
				partition = cur + 1

				// skip if not a command
				if len(str) < 1 || str[0] != '-' {
					continue
				}

				if argTypeIs[str] == EXPLORATION {
					partition = cur
					break
				}
			}

			// parse parent/child construct
			// colon indicates a namespace prefix in any or all of the components
			prnt, match := SplitInTwoAt(visit, "/", RIGHT)

			// promote arguments parsed at this level
			return &Block{Visit: visit, Parent: prnt, Match: match, Parsed: args[0:partition], Working: args[partition:]}
		}

		cur := 0

		// search for positions of current exploration command

		for idx, txt := range arguments {
			if txt == lctag || txt == uctag {
				if idx == 0 {
					continue
				}

				blk := subsetCommands(arguments[cur:idx])
				parseCommands(blk, level-1)
				parent.Subtasks = append(parent.Subtasks, blk)

				cur = idx
			}
		}

		if cur < len(arguments) {
			blk := subsetCommands(arguments[cur:])
			parseCommands(blk, level-1)
			parent.Subtasks = append(parent.Subtasks, blk)
		}

		// clear execution arguments from parent after subsetting
		parent.Working = nil
	}

	parseConditionals := func(cmds *Block, arguments []string) []*Operation {

		max := len(arguments)
		if max < 1 {
			return nil
		}

		// check for missing condition command
		txt := arguments[0]
		if txt != "-if" && txt != "-unless" && txt != "-match" && txt != "-avoid" && txt != "-position" {
			fmt.Fprintf(os.Stderr, "\nERROR: Missing -if command before '%s'\n", txt)
			os.Exit(1)
		}
		if txt == "-position" && max > 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: Cannot combine -position with -if or -unless commands\n")
			os.Exit(1)
		}
		// check for missing argument after last condition
		txt = arguments[max-1]
		if len(txt) > 0 && txt[0] == '-' {
			fmt.Fprintf(os.Stderr, "\nERROR: Item missing after %s command\n", txt)
			os.Exit(1)
		}

		cond := make([]*Operation, 0, max)

		status := UNSET

		// function to parse conditional clause into execution step
		parseStep := func(op *Operation, elementColonValue bool) {

			if op == nil {
				return
			}

			str := op.Value

			status := ELEMENT

			// check for pound, percent, or caret character at beginning of name
			if len(str) > 1 {
				switch str[0] {
				case '&':
					if IsAllCapsOrDigits(str[1:]) {
						status = VARIABLE
						str = str[1:]
					} else if strings.Contains(str, ":") {
						fmt.Fprintf(os.Stderr, "\nERROR: Unsupported construct '%s', use -if &VARIABLE -equals VALUE instead\n", str)
						os.Exit(1)
					} else {
						fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized variable '%s'\n", str)
						os.Exit(1)
					}
				case '#':
					status = COUNT
					str = str[1:]
				case '%':
					status = LENGTH
					str = str[1:]
				case '^':
					status = DEPTH
					str = str[1:]
				default:
				}
			} else if str == "+" {
				status = INDEX
			}

			// parse parent/element@attribute construct
			// colon indicates a namespace prefix in any or all of the components
			prnt, match := SplitInTwoAt(str, "/", RIGHT)
			match, attrib := SplitInTwoAt(match, "@", LEFT)
			val := ""

			// leading colon indicates namespace prefix wildcard
			wildcard := false
			if strings.HasPrefix(prnt, ":") || strings.HasPrefix(match, ":") || strings.HasPrefix(attrib, ":") {
				wildcard = true
			}

			if elementColonValue {

				// allow parent/element@attribute:value construct for deprecated -match and -avoid, and for subsequent -and and -or commands
				match, val = SplitInTwoAt(str, ":", LEFT)
				prnt, match = SplitInTwoAt(match, "/", RIGHT)
				match, attrib = SplitInTwoAt(match, "@", LEFT)
			}

			tsk := &Step{Type: status, Value: str, Parent: prnt, Match: match, Attrib: attrib, Wild: wildcard}

			op.Stages = append(op.Stages, tsk)

			// transform old -match "element:value" to -match element -equals value
			if val != "" {
				tsk := &Step{Type: EQUALS, Value: val}
				op.Stages = append(op.Stages, tsk)
			}
		}

		idx := 0

		// conditionals should alternate between command and object/value
		expectDash := true
		last := ""

		var op *Operation

		// flag to allow element-colon-value for deprecated -match and -avoid commands, otherwise colon is for namespace prefixes
		elementColonValue := false

		// parse command strings into operation structure
		for idx < max {
			str := arguments[idx]
			idx++

			// conditionals should alternate between command and object/value
			if expectDash {
				if len(str) < 1 || str[0] != '-' {
					fmt.Fprintf(os.Stderr, "\nERROR: Unexpected '%s' argument after '%s'\n", str, last)
					os.Exit(1)
				}
				expectDash = false
			} else {
				if len(str) > 0 && str[0] == '-' {
					fmt.Fprintf(os.Stderr, "\nERROR: Unexpected '%s' command after '%s'\n", str, last)
					os.Exit(1)
				}
				expectDash = true
			}
			last = str

			switch status {
			case UNSET:
				status = ParseFlag(str)
			case POSITION:
				cmds.Position = str
				status = UNSET
			case MATCH, AVOID:
				elementColonValue = true
				fallthrough
			case IF, UNLESS, AND, OR:
				op = &Operation{Type: status, Value: str}
				cond = append(cond, op)
				parseStep(op, elementColonValue)
				status = UNSET
			case EQUALS, CONTAINS, STARTSWITH, ENDSWITH, ISNOT:
				if op != nil {
					if len(str) > 1 && str[0] == '\\' {
						// first character may be backslash protecting dash (undocumented)
						str = str[1:]
					}
					tsk := &Step{Type: status, Value: str}
					op.Stages = append(op.Stages, tsk)
					op = nil
				} else {
					fmt.Fprintf(os.Stderr, "\nERROR: Unexpected adjacent string match constraints\n")
					os.Exit(1)
				}
				status = UNSET
			case GT, GE, LT, LE, EQ, NE:
				if op != nil {
					if len(str) > 1 && str[0] == '\\' {
						// first character may be backslash protecting minus sign (undocumented)
						str = str[1:]
					}
					if len(str) < 1 {
						fmt.Fprintf(os.Stderr, "\nERROR: Empty numeric match constraints\n")
						os.Exit(1)
					}
					ch := str[0]
					if (ch >= '0' && ch <= '9') || ch == '-' || ch == '+' {
						// literal numeric constant
						tsk := &Step{Type: status, Value: str}
						op.Stages = append(op.Stages, tsk)
					} else if (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') {
						// numeric test allows element as second argument
						prnt, match := SplitInTwoAt(str, "/", RIGHT)
						match, attrib := SplitInTwoAt(match, "@", LEFT)
						wildcard := false
						if strings.HasPrefix(prnt, ":") || strings.HasPrefix(match, ":") || strings.HasPrefix(attrib, ":") {
							wildcard = true
						}
						tsk := &Step{Type: status, Value: str, Parent: prnt, Match: match, Attrib: attrib, Wild: wildcard}
						op.Stages = append(op.Stages, tsk)
					} else {
						fmt.Fprintf(os.Stderr, "\nERROR: Unexpected numeric match constraints\n")
						os.Exit(1)
					}
					op = nil
				} else {
					fmt.Fprintf(os.Stderr, "\nERROR: Unexpected adjacent numeric match constraints\n")
					os.Exit(1)
				}
				status = UNSET
			case UNRECOGNIZED:
				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized argument '%s'\n", str)
				os.Exit(1)
			default:
				fmt.Fprintf(os.Stderr, "\nERROR: Unexpected argument '%s'\n", str)
				os.Exit(1)
			}
		}

		return cond
	}

	parseExtractions := func(cmds *Block, arguments []string) []*Operation {

		max := len(arguments)
		if max < 1 {
			return nil
		}

		// check for missing -element (or -first, etc.) command
		txt := arguments[0]
		if len(txt) < 1 || txt[0] != '-' {
			fmt.Fprintf(os.Stderr, "\nERROR: Missing -element command before '%s'\n", txt)
			os.Exit(1)
		}
		// check for missing argument after last -element (or -first, etc.) command
		txt = arguments[max-1]
		if len(txt) > 0 && txt[0] == '-' {
			if txt == "-rst" {
				fmt.Fprintf(os.Stderr, "\nERROR: Unexpected position for %s command\n", txt)
				os.Exit(1)
			} else if txt == "-clr" {
			} else if max < 2 || arguments[max-2] != "-lbl" {
				fmt.Fprintf(os.Stderr, "\nERROR: Item missing after %s command\n", txt)
				os.Exit(1)
			}
		}

		comm := make([]*Operation, 0, max)

		status := UNSET

		// function to parse next argument
		nextStatus := func(str string) OpType {

			status = ParseFlag(str)

			switch status {
			case VARIABLE:
				op := &Operation{Type: status, Value: str[1:]}
				comm = append(comm, op)
				status = VALUE
			case CLR, RST:
				op := &Operation{Type: status, Value: ""}
				comm = append(comm, op)
				status = UNSET
			case ELEMENT, FIRST, LAST, ENCODE, UPPER, LOWER, TERMS, WORDS, PAIRS, PHRASE:
			case NUM, LEN, SUM, MIN, MAX, INC, DEC, SUB, AVG, DEV, ZEROBASED, ONEBASED, UCSC:
			case TAB, RET, PFX, SFX, SEP, LBL, PFC, DEF:
			case UNSET:
				fmt.Fprintf(os.Stderr, "\nERROR: No -element before '%s'\n", str)
				os.Exit(1)
			case UNRECOGNIZED:
				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized argument '%s'\n", str)
				os.Exit(1)
			default:
				fmt.Fprintf(os.Stderr, "\nERROR: Misplaced %s command\n", str)
				os.Exit(1)
			}

			return status
		}

		// function to parse extraction clause into individual steps
		parseSteps := func(op *Operation, pttrn string) {

			if op == nil {
				return
			}

			stat := op.Type
			str := op.Value

			// element names combined with commas are treated as a prefix-separator-suffix group
			comma := strings.Split(str, ",")

			for _, item := range comma {
				status := stat

				// check for special character at beginning of name
				if len(item) > 1 {
					switch item[0] {
					case '&':
						if IsAllCapsOrDigits(item[1:]) {
							status = VARIABLE
							item = item[1:]
						} else {
							fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized variable '%s'\n", item)
							os.Exit(1)
						}
					case '#':
						status = COUNT
						item = item[1:]
					case '%':
						status = LENGTH
						item = item[1:]
					case '^':
						status = DEPTH
						item = item[1:]
					case '*':
						for _, ch := range item {
							if ch != '*' {
								break
							}
						}
						status = STAR
					default:
					}
				} else if item == "*" {
					status = STAR
				} else if item == "+" {
					status = INDEX
				}

				// parse parent/element@attribute construct
				// colon indicates a namespace prefix in any or all of the components
				prnt, match := SplitInTwoAt(item, "/", RIGHT)
				match, attrib := SplitInTwoAt(match, "@", LEFT)

				// leading colon indicates namespace prefix wildcard
				wildcard := false
				if strings.HasPrefix(prnt, ":") || strings.HasPrefix(match, ":") || strings.HasPrefix(attrib, ":") {
					wildcard = true
				}

				// sequence coordinate adjustments
				switch status {
				case ZEROBASED, ONEBASED, UCSC:
					seq := pttrn + ":"
					if attrib != "" {
						seq += "@"
						seq += attrib
					} else if match != "" {
						seq += match
					}
					// confirm -0-based or -1-based arguments are known sequence position elements or attributes
					seqtype, ok := sequenceTypeIs[seq]
					if !ok {
						fmt.Fprintf(os.Stderr, "\nERROR: Element '%s' is not suitable for sequence coordinate conversion\n", item)
						os.Exit(1)
					}
					switch status {
					case ZEROBASED:
						status = ELEMENT
						// if 1-based coordinates, decrement to get 0-based value
						if seqtype.Based == 1 {
							status = DEC
						}
					case ONEBASED:
						status = ELEMENT
						// if 0-based coordinates, increment to get 1-based value
						if seqtype.Based == 0 {
							status = INC
						}
					case UCSC:
						status = ELEMENT
						// half-open intervals, start is 0-based, stop is 1-based
						if seqtype.Based == 0 && seqtype.Which == ISSTOP {
							status = INC
						} else if seqtype.Based == 1 && seqtype.Which == ISSTART {
							status = DEC
						}
					default:
						status = ELEMENT
					}
				default:
				}

				tsk := &Step{Type: status, Value: item, Parent: prnt, Match: match, Attrib: attrib, Wild: wildcard}

				op.Stages = append(op.Stages, tsk)
			}
		}

		idx := 0

		// parse command strings into operation structure
		for idx < max {
			str := arguments[idx]
			idx++

			if argTypeIs[str] == CONDITIONAL {
				fmt.Fprintf(os.Stderr, "\nERROR: Misplaced %s command\n", str)
				os.Exit(1)
			}

			switch status {
			case UNSET:
				status = nextStatus(str)
			case ELEMENT, FIRST, LAST, ENCODE, UPPER, LOWER, TERMS, WORDS, PAIRS, PHRASE, NUM, LEN, SUM, MIN, MAX, INC, DEC, SUB, AVG, DEV, ZEROBASED, ONEBASED, UCSC:
				for !strings.HasPrefix(str, "-") {
					// create one operation per argument, even if under a single -element statement
					op := &Operation{Type: status, Value: str}
					comm = append(comm, op)
					parseSteps(op, pttrn)
					if idx >= max {
						break
					}
					str = arguments[idx]
					idx++
				}
				status = UNSET
				if idx < max {
					status = nextStatus(str)
				}
			case TAB, RET, PFX, SFX, SEP, LBL, PFC, DEF:
				op := &Operation{Type: status, Value: ConvertSlash(str)}
				comm = append(comm, op)
				status = UNSET
			case VARIABLE:
				op := &Operation{Type: status, Value: str[1:]}
				comm = append(comm, op)
				status = VALUE
			case VALUE:
				op := &Operation{Type: status, Value: str}
				comm = append(comm, op)
				parseSteps(op, pttrn)
				status = UNSET
			case UNRECOGNIZED:
				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized argument '%s'\n", str)
				os.Exit(1)
			default:
			}
		}

		return comm
	}

	// parseOperations recursive definition
	var parseOperations func(parent *Block)

	// parseOperations converts parsed arguments to operations lists
	parseOperations = func(parent *Block) {

		args := parent.Parsed

		partition := 0
		for cur, str := range args {

			// record junction between conditional and extraction commands
			partition = cur + 1

			// skip if not a command
			if len(str) < 1 || str[0] != '-' {
				continue
			}

			if argTypeIs[str] != CONDITIONAL {
				partition = cur
				break
			}
		}

		// split arguments into conditional tests and extraction or customization commands
		conditionals := args[0:partition]
		args = args[partition:]

		partition = 0
		foundElse := false
		for cur, str := range args {

			// record junction at -else command
			partition = cur + 1

			// skip if not a command
			if len(str) < 1 || str[0] != '-' {
				continue
			}

			if str == "-else" {
				partition = cur
				foundElse = true
				break
			}
		}

		extractions := args[0:partition]
		alternative := args[partition:]

		if len(alternative) > 0 && alternative[0] == "-else" {
			alternative = alternative[1:]
		}

		// validate argument structure and convert to operations lists
		parent.Conditions = parseConditionals(parent, conditionals)
		parent.Commands = parseExtractions(parent, extractions)
		parent.Failure = parseExtractions(parent, alternative)

		// reality checks on placement of -else command
		if foundElse {
			if len(conditionals) < 1 {
				fmt.Fprintf(os.Stderr, "\nERROR: Misplaced -else command\n")
				os.Exit(1)
			}
			if len(alternative) < 1 {
				fmt.Fprintf(os.Stderr, "\nERROR: Misplaced -else command\n")
				os.Exit(1)
			}
			if len(parent.Subtasks) > 0 {
				fmt.Fprintf(os.Stderr, "\nERROR: Misplaced -else command\n")
				os.Exit(1)
			}
		}

		for _, sub := range parent.Subtasks {
			parseOperations(sub)
		}
	}

	// ParseArguments

	head := &Block{}

	for _, txt := range args {
		head.Working = append(head.Working, txt)
	}

	// initial parsing of exploration command structure
	parseCommands(head, PATTERN)

	if len(head.Subtasks) != 1 {
		return nil
	}

	// skip past empty placeholder
	head = head.Subtasks[0]

	// convert command strings to array of operations for faster processing
	parseOperations(head)

	// check for no -element or multiple -pattern commands
	noElement := true
	numPatterns := 0
	for _, txt := range args {
		if argTypeIs[txt] == EXTRACTION {
			noElement = false
		}
		if txt == "-pattern" || txt == "-Pattern" {
			numPatterns++
		}
	}

	if numPatterns < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: No -pattern in command-line arguments\n")
		os.Exit(1)
	}

	if numPatterns > 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Only one -pattern command is permitted\n")
		os.Exit(1)
	}

	if noElement {
		fmt.Fprintf(os.Stderr, "\nERROR: No -element statement in argument list\n")
		os.Exit(1)
	}

	return head
}

// READ XML INPUT FILE INTO SET OF BLOCKS

type XMLReader struct {
	Reader     io.Reader
	Buffer     []byte
	Remainder  string
	Closed     bool
	Docompress bool
	Docleanup  bool
}

func NewXMLReader(in io.Reader, doCompress, doCleanup bool) *XMLReader {

	if in == nil {
		return nil
	}

	rdr := &XMLReader{Reader: in, Docompress: doCompress, Docleanup: doCleanup}

	// 65536 appears to be the maximum number of characters presented to io.Reader when input is piped from stdin
	// increasing size of buffer when input is from a file does not improve program performance
	// additional 16384 bytes are reserved for copying previous remainder to start of buffer before next read
	const XMLBUFSIZE = 65536 + 16384

	rdr.Buffer = make([]byte, XMLBUFSIZE)

	return rdr
}

// NextBlock reads buffer, concatenates if necessary to place long element content into a single string
// all result strings end in > character that is used as a sentinel in subsequent code
func (rdr *XMLReader) NextBlock() string {

	if rdr == nil {
		return ""
	}

	// read one buffer, trim at last > and retain remainder for next call, signal if no > character
	nextBuffer := func() (string, bool, bool) {

		if rdr.Closed {
			return "", false, true
		}

		// prepend previous remainder to beginning of buffer
		m := copy(rdr.Buffer, rdr.Remainder)
		rdr.Remainder = ""
		if m > 16384 {
			// previous remainder is larger than reserved section, write and signal need to continue reading
			return string(rdr.Buffer[:]), true, false
		}

		// read next block, append behind copied remainder from previous read
		n, err := rdr.Reader.Read(rdr.Buffer[m:])
		// with data piped through stdin, read function may not always return the same number of bytes each time
		if err != nil {
			// end of file
			rdr.Closed = true
			// do not send final remainder (not terminated by right angle bracket that is used as a sentinel)
			return "", false, true
		}

		// slice of actual characters read
		bufr := rdr.Buffer[:n+m]

		// look for last > character
		pos := -1
		for pos = len(bufr) - 1; pos >= 0; pos-- {
			if bufr[pos] == '>' {
				break
			}
		}

		// trim back to last > character, save remainder for next buffer
		if pos > -1 {
			pos++
			rdr.Remainder = string(bufr[pos:])
			return string(bufr[:pos]), false, false
		}

		// no > found, signal need to continue reading long content
		return string(bufr[:]), true, false
	}

	// read next buffer
	line, cont, closed := nextBuffer()

	if closed {
		// no sentinel in remainder at end of file
		return ""
	}

	// if buffer does not end with > character
	if cont {
		var buff bytes.Buffer

		// keep reading long content blocks
		for {
			if line != "" {
				buff.WriteString(line)
			}
			if !cont {
				// last buffer ended with sentinel
				break
			}
			line, cont, closed = nextBuffer()
			if closed {
				// no sentinel in multi-block buffer at end of file
				return ""
			}
		}

		// concatenate blocks
		line = buff.String()
	}

	// trimming spaces here would throw off line tracking

	// optionally compress/cleanup tags/attributes and contents
	if rdr.Docompress {
		line = CompressRunsOfSpaces(line)
	}
	if rdr.Docleanup {
		if HasBadSpace(line) {
			line = CleanupBadSpaces(line)
		}
	}

	return line
}

// PARSE XML BLOCK STREAM INTO STRINGS FROM <PATTERN> TO </PATTERN>

// PartitionPattern splits XML input by pattern and sends individual records to a callback
func PartitionPattern(pat, star string, rdr *XMLReader, proc func(int, string)) {

	if pat == "" || rdr == nil || proc == nil {
		return
	}

	type Scanner struct {
		Pattern   string
		PatLength int
		CharSkip  [256]int
	}

	// function to initialize <pattern> to </pattern> scanner
	newScanner := func(pattern string) *Scanner {

		if pattern == "" {
			return nil
		}

		scr := &Scanner{Pattern: pattern}

		patlen := len(pattern)
		scr.PatLength = patlen

		// position of last character in pattern
		last := patlen - 1

		// initialize bad character displacement table
		for i := range scr.CharSkip {
			scr.CharSkip[i] = patlen
		}
		for i := 0; i < last; i++ {
			ch := pattern[i]
			scr.CharSkip[ch] = last - i
		}

		return scr
	}

	// function check surroundings of match candidate
	isAnElement := func(text string, lf, rt, mx int) bool {

		if (lf >= 0 && text[lf] == '<') || (lf > 0 && text[lf] == '/' && text[lf-1] == '<') {
			if (rt < mx && (text[rt] == '>' || text[rt] == ' ')) || (rt+1 < mx && text[rt] == '/' && text[rt+1] == '>') {
				return true
			}
		}

		return false
	}

	// modified Boyer-Moore-Horspool search function
	findNextMatch := func(scr *Scanner, text string, offset int) (int, int, int) {

		if scr == nil || text == "" {
			return -1, -1, -1
		}

		// copy values into local variables for speed
		txtlen := len(text)
		pattern := scr.Pattern[:]
		patlen := scr.PatLength
		max := txtlen - patlen
		last := patlen - 1
		skip := scr.CharSkip[:]

		i := offset

		for i <= max {
			j := last
			k := i + last
			for j >= 0 && text[k] == pattern[j] {
				j--
				k--
			}
			// require match candidate to be element name, i.e., <pattern ... >, </pattern ... >, or <pattern ... />
			if j < 0 && isAnElement(text, i-1, i+patlen, txtlen) {
				// find positions of flanking brackets
				lf := i - 1
				for lf > 0 && text[lf] != '<' {
					lf--
				}
				rt := i + patlen
				for rt < txtlen && text[rt] != '>' {
					rt++
				}
				return i + 1, lf, rt + 1
			}
			// find character in text above last character in pattern
			ch := text[i+last]
			// displacement table can shift pattern by one or more positions
			i += skip[ch]
		}

		return -1, -1, -1
	}

	type PatternType int

	const (
		NOPATTERN PatternType = iota
		STARTPATTERN
		SELFPATTERN
		STOPPATTERN
	)

	// function to find next element with pattern name
	nextPattern := func(scr *Scanner, text string, pos int) (PatternType, int, int) {

		if scr == nil || text == "" {
			return NOPATTERN, 0, 0
		}

		prev := pos

		for {
			next, start, stop := findNextMatch(scr, text, prev)
			if next < 0 {
				return NOPATTERN, 0, 0
			}

			prev = next + 1

			if text[start+1] == '/' {
				return STOPPATTERN, stop, prev
			} else if text[stop-2] == '/' {
				return SELFPATTERN, start, prev
			} else {
				return STARTPATTERN, start, prev
			}
		}
	}

	// -pattern Object construct

	doNormal := func() {

		// current depth of -pattern objects
		level := 0

		begin := 0
		inPattern := false

		line := ""
		var accumulator bytes.Buffer

		match := NOPATTERN
		pos := 0
		next := 0

		rec := 0

		scr := newScanner(pat)
		if scr == nil {
			return
		}

		for {

			begin = 0
			next = 0

			line = rdr.NextBlock()
			if line == "" {
				return
			}

			for {
				match, pos, next = nextPattern(scr, line, next)
				if match == STARTPATTERN {
					if level == 0 {
						inPattern = true
						begin = pos
					}
					level++
				} else if match == STOPPATTERN {
					level--
					if level == 0 {
						inPattern = false
						accumulator.WriteString(line[begin:pos])
						// read and process one -pattern object at a time
						str := accumulator.String()
						if str != "" {
							rec++
							proc(rec, str[:])
						}
						// reset accumulator
						accumulator.Reset()
					}
				} else if match == SELFPATTERN {
				} else {
					if inPattern {
						accumulator.WriteString(line[begin:])
					}
					break
				}
			}
		}
	}

	// -pattern Parent/* construct now works with catenated files, but not if components
	// are recursive or self-closing objects, process those through xtract -format first

	doStar := func() {

		// current depth of -pattern objects
		level := 0

		begin := 0
		inPattern := false

		line := ""
		var accumulator bytes.Buffer

		match := NOPATTERN
		pos := 0
		next := 0

		rec := 0

		scr := newScanner(pat)
		if scr == nil {
			return
		}

		last := pat

		// read to first <pattern> element
		for {

			next = 0

			line = rdr.NextBlock()
			if line == "" {
				break
			}

			match, pos, next = nextPattern(scr, line, next)
			if match == STARTPATTERN {
				break
			}
		}

		if match != STARTPATTERN {
			return
		}

		// function to find next element in XML
		nextElement := func(text string, pos int) string {

			txtlen := len(text)

			tag := ""
			for i := pos; i < txtlen; i++ {
				if text[i] == '<' {
					tag = text[i+1:]
					break
				}
			}
			if tag == "" {
				return ""
			}
			if tag[0] == '/' {
				if strings.HasPrefix(tag[1:], pat) {
					//should be </pattern> at end, want to continue if catenated files
					return "/"
				}
				return ""
			}
			for i, ch := range tag {
				if ch == '>' || ch == ' ' || ch == '/' {
					return tag[0:i]
				}
			}

			return ""
		}

		// read and process heterogeneous objects immediately below <pattern> parent
		for {
			tag := nextElement(line, next)
			if tag == "" {

				begin = 0
				next = 0

				line = rdr.NextBlock()
				if line == "" {
					break
				}

				tag = nextElement(line, next)
			}
			if tag == "" {
				return
			}

			// check for catenated parent set files
			if tag[0] == '/' {
				scr = newScanner(pat)
				if scr == nil {
					return
				}
				last = pat
				// confirm end </pattern> just found
				match, pos, next = nextPattern(scr, line, next)
				if match != STOPPATTERN {
					return
				}
				// now look for a new start <pattern> tag
				for {
					match, pos, next = nextPattern(scr, line, next)
					if match == STARTPATTERN {
						break
					}
					next = 0
					line = rdr.NextBlock()
					if line == "" {
						break
					}
				}
				if match != STARTPATTERN {
					return
				}
				// continue with processing loop
				continue
			}

			if tag != last {
				scr = newScanner(tag)
				if scr == nil {
					return
				}
				last = tag
			}

			for {
				match, pos, next = nextPattern(scr, line, next)
				if match == STARTPATTERN {
					if level == 0 {
						inPattern = true
						begin = pos
					}
					level++
				} else if match == STOPPATTERN {
					level--
					if level == 0 {
						inPattern = false
						accumulator.WriteString(line[begin:pos])
						// read and process one -pattern/* object at a time
						str := accumulator.String()
						if str != "" {
							rec++
							proc(rec, str[:])
						}
						// reset accumulator
						accumulator.Reset()
						break
					}
				} else {
					if inPattern {
						accumulator.WriteString(line[begin:])
					}

					begin = 0
					next = 0

					line = rdr.NextBlock()
					if line == "" {
						break
					}
				}
			}
		}
	}

	// call appropriate handler
	if star == "" {
		doNormal()
	} else if star == "*" {
		doStar()
	}
}

// XML VALIDATION AND FORMATTING FUNCTIONS

// ProcessXMLStream tokenizes and runs designated operations on an entire XML file
func ProcessXMLStream(in *XMLReader, tbls *Tables, args []string, action SpecialType) (int, int) {

	if in == nil || tbls == nil {
		return 0, 0
	}

	blockCount := 0

	// token parser variables
	Text := ""
	Txtlen := 0
	Idx := 0
	Line := 1

	// variables to track comments or CDATA sections that span reader blocks
	Which := NOTAG
	SkipTo := ""

	nextToken := func(idx int) (TagType, string, string, int, int) {

		if Text == "" {
			// if buffer is empty, read next block
			Text = in.NextBlock()
			Txtlen = len(Text)
			Idx = 0
			idx = 0
			blockCount++
		}

		if Text == "" {
			return ISCLOSED, "", "", Line, 0
		}

		// lookup table array pointers
		inBlank := &tbls.AltBlank
		inFirst := &tbls.InFirst
		inElement := &tbls.InElement

		text := Text[:]
		txtlen := Txtlen
		line := Line

		if Which != NOTAG && SkipTo != "" {
			which := Which
			// previous block ended inside CDATA object or comment
			start := idx
			found := strings.Index(text[:], SkipTo)
			if found < 0 {
				// no stop signal found in next block
				// count lines
				for i := 0; i < txtlen; i++ {
					if text[i] == '\n' {
						line++
					}
				}
				Line = line
				str := text[:]
				if HasFlankingSpace(str) {
					str = strings.TrimSpace(str)
				}
				// signal end of current block
				Text = ""
				// leave Which and SkipTo values unchanged as another continuation signal
				// send CDATA or comment contents
				return which, str[:], "", Line, 0
			}
			// otherwise adjust position past end of skipTo string and return to normal processing
			idx += found
			// count lines
			for i := 0; i < idx; i++ {
				if text[i] == '\n' {
					line++
				}
			}
			Line = line
			str := text[start:idx]
			if HasFlankingSpace(str) {
				str = strings.TrimSpace(str)
			}
			idx += len(SkipTo)
			// clear tracking variables
			Which = NOTAG
			SkipTo = ""
			// send CDATA or comment contents
			return which, str[:], "", Line, idx
		}

		// all blocks end with > character, acts as sentinel to check if past end of text
		if idx >= txtlen {
			// signal end of current block, will read next block on next call
			Text = ""
			Line = line
			return NOTAG, "", "", Line, 0
		}

		// skip past leading blanks
		ch := text[idx]
		for {
			for inBlank[ch] {
				idx++
				ch = text[idx]
			}
			if ch != '\n' {
				break
			}
			line++
			idx++
			ch = text[idx]
		}
		Line = line

		start := idx

		if ch == '<' {

			// at start of element
			idx++
			ch = text[idx]

			// check for legal first character of element
			if inFirst[ch] {

				// read element name
				start = idx
				idx++

				ch = text[idx]
				for inElement[ch] {
					idx++
					ch = text[idx]
				}

				str := text[start:idx]

				switch ch {
				case '>':
					// end of element
					idx++

					return STARTTAG, str[:], "", Line, idx
				case '/':
					// self-closing element without attributes
					idx++
					ch = text[idx]
					if ch != '>' {
						fmt.Fprintf(os.Stderr, "\nSelf-closing element missing right angle bracket\n")
					}
					idx++

					return SELFTAG, str[:], "", Line, idx
				case '\n':
					line++
					fallthrough
				case ' ', '\t', '\r', '\f':
					// attributes
					idx++
					start = idx
					ch = text[idx]
					for {
						for ch != '<' && ch != '>' && ch != '\n' {
							idx++
							ch = text[idx]
						}
						if ch != '\n' {
							break
						}
						line++
						idx++
						ch = text[idx]
					}
					Line = line
					if ch != '>' {
						fmt.Fprintf(os.Stderr, "\nAttributes not followed by right angle bracket\n")
					}
					if text[idx-1] == '/' {
						// self-closing
						atr := text[start : idx-1]
						idx++
						return SELFTAG, str[:], atr[:], Line, idx
					}
					atr := text[start:idx]
					idx++
					return STARTTAG, str[:], atr[:], Line, idx
				default:
					fmt.Fprintf(os.Stderr, "\nUnexpected punctuation in XML element\n")
					return STARTTAG, str[:], "", Line, idx
				}

			} else {

				// punctuation character immediately after first angle bracket
				switch ch {
				case '/':
					// at start of end tag
					idx++
					start = idx
					ch = text[idx]
					// expect legal first character of element
					if inFirst[ch] {
						idx++
						ch = text[idx]
						for inElement[ch] {
							idx++
							ch = text[idx]
						}
						str := text[start:idx]
						if ch != '>' {
							fmt.Fprintf(os.Stderr, "\nUnexpected characters after end element name\n")
						}
						idx++

						return STOPTAG, str[:], "", Line, idx
					}
					fmt.Fprintf(os.Stderr, "\nUnexpected punctuation in XML element\n")
				case '?':
					// skip ?xml and ?processing instructions
					idx++
					ch = text[idx]
					for ch != '>' {
						idx++
						ch = text[idx]
					}
					idx++
					return NOTAG, "", "", Line, idx
				case '!':
					// skip !DOCTYPE, !comment, and ![CDATA[
					idx++
					start = idx
					ch = text[idx]
					Which := NOTAG
					SkipTo := ""
					if ch == '[' && strings.HasPrefix(text[idx:], "[CDATA[") {
						Which = CDATATAG
						SkipTo = "]]>"
						start += 7
					} else if ch == '-' && strings.HasPrefix(text[idx:], "--") {
						Which = COMMENTTAG
						SkipTo = "-->"
						start += 2
					}
					if Which != NOTAG && SkipTo != "" {
						which := Which
						// CDATA or comment block may contain internal angle brackets
						found := strings.Index(text[idx:], SkipTo)
						if found < 0 {
							// string stops in middle of CDATA or comment
							// count lines
							for i := start; i < txtlen; i++ {
								if text[i] == '\n' {
									line++
								}
							}
							Line = line
							str := text[start:]
							if HasFlankingSpace(str) {
								str = strings.TrimSpace(str)
							}
							// signal end of current block
							Text = ""
							// leave Which and SkipTo values unchanged as another continuation signal
							// send CDATA or comment contents
							return which, str[:], "", Line, 0
						}
						// adjust position past end of CDATA or comment
						idx += found
						// count lines
						for i := start; i < idx; i++ {
							if text[i] == '\n' {
								line++
							}
						}
						Line = line
						str := text[start:idx]
						if HasFlankingSpace(str) {
							str = strings.TrimSpace(str)
						}
						idx += len(SkipTo)
						// clear tracking variables
						Which = NOTAG
						SkipTo = ""
						// send CDATA or comment contents
						return which, str[:], "", Line, idx
					}
					// otherwise just skip to next right angle bracket
					for ch != '>' {
						if ch == '\n' {
							line++
						}
						idx++
						ch = text[idx]
					}
					Line = line
					idx++
					return NOTAG, "", "", Line, idx
				default:
					fmt.Fprintf(os.Stderr, "\nUnexpected punctuation in XML element\n")
				}
			}

		} else if ch != '>' {

			// at start of contents
			start = idx

			// find end of contents
			for {
				for ch != '<' && ch != '>' && ch != '\n' {
					idx++
					ch = text[idx]
				}
				if ch != '\n' {
					break
				}
				line++
				idx++
				ch = text[idx]
			}
			Line = line

			// trim back past trailing blanks
			lst := idx - 1
			ch = text[lst]
			for inBlank[ch] && lst > start {
				lst--
				ch = text[lst]
			}

			str := text[start : lst+1]

			return CONTENTTAG, str[:], "", Line, idx
		}

		// signal end of current block, will read next block on next call
		Text = ""
		Line = line
		return NOTAG, "", "", Line, 0
	}

	// common output buffer
	var buffer bytes.Buffer
	count := 0

	// processOutline displays outline of XML structure
	processOutline := func() {

		indent := 0

		for {
			tag, name, _, _, idx := nextToken(Idx)
			Idx = idx

			switch tag {
			case STARTTAG:
				if name == "eSummaryResult" ||
					name == "eLinkResult" ||
					name == "eInfoResult" ||
					name == "PubmedArticleSet" ||
					name == "DocumentSummarySet" ||
					name == "INSDSet" ||
					name == "Entrezgene-Set" ||
					name == "TaxaSet" {
					break
				}
				for i := 0; i < indent; i++ {
					buffer.WriteString("  ")
				}
				buffer.WriteString(name)
				buffer.WriteString("\n")
				indent++
			case SELFTAG:
				for i := 0; i < indent; i++ {
					buffer.WriteString("  ")
				}
				buffer.WriteString(name)
				buffer.WriteString("\n")
			case STOPTAG:
				indent--
			case NOTAG:
			case ISCLOSED:
				txt := buffer.String()
				if txt != "" {
					// print final buffer
					fmt.Fprintf(os.Stdout, "%s", txt)
				}
				return
			default:
			}

			count++
			if count > 1000 {
				count = 0
				txt := buffer.String()
				if txt != "" {
					// print current buffered output
					fmt.Fprintf(os.Stdout, "%s", txt)
				}
				buffer.Reset()
			}
		}
	}

	// processSynopsis displays paths to XML elements
	processSynopsis := func() {

		// synopsisLevel recursive definition
		var synopsisLevel func(string) bool

		synopsisLevel = func(parent string) bool {

			for {
				tag, name, _, _, idx := nextToken(Idx)
				Idx = idx

				switch tag {
				case STARTTAG:
					if name == "eSummaryResult" ||
						name == "eLinkResult" ||
						name == "eInfoResult" ||
						name == "PubmedArticleSet" ||
						name == "DocumentSummarySet" ||
						name == "INSDSet" ||
						name == "Entrezgene-Set" ||
						name == "TaxaSet" {
						break
					}
					if parent != "" {
						buffer.WriteString(parent)
						buffer.WriteString("/")
					}
					buffer.WriteString(name)
					buffer.WriteString("\n")
					path := parent
					if path != "" {
						path += "/"
					}
					path += name
					if synopsisLevel(path) {
						return true
					}
				case SELFTAG:
					if parent != "" {
						buffer.WriteString(parent)
						buffer.WriteString("/")
					}
					buffer.WriteString(name)
					buffer.WriteString("\n")
				case STOPTAG:
					// break recursion
					return false
				case NOTAG:
				case ISCLOSED:
					txt := buffer.String()
					if txt != "" {
						// print final buffer
						fmt.Fprintf(os.Stdout, "%s", txt)
					}
					return true
				default:
				}

				count++
				if count > 1000 {
					count = 0
					txt := buffer.String()
					if txt != "" {
						// print current buffered output
						fmt.Fprintf(os.Stdout, "%s", txt)
					}
					buffer.Reset()
				}
			}
		}

		for {
			// may have concatenated XMLs, loop through all
			if synopsisLevel("") {
				return
			}
		}
	}

	// processVerify checks for well-formed XML
	processVerify := func() {

		type VerifyType int

		const (
			_ VerifyType = iota
			START
			STOP
			CHAR
			OTHER
		)

		// skip past command name
		args = args[1:]

		pttrn := ""

		if len(args) > 0 {
			pttrn = args[0]
			args = args[1:]
		}

		// if pattern supplied, report maximum nesting depth and record spanning the most blocks (undocumented)
		maxDepth := 0
		depthLine := 0
		maxBlocks := 0
		blockLine := 0
		startLine := 0

		// verifyLevel recursive definition
		var verifyLevel func(string, int)

		// verify integrity of XML object nesting (well-formed)
		verifyLevel = func(parent string, level int) {

			status := START
			for {
				// use alternative low-level tokenizer
				tag, name, _, line, idx := nextToken(Idx)
				Idx = idx

				if level > maxDepth {
					maxDepth = level
					depthLine = line
				}

				switch tag {
				case STARTTAG:
					if status == CHAR {
						fmt.Fprintf(os.Stdout, "<%s> not expected after contents, line %d\n", name, line)
					}
					if name == pttrn {
						blockCount = 1
						startLine = line
					}
					verifyLevel(name, level+1)
					// returns here after recursion
					status = STOP
				case SELFTAG:
					status = OTHER
				case STOPTAG:
					if name == pttrn {
						if blockCount > maxBlocks {
							maxBlocks = blockCount
							blockLine = startLine
						}
					}
					if parent != name && parent != "" {
						fmt.Fprintf(os.Stdout, "Expected </%s>, found </%s>, line %d\n", parent, name, line)
					}
					if level < 1 {
						fmt.Fprintf(os.Stdout, "Unexpected </%s> at end of XML, line %d\n", name, line)
					}
					// break recursion
					return
				case CONTENTTAG:
					if status != START {
						fmt.Fprintf(os.Stdout, "Contents not expected before </%s>, line %d - status %d\n", parent, line, status)
					}
					status = CHAR
				case CDATATAG, COMMENTTAG:
					status = OTHER
				case NOTAG:
				case ISCLOSED:
					if level > 0 {
						fmt.Fprintf(os.Stdout, "Unexpected end of data\n")
					}
					return
				default:
					status = OTHER
				}
			}
		}

		verifyLevel("", 0)

		if pttrn != "" {
			fmt.Fprintf(os.Stdout, "Maximum nesting (%d levels) at line %d\n", maxDepth, depthLine)
			fmt.Fprintf(os.Stdout, "Longest pattern (%d blocks) at line %d\n", maxBlocks, blockLine)
		}
	}

	// processFilter modifies XML content, comments, or CDATA
	processFilter := func() {

		// skip past command name
		args = args[1:]

		max := len(args)
		if max < 1 {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract -filter\n")
			os.Exit(1)
		}

		pttrn := args[0]

		args = args[1:]
		max--

		if max < 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: No object name supplied to xtract -filter\n")
			os.Exit(1)
		}

		type ActionType int

		const (
			NOACTION ActionType = iota
			DORETAIN
			DOREMOVE
			DOENCODE
			DODECODE
			DOSHRINK
		)

		action := args[0]

		what := NOACTION
		switch action {
		case "retain":
			what = DORETAIN
		case "remove":
			what = DOREMOVE
		case "encode":
			what = DOENCODE
		case "decode":
			what = DODECODE
		case "shrink":
			what = DOSHRINK
		default:
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized action '%s' supplied to xtract -filter\n", action)
			os.Exit(1)
		}

		trget := args[1]

		which := NOTAG
		switch trget {
		case "attribute", "attributes":
			which = ATTRIBTAG
		case "content", "contents":
			which = CONTENTTAG
		case "cdata", "CDATA":
			which = CDATATAG
		case "comment", "comments":
			which = COMMENTTAG
		case "object":
			// object normally retained
			which = OBJECTTAG
		default:
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized target '%s' supplied to xtract -filter\n", trget)
			os.Exit(1)
		}

		inPattern := false

		for {
			tag, name, attr, _, idx := nextToken(Idx)
			Idx = idx

			switch tag {
			case STARTTAG:
				if name == pttrn {
					inPattern = true
				}
				if inPattern && which == OBJECTTAG && what == DOREMOVE {
					continue
				}
				buffer.WriteString("<")
				buffer.WriteString(name)
				if attr != "" {
					if which != ATTRIBTAG || what != DOREMOVE {
						attr = strings.TrimSpace(attr)
						attr = CompressRunsOfSpaces(attr)
						buffer.WriteString(" ")
						buffer.WriteString(attr)
					}
				}
				buffer.WriteString(">\n")
			case SELFTAG:
				if inPattern && which == OBJECTTAG && what == DOREMOVE {
					continue
				}
				buffer.WriteString("<")
				buffer.WriteString(name)
				if attr != "" {
					if which != ATTRIBTAG || what != DOREMOVE {
						attr = strings.TrimSpace(attr)
						attr = CompressRunsOfSpaces(attr)
						buffer.WriteString(" ")
						buffer.WriteString(attr)
					}
				}
				buffer.WriteString("/>\n")
			case STOPTAG:
				if name == pttrn {
					inPattern = false
					if which == OBJECTTAG && what == DOREMOVE {
						continue
					}
				}
				if inPattern && which == OBJECTTAG && what == DOREMOVE {
					continue
				}
				buffer.WriteString("</")
				buffer.WriteString(name)
				buffer.WriteString(">\n")
			case CONTENTTAG:
				if inPattern && which == OBJECTTAG && what == DOREMOVE {
					continue
				}
				if inPattern && which == tag {
					switch what {
					case DORETAIN:
						// default behavior for content - can use -filter X retain content as a no-op
					case DOREMOVE:
						continue
					case DOENCODE:
						name = html.EscapeString(name)
					case DODECODE:
						name = html.UnescapeString(name)
					case DOSHRINK:
						name = CompressRunsOfSpaces(name)
					default:
						continue
					}
				}
				// content normally printed
				if HasFlankingSpace(name) {
					name = strings.TrimSpace(name)
				}
				buffer.WriteString(name)
				buffer.WriteString("\n")
			case CDATATAG, COMMENTTAG:
				if inPattern && which == OBJECTTAG && what == DOREMOVE {
					continue
				}
				if inPattern && which == tag {
					switch what {
					case DORETAIN:
						// cdata and comment require explicit retain command
					case DOREMOVE:
						continue
					case DOENCODE:
						name = html.EscapeString(name)
					case DODECODE:
						name = html.UnescapeString(name)
					case DOSHRINK:
						name = CompressRunsOfSpaces(name)
					default:
						continue
					}
					// cdata and comment normally removed
					if HasFlankingSpace(name) {
						name = strings.TrimSpace(name)
					}
					buffer.WriteString(name)
					buffer.WriteString("\n")
				}
			case NOTAG:
			case ISCLOSED:
				txt := buffer.String()
				if txt != "" {
					// print final buffer
					fmt.Fprintf(os.Stdout, "%s", txt)
				}
				return
			default:
			}

			count++
			if count > 1000 {
				count = 0
				txt := buffer.String()
				if txt != "" {
					// print current buffered output
					fmt.Fprintf(os.Stdout, "%s", txt)
				}
				buffer.Reset()
			}
		}
	}

	// processFormat reformats XML for ease of reading
	processFormat := func() {

		// skip past command name
		args = args[1:]

		compRecrd := false
		wrapAttrs := false
		ret := "\n"
		frst := true

		if len(args) > 0 {
			switch args[0] {
			case "compact", "compacted", "compress", "compressed", "terse", "*":
				// compress to one record per line
				compRecrd = true
				ret = ""
			case "expand", "expanded", "verbose", "@":
				// each attribute on its own line
				wrapAttrs = true
			case "indent", "indented", "normal":
			default:
				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized option after -format command\n")
				os.Exit(1)
			}
		}

		type FormatType int

		const (
			NOTSET FormatType = iota
			START
			STOP
			CHAR
			OTHER
		)

		// array to speed up indentation
		indentSpaces := []string{
			"",
			"  ",
			"    ",
			"      ",
			"        ",
			"          ",
			"            ",
			"              ",
			"                ",
			"                  ",
		}

		indent := 0

		// parent used to detect first start tag, will place in doctype line
		parent := ""

		status := NOTSET

		// delay printing right bracket of start tag to support self-closing tag style
		needsRightBracket := ""

		// delay printing start tag if no attributes, suppress empty start-end pair if followed by end
		justStartName := ""
		justStartIndent := 0

		// function to indent a specified number of spaces
		doIndent := func(indt int) {
			if compRecrd {
				return
			}
			i := indt
			for i > 9 {
				buffer.WriteString("                    ")
				i -= 10
			}
			if i < 0 {
				return
			}
			buffer.WriteString(indentSpaces[i])
		}

		// function to handle delayed start tag
		doDelayedName := func() {
			if needsRightBracket != "" {
				buffer.WriteString(">")
				needsRightBracket = ""
			}
			if justStartName != "" {
				doIndent(justStartIndent)
				buffer.WriteString("<")
				buffer.WriteString(justStartName)
				buffer.WriteString(">")
				justStartName = ""
			}
		}

		closingTag := ""

		// function to print attributes
		printAttributes := func(attr string) {

			attr = strings.TrimSpace(attr)
			attr = CompressRunsOfSpaces(attr)

			if wrapAttrs {

				start := 0
				idx := 0

				attlen := len(attr)

				for idx < attlen {
					ch := attr[idx]
					if ch == '=' {
						str := attr[start:idx]
						buffer.WriteString("\n")
						doIndent(indent)
						buffer.WriteString(" ")
						buffer.WriteString(str)
						// skip past equal sign and leading double quote
						idx += 2
						start = idx
					} else if ch == '"' {
						str := attr[start:idx]
						buffer.WriteString("=\"")
						buffer.WriteString(str)
						buffer.WriteString("\"")
						// skip past trailing double quote and (possible) space
						idx += 2
						start = idx
					} else {
						idx++
					}
				}

				buffer.WriteString("\n")
				doIndent(indent)

			} else {

				buffer.WriteString(" ")
				buffer.WriteString(attr)
			}
		}

		for {
			tag, name, attr, _, idx := nextToken(Idx)
			Idx = idx

			switch tag {
			case STARTTAG:
				doDelayedName()
				if status == START {
					buffer.WriteString(ret)
				}
				// remove internal copies of </parent><parent> tags
				if parent != "" && name == parent && indent == 1 {
					continue
				}

				// detect first start tag, print xml and doctype parent
				if indent == 0 && parent == "" {
					parent = name
					buffer.WriteString("<?xml version=\"1.0\"?>\n")
					buffer.WriteString("<!DOCTYPE ")
					buffer.WriteString(parent)
					buffer.WriteString(">\n")
					// now filtering internal </parent><parent> tags, so queue printing of closing tag
					closingTag = fmt.Sprintf("</%s>\n", parent)
					// already past </parent><parent> test, so opening tag will print normally
				}

				// check for attributes
				if attr != "" {
					doIndent(indent)

					buffer.WriteString("<")
					buffer.WriteString(name)

					printAttributes(attr)

					needsRightBracket = name

				} else {
					justStartName = name
					justStartIndent = indent
				}

				if compRecrd && frst && indent == 0 {
					frst = false
					doDelayedName()
					buffer.WriteString("\n")
				}

				indent++

				status = START
			case SELFTAG:
				doDelayedName()
				if status == START {
					buffer.WriteString(ret)
				}

				// suppress self-closing tag without attributes attributes
				if attr != "" {
					doIndent(indent)

					buffer.WriteString("<")
					buffer.WriteString(name)

					printAttributes(attr)

					buffer.WriteString("/>")
					buffer.WriteString(ret)
				}

				status = STOP
			case STOPTAG:
				// if end immediately follows start, turn into self-closing tag if there were attributes, otherwise suppress empty tag
				if needsRightBracket != "" {
					if status == START && name == needsRightBracket {
						// end immediately follows start, produce self-closing tag
						buffer.WriteString("/>")
						buffer.WriteString(ret)
						needsRightBracket = ""
						indent--
						status = STOP
						break
					}
					buffer.WriteString(">")
					needsRightBracket = ""
				}
				if justStartName != "" {
					if status == START && name == justStartName {
						// end immediately follows delayed start with no attributes, suppress
						justStartName = ""
						indent--
						status = STOP
						break
					}
					doIndent(justStartIndent)
					buffer.WriteString("<")
					buffer.WriteString(justStartName)
					buffer.WriteString(">")
					justStartName = ""
				}

				// remove internal copies of </parent><parent> tags
				if parent != "" && name == parent && indent == 1 {
					continue
				}
				indent--
				if status == CHAR {
					buffer.WriteString("</")
					buffer.WriteString(name)
					buffer.WriteString(">")
					buffer.WriteString(ret)
				} else if status == START {
					buffer.WriteString("</")
					buffer.WriteString(name)
					buffer.WriteString(">")
					buffer.WriteString(ret)
				} else {
					doIndent(indent)

					buffer.WriteString("</")
					buffer.WriteString(name)
					buffer.WriteString(">")
					buffer.WriteString(ret)
				}
				status = STOP
				if compRecrd && indent == 1 {
					buffer.WriteString("\n")
				}
			case CONTENTTAG:
				doDelayedName()
				if len(name) > 0 && IsNotJustWhitespace(name) {
					if HasFlankingSpace(name) {
						name = strings.TrimSpace(name)
					}
					buffer.WriteString(name)
					status = CHAR
				}
			case CDATATAG, COMMENTTAG:
				// ignore
			case NOTAG:
			case ISCLOSED:
				doDelayedName()
				if closingTag != "" {
					buffer.WriteString(closingTag)
				}
				txt := buffer.String()
				if txt != "" {
					// print final buffer
					fmt.Fprintf(os.Stdout, "%s", txt)
				}
				return
			default:
				doDelayedName()
				status = OTHER
			}

			count++
			if count > 1000 {
				count = 0
				txt := buffer.String()
				if txt != "" {
					// print current buffered output
					fmt.Fprintf(os.Stdout, "%s", txt)
				}
				buffer.Reset()
			}
		}
	}

	// return variables set by performance test functions
	recordCount := 0
	byteCount := 0

	// processChunk reads a set of XML blocks
	processChunk := func() {

		for {
			str := in.NextBlock()
			if str == "" {
				break
			}
			recordCount++
			byteCount += len(str)
		}
	}

	// processSplit partitions XML by pattern
	processSplit := func() {

		if len(args) > 1 {
			if args[1] == "-pattern" || args[1] == "-Pattern" {
				// skip past -split if followed by -pattern
				args = args[1:]
			}
		}
		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -split command\n")
			os.Exit(1)
		}
		pat := args[1]

		PartitionPattern(pat, "", in,
			func(rec int, str string) {
				recordCount++
				byteCount += len(str)
			})
	}

	// processDrain partitions XML by pattern and sends them down a channel
	processDrain := func() {

		if len(args) > 1 {
			if args[1] == "-pattern" || args[1] == "-Pattern" {
				// skip past -drain if followed by -pattern
				args = args[1:]
			}
		}
		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -drain command\n")
			os.Exit(1)
		}
		pat := args[1]

		chn := make(chan string, tbls.ChanDepth)
		if chn == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create output channel\n")
			os.Exit(1)
		}

		sendPatterns := func(pat string, out chan<- string) {
			defer close(out)
			PartitionPattern(pat, "", in,
				func(rec int, str string) {
					out <- str
				})
		}

		go sendPatterns(pat, chn)

		for str := range chn {
			recordCount++
			byteCount += len(str)
		}
	}

	// processToken tokenizes the XML block stream
	processToken := func() {

		for {
			tag, name, attr, _, idx := nextToken(Idx)
			Idx = idx

			if tag == ISCLOSED {
				break
			}
			recordCount++
			byteCount += len(name) + len(attr)
		}
	}

	// ProcessXMLStream

	// call specific function
	switch action {
	case DOFORMAT:
		processFormat()
	case DOOUTLINE:
		processOutline()
	case DOSYNOPSIS:
		processSynopsis()
	case DOVERIFY:
		processVerify()
	case DOFILTER:
		processFilter()
	case DOCHUNK:
		processChunk()
	case DOSPLIT:
		processSplit()
	case DODRAIN:
		processDrain()
	case DOTOKEN:
		processToken()
	default:
	}

	return recordCount, byteCount
}

// INSDSEQ EXTRACTION COMMAND GENERATOR

// e.g., xtract -insd complete mat_peptide "%peptide" product peptide

// ProcessINSD generates extraction commands for GenBank/RefSeq records in INSDSet format
func ProcessINSD(args []string, isPipe, addDash bool) []string {

	// legal GenBank / GenPept / RefSeq features

	features := []string{
		"-10_signal",
		"-35_signal",
		"3'clip",
		"3'UTR",
		"5'clip",
		"5'UTR",
		"allele",
		"assembly_gap",
		"attenuator",
		"Bond",
		"C_region",
		"CAAT_signal",
		"CDS",
		"centromere",
		"conflict",
		"D_segment",
		"D-loop",
		"enhancer",
		"exon",
		"gap",
		"GC_signal",
		"gene",
		"iDNA",
		"intron",
		"J_segment",
		"LTR",
		"mat_peptide",
		"misc_binding",
		"misc_difference",
		"misc_feature",
		"misc_recomb",
		"misc_RNA",
		"misc_signal",
		"misc_structure",
		"mobile_element",
		"modified_base",
		"mRNA",
		"mutation",
		"N_region",
		"ncRNA",
		"old_sequence",
		"operon",
		"oriT",
		"polyA_signal",
		"polyA_site",
		"precursor_RNA",
		"prim_transcript",
		"primer_bind",
		"promoter",
		"propeptide",
		"protein_bind",
		"Protein",
		"RBS",
		"Region",
		"regulatory",
		"rep_origin",
		"repeat_region",
		"repeat_unit",
		"rRNA",
		"S_region",
		"satellite",
		"scRNA",
		"sig_peptide",
		"Site",
		"snoRNA",
		"snRNA",
		"source",
		"stem_loop",
		"STS",
		"TATA_signal",
		"telomere",
		"terminator",
		"tmRNA",
		"transit_peptide",
		"tRNA",
		"unsure",
		"V_region",
		"V_segment",
		"variation",
	}

	// legal GenBank / GenPept / RefSeq qualifiers

	qualifiers := []string{
		"allele",
		"altitude",
		"anticodon",
		"artificial_location",
		"bio_material",
		"bond_type",
		"bound_moiety",
		"breed",
		"calculated_mol_wt",
		"cell_line",
		"cell_type",
		"chloroplast",
		"chromoplast",
		"chromosome",
		"citation",
		"clone_lib",
		"clone",
		"coded_by",
		"codon_start",
		"codon",
		"collected_by",
		"collection_date",
		"compare",
		"cons_splice",
		"country",
		"cultivar",
		"culture_collection",
		"cyanelle",
		"db_xref",
		"derived_from",
		"dev_stage",
		"direction",
		"EC_number",
		"ecotype",
		"encodes",
		"endogenous_virus",
		"environmental_sample",
		"estimated_length",
		"evidence",
		"exception",
		"experiment",
		"focus",
		"frequency",
		"function",
		"gap_type",
		"gdb_xref",
		"gene_synonym",
		"gene",
		"germline",
		"haplogroup",
		"haplotype",
		"host",
		"identified_by",
		"inference",
		"insertion_seq",
		"isolate",
		"isolation_source",
		"kinetoplast",
		"lab_host",
		"label",
		"lat_lon",
		"linkage_evidence",
		"locus_tag",
		"macronuclear",
		"map",
		"mating_type",
		"metagenome_source",
		"metagenomic",
		"mitochondrion",
		"mobile_element_type",
		"mobile_element",
		"mod_base",
		"mol_type",
		"name",
		"nat_host",
		"ncRNA_class",
		"non_functional",
		"note",
		"number",
		"old_locus_tag",
		"operon",
		"organelle",
		"organism",
		"partial",
		"PCR_conditions",
		"PCR_primers",
		"peptide",
		"phenotype",
		"plasmid",
		"pop_variant",
		"product",
		"protein_id",
		"proviral",
		"pseudo",
		"pseudogene",
		"rearranged",
		"recombination_class",
		"region_name",
		"regulatory_class",
		"replace",
		"ribosomal_slippage",
		"rpt_family",
		"rpt_type",
		"rpt_unit_range",
		"rpt_unit_seq",
		"rpt_unit",
		"satellite",
		"segment",
		"sequenced_mol",
		"serotype",
		"serovar",
		"sex",
		"site_type",
		"specific_host",
		"specimen_voucher",
		"standard_name",
		"strain",
		"structural_class",
		"sub_clone",
		"sub_species",
		"sub_strain",
		"tag_peptide",
		"tissue_lib",
		"tissue_type",
		"trans_splicing",
		"transcript_id",
		"transcription",
		"transgenic",
		"transl_except",
		"transl_table",
		"translation",
		"transposon",
		"type_material",
		"UniProtKB_evidence",
		"usedin",
		"variety",
		"virion",
	}

	// legal INSDSeq XML fields

	insdtags := []string{
		"INSDAltSeqData_items",
		"INSDAltSeqData",
		"INSDAltSeqItem_first-accn",
		"INSDAltSeqItem_gap-comment",
		"INSDAltSeqItem_gap-length",
		"INSDAltSeqItem_gap-linkage",
		"INSDAltSeqItem_gap-type",
		"INSDAltSeqItem_interval",
		"INSDAltSeqItem_isgap",
		"INSDAltSeqItem_isgap@value",
		"INSDAltSeqItem_last-accn",
		"INSDAltSeqItem_value",
		"INSDAltSeqItem",
		"INSDAuthor",
		"INSDComment_paragraphs",
		"INSDComment_type",
		"INSDComment",
		"INSDCommentParagraph",
		"INSDFeature_intervals",
		"INSDFeature_key",
		"INSDFeature_location",
		"INSDFeature_operator",
		"INSDFeature_partial3",
		"INSDFeature_partial3@value",
		"INSDFeature_partial5",
		"INSDFeature_partial5@value",
		"INSDFeature_quals",
		"INSDFeature_xrefs",
		"INSDFeature",
		"INSDFeatureSet_annot-source",
		"INSDFeatureSet_features",
		"INSDFeatureSet",
		"INSDInterval_accession",
		"INSDInterval_from",
		"INSDInterval_interbp",
		"INSDInterval_interbp@value",
		"INSDInterval_iscomp",
		"INSDInterval_iscomp@value",
		"INSDInterval_point",
		"INSDInterval_to",
		"INSDInterval",
		"INSDKeyword",
		"INSDQualifier_name",
		"INSDQualifier_value",
		"INSDQualifier",
		"INSDReference_authors",
		"INSDReference_consortium",
		"INSDReference_journal",
		"INSDReference_position",
		"INSDReference_pubmed",
		"INSDReference_reference",
		"INSDReference_remark",
		"INSDReference_title",
		"INSDReference_xref",
		"INSDReference",
		"INSDSecondary-accn",
		"INSDSeq_accession-version",
		"INSDSeq_alt-seq",
		"INSDSeq_comment-set",
		"INSDSeq_comment",
		"INSDSeq_contig",
		"INSDSeq_create-date",
		"INSDSeq_create-release",
		"INSDSeq_database-reference",
		"INSDSeq_definition",
		"INSDSeq_division",
		"INSDSeq_entry-version",
		"INSDSeq_feature-set",
		"INSDSeq_feature-table",
		"INSDSeq_keywords",
		"INSDSeq_length",
		"INSDSeq_locus",
		"INSDSeq_moltype",
		"INSDSeq_organism",
		"INSDSeq_other-seqids",
		"INSDSeq_primary-accession",
		"INSDSeq_primary",
		"INSDSeq_project",
		"INSDSeq_references",
		"INSDSeq_secondary-accessions",
		"INSDSeq_segment",
		"INSDSeq_sequence",
		"INSDSeq_source-db",
		"INSDSeq_source",
		"INSDSeq_strandedness",
		"INSDSeq_struc-comments",
		"INSDSeq_taxonomy",
		"INSDSeq_topology",
		"INSDSeq_update-date",
		"INSDSeq_update-release",
		"INSDSeq_xrefs",
		"INSDSeq",
		"INSDSeqid",
		"INSDSet",
		"INSDStrucComment_items",
		"INSDStrucComment_name",
		"INSDStrucComment",
		"INSDStrucCommentItem_tag",
		"INSDStrucCommentItem_url",
		"INSDStrucCommentItem_value",
		"INSDStrucCommentItem",
		"INSDXref_dbname",
		"INSDXref_id",
		"INSDXref",
	}

	checkAgainstVocabulary := func(str, objtype string, arry []string) {

		if str == "" || arry == nil {
			return
		}

		// skip past pound, percent, or caret character at beginning of string
		if len(str) > 1 {
			switch str[0] {
			case '#', '%', '^':
				str = str[1:]
			default:
			}
		}

		for _, txt := range arry {
			if str == txt {
				return
			}
			if strings.ToUpper(str) == strings.ToUpper(txt) {
				fmt.Fprintf(os.Stderr, "\nERROR: Incorrect capitalization of '%s' %s, change to '%s'\n", str, objtype, txt)
				os.Exit(1)
			}
		}

		fmt.Fprintf(os.Stderr, "\nERROR: Item '%s' is not a legal -insd %s\n", str, objtype)
		os.Exit(1)
	}

	var acc []string

	max := len(args)
	if max < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract -insd\n")
		os.Exit(1)
	}

	acc = append(acc, "-pattern", "INSDSeq", "-ACCN", "INSDSeq_accession-version")
	printAccn := true

	// collect descriptors

	if strings.HasPrefix(args[0], "INSD") {

		if isPipe {
			acc = append(acc, "-clr", "-pfx", "\\n", "-element", "&ACCN")
			acc = append(acc, "-group", "INSDSeq", "-sep", "|", "-element")
		} else {
			acc = append(acc, "-clr", "-pfx", "\"\\n\"", "-element", "\"&ACCN\"")
			acc = append(acc, "-group", "INSDSeq", "-sep", "\"|\"", "-element")
		}
		printAccn = false

		for {
			if len(args) < 1 {
				return acc
			}
			str := args[0]
			if !strings.HasPrefix(args[0], "INSD") {
				break
			}
			checkAgainstVocabulary(str, "element", insdtags)
			acc = append(acc, str)
			args = args[1:]
		}

	} else if strings.HasPrefix(strings.ToUpper(args[0]), "INSD") {

		// report capitalization or vocabulary failure
		checkAgainstVocabulary(args[0], "element", insdtags)

		// program should not get to this point, but warn and exit anyway
		fmt.Fprintf(os.Stderr, "\nERROR: Item '%s' is not a legal -insd %s\n", args[0], "element")
		os.Exit(1)
	}

	// collect qualifiers

	partial := false
	complete := false

	if args[0] == "+" || args[0] == "complete" {
		complete = true
		args = args[1:]
		max--
	} else if args[0] == "-" || args[0] == "partial" {
		partial = true
		args = args[1:]
		max--
	}

	if max < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: No feature key supplied to xtract -insd\n")
		os.Exit(1)
	}

	acc = append(acc, "-group", "INSDFeature")

	// limit to designated features

	feature := args[0]

	fcmd := "-if"

	// can specify multiple features separated by plus sign (e.g., CDS+mRNA) or comma (e.g., CDS,mRNA)
	plus := strings.Split(feature, "+")
	for _, pls := range plus {
		comma := strings.Split(pls, ",")
		for _, cma := range comma {

			checkAgainstVocabulary(cma, "feature", features)
			acc = append(acc, fcmd, "INSDFeature_key", "-equals", cma)

			fcmd = "-or"
		}
	}

	if max < 2 {
		// still need at least one qualifier even on legal feature
		fmt.Fprintf(os.Stderr, "\nERROR: Feature '%s' must be followed by at least one qualifier\n", feature)
		os.Exit(1)
	}

	args = args[1:]

	if complete {
		acc = append(acc, "-unless", "INSDFeature_partial5", "-or", "INSDFeature_partial3")
	} else if partial {
		acc = append(acc, "-if", "INSDFeature_partial5", "-or", "INSDFeature_partial3")
	}

	if printAccn {
		if isPipe {
			acc = append(acc, "-clr", "-pfx", "\\n", "-element", "&ACCN")
		} else {
			acc = append(acc, "-clr", "-pfx", "\"\\n\"", "-element", "\"&ACCN\"")
		}
	}

	for _, str := range args {
		if strings.HasPrefix(str, "INSD") {

			checkAgainstVocabulary(str, "element", insdtags)
			if isPipe {
				acc = append(acc, "-block", "INSDFeature", "-sep", "|", "-element")
			} else {
				acc = append(acc, "-block", "INSDFeature", "-sep", "\"|\"", "-element")
			}
			acc = append(acc, str)
			if addDash {
				acc = append(acc, "-block", "INSDFeature", "-unless", str)
				if strings.HasSuffix(str, "@value") {
					if isPipe {
						acc = append(acc, "-lbl", "false")
					} else {
						acc = append(acc, "-lbl", "\"false\"")
					}
				} else {
					if isPipe {
						acc = append(acc, "-lbl", "\\-")
					} else {
						acc = append(acc, "-lbl", "\"\\-\"")
					}
				}
			}

		} else if strings.HasPrefix(str, "#INSD") {

			checkAgainstVocabulary(str, "element", insdtags)
			if isPipe {
				acc = append(acc, "-block", "INSDFeature", "-sep", "|", "-element")
				acc = append(acc, str)
			} else {
				acc = append(acc, "-block", "INSDFeature", "-sep", "\"|\"", "-element")
				ql := fmt.Sprintf("\"%s\"", str)
				acc = append(acc, ql)
			}

		} else if strings.HasPrefix(strings.ToUpper(str), "#INSD") || strings.HasPrefix(strings.ToUpper(str), "#INSD") {

			// report capitalization or vocabulary failure
			checkAgainstVocabulary(str, "element", insdtags)

		} else {

			acc = append(acc, "-block", "INSDQualifier")

			checkAgainstVocabulary(str, "qualifier", qualifiers)
			if len(str) > 2 && str[0] == '%' {
				acc = append(acc, "-if", "INSDQualifier_name", "-equals", str[1:])
				if isPipe {
					acc = append(acc, "-element", "%INSDQualifier_value")
				} else {
					acc = append(acc, "-element", "\"%INSDQualifier_value\"")
				}
				if addDash {
					acc = append(acc, "-block", "INSDFeature", "-unless", "INSDQualifier_name", "-equals", str[1:])
					if isPipe {
						acc = append(acc, "-lbl", "\\-")
					} else {
						acc = append(acc, "-lbl", "\"\\-\"")
					}
				}
			} else {
				acc = append(acc, "-if", "INSDQualifier_name", "-equals", str)
				acc = append(acc, "-element", "INSDQualifier_value")
				if addDash {
					acc = append(acc, "-block", "INSDFeature", "-unless", "INSDQualifier_name", "-equals", str)
					if isPipe {
						acc = append(acc, "-lbl", "\\-")
					} else {
						acc = append(acc, "-lbl", "\"\\-\"")
					}
				}
			}
		}
	}

	return acc
}

// HYDRA CITATION MATCHER COMMAND GENERATOR

// ProcessHydra generates extraction commands for NCBI's in-house citation matcher (undocumented)
func ProcessHydra(isPipe bool) []string {

	var acc []string

	// acceptable scores are 0.8 or higher, exact match on "1" rejects low value in scientific notation with minus sign present

	acc = append(acc, "-pattern", "Id")
	acc = append(acc, "-if", "@score", "-equals", "1")
	acc = append(acc, "-or", "@score", "-starts-with", "0.9")
	acc = append(acc, "-or", "@score", "-starts-with", "0.8")
	acc = append(acc, "-element", "Id")

	return acc
}

// COLLECT AND FORMAT REQUESTED XML VALUES

// ExploreElements returns matching element values to callback
func ExploreElements(curr *Node, mask, prnt, match, attrib string, wildcard bool, level int, proc func(string, int)) {

	if curr == nil || proc == nil {
		return
	}

	// **/Object performs deep exploration of recursive data (*/Object also supported)
	deep := false
	if prnt == "**" || prnt == "*" {
		prnt = ""
		deep = true
	}

	// exploreElements recursive definition
	var exploreElements func(curr *Node, skip string, lev int)

	exploreElements = func(curr *Node, skip string, lev int) {

		if !deep && curr.Name == skip {
			// do not explore within recursive object
			return
		}

		// parseAttributes is only run if attribute values are requested in element statements
		parseAttributes := func(attrb string) []string {

			if attrb == "" {
				return nil
			}

			attlen := len(attrb)

			// count equal signs
			num := 0
			for i := 0; i < attlen; i++ {
				if attrb[i] == '=' {
					num += 2
				}
			}
			if num < 1 {
				return nil
			}

			// allocate array of proper size
			arry := make([]string, num)
			if arry == nil {
				return nil
			}

			start := 0
			idx := 0
			itm := 0

			// place tag and value in successive array slots
			for idx < attlen && itm < num {
				ch := attrb[idx]
				if ch == '=' {
					// skip past possible leading blanks
					for start < attlen {
						ch = attrb[start]
						if ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r' || ch == '\f' {
							start++
						} else {
							break
						}
					}
					// =
					arry[itm] = attrb[start:idx]
					itm++
					// skip past equal sign and leading double quote
					idx += 2
					start = idx
				} else if ch == '"' {
					// "
					arry[itm] = attrb[start:idx]
					itm++
					// skip past trailing double quote and (possible) space
					idx += 2
					start = idx
				} else {
					idx++
				}
			}

			return arry
		}

		// wildcard matches any namespace prefix
		if curr.Name == match ||
			(wildcard && strings.HasPrefix(match, ":") && strings.HasSuffix(curr.Name, match)) ||
			(match == "" && attrib != "") {

			if prnt == "" ||
				curr.Parent == prnt ||
				(wildcard && strings.HasPrefix(prnt, ":") && strings.HasSuffix(curr.Parent, prnt)) {

				if attrib != "" {
					if curr.Attributes != "" && curr.Attribs == nil {
						// parse attributes on-the-fly if queried
						curr.Attribs = parseAttributes(curr.Attributes)
					}
					for i := 0; i < len(curr.Attribs)-1; i += 2 {
						// attributes now parsed into array as [ tag, value, tag, value, tag, value, ... ]
						if curr.Attribs[i] == attrib ||
							(wildcard && strings.HasPrefix(attrib, ":") && strings.HasSuffix(curr.Attribs[i], attrib)) {
							proc(curr.Attribs[i+1], level)
							return
						}
					}

				} else if curr.Contents != "" {

					str := curr.Contents[:]

					if HasAmpOrNotASCII(str) {
						// processing of <, >, &, ", and ' characters is now delayed until element contents is requested
						str = html.UnescapeString(str)
					}

					proc(str, level)
					return

				} else if curr.Children != nil {

					// for XML container object, send empty string to callback to increment count
					proc("", level)
					// and continue exploring
				}
			}
		}

		for chld := curr.Children; chld != nil; chld = chld.Next {
			// inner exploration is subject to recursive object exclusion
			exploreElements(chld, mask, lev+1)
		}
	}

	exploreElements(curr, "", level)
}

// PrintSubtree supports compression styles selected by -element "*" through "****"
func PrintSubtree(node *Node, style IndentType, printAttrs bool, proc func(string)) {

	if node == nil || proc == nil {
		return
	}

	// WRAPPED is SUBTREE plus each attribute on its own line
	wrapped := false
	if style == WRAPPED {
		style = SUBTREE
		wrapped = true
	}

	// INDENT is offset by two spaces to allow for parent tag, SUBTREE is not offset
	initial := 1
	if style == SUBTREE {
		style = INDENT
		initial = 0
	}

	// array to speed up indentation
	indentSpaces := []string{
		"",
		"  ",
		"    ",
		"      ",
		"        ",
		"          ",
		"            ",
		"              ",
		"                ",
		"                  ",
	}

	// function to indent a specified number of spaces
	doIndent := func(indt int) {
		i := indt
		for i > 9 {
			proc("                    ")
			i -= 10
		}
		if i < 0 {
			return
		}
		proc(indentSpaces[i])
	}

	// doSubtree recursive definition
	var doSubtree func(*Node, int)

	doSubtree = func(curr *Node, depth int) {

		// suppress if it would be an empty self-closing tag
		if !IsNotJustWhitespace(curr.Attributes) && curr.Contents == "" && curr.Children == nil {
			return
		}

		if style == INDENT {
			doIndent(depth)
		}

		proc("<")
		proc(curr.Name)

		if printAttrs {

			attr := strings.TrimSpace(curr.Attributes)
			attr = CompressRunsOfSpaces(attr)

			if attr != "" {

				if wrapped {

					start := 0
					idx := 0

					attlen := len(attr)

					for idx < attlen {
						ch := attr[idx]
						if ch == '=' {
							str := attr[start:idx]
							proc("\n")
							doIndent(depth)
							proc(" ")
							proc(str)
							// skip past equal sign and leading double quote
							idx += 2
							start = idx
						} else if ch == '"' {
							str := attr[start:idx]
							proc("=\"")
							proc(str)
							proc("\"")
							// skip past trailing double quote and (possible) space
							idx += 2
							start = idx
						} else {
							idx++
						}
					}

					proc("\n")
					doIndent(depth)

				} else {

					proc(" ")
					proc(attr)
				}
			}
		}

		// see if suitable for for self-closing tag
		if curr.Contents == "" && curr.Children == nil {
			proc("/>")
			if style != COMPACT {
				proc("\n")
			}
			return
		}

		proc(">")

		if curr.Contents != "" {

			proc(curr.Contents[:])

		} else {

			if style != COMPACT {
				proc("\n")
			}

			for chld := curr.Children; chld != nil; chld = chld.Next {
				doSubtree(chld, depth+1)
			}

			if style == INDENT {
				i := depth
				for i > 9 {
					proc("                    ")
					i -= 10
				}
				proc(indentSpaces[i])
			}
		}

		proc("<")
		proc("/")
		proc(curr.Name)
		proc(">")

		if style != COMPACT {
			proc("\n")
		}
	}

	doSubtree(node, initial)
}

// ProcessClause handles comma-separated -element arguments
func ProcessClause(curr *Node, stages []*Step, mask, prev, pfx, sfx, sep, def string, status OpType, index, level int, variables map[string]string) (string, bool) {

	if curr == nil || stages == nil {
		return "", false
	}

	// processElement handles individual -element constructs
	processElement := func(acc func(string)) {

		if acc == nil {
			return
		}

		// element names combined with commas are treated as a prefix-separator-suffix group
		for _, stage := range stages {

			stat := stage.Type
			item := stage.Value
			prnt := stage.Parent
			match := stage.Match
			attrib := stage.Attrib
			wildcard := stage.Wild

			// exploreElements is a wrapper for ExploreElements, obtaining most arguments as closures
			exploreElements := func(proc func(string, int)) {
				ExploreElements(curr, mask, prnt, match, attrib, wildcard, level, proc)
			}

			switch stat {
			case ELEMENT, TERMS, WORDS, PAIRS, PHRASE, VALUE, LEN, SUM, MIN, MAX, SUB, AVG, DEV:
				exploreElements(func(str string, lvl int) {
					if str != "" {
						acc(str)
					}
				})
			case FIRST:
				single := ""

				exploreElements(func(str string, lvl int) {
					if single == "" {
						single = str
					}
				})

				if single != "" {
					acc(single)
				}
			case LAST:
				single := ""

				exploreElements(func(str string, lvl int) {
					single = str
				})

				if single != "" {
					acc(single)
				}
			case ENCODE:
				exploreElements(func(str string, lvl int) {
					if str != "" {
						str = html.EscapeString(str)
						acc(str)
					}
				})
			case UPPER:
				exploreElements(func(str string, lvl int) {
					if str != "" {
						str = strings.ToUpper(str)
						acc(str)
					}
				})
			case LOWER:
				exploreElements(func(str string, lvl int) {
					if str != "" {
						str = strings.ToLower(str)
						acc(str)
					}
				})
			case VARIABLE:
				// use value of stored variable
				val, ok := variables[match]
				if ok {
					acc(val)
				}
			case NUM, COUNT:
				count := 0

				exploreElements(func(str string, lvl int) {
					count++
				})

				// number of element objects
				val := strconv.Itoa(count)
				acc(val)
			case LENGTH:
				length := 0

				exploreElements(func(str string, lvl int) {
					length += len(str)
				})

				// length of element strings
				val := strconv.Itoa(length)
				acc(val)
			case DEPTH:
				exploreElements(func(str string, lvl int) {
					// depth of each element in scope
					val := strconv.Itoa(lvl)
					acc(val)
				})
			case INDEX:
				// -element "+" prints index of current XML object
				val := strconv.Itoa(index)
				acc(val)
			case INC:
				// -inc, or component of -0-based, -1-based, or -ucsc
				exploreElements(func(str string, lvl int) {
					if str != "" {
						num, err := strconv.Atoi(str)
						if err == nil {
							// increment value
							num++
							val := strconv.Itoa(num)
							acc(val)
						}
					}
				})
			case DEC:
				// -dec, or component of -0-based, -1-based, or -ucsc
				exploreElements(func(str string, lvl int) {
					if str != "" {
						num, err := strconv.Atoi(str)
						if err == nil {
							// decrement value
							num--
							val := strconv.Itoa(num)
							acc(val)
						}
					}
				})
			case STAR:
				// -element "*" prints current XML subtree on a single line
				style := SINGULARITY
				printAttrs := true

				for _, ch := range item {
					if ch == '*' {
						style++
					} else if ch == '@' {
						printAttrs = false
					}
				}
				if style > WRAPPED {
					style = WRAPPED
				}
				if style < COMPACT {
					style = COMPACT
				}

				var buffer bytes.Buffer

				PrintSubtree(curr, style, printAttrs,
					func(str string) {
						if str != "" {
							buffer.WriteString(str)
						}
					})

				txt := buffer.String()
				if txt != "" {
					acc(txt)
				}
			default:
			}
		}
	}

	ok := false

	// format results in buffer
	var buffer bytes.Buffer

	buffer.WriteString(prev)
	buffer.WriteString(pfx)
	between := ""

	switch status {
	case ELEMENT, ENCODE, UPPER, LOWER, VALUE, NUM, INC, DEC, ZEROBASED, ONEBASED, UCSC:
		processElement(func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(between)
				buffer.WriteString(str)
				between = sep
			}
		})
	case FIRST:
		single := ""

		processElement(func(str string) {
			ok = true
			if single == "" {
				single = str
			}
		})

		if single != "" {
			buffer.WriteString(between)
			buffer.WriteString(single)
			between = sep
		}
	case LAST:
		single := ""

		processElement(func(str string) {
			ok = true
			single = str
		})

		if single != "" {
			buffer.WriteString(between)
			buffer.WriteString(single)
			between = sep
		}
	case TERMS:
		processElement(func(str string) {
			if str != "" {
				words := strings.Fields(str)
				for _, item := range words {
					max := len(item)
					for max > 1 {
						ch := item[max-1]
						if ch != '.' && ch != ',' && ch != ':' && ch != ';' {
							break
						}
						// trim trailing period, comma, colon, and semicolon
						item = item[:max-1]
						// continue checking for runs of punctuation at end
						max--
					}
					ok = true
					buffer.WriteString(between)
					buffer.WriteString(item)
					between = sep
				}
			}
		})
	case WORDS:
		processElement(func(str string) {
			if str != "" {
				words := strings.FieldsFunc(str, func(c rune) bool {
					return !unicode.IsLetter(c) && !unicode.IsNumber(c)
				})
				for _, item := range words {
					item = strings.ToLower(item)
					ok = true
					buffer.WriteString(between)
					buffer.WriteString(item)
					between = sep
				}
			}
		})
	case PAIRS:
		processElement(func(str string) {
			if str != "" {
				words := strings.FieldsFunc(str, func(c rune) bool {
					return !unicode.IsLetter(c) && !unicode.IsNumber(c)
				})
				if len(words) > 1 {
					past := ""
					for _, item := range words {
						item = strings.ToLower(item)
						if isStopWord[item] {
							past = ""
							continue
						}
						if past != "" {
							ok = true
							buffer.WriteString(between)
							buffer.WriteString(past + " " + item)
							between = sep
						}
						past = item
					}
				}
			}
		})
	case PHRASE:
		processElement(func(str string) {
			if str != "" {
				words := strings.FieldsFunc(str, func(c rune) bool {
					return !unicode.IsLetter(c) && !unicode.IsNumber(c)
				})
				buffer.WriteString("<Phrase>")
				if len(words) > 0 {
					buffer.WriteString("<Terms>")
					for _, item := range words {
						item = strings.ToLower(item)
						if isStopWord[item] {
							continue
						}
						ok = true
						buffer.WriteString("<Term>")
						buffer.WriteString(item)
						buffer.WriteString("</Term>")
					}
					buffer.WriteString("</Terms>")
				}
				if len(words) > 1 {
					buffer.WriteString("<Pairs>")
					past := ""
					for _, item := range words {
						item = strings.ToLower(item)
						if isStopWord[item] {
							past = ""
							continue
						}
						if past != "" {
							ok = true
							buffer.WriteString("<Pair>")
							buffer.WriteString(past + " " + item)
							buffer.WriteString("</Pair>")
							between = sep
						}
						past = item
					}
					buffer.WriteString("</Pairs>")
				}
				buffer.WriteString("</Phrase>")
			}
		})
	case LEN:
		length := 0

		processElement(func(str string) {
			ok = true
			length += len(str)
		})

		// length of element strings
		val := strconv.Itoa(length)
		buffer.WriteString(between)
		buffer.WriteString(val)
		between = sep
	case SUM:
		sum := 0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				sum += value
				ok = true
			}
		})

		if ok {
			// sum of element values
			val := strconv.Itoa(sum)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}
	case MIN:
		min := 0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				if !ok || value < min {
					min = value
				}
				ok = true
			}
		})

		if ok {
			// minimum of element values
			val := strconv.Itoa(min)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}
	case MAX:
		max := 0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				if !ok || value > max {
					max = value
				}
				ok = true
			}
		})

		if ok {
			// maximum of element values
			val := strconv.Itoa(max)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}
	case SUB:
		first := 0
		second := 0
		count := 0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				count++
				if count == 1 {
					first = value
				} else if count == 2 {
					second = value
				}
			}
		})

		if count == 2 {
			// must have exactly 2 elements
			ok = true
			// difference of element values
			val := strconv.Itoa(first - second)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}
	case AVG:
		sum := 0
		count := 0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				sum += value
				count++
				ok = true
			}
		})

		if ok {
			// average of element values
			avg := int(float64(sum) / float64(count))
			val := strconv.Itoa(avg)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}
	case DEV:
		count := 0
		mean := 0.0
		m2 := 0.0

		processElement(func(str string) {
			value, err := strconv.Atoi(str)
			if err == nil {
				// Welford algorithm for one-pass standard deviation
				count++
				x := float64(value)
				delta := x - mean
				mean += delta / float64(count)
				m2 += delta * (x - mean)
			}
		})

		if count > 1 {
			// must have at least 2 elements
			ok = true
			// standard deviation of element values
			vrc := m2 / float64(count-1)
			dev := int(math.Sqrt(vrc))
			val := strconv.Itoa(dev)
			buffer.WriteString(between)
			buffer.WriteString(val)
			between = sep
		}
	default:
	}

	// use default value if nothing written
	if !ok && def != "" {
		ok = true
		buffer.WriteString(def)
	}

	buffer.WriteString(sfx)

	if !ok {
		return "", false
	}

	txt := buffer.String()

	return txt, true
}

// ProcessInstructions performs extraction commands on a subset of XML
func ProcessInstructions(commands []*Operation, curr *Node, mask, tab, ret string, index, level int, variables map[string]string, accum func(string)) (string, string) {

	if accum == nil {
		return tab, ret
	}

	sep := "\t"
	pfx := ""
	sfx := ""

	def := ""

	col := "\t"
	lin := "\n"

	varname := ""

	// process commands
	for _, op := range commands {

		str := op.Value

		switch op.Type {
		case ELEMENT, FIRST, LAST, ENCODE, UPPER, LOWER, TERMS, WORDS, PAIRS, PHRASE, NUM, LEN, SUM, MIN, MAX, INC, DEC, SUB, AVG, DEV, ZEROBASED, ONEBASED, UCSC:
			txt, ok := ProcessClause(curr, op.Stages, mask, tab, pfx, sfx, sep, def, op.Type, index, level, variables)
			if ok {
				tab = col
				ret = lin
				accum(txt)
			}
		case TAB:
			col = str
		case RET:
			lin = str
		case PFX:
			pfx = str
		case SFX:
			sfx = str
		case SEP:
			sep = str
		case LBL:
			lbl := str
			accum(tab)
			accum(lbl)
			tab = col
			ret = lin
		case PFC:
			// preface clears previous tab and sets prefix in one command
			pfx = str
			fallthrough
		case CLR:
			// clear previous tab after the fact
			tab = ""
		case RST:
			pfx = ""
			sfx = ""
			sep = "\t"
			def = ""
		case DEF:
			def = str
		case VARIABLE:
			varname = str
		case VALUE:
			length := len(str)
			if length > 1 && str[0] == '(' && str[length-1] == ')' {
				// set variable from literal text inside parentheses, e.g., -COM "(, )"
				variables[varname] = str[1 : length-1]
				// -if "&VARIABLE" will succeed if set to blank with empty parentheses "()"
			} else if str == "" {
				// -if "&VARIABLE" will fail if initialized with empty string ""
				delete(variables, varname)
			} else {
				txt, ok := ProcessClause(curr, op.Stages, mask, "", pfx, sfx, sep, def, op.Type, index, level, variables)
				if ok {
					variables[varname] = txt
				}
			}
			varname = ""
		default:
		}
	}

	return tab, ret
}

// CONDITIONAL EXECUTION USES -if AND -unless STATEMENT, WITH SUPPORT FOR DEPRECATED -match AND -avoid STATEMENTS

// ConditionsAreSatisfied tests a set of conditions to determine if extraction should proceed
func ConditionsAreSatisfied(conditions []*Operation, curr *Node, mask string, index, level int, variables map[string]string) bool {

	if curr == nil {
		return false
	}

	required := 0
	observed := 0
	forbidden := 0
	isMatch := false
	isAvoid := false

	// function to test string or numeric constraints
	testConstraint := func(str string, constraint *Step) bool {

		if str == "" || constraint == nil {
			return false
		}

		val := constraint.Value
		stat := constraint.Type

		switch stat {
		case EQUALS, CONTAINS, STARTSWITH, ENDSWITH, ISNOT:
			// substring test on element values
			str = strings.ToUpper(str)
			val = strings.ToUpper(val)

			switch stat {
			case EQUALS:
				if str == val {
					return true
				}
			case CONTAINS:
				if strings.Contains(str, val) {
					return true
				}
			case STARTSWITH:
				if strings.HasPrefix(str, val) {
					return true
				}
			case ENDSWITH:
				if strings.HasSuffix(str, val) {
					return true
				}
			case ISNOT:
				if str != val {
					return true
				}
			default:
			}
		case GT, GE, LT, LE, EQ, NE:
			// second argument of numeric test can be element specifier
			if constraint.Parent != "" || constraint.Match != "" || constraint.Attrib != "" {
				ExploreElements(curr, mask, constraint.Parent, constraint.Match, constraint.Attrib, constraint.Wild, level, func(str string, lvl int) {
					if str != "" {
						_, errz := strconv.Atoi(str)
						if errz == nil {
							val = str
						}
					}
				})
			}

			// numeric tests on element values
			x, errx := strconv.Atoi(str)
			y, erry := strconv.Atoi(val)

			// both arguments must resolve to integers
			if errx != nil || erry != nil {
				return false
			}

			switch stat {
			case GT:
				if x > y {
					return true
				}
			case GE:
				if x >= y {
					return true
				}
			case LT:
				if x < y {
					return true
				}
			case LE:
				if x <= y {
					return true
				}
			case EQ:
				if x == y {
					return true
				}
			case NE:
				if x != y {
					return true
				}
			default:
			}
		default:
		}

		return false
	}

	// matchFound tests individual conditions
	matchFound := func(stages []*Step) bool {

		if stages == nil || len(stages) < 1 {
			return false
		}

		stage := stages[0]

		var constraint *Step

		if len(stages) > 1 {
			constraint = stages[1]
		}

		status := stage.Type
		prnt := stage.Parent
		match := stage.Match
		attrib := stage.Attrib
		wildcard := stage.Wild

		found := false
		number := ""

		// exploreElements is a wrapper for ExploreElements, obtaining most arguments as closures
		exploreElements := func(proc func(string, int)) {
			ExploreElements(curr, mask, prnt, match, attrib, wildcard, level, proc)
		}

		switch status {
		case ELEMENT:
			exploreElements(func(str string, lvl int) {
				// match to XML container object sends empty string, so do not check for str != "" here
				// test every selected element individually if value is specified
				if constraint == nil || testConstraint(str, constraint) {
					found = true
				}
			})
		case VARIABLE:
			// use value of stored variable
			str, ok := variables[match]
			if ok {
				//  -if &VARIABLE -equals VALUE is the supported construct
				if constraint == nil || testConstraint(str, constraint) {
					found = true
				}
			}
		case COUNT:
			count := 0

			exploreElements(func(str string, lvl int) {
				count++
				found = true
			})

			// number of element objects
			number = strconv.Itoa(count)
		case LENGTH:
			length := 0

			exploreElements(func(str string, lvl int) {
				length += len(str)
				found = true
			})

			// length of element strings
			number = strconv.Itoa(length)
		case DEPTH:
			depth := 0

			exploreElements(func(str string, lvl int) {
				depth = lvl
				found = true
			})

			// depth of last element in scope
			number = strconv.Itoa(depth)
		case INDEX:
			// index of explored parent object
			number = strconv.Itoa(index)
			found = true
		default:
		}

		if number == "" {
			return found
		}

		if constraint == nil || testConstraint(number, constraint) {
			return true
		}

		return false
	}

	// test conditional arguments
	for _, op := range conditions {

		switch op.Type {
		// -if tests for presence of element (deprecated -match can test element:value)
		case IF, MATCH:
			// checking for failure here allows for multiple -if [ -and / -or ] clauses
			if isMatch && observed < required {
				return false
			}
			if isAvoid && forbidden > 0 {
				return false
			}
			required = 0
			observed = 0
			forbidden = 0
			isMatch = true
			isAvoid = false
			// continue on to next two cases
			fallthrough
		case AND:
			required++
			// continue on to next case
			fallthrough
		case OR:
			if matchFound(op.Stages) {
				observed++
				// record presence of forbidden element if in -unless clause
				forbidden++
			}
		// -unless tests for absence of element, or presence but with failure of subsequent value test (deprecated -avoid can test element:value)
		case UNLESS, AVOID:
			if isMatch && observed < required {
				return false
			}
			if isAvoid && forbidden > 0 {
				return false
			}
			required = 0
			observed = 0
			forbidden = 0
			isMatch = false
			isAvoid = true
			if matchFound(op.Stages) {
				forbidden++
			}
		default:
		}
	}

	if isMatch && observed < required {
		return false
	}
	if isAvoid && forbidden > 0 {
		return false
	}

	return true
}

// RECURSIVELY PROCESS EXPLORATION COMMANDS AND XML DATA STRUCTURE

// ProcessCommands visits XML nodes, performs conditional tests, and executes data extraction instructions
func ProcessCommands(cmds *Block, curr *Node, tab, ret string, index, level int, variables map[string]string, accum func(string)) (string, string) {

	if accum == nil {
		return tab, ret
	}

	prnt := cmds.Parent
	match := cmds.Match

	// leading colon indicates namespace prefix wildcard
	wildcard := false
	if strings.HasPrefix(prnt, ":") || strings.HasPrefix(match, ":") {
		wildcard = true
	}

	// **/Object performs deep exploration of recursive data
	deep := false
	if prnt == "**" {
		prnt = "*"
		deep = true
	}

	// closure passes local variables to callback, which can modify caller tab and ret values
	processNode := func(node *Node, idx, lvl int) {

		// apply -if or -unless tests
		if ConditionsAreSatisfied(cmds.Conditions, node, match, idx, lvl, variables) {

			// execute data extraction commands
			if len(cmds.Commands) > 0 {
				tab, ret = ProcessInstructions(cmds.Commands, node, match, tab, ret, idx, lvl, variables, accum)
			}

			// process sub commands on child node
			for _, sub := range cmds.Subtasks {
				tab, ret = ProcessCommands(sub, node, tab, ret, 1, lvl, variables, accum)
			}

		} else {

			// execute commands after -else statement
			if len(cmds.Failure) > 0 {
				tab, ret = ProcessInstructions(cmds.Failure, node, match, tab, ret, idx, lvl, variables, accum)
			}
		}
	}

	// exploreNodes recursive definition
	var exploreNodes func(*Node, int, int, func(*Node, int, int)) int

	// exploreNodes visits all nodes that match the selection criteria
	exploreNodes = func(curr *Node, indx, levl int, proc func(*Node, int, int)) int {

		if curr == nil || proc == nil {
			return indx
		}

		// match is "*" for heterogeneous data constructs, e.g., -group PubmedArticleSet/*
		// wildcard matches any namespace prefix
		if curr.Name == match ||
			match == "*" ||
			(wildcard && strings.HasPrefix(match, ":") && strings.HasSuffix(curr.Name, match)) {

			if prnt == "" ||
				curr.Parent == prnt ||
				(wildcard && strings.HasPrefix(prnt, ":") && strings.HasSuffix(curr.Parent, prnt)) {

				proc(curr, indx, levl)
				indx++

				if !deep {
					// do not explore within recursive object
					return indx
				}
			}
		}

		// clearing prnt "*" now allows nested exploration within recursive data, e.g., -pattern Taxon -block */Taxon
		if prnt == "*" {
			prnt = ""
		}

		// explore child nodes
		for chld := curr.Children; chld != nil; chld = chld.Next {
			indx = exploreNodes(chld, indx, levl+1, proc)
		}

		return indx
	}

	// apply -position test

	if cmds.Position == "" {

		exploreNodes(curr, index, level, processNode)

	} else {

		var single *Node
		lev := 0
		ind := 0

		if cmds.Position == "first" {

			exploreNodes(curr, index, level,
				func(node *Node, idx, lvl int) {
					if single == nil {
						single = node
						ind = idx
						lev = lvl
					}
				})

		} else if cmds.Position == "last" {

			exploreNodes(curr, index, level,
				func(node *Node, idx, lvl int) {
					single = node
					ind = idx
					lev = lvl
				})

		} else {

			// use numeric position
			number, err := strconv.Atoi(cmds.Position)
			if err == nil {

				pos := 0

				exploreNodes(curr, index, level,
					func(node *Node, idx, lvl int) {
						pos++
						if pos == number {
							single = node
							ind = idx
							lev = lvl
						}
					})

			} else {

				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized position '%s'\n", cmds.Position)
				os.Exit(1)
			}
		}

		if single != nil {
			processNode(single, ind, lev)
		}
	}

	return tab, ret
}

// PROCESS ONE XML COMPONENT RECORD

// ParseXML combines the tokenizer and parser to work on a single partitioned XML string
func ParseXML(Text, parent string, tbls *Tables) (*Node, bool) {

	if Text == "" || tbls == nil {
		return nil, false
	}

	// node farm variables
	FarmPos := 0
	FarmMax := tbls.FarmSize
	FarmItems := make([]Node, FarmMax)

	// function to allocate multiple nodes in a large array for memory management efficiency
	nextNode := func(strt, attr, prnt string) *Node {

		// if farm array slots used up, allocate new array
		if FarmPos >= FarmMax {
			FarmItems = make([]Node, FarmMax)
			FarmPos = 0
		}

		if FarmItems == nil {
			return nil
		}

		// take node from next available slot in farm array
		node := &FarmItems[FarmPos]

		node.Name = strt[:]
		node.Attributes = attr[:]
		node.Parent = prnt[:]

		FarmPos++

		return node
	}

	// token parser variables
	Txtlen := len(Text)
	Idx := 0

	// function to get next XML token
	nextToken := func(idx int) (TagType, string, string, int) {

		// lookup table array pointers
		inBlank := &tbls.InBlank
		inFirst := &tbls.InFirst
		inElement := &tbls.InElement

		text := Text[:]
		txtlen := Txtlen

		// XML string ends with > character, acts as sentinel to check if past end of text
		if idx >= txtlen {
			// signal end of XML string
			return ISCLOSED, "", "", 0
		}

		// skip past leading blanks
		ch := text[idx]
		for inBlank[ch] {
			idx++
			ch = text[idx]
		}

		start := idx

		if ch == '<' {

			// at start of element
			idx++
			ch = text[idx]

			// check for legal first character of element
			if inFirst[ch] {

				// read element name
				start = idx
				idx++

				ch = text[idx]
				for inElement[ch] {
					idx++
					ch = text[idx]
				}

				str := text[start:idx]

				switch ch {
				case '>':
					// end of element
					idx++

					return STARTTAG, str[:], "", idx
				case '/':
					// self-closing element without attributes
					idx++
					ch = text[idx]
					if ch != '>' {
						fmt.Fprintf(os.Stderr, "\nSelf-closing element missing right angle bracket\n")
					}
					idx++

					return SELFTAG, str[:], "", idx
				case ' ', '\t', '\n', '\r', '\f':
					// attributes
					idx++
					start = idx
					ch = text[idx]
					for ch != '<' && ch != '>' {
						idx++
						ch = text[idx]
					}
					if ch != '>' {
						fmt.Fprintf(os.Stderr, "\nAttributes not followed by right angle bracket\n")
					}
					if text[idx-1] == '/' {
						// self-closing
						atr := text[start : idx-1]
						idx++
						return SELFTAG, str[:], atr[:], idx
					}
					atr := text[start:idx]
					idx++
					return STARTTAG, str[:], atr[:], idx
				default:
					fmt.Fprintf(os.Stderr, "\nUnexpected punctuation in XML element\n")
					return STARTTAG, str[:], "", idx
				}

			} else {

				// punctuation character immediately after first angle bracket
				switch ch {
				case '/':
					// at start of end tag
					idx++
					start = idx
					ch = text[idx]
					// expect legal first character of element
					if inFirst[ch] {
						idx++
						ch = text[idx]
						for inElement[ch] {
							idx++
							ch = text[idx]
						}
						str := text[start:idx]
						if ch != '>' {
							fmt.Fprintf(os.Stderr, "\nUnexpected characters after end element name\n")
						}
						idx++

						return STOPTAG, str[:], "", idx
					}
					fmt.Fprintf(os.Stderr, "\nUnexpected punctuation in XML element\n")
				case '?':
					// skip ?xml and ?processing instructions
					idx++
					ch = text[idx]
					for ch != '>' {
						idx++
						ch = text[idx]
					}
					idx++
				case '!':
					// skip !DOCTYPE, !comment, and ![CDATA[
					idx++
					start = idx
					ch = text[idx]
					which := NOTAG
					skipTo := ""
					if ch == '[' && strings.HasPrefix(text[idx:], "[CDATA[") {
						which = CDATATAG
						skipTo = "]]>"
						start += 7
					} else if ch == '-' && strings.HasPrefix(text[idx:], "--") {
						which = COMMENTTAG
						skipTo = "-->"
						start += 2
					}
					if which != NOTAG && skipTo != "" {
						// CDATA or comment block may contain internal angle brackets
						found := strings.Index(text[idx:], skipTo)
						if found < 0 {
							// string stops in middle of CDATA or comment
							return ISCLOSED, "", "", idx
						}
						// adjust position past end of CDATA or comment
						idx += found + len(skipTo)
					} else {
						// otherwise just skip to next right angle bracket
						for ch != '>' {
							idx++
							ch = text[idx]
						}
						idx++
					}
				default:
					fmt.Fprintf(os.Stderr, "\nUnexpected punctuation in XML element\n")
				}
			}

		} else if ch != '>' {

			// at start of contents
			start = idx

			// find end of contents
			for ch != '<' && ch != '>' {
				idx++
				ch = text[idx]
			}

			// trim back past trailing blanks
			lst := idx - 1
			ch = text[lst]
			for inBlank[ch] && lst > start {
				lst--
				ch = text[lst]
			}

			str := text[start : lst+1]

			return CONTENTTAG, str[:], "", idx
		}

		return NOTAG, "", "", idx
	}

	// parseLevel recursive definition
	var parseLevel func(string, string, string) (*Node, bool)

	parseLevel = func(strt, attr, prnt string) (*Node, bool) {

		ok := true

		// obtain next node from farm
		node := nextNode(strt, attr, prnt)
		if node == nil {
			return nil, false
		}

		var lastNode *Node

		for {
			tag, name, attr, idx := nextToken(Idx)
			if tag == ISCLOSED {
				break
			}
			Idx = idx

			switch tag {
			case STARTTAG:
				// read sub tree
				obj, ok := parseLevel(name, attr, node.Name)
				if !ok {
					break
				}

				// adding next child to end of linked list gives better performance than appending to slice of nodes
				if node.Children == nil {
					node.Children = obj
				}
				if lastNode != nil {
					lastNode.Next = obj
				}
				lastNode = obj
			case STOPTAG:
				// pop out of recursive call
				return node, ok
			case CONTENTTAG:
				node.Contents = name
			case SELFTAG:
				if attr == "" {
					// ignore if self-closing tag has no attributes
					continue
				}

				// self-closing tag has no contents, just create child node
				obj := nextNode(name, attr, node.Name)

				if node.Children == nil {
					node.Children = obj
				}
				if lastNode != nil {
					lastNode.Next = obj
				}
				lastNode = obj
				// continue on same level
			default:
			}
		}

		return node, ok
	}

	for {
		tag, name, attr, idx := nextToken(Idx)
		if tag == ISCLOSED {
			break
		}

		Idx = idx

		if tag != STARTTAG {
			continue
		}

		// call recursive function from beginning of XML
		return parseLevel(name, attr, parent)
	}

	return nil, false
}

// ProcessQuery calls XML combined tokenizer parser on a partitioned string
func ProcessQuery(text, parent, hd, tl string, index int, cmds *Block, tbls *Tables) string {

	if text == "" || cmds == nil || tbls == nil {
		return ""
	}

	// exit from function will collect garbage of node structure for current XML object
	pat, ok := ParseXML(text, parent, tbls)

	if !ok {
		return ""
	}

	// exit from function will also free map of recorded variables for current -pattern
	variables := make(map[string]string)

	var buffer bytes.Buffer

	ok = false

	if hd != "" {
		buffer.WriteString(hd[:])
	}

	// start processing at top of command tree and top of XML subregion selected by -pattern
	_, ret := ProcessCommands(cmds, pat, "", "", index, 1, variables,
		func(str string) {
			if str != "" {
				ok = true
				buffer.WriteString(str)
			}
		})

	if tl != "" {
		buffer.WriteString(tl[:])
	}

	if ret != "" {
		ok = true
		buffer.WriteString(ret)
	}

	txt := buffer.String()

	// remove leading newline (-insd -pfx artifact)
	if txt != "" && txt[0] == '\n' {
		txt = txt[1:]
	}

	if !ok {
		return ""
	}

	// return consolidated result string
	return txt
}

// UNSHUFFLER USES HEAP TO RESTORE OUTPUT OF MULTIPLE CONSUMERS TO ORIGINAL RECORD ORDER

type Extract struct {
	Index int
	Text  string
}

type ExtractHeap []Extract

// methods that satisfy heap.Interface
func (h ExtractHeap) Len() int {
	return len(h)
}
func (h ExtractHeap) Less(i, j int) bool {
	return h[i].Index < h[j].Index
}
func (h ExtractHeap) Swap(i, j int) {
	h[i], h[j] = h[j], h[i]
}
func (h *ExtractHeap) Push(x interface{}) {
	*h = append(*h, x.(Extract))
}
func (h *ExtractHeap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

// CONCURRENT CONSUMER GOROUTINES PARSE AND PROCESS PARTITIONED XML OBJECTS

// ReadBlocks -> SplitPattern => StreamTokens => ParseXML => ProcessQuery -> MergeResults

// XMLProducer sends partitioned XML strings through channel
func XMLProducer(pat, star string, rdr *XMLReader, out chan<- Extract) {

	// close channel when all records have been processed, so consumers can range over channel
	defer close(out)

	// partition all input by pattern and send XML substring to available consumer through channel
	PartitionPattern(pat, star, rdr,
		func(rec int, str string) {
			out <- Extract{rec, str}
		})
}

func CreateProducer(pat, star string, rdr *XMLReader, tbls *Tables) <-chan Extract {

	if tbls == nil {
		return nil
	}

	out := make(chan Extract, tbls.ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create producer channel\n")
		os.Exit(1)
	}

	// launch single producer goroutine
	go XMLProducer(pat, star, rdr, out)

	return out
}

// XMLConsumer reads partitioned XML from channel and calls parser for processing
func XMLConsumer(cmds *Block, tbls *Tables, parent, hd, tl string, wg *sync.WaitGroup, inp <-chan Extract, out chan<- Extract) {

	// report when this consumer has no more records to process
	defer wg.Done()

	// read partitioned XML from producer channel
	for ext := range inp {

		idx := ext.Index
		text := ext.Text

		if text == "" {
			// should never see empty input data
			out <- Extract{idx, text}
			continue
		}

		str := ProcessQuery(text[:], parent, hd, tl, idx, cmds, tbls)

		// send even if empty to get all record counts for reordering
		out <- Extract{idx, str}
	}
}

func CreateConsumers(cmds *Block, tbls *Tables, parent, hd, tl string, numServers int, inp <-chan Extract) <-chan Extract {

	if tbls == nil {
		return nil
	}

	out := make(chan Extract, tbls.ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create consumer channel\n")
		os.Exit(1)
	}

	var wg sync.WaitGroup

	// launch multiple consumer goroutines
	for i := 0; i < numServers; i++ {
		wg.Add(1)
		go XMLConsumer(cmds, tbls, parent, hd, tl, &wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all consumers are done, then close single output channel, so unshuffler can range over channel
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

// MAIN FUNCTION

// e.g., xtract -pattern PubmedArticle -element MedlineCitation/PMID -block Author -sep " " -element Initials,LastName

func main() {

	// skip past executable name
	args := os.Args[1:]

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: No command-line arguments supplied to xtract\n")
		os.Exit(1)
	}

	// CONCURRENCY, CLEANUP, AND DEBUGGING FLAGS

	// do these first because -defcpu and -maxcpu can be sent from wrapper before other arguments

	ncpu := runtime.NumCPU()
	if ncpu < 1 {
		ncpu = 1
	}

	// wrapper can limit maximum number of processors to use (undocumented)
	maxProcs := ncpu
	defProcs := 0

	// concurrent performance tuning parameters, can be overridden by -proc and -cons
	numProcs := 0
	serverRatio := 4

	// number of servers usually calculated by -cons server ratio, but can be overridden by -serv
	numServers := 0

	// number of channels usually equals number of servers, but can be overridden by -chan
	chanDepth := 0

	// miscellaneous tuning parameters
	heapSize := 16
	farmSize := 64

	// garbage collector control can be set by environment variable or default value with -gogc 0
	goGc := 600

	// XML data cleanup
	doCompress := false
	doCleanup := false

	// read data from file instead of stdin
	fileName := ""

	// debugging
	dbug := false
	mpty := false
	indx := false
	stts := false
	timr := false

	// profiling
	prfl := false

	// alternative source of sample record, processed a designated number of times, looping for each -proc from 1 to nCPU (undocumented)
	testCount := 0
	testType := ""
	testString := ""

	// repeat the specified extraction 5 times for each -proc from 1 to nCPU
	trial := false

	// function to get numeric value
	getNumericArg := func(name string, zer, min, max int) int {

		if len(args) < 2 {
			fmt.Fprintf(os.Stderr, "\nERROR: %s is missing\n", name)
			os.Exit(1)
		}
		value, err := strconv.Atoi(args[1])
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: %s (%s) is not an integer\n", name, args[1])
			os.Exit(1)
		}
		// skip past first of two arguments
		args = args[1:]

		// special case for argument value of 0
		if value < 1 {
			return zer
		}
		// limit value to between specified minimum and maximum
		if value < min {
			return min
		}
		if value > max {
			return max
		}
		return value
	}

	inSwitch := true

	// get concurrency, cleanup, and debugging flags in any order
	for {

		inSwitch = true

		switch args[0] {
		// concurrency override arguments can be passed in by local wrapper script (undocumented)
		case "-maxcpu":
			maxProcs = getNumericArg("Maximum number of processors", 1, 1, ncpu)
		case "-defcpu":
			defProcs = getNumericArg("Default number of processors", ncpu, 1, ncpu)
		// performance tuning flags
		case "-proc":
			numProcs = getNumericArg("Number of processors", ncpu, 1, ncpu)
		case "-cons":
			serverRatio = getNumericArg("Parser to processor ratio", 4, 1, 32)
		case "-serv":
			numServers = getNumericArg("Concurrent parser count", 0, ncpu, 128)
		case "-chan":
			chanDepth = getNumericArg("Communication channel depth", 0, ncpu, 128)
		case "-heap":
			heapSize = getNumericArg("Unshuffler heap size", 8, 8, 64)
		case "-farm":
			farmSize = getNumericArg("Node buffer length", 4, 4, 2048)
		case "-gogc":
			goGc = getNumericArg("Garbage collection percentage", 0, 100, 1000)
		// read data from file
		case "-input":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Input file name is missing\n")
				os.Exit(1)
			}
			fileName = args[1]
			// skip past first of two arguments
			args = args[1:]
		// data cleanup flags
		case "-compress":
			doCompress = true
		case "-cleanup":
			doCleanup = true
		// debugging flags
		case "-debug":
			dbug = true
		case "-empty":
			mpty = true
		case "-index":
			indx = true
		case "-stats", "-stat":
			stts = true
		case "-timer":
			timr = true
		case "-profile":
			prfl = true
		case "-trial":
			trial = true
		case "-test":
			testCount = getNumericArg("Test data counter", 0, 0, 1000000)
			if len(args) > 1 {
				next := args[1]
				// if next argument is not another flag
				if next != "" && next[0] != '-' {
					// get optional data source specifier
					testType = next
					// skip past second of three arguments
					args = args[1:]
				}
			}
		default:
			// if not any of the controls, set flag to break out of for loop
			inSwitch = false
		}

		if !inSwitch {
			break
		}

		// skip past argument
		args = args[1:]

		if len(args) < 1 {
			break
		}
	}

	// reality checks on number of processors to use
	// performance degrades if capacity is above maximum number of partitions per second (context switching?)
	if numProcs == 0 {
		if defProcs > 0 {
			numProcs = defProcs
		} else {
			// best performance measurement with current code is obtained when 4 to 6 processors are assigned,
			// varying slightly among queries on PubmedArticle, gene DocumentSummary, and INSDSeq sequence records
			numProcs = 4
		}
	}
	if numProcs > ncpu {
		numProcs = ncpu
	}
	if numProcs > maxProcs {
		numProcs = maxProcs
	}

	// allow simultaneous threads for multiplexed go routines
	runtime.GOMAXPROCS(numProcs)

	// adjust garbage collection target percentage
	if goGc >= 100 {
		debug.SetGCPercent(goGc)
	}

	// explicit -serv argument overrides -cons ratio
	if numServers > 0 {
		serverRatio = numServers / numProcs
		// if numServers / numProcs is not a whole number, do not print serverRatio in -stats
		if numServers != numProcs*serverRatio {
			serverRatio = 0
		}
	} else {
		numServers = numProcs * serverRatio
	}
	// server limits
	if numServers > 128 {
		numServers = 128
	} else if numServers < 1 {
		numServers = numProcs
	}

	// explicit -chan argument overrides default to number of servers
	if chanDepth == 0 {
		chanDepth = numServers
	}

	// -stats prints number of CPUs and performance tuning values if no other arguments (undocumented)
	if stts && len(args) < 1 {

		fmt.Fprintf(os.Stderr, "CPUs %d\n", ncpu)
		fmt.Fprintf(os.Stderr, "Proc %d\n", numProcs)
		if serverRatio > 0 {
			fmt.Fprintf(os.Stderr, "Cons %d\n", serverRatio)
		}
		fmt.Fprintf(os.Stderr, "Serv %d\n", numServers)
		fmt.Fprintf(os.Stderr, "Chan %d\n", chanDepth)
		fmt.Fprintf(os.Stderr, "Heap %d\n", heapSize)
		fmt.Fprintf(os.Stderr, "Farm %d\n", farmSize)
		if goGc >= 100 {
			fmt.Fprintf(os.Stderr, "Gogc %d\n", goGc)
		}
		fmt.Fprintf(os.Stderr, "\n")

		return
	}

	// -test N [pubmed|protein|insd|gene] repeats simple query on local XML to measure performance independent of stdin (undocumented)
	if testCount > 0 {

		var acc []string

		// select internal XML data source
		switch testType {
		case "pubmed":
			testString = pubMedArtSample
		case "protein", "sequence":
			testString = insdSeqSample
		case "insd":
			testString = insdSeqSample
		case "gene", "docsum":
			testString = geneDocSumSample
		default:
			testString = pubMedArtSample
		}

		// default commands if no other arguments
		if len(args) < 1 {
			switch testType {
			case "pubmed":
				acc = append(acc, "-pattern", "PubmedArticle", "-element", "LastName")
			case "protein", "sequence":
				acc = append(acc, "-pattern", "INSDSeq", "-element", "INSDSeq_accession-version")
			case "insd":
				acc = append(acc, "-insd", "mat_peptide", "%peptide", "product", "peptide")
			case "gene", "docsum":
				acc = append(acc, "-pattern", "DocumentSummary", "-element", "Name")
			default:
				acc = append(acc, "-pattern", "PubmedArticle", "-element", "LastName")
			}
		}

		// otherwise use remaining arguments for extraction commands
		for len(args) > 0 {
			acc = append(acc, args[0])
			args = args[1:]
		}

		args = acc
	}

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract\n")
		os.Exit(1)
	}

	// DOCUMENTATION COMMANDS

	inSwitch = true

	switch args[0] {
	case "-version":
		fmt.Printf("%s\n", xtractVersion)
	case "-help", "-extras", "-extra":
		fmt.Printf("xtract %s\n%s\n", xtractVersion, xtractHelp)
	case "-examples", "-example", "-scripts", "-script":
		fmt.Printf("xtract %s\n%s\n", xtractVersion, xtractExamples)
	case "-internal", "-internals":
		fmt.Printf("xtract %s\n%s\n", xtractVersion, xtractInternal)
	case "-sample", "-samples":
		// -sample [pubmed|protein|gene] sends specified sample record to stdout (undocumented)
		if len(args) > 1 {
			testType = args[1]
		}
		switch testType {
		case "pubmed":
			fmt.Printf("%s\n", pubMedArtSample)
		case "protein", "sequence", "insd":
			fmt.Printf("%s\n", insdSeqSample)
		case "gene", "docsum":
			fmt.Printf("%s\n", geneDocSumSample)
		default:
			fmt.Printf("%s\n", pubMedArtSample)
		}
	case "-keys":
		fmt.Printf("%s\n", keyboardShortcuts)
	case "-unix":
		fmt.Printf("%s\n", unixCommands)
	default:
		// if not any of the documentation commands, keep going
		inSwitch = false
	}

	if inSwitch {
		return
	}

	// INITIALIZE TABLES

	tbls := InitTables()
	if tbls == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Problem creating token streamer lookup tables\n")
		os.Exit(1)
	}

	// additional fields passed in master table
	tbls.ChanDepth = chanDepth
	tbls.FarmSize = farmSize

	// FILE NAME CAN BE SUPPLIED WITH -input COMMAND

	in := os.Stdin

	// check for data being piped into stdin
	fi, _ := os.Stdin.Stat()
	isPipe := bool((fi.Mode() & os.ModeNamedPipe) != 0)

	usingFile := false

	if fileName != "" {

		inFile, err := os.Open(fileName)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
			os.Exit(1)
		}

		defer inFile.Close()

		// use indicated file instead of stdin
		in = inFile
		usingFile = true

		if isPipe && runtime.GOOS != "windows" {
			mode := fi.Mode().String()
			fmt.Fprintf(os.Stderr, "\nERROR: Input data from both stdin and file '%s', mode is '%s'\n", fileName, mode)
			os.Exit(1)
		}

	} else if runtime.GOOS != "windows" {

		fromStdin := bool((fi.Mode() & os.ModeCharDevice) == 0)
		if !isPipe || !fromStdin {
			mode := fi.Mode().String()
			fmt.Fprintf(os.Stderr, "\nERROR: No data supplied to xtract from stdin or file, mode is '%s'\n", mode)
			os.Exit(1)
		}
	}

	// check for -input command after extraction arguments
	for _, str := range args {
		if str == "-input" {
			fmt.Fprintf(os.Stderr, "\nERROR: Misplaced -input command\n")
			os.Exit(1)
		}
	}

	// CREATE XML BLOCK READER FROM STDIN OR FILE

	rdr := NewXMLReader(in, doCompress, doCleanup)
	if rdr == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
		os.Exit(1)
	}

	// SEQUENCE RECORD EXTRACTION COMMAND GENERATOR

	// -insd simplifies extraction of INSDSeq qualifiers
	if args[0] == "-insd" || args[0] == "-insd-" {

		addDash := true
		// -insd- variant suppresses use of dash as placeholder for missing qualifiers (undocumented)
		if args[0] == "-insd-" {
			addDash = false
		}

		args = args[1:]

		insd := ProcessINSD(args, isPipe || usingFile || testCount > 0, addDash)

		if !isPipe && !usingFile && testCount < 1 {
			// no piped input, so write output instructions
			fmt.Printf("xtract")
			for _, str := range insd {
				fmt.Printf(" %s", str)
			}
			fmt.Printf("\n")
			return
		}

		// data in pipe, so replace arguments, execute dynamically
		args = insd
	}

	// CITATION MATCHER EXTRACTION COMMAND GENERATOR

	// -hydra filters HydraResponse output by relevance score (undocumented)
	if args[0] == "-hydra" {

		hydra := ProcessHydra(isPipe || usingFile)

		if !isPipe && !usingFile {
			// no piped input, so write output instructions
			fmt.Printf("xtract")
			for _, str := range hydra {
				fmt.Printf(" %s", str)
			}
			fmt.Printf("\n")
			return
		}

		// data in pipe, so replace arguments, execute dynamically
		args = hydra
	}

	// XML DATA FORMATTING/COMPRESSION COMMAND GENERATOR

	// -reformat takes a parent pattern and compresses each object for fastest processing (undocumented)
	if args[0] == "-reformat" {

		args = args[1:]

		max := len(args)
		if max < 1 {
			fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -reformat command\n")
			os.Exit(1)
		}

		// required first argument is parent pattern, will explore using Parent/* construct, write component with -element "*"
		prnt := args[0]
		if prnt == "" {
			fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -reformat command\n")
			os.Exit(1)
		}
		if prnt == "-xml" || prnt == "-doctype" || prnt == "-pfx" || prnt == "-sfx" {
			fmt.Fprintf(os.Stderr, "\nERROR: Deprecated argument '%s' used in -reformat command\n", prnt)
			os.Exit(1)
		}

		// optional second argument controls XML expansion or compression level
		// asterisks MUST be quoted to avoid interpration as file wildcard by Unix shell
		elm := "*"
		addRet := false
		hideDoctype := false
		hideWrapper := false
		if max > 1 {
			// * = compact, ** = flush, *** = indented, **** = subtree, ***** = attributes on separate lines
			// @ = remove attributes
			// ^ = suppress xml and doctype, ^^ = also suppress parent set wrapper
			numStars := 0
			numCarets := 0
			hideAttrs := false
			for _, ch := range args[1] {
				if ch == '*' {
					numStars++
				} else if ch == '@' {
					hideAttrs = true
				} else if ch == '^' {
					numCarets++
				}
			}
			if numStars > 1 {
				addRet = true
			}
			if numCarets > 0 {
				hideDoctype = true
				if numCarets > 1 {
					hideWrapper = true
				}
			}
			// construct legal element argument for PrintSubtree
			switch numStars {
			case 1:
				elm = "*"
			case 2:
				elm = "**"
			case 3:
				elm = "***"
			case 4:
				elm = "****"
			case 5:
				elm = "*****"
			default:
				elm = "*"
			}
			if hideAttrs {
				elm += "@"
			}
		}

		// optional third argument provides detailed DOCTYPE construct
		doctype := ""
		if max > 2 {
			str := ConvertSlash(args[2])
			if strings.HasPrefix(str, "<!DOCTYPE ") && strings.HasSuffix(str, ">") {
				doctype = str
			}
		}

		if !isPipe && !usingFile {
			// no piped input, so write output instructions (without -head and -tail arguments)
			if addRet {
				fmt.Printf("xtract -pattern %s/* -ret \"\" -element \"%s\"\n", prnt, elm)
			} else {
				fmt.Printf("xtract -pattern %s/* -element \"%s\"\n", prnt, elm)
			}
			return
		}

		// add xml, DOCTYPE, and <Parent> lines at the beginning
		hd := fmt.Sprintf("<?xml version=\"1.0\"?>\n<!DOCTYPE %s>\n<%s>", prnt, prnt)
		if doctype != "" {
			// use supplied DOCTYPE argument
			hd = fmt.Sprintf("<?xml version=\"1.0\"?>\n%s\n<%s>", doctype, prnt)
		}
		if hideDoctype {
			// or just <Parent> line
			hd = fmt.Sprintf("<%s>", prnt)
		}
		// add </Parent> line at the end
		tl := fmt.Sprintf("</%s>", prnt)

		// use -pattern Parent/* construct
		prnt += "/*"

		var acc []string

		if !hideWrapper {
			acc = append(acc, "-head", hd, "-tail", tl)
		}
		acc = append(acc, "-pattern", prnt)
		if addRet {
			acc = append(acc, "-ret", "")
		}
		acc = append(acc, "-element", elm)

		// data in pipe, so replace arguments, execute dynamically
		args = acc
	}

	// CONFIRM INPUT DATA AVAILABILITY AFTER RUNNING COMMAND GENERATORS

	if testCount < 1 && !usingFile && !isPipe {

		fmt.Fprintf(os.Stderr, "\nERROR: No XML input data supplied to xtract\n")
		os.Exit(1)
	}

	// START PROFILING IF REQUESTED

	if prfl {

		dbug = true

		f, err := os.Create("cpu.pprof")
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create profile output file\n")
			os.Exit(1)
		}

		pprof.StartCPUProfile(f)

		defer pprof.StopCPUProfile()
	}

	// INITIALIZE PROCESS TIMER AND RECORD COUNT

	startTime := time.Now()
	recordCount := 0
	byteCount := 0

	// function to print processing rate and program duration
	printDuration := func(name string) {
		stopTime := time.Now()
		duration := stopTime.Sub(startTime)
		seconds := float64(duration.Nanoseconds()) / 1e9

		if recordCount >= 1000000 {
			fmt.Fprintf(os.Stderr, "\nXtract processed %d million %s in %.3f seconds", recordCount/1000000, name, seconds)
		} else {
			fmt.Fprintf(os.Stderr, "\nXtract processed %d %s in %.3f seconds", recordCount, name, seconds)
		}
		if seconds >= 0.001 && recordCount > 0 {
			rate := int(float64(recordCount) / seconds)
			if rate >= 1000000 {
				fmt.Fprintf(os.Stderr, " (%d mega%s/second", rate/1000000, name)
			} else {
				fmt.Fprintf(os.Stderr, " (%d %s/second", rate, name)
			}
			if byteCount > 0 {
				rate := int(float64(byteCount) / seconds)
				if rate >= 1000000 {
					fmt.Fprintf(os.Stderr, ", %d megabytes/second", rate/1000000)
				} else if rate >= 1000 {
					fmt.Fprintf(os.Stderr, ", %d kilobytes/second", rate/1000)
				} else {
					fmt.Fprintf(os.Stderr, ", %d bytes/second", rate)
				}
			}
			fmt.Fprintf(os.Stderr, ")")
		}
		fmt.Fprintf(os.Stderr, "\n\n")
	}

	// PERFORMANCE TIMING COMMANDS

	inSwitch = true
	action := NOPROCESS
	recordType := ""

	switch args[0] {
	case "-chunk":
		action = DOCHUNK
		recordType = "blocks"
	case "-split":
		action = DOSPLIT
		recordType = "patterns"
	case "-drain":
		action = DODRAIN
		recordType = "patterns"
	case "-token":
		action = DOTOKEN
		recordType = "tokens"
	default:
		// if not any of the formatting commands, keep going
		inSwitch = false
	}

	if inSwitch {
		recordCount, byteCount = ProcessXMLStream(rdr, tbls, args, action)
		printDuration(recordType)
		return
	}

	// SPECIAL FORMATTING COMMANDS

	inSwitch = true
	action = NOPROCESS

	switch args[0] {
	case "-format":
		action = DOFORMAT
	case "-outline":
		action = DOOUTLINE
	case "-synopsis":
		action = DOSYNOPSIS
	case "-verify", "-validate":
		action = DOVERIFY
	case "-filter":
		action = DOFILTER
	default:
		// if not any of the formatting commands, keep going
		inSwitch = false
	}

	if inSwitch {
		ProcessXMLStream(rdr, tbls, args, action)
		return
	}

	// SPECIFY STRINGS TO GO BEFORE AND AFTER ENTIRE OUTPUT OR EACH RECORD

	head := ""
	tail := ""

	hd := ""
	tl := ""

	for {

		inSwitch = true

		switch args[0] {
		case "-head":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -head command\n")
				os.Exit(1)
			}
			head = ConvertSlash(args[1])
		case "-tail":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -tail command\n")
				os.Exit(1)
			}
			tail = ConvertSlash(args[1])
		case "-hd":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -hd command\n")
				os.Exit(1)
			}
			hd = ConvertSlash(args[1])
		case "-tl":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Pattern missing after -tl command\n")
				os.Exit(1)
			}
			tl = ConvertSlash(args[1])
		default:
			// if not any of the controls, set flag to break out of for loop
			inSwitch = false
		}

		if !inSwitch {
			break
		}

		// skip past arguments
		args = args[2:]

		if len(args) < 1 {
			fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract\n")
			os.Exit(1)
		}
	}

	// ENSURE PRESENCE OF PATTERN ARGUMENT

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to xtract\n")
		os.Exit(1)
	}

	// make sure top-level -pattern command is next
	if args[0] != "-pattern" && args[0] != "-Pattern" {
		fmt.Fprintf(os.Stderr, "\nERROR: No -pattern in command-line arguments\n")
		os.Exit(1)
	}
	if len(args) < 2 {
		fmt.Fprintf(os.Stderr, "\nERROR: Item missing after -pattern command\n")
		os.Exit(1)
	}

	topPat := args[1]
	if topPat == "" {
		fmt.Fprintf(os.Stderr, "\nERROR: Item missing after -pattern command\n")
		os.Exit(1)
	}
	if strings.HasPrefix(topPat, "-") {
		fmt.Fprintf(os.Stderr, "\nERROR: Misplaced %s command\n", topPat)
		os.Exit(1)
	}

	// look for -pattern Parent/* construct for heterogeneous data, e.g., -pattern PubmedArticleSet/*
	topPattern, star := SplitInTwoAt(topPat, "/", LEFT)
	if topPattern == "" {
		return
	}

	parent := ""
	if star == "*" {
		parent = topPattern
	} else if star != "" {
		fmt.Fprintf(os.Stderr, "\nERROR: -pattern Parent/Child construct is not supported\n")
		os.Exit(1)
	}

	// PARSE AND VALIDATE EXTRACTION ARGUMENTS

	// parse nested exploration instruction from command-line arguments
	cmds := ParseArguments(args, topPattern)
	if cmds == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Problem parsing command-line arguments\n")
		os.Exit(1)
	}

	// PERFORMANCE TIMING COMMAND

	// -stats with an extraction command prints XML size and processing time for each record
	if stts {

		legend := "REC\tSIZE\tTIME"

		PartitionPattern(topPattern, star, rdr,
			func(rec int, str string) {
				beginTime := time.Now()
				ProcessQuery(str[:], parent, "", "", rec, cmds, tbls)
				endTime := time.Now()
				duration := endTime.Sub(beginTime)
				micro := int(float64(duration.Nanoseconds()) / 1e3)
				if legend != "" {
					fmt.Printf("%s\n", legend)
					legend = ""
				}
				fmt.Printf("%d\t%d\t%d\n", rec, len(str), micro)
			})

		return
	}

	// PERFORMANCE OPTIMIZATION FUNCTIONS

	// -test N runs a test extraction N times for each -proc from 1 to nCPU (undocumented)
	if testCount > 0 && testString != "" {

		// clean up copy of sample string included in source code
		sample := strings.TrimSpace(testString)
		sample = CleanupBadSpaces(sample)

		legend := "CPU\tTIME\tRATE"

		for numServ := 1; numServ <= ncpu; numServ++ {

			runtime.GOMAXPROCS(numServ)

			// alternative producer sends sample XML through channel N times
			xmlq := make(chan Extract, chanDepth)
			go func(out chan<- Extract) {
				for rec := 1; rec <= testCount; rec++ {
					out <- Extract{rec, sample}
				}
				close(out)
			}(xmlq)
			tblq := CreateConsumers(cmds, tbls, parent, "", "", numServ, xmlq)

			if xmlq == nil || tblq == nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create servers\n")
				os.Exit(1)
			}

			begTime := time.Now()
			recordCount = 0

			for _ = range tblq {
				recordCount++
				runtime.Gosched()
			}

			debug.FreeOSMemory()

			endTime := time.Now()
			expended := endTime.Sub(begTime)
			secs := float64(expended.Nanoseconds()) / 1e9

			if secs >= 0.000001 && recordCount > 0 {
				speed := int(float64(recordCount) / secs)
				if legend != "" {
					fmt.Printf("%s\n", legend)
					legend = ""
				}
				fmt.Printf("%d\t%.3f\t%d\n", numServ, secs, speed)
			}
		}

		return
	}

	// -trial -input fileName runs the specified extraction for each -proc from 1 to nCPU
	if trial && fileName != "" {

		legend := "CPU\tRATE\tDEV"

		for numServ := 1; numServ <= ncpu; numServ++ {

			runtime.GOMAXPROCS(numServ)

			sum := 0
			count := 0
			mean := 0.0
			m2 := 0.0

			// calculate mean and standard deviation of processing rate
			for trials := 0; trials < 5; trials++ {

				inFile, err := os.Open(fileName)
				if err != nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
					os.Exit(1)
				}

				rdr := NewXMLReader(inFile, doCompress, doCleanup)
				if rdr == nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to read input file\n")
					os.Exit(1)
				}

				xmlq := CreateProducer(topPattern, star, rdr, tbls)
				tblq := CreateConsumers(cmds, tbls, parent, hd, tl, numServ, xmlq)

				if xmlq == nil || tblq == nil {
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to create servers\n")
					os.Exit(1)
				}

				begTime := time.Now()
				recordCount = 0

				for _ = range tblq {
					recordCount++
					runtime.Gosched()
				}

				inFile.Close()

				debug.FreeOSMemory()

				endTime := time.Now()
				expended := endTime.Sub(begTime)
				secs := float64(expended.Nanoseconds()) / 1e9

				if secs >= 0.000001 && recordCount > 0 {
					speed := int(float64(recordCount) / secs)
					sum += speed
					count++
					x := float64(speed)
					delta := x - mean
					mean += delta / float64(count)
					m2 += delta * (x - mean)
				}
			}

			if legend != "" {
				fmt.Printf("%s\n", legend)
				legend = ""
			}
			if count > 1 {
				vrc := m2 / float64(count-1)
				dev := int(math.Sqrt(vrc))
				fmt.Printf("%d\t%d\t%d\n", numServ, sum/count, dev)
			}
		}

		return
	}

	// PROCESS SINGLE SELECTED RECORD IF -pattern ARGUMENT IS IMMEDIATELY FOLLOWED BY -position COMMAND

	if cmds.Visit == topPat && cmds.Position != "" {

		qry := ""
		idx := 0

		if cmds.Position == "first" {

			PartitionPattern(topPattern, star, rdr,
				func(rec int, str string) {
					if rec == 1 {
						qry = str
						idx = rec
					}
				})

		} else if cmds.Position == "last" {

			PartitionPattern(topPattern, star, rdr,
				func(rec int, str string) {
					qry = str
					idx = rec
				})

		} else {

			// use numeric position
			number, err := strconv.Atoi(cmds.Position)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized position '%s'\n", cmds.Position)
				os.Exit(1)
			}

			PartitionPattern(topPattern, star, rdr,
				func(rec int, str string) {
					if rec == number {
						qry = str
						idx = rec
					}
				})
		}

		if qry == "" {
			return
		}

		// clear position on top node to prevent condition test failure
		cmds.Position = ""

		// process single selected record
		res := ProcessQuery(qry[:], parent, "", "", idx, cmds, tbls)
		if res != "" {
			fmt.Printf("%s\n", res)
		}

		return
	}

	// LAUNCH PRODUCER AND CONSUMER SERVERS

	// launch producer goroutine to partition XML by pattern
	xmlq := CreateProducer(topPattern, star, rdr, tbls)
	if xmlq == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create producer\n")
		os.Exit(1)
	}

	// launch consumer goroutines to parse and explore partitioned XML objects
	tblq := CreateConsumers(cmds, tbls, parent, hd, tl, numServers, xmlq)
	if tblq == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create consumers\n")
		os.Exit(1)
	}

	// PERFORMANCE SUMMARY

	if dbug {

		// drain results, but suppress extraction output
		for ext := range tblq {
			byteCount += len(ext.Text)
			recordCount++
			runtime.Gosched()
		}

		// force garbage collection, return memory to operating system
		debug.FreeOSMemory()

		// print processing parameters as XML object
		stopTime := time.Now()
		duration := stopTime.Sub(startTime)
		seconds := float64(duration.Nanoseconds()) / 1e9

		// Threads is a more easily explained concept than GOMAXPROCS
		fmt.Printf("<Xtract>\n")
		fmt.Printf("  <Threads>%d</Threads>\n", numProcs)
		fmt.Printf("  <Parsers>%d</Parsers>\n", numServers)
		fmt.Printf("  <Time>%.3f</Time>\n", seconds)
		if seconds >= 0.001 && recordCount > 0 {
			rate := int(float64(recordCount) / seconds)
			fmt.Printf("  <Rate>%d</Rate>\n", rate)
		}
		fmt.Printf("</Xtract>\n")

		return
	}

	// DRAIN OUTPUT CHANNEL TO EXECUTE EXTRACTION COMMANDS, RESTORE OUTPUT ORDER WITH HEAP

	// initialize empty heap
	hp := &ExtractHeap{}
	heap.Init(hp)

	// index of next desired result
	next := 1

	delay := 0

	var buffer bytes.Buffer
	count := 0
	okay := false

	if head != "" {
		buffer.WriteString(head[:])
		buffer.WriteString("\n")
	}

	// printResult prints output for current pattern, handles -empty and -index flags, and periodically flushes buffer
	printResult := func(curr Extract) {

		str := curr.Text

		if mpty {

			if str == "" {

				okay = true

				idx := curr.Index
				val := strconv.Itoa(idx)
				buffer.WriteString(val[:])
				buffer.WriteString("\n")

				count++
			}

		} else if str != "" {

			okay = true

			if indx {
				idx := curr.Index
				val := strconv.Itoa(idx)
				buffer.WriteString(val[:])
				buffer.WriteString("\t")
			}

			// save output to byte buffer
			buffer.WriteString(str[:])

			count++
		}

		if count > 1000 {
			count = 0
			txt := buffer.String()
			if txt != "" {
				// print current buffer
				os.Stdout.WriteString(txt[:])
			}
			buffer.Reset()
		}
	}

	for ext := range tblq {

		// push result onto heap
		heap.Push(hp, ext)

		// read several values before checking to see if next record to print has been processed
		if delay < heapSize {
			delay++
			continue
		}

		delay = 0

		for hp.Len() > 0 {
			// remove lowest item from heap, use interface type assertion
			curr := heap.Pop(hp).(Extract)

			if curr.Index == next {

				// if this is the desired item, send to output
				printResult(curr)

				recordCount++

				// increment index
				next++
				// and keep checking heap to see if next result is already available
			} else {
				// otherwise push back onto heap
				heap.Push(hp, curr)
				// and go back to waiting on input channel
				break
			}
		}
	}

	// send remainder of heap to output
	for hp.Len() > 0 {
		curr := heap.Pop(hp).(Extract)

		printResult(curr)

		recordCount++
	}

	if tail != "" {
		buffer.WriteString(tail[:])
		buffer.WriteString("\n")
	}

	// do not print head or tail if no extraction output
	if okay {
		txt := buffer.String()
		if txt != "" {
			// print final buffer
			os.Stdout.WriteString(txt[:])
		}
	}
	buffer.Reset()

	// force garbage collection and return memory before calculating processing rate
	debug.FreeOSMemory()

	if timr {
		printDuration("records")
	}
}
