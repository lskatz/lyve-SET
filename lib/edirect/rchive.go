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
// File Name:  rchive.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"container/heap"
	"encoding/binary"
	"fmt"
	"github.com/fiam/gounidecode/unidecode"
	"hash/crc32"
	"html"
	"io"
	"io/ioutil"
	"os"
	"os/user"
	"path"
	"regexp"
	"runtime"
	"runtime/debug"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
	"unicode"
)

// RCHIVE VERSION AND HELP MESSAGE TEXT

const rchiveVersion = "8.60"

const rchiveHelp = `
Processing Flags

  -strict     Remove HTML and MathML tags
  -mixed      Allow PubMed mixed content

Data Source

  -input      Read XML from file instead of stdin

Local Record Cache

  -archive    Base path for saving individual XML files
  -index      Use [parent/element@attribute^version] for identifier

  -fetch      Base path for retrieving XML files
  -stream     Path for retrieving compressed XML

  -flag       [strict|mixed|none]
  -gzip       Use compression for local XML files
  -hash       Print UIDs and checksum values to stdout

  -trie       Print archive trie

Local Record Index

  -e2index    Create Entrez index XML (in xtract)
  -invert     Generate inverted index on specified field
  -merge      Combine inverted indices, divide by term prefix
  -deleted    File of deleted and versioned UIDs
  -promote    Create term lists and posting files

  -query      Search on individual words or overlapping word pairs

  -count      Print terms and counts, merging wildcards
  -counts     Expand wildcards, print individual term counts

  -stems      Apply Porter2 stemming to queries
  -stops      Retain stop words in queries

  -flag       [stems|none]

Documentation

  -help       Print this document
  -version    Print version number

Sample File Download

  ftp-cp ftp.ncbi.nlm.nih.gov /entrez/entrezdirect/samples carotene.xml.zip
  unzip carotene.xml.zip
  rm carotene.xml.zip

Mammalian Sequence Download

  download-sequence gbmam gbpri gbrod

Human Subset Extraction

  #!/bin/sh

  for fl in gbpri?.aso.gz gbpri??.aso.gz
  do
    run-ncbi-converter asn2all -i "$fl" -a t -b -c -O 9606 -f s > ${fl%.aso.gz}.xml
  done

Populate PubMed Archive

  export EDIRECT_PUBMED_MASTER=/Volumes/alexandria

  for dir in Archive Pubmed Indexed Inverted Merged Postings
  do
    mkdir -p "$MASTER/$dir"
  done

  pm-prepare "$MASTER/Archive"

  cd "$MASTER/Pubmed"

  download-pubmed baseline updatefiles

  pm-clean "$MASTER/Cleaned"

  cd "$MASTER/Cleaned"

  pm-stash "$MASTER/Archive"

  pm-refresh "$MASTER/Archive"

Retrieve from PubMed Archive

  cat subset.uid |
  fetch-pubmed "$MASTER/Archive" > subset.xml

Record Counts

  rchive -count "$MASTER/Postings" "protease inhibit* AND catabolite repress*"

Wildcard Expansion

  rchive -counts "$MASTER/Postings" "protease inhibit* AND catabolite repress*"

Query Processing

  rchive -query "$MASTER/Postings" "protease inhibit* AND catabolite repress*"

DISABLE ANTI-VIRUS FILE SCANNING FOR LOCAL ARCHIVES OR DESIGNATE AS TRUSTED FILES

DISABLE SPOTLIGHT INDEXING FOR EXTERNAL DISKS CONTAINING LOCAL ARCHIVES
`

const rchiveExtras = `
Maintenance Commands

  -prepare    [release|report] Compare daily update to archive
  -ignore     Ignore contents of object in -prepare comparisons
  -damaged    Report UIDs containing damaged embedded HTML tags
  -missing    Print list of missing identifiers
  -unique     File of UIDs for skipping all but last version

Miscellaneous

  -head       Print before everything else
  -tail       Print after everything else
  -hd         Print before each record
  -tl         Print after each record

Update Candidate Report

  cd "$MASTER/Pubmed"
  gunzip -c *.xml.gz | xtract -strict -compress -format flush |
  rchive -prepare report -ignore DateRevised -archive "$MASTER/Archive" \
    -index MedlineCitation/PMID -pattern PubmedArticle

Unnecessary Update Removal

  cd "$MASTER/Pubmed"
  gunzip -c *.xml.gz | xtract -strict -compress -format flush |
  rchive -prepare release -ignore DateRevised -archive "$MASTER/Archive" -index MedlineCitation/PMID \
    -head "<PubmedArticleSet>" -tail "</PubmedArticleSet>" -pattern PubmedArticle |
  xtract -format indent -xml '<?xml version="1.0" encoding="UTF-8"?>' \
    -doctype '<!DOCTYPE PubmedArticleSet SYSTEM "http://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_180101.dtd">' |
  gzip > newupdate.xml.gz

Get Archive UID List

  pm-uids "$MASTER/Archive" > complete.uid

Reconstruct List of Versioned PMIDs

  cd "$MASTER/Pubmed"
  rm -f "$MASTER/Archive/versioned.uid"
  gunzip -c *.xml.gz |
  xtract -strict -pattern PubmedArticle -if MedlineCitation/PMID@Version -gt 1 \
    -element MedlineCitation/PMID > "$MASTER/Archive/versioned.uid"

Reconstruct Release Files

  split -a 3 -l 30000 release.uid uids-
  n=1
  for x in uids-???
  do
    xmlfile=$(printf "pubmed18n%04d.xml.gz" "$n")
    n=$((n+1))
    echo "$xmlfile"
    cat "$x" |
    rchive -fetch "$MASTER/Archive" -head "<PubmedArticleSet>" -tail "</PubmedArticleSet>" |
    xtract -strict -format indent -xml '<?xml version="1.0" encoding="UTF-8"?>' \
      -doctype '<!DOCTYPE PubmedArticleSet SYSTEM "http://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_180101.dtd">' |
    gzip > "$xmlfile"
  done
  rm -rf uids-???

Damaged Embedded HTML Tag Search

  for fl in *.xml.gz
  do
    echo "$fl"
    gunzip -c "$fl" | rchive -mixed -damaged -index MedlineCitation/PMID^Version -pattern PubmedArticle
  done

  grep -v pubmed18n | grep AMPER | cut -f 1,6

Entrez Indexing

  cat carotene.xml | xtract -strict -stems -e2index PubmedArticle \
    MedlineCitation/PMID ArticleTitle,Abstract/AbstractText > carotene.e2x

Index Inversion

  rchive -invert PAIR carotene.e2x > carotene.inv

Merge Indices

  rchive -deleted "$MASTER/Archive/deleted.uid" -merge "$MASTER/Merged" PAIR carotene.inv

Create Postings

  rchive -promote "$MASTER/Postings" PAIR carotene.mrg

Reconstruct Term List Keys

  for fld in PAIR NORM STEM GRFT
  do
    rm -f "$MASTER/Postings/$fld/sections.txt"
    find "$MASTER/Postings/$fld" -name "*.mst" |
    sed -e 's,.*/\(.*\)\.mst,\1,' |
    sort | uniq > "$MASTER/Postings/$fld/sections.txt"
  done

Generate Term List Paths

  for fld in PAIR NORM STEM GRFT
  do
    find "$MASTER/Postings/$fld" -name "*.trm" |
    sed -e 's,\(.*/\)\(.*\.trm\),\1 \2,' |
    sort -k 2 | uniq | tr -d ' '
  done

DISABLE ANTI-VIRUS FILE SCANNING FOR LOCAL ARCHIVES OR DESIGNATE AS TRUSTED FILES

DISABLE SPOTLIGHT INDEXING FOR EXTERNAL DISKS CONTAINING LOCAL ARCHIVES
`

const rchiveInternal = `
Performance Default Overrides

  -proc     Number of CPU processors used
  -cons     Ratio of parsers to processors
  -serv     Concurrent parser instances
  -chan     Communication channel depth
  -heap     Order restoration heap size
  -farm     Node allocation buffer length
  -gogc     Garbage collection tuning knob

Debugging

  -timer    Report processing duration and rate

Execution Profiling

  ./rchive -profile -invert PAIR carotene.e2x  > /dev/null
  go tool pprof --pdf ./rchive ./cpu.pprof > ./callgraph.pdf
`

// DATA OBJECTS

type Master struct {
	TermOffset int32
	PostOffset int32
}

// UTILITIES

func ReportEncodedMarkup(typ, id, str string) {

	var buffer strings.Builder

	max := len(str)

	lookAhead := func(txt string, to int) string {

		mx := len(txt)
		if to > mx {
			to = mx
		}
		pos := strings.Index(txt[:to], "gt;")
		if pos > 0 {
			to = pos + 3
		}
		return txt[:to]
	}

	findContext := func(fr, to int) string {

		numSpaces := 0

		for fr > 0 {
			ch := str[fr]
			if ch == ' ' {
				numSpaces++
				if numSpaces > 1 {
					fr++
					break
				}
			} else if ch == '\n' || ch == '>' {
				fr++
				break
			}
			fr--
		}

		numSpaces = 0

		for to < max {
			ch := str[to]
			if ch == ' ' {
				numSpaces++
				if numSpaces > 1 {
					break
				}
			} else if ch == '\n' || ch == '<' {
				break
			}
			to++
		}

		return str[fr:to]
	}

	reportMarkup := func(lbl string, fr, to int, txt string) {

		if lbl == typ || typ == "ALL" {
			// extract XML of SELF, SINGLE, DOUBLE, or AMPER types, or ALL
			buffer.WriteString(str)
			buffer.WriteString("\n")
		} else if typ == "" {
			// print report
			buffer.WriteString(id)
			buffer.WriteString("\t")
			buffer.WriteString(lbl)
			buffer.WriteString("\t")
			buffer.WriteString(txt)
			buffer.WriteString("\t| ")
			ctx := findContext(fr, to)
			buffer.WriteString(ctx)
			if HasUnicodeMarkup(ctx) {
				ctx = RepairUnicodeMarkup(ctx, SPACE)
			}
			ctx = RepairEncodedMarkup(ctx)
			buffer.WriteString("\t| ")
			buffer.WriteString(ctx)
			if HasAmpOrNotASCII(ctx) {
				ctx = html.UnescapeString(ctx)
			}
			buffer.WriteString("\t| ")
			buffer.WriteString(ctx)
			buffer.WriteString("\n")
		}
	}

	/*
		badTags := [10]string{
			"<i/>",
			"<i />",
			"<b/>",
			"<b />",
			"<u/>",
			"<u />",
			"<sup/>",
			"<sup />",
			"<sub/>",
			"<sub />",
		}
	*/

	skip := 0

	/*
		var prev rune
	*/

	for i, ch := range str {
		if skip > 0 {
			skip--
			continue
		}
		/*
			if ch > 127 {
				if IsUnicodeSuper(ch) {
					if IsUnicodeSubsc(prev) {
						// reportMarkup("UNIUP", i, i+2, string(ch))
					}
				} else if IsUnicodeSubsc(ch) {
					if IsUnicodeSuper(prev) {
						// reportMarkup("UNIDN", i, i+2, string(ch))
					}
				} else if ch == '\u0038' || ch == '\u0039' {
					// reportMarkup("ANGLE", i, i+2, string(ch))
				}
				prev = ch
				continue
			} else {
				prev = ' '
			}
		*/
		if ch == '<' {
			/*
				j := i + 1
				if j < max {
					nxt := str[j]
					if nxt == 'i' || nxt == 'b' || nxt == 'u' || nxt == 's' {
						for _, tag := range badTags {
							if strings.HasPrefix(str, tag) {
								k := len(tag)
								reportMarkup("SELF", i, i+k, tag)
								break
							}
						}
					}
				}
				if strings.HasPrefix(str[i:], "</sup><sub>") {
					// reportMarkup("SUPSUB", i, i+11, "</sup><sub>")
				} else if strings.HasPrefix(str[i:], "</sub><sup>") {
					// reportMarkup("SUBSUP", i, i+11, "</sub><sup>")
				}
			*/
			continue
		} else if ch != '&' {
			continue
		} else if strings.HasPrefix(str[i:], "&lt;") {
			sub := lookAhead(str[i:], 14)
			_, ok := htmlRepair[sub]
			if ok {
				skip = len(sub) - 1
				reportMarkup("SINGLE", i, i+skip+1, sub)
				continue
			}
		} else if strings.HasPrefix(str[i:], "&amp;lt;") {
			sub := lookAhead(str[i:], 22)
			_, ok := htmlRepair[sub]
			if ok {
				skip = len(sub) - 1
				reportMarkup("DOUBLE", i, i+skip+1, sub)
				continue
			}
		} else if strings.HasPrefix(str[i:], "&amp;amp;") {
			reportMarkup("AMPER", i, i+9, "&amp;amp;")
			skip = 8
			continue
		}
	}

	res := buffer.String()

	os.Stdout.WriteString(res)
}

func SplitMarkedClause(str string) ([]string, []string) {

	if str == "" {
		return nil, nil
	}

	var norm []string
	var pair []string

	parts := strings.Split(str, "+")

	for _, segment := range parts {

		segment = strings.TrimSpace(segment)

		if segment == "" {
			continue
		}

		words := strings.Fields(segment)

		if len(words) < 2 {
			// index single term
			norm = append(norm, segment)
			continue
		}

		past := ""
		for _, item := range words {

			if past != "" {
				// index informative adjacent word pair
				adjacent := past + " " + item
				adjacent = CompressRunsOfSpaces(adjacent)
				adjacent = strings.TrimSpace(adjacent)
				pair = append(pair, adjacent)
			}
			past = item

		}
	}

	return norm, pair
}

// DIRECTORY PATH UTILITIES

// MakeArchiveTrie allows a short prefix of letters with an optional underscore, and splits the remainder into character pairs
func MakeArchiveTrie(str string, arry [132]rune) string {

	if len(str) > 64 {
		return ""
	}

	if IsAllDigitsOrPeriod(str) {

		// limit trie to first 6 characters
		if len(str) > 6 {
			str = str[:6]
		}
	}

	max := 4
	k := 0
	for _, ch := range str {
		if unicode.IsLetter(ch) {
			k++
			continue
		}
		if ch == '_' {
			k++
			max = 6
		}
		break
	}

	// prefix is up to three letters if followed by digits, or up to four letters if followed by an underscore
	pfx := str[:k]
	if len(pfx) < max {
		str = str[k:]
	} else {
		pfx = ""
	}

	i := 0

	if pfx != "" {
		for _, ch := range pfx {
			arry[i] = ch
			i++
		}
		arry[i] = '/'
		i++
	}

	between := 0
	doSlash := false

	// remainder is divided in character pairs, e.g., NP_/06/00/51 for NP_060051.2
	for _, ch := range str {
		// break at period separating accession from version
		if ch == '.' {
			break
		}
		if doSlash {
			arry[i] = '/'
			i++
			doSlash = false
		}
		if ch == ' ' {
			ch = '_'
		}
		if !unicode.IsLetter(ch) && !unicode.IsDigit(ch) {
			ch = '_'
		}
		arry[i] = ch
		i++
		between++
		if between > 1 {
			doSlash = true
			between = 0
		}
	}

	res := string(arry[:i])

	if !strings.HasSuffix(res, "/") {
		arry[i] = '/'
		i++
		res = string(arry[:i])
	}

	return strings.ToUpper(res)
}

// MakePostingsTrie splits a string into characters, separated by path delimiting slashes
func MakePostingsTrie(str string, arry [516]rune) string {

	if len(str) > 256 {
		return ""
	}

	// expand Greek letters, anglicize characters in other alphabets
	if IsNotASCII(str) {
		if HasGreek(str) {
			str = SpellGreek(str)
			str = CompressRunsOfSpaces(str)
		}
		str = unidecode.Unidecode(str)
		str = strings.TrimSpace(str)
	}

	i := 0
	doSlash := false

	for _, ch := range str {
		if doSlash {
			arry[i] = '/'
			i++
		}
		if ch == ' ' {
			ch = '_'
		}
		if !unicode.IsLetter(ch) && !unicode.IsDigit(ch) {
			ch = '_'
		}
		arry[i] = ch
		i++
		doSlash = true
	}

	return strings.ToLower(string(arry[:i]))
}

// POSTINGS FILE UTILITIES

// trieLen directory depth parameters are based on the observed size distribution of PubMed PAIR indices
var trieLen = map[string]int{
	"ac": 4,
	"af": 4,
	"an": 4,
	"ca": 4,
	"ce": 4,
	"cl": 4,
	"co": 4,
	"di": 4,
	"ex": 4,
	"ge": 4,
	"gr": 4,
	"he": 4,
	"hi": 4,
	"in": 4,
	"me": 4,
	"mo": 4,
	"no": 4,
	"pa": 4,
	"pe": 4,
	"pl": 4,
	"po": 4,
	"pr": 4,
	"re": 4,
	"si": 4,
	"sp": 4,
	"st": 4,
	"su": 4,
	"tr": 4,
	"tw": 4,
	"un": 3,
	"va": 3,
	"ve": 3,
	"vi": 3,
	"wh": 3,
}

func PostingDir(term string) string {

	if len(term) < 3 {
		return term
	}

	key := term[:2]

	num, ok := trieLen[key]
	if ok && len(term) >= num {
		return term[:num]
	}

	switch term[0] {
	case 'u', 'v', 'w', 'x', 'y', 'z':
		return term[:2]
	}

	return term[:3]
}

func IdentifierKey(term string) string {

	// remove punctuation from term
	key := strings.Map(func(c rune) rune {
		if !unicode.IsLetter(c) && !unicode.IsDigit(c) && c != ' ' && c != '-' && c != '_' {
			return -1
		}
		return c
	}, term)

	key = strings.Replace(key, " ", "_", -1)
	key = strings.Replace(key, "-", "_", -1)

	// use first 2, 3, or 4 characters of identifier for directory
	key = PostingDir(key)

	return key
}

func PostingPath(prom, field, term string, arry [516]rune) (string, string) {

	// use first few characters of identifier for directory
	dir := IdentifierKey(term)

	trie := MakePostingsTrie(dir, arry)
	if trie == "" {
		return "", ""
	}

	dpath := path.Join(prom, field, trie)

	return dpath, dir
}

func CommonOpenFile(dpath, fname string) (*os.File, int64) {

	fpath := path.Join(dpath, fname)
	if fpath == "" {
		return nil, 0
	}

	inFile, err := os.Open(fpath)
	if err != nil && os.IsNotExist(err) {
		return nil, 0
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil, 0
	}

	fi, err := inFile.Stat()
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil, 0
	}

	size := fi.Size()

	return inFile, size
}

func ReadMasterIndex(dpath, key string) []Master {

	inFile, size := CommonOpenFile(dpath, key+".mst")
	if inFile == nil {
		return nil
	}

	defer inFile.Close()

	data := make([]Master, size/8)
	if data == nil {
		return nil
	}

	err := binary.Read(inFile, binary.LittleEndian, &data)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	return data
}

func ReadTermList(dpath, key string) []byte {

	inFile, size := CommonOpenFile(dpath, key+".trm")
	if inFile == nil {
		return nil
	}

	defer inFile.Close()

	data := make([]byte, size)
	if data == nil {
		return nil
	}

	err := binary.Read(inFile, binary.LittleEndian, &data)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	return data
}

func ReadPostingData(dpath, key string, offset int32, size int32) []int32 {

	inFile, _ := CommonOpenFile(dpath, key+".pst")
	if inFile == nil {
		return nil
	}

	defer inFile.Close()

	data := make([]int32, size/4)
	if data == nil {
		return nil
	}

	_, err := inFile.Seek(int64(offset), io.SeekStart)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	err = binary.Read(inFile, binary.LittleEndian, data)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return nil
	}

	return data
}

func GetPostingIDs(prom, field, term string) []int32 {

	var arry [516]rune
	dpath, key := PostingPath(prom, field, term, arry)
	if dpath == "" {
		return nil
	}

	indx := ReadMasterIndex(dpath, key)
	if indx == nil {
		return nil
	}

	trms := ReadTermList(dpath, key)
	if trms == nil {
		return nil
	}

	strs := make([]string, len(indx)-1)
	if strs == nil {
		return nil
	}

	// populate array of strings from term list
	for i, j := 0, 1; i < len(indx)-1; i++ {
		from := indx[i].TermOffset
		to := indx[j].TermOffset - 1
		j++
		strs[i] = string(trms[from:to])
	}

	// change protecting underscore to space
	term = strings.Replace(term, "_", " ", -1)

	isWildCard := false
	if strings.HasSuffix(term, "*") {
		tlen := len(term)
		if tlen > 1 {
			isWildCard = true
			term = term[:tlen-1]
		}
		if tlen < 4 {
			fmt.Fprintf(os.Stderr, "Wildcard term '%s' must be at least 4 characters long - ignoring this word\n", term)
			return nil
		}
	}

	// binary search in term list
	num := len(strs)
	L, R := 0, num-1
	for L < R {
		mid := (L + R) / 2
		if strs[mid] < term {
			L = mid + 1
		} else {
			R = mid
		}
	}

	// wild card search scans term lists, fuses adjacent postings lists
	if isWildCard {
		if R < num && strings.HasPrefix(strs[R], term) {
			offset := indx[R].PostOffset
			for R < num && strings.HasPrefix(strs[R], term) {
				R++
			}
			size := indx[R].PostOffset - offset

			// term match, read relevant postings list section
			data := ReadPostingData(dpath, key, offset, size)
			if data == nil {
				return nil
			}

			merged := make(map[int32]bool)

			// combine all postings in term range
			for _, val := range data {
				merged[val] = true
			}

			fused := make([]int32, len(merged))

			// convert map to slice
			i := 0
			for num := range merged {
				fused[i] = num
				i++
			}

			sort.Slice(fused, func(i, j int) bool { return fused[i] < fused[j] })

			return fused
		}
		return nil
	}

	// regular search requires exact match from binary search
	if R < num && strs[R] == term {

		offset := indx[R].PostOffset
		size := indx[R+1].PostOffset - offset

		// term match, read relevant postings list section
		data := ReadPostingData(dpath, key, offset, size)
		if data == nil {
			return nil
		}
		return data
	}

	return nil
}

func PrintTermCounts(base, field, term string) int {

	if len(term) < 4 {
		fmt.Fprintf(os.Stderr, "\nERROR: Term count argument must be at least 4 characters\n")
		os.Exit(1)
	}

	if strings.Contains(term[:4], "*") {
		fmt.Fprintf(os.Stderr, "\nERROR: Wildcard asterisk must not be in first 4 characters\n")
		os.Exit(1)
	}

	var arry [516]rune
	dpath, key := PostingPath(base, field, term, arry)
	if dpath == "" {
		return 0
	}

	indx := ReadMasterIndex(dpath, key)
	if indx == nil {
		return 0
	}

	trms := ReadTermList(dpath, key)
	if trms == nil {
		return 0
	}

	strs := make([]string, len(indx)-1)
	if strs == nil {
		return 0
	}

	// populate array of strings from term list
	for i, j := 0, 1; i < len(indx)-1; i++ {
		from := indx[i].TermOffset
		to := indx[j].TermOffset - 1
		j++
		strs[i] = string(trms[from:to])
	}

	// change protecting underscore to space
	term = strings.Replace(term, "_", " ", -1)

	// flank pattern with start-of-string and end-of-string symbols
	pat := "^" + term + "$"

	// change asterisk in query to dot + star for regular expression
	pat = strings.Replace(pat, "*", ".*", -1)

	re, err := regexp.Compile(pat)

	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err.Error())
		return 0
	}

	count := 0

	for R, str := range strs {
		if re.MatchString(str) {
			offset := indx[R].PostOffset
			size := indx[R+1].PostOffset - offset
			fmt.Fprintf(os.Stdout, "%d\t%s\n", size/4, str)
			count++
		}
	}

	return count
}

// BOOLEAN OPERATIONS FOR POSTINGS LISTS

func IntersectIDs(N, M []int32) []int32 {

	if N == nil {
		return M
	}
	if M == nil {
		return N
	}

	n, m := len(N), len(M)

	// swap to make M the smaller list
	if n < m {
		N, M = M, N
		n, m = m, n
	}

	if m < 1 {
		return N
	}

	res := make([]int32, m)

	i, j, k := 0, 0, 0

	// use local variables for speed
	en := N[i]
	em := M[j]

	for {
		// do inequality tests first
		if en < em {
			// index to larger list most likely to be advanced
			i++
			if i == n {
				break
			}
			en = N[i]
		} else if en > em {
			j++
			if j == m {
				break
			}
			em = M[j]
		} else {
			// equality (intersection match) least likely
			res[k] = en
			k++
			i++
			j++
			if i == n || j == m {
				break
			}
			en = N[i]
			em = M[j]
		}
	}

	// truncate output array to actual size of intersection
	res = res[:k]

	return res
}

// if m * log(n) < m + n, binary search has fewer comparisons, but cache makes linear algorithm faster
/*
func IntersectBinary(N, M []int32) []int32 {

	if N == nil {
		return M
	}
	if M == nil {
		return N
	}

	n, m := len(N), len(M)

	// swap to make M the smaller list
	if n < m {
		N, M = M, N
		n, m = m, n
	}

	if m < 1 {
		return N
	}

	k := 0

	res := make([]int32, m)

	for _, uid := range M {
		// inline binary search is faster than sort.Search
		L, R := 0, n-1
		for L < R {
			mid := (L + R) / 2
			if N[mid] < uid {
				L = mid + 1
			} else {
				R = mid
			}
		}
		// R := sort.Search(len(N), func(i int) bool { return N[i] >= uid })
		if R < n && N[R] == uid {
			res[k] = uid
			k++
			// remove leading part of N for slight speed gain
			N = N[R:]
			n = len(N)
		}
	}

	res = res[:k]

	return res
}
*/

func CombineIDs(N, M []int32) []int32 {

	if N == nil {
		return M
	}
	if M == nil {
		return N
	}

	n, m := len(N), len(M)

	// swap to make M the smaller list
	if n < m {
		N, M = M, N
		n, m = m, n
	}

	if m < 1 {
		return N
	}

	i, j, k := 0, 0, 0

	res := make([]int32, n+m)

	for i < n && j < m {
		if N[i] < M[j] {
			res[k] = N[i]
			k++
			i++
		} else if N[i] > M[j] {
			res[k] = M[j]
			k++
			j++
		} else {
			res[k] = N[i]
			k++
			i++
			j++
		}
	}
	for i < n {
		res[k] = N[i]
		k++
		i++
	}
	for j < m {
		res[k] = M[j]
		k++
		j++
	}

	res = res[:k]

	return res
}

func ExcludeIDs(N, M []int32) []int32 {

	if N == nil {
		return nil
	}
	if M == nil {
		return N
	}

	combo := make(map[int32]bool)

	for _, uid := range N {
		combo[uid] = true
	}

	// delete postings that appear in both lists
	for _, uid := range M {
		if combo[uid] {
			delete(combo, uid)
		}
	}

	fused := make([]int32, len(combo))

	// convert map to slice
	i := 0
	for num := range combo {
		fused[i] = num
		i++
	}

	// need to sort results
	sort.Slice(fused, func(i, j int) bool { return fused[i] < fused[j] })

	return fused
}

// SEARCH TERM LISTS FOR PHRASES OR NORMALIZED TERMS, OR MATCH BY PATTERN

type Arrays struct {
	Data []int32
}

func ConvertTildes(phrase string) string {

	phrase = strings.Replace(phrase, "~ ~", "~~", -1)
	phrase = strings.Replace(phrase, "~ ~", "~~", -1)
	phrase = strings.Replace(phrase, "~~", "~", -1)
	phrase = strings.Replace(phrase, "~~", "~", -1)
	phrase = strings.Replace(phrase, "~", " AND ", -1)
	phrase = CompressRunsOfSpaces(phrase)

	return phrase
}

func ProcessQuery(base, phrase string) int {

	if phrase == "" {
		return 0
	}

	phrase = ConvertTildes(phrase)

	phrase = PrepareQuery(phrase)

	clauses := PartitionQuery(phrase)

	// to do wildcard on first term of pair, use "term **" to circumvent single character exclusion
	clauses = MarkClauses(clauses)

	if clauses == nil {
		return 0
	}

	nfield := "NORM"
	pfield := "PAIR"
	if DoStem {
		nfield = "STEM"
		pfield = "GRFT"
	}

	count := 0

	// check each phrase against record
	testPhrase := func(tokens []string) []int32 {

		eval := func(str string) []int32 {

			norm, pair := SplitMarkedClause(str)

			if norm == nil && pair == nil {
				return nil
			}

			var intersect []Arrays

			sort.Slice(pair, func(i, j int) bool { return pair[i] < pair[j] })

			prev := ""
			for _, term := range pair {

				if prev == term {
					continue
				}
				prev = term

				data := GetPostingIDs(base, pfield, term)
				if data == nil {
					// if a pair is not indexed, add individual words
					// if interested in closest match, it will have better selectivity
					// if using subsequent phrase confirmation, it will have fewer false positives to check
					words := strings.Fields(term)
					for _, item := range words {
						norm = append(norm, item)
					}
					// skip pairs that are not present in the index
					continue
				}

				one := Arrays{Data: data}
				intersect = append(intersect, one)
			}

			sort.Slice(norm, func(i, j int) bool { return norm[i] < norm[j] })

			prev = ""
			for _, term := range norm {

				if prev == term {
					continue
				}
				prev = term

				data := GetPostingIDs(base, nfield, term)
				if data == nil {
					// skip words that are not present in the index
					continue
				}

				one := Arrays{Data: data}
				intersect = append(intersect, one)
			}

			if len(intersect) < 1 {
				return nil
			}

			result := intersect[0].Data

			for i := 1; i < len(intersect); i++ {

				curr := intersect[i].Data

				result = IntersectIDs(result, curr)
			}

			count += len(intersect)

			return result
		}

		nextToken := func() string {

			if len(tokens) < 1 {
				return ""
			}

			tkn := tokens[0]
			tokens = tokens[1:]

			return tkn
		}

		// recursive definitions
		var excl func() ([]int32, string)
		var expr func() ([]int32, string)
		var fact func() ([]int32, string)
		var term func() ([]int32, string)

		fact = func() ([]int32, string) {

			var (
				data []int32
				tkn  string
			)

			tkn = nextToken()
			if tkn == "(" {
				data, tkn = expr()
				if tkn == ")" {
					tkn = nextToken()
				}
			} else {
				data = eval(tkn)
				tkn = nextToken()
			}

			return data, tkn
		}

		excl = func() ([]int32, string) {

			var next []int32

			data, tkn := fact()
			for tkn == "!" {
				next, tkn = fact()
				data = ExcludeIDs(data, next)
			}

			return data, tkn
		}

		term = func() ([]int32, string) {

			var next []int32

			data, tkn := excl()
			for tkn == "&" {
				next, tkn = excl()
				data = IntersectIDs(data, next)
			}

			return data, tkn
		}

		expr = func() ([]int32, string) {

			var next []int32

			data, tkn := term()
			for tkn == "|" {
				next, tkn = term()
				data = CombineIDs(data, next)
			}

			return data, tkn
		}

		// enter recursive descent parser
		result, _ := expr()

		return result
	}

	result := testPhrase(clauses)

	// sort final result
	sort.Slice(result, func(i, j int) bool { return result[i] < result[j] })

	for _, pmid := range result {
		fmt.Fprintf(os.Stdout, "%d\n", pmid)
	}

	return count
}

func ProcessCount(base, phrase string) int {

	if phrase == "" {
		return 0
	}

	phrase = ConvertTildes(phrase)

	phrase = PrepareQuery(phrase)

	clauses := PartitionQuery(phrase)

	// to do wildcard on first term of pair, use "term **" to circumvent single character exclusion
	clauses = MarkClauses(clauses)

	if clauses == nil {
		return 0
	}

	nfield := "NORM"
	pfield := "PAIR"
	if DoStem {
		nfield = "STEM"
		pfield = "GRFT"
	}

	count := 0

	checkTermCounts := func(str string) {

		norm, pair := SplitMarkedClause(str)

		if norm == nil && pair == nil {
			return
		}

		for _, term := range pair {

			data := GetPostingIDs(base, pfield, term)
			if data == nil {
				// if a pair is not indexed, add individual words
				words := strings.Fields(term)
				for _, item := range words {
					norm = append(norm, item)
				}
				// skip pairs that are not present in the index
				continue
			}
			size := len(data)
			fmt.Fprintf(os.Stdout, "%d\t%s\n", size, term)
			count += size
		}

		for _, term := range norm {

			data := GetPostingIDs(base, nfield, term)
			size := len(data)
			fmt.Fprintf(os.Stdout, "%d\t%s\n", size, term)
			count += size
		}
	}

	for _, item := range clauses {

		// skip control symbols
		if item == "(" || item == ")" || item == "&" || item == "|" || item == "!" {
			continue
		}

		checkTermCounts(item)
	}

	return count
}

func ProcessCounts(base, phrase string) int {

	if phrase == "" {
		return 0
	}

	phrase = ConvertTildes(phrase)

	phrase = PrepareQuery(phrase)

	clauses := PartitionQuery(phrase)

	// to do wildcard on first term of pair, use "term **" to circumvent single character exclusion
	clauses = MarkClauses(clauses)

	if clauses == nil {
		return 0
	}

	nfield := "NORM"
	pfield := "PAIR"
	if DoStem {
		nfield = "STEM"
		pfield = "GRFT"
	}

	count := 0

	checkTermCounts := func(str string) {

		norm, pair := SplitMarkedClause(str)

		if norm == nil && pair == nil {
			return
		}

		for _, term := range pair {

			count += PrintTermCounts(base, pfield, term)
		}

		for _, term := range norm {

			count += PrintTermCounts(base, nfield, term)
		}
	}

	for _, item := range clauses {

		// skip control symbols
		if item == "(" || item == ")" || item == "&" || item == "|" || item == "!" {
			continue
		}

		checkTermCounts(item)
	}

	return count
}

// READ FILE IN ARCHIVE THAT STORES UIDS WITH NON-DEFAULT VERSIONS

// CONCURRENT GOROUTINE SERVERS

// processes with single goroutine call defer close(out) so consumer(s) can range over channel
// processes with multiple instances call defer wg.Done(), separate goroutine uses wg.Wait() to delay close(out)

func CreateUIDReader(in io.Reader) <-chan Extract {

	if in == nil {
		return nil
	}

	out := make(chan Extract, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create uid reader channel\n")
		os.Exit(1)
	}

	// uidReader reads uids from input stream and sends through channel
	uidReader := func(in io.Reader, out chan<- Extract) {

		// close channel when all records have been processed
		defer close(out)

		scanr := bufio.NewScanner(in)

		idx := 0
		for scanr.Scan() {

			// read lines of identifiers
			file := scanr.Text()
			idx++

			pos := strings.Index(file, ".")
			if pos >= 0 {
				// remove version suffix
				file = file[:pos]
			}

			out <- Extract{idx, "", file, nil}
		}
	}

	// launch single uid reader goroutine
	go uidReader(in, out)

	return out
}

func CreateStashers(stash, parent, indx string, hash, zipp bool, inp <-chan Extract) <-chan string {

	if inp == nil {
		return nil
	}

	out := make(chan string, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create stasher channel\n")
		os.Exit(1)
	}

	find := ParseIndex(indx)

	sfx := ".xml"
	if zipp {
		sfx = ".xml.gz"
	}

	type StasherType int

	const (
		OKAY StasherType = iota
		WAIT
		BAIL
	)

	// mutex to protect access to inUse map
	var wlock sync.Mutex

	// map to track files currently being written
	inUse := make(map[string]int)

	// lockFile function prevents colliding writes
	lockFile := func(id string, index int) StasherType {

		// map is non-reentrant, protect with mutex
		wlock.Lock()

		// multiple return paths, schedule the unlock command up front
		defer wlock.Unlock()

		idx, ok := inUse[id]

		if ok {
			if idx < index {
				// later version is being written by another goroutine, skip this
				return BAIL
			}
			// earlier version is being written by another goroutine, wait
			return WAIT
		}

		// okay to write file, mark in use to prevent collision
		inUse[id] = index
		return OKAY
	}

	// freeFile function removes entry from inUse map
	freeFile := func(id string) {

		wlock.Lock()

		defer wlock.Unlock()

		// free entry in map, later versions of same record can now be written
		delete(inUse, id)
	}

	// stashRecord saves individual XML record to archive file accessed by trie
	stashRecord := func(str, id string, index int) string {

		pos := strings.Index(id, ".")
		if pos >= 0 {
			// remove version from UID
			id = id[:pos]
		}

		var arry [132]rune
		trie := MakeArchiveTrie(id, arry)
		if trie == "" {
			return ""
		}

		attempts := 5
		keepChecking := true

		for keepChecking {
			// check if file is not being written by another goroutine
			switch lockFile(id, index) {
			case OKAY:
				// okay to save this record now
				keepChecking = false
			case WAIT:
				// earlier version is being saved, wait one second and try again
				time.Sleep(time.Second)
				attempts--
				if attempts < 1 {
					// cannot get lock after several attempts
					fmt.Fprintf(os.Stderr, "\nERROR: Unable to save '%s'\n", id)
					return ""
				}
			case BAIL:
				// later version is being saved, skip this one
				return ""
			default:
			}
		}

		// delete lock after writing file
		defer freeFile(id)

		dpath := path.Join(stash, trie)
		if dpath == "" {
			return ""
		}
		_, err := os.Stat(dpath)
		if err != nil && os.IsNotExist(err) {
			err = os.MkdirAll(dpath, os.ModePerm)
		}
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return ""
		}
		fpath := path.Join(dpath, id+sfx)
		if fpath == "" {
			return ""
		}

		// overwrites and truncates existing file
		fl, err := os.Create(fpath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return ""
		}

		defer fl.Close()

		defer fl.Sync()

		res := ""

		if hash {
			// calculate hash code for verification table
			hsh := crc32.NewIEEE()
			hsh.Write([]byte(str))
			val := hsh.Sum32()
			res = strconv.FormatUint(uint64(val), 10)
		}

		if zipp {

			zpr, err := gzip.NewWriterLevel(fl, gzip.DefaultCompression)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
				return ""
			}

			defer zpr.Close()

			wrtr := bufio.NewWriter(zpr)

			defer wrtr.Flush()

			// compress and copy record to file
			wrtr.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				wrtr.WriteString("\n")
			}

		} else {

			// copy record to file
			fl.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				fl.WriteString("\n")
			}
		}

		return res
	}

	// xmlStasher reads from channel and calls stashRecord
	xmlStasher := func(wg *sync.WaitGroup, inp <-chan Extract, out chan<- string) {

		defer wg.Done()

		for ext := range inp {

			ext.Ident = FindIdentifier(ext.Text, parent, find)

			hsh := stashRecord(ext.Text, ext.Ident, ext.Index)

			res := ext.Ident
			if hash {
				res += "\t" + hsh
			}
			res += "\n"

			out <- res
		}
	}

	var wg sync.WaitGroup

	// launch multiple stasher goroutines
	for i := 0; i < NumServe; i++ {
		wg.Add(1)
		go xmlStasher(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all stashers are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

func CreateFetchers(stash string, zipp bool, inp <-chan Extract) <-chan Extract {

	if inp == nil {
		return nil
	}

	out := make(chan Extract, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create fetcher channel\n")
		os.Exit(1)
	}

	sfx := ".xml"
	if zipp {
		sfx = ".xml.gz"
	}

	fetchRecord := func(file string, buf bytes.Buffer) string {

		var arry [132]rune
		trie := MakeArchiveTrie(file, arry)

		if file == "" || trie == "" {
			return ""
		}

		fpath := path.Join(stash, trie, file+sfx)
		if fpath == "" {
			return ""
		}

		iszip := zipp

		inFile, err := os.Open(fpath)

		// if failed to find ".xml" file, try ".xml.gz" without requiring -gzip
		if err != nil && os.IsNotExist(err) && !zipp {
			iszip = true
			fpath := path.Join(stash, trie, file+".xml.gz")
			if fpath == "" {
				return ""
			}
			inFile, err = os.Open(fpath)
		}
		if err != nil {
			msg := err.Error()
			if !strings.HasSuffix(msg, "no such file or directory") && !strings.HasSuffix(msg, "cannot find the path specified.") {
				fmt.Fprintf(os.Stderr, "%s\n", msg)
			}
			return ""
		}

		defer inFile.Close()

		brd := bufio.NewReader(inFile)

		if iszip {

			zpr, err := gzip.NewReader(brd)

			defer zpr.Close()

			if err == nil {
				// copy and decompress cached file contents
				buf.ReadFrom(zpr)
			}

		} else {

			// copy cached file contents
			buf.ReadFrom(brd)
		}

		str := buf.String()

		return str
	}

	// xmlFetcher reads XML from file
	xmlFetcher := func(wg *sync.WaitGroup, inp <-chan Extract, out chan<- Extract) {

		// report when more records to process
		defer wg.Done()

		var buf bytes.Buffer

		for ext := range inp {

			buf.Reset()

			str := fetchRecord(ext.Text, buf)

			out <- Extract{ext.Index, "", str, nil}
		}
	}

	var wg sync.WaitGroup

	// launch multiple fetcher goroutines
	for i := 0; i < NumServe; i++ {
		wg.Add(1)
		go xmlFetcher(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all fetchers are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

func CreateStreamers(stash string, inp <-chan Extract) <-chan Extract {

	if inp == nil {
		return nil
	}

	out := make(chan Extract, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create streamer channel\n")
		os.Exit(1)
	}

	sfx := ".xml.gz"

	getRecord := func(file string, buf bytes.Buffer) []byte {

		var arry [132]rune
		trie := MakeArchiveTrie(file, arry)

		if file == "" || trie == "" {
			return nil
		}

		fpath := path.Join(stash, trie, file+sfx)
		if fpath == "" {
			return nil
		}

		inFile, err := os.Open(fpath)

		if err != nil {
			msg := err.Error()
			if !strings.HasSuffix(msg, "no such file or directory") && !strings.HasSuffix(msg, "cannot find the path specified.") {
				fmt.Fprintf(os.Stderr, "%s\n", msg)
			}
			return nil
		}

		defer inFile.Close()

		brd := bufio.NewReader(inFile)

		// copy cached file contents
		buf.ReadFrom(brd)

		data := buf.Bytes()

		return data
	}

	// xmlStreamer reads compressed XML from file
	xmlStreamer := func(wg *sync.WaitGroup, inp <-chan Extract, out chan<- Extract) {

		// report when more records to process
		defer wg.Done()

		var buf bytes.Buffer

		for ext := range inp {

			buf.Reset()

			data := getRecord(ext.Text, buf)

			out <- Extract{ext.Index, "", "", data}
		}
	}

	var wg sync.WaitGroup

	// launch multiple streamer goroutines
	for i := 0; i < NumServe; i++ {
		wg.Add(1)
		go xmlStreamer(&wg, inp, out)
	}

	// launch separate anonymous goroutine to wait until all streamers are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

func CreateInverters(args []string, nvrt string) <-chan []string {

	if args == nil {
		return nil
	}

	out := make(chan []string, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create inverter channel array\n")
		os.Exit(1)
	}

	// mutex for inverted index
	var ilock sync.Mutex

	// map for inverted index
	inverted := make(map[string][]string)

	// add single posting
	addPost := func(term, uid string) {

		// protect with mutex
		ilock.Lock()

		defer ilock.Unlock()

		data, ok := inverted[term]
		if !ok {
			data = make([]string, 0, 2)
			// first entry on new slice is term
			data = append(data, term)
		}
		data = append(data, uid)
		// always need to update inverted, since data may be reallocated
		inverted[term] = data
	}

	// xmlInverter sends UID and term strings to a callback for inversion
	xmlInverter := func(wg *sync.WaitGroup, fileName string) {

		defer wg.Done()

		f, err := os.Open(fileName)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
			os.Exit(1)
		}

		// close input file when all records have been processed
		defer f.Close()

		var in io.Reader

		in = f

		// if suffix is ".gz", use decompressor
		iszip := false
		if strings.HasSuffix(fileName, ".gz") {
			iszip = true
		}

		if iszip {
			brd := bufio.NewReader(f)
			if brd == nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create buffered reader on '%s'\n", fileName)
				os.Exit(1)
			}
			zpr, err := gzip.NewReader(brd)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create decompressor on '%s'\n", fileName)
				os.Exit(1)
			}

			// close decompressor when all records have been processed
			defer zpr.Close()

			// use decompressor for reading file
			in = zpr
		}

		rdr := NewXMLReader(in)

		if rdr == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
			os.Exit(1)
		}

		currUID := ""

		doInvert := func(tag, content string) {

			if tag == "IdxUid" {
				currUID = content
			} else if tag == nvrt {

				// expand Greek letters, anglicize characters in other alphabets
				if IsNotASCII(content) {

					if HasGreek(content) {
						content = SpellGreek(content)
						content = CompressRunsOfSpaces(content)
					}
					content = unidecode.Unidecode(content)

					content = DoAccentTransform(content)
					content = UnicodeToASCII(content)

					content = strings.TrimSpace(content)
				}

				// remove punctuation from term
				content := strings.Map(func(c rune) rune {
					if !unicode.IsLetter(c) && !unicode.IsDigit(c) && c != ' ' && c != '-' && c != '_' {
						return -1
					}
					return c
				}, content)

				content = strings.Replace(content, "_", " ", -1)
				content = strings.Replace(content, "-", " ", -1)

				content = CompressRunsOfSpaces(content)
				content = strings.TrimSpace(content)

				if content != "" && currUID != "" {
					addPost(content, currUID)
				}
			}
		}

		// partition all input by pattern and send XML substring through channel
		PartitionPattern("IdxDocument", "", rdr,
			func(rec int, ofs int64, str string) {
				StreamValues(str[:], "IdxDocument", doInvert)
			})
	}

	var wg sync.WaitGroup

	// launch multiple inverter goroutines
	for _, str := range args {
		wg.Add(1)
		go xmlInverter(&wg, str)
	}

	// launch separate anonymous goroutine to wait until all inverters are done
	go func() {
		wg.Wait()

		// send results to sorters
		for _, data := range inverted {
			out <- data
		}

		close(out)
	}()

	return out
}

func CreateSorters(nvrt string, inp <-chan []string) <-chan Extract {

	if inp == nil {
		return nil
	}

	out := make(chan Extract, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create sorter channel\n")
		os.Exit(1)
	}

	// xmlSorter reads from channel and prints one sorted posting list
	xmlSorter := func(wg *sync.WaitGroup, nvrt string, inp <-chan []string, out chan<- Extract) {

		defer wg.Done()

		var buffer strings.Builder

		printPosting := func(key string, data []string) string {

			if len(data) > 1 {
				sort.Slice(data, func(i, j int) bool {
					// numeric sort on strings checks lengths first
					lni := len(data[i])
					lnj := len(data[j])
					// shorter string is numerically less, assuming no leading zeros
					if lni < lnj {
						return true
					}
					if lni > lnj {
						return false
					}
					// same length, can now do string comparison on contents
					return data[i] < data[j]
				})
			}

			buffer.Reset()

			buffer.WriteString("  <InvDocument>\n")
			buffer.WriteString("    <InvKey>")
			buffer.WriteString(key)
			buffer.WriteString("</InvKey>\n")

			// print list of UIDs, skipping duplicates
			buffer.WriteString("    <InvIDs>\n")
			prev := ""
			for _, uid := range data {
				if uid == prev {
					continue
				}

				buffer.WriteString("      <")
				buffer.WriteString(nvrt)
				buffer.WriteString(">")
				buffer.WriteString(uid)
				buffer.WriteString("</")
				buffer.WriteString(nvrt)
				buffer.WriteString(">\n")

				prev = uid
			}
			buffer.WriteString("    </InvIDs>\n")
			buffer.WriteString("  </InvDocument>\n")

			str := buffer.String()

			return str
		}

		for inv := range inp {

			key := inv[0]
			data := inv[1:]

			str := printPosting(key, data)

			out <- Extract{0, key, str, nil}
		}
	}

	var wg sync.WaitGroup

	// launch multiple sorter goroutines
	for i := 0; i < NumServe; i++ {
		wg.Add(1)
		go xmlSorter(&wg, nvrt, inp, out)
	}

	// launch separate anonymous goroutine to wait until all sorters are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

func CreateResolver(inp <-chan Extract) <-chan string {

	if inp == nil {
		return nil
	}

	out := make(chan string, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create resolver channel\n")
		os.Exit(1)
	}

	// xmlResolver distributes adjacent records with the same identifier prefix
	xmlResolver := func(inp <-chan Extract, out chan<- string) {

		// close channel when all records have been processed
		defer close(out)

		// map for inverted index
		inverted := make(map[string]string)

		// drain channel, populate map for alphabetizing
		for curr := range inp {

			inverted[curr.Ident] = curr.Text
		}

		// sort map keys in alphabetical order
		ordered := make([]Extract, len(inverted))

		i := 0
		for key, str := range inverted {
			ordered[i].Ident = key
			ordered[i].Text = str
			i++
		}
		sort.Slice(ordered, func(i, j int) bool { return ordered[i].Ident < ordered[j].Ident })

		// iterate through alphabetized results
		for _, curr := range ordered {

			// send result to output
			out <- curr.Text

			runtime.Gosched()
		}
	}

	// launch single resolver goroutine
	go xmlResolver(inp, out)

	return out
}

type Plex struct {
	Which int
	Ident string
	Text  string
	Versd bool
}

type PlexHeap []Plex

// methods that satisfy heap.Interface
func (h PlexHeap) Len() int {
	return len(h)
}
func (h PlexHeap) Less(i, j int) bool {
	return h[i].Ident < h[j].Ident
}
func (h PlexHeap) Swap(i, j int) {
	h[i], h[j] = h[j], h[i]
}
func (h *PlexHeap) Push(x interface{}) {
	*h = append(*h, x.(Plex))
}
func (h *PlexHeap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

func CreatePresenters(args []string) []<-chan Plex {

	if args == nil {
		return nil
	}

	numFiles := len(args)
	if numFiles < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Not enough inverted files to merge\n")
		os.Exit(1)
	}

	chns := make([]<-chan Plex, numFiles)
	if chns == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create presenter channel array\n")
		os.Exit(1)
	}

	// xmlPresenter sends partitioned XML strings through channel
	xmlPresenter := func(fileNum int, fileName string, out chan<- Plex) {

		// close channel when all records have been processed
		defer close(out)

		f, err := os.Open(fileName)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
			os.Exit(1)
		}

		// close input file when all records have been processed
		defer f.Close()

		var in io.Reader

		in = f

		isVersioned := false
		if strings.Contains(fileName, "versioned") {
			isVersioned = true
		}

		// if suffix is ".gz", use decompressor
		iszip := false
		if strings.HasSuffix(fileName, ".gz") {
			iszip = true
		}

		if iszip {
			brd := bufio.NewReader(f)
			if brd == nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create buffered reader on '%s'\n", fileName)
				os.Exit(1)
			}
			zpr, err := gzip.NewReader(brd)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create decompressor on '%s'\n", fileName)
				os.Exit(1)
			}

			// close decompressor when all records have been processed
			defer zpr.Close()

			// use decompressor for reading file
			in = zpr
		}

		rdr := NewXMLReader(in)

		if rdr == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
			os.Exit(1)
		}

		find := ParseIndex("InvKey")

		// partition all input by pattern and send XML substring through channel
		PartitionPattern("InvDocument", "", rdr,
			func(rec int, ofs int64, str string) {
				id := FindIdentifier(str[:], "InvDocument", find)

				out <- Plex{fileNum, id, str, isVersioned}
			})
	}

	// launch multiple presenter goroutines
	for i, str := range args {

		chn := make(chan Plex, ChanDepth)
		if chn == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create presenter channel\n")
			os.Exit(1)
		}

		go xmlPresenter(i, str, chn)

		chns[i] = chn
	}

	// no need for separate anonymous goroutine to wait until all presenters are done

	return chns
}

func CreateManifold(inp []<-chan Plex) <-chan Plex {

	if inp == nil {
		return nil
	}

	out := make(chan Plex, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create manifold channel\n")
		os.Exit(1)
	}

	// xmlManifold restores original order with heap
	xmlManifold := func(inp []<-chan Plex, out chan<- Plex) {

		// close channel when all records have been processed
		defer close(out)

		// initialize empty heap
		hp := &PlexHeap{}
		heap.Init(hp)

		// read first object from all input channels in turn
		for _, chn := range inp {
			plx, ok := <-chn
			if ok {
				heap.Push(hp, plx)
			}
		}

		// reading from heap returns objects in alphabetical order
		for hp.Len() > 0 {

			// remove lowest item from heap, use interface type assertion
			curr := heap.Pop(hp).(Plex)

			out <- Plex{curr.Which, curr.Ident, curr.Text, curr.Versd}

			// read next object from channel that just supplied lowest item
			chn := inp[curr.Which]
			plx, ok := <-chn
			if ok {
				heap.Push(hp, plx)
			}
		}
	}

	// launch single manifold goroutine
	go xmlManifold(inp, out)

	return out
}

func CreateMerger(field, dltd string, inp <-chan Plex) <-chan Plex {

	if inp == nil {
		return nil
	}

	out := make(chan Plex, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create merger channel\n")
		os.Exit(1)
	}

	deletedUIDs := make(map[string]bool)

	readDeletedUIDs := func() {
		in, err := os.Open(dltd)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return
		}
		defer in.Close()

		scanr := bufio.NewScanner(in)
		for scanr.Scan() {

			// read lines of identifiers
			uid := scanr.Text()
			// store in map of deleted UIDs
			deletedUIDs[uid] = true
		}
	}

	checkDeleted := false

	if dltd != "" {
		readDeletedUIDs()
		checkDeleted = true
	}

	fuseInverts := func(key string, arry []string) string {

		var data []string

		addUID := func(tag, content string) {

			if tag == field {

				data = append(data, content)
			}
		}

		filterUID := func(tag, content string) {

			if tag == field {

				if deletedUIDs[content] {
					// skip UID if in deleted list
					return
				}
				data = append(data, content)
			}
		}
		for _, str := range arry {
			// first character is Y for versioned.inv, N for all others
			prefix := str[0]
			str = str[1:]
			if checkDeleted && prefix == 'N' {
				StreamValues(str[:], "InvDocument", filterUID)
			} else {
				StreamValues(str[:], "InvDocument", addUID)
			}
		}

		if len(data) > 1 {
			// sort each posting list in numeric order
			sort.Slice(data, func(i, j int) bool {
				// numeric sort on strings checks lengths first
				lni := len(data[i])
				lnj := len(data[j])
				// shorter string is numerically less, assuming no leading zeros
				if lni < lnj {
					return true
				}
				if lni > lnj {
					return false
				}
				// same length, can now do string comparison on contents
				return data[i] < data[j]
			})
		}

		var buffer strings.Builder

		buffer.WriteString("  <InvDocument>\n")
		buffer.WriteString("    <InvKey>")
		buffer.WriteString(key)
		buffer.WriteString("</InvKey>\n")

		// print list of UIDs, skipping duplicates
		buffer.WriteString("    <InvIDs>\n")
		last := ""
		for _, uid := range data {
			// detect duplicate UIDs, now in same list after conversion of one term entry from foreign alphabet
			if uid == last {
				continue
			}
			buffer.WriteString("      <")
			buffer.WriteString(field)
			buffer.WriteString(">")
			buffer.WriteString(uid)
			buffer.WriteString("</")
			buffer.WriteString(field)
			buffer.WriteString(">\n")

			last = uid
		}
		buffer.WriteString("    </InvIDs>\n")
		buffer.WriteString("  </InvDocument>\n")

		txt := buffer.String()

		return txt
	}

	// xmlMerger fuses adjacent inverted records with the same identifier
	xmlMerger := func(inp <-chan Plex, out chan<- Plex) {

		// close channel when all records have been processed
		defer close(out)

		// remember previous record
		prev := Plex{}

		// array to collect strings with same identifier
		var arry []string

		for curr := range inp {

			if curr.Text == "" {
				continue
			}

			prefix := "N"
			if curr.Versd {
				prefix = "Y"
			}

			// compare adjacent record identifiers
			if prev.Ident == curr.Ident {

				// save next string in slice
				arry = append(arry, prefix+curr.Text)

			} else {

				if len(arry) > 0 {

					// if next identifier is different, fuse all saved components
					prev.Text = fuseInverts(prev.Ident, arry)

					// and send combined results to output channel
					out <- prev

					// empty the slice
					arry = nil
				}

				// now remember this record
				prev = curr

				// and save its string as the first entry in the slice
				arry = append(arry, prefix+curr.Text)
			}
		}

		if len(arry) > 0 {

			// fuse remaining saved components
			prev.Text = fuseInverts(prev.Ident, arry)

			// send last record
			out <- prev

			arry = nil
		}
	}

	// launch single merger goroutine
	go xmlMerger(inp, out)

	return out
}

func CreateSplitter(merg, field string, zipp bool, inp <-chan Plex) <-chan string {

	if inp == nil {
		return nil
	}

	out := make(chan string, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create splitter channel\n")
		os.Exit(1)
	}

	openSaver := func(merg, field, key string, zipp bool) (*os.File, *bufio.Writer, *gzip.Writer) {

		var (
			fl   *os.File
			wrtr *bufio.Writer
			zpr  *gzip.Writer
			err  error
		)

		sfx := ".mrg"
		if zipp {
			sfx = ".mrg.gz"
		}

		fpath := path.Join(merg, field, key+sfx)
		if fpath == "" {
			return nil, nil, nil
		}

		fl, err = os.OpenFile(fpath, os.O_CREATE|os.O_WRONLY, 0600)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			return nil, nil, nil
		}

		var out io.Writer

		out = fl

		if zipp {

			zpr, err = gzip.NewWriterLevel(fl, gzip.BestSpeed)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
				return nil, nil, nil
			}

			out = zpr
		}

		// create buffered writer layer
		wrtr = bufio.NewWriter(out)

		return fl, wrtr, zpr
	}

	closeSaver := func(fl *os.File, wrtr *bufio.Writer, zpr *gzip.Writer) {

		wrtr.Flush()
		if zpr != nil {
			zpr.Close()
		}
		fl.Sync()
	}

	// xmlSplitter distributes adjacent records with the same identifier prefix
	xmlSplitter := func(inp <-chan Plex, out chan<- string) {

		// close channel when all records have been processed
		defer close(out)

		var (
			fl   *os.File
			wrtr *bufio.Writer
			zpr  *gzip.Writer
		)

		currTag := ""
		prevTag := ""

		// remember previous record
		prev := Plex{}

		for curr := range inp {

			// use first few characters of identifier
			currTag = IdentifierKey(curr.Ident)
			if currTag == "" {
				continue
			}

			// then truncate to 23 character prefix
			if len(currTag) > 2 {
				currTag = currTag[:2]
			}

			if fl == nil {
				// open initial file
				fl, wrtr, zpr = openSaver(merg, field, currTag, zipp)

				// send first opening tag and indent
				wrtr.WriteString("<InvDocumentSet>\n  ")
			}

			// compare keys from adjacent term lists
			if prev.Text != "" && prevTag != currTag {

				// after IdentifierKey converts space to underscore,
				// okay that x_ and x0 will be out of alphabetical order

				// send closing tag
				wrtr.WriteString("</InvDocumentSet>\n")

				closeSaver(fl, wrtr, zpr)

				out <- currTag

				// open next file
				fl, wrtr, zpr = openSaver(merg, field, currTag, zipp)

				// send opening tag and indent
				wrtr.WriteString("<InvDocumentSet>\n  ")
			}

			// send one InvDocument
			str := strings.TrimSpace(curr.Text)

			wrtr.WriteString(str)
			if !strings.HasSuffix(str, "\n") {
				wrtr.WriteString("\n")
			}

			// now remember this record
			prev = curr

			prevTag = currTag
		}

		if prev.Text != "" {

			// send last closing tag
			wrtr.WriteString("</InvDocumentSet>\n")

			closeSaver(fl, wrtr, zpr)

			out <- currTag
		}
	}

	// launch single splitter goroutine
	go xmlSplitter(inp, out)

	return out
}

func CreatePromoters(args []string, prom, field string) <-chan string {

	if args == nil {
		return nil
	}

	out := make(chan string, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create promoter channel\n")
		os.Exit(1)
	}

	// xmlPromoter splits inverted index groups into subdirectories for term lists and postings files
	xmlPromoter := func(wg *sync.WaitGroup, fileName string, out chan<- string) {

		defer wg.Done()

		f, err := os.Open(fileName)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open input file '%s'\n", fileName)
			os.Exit(1)
		}

		// close input file when all records have been processed
		defer f.Close()

		var in io.Reader

		in = f

		// if suffix is ".gz", use decompressor
		iszip := false
		if strings.HasSuffix(fileName, ".gz") {
			iszip = true
		}

		if iszip {
			brd := bufio.NewReader(f)
			if brd == nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create buffered reader on '%s'\n", fileName)
				os.Exit(1)
			}
			zpr, err := gzip.NewReader(brd)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create decompressor on '%s'\n", fileName)
				os.Exit(1)
			}

			// close decompressor when all records have been processed
			defer zpr.Close()

			// use decompressor for reading file
			in = zpr
		}

		rdr := NewXMLReader(in)

		if rdr == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
			os.Exit(1)
		}

		getOnePosting := func(text string) (string, []int32) {

			var data []int32

			term := ""

			doPromote := func(tag, content string) {

				if tag == "InvKey" {

					// term used for postings file name
					term = content

					term = strings.ToLower(term)

				} else if tag == field {

					// convert UID string to integer
					value, err := strconv.ParseInt(content, 10, 32)
					if err != nil {
						fmt.Fprintf(os.Stderr, "%s\n", err.Error())
						return
					}
					data = append(data, int32(value))
				}
			}

			// explore data fields
			StreamValues(text[:], "InvDocument", doPromote)

			if term == "" || len(data) < 1 {
				return "", nil
			}

			return term, data
		}

		var termPos int32
		var postPos int32

		var indxList bytes.Buffer
		var termList bytes.Buffer
		var postList bytes.Buffer

		retlength := len("\n")

		addOnePosting := func(term string, data []int32) {

			dlength := len(data)
			tlength := len(term)

			// write to term list buffer
			termList.WriteString(term[:])
			termList.WriteString("\n")

			// write to postings buffer
			binary.Write(&postList, binary.LittleEndian, data)

			// write to master index buffer
			binary.Write(&indxList, binary.LittleEndian, termPos)
			binary.Write(&indxList, binary.LittleEndian, postPos)

			postPos += int32(dlength * 4)
			termPos += int32(tlength + retlength)
		}

		topOffMaster := func() {

			// phantom term and postings positions eliminates special case calculation at end
			binary.Write(&indxList, binary.LittleEndian, termPos)
			binary.Write(&indxList, binary.LittleEndian, postPos)
		}

		writeFile := func(dpath, fname string, bfr bytes.Buffer) {

			fpath := path.Join(dpath, fname)
			if fpath == "" {
				return
			}

			fl, err := os.OpenFile(fpath, os.O_CREATE|os.O_WRONLY, 0600)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
				return
			}

			defer fl.Close()

			defer fl.Sync()

			data := bfr.Bytes()

			wrtr := bufio.NewWriter(fl)

			_, err = wrtr.Write(data)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			}

			wrtr.Flush()
		}

		writeThreeFiles := func(key string) {

			var arry [516]rune
			dpath, key := PostingPath(prom, field, key, arry)
			if dpath == "" {
				return
			}

			// make subdirectories, if necessary
			_, err := os.Stat(dpath)
			if err != nil && os.IsNotExist(err) {
				err = os.MkdirAll(dpath, os.ModePerm)
			}
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
				return
			}

			writeFile(dpath, key+".trm", termList)

			writeFile(dpath, key+".pst", postList)

			writeFile(dpath, key+".mst", indxList)
		}

		currTag := ""
		prevTag := ""

		ok := false

		PartitionPattern("InvDocument", "", rdr,
			func(rec int, ofs int64, str string) {

				term, data := getOnePosting(str)

				if term == "" || data == nil {
					return
				}

				ok = true

				// use first few characters of identifier
				currTag = IdentifierKey(term)

				if prevTag != currTag {

					// after IdentifierKey converts space to underscore,
					// okay that xxx_ and xxx0 will be out of alphabetical order

					// directory prefix changed from last posting
					if prevTag != "" {

						topOffMaster()
						writeThreeFiles(prevTag)
						out <- prevTag
					}

					// reset buffers and position counters
					termPos = 0
					postPos = 0

					indxList.Reset()
					termList.Reset()
					postList.Reset()

				}

				addOnePosting(term, data)

				prevTag = currTag
			})

		if ok {

			// write last set of files
			topOffMaster()
			writeThreeFiles(prevTag)
			out <- prevTag
		}
	}

	var wg sync.WaitGroup

	// launch multiple promoter goroutines
	for _, str := range args {
		wg.Add(1)
		go xmlPromoter(&wg, str, out)
	}

	// launch separate anonymous goroutine to wait until all promoters are done
	go func() {
		wg.Wait()
		close(out)
	}()

	return out
}

// MAIN FUNCTION

func main() {

	// skip past executable name
	args := os.Args[1:]

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: No command-line arguments supplied to rchive\n")
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
	NumServe = 0

	// number of channels usually equals number of servers, but can be overridden by -chan
	ChanDepth = 0

	// miscellaneous tuning parameters
	HeapSize = 16
	FarmSize = 64

	// garbage collector control can be set by environment variable or default value with -gogc 0
	goGc := 600

	// -flag sets -strict or -mixed cleanup flags from argument
	flgs := ""

	DeStop = true

	// read data from file instead of stdin
	fileName := ""

	// debugging
	timr := false

	// profiling
	prfl := false

	// element to use as local data index
	indx := ""

	// file of index values for removing duplicates
	unqe := ""

	// path for local data indexed as trie
	stsh := ""
	ftch := ""
	strm := ""

	// field for inverted index
	nvrt := ""

	// field for merging index files
	merg := ""

	// file of deleted and versioned UIDs
	dltd := ""

	// base for promoting inverted index to retrieval indices
	prom := ""

	// field for promoting inverted index files
	fild := ""

	// base for queries
	base := ""

	// query by phrase, normalized terms (with truncation wildcarding), optionally using stemmed terms
	phrs := ""

	// print term list with counts
	trms := ""
	plrl := false

	// use gzip compression on local data files
	zipp := false

	// print UIDs and hash values
	hshv := false

	// convert UIDs to archive trie
	trei := false

	// compare input record against stash
	cmpr := false
	cmprType := ""
	ignr := ""

	// flag missing identifiers
	msng := false

	// flag records with damaged embedded HTML tags
	dmgd := false
	dmgdType := ""

	// kludge to use non-threaded fetching for windows
	windows := false

	// get numeric value
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
			NumServe = getNumericArg("Concurrent parser count", 0, ncpu, 128)
		case "-chan":
			ChanDepth = getNumericArg("Communication channel depth", 0, ncpu, 128)
		case "-heap":
			HeapSize = getNumericArg("Unshuffler heap size", 8, 8, 64)
		case "-farm":
			FarmSize = getNumericArg("Node buffer length", 4, 4, 2048)
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

		// file with selected indexes for removing duplicates
		case "-unique":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Unique identifier file is missing\n")
				os.Exit(1)
			}
			unqe = args[1]
			// skip past first of two arguments
			args = args[1:]

		// local directory path for indexing
		case "-archive", "-stash":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Archive path is missing\n")
				os.Exit(1)
			}
			stsh = args[1]
			if stsh != "" && !strings.HasSuffix(stsh, "/") {
				stsh += "/"
			}
			// skip past first of two arguments
			args = args[1:]
		// local directory path for retrieval
		case "-fetch":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Fetch path is missing\n")
				os.Exit(1)
			}
			ftch = args[1]
			if ftch != "" && !strings.HasSuffix(ftch, "/") {
				ftch += "/"
			}
			// skip past first of two arguments
			args = args[1:]
		// local directory path for retrieval of compressed XML
		case "-stream":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Stream path is missing\n")
				os.Exit(1)
			}
			strm = args[1]
			if strm != "" && !strings.HasSuffix(strm, "/") {
				strm += "/"
			}
			// skip past first of two arguments
			args = args[1:]

		// data element for indexing
		case "-index":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Index element is missing\n")
				os.Exit(1)
			}
			indx = args[1]
			// skip past first of two arguments
			args = args[1:]

		// build inverted index
		case "-invert":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Invert field is missing\n")
				os.Exit(1)
			}
			nvrt = args[1]
			// skip past first of two arguments
			args = args[1:]

		// merge inverted index files, distribute by prefix
		case "-merge":
			if len(args) < 3 {
				fmt.Fprintf(os.Stderr, "\nERROR: Merge field is missing\n")
				os.Exit(1)
			}
			merg = args[1]
			fild = args[2]
			// skip past first and second arguments
			args = args[2:]

		// file of deleted and versioned UIDs for merging
		case "-deleted":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Deleted and versioned PMID file is missing\n")
				os.Exit(1)
			}
			dltd = args[1]
			// skip past first of two arguments
			args = args[1:]

		// promote inverted index to term-specific postings files
		case "-promote":
			if len(args) < 3 {
				fmt.Fprintf(os.Stderr, "\nERROR: Promote path is missing\n")
				os.Exit(1)
			}
			prom = args[1]
			fild = args[2]
			// skip past first and second arguments
			args = args[2:]

		case "-query":
			if len(args) < 3 {
				fmt.Fprintf(os.Stderr, "\nERROR: Postings path is missing\n")
				os.Exit(1)
			}
			base = args[1]
			phrs = args[2]
			// skip past first and second arguments
			args = args[2:]

		case "-counts":
			plrl = true
			fallthrough
		case "-count":
			if len(args) < 3 {
				fmt.Fprintf(os.Stderr, "\nERROR: Postings path is missing\n")
				os.Exit(1)
			}
			base = args[1]
			trms = args[2]
			// skip past first and second arguments
			args = args[2:]

		case "-gzip":
			zipp = true
		case "-hash":
			hshv = true
		case "-trie":
			trei = true
		// check for missing records
		case "-missing":
			msng = true

		// use non-threaded fetch function for windows (undocumented)
		case "-windows":
			windows = true

		// data cleanup flags
		case "-compress", "-compressed":
			DoCompress = true
		case "-spaces", "-cleanup":
			DoCleanup = true
		case "-strict":
			DoStrict = true
		case "-mixed", "-relaxed":
			DoMixed = true
		case "-accent", "-plain":
			DeAccent = true
		case "-ascii":
			DoASCII = true

		case "-stems", "-stem":
			DoStem = true
		case "-stops", "-stop":
			DeStop = false

		case "-unicode", "-repair":
			DoUnicode = true
		case "-script":
			DoScript = true
		case "-mathml":
			DoMathML = true

		case "-flag", "-flags":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: Flags argument is missing\n")
				os.Exit(1)
			}
			flgs = args[1]
			// skip past first of two arguments
			args = args[1:]

		// debugging flags
		case "-damaged", "-damage", "-broken":
			dmgd = true
			if len(args) > 1 {
				next := args[1]
				// if next argument is not another flag
				if next != "" && next[0] != '-' {
					// get optional extraction class (SELF, SINGLE, DOUBLE, AMPER, or ALL)
					dmgdType = next
					// skip past first of two arguments
					args = args[1:]
				}
			}
		case "-prepare":
			cmpr = true
			if len(args) > 1 {
				next := args[1]
				// if next argument is not another flag
				if next != "" && next[0] != '-' {
					// get optional data source specifier
					cmprType = next
					// skip past first of two arguments
					args = args[1:]
				}
			}
		case "-ignore":
			if len(args) < 2 {
				fmt.Fprintf(os.Stderr, "\nERROR: -ignore value is missing\n")
				os.Exit(1)
			}
			ignr = args[1]
			// skip past first of two arguments
			args = args[1:]

		case "-timer":
			timr = true
		case "-profile":
			prfl = true

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

	// -flag allows script to set -strict or -mixed (or -stems, or -stops) from argument
	switch flgs {
	case "strict":
		DoStrict = true
	case "mixed":
		DoMixed = true
	case "stems", "stem":
		DoStem = true
	case "stops", "stop":
		DeStop = false
	case "none", "default":
	default:
		if flgs != "" {
			fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized -flag value '%s'\n", flgs)
			os.Exit(1)
		}
	}

	AllowEmbed = DoStrict || DoMixed
	ContentMods = AllowEmbed || DoCompress || DoUnicode || DoScript || DoMathML || DeAccent || DoASCII

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

	// allow simultaneous threads for multiplexed goroutines
	runtime.GOMAXPROCS(numProcs)

	// adjust garbage collection target percentage
	if goGc >= 100 {
		debug.SetGCPercent(goGc)
	}

	// explicit -serv argument overrides -cons ratio
	if NumServe > 0 {
		serverRatio = NumServe / numProcs
		// if numServers / numProcs is not a whole number, do not print serverRatio in -stats
		if NumServe != numProcs*serverRatio {
			serverRatio = 0
		}
	} else {
		NumServe = numProcs * serverRatio
	}
	// server limits
	if NumServe > 128 {
		NumServe = 128
	} else if NumServe < 1 {
		NumServe = numProcs
	}

	// explicit -chan argument overrides default to number of servers
	if ChanDepth == 0 {
		ChanDepth = NumServe
	}

	// if copying from local files accessed by identifier, add dummy argument to bypass length tests
	if stsh != "" && indx == "" {
		args = append(args, "-dummy")
	} else if ftch != "" || strm != "" {
		args = append(args, "-dummy")
	} else if base != "" {
		args = append(args, "-dummy")
	} else if trei || dmgd || cmpr {
		args = append(args, "-dummy")
	}

	// expand -archive ~/ to home directory path
	if stsh != "" {

		if stsh[:2] == "~/" {
			cur, err := user.Current()
			if err == nil {
				hom := cur.HomeDir
				stsh = strings.Replace(stsh, "~/", hom+"/", 1)
			}
		}
	}

	// expand -fetch ~/ to home directory path
	if ftch != "" {

		if ftch[:2] == "~/" {
			cur, err := user.Current()
			if err == nil {
				hom := cur.HomeDir
				ftch = strings.Replace(ftch, "~/", hom+"/", 1)
			}
		}
	}

	// expand -stream ~/ to home directory path
	if strm != "" {

		if strm[:2] == "~/" {
			cur, err := user.Current()
			if err == nil {
				hom := cur.HomeDir
				strm = strings.Replace(strm, "~/", hom+"/", 1)
			}
		}
	}

	// expand -promote ~/ to home directory path
	if prom != "" {

		if prom[:2] == "~/" {
			cur, err := user.Current()
			if err == nil {
				hom := cur.HomeDir
				prom = strings.Replace(prom, "~/", hom+"/", 1)
			}
		}
	}

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to rchive\n")
		os.Exit(1)
	}

	// DOCUMENTATION COMMANDS

	inSwitch = true

	switch args[0] {
	case "-version":
		fmt.Printf("%s\n", rchiveVersion)
	case "-help":
		fmt.Printf("rchive %s\n%s\n", rchiveVersion, rchiveHelp)
	case "-extras", "-extra", "-advanced":
		fmt.Printf("rchive %s\n%s\n", rchiveVersion, rchiveExtras)
	case "-internal":
		fmt.Printf("rchive %s\n%s\n", rchiveVersion, rchiveInternal)
	default:
		// if not any of the documentation commands, keep going
		inSwitch = false
	}

	if inSwitch {
		return
	}

	// FILE NAME CAN BE SUPPLIED WITH -input COMMAND

	in := os.Stdin

	// check for data being piped into stdin
	isPipe := false
	fi, err := os.Stdin.Stat()
	if err == nil {
		isPipe = bool((fi.Mode() & os.ModeNamedPipe) != 0)
	}

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
	}

	// check for -input command after extraction arguments
	for _, str := range args {
		if str == "-input" {
			fmt.Fprintf(os.Stderr, "\nERROR: Misplaced -input command\n")
			os.Exit(1)
		}
	}

	// START PROFILING IF REQUESTED

	if prfl {

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

	// print processing rate and program duration
	printDuration := func(name string) {

		stopTime := time.Now()
		duration := stopTime.Sub(startTime)
		seconds := float64(duration.Nanoseconds()) / 1e9

		if recordCount >= 1000000 {
			throughput := float64(recordCount/100000) / 10.0
			fmt.Fprintf(os.Stderr, "\nRchive processed %.1f million %s in %.3f seconds", throughput, name, seconds)
		} else {
			fmt.Fprintf(os.Stderr, "\nRchive processed %d %s in %.3f seconds", recordCount, name, seconds)
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

	// ENTREZ INDEX INVERSION

	// -invert PAIR *.e2x.gz reads IdxDocumentSet XML and creates an inverted index
	if nvrt != "" {

		invq := CreateInverters(args, nvrt)
		srtq := CreateSorters(nvrt, invq)
		rslq := CreateResolver(srtq)

		if invq == nil || srtq == nil || rslq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create inverter\n")
			os.Exit(1)
		}

		var out io.Writer

		out = os.Stdout

		if zipp {

			zpr, err := gzip.NewWriterLevel(out, gzip.BestSpeed)
			if err != nil {
				fmt.Fprintf(os.Stderr, "\nERROR: Unable to create compressor\n")
				os.Exit(1)
			}

			// close decompressor when all records have been processed
			defer zpr.Close()

			// use compressor for writing file
			out = zpr
		}

		// create buffered writer layer
		wrtr := bufio.NewWriter(out)

		wrtr.WriteString("<InvDocumentSet>\n")

		// drain channel of alphabetized results
		for str := range rslq {

			// send result to output
			wrtr.WriteString(str)

			recordCount++
			runtime.Gosched()
		}

		wrtr.WriteString("</InvDocumentSet>\n\n")

		wrtr.Flush()

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	// MERGE INVERTED INDEX FILES AND GROUP BY TERM

	// -merge combines inverted files, distributes by prefix
	if merg != "" && fild != "" {

		chns := CreatePresenters(args)
		mfld := CreateManifold(chns)
		mrgr := CreateMerger(fild, dltd, mfld)
		sptr := CreateSplitter(merg, fild, zipp, mrgr)

		if chns == nil || mfld == nil || mrgr == nil || sptr == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create inverted index merger\n")
			os.Exit(1)
		}

		// drain channel, print two-character index name
		for str := range sptr {

			fmt.Fprintf(os.Stdout, "%s\n", str)

			recordCount++
			runtime.Gosched()
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	// PROMOTE MERGED INVERTED INDEX TO TERM LIST AND POSTINGS FILES

	if prom != "" && fild != "" {

		prmq := CreatePromoters(args, prom, fild)

		if prmq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create new postings file generator\n")
			os.Exit(1)
		}

		// drain channel, print 2-4 character file prefix
		for str := range prmq {

			fmt.Fprintf(os.Stdout, "%s\n", str)

			recordCount++
			runtime.Gosched()
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	// TEST POSTINGS FILES

	if base != "" && phrs != "" {

		// deStop should match value used in building the NORM, PAIR, STEM, and GRFT indices
		recordCount = ProcessQuery(base, phrs)

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	if base != "" && trms != "" {

		// deStop should match value used in building the NORM, PAIR, STEM, and GRFT indices
		if plrl {
			recordCount = ProcessCounts(base, trms)
		} else {
			recordCount = ProcessCount(base, trms)
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("terms")
		}

		return
	}

	// CREATE XML BLOCK READER FROM STDIN OR FILE

	rdr := NewXMLReader(in)
	if rdr == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create XML Block Reader\n")
		os.Exit(1)
	}

	// CONFIRM INPUT DATA AVAILABILITY AFTER RUNNING COMMAND GENERATORS

	if fileName == "" && runtime.GOOS != "windows" {

		fromStdin := bool((fi.Mode() & os.ModeCharDevice) == 0)
		if !isPipe || !fromStdin {
			mode := fi.Mode().String()
			fmt.Fprintf(os.Stderr, "\nERROR: No data supplied to rchive from stdin or file, mode is '%s'\n", mode)
			os.Exit(1)
		}
	}

	if !usingFile && !isPipe {

		fmt.Fprintf(os.Stderr, "\nERROR: No XML input data supplied to rchive\n")
		os.Exit(1)
	}

	// SPECIFY STRINGS TO GO BEFORE AND AFTER ENTIRE OUTPUT OR EACH RECORD

	head := ""
	tail := ""

	hd := ""
	tl := ""

	for {

		if len(args) < 1 {
			break
		}

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
	}

	// PRODUCE ARCHIVE SUBPATH FROM IDENTIFIER

	// -trie converts identifier to directory subpath plus file name (undocumented)
	if trei {

		scanr := bufio.NewScanner(rdr.Reader)

		sfx := ".xml"
		if zipp {
			sfx = ".xml.gz"
		}

		// read lines of identifiers
		for scanr.Scan() {

			file := scanr.Text()

			var arry [132]rune
			trie := MakeArchiveTrie(file, arry)
			if trie == "" || file == "" {
				continue
			}

			fpath := path.Join(trie, file+sfx)
			if fpath == "" {
				continue
			}

			os.Stdout.WriteString(fpath)
			os.Stdout.WriteString("\n")
		}

		return
	}

	// CHECK FOR MISSING RECORDS IN LOCAL DIRECTORY INDEXED BY TRIE ON IDENTIFIER

	// -archive plus -missing checks for missing records
	if stsh != "" && msng {

		scanr := bufio.NewScanner(rdr.Reader)

		sfx := ".xml"
		if zipp {
			sfx = ".xml.gz"
		}

		// read lines of identifiers
		for scanr.Scan() {

			file := scanr.Text()

			pos := strings.Index(file, ".")
			if pos >= 0 {
				// remove version suffix
				file = file[:pos]
			}

			var arry [132]rune
			trie := MakeArchiveTrie(file, arry)

			if file == "" || trie == "" {
				continue
			}

			fpath := path.Join(stsh, trie, file+sfx)
			if fpath == "" {
				continue
			}

			_, err := os.Stat(fpath)

			// if failed to find ".xml" file, try ".xml.gz" without requiring -gzip
			if err != nil && os.IsNotExist(err) && !zipp {
				fpath := path.Join(stsh, trie, file+".xml.gz")
				if fpath == "" {
					continue
				}
				_, err = os.Stat(fpath)
			}
			if err != nil && os.IsNotExist(err) {
				// record is missing from local file cache
				os.Stdout.WriteString(file)
				os.Stdout.WriteString("\n")
			}
		}

		return
	}

	// RETRIEVE XML COMPONENT RECORDS FROM LOCAL DIRECTORY INDEXED BY TRIE ON IDENTIFIER

	// alternative windows version limits memory by not using goroutines
	if ftch != "" && indx == "" && runtime.GOOS == "windows" && windows {

		scanr := bufio.NewScanner(rdr.Reader)
		if scanr == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create UID scanner\n")
			os.Exit(1)
		}

		sfx := ".xml"
		if zipp {
			sfx = ".xml.gz"
		}

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		var buf bytes.Buffer

		for scanr.Scan() {

			// read next identifier
			file := scanr.Text()

			pos := strings.Index(file, ".")
			if pos >= 0 {
				// remove version suffix
				file = file[:pos]
			}

			var arry [132]rune
			trie := MakeArchiveTrie(file, arry)

			if file == "" || trie == "" {
				continue
			}

			fpath := path.Join(ftch, trie, file+sfx)
			if fpath == "" {
				continue
			}

			iszip := zipp

			inFile, err := os.Open(fpath)

			// if failed to find ".xml" file, try ".xml.gz" without requiring -gzip
			if err != nil && os.IsNotExist(err) && !zipp {
				iszip = true
				fpath := path.Join(ftch, trie, file+".xml.gz")
				if fpath == "" {
					continue
				}
				inFile, err = os.Open(fpath)
			}
			if err != nil {
				continue
			}

			buf.Reset()

			brd := bufio.NewReader(inFile)

			if iszip {

				zpr, err := gzip.NewReader(brd)

				if err == nil {
					// copy and decompress cached file contents
					buf.ReadFrom(zpr)
				}

				zpr.Close()

			} else {

				// copy cached file contents
				buf.ReadFrom(brd)
			}

			inFile.Close()

			str := buf.String()

			if str == "" {
				continue
			}

			recordCount++

			if hd != "" {
				os.Stdout.WriteString(hd)
				os.Stdout.WriteString("\n")
			}

			if hshv {
				// calculate hash code for verification table
				hsh := crc32.NewIEEE()
				hsh.Write([]byte(str))
				val := hsh.Sum32()
				res := strconv.FormatUint(uint64(val), 10)
				txt := file + "\t" + res + "\n"
				os.Stdout.WriteString(txt)
			} else {
				// send result to output
				os.Stdout.WriteString(str)
				if !strings.HasSuffix(str, "\n") {
					os.Stdout.WriteString("\n")
				}
			}

			if tl != "" {
				os.Stdout.WriteString(tl)
				os.Stdout.WriteString("\n")
			}

			debug.FreeOSMemory()
		}

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// -fetch without -index retrieves XML files in trie-based directory structure
	if ftch != "" && indx == "" {

		uidq := CreateUIDReader(rdr.Reader)
		strq := CreateFetchers(ftch, zipp, uidq)
		unsq := CreateUnshuffler(strq)

		if uidq == nil || strq == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create archive reader\n")
			os.Exit(1)
		}

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		// drain output channel
		for curr := range unsq {

			str := curr.Text

			if str == "" {
				continue
			}

			recordCount++

			if hd != "" {
				os.Stdout.WriteString(hd)
				os.Stdout.WriteString("\n")
			}

			if hshv {
				// calculate hash code for verification table
				hsh := crc32.NewIEEE()
				hsh.Write([]byte(str))
				val := hsh.Sum32()
				res := strconv.FormatUint(uint64(val), 10)
				txt := curr.Ident + "\t" + res + "\n"
				os.Stdout.WriteString(txt)
			} else {
				// send result to output
				os.Stdout.WriteString(str)
				if !strings.HasSuffix(str, "\n") {
					os.Stdout.WriteString("\n")
				}
			}

			if tl != "" {
				os.Stdout.WriteString(tl)
				os.Stdout.WriteString("\n")
			}
		}

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// -stream without -index retrieves compressed XML files in trie-based directory structure
	if strm != "" && indx == "" {

		uidq := CreateUIDReader(rdr.Reader)
		strq := CreateStreamers(strm, uidq)
		unsq := CreateUnshuffler(strq)

		if uidq == nil || strq == nil || unsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create archive reader\n")
			os.Exit(1)
		}

		// drain output channel
		for curr := range unsq {

			data := curr.Data

			if data == nil {
				continue
			}

			recordCount++

			_, err = os.Stdout.Write(data)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			}
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// ENSURE PRESENCE OF PATTERN ARGUMENT

	if len(args) < 1 {
		fmt.Fprintf(os.Stderr, "\nERROR: Insufficient command-line arguments supplied to rchive\n")
		os.Exit(1)
	}

	// allow -record as synonym of -pattern (undocumented)
	if args[0] == "-record" || args[0] == "-Record" {
		args[0] = "-pattern"
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

	// REPORT RECORDS THAT CONTAIN DAMAGED EMBEDDED HTML TAGS

	// -damaged plus -index plus -pattern reports records with multiply-encoded HTML tags
	if dmgd && indx != "" {

		find := ParseIndex(indx)

		PartitionPattern(topPattern, star, rdr,
			func(rec int, ofs int64, str string) {
				recordCount++

				id := FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}

				// remove default version suffix
				if strings.HasSuffix(id, ".1") {
					idlen := len(id)
					id = id[:idlen-2]
				}

				ReportEncodedMarkup(dmgdType, id, str)
			})

		if timr {
			printDuration("records")
		}

		return
	}

	// COMPARE XML UPDATES TO LOCAL DIRECTORY, RETAIN NEW OR SUBSTANTIVELY CHANGED RECORDS

	// -prepare plus -archive plus -index plus -pattern compares XML files against stash
	if stsh != "" && indx != "" && cmpr {

		doReport := false
		if cmprType == "" || cmprType == "report" {
			doReport = true
		} else if cmprType != "release" {
			fmt.Fprintf(os.Stderr, "\nERROR: -prepare argument must be release or report\n")
			os.Exit(1)
		}

		find := ParseIndex(indx)

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		PartitionPattern(topPattern, star, rdr,
			func(rec int, ofs int64, str string) {
				recordCount++

				id := FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}

				pos := strings.Index(id, ".")
				if pos >= 0 {
					// remove version suffix
					id = id[:pos]
				}

				var arry [132]rune
				trie := MakeArchiveTrie(id, arry)

				if id == "" || trie == "" {
					return
				}

				fpath := path.Join(stsh, trie, id+".xml")
				if fpath == "" {
					return
				}

				// print new or updated XML record
				printRecord := func(stn string, isNew bool) {

					if stn == "" {
						return
					}

					if doReport {
						if isNew {
							os.Stdout.WriteString("NW ")
							os.Stdout.WriteString(id)
							os.Stdout.WriteString("\n")
						} else {
							os.Stdout.WriteString("UP ")
							os.Stdout.WriteString(id)
							os.Stdout.WriteString("\n")
						}
						return
					}

					if hd != "" {
						os.Stdout.WriteString(hd)
						os.Stdout.WriteString("\n")
					}

					os.Stdout.WriteString(stn)
					os.Stdout.WriteString("\n")

					if tl != "" {
						os.Stdout.WriteString(tl)
						os.Stdout.WriteString("\n")
					}
				}

				_, err := os.Stat(fpath)
				if err != nil && os.IsNotExist(err) {
					// new record
					printRecord(str, true)
					return
				}
				if err != nil {
					return
				}

				buf, err := ioutil.ReadFile(fpath)
				if err != nil {
					return
				}

				txt := string(buf[:])
				if strings.HasSuffix(txt, "\n") {
					tlen := len(txt)
					txt = txt[:tlen-1]
				}

				// check for optional -ignore argument
				if ignr != "" {

					// ignore differences inside specified object
					ltag := "<" + ignr + ">"
					sleft, _ := SplitInTwoAt(str, ltag, LEFT)
					tleft, _ := SplitInTwoAt(txt, ltag, LEFT)

					rtag := "</" + ignr + ">"
					_, srght := SplitInTwoAt(str, rtag, RIGHT)
					_, trght := SplitInTwoAt(txt, rtag, RIGHT)

					if sleft == tleft && srght == trght {
						if doReport {
							os.Stdout.WriteString("NO ")
							os.Stdout.WriteString(id)
							os.Stdout.WriteString("\n")
						}
						return
					}

				} else {

					// compare entirety of objects
					if str == txt {
						if doReport {
							os.Stdout.WriteString("NO ")
							os.Stdout.WriteString(id)
							os.Stdout.WriteString("\n")
						}
						return
					}
				}

				// substantively modified record
				printRecord(str, false)
			})

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		if timr {
			printDuration("records")
		}

		return
	}

	// SAVE XML COMPONENT RECORDS TO LOCAL DIRECTORY INDEXED BY TRIE ON IDENTIFIER

	// -archive plus -index plus -pattern saves XML files in trie-based directory structure
	if stsh != "" && indx != "" {

		xmlq := CreateProducer(topPattern, star, rdr)
		stsq := CreateStashers(stsh, parent, indx, hshv, zipp, xmlq)

		if xmlq == nil || stsq == nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to create stash generator\n")
			os.Exit(1)
		}

		// drain output channel
		for str := range stsq {

			if hshv {
				// print table of UIDs and hash values
				os.Stdout.WriteString(str)
			}

			recordCount++
			runtime.Gosched()
		}

		debug.FreeOSMemory()

		if timr {
			printDuration("records")
		}

		return
	}

	// READ FILE OF IDENTIFIERS AND EXTRACT SELECTED RECORDS FROM XML INPUT FILE

	// -index plus -unique [plus -head/-tail/-hd/-tl] plus -pattern with no other extraction arguments
	// takes an XML input file and a file of its UIDs and keeps only the last version of each record
	if indx != "" && unqe != "" && len(args) == 2 {

		// read file of identifiers to use for filtering
		fl, err := os.Open(unqe)
		if err != nil {
			fmt.Fprintf(os.Stderr, "\nERROR: Unable to open identifier file '%s'\n", unqe)
			os.Exit(1)
		}

		// create map that counts instances of each UID
		order := make(map[string]int)

		scanr := bufio.NewScanner(fl)

		// read lines of identifiers
		for scanr.Scan() {

			id := scanr.Text()

			// map records count for given identifier
			val := order[id]
			val++
			order[id] = val
		}

		fl.Close()

		find := ParseIndex(indx)

		if head != "" {
			os.Stdout.WriteString(head)
			os.Stdout.WriteString("\n")
		}

		PartitionPattern(topPattern, star, rdr,
			func(rec int, ofs int64, str string) {
				recordCount++

				id := FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}

				val, ok := order[id]
				if !ok {
					// not in identifier list, skip
					return
				}
				// decrement count in map
				val--
				order[id] = val
				if val > 0 {
					// only write last record with a given identifier
					return
				}

				if hd != "" {
					os.Stdout.WriteString(hd)
					os.Stdout.WriteString("\n")
				}

				// write selected record
				os.Stdout.WriteString(str[:])
				os.Stdout.WriteString("\n")

				if tl != "" {
					os.Stdout.WriteString(tl)
					os.Stdout.WriteString("\n")
				}
			})

		if tail != "" {
			os.Stdout.WriteString(tail)
			os.Stdout.WriteString("\n")
		}

		if timr {
			printDuration("records")
		}

		return
	}

	// GENERATE RECORD INDEX ON XML INPUT FILE

	// -index plus -pattern prints record identifier, file offset, and XML size
	if indx != "" {

		lbl := ""
		// check for optional filename label after -pattern argument (undocumented)
		if len(args) > 3 && args[2] == "-lbl" {
			lbl = args[3]

			lbl = strings.TrimSpace(lbl)
			if strings.HasPrefix(lbl, "pubmed") {
				lbl = lbl[7:]
			}
			if strings.HasSuffix(lbl, ".xml.gz") {
				xlen := len(lbl)
				lbl = lbl[:xlen-7]
			}
			lbl = strings.TrimSpace(lbl)
		}

		// legend := "ID\tREC\tOFST\tSIZE"

		find := ParseIndex(indx)

		PartitionPattern(topPattern, star, rdr,
			func(rec int, ofs int64, str string) {
				recordCount++

				id := FindIdentifier(str[:], parent, find)
				if id == "" {
					return
				}
				if lbl != "" {
					fmt.Printf("%s\t%d\t%d\t%d\t%s\n", id, rec, ofs, len(str), lbl)
				} else {
					fmt.Printf("%s\t%d\t%d\t%d\n", id, rec, ofs, len(str))
				}
			})

		if timr {
			printDuration("records")
		}

		return
	}

	// REPORT UNRECOGNIZED COMMAND

	fmt.Fprintf(os.Stderr, "\nERROR: Unrecognized rchive command\n")
	os.Exit(1)
}
