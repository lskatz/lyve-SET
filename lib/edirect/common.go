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
// File Name:  common.go
//
// Author:  Jonathan Kans
//
// ==========================================================================

/*
  Download external Go libraries by running:

  cd "$GOPATH"
  go get -u github.com/antzucaro/matchr
  go get -u github.com/fiam/gounidecode/unidecode
  go get -u github.com/surgebase/porter2
  go get -u golang.org/x/text/runes
  go get -u golang.org/x/text/transform
  go get -u golang.org/x/text/unicode/norm

  Test for presence of Go compiler, and cross-compile xtract executable, by running:

  if hash go 2>/dev/null
  then
    osname=`uname -s`
    cputype=`uname -m`
    case "$osname-$cputype" in
      Linux-x86_64 )
        platform=Linux
        goos=linux
        goarch=amd64
        ;;
      Darwin-x86_64 )
        platform=Darwin
        goos=darwin
        goarch=amd64
        ;;
      CYGWIN_NT-* | MINGW*-* )
        platform=CYGWIN_NT
        goos=windows
        goarch=386
        ;;
      Linux-*arm* )
        platform=ARM
        goos=linux
        goarch=arm
        ;;
    esac
    if [ -n "$platform" ]
    then
      env GOOS="$goos" GOARCH="$goarch" go build -o xtract."$platform" xtract.go common.go
      env GOOS="$goos" GOARCH="$goarch" go build -o rchive."$platform" rchive.go common.go
    fi
  fi
*/

package main

import (
	"container/heap"
	"fmt"
	"github.com/surgebase/porter2"
	"golang.org/x/text/runes"
	"golang.org/x/text/transform"
	"golang.org/x/text/unicode/norm"
	"html"
	"io"
	"os"
	"strconv"
	"strings"
	"sync"
	"unicode"
)

// TYPED CONSTANTS

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
	DOCTYPETAG
	OBJECTTAG
	CONTAINERTAG
	ISCLOSED
	BADTAG
)

type ContentType int

const (
	NONE  ContentType = iota
	PLAIN ContentType = 1 << iota
	MIXED
	AMPER
	ASCII
)

type VerifyType int

const (
	_ VerifyType = iota
	START
	STOP
	CHAR
	OTHER
)

// ARGUMENT MAPS

var scriptRunes = map[rune]rune{
	'\u00B2': '2',
	'\u00B3': '3',
	'\u00B9': '1',
	'\u2070': '0',
	'\u2071': '1',
	'\u2074': '4',
	'\u2075': '5',
	'\u2076': '6',
	'\u2077': '7',
	'\u2078': '8',
	'\u2079': '9',
	'\u207A': '+',
	'\u207B': '-',
	'\u207C': '=',
	'\u207D': '(',
	'\u207E': ')',
	'\u207F': 'n',
	'\u2080': '0',
	'\u2081': '1',
	'\u2082': '2',
	'\u2083': '3',
	'\u2084': '4',
	'\u2085': '5',
	'\u2086': '6',
	'\u2087': '7',
	'\u2088': '8',
	'\u2089': '9',
	'\u208A': '+',
	'\u208B': '-',
	'\u208C': '=',
	'\u208D': '(',
	'\u208E': ')',
}

var accentRunes = map[rune]rune{
	'\u00D8': 'O',
	'\u00F0': 'd',
	'\u00F8': 'o',
	'\u0111': 'd',
	'\u0131': 'i',
	'\u0141': 'L',
	'\u0142': 'l',
	'\u02BC': '\'',
}

var ligatureRunes = map[rune]string{
	'\u00DF': "ss",
	'\u00E6': "ae",
	'\uFB00': "ff",
	'\uFB01': "fi",
	'\uFB02': "fl",
	'\uFB03': "ffi",
	'\uFB04': "ffl",
	'\uFB05': "ft",
	'\uFB06': "st",
}

var greekRunes = map[rune]string{
	'\u0190': "epsilon",
	'\u025B': "epsilon",
	'\u03B1': "alpha",
	'\u03B2': "beta",
	'\u03B3': "gamma",
	'\u03B4': "delta",
	'\u03B5': "epsilon",
	'\u03B6': "zeta",
	'\u03B7': "eta",
	'\u03B8': "theta",
	'\u03B9': "iota",
	'\u03BA': "kappa",
	'\u03BB': "lambda",
	'\u03BC': "mu",
	'\u03BD': "nu",
	'\u03BE': "xi",
	'\u03BF': "omicron",
	'\u03C0': "pi",
	'\u03C1': "rho",
	'\u03C3': "sigma",
	'\u03C4': "tau",
	'\u03C5': "upsilon",
	'\u03C6': "phi",
	'\u03C7': "chi",
	'\u03C8': "psi",
	'\u03C9': "omega",
	'\u0391': "alpha",
	'\u0392': "beta",
	'\u0393': "gamma",
	'\u0394': "delta",
	'\u0395': "epsilon",
	'\u0396': "zeta",
	'\u0397': "eta",
	'\u0398': "theta",
	'\u0399': "iota",
	'\u039A': "kappa",
	'\u039B': "lambda",
	'\u039C': "mu",
	'\u039D': "nu",
	'\u039E': "xi",
	'\u039F': "omicron",
	'\u03A0': "pi",
	'\u03A1': "rho",
	'\u03A3': "sigma",
	'\u03A4': "tau",
	'\u03A5': "upsilon",
	'\u03A6': "phi",
	'\u03A7': "chi",
	'\u03A8': "psi",
	'\u03A9': "omega",
	'\u03D1': "theta",
	'\u03D5': "phi",
	'\u03D6': "pi",
	'\u03F0': "kappa",
	'\u03F1': "rho",
	'\u03F5': "epsilon",
}

var isStopWord = map[string]bool{
	"!":             true,
	"\"":            true,
	"#":             true,
	"$":             true,
	"%":             true,
	"&":             true,
	"'":             true,
	"(":             true,
	")":             true,
	"*":             true,
	"+":             true,
	",":             true,
	"-":             true,
	".":             true,
	"/":             true,
	":":             true,
	";":             true,
	"<":             true,
	"=":             true,
	">":             true,
	"?":             true,
	"@":             true,
	"[":             true,
	"\\":            true,
	"]":             true,
	"^":             true,
	"_":             true,
	"`":             true,
	"{":             true,
	"|":             true,
	"}":             true,
	"~":             true,
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

var htmlRepair = map[string]string{
	"&amp;lt;b&amp;gt;":     "<b>",
	"&amp;lt;i&amp;gt;":     "<i>",
	"&amp;lt;u&amp;gt;":     "<u>",
	"&amp;lt;/b&amp;gt;":    "</b>",
	"&amp;lt;/i&amp;gt;":    "</i>",
	"&amp;lt;/u&amp;gt;":    "</u>",
	"&amp;lt;b/&amp;gt;":    "<b/>",
	"&amp;lt;i/&amp;gt;":    "<i/>",
	"&amp;lt;u/&amp;gt;":    "<u/>",
	"&amp;lt;b /&amp;gt;":   "<b/>",
	"&amp;lt;i /&amp;gt;":   "<i/>",
	"&amp;lt;u /&amp;gt;":   "<u/>",
	"&amp;lt;sub&amp;gt;":   "<sub>",
	"&amp;lt;sup&amp;gt;":   "<sup>",
	"&amp;lt;/sub&amp;gt;":  "</sub>",
	"&amp;lt;/sup&amp;gt;":  "</sup>",
	"&amp;lt;sub/&amp;gt;":  "<sub/>",
	"&amp;lt;sup/&amp;gt;":  "<sup/>",
	"&amp;lt;sub /&amp;gt;": "<sub/>",
	"&amp;lt;sup /&amp;gt;": "<sup/>",
	"&lt;b&gt;":             "<b>",
	"&lt;i&gt;":             "<i>",
	"&lt;u&gt;":             "<u>",
	"&lt;/b&gt;":            "</b>",
	"&lt;/i&gt;":            "</i>",
	"&lt;/u&gt;":            "</u>",
	"&lt;b/&gt;":            "<b/>",
	"&lt;i/&gt;":            "<i/>",
	"&lt;u/&gt;":            "<u/>",
	"&lt;b /&gt;":           "<b/>",
	"&lt;i /&gt;":           "<i/>",
	"&lt;u /&gt;":           "<u/>",
	"&lt;sub&gt;":           "<sub>",
	"&lt;sup&gt;":           "<sup>",
	"&lt;/sub&gt;":          "</sub>",
	"&lt;/sup&gt;":          "</sup>",
	"&lt;sub/&gt;":          "<sub/>",
	"&lt;sup/&gt;":          "<sup/>",
	"&lt;sub /&gt;":         "<sub/>",
	"&lt;sup /&gt;":         "<sup/>",
}

var hyphenatedPrefixes = map[string]bool{
	"anti":    true,
	"bi":      true,
	"co":      true,
	"contra":  true,
	"counter": true,
	"de":      true,
	"di":      true,
	"extra":   true,
	"infra":   true,
	"inter":   true,
	"intra":   true,
	"micro":   true,
	"mid":     true,
	"mono":    true,
	"multi":   true,
	"non":     true,
	"over":    true,
	"peri":    true,
	"post":    true,
	"pre":     true,
	"pro":     true,
	"proto":   true,
	"pseudo":  true,
	"re":      true,
	"semi":    true,
	"sub":     true,
	"super":   true,
	"supra":   true,
	"tetra":   true,
	"trans":   true,
	"tri":     true,
	"ultra":   true,
	"un":      true,
	"under":   true,
	"whole":   true,
}

// DATA OBJECTS

type Node struct {
	Name       string
	Parent     string
	Contents   string
	Attributes string
	Attribs    []string
	Children   *Node
	Next       *Node
}

type Find struct {
	Index  string
	Parent string
	Match  string
	Attrib string
	Versn  string
}

type Token struct {
	Tag   TagType
	Cont  ContentType
	Name  string
	Attr  string
	Index int
	Line  int
}

type MarkupType int

const (
	NOSCRIPT MarkupType = iota
	SUPSCRIPT
	SUBSCRIPT
	PLAINDIGIT
)

type MarkupPolicy int

const (
	NOMARKUP MarkupPolicy = iota
	FUSE
	SPACE
	PERIOD
	BRACKETS
	MARKDOWN
	SLASH
	TAGS
	TERSE
)

// GLOBAL VARIABLES

var (
	InBlank   [256]bool
	InFirst   [256]bool
	InElement [256]bool
	InLower   [256]bool
	InContent [256]bool

	ChanDepth int
	FarmSize  int
	HeapSize  int
	NumServe  int

	DoCompress  bool
	DoCleanup   bool
	DoStrict    bool
	DoMixed     bool
	DoUnicode   bool
	DoScript    bool
	DoMathML    bool
	DeAccent    bool
	DoASCII     bool
	DoStem      bool
	DeStop      bool
	AllowEmbed  bool
	ContentMods bool
	CountLines  bool

	UnicodeFix = NOMARKUP
	ScriptFix  = NOMARKUP
	MathMLFix  = NOMARKUP
)

// UTILITIES

func IsNotJustWhitespace(str string) bool {

	for _, ch := range str {
		if ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r' && ch != '\f' {
			return true
		}
	}

	return false
}

func IsNotASCII(str string) bool {

	for _, ch := range str {
		if ch > 127 {
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

	for _, ch := range str {
		if !unicode.IsUpper(ch) && !unicode.IsDigit(ch) {
			return false
		}
	}

	return true
}

func IsAllDigitsOrPeriod(str string) bool {

	for _, ch := range str {
		if !unicode.IsDigit(ch) && ch != '.' {
			return false
		}
	}

	return true
}

func IsAllNumeric(str string) bool {

	for _, ch := range str {
		if !unicode.IsDigit(ch) &&
			ch != '.' &&
			ch != '+' &&
			ch != '-' &&
			ch != '*' &&
			ch != '/' &&
			ch != ',' &&
			ch != '$' &&
			ch != '#' &&
			ch != '%' &&
			ch != '(' &&
			ch != ')' {
			return false
		}
	}

	return true
}

func HasAngleBracket(str string) bool {

	hasAmp := false
	hasSemi := false

	for _, ch := range str {
		if ch == '<' || ch == '>' {
			return true
		} else if ch == '&' {
			hasAmp = true
		} else if ch == ';' {
			hasSemi = true
		}
	}

	if hasAmp && hasSemi {
		if strings.Contains(str, "&lt;") ||
			strings.Contains(str, "&gt;") ||
			strings.Contains(str, "&amp;") {
			return true
		}
	}

	return false
}

func HasAdjacentSpaces(str string) bool {

	whiteSpace := false

	for _, ch := range str {
		if ch == ' ' || ch == '\n' {
			if whiteSpace {
				return true
			}
			whiteSpace = true
		} else {
			whiteSpace = false
		}
	}

	return false
}

func HasAdjacentSpacesOrNewline(str string) bool {

	whiteSpace := false

	for _, ch := range str {
		if ch == '\n' {
			return true
		}
		if ch == ' ' {
			if whiteSpace {
				return true
			}
			whiteSpace = true
		} else {
			whiteSpace = false
		}
	}

	return false
}

func CompressRunsOfSpaces(str string) string {

	whiteSpace := false
	var buffer strings.Builder

	for _, ch := range str {
		if ch < 127 && InBlank[ch] {
			if !whiteSpace {
				buffer.WriteRune(' ')
			}
			whiteSpace = true
		} else {
			buffer.WriteRune(ch)
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
	if ch < 127 && InBlank[ch] {
		return true
	}

	strlen := len(str)
	ch = str[strlen-1]
	if ch < 127 && InBlank[ch] {
		return true
	}

	return false
}

func HasBadSpace(str string) bool {

	for _, ch := range str {
		if ch > 127 && unicode.IsSpace(ch) {
			return true
		}
	}

	return false
}

func CleanupBadSpaces(str string) string {

	var buffer strings.Builder

	for _, ch := range str {
		if ch > 127 && unicode.IsSpace(ch) {
			buffer.WriteRune(' ')
		} else {
			buffer.WriteRune(ch)
		}
	}

	return buffer.String()
}

func RepairEncodedMarkup(str string) string {

	var buffer strings.Builder

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

	skip := 0

	for i, ch := range str {
		if skip > 0 {
			skip--
			continue
		}
		if ch == '<' {
			// remove internal tags in runs of subscripts or superscripts
			if strings.HasPrefix(str[i:], "</sub><sub>") || strings.HasPrefix(str[i:], "</sup><sup>") {
				skip = 10
				continue
			}
			buffer.WriteRune(ch)
			continue
		} else if ch != '&' {
			buffer.WriteRune(ch)
			continue
		} else if strings.HasPrefix(str[i:], "&lt;") {
			sub := lookAhead(str[i:], 14)
			txt, ok := htmlRepair[sub]
			if ok {
				adv := len(sub) - 1
				// do not convert if flanked by spaces - it may be a scientific symbol,
				// e.g., fragments <i> in PMID 9698410, or escaped <b> and <d> tags used
				// to indicate stem position in letters in PMID 21892341
				if i < 1 || str[i-1] != ' ' || !strings.HasPrefix(str[i+adv:], "; ") {
					buffer.WriteString(txt)
					skip = adv
					continue
				}
			}
		} else if strings.HasPrefix(str[i:], "&amp;") {
			if strings.HasPrefix(str[i:], "&amp;lt;") {
				sub := lookAhead(str[i:], 22)
				txt, ok := htmlRepair[sub]
				if ok {
					buffer.WriteString(txt)
					skip = len(sub) - 1
					continue
				} else {
					buffer.WriteString("&lt;")
					skip = 7
					continue
				}
			} else if strings.HasPrefix(str[i:], "&amp;gt;") {
				buffer.WriteString("&gt;")
				skip = 7
				continue
			} else {
				skip = 4
				j := i + 5
				// remove runs of multiply-encoded ampersands
				for strings.HasPrefix(str[j:], "amp;") {
					skip += 4
					j += 4
				}
				// then look for special symbols used in PubMed records
				if strings.HasPrefix(str[j:], "lt;") {
					buffer.WriteString("&lt;")
					skip += 3
				} else if strings.HasPrefix(str[j:], "gt;") {
					buffer.WriteString("&gt;")
					skip += 3
				} else if strings.HasPrefix(str[j:], "frac") {
					buffer.WriteString("&frac")
					skip += 4
				} else if strings.HasPrefix(str[j:], "plusmn") {
					buffer.WriteString("&plusmn")
					skip += 6
				} else if strings.HasPrefix(str[j:], "acute") {
					buffer.WriteString("&acute")
					skip += 5
				} else if strings.HasPrefix(str[j:], "aacute") {
					buffer.WriteString("&aacute")
					skip += 6
				} else if strings.HasPrefix(str[j:], "rsquo") {
					buffer.WriteString("&rsquo")
					skip += 5
				} else if strings.HasPrefix(str[j:], "lsquo") {
					buffer.WriteString("&lsquo")
					skip += 5
				} else if strings.HasPrefix(str[j:], "micro") {
					buffer.WriteString("&micro")
					skip += 5
				} else if strings.HasPrefix(str[j:], "oslash") {
					buffer.WriteString("&oslash")
					skip += 6
				} else if strings.HasPrefix(str[j:], "kgr") {
					buffer.WriteString("&kgr")
					skip += 3
				} else if strings.HasPrefix(str[j:], "apos") {
					buffer.WriteString("&apos")
					skip += 4
				} else if strings.HasPrefix(str[j:], "quot") {
					buffer.WriteString("&quot")
					skip += 4
				} else if strings.HasPrefix(str[j:], "alpha") {
					buffer.WriteString("&alpha")
					skip += 5
				} else if strings.HasPrefix(str[j:], "beta") {
					buffer.WriteString("&beta")
					skip += 4
				} else if strings.HasPrefix(str[j:], "gamma") {
					buffer.WriteString("&gamma")
					skip += 5
				} else if strings.HasPrefix(str[j:], "Delta") {
					buffer.WriteString("&Delta")
					skip += 5
				} else if strings.HasPrefix(str[j:], "phi") {
					buffer.WriteString("&phi")
					skip += 3
				} else if strings.HasPrefix(str[j:], "ge") {
					buffer.WriteString("&ge")
					skip += 2
				} else if strings.HasPrefix(str[j:], "sup2") {
					buffer.WriteString("&sup2")
					skip += 4
				} else if strings.HasPrefix(str[j:], "#") {
					buffer.WriteString("&")
				} else {
					buffer.WriteString("&amp;")
				}
				continue
			}
		}

		buffer.WriteRune(ch)
	}

	return buffer.String()
}

func RemoveEmbeddedMarkup(str string) string {

	inContent := true
	var buffer strings.Builder

	for _, ch := range str {
		if ch == '<' {
			inContent = false
		} else if ch == '>' {
			inContent = true
		} else if inContent {
			buffer.WriteRune(ch)
		}
	}

	return buffer.String()
}

func HTMLAhead(text string, idx int) int {

	// record position of < character
	start := idx

	// at start of element
	idx++
	ch := text[idx]

	if ch == '/' {
		// skip past end tag symbol
		idx++
		ch = text[idx]
	}

	// all embedded markup tags start with a lower-case letter
	if ch < 'a' || ch > 'z' {
		// except for DispFormula in PubmedArticle
		if ch == 'D' && strings.HasPrefix(text[idx:], "DispFormula") {
			for ch != '>' {
				idx++
				ch = text[idx]
			}
			return idx + 1 - start
		}

		// otherwise not a recognized markup tag
		return 0
	}

	idx++
	ch = text[idx]
	for InLower[ch] {
		idx++
		ch = text[idx]
	}

	// if tag name was not all lower-case, then exit
	if ch >= 'A' && ch <= 'Z' {
		return 0
	}

	// skip to end of element, past any attributes or slash character
	for ch != '>' {
		idx++
		ch = text[idx]
	}

	// return number of characters to advance to skip this markup tag
	return idx + 1 - start
}

func HTMLBehind(bufr []byte, pos int) bool {

	for pos >= 0 {
		if bufr[pos] == '<' {
			return HTMLAhead(string(bufr), pos) != 0
		}
		pos--
	}

	return false
}

func HasUnicodeMarkup(str string) bool {

	for _, ch := range str {
		if ch <= 127 {
			continue
		}
		// check for Unicode superscript or subscript characters
		if ch == '\u00B2' || ch == '\u00B3' || ch == '\u00B9' || (ch >= '\u2070' && ch <= '\u208E') {
			return true
		}
	}

	return false
}

func IsUnicodeSuper(ch rune) bool {
	return ch == '\u00B2' || ch == '\u00B3' || ch == '\u00B9' || (ch >= '\u2070' && ch <= '\u207F')
}

func IsUnicodeSubsc(ch rune) bool {
	return ch >= '\u2080' && ch <= '\u208E'
}

func RepairUnicodeMarkup(str string, policy MarkupPolicy) string {

	type MarkupType int

	const (
		NOSCRIPT MarkupType = iota
		SUPSCRIPT
		SUBSCRIPT
		PLAINDIGIT
	)

	var buffer strings.Builder

	// to improve readability, keep track of switches between numeric types, add period at transitions when converting to plain ASCII
	level := NOSCRIPT

	for _, ch := range str {
		if ch > 127 {
			if IsUnicodeSuper(ch) {
				rn, ok := scriptRunes[ch]
				if ok {
					ch = rn
					switch level {
					case NOSCRIPT:
						switch policy {
						case PERIOD:
						case SPACE:
						case BRACKETS:
							buffer.WriteRune('[')
						case MARKDOWN:
							buffer.WriteRune('^')
						case SLASH:
						case TAGS:
							buffer.WriteString("<sup>")
						}
					case SUPSCRIPT:
						switch policy {
						case PERIOD:
						case SPACE:
						case BRACKETS:
						case MARKDOWN:
						case SLASH:
						case TAGS:
						}
					case SUBSCRIPT:
						switch policy {
						case PERIOD:
							buffer.WriteRune('.')
						case SPACE:
							buffer.WriteRune(' ')
						case BRACKETS:
							buffer.WriteRune(')')
							buffer.WriteRune('[')
						case MARKDOWN:
							buffer.WriteRune('~')
							buffer.WriteRune('^')
						case SLASH:
							buffer.WriteRune('\\')
						case TAGS:
							buffer.WriteString("</sub>")
							buffer.WriteString("<sup>")
						}
					case PLAINDIGIT:
						switch policy {
						case PERIOD:
							buffer.WriteRune('.')
						case SPACE:
							buffer.WriteRune(' ')
						case BRACKETS:
							buffer.WriteRune('[')
						case MARKDOWN:
							buffer.WriteRune('^')
						case SLASH:
							buffer.WriteRune('\\')
						case TAGS:
							buffer.WriteString("<sup>")
						}
					}
					level = SUPSCRIPT
				}
			} else if IsUnicodeSubsc(ch) {
				rn, ok := scriptRunes[ch]
				if ok {
					ch = rn
					switch level {
					case NOSCRIPT:
						switch policy {
						case PERIOD:
						case SPACE:
						case BRACKETS:
							buffer.WriteRune('(')
						case MARKDOWN:
							buffer.WriteRune('~')
						case SLASH:
						case TAGS:
							buffer.WriteString("<sub>")
						}
					case SUPSCRIPT:
						switch policy {
						case PERIOD:
							buffer.WriteRune('.')
						case SPACE:
							buffer.WriteRune(' ')
						case BRACKETS:
							buffer.WriteRune(']')
							buffer.WriteRune('(')
						case MARKDOWN:
							buffer.WriteRune('^')
							buffer.WriteRune('~')
						case SLASH:
							buffer.WriteRune('/')
						case TAGS:
							buffer.WriteString("</sup>")
							buffer.WriteString("<sub>")
						}
					case SUBSCRIPT:
						switch policy {
						case PERIOD:
						case SPACE:
						case BRACKETS:
						case MARKDOWN:
						case SLASH:
						case TAGS:
						}
					case PLAINDIGIT:
						switch policy {
						case PERIOD:
							buffer.WriteRune('.')
						case SPACE:
							buffer.WriteRune(' ')
						case BRACKETS:
							buffer.WriteRune('(')
						case MARKDOWN:
							buffer.WriteRune('~')
						case SLASH:
							buffer.WriteRune('/')
						case TAGS:
							buffer.WriteString("<sub>")
						}
					}
					level = SUBSCRIPT
				}
			} else {
				level = NOSCRIPT
			}
		} else if ch >= '0' && ch <= '9' {
			switch level {
			case NOSCRIPT:
				switch policy {
				case PERIOD:
				case SPACE:
				case BRACKETS:
				case MARKDOWN:
				case SLASH:
				case TAGS:
				}
			case SUPSCRIPT:
				switch policy {
				case PERIOD:
					buffer.WriteRune('.')
				case SPACE:
					buffer.WriteRune(' ')
				case BRACKETS:
					buffer.WriteRune(']')
				case MARKDOWN:
					buffer.WriteRune('^')
				case SLASH:
					buffer.WriteRune('/')
				case TAGS:
					buffer.WriteString("</sup>")
				}
			case SUBSCRIPT:
				switch policy {
				case PERIOD:
					buffer.WriteRune('.')
				case SPACE:
					buffer.WriteRune(' ')
				case BRACKETS:
					buffer.WriteRune(')')
				case MARKDOWN:
					buffer.WriteRune('~')
				case SLASH:
					buffer.WriteRune('\\')
				case TAGS:
					buffer.WriteString("</sub>")
				}
			case PLAINDIGIT:
				switch policy {
				case PERIOD:
				case SPACE:
				case BRACKETS:
				case MARKDOWN:
				case SLASH:
				case TAGS:
				}
			}
			level = PLAINDIGIT
		} else {
			switch level {
			case NOSCRIPT:
				switch policy {
				case PERIOD:
				case SPACE:
				case BRACKETS:
				case MARKDOWN:
				case SLASH:
				case TAGS:
				}
			case SUPSCRIPT:
				switch policy {
				case PERIOD:
				case SPACE:
				case BRACKETS:
					buffer.WriteRune(']')
				case MARKDOWN:
					buffer.WriteRune('^')
				case SLASH:
				case TAGS:
					buffer.WriteString("</sup>")
				}
			case SUBSCRIPT:
				switch policy {
				case PERIOD:
				case SPACE:
				case BRACKETS:
					buffer.WriteRune(')')
				case MARKDOWN:
					buffer.WriteRune('~')
				case SLASH:
				case TAGS:
					buffer.WriteString("</sub>")
				}
			case PLAINDIGIT:
				switch policy {
				case PERIOD:
				case SPACE:
				case BRACKETS:
				case MARKDOWN:
				case SLASH:
				case TAGS:
				}
			}
			level = NOSCRIPT
		}
		buffer.WriteRune(ch)
	}

	switch level {
	case NOSCRIPT:
		switch policy {
		case PERIOD:
		case SPACE:
		case BRACKETS:
		case MARKDOWN:
		case SLASH:
		case TAGS:
		}
	case SUPSCRIPT:
		switch policy {
		case PERIOD:
		case SPACE:
		case BRACKETS:
			buffer.WriteRune(']')
		case MARKDOWN:
			buffer.WriteRune('^')
		case SLASH:
		case TAGS:
			buffer.WriteString("</sup>")
		}
	case SUBSCRIPT:
		switch policy {
		case PERIOD:
		case SPACE:
		case BRACKETS:
			buffer.WriteRune(')')
		case MARKDOWN:
			buffer.WriteRune('~')
		case SLASH:
		case TAGS:
			buffer.WriteString("</sub>")
		}
	case PLAINDIGIT:
		switch policy {
		case PERIOD:
		case SPACE:
		case BRACKETS:
		case MARKDOWN:
		case SLASH:
		case TAGS:
		}
	}

	return buffer.String()
}

func FlattenMathML(str string, policy MarkupPolicy) string {

	findNextXMLBlock := func(txt string) (int, int, bool) {

		beg := strings.Index(txt, "<")
		if beg < 0 {
			return -1, -1, false
		}
		end := strings.Index(txt, ">")
		if end < 0 {
			return -1, -1, false
		}
		end++
		return beg, end, true
	}

	var arry []string

	for {
		beg, end, ok := findNextXMLBlock(str)
		if !ok {
			break
		}
		pfx := str[:beg]
		pfx = strings.TrimSpace(pfx)
		if pfx != "" {
			arry = append(arry, pfx)
		}
		tmp := str[beg:end]
		tmp = strings.TrimSpace(tmp)
		str = str[end:]
	}

	switch policy {
	case PERIOD:
	case SPACE:
		str = strings.Join(arry, " ")
	case BRACKETS:
	case MARKDOWN:
	case SLASH:
	case TAGS:
	case TERSE:
		str = strings.Join(arry, "")
	}

	str = strings.TrimSpace(str)

	// str = RemoveEmbeddedMarkup(str)

	return str
}

func RepairMathMLMarkup(str string, policy MarkupPolicy) string {

	str = strings.Replace(str, "> <mml:", "><mml:", -1)
	str = strings.Replace(str, "> </mml:", "></mml:", -1)

	findNextMathBlock := func(txt string) (int, int, bool) {

		beg := strings.Index(txt, "<DispFormula")
		if beg < 0 {
			return -1, -1, false
		}
		end := strings.Index(txt, "</DispFormula>")
		if end < 0 {
			return -1, -1, false
		}
		end += 14
		return beg, end, true
	}

	var arry []string

	for {
		beg, end, ok := findNextMathBlock(str)
		if !ok {
			break
		}
		pfx := str[:beg]
		pfx = strings.TrimSpace(pfx)
		arry = append(arry, pfx)
		tmp := str[beg:end]
		if strings.HasPrefix(tmp, "<DispFormula") {
			tmp = FlattenMathML(tmp, policy)
		}
		tmp = strings.TrimSpace(tmp)
		arry = append(arry, tmp)
		str = str[end:]
	}

	str = strings.TrimSpace(str)
	arry = append(arry, str)

	return strings.Join(arry, " ")
}

func RepairScriptMarkup(str string, policy MarkupPolicy) string {

	var buffer strings.Builder

	skip := 0

	for i, ch := range str {
		if skip > 0 {
			skip--
			continue
		}
		if ch == '<' {
			if strings.HasPrefix(str[i:], "<sub>") {
				switch policy {
				case PERIOD:
				case SPACE:
				case BRACKETS:
					buffer.WriteRune('(')
				case MARKDOWN:
					buffer.WriteRune('~')
				}
				skip = 4
				continue
			}
			if strings.HasPrefix(str[i:], "<sup>") {
				switch policy {
				case PERIOD:
				case SPACE:
				case BRACKETS:
					buffer.WriteRune('[')
				case MARKDOWN:
					buffer.WriteRune('^')
				}
				skip = 4
				continue
			}
			if strings.HasPrefix(str[i:], "</sub>") {
				if strings.HasPrefix(str[i+6:], "<sup>") {
					switch policy {
					case PERIOD:
						buffer.WriteRune('.')
					case SPACE:
						buffer.WriteRune(' ')
					case BRACKETS:
						buffer.WriteRune(')')
						buffer.WriteRune('[')
					case MARKDOWN:
						buffer.WriteRune('~')
						buffer.WriteRune('^')
					}
					skip = 10
					continue
				}
				switch policy {
				case PERIOD:
				case SPACE:
				case BRACKETS:
					buffer.WriteRune(')')
				case MARKDOWN:
					buffer.WriteRune('~')
				}
				skip = 5
				continue
			}
			if strings.HasPrefix(str[i:], "</sup>") {
				if strings.HasPrefix(str[i+6:], "<sub>") {
					switch policy {
					case PERIOD:
						buffer.WriteRune('.')
					case SPACE:
						buffer.WriteRune(' ')
					case BRACKETS:
						buffer.WriteRune(']')
						buffer.WriteRune('(')
					case MARKDOWN:
						buffer.WriteRune('^')
						buffer.WriteRune('~')
					}
					skip = 10
					continue
				}
				switch policy {
				case PERIOD:
				case SPACE:
				case BRACKETS:
					buffer.WriteRune(']')
				case MARKDOWN:
					buffer.WriteRune('^')
				}
				skip = 5
				continue
			}
		}

		buffer.WriteRune(ch)
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
	for _, ch := range str {
		if isSlash {
			switch ch {
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
				res[idx] = byte(ch)
			}
			idx++
			isSlash = false
		} else if ch == '\\' {
			isSlash = true
		} else {
			res[idx] = byte(ch)
			idx++
		}
	}

	res = res[0:idx]

	return string(res)
}

func DoTrimFlankingHTML(str string) string {

	badPrefix := [10]string{
		"<i></i>",
		"<b></b>",
		"<u></u>",
		"<sup></sup>",
		"<sub></sub>",
		"</i>",
		"</b>",
		"</u>",
		"</sup>",
		"</sub>",
	}

	badSuffix := [10]string{
		"<i></i>",
		"<b></b>",
		"<u></u>",
		"<sup></sup>",
		"<sub></sub>",
		"<i>",
		"<b>",
		"<u>",
		"<sup>",
		"<sub>",
	}

	if strings.Contains(str, "<") {
		goOn := true
		for goOn {
			goOn = false
			for _, tag := range badPrefix {
				if strings.HasPrefix(str, tag) {
					str = str[len(tag):]
					goOn = true
				}
			}
			for _, tag := range badSuffix {
				if strings.HasSuffix(str, tag) {
					str = str[:len(str)-len(tag)]
					goOn = true
				}
			}
		}
	}

	return str
}

func HasBadAccent(str string) bool {

	for _, ch := range str {
		if ch <= 127 {
			continue
		}
		// quick min-to-max check for additional characters to treat as accents
		if ch >= '\u00D8' && ch <= '\u02BC' {
			return true
		} else if ch >= '\uFB00' && ch <= '\uFB06' {
			return true
		}
	}

	return false
}

func FixBadAccent(str string) string {

	var buffer strings.Builder

	for _, ch := range str {
		if ch > 127 {
			if ch >= '\u00D8' && ch <= '\u02BC' {
				rn, ok := accentRunes[ch]
				if ok {
					buffer.WriteRune(rn)
					continue
				}
				st, ok := ligatureRunes[ch]
				if ok {
					buffer.WriteString(st)
					continue
				}
			}
			if ch >= '\uFB00' && ch <= '\uFB06' {
				st, ok := ligatureRunes[ch]
				if ok {
					buffer.WriteString(st)
					continue
				}
			}
		}
		buffer.WriteRune(ch)
	}

	return buffer.String()
}

var (
	tlock sync.Mutex
	tform transform.Transformer
)

func DoAccentTransform(str string) string {

	// transformer not reentrant, protected by mutex
	tlock.Lock()

	defer tlock.Unlock()

	if tform == nil {
		tform = transform.Chain(norm.NFD, runes.Remove(runes.In(unicode.Mn)), norm.NFC)
	}

	if tform != nil {

		var arry []string

		// split long string into words to avoid transform short internal buffer error
		terms := strings.Fields(str)

		for _, item := range terms {

			// remove accents from single word
			tmp, _, err := transform.String(tform, item)
			if err == nil {
				// collect transformed result
				arry = append(arry, tmp)
			} else {
				fmt.Fprintf(os.Stderr, "%s\n", err.Error())
			}
		}

		// reconstruct string from transformed words
		str = strings.Join(arry, " ")
	}

	// look for characters not in current external runes conversion table
	if HasBadAccent(str) {
		str = FixBadAccent(str)
	}

	return str
}

func UnicodeToASCII(str string) string {

	var buffer strings.Builder

	for _, ch := range str {
		if ch > 127 {
			s := strconv.QuoteToASCII(string(ch))
			s = strings.ToUpper(s[3:7])
			for {
				if !strings.HasPrefix(s, "0") {
					break
				}
				s = s[1:]
			}
			buffer.WriteString("&#x")
			buffer.WriteString(s)
			buffer.WriteRune(';')
			continue
		}
		buffer.WriteRune(ch)
	}

	return buffer.String()
}

func HasGreek(str string) bool {

	for _, ch := range str {
		if ch <= 127 {
			continue
		}
		// quick min-to-max check for Greek characters to convert to english words
		if ch >= '\u03B1' && ch <= '\u03C9' {
			return true
		} else if ch >= '\u0391' && ch <= '\u03A9' {
			return true
		} else if ch >= '\u03D1' && ch <= '\u03D6' {
			return true
		} else if ch >= '\u03F0' && ch <= '\u03F5' {
			return true
		} else if ch == '\u0190' || ch == '\u025B' {
			return true
		}
	}

	return false
}

func SpellGreek(str string) string {

	var buffer strings.Builder

	for _, ch := range str {
		st := ""
		ok := false
		if ch > 127 {
			if (ch >= '\u03B1' && ch <= '\u03C9') || (ch >= '\u0391' && ch <= '\u03A9') {
				st, ok = greekRunes[ch]
			} else if (ch >= '\u03D1' && ch <= '\u03D6') || (ch >= '\u03F0' && ch <= '\u03F5') {
				// glyph variants of Greek letters
				st, ok = greekRunes[ch]
			} else if ch == '\u0190' || ch == '\u025B' {
				// treat latin letter open E as epsilon
				st, ok = greekRunes[ch]
			}
		}
		if ok {
			buffer.WriteString(" ")
			buffer.WriteString(st)
			buffer.WriteString(" ")
			continue
		}
		buffer.WriteRune(ch)
	}

	return buffer.String()
}

func RemoveHyphenFromPrefix(str string) string {

	var buffer strings.Builder

	for i, ch := range str {
		if ch == '-' {
			_, ok := hyphenatedPrefixes[str[0:i]]
			if ok {
				continue
			}
		}
		buffer.WriteRune(ch)
	}

	return buffer.String()
}

func RemoveAllPrefixHyphens(str string) string {

	var arry []string

	terms := strings.Fields(str)

	for _, item := range terms {

		item = RemoveHyphenFromPrefix(item)
		arry = append(arry, item)
	}

	// reconstruct string from transformed words
	str = strings.Join(arry, " ")

	return str
}

var plock sync.RWMutex

func IsStopWord(str string) bool {

	plock.RLock()
	isSW := isStopWord[str]
	plock.RUnlock()

	return isSW
}

func ParseIndex(indx string) *Find {

	if indx == "" {
		return &Find{}
	}

	// parse parent/element@attribute^version index
	prnt, match := SplitInTwoAt(indx, "/", RIGHT)
	match, versn := SplitInTwoAt(match, "^", LEFT)
	match, attrib := SplitInTwoAt(match, "@", LEFT)

	return &Find{Index: indx, Parent: prnt, Match: match, Attrib: attrib, Versn: versn}
}

func CleanupContents(str string, ascii, amper, mixed bool) string {

	if DoCompress {
		if !AllowEmbed {
			if ascii && HasBadSpace(str) {
				str = CleanupBadSpaces(str)
			}
		}
		if HasAdjacentSpacesOrNewline(str) {
			str = CompressRunsOfSpaces(str)
		}
	}
	if DoUnicode {
		if ascii && HasUnicodeMarkup(str) {
			str = RepairUnicodeMarkup(str, UnicodeFix)
		}
	}
	if AllowEmbed {
		if amper {
			str = RepairEncodedMarkup(str)
		}
	}
	if DoScript {
		if mixed && HasAngleBracket(str) {
			str = RepairScriptMarkup(str, ScriptFix)
		}
	}
	if DoMathML {
		if mixed && HasAngleBracket(str) {
			str = RepairMathMLMarkup(str, MathMLFix)
		}
	}
	if DoStrict {
		if mixed || amper {
			if HasAngleBracket(str) {
				str = RemoveEmbeddedMarkup(str)
			}
		}
		if ascii && HasBadSpace(str) {
			str = CleanupBadSpaces(str)
		}
		if HasAdjacentSpaces(str) {
			str = CompressRunsOfSpaces(str)
		}
	}
	if DoMixed {
		if mixed {
			str = DoTrimFlankingHTML(str)
		}
		if ascii && HasBadSpace(str) {
			str = CleanupBadSpaces(str)
		}
		if HasAdjacentSpaces(str) {
			str = CompressRunsOfSpaces(str)
		}
	}
	if DeAccent {
		if ascii {
			str = DoAccentTransform(str)
		}
	}
	if DoASCII {
		if ascii {
			str = UnicodeToASCII(str)
		}
	}

	if HasFlankingSpace(str) {
		str = strings.TrimSpace(str)
	}

	return str
}

var (
	rlock sync.Mutex
	qfix  *strings.Replacer
)

func ProtectSpecialTerms(str string) string {

	// NewReplacer not reentrant (?), protected by mutex
	rlock.Lock()

	defer rlock.Unlock()

	if qfix == nil {
		// handles biologically important terms that would otherwise be dropped from index
		qfix = strings.NewReplacer(
			" 5' ", " 5_prime ",
			" 5'-", " 5_prime ",
			" 3' ", " 3_prime ",
			" 3'-", " 3_prime ",
			" b cell ", " b_cell ",
			" b-cell ", " b_cell ",
			" t cell ", " t_cell ",
			" t-cell ", " t_cell ",
			" il 1", " il_1",
			" il 2", " il_2",
			" il 3", " il_3",
			" il 4", " il_4",
			" il 5", " il_5",
			" il 6", " il_6",
			" il 7", " il_7",
			" il 8", " il_8",
			" il 9", " il_9",
			" il-1", " il_1",
			" il-2", " il_2",
			" il-3", " il_3",
			" il-4", " il_4",
			" il-5", " il_5",
			" il-6", " il_6",
			" il-7", " il_7",
			" il-8", " il_8",
			" il-9", " il_9",
			" interleukin 1", " interleukin_1",
			" interleukin 2", " interleukin_2",
			" interleukin 3", " interleukin_3",
			" interleukin 4", " interleukin_4",
			" interleukin 5", " interleukin_5",
			" interleukin 6", " interleukin_6",
			" interleukin 7", " interleukin_7",
			" interleukin 8", " interleukin_8",
			" interleukin 9", " interleukin_9",
			" interleukin-1", " interleukin_1",
			" interleukin-2", " interleukin_2",
			" interleukin-3", " interleukin_3",
			" interleukin-4", " interleukin_4",
			" interleukin-5", " interleukin_5",
			" interleukin-6", " interleukin_6",
			" interleukin-7", " interleukin_7",
			" interleukin-8", " interleukin_8",
			" interleukin-9", " interleukin_9",
		)
	}

	// must do after lower casing and removing underscores, but before removing hyphens
	if qfix != nil {

		str = CompressRunsOfSpaces(str)
		str = strings.TrimSpace(str)

		str = " " + str + " "

		str = qfix.Replace(str)

		str = CompressRunsOfSpaces(str)
		str = strings.TrimSpace(str)
	}

	return str
}

func PrepareQuery(str string) string {

	if str == "" {
		return ""
	}

	// cleanup string
	if IsNotASCII(str) {
		if HasGreek(str) {
			str = SpellGreek(str)
		}
		str = DoAccentTransform(str)
		if HasUnicodeMarkup(str) {
			str = RepairUnicodeMarkup(str, SPACE)
		}
	}

	if HasAmpOrNotASCII(str) {
		str = html.UnescapeString(str)
	}

	if HasBadSpace(str) {
		str = CleanupBadSpaces(str)
	}
	if HasAngleBracket(str) {
		str = RepairEncodedMarkup(str)
		str = RepairScriptMarkup(str, SPACE)
		str = RepairMathMLMarkup(str, SPACE)
		str = RemoveEmbeddedMarkup(str)
	}

	str = strings.Replace(str, " AND ", " & ", -1)
	str = strings.Replace(str, " OR ", " | ", -1)
	str = strings.Replace(str, " NOT ", " ! ", -1)

	str = strings.Replace(str, "(", " ( ", -1)
	str = strings.Replace(str, ")", " ) ", -1)
	str = strings.Replace(str, "&", " & ", -1)
	str = strings.Replace(str, "|", " | ", -1)
	str = strings.Replace(str, "!", " ! ", -1)

	str = strings.ToLower(str)

	str = strings.Replace(str, "_", " ", -1)

	// must do after lower casing and removing underscores, but before removing hyphens
	str = ProtectSpecialTerms(str)

	str = RemoveAllPrefixHyphens(str)

	str = strings.Replace(str, "-", " ", -1)

	// break terms at punctuation, and at non-ASCII characters, allowing Boolean control symbols, along with
	// asterisk to indicate truncation wildcard, underscore for protected terms, and tilde for proximity
	terms := strings.FieldsFunc(str, func(c rune) bool {
		return (!unicode.IsLetter(c) && !unicode.IsDigit(c) && c != '*' && c != '_' && c != '~' &&
			c != '&' && c != '|' && c != '!' && c != '(' && c != ')') || c > 127
	})

	// rejoin into processed sentence
	tmp := strings.Join(terms, " ")

	tmp = CompressRunsOfSpaces(tmp)
	tmp = strings.TrimSpace(tmp)

	return tmp
}

func PartitionQuery(str string) []string {

	if str == "" {
		return nil
	}

	str = CompressRunsOfSpaces(str)
	str = strings.TrimSpace(str)

	// fuse tildes
	str = strings.Replace(str, "~ ~", "~~", -1)
	str = strings.Replace(str, "~ ~", "~~", -1)

	str = " " + str + " "

	// flank all operators with caret
	str = strings.Replace(str, " ( ", " ^ ( ^ ", -1)
	str = strings.Replace(str, " ) ", " ^ ) ^ ", -1)
	str = strings.Replace(str, " & ", " ^ & ^ ", -1)
	str = strings.Replace(str, " | ", " ^ | ^ ", -1)
	str = strings.Replace(str, " ! ", " ^ ! ^ ", -1)

	str = CompressRunsOfSpaces(str)
	str = strings.TrimSpace(str)

	str = strings.Replace(str, "^ ^", "^", -1)

	if strings.HasPrefix(str, "^ ") {
		str = str[2:]
	}
	if strings.HasSuffix(str, " ^") {
		max := len(str)
		str = str[:max-2]
	}

	// split into non-broken phrase segments or operator symbols
	tmp := strings.Split(str, " ^ ")

	return tmp
}

func MarkStopWords(str string) string {

	if str == "" {
		return ""
	}

	var chain []string

	terms := strings.Fields(str)

	// replace unwanted and stop words with plus sign
	for _, item := range terms {

		if item == "~" {
			// allow tilde proximity indicator
			chain = append(chain, item)
			continue
		}

		if len(item) < 2 {
			// skip a single character
			chain = append(chain, "+")
			continue
		}

		if IsAllDigitsOrPeriod(item) {
			// skip terms that are all digits
			chain = append(chain, "+")
			continue
		}

		if DeStop && IsStopWord(item) {
			// skip if stop word, breaking word pair chain
			chain = append(chain, "+")
			continue
		}
		if DoStem && !strings.HasSuffix(item, "*") {
			// apply stemming algorithm
			item = porter2.Stem(item)
			item = strings.TrimSpace(item)
		}

		// record single term
		chain = append(chain, item)
	}

	// rejoin into processed sentence
	tmp := strings.Join(chain, " ")

	tmp = CompressRunsOfSpaces(tmp)
	tmp = strings.TrimSpace(tmp)

	return tmp
}

func MarkClauses(clauses []string) []string {

	var res []string

	if clauses == nil {
		return nil
	}

	for _, str := range clauses {

		// pass control symbols or angle bracket content delimiters unchanged
		if str == "(" || str == ")" || str == "&" || str == "|" || str == "!" || str == "<" || str == ">" {
			res = append(res, str)
			continue
		}

		// process clause, using plus sign to break runs of words
		tmp := MarkStopWords(str)
		res = append(res, tmp)
	}

	return res
}

// READ XML INPUT FILE INTO SET OF BLOCKS

type XMLReader struct {
	Reader    io.Reader
	Buffer    []byte
	Remainder string
	Position  int64
	Delta     int
	Closed    bool
}

func NewXMLReader(in io.Reader) *XMLReader {

	if in == nil {
		return nil
	}

	rdr := &XMLReader{Reader: in}

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
			return string(rdr.Buffer[:m]), true, false
		}

		// read next block, append behind copied remainder from previous read
		n, err := rdr.Reader.Read(rdr.Buffer[m:])
		// with data piped through stdin, read function may not always return the same number of bytes each time
		if err != nil {
			if err != io.EOF {
				// real error
				fmt.Fprintf(os.Stderr, "\nERROR: %s\n", err.Error())
				// Ignore bytes - non-conforming implementations of io.Reader may returned mangled data on non-EOF errors
				rdr.Closed = true
				return "", false, true
			}
			// end of file
			rdr.Closed = true
			if n == 0 {
				// if EOF and no more data, do not send final remainder (not terminated by right angle bracket that is used as a sentinel)
				return "", false, true
			}
		}
		if n < 0 {
			// Reality check - non-conforming implementations of io.Reader may return -1
			fmt.Fprintf(os.Stderr, "\nERROR: io.Reader returned negative count %d\n", n)
			// treat as n == 0 in order to update file offset and avoid losing previous remainder
			n = 0
		}

		// keep track of file offset
		rdr.Position += int64(rdr.Delta)
		rdr.Delta = n

		// slice of actual characters read
		bufr := rdr.Buffer[:n+m]

		// look for last > character
		// safe to back up on UTF-8 rune array when looking for 7-bit ASCII character
		pos := -1
		for pos = len(bufr) - 1; pos >= 0; pos-- {
			if bufr[pos] == '>' {
				if AllowEmbed {
					// optionally skip backwards past embedded i, b, u, sub, and sup HTML open, close, and empty tags, and MathML
					if HTMLBehind(bufr, pos) {
						continue
					}
				}
				// found end of XML tag, break
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
		var buff strings.Builder

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

	// optionally compress/cleanup tags/attributes and contents (undocumented)
	if DoCleanup {
		if HasBadSpace(line) {
			line = CleanupBadSpaces(line)
		}
		if HasAdjacentSpaces(line) {
			line = CompressRunsOfSpaces(line)
		}
	}

	return line
}

// PARSE XML BLOCK STREAM INTO STRINGS FROM <PATTERN> TO </PATTERN>

// PartitionPattern splits XML input by pattern and sends individual records to a callback
func PartitionPattern(pat, star string, rdr *XMLReader, proc func(int, int64, string)) {

	if pat == "" || rdr == nil || proc == nil {
		return
	}

	type Scanner struct {
		Pattern   string
		PatLength int
		CharSkip  [256]int
	}

	// initialize <pattern> to </pattern> scanner
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

	// check surroundings of match candidate
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

	// find next element with pattern name
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
		var accumulator strings.Builder

		match := NOPATTERN
		pos := 0
		next := 0

		offset := int64(0)

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
						offset = rdr.Position + int64(pos)
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
							proc(rec, offset, str[:])
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
	// are recursive or self-closing objects, process those through -format first

	doStar := func() {

		// current depth of -pattern objects
		level := 0

		begin := 0
		inPattern := false

		line := ""
		var accumulator strings.Builder

		match := NOPATTERN
		pos := 0
		next := 0

		offset := int64(0)

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

		// find next element in XML
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
						offset = rdr.Position + int64(pos)
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
							proc(rec, offset, str[:])
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

// PARSE ATTRIBUTES INTO TAG/VALUE PAIRS

// ParseAttributes is only run if attribute values are requested in element statements
func ParseAttributes(attrb string) []string {

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

// ParseXML calls XML parser on a partitioned string, optimized for maximum processing speed
func ParseXML(Text, parent string, tokens func(Token), find *Find) (*Node, string) {

	if Text == "" {
		return nil, ""
	}

	// node farm variables
	FarmPos := 0
	FarmMax := FarmSize
	FarmItems := make([]Node, FarmMax)

	// allocate multiple nodes in a large array for memory management efficiency
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

	// get next XML token
	nextToken := func(idx int) (TagType, ContentType, string, string, int) {

		if Text == "" {
			// signal end of XML data
			return ISCLOSED, NONE, "", "", 0
		}

		text := Text[:]
		txtlen := Txtlen

		// XML string ends with > character, acts as sentinel to check if past end of text
		if idx >= txtlen {
			// signal end of XML string
			return ISCLOSED, NONE, "", "", 0
		}

		// skip past leading blanks
		ch := text[idx]
		for InBlank[ch] {
			idx++
			ch = text[idx]
		}

		start := idx

		plainContent := true

		if AllowEmbed && ch == '<' {
			// check to see if an HTML or MathML element is at the beginning of a content string
			if HTMLAhead(text, idx) != 0 {
				plainContent = false
			}
		}

		if plainContent && ch == '<' {

			// at start of element
			idx++
			ch = text[idx]

			// check for legal first character of element
			if InFirst[ch] {

				// read element name
				start = idx
				idx++

				ch = text[idx]
				for InElement[ch] {
					idx++
					ch = text[idx]
				}

				str := text[start:idx]

				switch ch {
				case '>':
					// end of element
					idx++

					return STARTTAG, NONE, str[:], "", idx
				case '/':
					// self-closing element without attributes
					idx++
					ch = text[idx]
					if ch != '>' {
						fmt.Fprintf(os.Stderr, "\nSelf-closing element missing right angle bracket\n")
					}
					idx++

					return SELFTAG, NONE, str[:], "", idx
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
						return SELFTAG, NONE, str[:], atr[:], idx
					}
					atr := text[start:idx]
					idx++
					return STARTTAG, NONE, str[:], atr[:], idx
				default:
					fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element\n", ch)
					return STARTTAG, NONE, str[:], "", idx
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
					if InFirst[ch] {
						idx++
						ch = text[idx]
						for InElement[ch] {
							idx++
							ch = text[idx]
						}
						str := text[start:idx]
						if ch != '>' {
							fmt.Fprintf(os.Stderr, "\nUnexpected characters after end element name\n")
						}
						idx++

						return STOPTAG, NONE, str[:], "", idx
					}
					fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element\n", ch)
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
					// skip !DOCTYPE, !COMMENT, and ![CDATA[
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
						// CDATA or COMMENT block may contain internal angle brackets
						found := strings.Index(text[idx:], skipTo)
						if found < 0 {
							// string stops in middle of CDATA or COMMENT
							return ISCLOSED, NONE, "", "", idx
						}
						// adjust position past end of CDATA or COMMENT
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
					fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element\n", ch)
				}
			}

		} else if ch != '>' {

			// at start of contents
			start = idx

			hasMarkup := false
			hasNonASCII := false

			// find end of contents
			if AllowEmbed {

				for {
					for InContent[ch] {
						idx++
						ch = text[idx]
					}
					// set flags to speed up conditional content processing
					if ch == '&' {
						idx++
						ch = text[idx]
						if ch == 'a' {
							if strings.HasPrefix(text[idx:], "amp;") {
								hasMarkup = true
							}
						} else if ch == 'g' {
							if strings.HasPrefix(text[idx:], "gt;") {
								hasMarkup = true
							}
						} else if ch == 'l' {
							if strings.HasPrefix(text[idx:], "lt;") {
								hasMarkup = true
							}
						}
						continue
					}
					if ch > 127 {
						hasNonASCII = true
						idx++
						ch = text[idx]
						continue
					}
					if ch == '<' {
						// optionally allow HTML text formatting elements and super/subscripts
						advance := HTMLAhead(text, idx)
						if advance > 0 {
							idx += advance
							ch = text[idx]
							plainContent = false
							continue
						}
					}
					break
				}
			} else {
				for ch != '<' && ch != '>' {
					idx++
					ch = text[idx]
				}
			}

			// trim back past trailing blanks
			lst := idx - 1
			ch = text[lst]
			for InBlank[ch] && lst > start {
				lst--
				ch = text[lst]
			}

			str := text[start : lst+1]

			ctype := PLAIN

			if AllowEmbed {
				if !plainContent {
					ctype |= MIXED
				}
				if hasMarkup {
					ctype |= AMPER
				}
				if hasNonASCII {
					ctype |= ASCII
				}
			}

			return CONTENTTAG, ctype, str[:], "", idx
		}

		return BADTAG, NONE, "", "", idx
	}

	// Parse tokens into tree structure for exploration

	// parseLevel recursive definition
	var parseLevel func(string, string, string) (*Node, bool)

	// parse XML tags into tree structure for searching
	parseLevel = func(strt, attr, prnt string) (*Node, bool) {

		ok := true

		// obtain next node from farm
		node := nextNode(strt, attr, prnt)
		if node == nil {
			return nil, false
		}

		var lastNode *Node

		status := START
		for {
			tag, ctype, name, attr, idx := nextToken(Idx)
			Idx = idx

			if tag == BADTAG {
				fmt.Fprintf(os.Stderr, "\nERROR: Unparsable XML element\n")
				break
			}
			if tag == ISCLOSED {
				break
			}

			switch tag {
			case STARTTAG:
				if status == CHAR {
					if AllowEmbed {
						fmt.Fprintf(os.Stdout, "ERROR: UNRECOGNIZED MIXED CONTENT <%s> in <%s>\n", name, prnt)
					} else {
						fmt.Fprintf(os.Stdout, "ERROR: UNEXPECTED MIXED CONTENT <%s> in <%s>\n", name, prnt)
					}
				}
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
				status = STOP
			case STOPTAG:
				// pop out of recursive call
				return node, ok
			case CONTENTTAG:
				if ContentMods {
					node.Contents = CleanupContents(name, (ctype&ASCII) != 0, (ctype&AMPER) != 0, (ctype&MIXED) != 0)
				} else {
					node.Contents = name
				}
				status = CHAR
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
				status = OTHER
				// continue on same level
			default:
				status = OTHER
			}
		}

		return node, ok
	}

	// parseIndex recursive definition
	var parseIndex func(string, string, string) string

	// parse XML tags looking for trie index element
	parseIndex = func(strt, attr, prnt string) string {

		versn := ""

		// check for version attribute match
		if attr != "" && find.Versn != "" && strings.Contains(attr, find.Versn) {
			if strt == find.Match || find.Match == "" {
				if find.Parent == "" || prnt == find.Parent {
					attribs := ParseAttributes(attr)
					for i := 0; i < len(attribs)-1; i += 2 {
						if attribs[i] == find.Versn {
							versn = attribs[i+1]
						}
					}
				}
			}
		}

		// check for attribute index match
		if attr != "" && find.Attrib != "" && strings.Contains(attr, find.Attrib) {
			if strt == find.Match || find.Match == "" {
				if find.Parent == "" || prnt == find.Parent {
					attribs := ParseAttributes(attr)
					for i := 0; i < len(attribs)-1; i += 2 {
						if attribs[i] == find.Attrib {
							return attribs[i+1]
						}
					}
				}
			}
		}

		for {
			tag, _, name, attr, idx := nextToken(Idx)
			Idx = idx

			if tag == BADTAG {
				fmt.Fprintf(os.Stderr, "\nERROR: Unparsable XML element\n")
				break
			}
			if tag == ISCLOSED {
				break
			}

			switch tag {
			case STARTTAG:
				id := parseIndex(name, attr, strt)
				if id != "" {
					return id
				}
			case SELFTAG:
			case STOPTAG:
				// break recursion
				return ""
			case CONTENTTAG:
				// check for content index match
				if strt == find.Match || find.Match == "" {
					if find.Parent == "" || prnt == find.Parent {
						// append version if specified as parent/element@attribute^version
						if versn != "" {
							name += "."
							name += versn
						}
						return name
					}
				}
			default:
			}
		}

		return ""
	}

	// ParseXML

	// stream all tokens through callback
	if tokens != nil {

		for {
			tag, ctype, name, attr, idx := nextToken(Idx)
			Idx = idx

			if tag == BADTAG {
				fmt.Fprintf(os.Stderr, "\nERROR: Unparsable XML element\n")
				break
			}
			if tag == ISCLOSED {
				break
			}

			tkn := Token{tag, ctype, name, attr, idx, 0}

			tokens(tkn)

			if tag == ISCLOSED {
				break
			}
		}

		return nil, ""
	}

	// find value of index element
	if find != nil && find.Index != "" {

		// return indexed identifier

		tag, _, name, attr, idx := nextToken(Idx)

		// loop until start tag
		for {
			Idx = idx

			if tag == BADTAG {
				fmt.Fprintf(os.Stderr, "\nERROR: Unparsable XML element\n")
				break
			}
			if tag == ISCLOSED {
				break
			}

			if tag == STARTTAG {
				break
			}

			tag, _, name, attr, idx = nextToken(Idx)
		}

		return nil, parseIndex(name, attr, parent)
	}

	// otherwise create node tree for general data extraction
	tag, _, name, attr, idx := nextToken(Idx)

	// loop until start tag
	for {
		Idx = idx

		if tag == BADTAG {
			fmt.Fprintf(os.Stderr, "\nERROR: Unparsable XML element\n")
			break
		}
		if tag == ISCLOSED {
			break
		}

		if tag == STARTTAG {
			break
		}

		tag, _, name, attr, idx = nextToken(Idx)
	}

	top, ok := parseLevel(name, attr, parent)

	if !ok {
		return nil, ""
	}

	return top, ""
}

// StreamXML calls XML tokenizer parser on an XML reader, sends tokens for CDATA and COMMENT sections, and optionally tracks line numbers
func StreamXML(in *XMLReader, tokens func(Token)) {

	if in == nil || tokens == nil {
		return
	}

	// token parser variables
	Text := ""
	Txtlen := len(Text)
	Idx := 0

	// line tracking variables
	Line := 1
	Lag := 0

	// variables to track COMMENT or CDATA sections that span reader blocks
	Which := NOTAG
	SkipTo := ""

	// get next XML token
	nextToken := func(idx int) (TagType, ContentType, string, string, int, int) {

		if Text == "" {
			// buffer is empty, read next block if available
			if in != nil {
				Text = in.NextBlock()
				Txtlen = len(Text)
				Idx = 0
				idx = 0
				Lag = 0
			}

			if Text == "" {
				// signal end of XML data
				return ISCLOSED, NONE, "", "", 0, Line
			}
		}

		text := Text[:]
		txtlen := Txtlen

		updateLineCount := func(max int) {
			// count lines
			for i := Lag; i < max; i++ {
				if text[i] == '\n' {
					Line++
				}
			}
			Lag = idx
		}

		if Which != NOTAG && SkipTo != "" {
			which := Which
			// previous block ended inside CDATA object or COMMENT
			start := idx
			found := strings.Index(text[:], SkipTo)
			if found < 0 {
				// no stop signal found in next block
				str := text[:]
				if HasFlankingSpace(str) {
					str = strings.TrimSpace(str)
				}
				// signal end of current block
				Text = ""

				if CountLines {
					updateLineCount(txtlen)
				}

				// leave Which and SkipTo values unchanged as another continuation signal
				// send CDATA or COMMENT contents
				return which, NONE, str[:], "", 0, Line
			}
			// otherwise adjust position past end of skipTo string and return to normal processing
			idx += found
			str := text[start:idx]
			if HasFlankingSpace(str) {
				str = strings.TrimSpace(str)
			}
			idx += len(SkipTo)
			// clear tracking variables
			Which = NOTAG
			SkipTo = ""
			// send CDATA or COMMENT contents
			return which, NONE, str[:], "", idx, Line
		}

		// XML string, and all blocks, end with > character, acts as sentinel to check if past end of text
		if idx >= txtlen {
			// signal end of XML string or current block, will read next block on next call
			Text = ""

			if CountLines {
				updateLineCount(txtlen)
			}

			return NOTAG, NONE, "", "", 0, Line
		}

		// skip past leading blanks
		ch := text[idx]
		for InBlank[ch] {
			idx++
			ch = text[idx]
		}

		start := idx

		plainContent := true

		if AllowEmbed && ch == '<' {
			// check to see if an HTML or MathML element is at the beginning of a content string
			if HTMLAhead(text, idx) != 0 {
				plainContent = false
			}
		}

		if plainContent && ch == '<' {

			// at start of element
			idx++
			ch = text[idx]

			// check for legal first character of element
			if InFirst[ch] {

				// read element name
				start = idx
				idx++

				ch = text[idx]
				for InElement[ch] {
					idx++
					ch = text[idx]
				}

				str := text[start:idx]

				switch ch {
				case '>':
					// end of element
					idx++

					if CountLines {
						updateLineCount(idx)
					}

					return STARTTAG, NONE, str[:], "", idx, Line
				case '/':
					// self-closing element without attributes
					idx++
					ch = text[idx]
					if ch != '>' {
						fmt.Fprintf(os.Stderr, "\nSelf-closing element missing right angle bracket\n")
					}
					idx++

					if CountLines {
						updateLineCount(idx)
					}

					return SELFTAG, NONE, str[:], "", idx, Line
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

						if CountLines {
							updateLineCount(idx)
						}

						return SELFTAG, NONE, str[:], atr[:], idx, Line
					}
					atr := text[start:idx]
					idx++

					if CountLines {
						updateLineCount(idx)
					}

					return STARTTAG, NONE, str[:], atr[:], idx, Line
				default:
					if CountLines {
						fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element, line %d\n", ch, Line)
					} else {
						fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element\n", ch)
					}

					if CountLines {
						updateLineCount(idx)
					}

					return STARTTAG, NONE, str[:], "", idx, Line
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
					if InFirst[ch] {
						idx++
						ch = text[idx]
						for InElement[ch] {
							idx++
							ch = text[idx]
						}
						str := text[start:idx]
						if ch != '>' {
							fmt.Fprintf(os.Stderr, "\nUnexpected characters after end element name\n")
						}
						idx++

						if CountLines {
							updateLineCount(idx)
						}

						return STOPTAG, NONE, str[:], "", idx, Line
					}
					if CountLines {
						fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element, line %d\n", ch, Line)
					} else {
						fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element\n", ch)
					}
				case '?':
					// skip ?xml and ?processing instructions
					idx++
					ch = text[idx]
					for ch != '>' {
						idx++
						ch = text[idx]
					}
					idx++
					return NOTAG, NONE, "", "", idx, Line
				case '!':
					// skip !DOCTYPE, !COMMENT, and ![CDATA[
					idx++
					start = idx
					ch = text[idx]
					Which = NOTAG
					SkipTo = ""
					if ch == '[' && strings.HasPrefix(text[idx:], "[CDATA[") {
						Which = CDATATAG
						SkipTo = "]]>"
						start += 7
					} else if ch == '-' && strings.HasPrefix(text[idx:], "--") {
						Which = COMMENTTAG
						SkipTo = "-->"
						start += 2
					} else if ch == 'D' && strings.HasPrefix(text[idx:], "DOCTYPE") {
						Which = DOCTYPETAG
						SkipTo = ">"
					}
					if Which != NOTAG && SkipTo != "" {
						which := Which
						// CDATA or COMMENT block may contain internal angle brackets
						found := strings.Index(text[idx:], SkipTo)
						if found < 0 {
							// string stops in middle of CDATA or COMMENT
							str := text[start:]
							if HasFlankingSpace(str) {
								str = strings.TrimSpace(str)
							}
							// signal end of current block
							Text = ""

							if CountLines {
								updateLineCount(txtlen)
							}

							// leave Which and SkipTo values unchanged as another continuation signal
							// send CDATA or COMMENT contents
							return which, NONE, str[:], "", 0, Line
						}
						// adjust position past end of CDATA or COMMENT
						idx += found
						str := text[start:idx]
						if HasFlankingSpace(str) {
							str = strings.TrimSpace(str)
						}
						idx += len(SkipTo)
						// clear tracking variables
						Which = NOTAG
						SkipTo = ""
						// send CDATA or COMMENT contents
						return which, NONE, str[:], "", idx, Line
					}
					// otherwise just skip to next right angle bracket
					for ch != '>' {
						idx++
						ch = text[idx]
					}
					idx++
					return NOTAG, NONE, "", "", idx, Line
				default:
					if CountLines {
						fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element, line %d\n", ch, Line)
					} else {
						fmt.Fprintf(os.Stderr, "\nUnexpected punctuation '%c' in XML element\n", ch)
					}
				}
			}

		} else if ch != '>' {

			// at start of contents
			start = idx

			hasMarkup := false
			hasNonASCII := false

			// find end of contents
			if AllowEmbed {

				for {
					for InContent[ch] {
						idx++
						ch = text[idx]
					}
					// set flags to speed up conditional content processing
					if ch == '&' {
						idx++
						ch = text[idx]
						if ch == 'a' {
							if strings.HasPrefix(text[idx:], "amp;") {
								hasMarkup = true
							}
						} else if ch == 'g' {
							if strings.HasPrefix(text[idx:], "gt;") {
								hasMarkup = true
							}
						} else if ch == 'l' {
							if strings.HasPrefix(text[idx:], "lt;") {
								hasMarkup = true
							}
						}
						continue
					}
					if ch > 127 {
						hasNonASCII = true
						idx++
						ch = text[idx]
						continue
					}
					if ch == '<' {
						// optionally allow HTML text formatting elements and super/subscripts
						advance := HTMLAhead(text, idx)
						if advance > 0 {
							idx += advance
							ch = text[idx]
							plainContent = false
							continue
						}
					}
					break
				}
			} else {
				for ch != '<' && ch != '>' {
					idx++
					ch = text[idx]
				}
			}

			// trim back past trailing blanks
			lst := idx - 1
			ch = text[lst]
			for InBlank[ch] && lst > start {
				lst--
				ch = text[lst]
			}

			str := text[start : lst+1]

			ctype := PLAIN

			if AllowEmbed {
				if !plainContent {
					ctype |= MIXED
				}
				if hasMarkup {
					ctype |= AMPER
				}
				if hasNonASCII {
					ctype |= ASCII
				}
			}

			if CountLines {
				updateLineCount(idx)
			}

			return CONTENTTAG, ctype, str[:], "", idx, Line
		}

		return BADTAG, NONE, "", "", idx, Line
	}

	// StreamXML

	// stream all tokens through callback

	for {
		tag, ctype, name, attr, idx, line := nextToken(Idx)
		Idx = idx

		if tag == BADTAG {
			if CountLines {
				fmt.Fprintf(os.Stderr, "\nERROR: Unparsable XML element, line %d\n", line)
			} else {
				fmt.Fprintf(os.Stderr, "\nERROR: Unparsable XML element\n")
			}
			break
		}

		tkn := Token{tag, ctype, name, attr, idx, line}

		tokens(tkn)

		if tag == ISCLOSED {
			break
		}
	}
}

// specialized public ParseXML / StreamXML shortcuts

func ParseRecord(Text, parent string) *Node {

	pat, _ := ParseXML(Text, parent, nil, nil)

	return pat
}

func FindIdentifier(Text, parent string, find *Find) string {

	_, id := ParseXML(Text, parent, nil, find)

	return id
}

func StreamValues(Text, parent string, stream func(string, string)) {

	elementName := ""

	streamer := func(tkn Token) {

		switch tkn.Tag {
		case STARTTAG:
			elementName = tkn.Name
		case CONTENTTAG:
			// send element name and content to callback
			stream(elementName, tkn.Name)
		default:
		}
	}

	ParseXML(Text, parent, streamer, nil)
}

func CreateTokenizer(in *XMLReader) <-chan Token {

	if in == nil {
		return nil
	}

	out := make(chan Token, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create tokenizer channel\n")
		os.Exit(1)
	}

	// xmlTokenizer sends XML tokens through channel
	xmlTokenizer := func(rdr *XMLReader, out chan<- Token) {

		// close channel when all records have been processed
		defer close(out)

		// parse XML and send tokens through channel
		StreamXML(rdr, func(tkn Token) { out <- tkn })
	}

	// launch single tokenizer goroutine
	go xmlTokenizer(in, out)

	return out
}

// UNSHUFFLER USES HEAP TO RESTORE OUTPUT OF MULTIPLE CONSUMERS TO ORIGINAL RECORD ORDER

type Extract struct {
	Index int
	Ident string
	Text  string
	Data  []byte
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

func CreateProducer(pat, star string, rdr *XMLReader) <-chan Extract {

	if rdr == nil {
		return nil
	}

	out := make(chan Extract, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create producer channel\n")
		os.Exit(1)
	}

	// xmlProducer sends partitioned XML strings through channel
	xmlProducer := func(pat, star string, rdr *XMLReader, out chan<- Extract) {

		// close channel when all records have been processed
		defer close(out)

		// partition all input by pattern and send XML substring to available consumer through channel
		PartitionPattern(pat, star, rdr,
			func(rec int, ofs int64, str string) {
				out <- Extract{rec, "", str, nil}
			})
	}

	// launch single producer goroutine
	go xmlProducer(pat, star, rdr, out)

	return out
}

func CreateUnshuffler(inp <-chan Extract) <-chan Extract {

	if inp == nil {
		return nil
	}

	out := make(chan Extract, ChanDepth)
	if out == nil {
		fmt.Fprintf(os.Stderr, "\nERROR: Unable to create unshuffler channel\n")
		os.Exit(1)
	}

	// xmlUnshuffler restores original order with heap
	xmlUnshuffler := func(inp <-chan Extract, out chan<- Extract) {

		// close channel when all records have been processed
		defer close(out)

		// initialize empty heap
		hp := &ExtractHeap{}
		heap.Init(hp)

		// index of next desired result
		next := 1

		delay := 0

		for ext := range inp {

			// push result onto heap
			heap.Push(hp, ext)

			// read several values before checking to see if next record to print has been processed
			if delay < HeapSize {
				delay++
				continue
			}

			delay = 0

			for hp.Len() > 0 {

				// remove lowest item from heap, use interface type assertion
				curr := heap.Pop(hp).(Extract)

				if curr.Index > next {

					// record should be printed later, push back onto heap
					heap.Push(hp, curr)
					// and go back to waiting on input channel
					break
				}

				// send even if empty to get all record counts for reordering
				out <- Extract{curr.Index, curr.Ident, curr.Text, curr.Data}

				// prevent ambiguous -limit filter from clogging heap (deprecated)
				if curr.Index == next {
					// increment index for next expected match
					next++
				}

				// keep checking heap to see if next result is already available
			}
		}

		// send remainder of heap to output
		for hp.Len() > 0 {
			curr := heap.Pop(hp).(Extract)

			out <- Extract{curr.Index, curr.Ident, curr.Text, curr.Data}
		}
	}

	// launch single unshuffler goroutine
	go xmlUnshuffler(inp, out)

	return out
}

// CREATE COMMON DRIVER TABLES

// InitTables creates lookup tables to simplify the tokenizer
func InitTables() {

	for i := range InBlank {
		InBlank[i] = false
	}
	InBlank[' '] = true
	InBlank['\t'] = true
	InBlank['\n'] = true
	InBlank['\r'] = true
	InBlank['\f'] = true

	// first character of element cannot be a digit, dash, or period
	for i := range InFirst {
		InFirst[i] = false
	}
	for ch := 'A'; ch <= 'Z'; ch++ {
		InFirst[ch] = true
	}
	for ch := 'a'; ch <= 'z'; ch++ {
		InFirst[ch] = true
	}
	InFirst['_'] = true

	// remaining characters also includes colon for namespace
	for i := range InElement {
		InElement[i] = false
	}
	for ch := 'A'; ch <= 'Z'; ch++ {
		InElement[ch] = true
	}
	for ch := 'a'; ch <= 'z'; ch++ {
		InElement[ch] = true
	}
	for ch := '0'; ch <= '9'; ch++ {
		InElement[ch] = true
	}
	InElement['_'] = true
	InElement['-'] = true
	InElement['.'] = true
	InElement[':'] = true

	// embedded markup and math tags are lower case
	for i := range InLower {
		InLower[i] = false
	}
	for ch := 'a'; ch <= 'z'; ch++ {
		InLower[ch] = true
	}
	for ch := '0'; ch <= '9'; ch++ {
		InLower[ch] = true
	}
	InLower['_'] = true
	InLower['-'] = true
	InLower['.'] = true
	InLower[':'] = true

	// shortcut to find <, >, or &, or non-ASCII
	for i := range InContent {
		InContent[i] = false
	}
	for i := 0; i <= 127; i++ {
		InContent[i] = true
	}
	InContent['<'] = false
	InContent['>'] = false
	InContent['&'] = false
}

// GO AUTOMATIC INITIALIZER

func init() {
	InitTables()
}
