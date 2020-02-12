# Fasta-Region-Inspector: A Somatic Hypermutation Analysis Tool

## Introduction 

**Fasta-Region-Inspector (FRI)** is a computational tool for analyzing somatic hypermutation.<br/>
This haskell script takes in variant information, corresponding region information, and ambiguity code string(s) to determine:
1. Whether the user-defined variants are within the transcription start site (TSS) of the corresponding gene.
2. All possible start locations of mapped ambiguity code strings within 2 Kb of TSS.
3. Final list of user-defined variants that lie within an mapped ambiguity code string inside of a 2 Kb window of TSS of corresponding gene.<br/><br/>

**FRI** outputs three files:<br/>
1. **variants.tsv** - Outputs list of user-defined variants and corresponding region information and whether the variant is within 2 Kb (Y/N).
2. **ambiguity_codes.tsv** - Outputs list of mapped ambiguity codes strings, along with corresponding region information, and all start locations where these mapped ambiguity code strings are found.
3. **variants_in_ambiguity_codes.tsv** - Outputs final list of user-defined variants that are within 2 Kb of TSS of corresponding gene that lie within a mapped ambiguity code string.

Downstream analysis is left to the user, as the output files can easily be filtered using a scripting language like **awk**.  Filtering on reference/alternate bases, strand orientation (1 vs. -1), and SYMBOL (gene) are some examples of ways to continue narrowing down the output of **FRI**.

## Theory and Implementation

String-searching plays a large role in this program, as it makes up an overwhelming proportion of the programs runtime.  An implementation of the `Knuth-Morris-Pratt` algorithm provided by the [stringsearch](https://hackage.haskell.org/package/stringsearch-0.3.6.6/docs/Data-ByteString-Search-DFA.html) is used to find all possible locations of mapped ambiguity code strings across the 2 Kb window of the TSS.<br/>
The `Knuth-Morris-Pratt` implementation was chosen as opposed to other string-searching algorithms like `Boyer-Moore` or `Rabin-Karp` due to its strength with low complexity search alphabets (ATGC). <br/><br/>

Generating all posssible strings that match a regular expression is needed for this program, as the user defined ambiguity code string(s) can be thought of as regular expressions.<br/>
[Thompson's construction algorithm](https://en.wikipedia.org/wiki/Thompson%27s_construction) can be used to create a nondeterministic finite automaton (NFA) to solve this problem.<br/>
An example of tranforming an ambiguity code string into a regular expression would be the following:<br/>
**ambiguity code**: `WRCY` <-> **regular expression**: `[A|T][A|G]C[C|T]`<br/>
**FRI** utilizes the [sbv](https://hackage.haskell.org/package/sbv) package to generate all possible strings from user-defined ambiguity code(s), which utilizes a theorem prover to construct a NFA.<br/>

## Prerequisites

**fri.hs** assumes you have a the [GHC](https://www.haskell.org/ghc/) compiler and packages installed that it imports.  The easiest way to do this is to download the [Haskell Platform](https://www.haskell.org/platform/).<br/><br/>

## Installing required packages

To install the peripheral packages **fri.hs** requires, you can call the following command assuming you have [cabal](https://www.haskell.org/cabal/), a package manager and build system for Haskell, installed on your system (it comes with the [Haskell Platform](https://www.haskell.org/platform/)).<br/><br/>
`$ cabal install [packagename]`<br/><br/>

**Required packages**
- Bio.Core.Sequence
- Bio.Sequence.Fasta 
- Control.DeepSeq 
- Data.ByteString 
- Data.ByteString.Char8 
- Data.ByteString.Lazy 
- Data.ByteString.Search.DFA 
- Data.Char 
- Data.Functor 
- Data.List 
- Data.List.Split 
- Data.Ord 
- Data.SBV 
- Data.SBV.String 
- Data.SBV.RegExp
- Data.Traversable 
- System.Console.GetOpt 
- System.Process 
- System.Environment 
- System.Exit 
- System.IO 
- System.IO.Temp 
- Text.PrettyPrint.Boxes 
- Text.Regex.TDFA

## Input

**FRI** requires three inputs:<br/><br/>

1. **Variant Input** - This input tsv file needs to have the following fields:<br/><br/>

`Sample\tSymbol\tChromosome\tStart\tStop\tRef\tAlt`<br/><br/>

This file hold all of variant information.<br/><br/>

2. **Region Input** - This input tsv file needs to have the following fields:<br/><br/>

`Chromosome\tTSS\tStrand\tGene_Name`<br/><br/>

This file holds all of the corresponding region information.<br/><br/>

To create this file, you will typically take each variant of interest's corresponding `ENST`, and query [bioMart](https://useast.ensembl.org/biomart/martview) to return the following fields:<br/>

`Gene stable ID`<br/>
`Transcript stable ID`<br/>
`Strand`<br/>
`Chromosome/scaffold name`<br/>
`Transcription start site (TSS)`<br/> 
`Gene name`<br/><br/>

You can then subset this file to contain only fields `Chromosome/scaffold name`, `Transcription start site (TSS)`, `Strand` and `Gene name`.<br/><br/>

**IMPORTANT: Remove the "chr" before the chromosome number.**<br/><br/>

3. **Ambiguity Codes String** - This string argument describes the ambiguity codes to search for within the TSS of each gene.<br/><br/>

The user may define as many ambiguity code strings as desired, but keep in mind, the more ambiguity codes you input, the more strings that need to be searched for within the TSS.<br/><br/>

The string input is **semicolon** delimited and the string needs to **start** and **end** with a **semicolon**.<br/><br/>

To following is an example string:<br/><br/>

Assume the user wants to search for ``, ``, `` and `` in the TSS.<br/><br/>

**Ambiguity Code String** -> `;WRCY;WRC;YYGG;CCGY;`<br/><br/>
