{-=FastaRegionInspector (FRI): A Haskell-based solution to=-}
{-=analyze specific regions of a fasta file for specific=-}
{-=ambiguity codes and variants.=-}
{-=Author: Matthew Mosior=-}
{-=Version: 1.0=-}
{-=Synopsis:  This Haskell Script will take in=-}
{-=will access user specified regions of a user specified fasta=-}
{-=file and search this region for user specified ambiguity codes=-}
{-=(i.e. WRCY) and will return wether user specified variants=-}
{-=are found in these specified ambiguity codes.=-}


{-Syntax extension by language pragma.-}

{-# LANGUAGE OverloadedStrings #-}

{--------------------------------------}

{-Imports-}

import Bio.Core.Sequence as BCS
import Bio.Sequence.Fasta as BSF
import Control.DeepSeq as CD
import Control.Parallel.Strategies as CPS
import Data.ByteString as DB
import Data.ByteString.Char8 as DBC
import Data.ByteString.Lazy as DBL
import Data.ByteString.Search.DFA as DBSDFA
import Data.Char as DC
import Data.Functor as DF
import Data.List as DL
import Data.List.Split as DLS
import Data.Ord as DO
import Data.SBV as DSBV
import qualified Data.SBV.String as DSBVS
import Data.SBV.RegExp as DSBVRE
import Data.Traversable as DT
import System.Console.GetOpt as SCG
import System.Process as SP
import System.Environment as SE
import System.Exit as SX
import System.IO as SIO
import System.IO.Temp as SIOT
import Text.PrettyPrint.Boxes as TPB
import Text.Regex.TDFA as TRTDFA

{---------}


{-Custom CML Option Datatype.-}

data Flag
    = Verbose                -- -v
    | Version                -- -V -?
    | OutputDirectory String -- -o
    | TSSWindowSize   String -- 
    | Help                   -- --help
    deriving (Eq,Ord,Show)

{-----------------------------}


{-Custom bool functions for Flag Datatype.-}

--isOutputDirectory -> This function will
--test for OutputFile flag.
isOutputDirectory :: Flag -> Bool
isOutputDirectory (OutputDirectory _) = True
isOutputDirectory _                   = False

--isTSSWindowSize -> This function will
--test for TSSWindowSize flag.
isTSSWindowSize :: Flag -> Bool
isTSSWindowSize (TSSWindowSize _) = True
isTSSWindowSize _                 = False

{------------------------------------------}


{-Custom extraction functions for Flag Datatype.-}

--extractOutputDirectory -> This function will
--extract the string associated with 
--OutputFile.
extractOutputDirectory :: Flag -> String
extractOutputDirectory (OutputDirectory x) = x

--extractTSSWindowSize -> This function will
--extract the string associated with 
--TSSWindowSize.
extractTSSWindowSize :: Flag -> String
extractTSSWindowSize (TSSWindowSize x) = x 

{------------------------------------------------}


{-Custom extraction functions for Sequence Datatype.-}

--extractSeqData -> This function will
--extract the SeqData field of Seq constructor 
--of the Sequence Datatype.
extractSeqData :: Sequence -> SeqData
extractSeqData (Seq _ x _) = x

--extractSeqLabel -> This function will
--extract the SeqLabel field of Seq constructor
--of the Sequence Datatype.
extractSeqLabel :: Sequence -> SeqLabel
extractSeqLabel (Seq x _ _) = x

--extractQualData -> This function will
--extract the QualData field of Seq constructor
--of the Sequence Datatype.
extractQualData :: Sequence -> Maybe QualData
extractQualData (Seq _ _ x) = x

--extractunSD -> This function will
--extract the unSD field of SeqData constructor.
extractunSD :: SeqData -> DBL.ByteString
extractunSD (SeqData unSD) = unSD

--extractunSL -> This function will
--extract the unSL field of SeqLabel constructor.
extractunSL :: SeqLabel -> DBL.ByteString
extractunSL (SeqLabel unSL) = unSL

--extractunQD -> This function will
--extract the unQD field of QualData constructor.
extractunQD :: QualData -> DBL.ByteString
extractunQD (QualData unQD) = unQD

{----------------------------------------------------}

{-Option Description function relating to datatype above.-}

--options -> This function will
--describe flags.
options :: [OptDescr Flag]
options =
    [ Option ['v']     ["verbose"]             (NoArg Verbose)                         "Output on stderr.",
      Option ['V','?'] ["version"]             (NoArg Version)                         "Show version number.",
      Option ['o']     ["outputdirectory"]     (ReqArg OutputDirectory "OUTDIRECTORY") "The directory path where output files will be printed.",
      Option []        ["TSSwindowsize"]       (ReqArg TSSWindowSize "TSSWINSIZE")     "The size of the window of which to search each region from the TSS.",
      Option []        ["help"]                (NoArg Help)                            "Print this help message."
    ]

--compilerOpts -> This function will
--parse incoming command line arguments.
compilerOpts :: [String] -> IO ([Flag],[String])
compilerOpts argv =
    case getOpt Permute Main.options argv of
        (args,files,[]) ->
            if DL.elem Help args
                then do SIO.hPutStrLn stderr (greeting ++ SCG.usageInfo header Main.options)
                        SX.exitWith SX.ExitSuccess
                else if DL.elem Version args
                    then do SIO.hPutStrLn stderr (greeting ++ version ++ SCG.usageInfo header Main.options)
                            SX.exitWith SX.ExitSuccess
                    else if (DL.length (DL.filter (isOutputDirectory) args) < 1)
                        then do SIO.hPutStrLn stderr (odmisserror ++ github ++ SCG.usageInfo header Main.options)
                                SX.exitWith (SX.ExitFailure 1)
                        else if (DL.length files > 4 || DL.length files < 4)
                            then do SIO.hPutStrLn stderr (flerror ++  acmisserror ++ github ++ SCG.usageInfo header Main.options)
                                    SX.exitWith (SX.ExitFailure 1)
                            else if (DL.length (files DL.!! 3) > 0) &&
                                    (not (ambiguityCodesCheck (files DL.!! 3)))
                                then do SIO.hPutStrLn stderr (acerror ++ github ++ SCG.usageInfo header Main.options)
                                        SX.exitWith (SX.ExitFailure 1)
                                else if (DL.length (DL.filter (isTSSWindowSize) args) > 0) &&
                                        (not (tssWindowSizeCheck (extractTSSWindowSize (DL.head (DL.filter (isTSSWindowSize) args))))) 
                                    then do SIO.hPutStrLn stderr (tsswinsizeerror ++ github ++ SCG.usageInfo header Main.options)
                                            SX.exitWith (SX.ExitFailure 1)
                                    else return (DL.nub args, files)
        (_,_,errors) -> do
            SIO.hPutStrLn stderr (DL.concat errors ++ SCG.usageInfo header Main.options)
            SX.exitWith (SX.ExitFailure 1)
        where
            greeting        = "Fasta Region Inspector, Copyright (c) 2020 Matthew Mosior.\n"
            header          = "Usage: fri [-vV?o] [Fasta File] [Region File] [Variant File] [Ambiguity Codes String]"
            version         = "Fasta Region Inspector (FRI), Version 1.0.\n"
            github          = "Please see https://github.com/Matthew-Mosior/Fasta-Region-Inspector/wiki for more information.\n"
            flerror         = "Incorrect number of input files:  Please provide exactly three files.\n\
                              \First input file  -> Fasta file\n\
                              \Second input file -> Region file\n\
                              \Third input file -> Variant File\n"
            odmisserror     = "Output directory argument missing.\n\
                              \Please provide the output directory."
            acmisserror     = "Ambiguity codes string missing.\n\
                              \Please provide ambiguity codes string with the following structure:\n\
                              \;[CODE1];[CODE2];...;[CODEN];\n\
                              \Please see https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html for more information on ambiguity codes.\n"
            acerror         = "Incorrect structure of ambiguity codes string.\n\
                              \Please provide ambiguity codes string with the following structure:\n\
                              \;[CODE1];[CODE2];...;[CODEN];\n\
                              \Please see https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html for more information on ambiguity codes.\n"
            tsswinsizeerror = "Incorrect structure of TSS window size string.\n\
                              \This string should only contain digits (0..9).\n"
          
{---------------------------------------------------------}


{-General Utility Functions.-}

--lineFeed -> This function will
--read the file in and split on
--whitespace, returning a list
--of lists.
lineFeed :: String -> [[String]]
lineFeed [] = []
lineFeed xs = DL.map DL.words (DL.lines xs)

--mapNotLast -> This function will
--work like the traditional map 
--function in Data.List, but not
--map to the last element of a list.
mapNotLast :: (a -> a) -> [a] -> [a]
mapNotLast fn []     = []
mapNotLast fn [x]    = [x]
mapNotLast fn (x:xs) = fn x : mapNotLast fn xs

--stringToTuple -> This function will
--annotate all mapped strings with
--directionality (tuple). 
stringToTuple :: [[String]] -> [[(String,String)]]
stringToTuple [] = []
stringToTuple xs = DL.map (DL.map (\y -> (y,"1"))) (DL.take ((DL.length xs) `div` 2) xs) 
                ++ DL.map (DL.map (\y -> (y,"-1"))) (DL.drop ((DL.length xs) `div` 2) xs)

--bslToStr -> This function will
--Convert from Bytestring (Lazy) to String.
bslToStr :: DBL.ByteString -> String
bslToStr = DL.map (DC.chr . fromEnum) . DBL.unpack

--strToBSC8 -> This function will
--convert Strings to Bytestring (Char8).
strToBSC8 :: String -> DBC.ByteString
strToBSC8 xs = DBC.pack xs

{----------------------------}

{-Ambiguity Codes functions.-}

--ambiguityCodesCheck -> This function will
--check that the ambiguity codes string 
--provided by the user has the correct structure.
ambiguityCodesCheck :: String -> Bool
ambiguityCodesCheck xs = if DL.head xs == ';' &&
                            DL.last xs == ';' &&
                            DL.all (\y -> y `DL.elem` nucleotideAmbiguityCodes) (DL.filter (\x -> x /= ';') xs) 
                             then True
                             else False 
    where
        --Local definition.--
        --nucleotideAmbiguityCodes
        nucleotideAmbiguityCodes = ['A','G','C','T','Y','R','W','S','K','M','D','V','H','B','X','N','-']
        ---------------------

{----------------------------}


{-TSSWindowSize Functions.-}

--tssWindowSizeCheck -> This function will
--check that the TSSWindowSize string
--provided by the user has the correct structure.
tssWindowSizeCheck :: String -> Bool
tssWindowSizeCheck xs = if (DL.all (DC.isDigit) xs)
                            then True
                            else False

{--------------------------}


{-Region file functions.-}

--variantWithinRegionCheck -> This function will
--check to see if the given variant is within genes
--2 kb (or custom TSS window size if provided)
--of the genes TSS.
variantWithinRegionCheck :: [Flag] -> [[String]] -> [[String]] -> [[([String],[String],Char)]]
variantWithinRegionCheck []    []     []    = []
variantWithinRegionCheck []    []     (_:_) = []
variantWithinRegionCheck (_:_) []     _     = []
variantWithinRegionCheck opts  (x:xs) ys    = if DL.length (DL.filter (isTSSWindowSize) opts) > 0
                                                  then do --Grab just "TSSWINSIZE".
                                                      let tsswinsize = DL.head (DL.filter (isTSSWindowSize) opts)    
                                                      --Grab the string from tsswinsize.
                                                      let tsswinsizestring = extractTSSWindowSize tsswinsize
                                                      --Check to see if variant is within tsswinsizestring range of region.
                                                      [variantRegionCheck x ys (read tsswinsizestring)] ++ (variantWithinRegionCheck opts xs ys)
                                                  else do --User did not provide TSSWindowSize argument.
                                                          --Use 2 kb as default TSS window size. 
                                                          --Call variantsWithinRegionCheckDefault and recurse.
                                                          [variantRegionCheck x ys 2000] ++ (variantWithinRegionCheck opts xs ys)
    where
        --Local function definition.--
        --variantRegionCheck -> This function will
        --check to see if current variant is within
        --region of gene.
        variantRegionCheck :: [String] -> [[String]] -> Int -> [([String],[String],Char)]
        variantRegionCheck [] [] _ = []
        variantRegionCheck x  ys z = windowChecker x (DL.filter (\a -> (a DL.!! 3) == (x DL.!! 1)) ys) z

        --windowChecker -> This function will
        --check variants position versus region
        --positions.
        windowChecker :: [String] -> [[String]] -> Int -> [([String],[String],Char)]
        windowChecker []     []     _  = []
        windowChecker (_:_)  []     _  = []
        windowChecker x      (y:ys) z  = if ((y !! 2) == "-1") &&
                                            (((read (y !! 1) :: Int) - z) <= (read (x !! 3) :: Int)) &&
                                            ((read (x !! 3) :: Int) <= (read (y !! 1) :: Int))
                                              then [(x,y,'Y')] ++ (windowChecker x ys z)
                                              else if ((y !! 2) == "1") &&
                                                      ((read (y !! 1) :: Int) <= (read (x !! 3) :: Int)) &&
                                                      ((read (x !! 3) :: Int) <= (read (y !! 1) :: Int) + z)
                                                        then [(x,y,'Y')] ++ (windowChecker x ys z)
                                                        else [(x,y,'N')] ++ (windowChecker x ys z)
        ------------------------------

--prepareWithinTSS -> This function will
--prepare withintss (see main) for
--printFile function.
prepareWithinTSS :: [[([String],[String],Char)]] -> [[String]]
prepareWithinTSS []     = []
prepareWithinTSS (x:xs) = [convertToList x] ++ (prepareWithinTSS xs)
    where
        --Local definitions.--
        --convertToList -> This function will
        --convert 4-tuple to list.
        convertToList :: [([String],[String],Char)] -> [String]
        convertToList [] = []
        convertToList xs = DL.concat (DL.concat (DL.map (\(a,b,c) -> [[DL.intercalate ":" a] ++ [DL.intercalate ":" b] ++ [[c]]]) xs))
        ----------------------

{------------------------}


{-Ambiguity Code Functions.-}

--ambiguityCodesWithinRegionCheck -> This function will
--check to see if the ambiguity codes are within genes
--2 kb (or custom TSS window size if provided)
--of the genes TSS.
ambiguityCodesWithinRegionCheck :: [(String,String)] -> [[(String,String)]] -> [[String]] -> [Sequence] -> [Flag] -> [[(String,[String],[String],[[Int]])]]
ambiguityCodesWithinRegionCheck []     []      []    []    []    = []
ambiguityCodesWithinRegionCheck []     []      []    []    (_:_) = []
ambiguityCodesWithinRegionCheck []     []      []    (_:_) _     = []
ambiguityCodesWithinRegionCheck []     (_:_)   []    _     _     = []
ambiguityCodesWithinRegionCheck (_:_)  _       []    _     _     = []
ambiguityCodesWithinRegionCheck []     []      (_:_) _     _     = []
ambiguityCodesWithinRegionCheck (_:_)  []      (_:_) _     _     = []
ambiguityCodesWithinRegionCheck []     (_:_)   (_:_) _     _     = []
ambiguityCodesWithinRegionCheck (x:xs) (r:rs)  ys    zs    opts  = [ambiguityCodesWithinRegionCheckSmall x r ys zs opts] ++ (ambiguityCodesWithinRegionCheck xs rs ys zs opts)
--ambiguityCodesWithinRegionCheckSmall -> This function will
--check to see if the ambiguity codes are within genes
--2 kb (or custom TSS window size if provided)
--of the genes TSS.
ambiguityCodesWithinRegionCheckSmall :: (String,String) -> [(String,String)] -> [[String]] -> [Sequence] -> [Flag] -> [(String,[String],[String],[[Int]])]
ambiguityCodesWithinRegionCheckSmall ([],[]) []    []       []    []    = []
ambiguityCodesWithinRegionCheckSmall ([],[]) []    []       []    (_:_) = []
ambiguityCodesWithinRegionCheckSmall ([],[]) []    []       (_:_) _     = []
ambiguityCodesWithinRegionCheckSmall ([],[]) (_:_) []       _     _     = []
ambiguityCodesWithinRegionCheckSmall (_,_)   _     []       _     _     = []
ambiguityCodesWithinRegionCheckSmall xs      rs    (y:ys)   zs    opts  = if y DL.!! 2 == snd xs
                                                                              then if DL.length (DL.filter (isTSSWindowSize) opts) > 0
                                                                                  then do --Grab just "TSSWINSIZE".
                                                                                          let tsswinsize = DL.head (DL.filter (isTSSWindowSize) opts)
                                                                                          --Grab the string from tsswinsize.
                                                                                          let tsswinsizestring = extractTSSWindowSize tsswinsize
                                                                                          --Grab locations of mapped am codes,
                                                                                          --and recurse.
                                                                                          [(fst xs,DL.map (fst) rs,y,subStrLocations (DL.map (fst) rs) y (read tsswinsizestring) zs)] ++ (ambiguityCodesWithinRegionCheckSmall xs rs ys zs opts)
                                                                                  else do --User did not provide TSSWindowSize argument.
                                                                                          --Use 2 kb as default TSS window size.
                                                                                          [(fst xs,DL.map (fst) rs,y,subStrLocations (DL.map (fst) rs) y 2000 zs)] ++ (ambiguityCodesWithinRegionCheckSmall xs rs ys zs opts)
                                                                          else do --Current ambiguity codes and mapped strings
                                                                                  --are not correct for region strand 
                                                                                  --(i.e. "-1" != "1" or "1" != "-1").  
                                                                                  ambiguityCodesWithinRegionCheckSmall xs rs ys zs opts       
    where
        --Local definitions.--
        --subStrLocations -> This function will
        --find the locations for all given substrings
        --found using allStrGeneration.
        subStrLocations :: [String] -> [String] -> Int -> [Sequence] -> [[Int]]
        subStrLocations []     []    _  []    = []
        subStrLocations []     []    _  (_:_) = []
        subStrLocations []     (_:_) _  _     = []     
        subStrLocations xs     ys    zs rs    = if (ys DL.!! 2 == "-1")
                                                    then DL.map (DL.map (\i -> ((((read (ys DL.!! 1)) - (zs)) + i) + 2))) (subStrLocationsSmallReverse xs ys zs rs)
                                                    else DL.map (DL.map (\i -> ((read (ys DL.!! 1)) + i))) (subStrLocationsSmallForward xs ys zs rs)
        --subStrLocationsSmallReverse -> This function will
        --find the locations for all given substrings
        --found using allStrGeneration.
        subStrLocationsSmallReverse :: [String] -> [String] -> Int -> [Sequence] -> [[Int]]
        subStrLocationsSmallReverse []     []    _  []    = []
        subStrLocationsSmallReverse []     []    _  (_:_) = []
        subStrLocationsSmallReverse []     (_:_) _  _     = []
        subStrLocationsSmallReverse (x:xs) ys    zs rs    = [DL.map (\a -> (DBC.length ((grabRegionSequence (grabFastaSequence (read (y !! 0) :: Int) rs) ys zs))) - a - 1) 
                                                                     (DBSDFA.indices (DBC.pack x) 
                                                                     (reverseComplementNucleotide 
                                                                     (grabRegionSequence 
                                                                     (grabFastaSequence (read (y !! 0) :: Int) rs) ys zs)))] 
                                                         ++ (subStrLocationsSmallReverse xs ys zs rs)
        --subStrLocationsSmallForward -> This function will
        --find the locations for all given substrings
        --found using allStrGeneration.
        subStrLocationsSmallForward :: [String] -> [String] -> Int -> [Sequence] -> [[Int]]
        subStrLocationsSmallForward []     []    _  []    = []
        subStrLocationsSmallForward []     []    _  (_:_) = []
        subStrLocationsSmallForward []     (_:_) _  _     = []
        subStrLocationsSmallForward (x:xs) ys    zs rs    = [DBSDFA.indices (DBC.pack x) (grabRegionSequence (grabFastaSequence (read (y !! 0) :: Int) rs) ys zs)] ++ (subStrLocationsSmallForward xs ys zs rs)
        --grabRegionSequence -> This function will
        --grab the region of the correct chr 
        --returned from grabFastaSequence.
        grabRegionSequence :: DBC.ByteString -> [String] -> Int -> DBC.ByteString 
        grabRegionSequence xs ys zs = if (ys DL.!! 2 == "-1")
                                          then DBC.take (zs) (DBC.drop ((read (ys DL.!! 1)) - (zs - 1)) (xs))
                                          else DBC.take (zs) (DBC.drop ((read (ys DL.!! 1)) - 1) (xs))
        --grabFastaSequence -> This function will
        --grab the correct fasta sequence
        --using chromosome information
        --in the region file.
        grabFastaSequence :: Int -> [Sequence] -> DBC.ByteString
        grabFastaSequence x ys = smallGrabFastaSequence x ys [0..(DL.length ys) - 1]
        --smallGrabFastaSequence -> This function will
        --grab the correct fasta sequence
        --using chromosome information
        --in the region file.
        smallGrabFastaSequence :: Int -> [Sequence] -> [Int] -> DBC.ByteString
        smallGrabFastaSequence _ _ [] = DBC.empty
        smallGrabFastaSequence x ys (z:zs) = if ((bslToStr (extractunSL (extractSeqLabel (ys !! z)))) == ("chr" ++ show x))
                                                 then strToBSC8 (bslToStr (extractunSD (extractSeqData (ys !! z))))
                                                 else smallGrabFastaSequence x ys zs
        --reverseComplementNucleotide -> This function will
        --return the reverse complement of a DBC.ByteString.
        reverseComplementNucleotide :: DBC.ByteString -> DBC.ByteString
        reverseComplementNucleotide xs = DBC.pack (DL.map (snd) (DL.concatMap (\y -> DL.filter (\(r,_) -> r == y) revcomplementmapping) (DBC.unpack (DBC.reverse xs))))
            where
                --Local definitions.--
                --revcomplementmapping -> List containing reverse
                --complement mapping for nucleotides.
                --(NUCLEOTIDE,NUCLEOTIDE_REV_COMPLEMENT).
                revcomplementmapping = [('A','T'),('T','A'),('G','C'),('C','G')]
                ----------------------

--allStrGeneration -> This function will
--return all substring locations with all
--ambiguity code mapped to true nucleotides.
allStrGeneration :: [String] -> IO [[String]]
allStrGeneration [] = return []
allStrGeneration xs = mapM (allStrGenerationSmall) xs
--allStrGenerationSmall -> This function will
--return all substring locations with current
--ambiguity code mapped to true nucleotides.
allStrGenerationSmall :: String -> IO [String]
allStrGenerationSmall [] = return []
allStrGenerationSmall xs = do --Need to create all possible strings created from nucleotidemapping.
                             --Use Data.SBV library to create all possible satisfiability predicates
                             --(Strings generated from a particular regular expression pattern).
                             generatedregexstrs <- DSBV.allSat $ \s -> (s :: SString) `DSBVRE.match` (DSBVRE.Conc (DL.map (DSBVRE.oneOf) nucleotidemappingfinal))
                             --Return only true results of generatedregexstrs.        
                             return (DL.filter (\z -> not (DL.null z)) (DL.filter (\y -> DL.all (DC.isUpper) y) (DLS.splitOneOf " " (DL.filter (\x -> not (x `DL.elem` ("\"" :: String))) (show generatedregexstrs))))) 
    where
        --charConversion -> This function will
        --compute the conversion to prepare for
        --the satisfiability predicates.
        charConversion :: String -> [(Char,Char)] -> [(Char,String)]
        charConversion []     [] = []
        charConversion _      [] = []
        charConversion []     _  = []
        charConversion (x:xs) ys = [(DL.head (DL.map (fst) (DL.filter (\(a,_) -> a == x) ys)),DL.map (snd) (DL.filter (\(a,_) -> a == x) ys))] ++ (charConversion xs ys)
        --Local function and variable definitions.--
        --nucleotidemappingfinal -> This returns only the mapped nucleotides
        --of nucleotidemapping that fit the user defined ambiguity code string.
        nucleotidemappingfinal    = DL.map (snd) (charConversion xs nucleotidemappingfiltered) 
        --nucleotidemappingfiltered -> This returns only the mapped tuples
        --of nucleotidemapping that fit the user defined ambiguity code string.
        nucleotidemappingfiltered = DL.filter (\(b,_) -> b `DL.elem` (DL.filter (\a -> a `DL.elem` DL.map (fst) nucleotidemapping) xs)) nucleotidemapping 
        --nucleotidemappings -> List containing ambiguity code mapping (AMBIGUITY_CODE,MAPPEDNUCLEOTIDE). 
        nucleotidemapping = [('A','A'),('T','T'),('G','G'),('C','C'),('Y','C'),('Y','T'),('R','A'),('R','G'),('W','A'),('W','T'),('S','G'),('S','C'),('K','T'),('K','G'),('M','C'),('M','A'),('D','A'),('D','G'),('D','T'),('V','A'),('V','C'),('V','G'),('H','A'),('H','C'),('H','T'),('B','C'),('B','G'),('B','T'),('N','A'),('N','T'),('N','G'),('N','C'),('X','A'),('X','T'),('X','G'),('X','C')]
        --------------------------------------------

--ambiguityCodesReverseComplement -> This function will
--calculate the reverse complement for user specified 
--ambiguity codes.
ambiguityCodesReverseComplement :: [String] -> [String]
ambiguityCodesReverseComplement []     = []
ambiguityCodesReverseComplement (x:xs) = [reversecomplementfinal] ++ (ambiguityCodesReverseComplement xs) 
    where
        --Local definitions.--
        --reversecomplementfinal -> This returns only the
        --mapped ambiguity codes of reversecomplementfiltered 
        --reversed and in the correct order.
        reversecomplementfinal = DL.map (snd) (DL.concatMap (\y -> DL.filter (\(r,_) -> r == y) reversecomplementfiltered) (DL.reverse x)) 
        --reversecomplementfiltered -> This returns only the mapped tuples
        --of reversecomplementmapping that fit the user defined ambiguity code string.
        reversecomplementfiltered = DL.filter (\(b,_) -> b `DL.elem` (DL.filter (\a -> a `DL.elem` DL.map (fst) revcomplementmapping) x)) revcomplementmapping
        --revcomplementmapping -> List containing reverse 
        --complement mapping for ambiguity codes 
        --(AMBIGUITY_CODE,AMBIGUITY_CODE_REV_COMPLEMENT).
        revcomplementmapping = [('A','T'),('T','A'),('G','C'),('C','G'),('Y','R'),('R','Y'),('W','W'),('S','S'),('K','M'),('M','K'),('D','H'),('H','D'),('V','B'),('B','V'),('X','X'),('N','N'),('-','-')] 
        ---------------------- 

--prepareAmbiguityCodesWithinTSS -> This function will
--prepare ambiguitycodeswithintss (see main) for
--printFile function.
prepareAmbiguityCodesWithinTSS :: [[(String,[String],[String],[[Int]])]] -> [[[String]]]
prepareAmbiguityCodesWithinTSS []     = []
prepareAmbiguityCodesWithinTSS (x:xs) = (convertToList x) ++ (prepareAmbiguityCodesWithinTSS xs)
    where
        --Local definitions.--
        --convertToList -> This function will
        --convert 4-tuple to list.
        convertToList :: [(String,[String],[String],[[Int]])] -> [[[String]]]
        convertToList [] = [] 
        convertToList xs = DL.map (\y -> tupleToList y) (DL.map (\x -> tupleConverter x) xs)
        --tupleToList -> This function will
        --convert a 4-tuple to a [
        tupleToList :: ([[String]],[[String]],[[String]],[[Int]]) -> [[String]]
        tupleToList ([],    _,     _,     _)      = []
        tupleToList ((_:_), [],    _,     _)      = []
        tupleToList ((_:_), (_:_), [],    _)      = []
        tupleToList ((_:_), (_:_), (_:_), [])     = []
        tupleToList ((a:as),(b:bs),(c:cs),(d:ds)) = [a ++ b ++ c ++ (DL.map (show) d)] ++ (tupleToList (as,bs,cs,ds)) 
        --tupleConverter -> This function will
        --convert a 4-tuple to the corrected
        --4-tuple.
        tupleConverter :: (String,[String],[String],[[Int]]) -> ([[String]],[[String]],[[String]],[[Int]])
        tupleConverter (a,b,c,d) = ((DL.map (\x -> [x]) (DL.replicate (DL.length b) a)),(DL.map (\x -> [x]) b),(DL.replicate (DL.length b) c),d)  
        ----------------------

{---------------------------}


{-Variant functions.-}

--variantsWithinAmbiguityCodesAndTSS ->  This function will
--identify all variants that are within ambiguity codes
--and corresponding genes TSS.
variantsWithinAmbiguityCodesAndTSS :: [[([String],[String],Char)]] -> [[[String]]] -> [[String]]
variantsWithinAmbiguityCodesAndTSS []     [] = []
variantsWithinAmbiguityCodesAndTSS _      [] = []
variantsWithinAmbiguityCodesAndTSS []     _  = []
variantsWithinAmbiguityCodesAndTSS (x:xs) ys = (variantsWithinAmbiguityCodesAndTSSSmall variantsfiltered ys) ++ (variantsWithinAmbiguityCodesAndTSS xs ys)
    where 
      --Local definitions.--
      --variantsfiltered -> List containing only variants within TSS
      --specified by user, or within 2 Kb of TSS otherwise.
      variantsfiltered = DL.filter (\(_,_,c) -> c == 'Y') x
      --variantsWithinAmbiguityCodesAndTSSSmall -> This function will
      --grab all filtered ambiguity codes for matching genes.
      variantsWithinAmbiguityCodesAndTSSSmall :: [([String],[String],Char)] -> [[[String]]] -> [[String]]
      variantsWithinAmbiguityCodesAndTSSSmall []     [] = []
      variantsWithinAmbiguityCodesAndTSSSmall _      [] = []
      variantsWithinAmbiguityCodesAndTSSSmall []     _  = []
      variantsWithinAmbiguityCodesAndTSSSmall (x:xs) ys = (variantsAmbiguityCodesChecker x ambiguitycodesregionsfiltered) ++ (variantsWithinAmbiguityCodesAndTSSSmall xs ys)  
          where 
              --Local definitions.--
              ambiguitycodesregionsfiltered = DL.map (DL.filter (\y -> (y DL.!! 5) == ((\(a,_,_) -> a) (x) DL.!! 1))) ys
      --variantsAmbiguityCodesChecker -> This function will
      --call variantsAmbiguityCodesCheckerSmall.
      variantsAmbiguityCodesChecker :: ([String],[String],Char) -> [[[String]]] -> [[String]]
      variantsAmbiguityCodesChecker _      []     = []
      variantsAmbiguityCodesChecker xs     (y:ys) = (variantsAmbiguityCodesCheckerSmall xs y) ++ (variantsAmbiguityCodesChecker xs ys)
      --variantsAmbiguityCodesCheckerSmall -> This function will
      --call variantsAmbiguityCodesCheckerSmaller.
      variantsAmbiguityCodesCheckerSmall :: ([String],[String],Char) -> [[String]] -> [[String]]
      variantsAmbiguityCodesCheckerSmall _  []     = []
      variantsAmbiguityCodesCheckerSmall xs (y:ys) = if not (DL.null (variantsAmbiguityCodesCheckerSmaller xs (DL.take ((DL.length y) - 6) (DL.drop 6 y)) (DL.length (y DL.!! 1))))
                                                         then [[DL.intercalate ":" ((\(a,_,_) -> a) xs)] 
                                                            ++ [DL.intercalate ":" ((\(_,b,_) -> b) xs)] 
                                                            ++ [[((\(_,_,c) -> c) xs)]]
                                                            ++ [y DL.!! 0]
                                                            ++ [y DL.!! 1]
                                                            ++ [DL.intercalate "," (variantsAmbiguityCodesCheckerSmaller xs (DL.take ((DL.length y) - 6) (DL.drop 6 y)) (DL.length (y DL.!! 1)))]]
                                                           ++ (variantsAmbiguityCodesCheckerSmall xs ys)
                                                         else variantsAmbiguityCodesCheckerSmall xs ys
      --variantsAmbiguityCodesCheckerSmaller -> This function will
      --check whether the variant in question lies within an
      --ambiguity code sequence.
      variantsAmbiguityCodesCheckerSmaller :: ([String],[String],Char) -> [String] -> Int -> [String]
      variantsAmbiguityCodesCheckerSmaller _  []     _ = []
      variantsAmbiguityCodesCheckerSmaller xs (y:ys) z = --TSS reads in reverse direction (-1).
                                                         if (((\(_,b,_) -> b) xs) DL.!! 2) == "-1"
                                                             then if (read y :: Int) >= (read (((\(a,_,_) -> a) xs) DL.!! 3) :: Int) &&
                                                                     (read (((\(a,_,_) -> a) xs) DL.!! 3) :: Int) >= ((((read y) - z) + 1) :: Int)
                                                                 then [y] ++ (variantsAmbiguityCodesCheckerSmaller xs ys z)
                                                                 else variantsAmbiguityCodesCheckerSmaller xs ys z
                                                         --TSS reads in forward direction (1).
                                                         else if (read y :: Int) <= (read (((\(a,_,_) -> a) xs) DL.!! 3) :: Int) &&
                                                                 (read (((\(a,_,_) -> a) xs) DL.!! 3) :: Int) <= ((((read y) + z) - 1) :: Int) 
                                                             then [y] ++ (variantsAmbiguityCodesCheckerSmaller xs ys z)
                                                             else variantsAmbiguityCodesCheckerSmaller xs ys z
              ---------------------- 
      ----------------------

{--------------------}

{-Printing function.-}

--printFile -> This function will
--print the file to either stdout
--or to a output file based on
--command-lines options provided.
printFile :: [Flag] -> String -> [[String]] -> IO ()
printFile []    []       []    = return ()
printFile []    []       (_:_) = return ()
printFile (_:_) []       _     = return ()
printFile opts  filename xs    = do
    --Grab just "OUTFILE".
    let outdir = DL.head (DL.filter (isOutputDirectory) opts)
    --Extract the string from OutputDirectory.
    let outdirstring = extractOutputDirectory outdir
    --mapNotLast tabs and newlines in xs. 
    let tabsandnewlinesadded = DL.map (mapNotLast (++ "\t")) xs
    --Write the output file to the user-specified output directory.
    SIO.writeFile (outdirstring ++ filename) $
                  (TPB.render $
                  (TPB.hsep 0 TPB.left . DL.map (TPB.vcat TPB.left) . DL.map (DL.map (TPB.text)))
                  (DL.transpose tabsandnewlinesadded))

{---------------------}


{-FRI Specific Functions.-}

--processArgsAndFiles -> This function will
--walk through each of the command-line
--arguments and files provided by the user.
processArgsAndFiles :: ([Flag],[String]) -> IO ()
processArgsAndFiles ([],[]) = return ()
processArgsAndFiles (options,files) = do
    --Process the fasta file first.
    --Read in the data file.
    readfastafile <- BSF.readFasta (files DL.!! 0)  
    --Process the region file second.
    --Read in the region file.
    readregionfile <- SIO.readFile (files DL.!! 1)
    --Process readregionfile with lineFeed.
    let regions = lineFeed readregionfile
    --Remove header line from regions.
    let regionsnoheader = DL.tail regions 
    --Process the variant file third.
    --Read in the variant file.
    readvariantfile <- SIO.readFile (files DL.!! 2)
    --Process readvariantfile with lineFeed.
    let variants = lineFeed readvariantfile
    --Remove header line from variants.
    let variantsnoheader = DL.tail variants 
    --Process the ambiguity codes fourth.
    --Read in the ambiguity codes string.
    let ambiguitycodesstring = files DL.!! 3
    --Grab just the ambiguity codes, no semicolons.
    let ambiguitycodesfinal = (DL.filter (\x -> not (DL.null x)) (DLS.splitOn ";" ambiguitycodesstring))
    --Create list of tuples defining directionality of ambiguitycodesfinal.
    let ambiguitycodesfinaltuple = DL.map (\x -> (x,"1")) ambiguitycodesfinal 
    --Determine whether each variant is within
    --the TSS of the genes each variant is located in.
    let withintss = variantWithinRegionCheck options variantsnoheader regionsnoheader 
    --Grab the reverse complement of the 
    --user defined ambiguity codes.
    let ambiguitycodesreversecomplements = ambiguityCodesReverseComplement ambiguitycodesfinal
    --Create list of tuples defining directionality of ambiguitycodesreversecomplements.
    let ambiguitycodesreversecomplementstuple = DL.map (\x -> (x,"-1")) ambiguitycodesreversecomplements 
    --Grab all possible strings created from each ambiguity codes.
    allmappedambiguitystrs <- allStrGeneration (ambiguitycodesfinal ++ ambiguitycodesreversecomplements) 
    --Prepare allmappedambiguitystrs for ambiguityCodesWithinRegionCheck.
    let allmappedambiguitystrstuple = stringToTuple allmappedambiguitystrs 
    --Determine whether there are ambiguity codes strings
    --present within the TSS of each region (PARALLELIZED). 
    let ambiguitycodeswithintss = (ambiguityCodesWithinRegionCheck (ambiguitycodesfinaltuple ++ ambiguitycodesreversecomplementstuple) allmappedambiguitystrstuple regionsnoheader readfastafile options) `CPS.using` (CPS.parList CPS.rdeepseq)
    --Prepare ambiguitycodeswithintss for printing.
    let analysisreadyambiguitycodeswithintss = prepareAmbiguityCodesWithinTSS ambiguitycodeswithintss 
    --Determine whether there are variants present
    --within ambiguity codes within corresponding regions.
    let variantsinambiguitycodesandtss = ambiguitycodeswithintss `CD.deepseq` variantsWithinAmbiguityCodesAndTSS (DL.filter (\x -> not (DL.null x)) withintss) analysisreadyambiguitycodeswithintss  
    --Prepare withintss, ambiguitycodeswithintss, and variantsinambiguitycodesandtss for printing. 
    let printreadywithintss = prepareWithinTSS withintss
    let printreadyambiguitycodeswithintss = ambiguitycodeswithintss `CD.deepseq` DL.concat (DL.map (DL.map (\xs -> (DL.take 6 xs) ++ [DL.intercalate "," (DL.drop 6 xs)])) (prepareAmbiguityCodesWithinTSS ambiguitycodeswithintss)) 
    let printreadyvariantsinambiguitycodesandtss = variantsinambiguitycodesandtss
    --Prepare final print ready files with headers.
    let finalprintreadywithintss = [["Variant","Region","Variant_Within_Region"]] ++ printreadywithintss
    let finalprintreadyambiguitycodeswithintss = [["Ambiguity_Code","Mapped_Nucleotide_String","Chromosome","TSS","Strand","SYMBOL","Ambiguity_Code_String_Locations_Within_TSS"]] ++ printreadyambiguitycodeswithintss
    let finalprintreadyvariantsinambiguitycodesandtss = [["Variant","Region","Variant_Within_Region","Ambiguity_Code","Mapped_Nucleotide_String","Ambiguity_Code_String_Locations_Within_TSS"]] ++ printreadyvariantsinambiguitycodesandtss
    --Print  withintss, ambiguitycodeswithintss, and variantsinambiguitycodesandtss to files. 
    finalprintreadywithintss `CD.deepseq` printFile options "variants.tsv" finalprintreadywithintss
    finalprintreadyambiguitycodeswithintss `CD.deepseq` printFile options "ambiguity_codes.tsv" finalprintreadyambiguitycodeswithintss
    finalprintreadyvariantsinambiguitycodesandtss `CD.deepseq` printFile options "variants_in_ambiguity_codes.tsv" finalprintreadyvariantsinambiguitycodesandtss 

{-------------------------}


{-Main function.-}

main :: IO ()
main = do
    --Get command line arguments.
    (args,files) <- SE.getArgs >>= compilerOpts
    --Run args and files through processArgsandFiles.
    processArgsAndFiles (args,files)

{----------------}
