module Main where

import Control.Monad
import qualified Data.Map.Strict as M
import qualified Data.List as L
import Data.List.Split
import qualified Data.Vector.Unboxed as V
import Debug.Trace

main :: IO()
main = rosalindSubs

rosalindSubs :: IO()
rosalindSubs = do
	content <- fmap (lines) $ readFile "rosalind_subs.txt"
	print $ unwordList $ stringIndexes (content!!1) (content!!0)

stringIndexes :: String -> String -> [Int]
stringIndexes target dna = reverse $ find dna 1 []
	where find [] _ found = found
	      find list index found
	      	|take len list == target = find (tail list) (index + 1) (index:found)
	      	|otherwise = find (tail list) (index + 1) found
	      len = length target

rosalindProt :: IO()
rosalindProt = do
	content <- readFile "rosalind_prot.txt"
	table <- codonTable
	print $ translation content table

translation :: String -> M.Map String String -> String
translation rna table = concat $ takeWhile (/= "Stop") $ fmap (\key -> M.findWithDefault "X" key table) (chunksOf 3 rna)

codonTable :: IO (M.Map String String)
codonTable = do
	content <- fmap (convert) $ readFile "codonTable.txt"
	return content
	where convert a = M.fromList $ fmap toTuple $ ((fmap words) . lines) a
	      toTuple x = (head x, last x)

rosalindIPrb :: IO()
rosalindIPrb = do
	content <- fmap convert $ readFile "rosalind_iprb.txt"
	print $ domProbability (content!!0) (content!!1) (content!!2)
	where convert a = ((fmap read) . words) a :: [Double]

domProbability :: Double -> Double -> Double -> Double
domProbability k m n = 1 - (0.25 * m * (m-1) + m * n + n * (n-1))/(t * (t-1))
	where t = k + m + n

rosalindHamm :: IO()
rosalindHamm = do
	content <- fmap lines $ readFile "rosalind_hamm.txt"
	print $ numPointMutations (content!!0) (content!!1)

numPointMutations :: String -> String -> Int
numPointMutations a b = (length . filter (== False)) $ L.zipWith (==) a b

formatRaw :: [String] -> [String]
formatRaw raw = map L.concat $  L.groupBy (\a b -> if head a /= '>'  && head b /= '>' then True else False) raw

rosalindGC :: IO()
rosalindGC = do
	content <- fmap convert $ readFile "rosalind_gc.txt"
	print $ (maxGCContent . packageData) content
	where convert a = (formatRaw . words) a :: [String]

packageData :: [String] -> [(String, String)]
packageData n = fmap toTuple $ L.groupBy (\a b -> if head a == '>'  && head b /= '>' then True else False) n
	where toTuple x = (head x, last x)

maxGCContent :: [(String, String)] -> (String, Double)
maxGCContent fasta = head $ L.sortBy sortData (fmap (\(dnaId, dna) -> (dnaId, gcContent dna)) fasta)

sortData :: (String, Double) -> (String, Double) -> Ordering
sortData a b = if (snd a >= snd b) then LT else GT

gcContent :: String -> Double
gcContent dna = (snd count / fst count) * 100
	where count = foldr (\x (len, gc) -> if x == 'G' || x == 'C' then (len + 1, gc + 1) else (len + 1, gc)) (0,0) dna

rosalindRabbits :: IO()
rosalindRabbits = do
	content <- fmap convert $ readFile "rosalind_fib.txt"
	print $ rabbits (head content) (last content)
	where convert a = ((fmap read) . words) a :: [Integer]

rabbits :: Integer -> Integer -> Integer
rabbits n k = rabbits' (0, 1) 1
    where rabbits' (old, young) r
            |r >= n = old + young
            |otherwise = rabbits' (old + young, old * k) (r+1)

rosalindIns :: IO()
rosalindIns = do
	content <- fmap convert $ readFile "rosalind_ins.txt"
	print $ snd $ insertionSort (tail content, 0)
	where convert a = ((fmap read) . words) a :: [Int]

insertionSort :: ([Int], Int) -> ([Int], Int)
insertionSort ([], _) = ([], 0)
insertionSort ((x:xs), n) = insert x (insertionSort (xs, n))

insert :: Int -> ([Int], Int) -> ([Int], Int)
insert x ([], n) = ([x], n)
insert x ((y:ys), n)
  | x <= y = ((x:y:ys), n)
  | otherwise = (y:(fst result), snd result)
  where result = insert x (ys, n+1)

rosalindGraphDeg :: IO()
rosalindGraphDeg = do
	content <- fmap (convert) $ readFile "rosalind_deg.txt"
	print . unwordList $ map length $ (L.group . L.sort) (tail $ tail content)
	where convert a = ((fmap read) . words) a :: [Int]

rosalindMajElem :: IO()
rosalindMajElem = do
	content <- fmap (convert . lines) $ readFile "rosalind_maj.txt"
	print . unwordList $ fmap (\list -> majElement list (content!!0!!1)) (tail content)
	where convert a = fmap ((fmap read) . words) a :: [[Int]]

majElement :: [Int] -> Int -> Int
majElement list size = if checkMaj sorted size point then point else -1
    where sorted = L.sort list
          point = sorted!!(size `div` 2)

checkMaj list size num = if (foldr (\x acc-> if x == num then acc + 1 else acc) 0 list) >  ceiling (fromIntegral size / 2) then True else False

unwordList :: Show a => [a] -> String
unwordList a = unwords $ map show a

rosalindSearch :: IO()
rosalindSearch = do
	content <- fmap (convert . lines) $ readFile "rosalind_bins.txt"
	print . unwordList $ fmap (binarySearch (V.fromList (content!!2)) (head $ head content)) (content!!3)
	where convert a = fmap ((fmap read) . words) a :: [[Int]]

binarySearch :: V.Vector Int -> Int -> Int -> Int
binarySearch list size num = find list 0 size num

find :: V.Vector Int -> Int -> Int -> Int -> Int
find list start end num
	|start == end = -1
	| pivotValue == num = pivot + 1 
	| pivotValue > num = find list start (pivot) num
	| pivotValue < num = find list (pivot+1) end num
	where pivot = (start + end) `div` 2
	      pivotValue = list V.! pivot

rosalindFib :: IO()
rosalindFib = (fmap (fib . read)) (readFile "rosalind_fibo.txt") >>= print

fib 1 = 1
fib 0 = 0
fib n = fib (n-1) + fib (n-2)
