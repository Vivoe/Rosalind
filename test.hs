module Main where

main :: IO()
main = print $ reading ["34", "2", "323"]

addDiv a b = (a+b)/2

reading list = map (read) list :: [Int]