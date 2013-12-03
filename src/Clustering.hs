{-------------------------------------------------------------------------------
 Clustering

 (c) 2009 Jan Snajder <jan.snajder@fer.hr>

-------------------------------------------------------------------------------}

module Clustering (DistanceMeasure) where

type DistanceMeasure a = a -> a -> Double

