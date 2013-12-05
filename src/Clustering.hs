{-------------------------------------------------------------------------------
 Clustering

 (c) 2009 Jan Snajder <jan.snajder@fer.hr>

-------------------------------------------------------------------------------}

module Clustering where

type DistMeasure a = a -> a -> Double
type DistMeasureM m a = a -> a -> m Double

