{-------------------------------------------------------------------------------
 Clustering.Hierarchical
 Hierarchical agglomerative clustering.

 (c) 2009 Jan Snajder <jan.snajder@fer.hr>

-------------------------------------------------------------------------------}

-- | Hierarchical agglomerative clustering.
module Clustering.Hierarchical
  (DistanceMeasure,
   Linkage(..),
   Dendrogram,
   cluster,
   distanceCut,
   numberCut,
   clusterAt,
   clusterInto,
   distances,
   scale,
   normalize,
   clusterPartition,
   clusterPartitionCb) where

import Clustering

import Data.List (insertBy,sortBy,partition)
import Data.Maybe
import Control.Monad (liftM,zipWithM)
import Data.Ord (comparing)
import qualified Data.Map as M

type Distances    = [((Int,Int),Double)]
type Clusters a   = M.Map Int (Int, Dendrogram a)
data Linkage      = Single | Complete | Average 
  deriving (Eq,Show)
data Dendrogram a = Item a | Merge Double (Dendrogram a) (Dendrogram a) 
  deriving (Show,Read,Eq)

------------------------------------------------------------------------------
-- Operations on dendrograms
------------------------------------------------------------------------------

-- | Converts cluster map into a dendrogram.
toDendrogram :: Clusters a -> Dendrogram a
toDendrogram = snd . head . M.elems

-- | Cuts a dendrogram at a specified distance level.
distanceCut :: Double -> Dendrogram a -> [[a]]
distanceCut _ (Item x) = [[x]]
distanceCut t (Merge d l r)
  | d<=t      = merge (distanceCut t l ++ distanceCut t r)
  | otherwise = distanceCut t l ++ distanceCut t r
  where merge = (:[]) . concat

infty = 1/0

-- | Returns a reverse-sorted list of finite dendrogram distance levels,
-- ignoring the infinite distances.
distances :: Dendrogram a -> [Double]
distances = dropWhile (==infty) . sortBy (flip compare) . dist
  where dist (Merge d l r) = d : dist l ++ dist r
        dist _ = []

-- | Scales dendrogram distances by a specified factor, 
-- thereby ignoring the infinite distances.
scale :: Double -> Dendrogram a -> Dendrogram a
scale _ i@(Item _ )      = i
scale x dg@(Merge d l r) = Merge (d/x) (scale x l) (scale x r)

-- | Normalizes dendrogram distances to [0,1] interval,
-- thereby ignoring the infinite distances.
-- Defined as:
-- @ normalize dg = scale (head $ distances dg) dg @
normalize :: Dendrogram a -> Dendrogram a
normalize dg = scale (head $ distances dg) dg

-- | Cuts a dendrogram at a specified number of clusters.
numberCut :: Int -> Dendrogram a -> [[a]]
numberCut n dg = distanceCut ((distances dg ++ repeat 0.0)!!n) dg

-- | Merges two dendrograms (merging distance is set to infinity).
merge :: Dendrogram a -> Dendrogram a -> Dendrogram a
merge dg1 dg2 = Merge infty dg1 dg2

-- | Concatenates a list of dendrograms (undefined for empty list).
concatDendrograms :: [Dendrogram a] -> Dendrogram a
concatDendrograms []  = error "empty list"
concatDendrograms dgs = foldr1 merge dgs

------------------------------------------------------------------------------
-- Clustering algorithm
------------------------------------------------------------------------------

-- | Creates initial clusters from a list of items.
initClusters :: [a] -> Clusters a
initClusters xs = M.fromList . zip [0..] $ [(1,Item x) | x <- xs]

-- | Merges two clusters in the map of clusters.
mergeClusters :: Clusters a -> ((Int,Int),Double) -> Clusters a
mergeClusters cs ((i,j),d) = 
  M.delete j . M.update (\(ni,ci) -> Just (ni+nj,Merge d ci cj)) i $ cs
  where Just (nj,cj) = M.lookup j cs

-- | Returns the initial distances list between singleton clusters.
initDist :: Clusters a -> DistanceMeasure a -> Distances
initDist cs dm = sortBy (comparing snd) dl
  where n = M.size cs
        dl = map calcDist [(i,j) | i <- [1..n-1], j <- [0..i-1]]
        calcDist ij@(i,j) = (ij,dm ei ej)
          where Just (_, Item ei) = M.lookup i cs
                Just (_, Item ej) = M.lookup j cs

-- | Returns the size of a specified cluster.
clusterSize :: Clusters a -> Int -> Int
clusterSize cs i = case M.lookup i cs of
                      Just (n,_) -> n

-- | Updates distance list after merging of clusters i and j.
-- (Distances to cluster i are updated, distances to cluster j are removed.)
updateDistances :: Clusters a -> Linkage -> Distances -> Distances
updateDistances cs l (((i,j),_):ds) = mergeBy (comparing snd) new rest
  where (aff,rest) = partition mustUpdate ds
        (del,upd) = partition mustDelete aff
        new = sortBy (comparing snd) . map update $ upd
        mustUpdate ((x,y),_) = x==i || x==j || y==i || y==j
        mustDelete ((x,y),_) = x==j || y==j
        update ((x,y),di) | l==Complete = ((x,y),max di dj)
                          | l==Single   = ((x,y),min di dj)
                          | l==Average  = ((x,y),(ni*di+nj*dj)/(ni+nj))
          where (Just dj) = lookup pair del
                pair = if (z>j) then (z,j) else (j,z)
                z = if x==i then y else x
                ni = fromInteger.toInteger $ clusterSize cs i
                nj = fromInteger.toInteger $ clusterSize cs j

-- | Merges two sorted lists using a user-specified comparison function.
mergeBy :: (Ord a) => (a -> a -> Ordering) -> [a] -> [a] -> [a]
mergeBy _ [] ys = ys
mergeBy _ xs [] = xs
mergeBy p (x:xs) (y:ys)
  | p x y == LT = x:mergeBy p xs (y:ys)
  | otherwise   = y:mergeBy p (x:xs) ys 

-- | Clusters the elements of a list using a given distance measure 
-- and linkage type. Returns a 'Dendrogram', which can be cut at 
-- arbitrary levels using 'distanceCut' and 'numberCut'.
cluster :: DistanceMeasure a -> Linkage -> [a] -> Dendrogram a
cluster dm l xs = toDendrogram $ cluster' cs ds
  where cs = initClusters xs
        ds = initDist cs dm
        cluster' cs [] = cs
        cluster' cs ds@(d:_) = 
          cluster' (mergeClusters cs d) (updateDistances cs l ds) 

------------------------------------------------------------------------------
-- Additional flavors
------------------------------------------------------------------------------

-- | Clusters until distances are above a specified threshold.
-- Defined as:
-- @ clusterAt t dm l = distanceCut t . cluster dm l @
clusterAt :: Double -> DistanceMeasure a -> Linkage -> [a] -> [[a]]
clusterAt t dm l = distanceCut t . cluster dm l

-- | Clusters into a specified number of clusters (if possible).
-- Defined as:
-- @ clusterInto n dm l = numberCut n . cluster dm l @
clusterInto :: Int -> DistanceMeasure a -> Linkage -> [a] -> [[a]]
clusterInto n dm l = numberCut n . cluster dm l

-- | Clusters pre-partitioned elements and produces a single
-- dendrogram.
clusterPartition :: DistanceMeasure a -> Linkage -> [[a]] -> Dendrogram a
clusterPartition dm l = concatDendrograms . map (cluster dm l)

-- | A call-back version of 'clusterPartition'. The call-back function is
-- invoked for each partition and the argument of the
-- function is the number of the partition being clustered.
clusterPartitionCb :: (Int -> IO t) -> DistanceMeasure a -> Linkage -> [[a]] -> IO (Dendrogram a)
clusterPartitionCb cb dm l xss = concatDendrograms `liftM` zipWithM step xss [1..]
  where step xs i = cb i >> let r = cluster dm l xs in r `seq` return r

