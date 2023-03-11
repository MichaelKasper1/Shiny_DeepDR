#Date: 20230210
#Function: Calculate Jaccard Similarity


ShinyDeepDR_jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  
  return (intersection/union)
}

