# A function to get the color codes of given subject id
# input: ids  Sample IDs
#        by   how to assign colour by day or subject
ids2colors = function(ids) {
  unids = unique(ids)
  col.unids = rainbow(length(unids))
  names(col.unids) = unids
  return(col.unids[ids])
}
