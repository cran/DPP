### This place is for miscellaneous functions we want
### to make available once the library is loaded

expectedNumberOfClusters = function(n,a) sum(a / (a + 1:n - 1))
