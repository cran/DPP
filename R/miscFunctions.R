### This place is for miscellaneous functions we want
### to make available once the library is loaded

#' Calculate the expected number of clusters from the number of individuals and a concentration parameter
#'
#' @param n the number of individuals
#' @param a the concentration parameter
#' @return the expected number of clusters
#' @examples
#' expectedNumberOfClusters(100,0.2)
#' expectedNumberOfClusters(100,0.15)
expectedNumberOfClusters = function(n,a) sum(a / (a + 1:n - 1))

#' Simulate a discrete distribution as in the chinese restaurant problem
#'
#' @param num_elements the number of elements to be grouped
#' @param chi the concentration parameter
#' @return The sum of \code{x} and \code{y}
#' @examples
#' simulateChineseRestaurant(100, 0.2)
#' simulateChineseRestaurant(100, 0.8)
simulateChineseRestaurant = function(num_elements,chi) {

  allocation <- numeric(num_elements)
  allocation[1] <- 1
  num_cats <- 1
  num_in_cats <- 1

  for(i in 2:num_elements) {

    # update the number of elements in each category
    num_in_cats <- table(allocation[1:(i-1)])

    u <- runif(1,0,1)
    if ( u < chi / (chi + i - 1) ) {
      # choose a new category
      num_cats <- num_cats + 1
      this_category <- num_cats
    } else {
      # choose an existing category
      this_category <- sample.int(n=num_cats, size=1, replace=TRUE, prob=num_in_cats / (chi + i - 1))
    }

    allocation[i] <- this_category

  }

  return(allocation)

}
