#--source("generics.R");

#'
#' @title Create a ragged array (an object of class "ragged_array")
#' @description Helper function to create a ragged array (an object of class "ragged_array").
#'
#' @param x - object coercible to a vector, or a ragged_array
#' @param dfrDims - tibble defining dimensions (optional if `x` is itself a ragged_array)
#'
#' @details A ragged_array object is essentially a vector with an attribute (a tibble)
#' that maps the index of the vector to a set of dimensions. `x` is coerced to a vector
#' (using `as.vector(x)`) and truncated or replicated enough to match the number of rows in
#' `dfrDims`. The columns of dfrDims then constitute the values of the dimension indices
#' for the resulting ragged_array object and a row indicates the index of each dimension
#' associated with the corresponding value of the vector.
#'
#' If `x` is a ragged_array object, then the tibble `dfrDims` is used to subset `x` based on
#' the latter's `dfrDims` attribute.
#'
#' This is a "helper" function (see https://adv-r.hadley.nz/s3.html#helpers) for package users.
#'
#' @export
#'
ragged_array<-function(x,
                       dfrDims=NULL){
  if (inherits(x,"ragged_array")){
    if (!is.null(dfrDims)){
      #--extract elements only in dfrDims
      dmnms = attr(dfrDims,"dmnms");
      if (is.null(dmnms)) dmnms=names(dfrDims);
      dfrDimsX = attr(x,"dfrDims");
      dfrDimsY = dplyr::right_join(dfrDimsX,
                                  dfrDims |> dplyr::select(dplyr::all_of(dmnms)),
                                  by=dmnms);
      rws = dfrDimsY[[1]];
      dfrDimsY[[1]] = 1:length(rws);
      y = new_ragged_array(x[rws],dfrDimsY);
    } else {
      y = x;
    }
  } else {
    if (!is.null(dfrDims)){
      if (is.vector(x)) {
        y = new_ragged_array(x,dfrDims);
      }
    } else {
      if (is.vector(x)) {
        dfrDims = createSparseDimsMap(i=1:length(x));
        y = new_ragged_array(x,dfrDims);
      }
      warning("Creating a ragged_array object from a ",class(x)[1]," object but dfrDims is NULL.")
    }
  }
  return(y);
}

#'
#' @title Create a new ragged_array object
#' @description Internal constructor to create a new ragged_array object.
#'
#' @param x - vector (or other object that can be converted to a vector using `as.vector(x)`)
#' @param dfrDims - dimensions dataframe (a DimsMap)
#'
#' @return a vector of class "ragged_array" with a "dfrDims" attribute
#'
#' @details A ragged_array object is essentially a vector with an attribute (a tibble)
#' that maps the index of the vector to a set of dimensions.
#' If `length(as.vector(x))` is < `nrow(dfrDims)`, it is replicated using `rep_len(as.vector(x),nrow(dfrDims))`
#' to a length of `nrow(dfrDims)`. If `length(x)` is > `nrow(dfrDims)`, it is truncated to the
#' latter.
#'
#' This function is not exported. Use `ragged_array(x,dfrDims)` instead.
#'
#' See https://adv-r.hadley.nz/s3.html#s3-constructor for details on constructors.
#'
new_ragged_array<-function(x,dfrDims){
  stopifnot(is.data.frame(dfrDims));
  n = nrow(dfrDims);
  y = as.vector(x);
  if (n<length(y)) y = y[1:n];
  if (n>length(y)) y = rep_len(y,length.out=n);
  ra = structure(y,dfrDims=dfrDims,class=unique(c("ragged_array",class(x))));
  return(ra);
}

#'
#' @title Subset a ragged_array object
#' @description Subset a ragged_array object by its vector index.
#' @param x - the ragged array object
#' @param i - a vector of indices by which to subset `x`
#' @return the subsetted ragged_array object
#' @details The "dfrDims" attribute for the returned ragged_array object is
#' the "dfrDims" attribute of `x`, but appropriately subsetted based on its row index
#' to match the subsetted vector underlying `x`.
#'
#' @examples
#' # simple subsetting
#' x = ragged_array(1:5);
#' x[2:3];
#'
#'
#' @export
#'
`[.ragged_array`<-function(x,i=NULL){
  if (is.null(i)) i = 1:length(x);
  dfrDims = attr(x,"dfrDims") |>
              dplyr::filter(dplyr::row_number() %in% i);
  dfrDims = dfrDims[i,];
  ra = new_ragged_array(unclass(x)[i],dfrDims);
  return(ra);
}

#'
#' @title Assign values to a ragged_array object
#' @description Assign values to a ragged_array object.
#'
#' @param x - the LHS ragged array object in the assignment
#' @param y - the indices into x, or a dimsMap
#' @param v - the values to assign
#'
#' @return x, the LHS object
#'
#' @details TODO!
#'
#' @export
#'
`[<-.ragged_array`<-function(v){
  cat("#--in [<-.ragged_array\n")
  if (is.DimsMap(y)) {
    #--replace values per dim values of y
    ##--TODO: complete this!!----
  } else if (is.vector(y)){
    cat("#--y is vector\n")
    `*tmp*` = v;
  }
}

#'
#' @title Test if an object is a ragged_array
#' @description Function to test if an object is a ragged_array.
#' @param x - the object
#' @return TRUE or FALSE
#' @details Returns TRUE of the object inherits from class "ragged_array".
#' @export
#'
`is.ragged_array`<-function(x){
  return(inherits(x,"ragged_array"));
}

#'
#' @title Display structure of a ragged_array
#' @description Function to display structure of a ragged_array.
#' @param x - the object
#' @return TRUE or FALSE
#' @details Returns TRUE of the object inherits from class "ragged_array".
#' @export
#'
`str.ragged_array`<-function(x,...){
  cat("# A ragged_array: \n")
  str(as_tibble.ragged_array(x));
}

#'
#' @title Convert a ragged_array to a vector
#' @description Function to convert a ragged_array to a vector.
#' @param x - the ragged_array
#' @param mode - mode for resulting vector (default="any")
#' @return vector of specified mode (default mode is the same as that of the underlying vector)
#' @details Returns the vector underlying the ragged_array without the "dfrDims" attribute.
#' This is an S3 generic method, so just ue `as.vector(x)`
#' @export
#'
as.vector.ragged_array<-function(x,mode="any"){
  x = unclass(x);
  attr(x,"dfrDims")<-NULL;
  return(x);
}

#'
#' @title print a ragged_array
#' @description S3 method to print a ragged_array (as a tibble).
#' @param x - the ragged array to print
#' @return none
#' @details Prints the ragged array as a tibble with the values of the vector
#' in the "val" column. This is an S3 generic method, so just use `print(x)`.
#' @export
#'
`print.ragged_array`<-function(x){
  print(as_tibble.ragged_array(x));
}

#'
#' @title Convert a ragged_array to a [tibble::tibble]
#' @description Function to convert a ragged_array to a [tibble::tibble].
#' @param x - the ragged array to convert
#' @return a tibble
#' @details The ragged array is converted to a tibble with the values of the vector
#' in the "val" column. Other columns correspond to the row index and the associated
#' dimension indices. This is *not* an S3 method, so use `as_tibble.ragged_array(x)`.
#' @export
#'
`as_tibble.ragged_array`<-function(x){
  if (!is.ragged_array(x)){
    stop("x is not a ragged array: has class(es) ",paste(class(x)), collapse=", ");
  }
  return(dplyr::bind_cols(attr(x,"dfrDims"),val=unclass(x)));
}

#'
#' @title Expand a ragged_array by a dimensions map
#' @description S3 method to expand a ragged_array using a dimensions map.
#' @param x - the ragged array to convert
#' @param y - the dimensions map used to expand `x`
#' @return a ragged_array
#' @details `x`'s dimension map (it's "dfrDims" attribute) is expanded by `y` using
#' expand.DimsMap, then the vector underlying `x` is suitably replicated across
#' the dimensions in `y`.
#'
#' This is an S3 generic method, so just use `expand(x,y)`.
#' @export
#'
`expand.ragged_array`<-function(x,y){
  if (!is.ragged_array(x)){
    stop("x is not a ragged array: has class(es) ",paste(class(x)), collapse=", ");
  }
  dfrDims = y;
  dfrDimsA = attr(x,"dfrDims");
  dmnmsA   = attr(dfrDimsA,"dmnms");
  dmnmsB   = attr(dfrDims,"dmnms");
  dmnmsC   = unique(c(dmnmsA,dmnmsB));
  dfrDimsC = expand.DimsMap(dfrDimsA,dfrDims);
  z = as.vector(x)[dfrDimsC$sparse_idx];
  dfrDimsC$sparse_idx = 1:nrow(dfrDimsC);
  z = new_ragged_array(unclass(z),dfrDimsC);
  return(z);
}

getDimIndices<-function(x){
  if (is.ragged_array(x)){
      dfrDimsX = attr(x,"dfrDims") |>
                   dplyr::select(!1);
      return(dfrDimsX);
  }else{
    warning("'x' is not a ragged_array: ",str(x))}
}




