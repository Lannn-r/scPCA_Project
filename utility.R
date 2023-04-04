#' Splice in default arguments into a function
#'
#' Arguments in ... supercede in case of collisions with `extra`
#' @param f `function`
#' @param ... key-value (named) arguments
#' @param extra named `list` of default arguments
#'
#' @return  value of `f` with supplied arguments
#' @export
#'
#' @examples
#' call_intercalate(sum, 3, 4, NA, extra = list(na.rm = TRUE))
#' call_intercalate_left(sum, 3, NA, na.rm = FALSE, extra = list(na.rm = TRUE))
#' call_intercalate_right(sum, 3, NA, na.rm = FALSE, extra = list(na.rm = TRUE))
#' meld_list_left(list(A=1, B=2), list(A = 0))
call_intercalate = function(f, ..., extra){
  nargs = meld_list_left(list(...), extra)
  if(length(nargs) != (length(list(...)) + length(extra))) warning("Duplicated arguments")
  do.call(f, nargs)
}

#' @describeIn call_intercalate don't warn with collision
#' @export
call_intercalate_left = function(f, ..., extra){
  nargs = meld_list_left(list(...), extra)
  do.call(f, nargs)
}

#' @describeIn call_intercalate arguments in `extra` take presidence
#' @export
call_intercalate_right = function(f, ..., extra){
  nargs = meld_list_left(extra, list(...))
  do.call(f, nargs)
}

#' @describeIn call_intercalate combine lists, preferentially taking elements from x if there are duplicate names
#' @param x list
#' @param y list
#' @export
meld_list_left = function(x, y){
  unite = c(x, y)
  dups = nchar(names(unite)) & duplicated(names(unite))
  unite[!dups]
}
