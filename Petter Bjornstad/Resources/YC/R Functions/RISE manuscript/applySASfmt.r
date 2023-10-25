applySASfmt <- function(x){
  if(!is.null(attr(x,'labels')) & !is.null(attr(attr(x, 'labels'),'names'))) x <- factor(x, levels = attr(x,'labels'), labels = attr(attr(x, 'labels'),'names'), exclude=NULL)
  return(x)
}
