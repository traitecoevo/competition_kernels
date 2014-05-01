## Hook to replace ```r -> ```S in generated output -- this renders
## better on github.  I don't have the patience to go through and work
## out which bits need to be escaped to have this work from within
## .knit.sh, but it seems like there's a bunch of trouble makers here
## (`, $, ').
knitr::knit_hooks$set(source=function(x, options)
                      paste0('\n\n```S\n', x, '\n```\n\n'))

## Hook to make more friendly (smaller) figure margins and set the
## hook to run by by default.  Pass small_mar=FALSE to disable.
knitr::knit_hooks$set(small_mar=function(before, options, envir) {
  if (before) par(mar=c(4, 4, .1, .1)) # smaller margin on top and right
})
knitr::opts_chunk$set(small_mar=TRUE)
