## Hook to replace ```r -> ```S in generated output -- this renders
## better on github.
knitr::knit_hooks$set(source=function(x, options)
                      paste0('\n\n```S\n', x, '\n```\n\n'))

## Hook to make more friendly (smaller) figure margins and set the
## hook to run by by default.  Pass small_mar=FALSE to disable.
knitr::knit_hooks$set(small_mar=function(before, options, envir) {
  if (before) par(mar=c(4, 4, .1, .1)) # smaller margin on top and right
})
knitr::opts_chunk$set(small_mar=TRUE)

local({
  knit_and_read <- function(filename) {
    readLines(knitr::knit(filename, tempfile(), quiet=TRUE))
  }
  document_with_footer <- function(x) {
    knitr::render_markdown()
    c(x, knit_and_read(".knitr_footer.Rmd"))
  }
  knitr::knit_hooks$set(document=document_with_footer)
})
