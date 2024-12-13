# guides fix for patchwork
# source: https://github.com/thomasp85/patchwork/issues/404
guides_build_mod <- function (guides, theme){
  legend.spacing.y <- calc_element("legend.spacing.y", theme)  # modified by me
  legend.spacing.x <- calc_element("legend.spacing.x", theme)  # modified by me
  legend.box.margin <- calc_element("legend.box.margin", theme) %||% 
    margin()
  widths <- exec(unit.c, !!!lapply(guides, gtable_width))
  heights <- exec(unit.c, !!!lapply(guides, gtable_height))
  just <- valid.just(calc_element("legend.box.just", theme))
  xjust <- just[1]
  yjust <- just[2]
  vert <- identical(calc_element("legend.box", theme), "horizontal")
  guides <- lapply(guides, function(g) {
    editGrob(g, vp = viewport(x = xjust, y = yjust, just = c(xjust, 
                                                             yjust), height = if (vert) 
                                                               heightDetails(g)
                              else 1, width = if (!vert) 
                                widthDetails(g)
                              else 1))
  })
  guide_ind <- seq(by = 2, length.out = length(guides))
  sep_ind <- seq(2, by = 2, length.out = length(guides) - 1)
  if (vert) {
    heights <- max(heights)
    if (length(widths) != 1) {
      w <- unit(rep_len(0, length(widths) * 2 - 1), "mm")
      w[guide_ind] <- widths
      w[sep_ind] <- legend.spacing.x
      widths <- w
    }
  }
  else {
    widths <- max(widths)
    if (length(heights) != 1) {
      h <- unit(rep_len(0, length(heights) * 2 - 1), "mm")
      h[guide_ind] <- heights
      h[sep_ind] <- legend.spacing.y
      heights <- h
    }
  }
  widths <- unit.c(legend.box.margin[4], widths, legend.box.margin[2])
  heights <- unit.c(legend.box.margin[1], heights, legend.box.margin[3])
  guides <- gtable_add_grob(gtable(widths, heights, name = "guide-box"), 
                            guides, t = 1 + if (!vert) 
                              guide_ind
                            else 1, l = 1 + if (vert) 
                              guide_ind
                            else 1, name = "guides")
  gtable_add_grob(guides, element_render(theme, "legend.box.background"), 
                  t = 1, l = 1, b = -1, r = -1, z = -Inf, clip = "off", 
                  name = "legend.box.background")
}

environment(guides_build_mod) <- asNamespace('patchwork')
assignInNamespace("guides_build", guides_build_mod, ns = "patchwork")