create_2d_funnel_plot <- function(xmin,xmax,ymin,ymax) {

  compute_U <- function(q1,q2) {
    y <- q1
    x <- q2
    (y^2)/18 + y/2 + (x^2)/(2*exp(y))
  }
  Uvec <- Vectorize(function(q1, q2) compute_U(q1,q2))
  expand.grid(q2 = seq(xmin,xmax,(xmax-xmin)/200), q1 = seq(ymin,ymax,(ymax-ymin)/200)) %>%
    as_tibble %>%
    mutate(U = Uvec(q1,q2)) %>%
    ggplot(aes(q2,q1)) +
    geom_raster(aes(fill = log(U+2))) +
    scale_fill_gradientn(colours = c("white", "orange", "black"))

}

plot_funnel_nuts_draw <- function(funnel_plot, nuts_draw_hist) {

  origin <- nuts_draw_hist %>% filter(is.na(depth))
  left_most_point <- head(nuts_draw_hist,1)
  right_most_point <- tail(nuts_draw_hist,1)

  # left_vector <- tibble(q0.1 = origin$q1, q0.2 = origin$q2, q1.1 = left.most.point$q1, q1.2 = left.most.point$q2)
  # right_vector <- tibble(q0.1 = origin$q1, q0.2 = origin$q2, q1.1 = right.most.point$q1, q1.2 = right.most.point$q2)
  # uturn <- nuts.draw.hist %>% filter(invalid == "U-Turn")
  # head.uturn <- head(uturn,1)
  # tail.uturn <- tail(uturn,1)

  funnel_plot +
    geom_path(data = nuts_draw_hist) +
    geom_point(size = 1.7, data = nuts_draw_hist %>% mutate(depth = as.integer(depth))) +
    geom_point(aes(color=as.factor(depth)), size = 1.2,
               data = nuts_draw_hist %>% mutate(depth = as.integer(depth)))
    # geom_point(color = "red", size = 3.0, data = nuts.draw.hist %>% filter(is.na(depth))) +
    # geom_point(color = "green",  size = 3.0,data = tibble(q1 = q[1], q2 = q[2]))
  # +
  # geom_segment(aes(xend = q2+p2, yend = q1+p1), color = "red",
  #              arrow = arrow(length = unit(0.01, "npc")),
  #              nuts.draw.hist %>% filter(is.na(depth)))
}
