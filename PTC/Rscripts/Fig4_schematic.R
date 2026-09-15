# Figure 4e: conceptual clonal and subclonal branching patterns.
# Geometry and colours illustrate the two architectures without inferring event order.
run_Fig4_schematic <- function(input_dir = NULL, output_dir) {
  ptc_font()
  # Convert drawing coordinates from a downward to an upward y axis.
  nodes <- data.frame(
    id = c('upper_root', 'upper_split', 'upper_long', 'upper_short',
           'lower_root', 'lower_split', 'lower_long', 'lower_short'),
    x = c(59.7434, 128.9574, 262.8008, 163.8916,
          56.7485, 200.6905, 267.2537, 239.6369),
    y = c(490.0670, 490.0670, 440.2462, 504.2843,
          556.6473, 556.6473, 531.0631, 572.9345),
    radius = c(rep(5.8551, 4), rep(6.0064, 4)))
  edges <- data.frame(
    from = c(1L, 2L, 2L, 5L, 6L, 6L),
    to = c(2L, 3L, 4L, 6L, 7L, 8L),
    width = c(7.1565, 3.2866, 4.9570, 7.3414, 3.4285, 5.1987),
    colour = rep(c('trunk', 'long', 'short'), 2))
  # Colours for trunk, subclonal branches and nodes.
  colours <- c(trunk = '#C82506', long = '#6CB9D2', short = '#479E88', node = '#A6AAA9')
  edge_polygons <- lapply(seq_len(nrow(edges)), function(i) {
    a <- nodes[edges$from[i], ]; b <- nodes[edges$to[i], ]
    dx <- b$x-a$x; dy <- b$y-a$y
    nx <- -dy/sqrt(dx^2+dy^2)*edges$width[i]/2
    ny <- dx/sqrt(dx^2+dy^2)*edges$width[i]/2
    data.frame(x=c(a$x+nx,b$x+nx,b$x-nx,a$x-nx),
               y=-c(a$y+ny,b$y+ny,b$y-ny,a$y-ny),
               group=paste0('edge',i), fill=unname(colours[edges$colour[i]]))
  })
  node_polygons <- lapply(seq_len(nrow(nodes)), function(i) {
    angle <- seq(0,2*pi,length.out=101)
    data.frame(x=nodes$x[i]+nodes$radius[i]*cos(angle),
               y=-nodes$y[i]+nodes$radius[i]*sin(angle),
               group=paste0('node',i),fill=colours[['node']])
  })
  labels <- data.frame(
    x=c(94,185,141,120), y=-c(479,458,541,573),
    text=c('BRAF','chr22q subclone deletion',
           'KRAS/TERT/Fusion/chr22q clonal deletion','BRAF+TERT'),
    angle=c(0,20,0,0), face=c('bold.italic','bold','bold.italic','bold.italic'))
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data=do.call(rbind,edge_polygons),
      ggplot2::aes(x,y,group=group,fill=fill),colour='black',linewidth=0.15) +
    ggplot2::geom_polygon(data=do.call(rbind,node_polygons),
      ggplot2::aes(x,y,group=group,fill=fill),colour='black',linewidth=0.15) +
    ggplot2::geom_text(data=labels,
      ggplot2::aes(x,y,label=text,angle=angle,fontface=face),
      family='Roboto Condensed',size=5.5) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_fixed(xlim=c(25,290),ylim=c(-590,-425),expand=FALSE,clip='off') +
    ggplot2::theme_void(base_family='Roboto Condensed',base_size=13) +
    ggplot2::theme(plot.margin=ggplot2::margin(10,10,10,10))
  paths <- ptc_save(p,output_dir,'Fig4e',width=10,height=6.5)
  list(panel_id='Fig. 4e',status='generated',
       result=paths,plot_type='conceptual schematic',
       reason='Illustrates the two branching architectures without establishing event order.',
       warnings=character())
}
