
mycolor <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
CCB_plot_style <- function(...)
{
  
  theme_bw()+
    
    
    theme(
      text = element_text(size=16,face="bold",color="black",family="sans"),
      axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
      axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
      axis.text.x=element_text(size=18,family="sans",face="bold",colour="#666666"),
      axis.text.y=element_text(size=18,family="sans",face="bold",colour="#666666"),
      plot.margin = margin(10, 10, 10, 10),
      axis.title.y = element_text(color = "black", size = 16, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
      legend.box.background = element_rect(colour = "black",linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank() ,
      panel.background = element_blank(),
      axis.line = element_line(color = '#666666',linewidth = 0.8),
      axis.ticks =element_line(color="#666666",linewidth =1),
      axis.ticks.length=unit(c(-0.3,0.3), "cm"))
  
  
}