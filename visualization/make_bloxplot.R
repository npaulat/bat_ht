#setwd("C:/Users/Nikki/OneDrive - Texas Tech University")
library(readxl)
library(ggplot2)
library(ggthemes)
library(dplyr)

te_prop <- read_excel("bat_te_proportions_summary.xlsx", sheet="Total DNA_RC")

boxplot(Total_Proportion~Group, data=te_prop, 
        ylab = "Genome Proportion", ylim=c(0,0.25), boxwex=0.2, col="grey")
#, subset = Group=="Bats"
#xlab=te_prop$Group,
max(te_prop$Total_Proportion)

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}
#te_plot <- ggplot(te_prop, aes(x=Group, y=Total_Proportion, fill=Group), position_dodge(1)) +
#  geom_boxplot(width=0.2) + scale_y_continuous(limits=(c(NA, 0.25))) +
#  labs(title= "Total DNA TE Content", x="", y="Genome Proportion")  
te_plot <- ggplot(te_prop, aes(x=Group, y=Total_Proportion), position_dodge(1)) +
    geom_boxplot(width=0.2) +
    coord_cartesian(ylim = c(0.01, 0.24)) +  
    labs(title= "Total DNA TE Content", x="", y="Genome Proportion") 
te_plot + geom_rangeframe() + scale_color_grey() + 
  theme(legend.position="none", panel.background = element_blank(), 
        plot.title = element_text(size = 8, hjust=0.5, color="black"),
        axis.text = element_text(size=8, colour="black"),
        axis.title.y = element_text(size=8, color="black"),
        panel.border = element_rect(color="black", fill=NA, size=0.5),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"))
#  + scale_fill_brewer(palette="Dark2")
ggsave("te_proportions_inset.png", dpi = 900, width = 2.25, height = 2.75, units ="in")
                           
