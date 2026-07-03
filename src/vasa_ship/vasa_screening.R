screening <- read.table("data/vasa_ship/screening.csv", sep = ",", header = T)

library(ggplot2)
library(scales)

hist(screening$GRR)
mean(screening$GRR)
sd(screening$GRR)

ggplot(screening, aes(
  x = GRR,
  y = reorder(Specimen.ID, GRR),
  color = GRR > 0.03
)) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "red"),
    guide = "none"
  ) +
  labs(
    x = "GRR",
    y = "Specimen ID"
  ) +
  scale_x_continuous(labels = label_percent(), breaks = seq(0, 0.10, by = 0.02)) +
  labs(
    x = "GRR (%)",
    y = "Specimen ID"
  ) +
  theme_bw()
