cols = readRDS("data-raw/clust2out_151673.RDS")
test1 = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(cols$col, levels = c(5,4,1,6,2,3,7))), alpha = cols$alpha,
            size = 2, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 1.5, family = "Lucida Sans Unicode", show.legend = F,alpha = 0.5) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F, col = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()
truth_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(sce$layer_guess_reordered)),
            size = 2, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 1.5, family = "Lucida Sans Unicode", show.legend = F, alpha = 0.5) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()

maynardclust_out = ggplot(data.frame(positions), aes(x, -y)) +
  geom_text(aes(color = factor(sce$HVG_PCA_spatial, levels = c(5,1,3,2,7,8,6,4))),
            size = 2, label = "\u2B22", family = "Lucida Sans Unicode") +
  geom_text(color = "black", label = "\u2B21", size = 1.5, family = "Lucida Sans Unicode", show.legend = F, alpha = 0.5) +
  labs(color = "Cluster", alpha = "Proportion", x = NULL, y = NULL) +
  guides(alpha = F, col = F) + 
  scale_color_viridis_d(option = "A") +
  scale_alpha_continuous(limits = c(0,1), breaks = seq(0.1,1,0.1), range = c(0,1))+
  theme_void() + coord_fixed()
(truth_out | maynardclust_out| test1)/test3 + plot_annotation(tag_levels = "A")
