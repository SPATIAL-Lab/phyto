pco2.out.925 <- read.csv("model_out/925out_cultonly.csv")
pco2.out.MD012392 <- read.csv("model_out/MD01out_cultonly.csv")
#pco2.out.688B <- read.csv("model_out/pco2_out688B.csv")

ggplot() + 
  geom_abline(slope=1, intercept=0, linetype=2,color="gray") +
  geom_errorbar(data = pco2.out.925, mapping = aes(x=ice, xmin =ice-6, xmax =ice+6, y=mean, ymin=low, ymax=high), color="black", linewidth=0.2, width=0) +
  geom_errorbar(data = pco2.out.925, mapping = aes(x=ice, xmin =ice-6, xmax =ice+6, y=mean), color="black", linewidth=0.2, width=0) +
  geom_errorbar(data = pco2.out.MD012392, mapping = aes(x=ice, xmin =ice-6, xmax =ice+6, y=mean, ymin=low, ymax=high), color="black", linewidth=0.2, width=0) +
  geom_errorbar(data = pco2.out.MD012392, mapping = aes(x=ice, xmin =ice-6, xmax =ice+6, y=mean), color="black", linewidth=0.2, width=0) +
  geom_point(data = pco2.out.MD012392, mapping = aes(x=ice, y=mean), size=2, color = "cornflowerblue") +
  geom_point(data = pco2.out.925, mapping = aes(x=ice, y=mean), size=2, color = "purple") +
  #geom_errorbar(data = pco2.out.688B, mapping = aes(x=ice, xmin =ice-6, xmax =ice+6, y=mean, ymin=X2.5., ymax=X97.5.), color="cyan", width=0) +
  #geom_errorbar(data = pco2.out.688B, mapping = aes(x=ice, xmin =ice-6, xmax =ice+6, y=mean), color="cyan", width=0) +
  #geom_point(data = pco2.out.688B, mapping = aes(x=ice, y=mean), color = "blue") +
  ylim(150,330) +
  xlim(150,330) +
  labs(x= expression(ice~core~CO[2]~(ppmv)), y = expression(alkenone~CO[2]~(ppmv))) +   
  theme_bw() 
