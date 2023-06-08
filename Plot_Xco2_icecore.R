pco2.out.925 <- read.csv("model_out/pco2_out925.csv")
pco2.out.MD012392 <- read.csv("model_out/pco2_outMD012392.csv")
pco2.out.688B <- read.csv("model_out/pco2_out688B.csv")

ggplot() + 
  geom_abline(slope=1, intercept=0, color="black") +
  geom_errorbar(data = pco2.out.925, mapping = aes(x=ice, xmin =ice-6, xmax =ice+6, y=mean, ymin=X2.5., ymax=X97.5.), color="pink", width=0) +
  geom_errorbar(data = pco2.out.925, mapping = aes(x=ice, xmin =ice-6, xmax =ice+6, y=mean), color="pink", width=0) +
  geom_point(data = pco2.out.925, mapping = aes(x=ice, y=mean), color = "red") +
  geom_errorbar(data = pco2.out.MD012392, mapping = aes(x=ice, xmin =ice-6, xmax =ice+6, y=mean, ymin=X2.5., ymax=X97.5.), color="gray", width=0) +
  geom_errorbar(data = pco2.out.MD012392, mapping = aes(x=ice, xmin =ice-6, xmax =ice+6, y=mean), color="gray", width=0) +
  geom_point(data = pco2.out.MD012392, mapping = aes(x=ice, y=mean), color = "purple") +
  geom_errorbar(data = pco2.out.688B, mapping = aes(x=ice, xmin =ice-6, xmax =ice+6, y=mean, ymin=X2.5., ymax=X97.5.), color="cyan", width=0) +
  geom_errorbar(data = pco2.out.688B, mapping = aes(x=ice, xmin =ice-6, xmax =ice+6, y=mean), color="cyan", width=0) +
  geom_point(data = pco2.out.688B, mapping = aes(x=ice, y=mean), color = "blue") +
  ylim(150,350) +
  xlim(150,350) +
  labs(x= "ice core pCO2 (uatm)", y = "pCO2 (uatm)") +   
  theme_bw()
