## schematic figures
library(ggplot2)
library(officer)
library(rvg)

## 'data' for these plots
bmirefpop <- data.table(k = 25.0, theta = 1.0) # for testing
## Saunders et al@ risk per one unit increase in BMI was 14.8% (95%CI: 13.3-16.3)
t <- -log(1.148) # risk function parameter
t1 <- -log(1.18) # risk function parameter
exp(-t1) # risk increase with 1 unit decreast
t2 <- -log(1.069) # risk function parameter
exp(-t2) # risk increase with 1 unit decreast
BL <- function(x, t1, t2) {
  ans <- (x - 25)
  less <- ans < 0
  ans[less] <- t1 * ans[less]
  ans[!less] <- t2 * ans[!less]
  ans
}
mean(exp(BL(bmi1, t1, t2))) / mean(exp(BL(bmi0, t1, t2)))



## ===== example dists: exaggerated
K <- 1e5
bmi1 <- rgamma(K, shape = 24, scale = 0.8)
bmi0 <- rgamma(K, shape = bmirefpop$k, scale = bmirefpop$theta)
## rrtxt <- mean(exp(t * (bmi1 - 30))) / mean(exp(t * (bmi0 - 30)))
rrtxt <- mean(exp(BL(bmi1, t1, t2))) / mean(exp(BL(bmi0, t1, t2)))
(rrtxt <- round(rrtxt,digits = 2))

GP <- ggplot() +
  xlim(10, 45) +
  geom_function(
    fun = function(x) dgamma(x, shape = 24, scale = 0.8),
    n = 500, col = 2
  ) +
  geom_function(
    fun = function(x) dgamma(x, shape = bmirefpop$k, scale = bmirefpop$theta),
    n = 500
  ) +
  annotate(geom = "text", label = paste0("RR = ", rrtxt), x = 30, y = 0.08, col = 2, size = 6) +
  xlab("BMI (kg/m^2)") +
  ylab("Density") +
  theme_classic() +
  ggpubr::grids()
GP

ggsave(GP, file = here("output/eg_dist.png"), w = 6, h = 5)


## ==== distribution of where TB comes from
wts <- exp(t * (bmi0 - 30))
bmi0tb <- bmi0[sample(1:K, size = K, replace = TRUE, prob = wts)]
bmiz <- data.table(
  population = c(rep("whole population", K), rep("TB", K)),
  BMI = c(bmi0, bmi0tb)
)

GPC <- ggplot(bmiz, aes(x = BMI, y = after_stat(density), fill = population)) +
  xlim(10, 45) +
  geom_density(alpha = 0.5, adjust = 2) +
  theme_classic() +
  theme(legend.position = "top", legend.title = element_blank()) +
  ggpubr::grids() +
  xlab("BMI (kg/m^2)") +
  ylab("Density")
GPC

ggsave(GPC, file = here("output/eg_tb_vs_pop.png"), w = 6, h = 5)

## ======= lopoffs

## truncated renormalized
lopoff17 <- function(x) {
  ifelse(x < 17, 0,
    dgamma(x, shape = bmirefpop$k, scale = bmirefpop$theta) /
      (1 - pgamma(17, shape = bmirefpop$k, scale = bmirefpop$theta))
  )
}

GP2 <- ggplot() +
  xlim(10, 45) +
  geom_function(
    fun = lopoff17,
    n = 1e3, col = 2
  ) +
  geom_function(
    fun = function(x) dgamma(x, shape = bmirefpop$k, scale = bmirefpop$theta),
    n = 1e3
  ) +
  geom_vline(xintercept = 17, col = 2, lty = 3) +
  xlab("BMI (kg/m^2)") +
  ylab("Density") +
  theme_classic() +
  ggpubr::grids()
GP2

ggsave(GP2, file = here("output/eg_lopoff.png"), w = 6, h = 5)


## ====== risk functions
GP3 <- ggplot() +
  xlim(10, 45) +
  geom_function(
    fun = function(x) BL(x, t1, t2),
    n = 1e3, col = 2
  ) +
  geom_vline(xintercept = 25, col = 2, lty = 2) +
  xlab("BMI (kg/m^2)") +
  ylab("log(Relative Risk)") +
  theme_classic() +
  ggpubr::grids()
GP3

ggsave(GP3, file = here("output/eg_blriskfun.png"), w = 6, h = 5)


## save out relevant plots as PPT
doc <- read_pptx()
doc <- add_slide(doc, layout = "Blank")
doc <- ph_with(doc, GP, location = ph_location_fullsize())
doc <- add_slide(doc, layout = "Blank")
doc <- ph_with(doc, GP2, location = ph_location_fullsize())
doc <- add_slide(doc, layout = "Blank")
doc <- ph_with(doc, GP3, location = ph_location_fullsize())
doc <- add_slide(doc, layout = "Blank")
doc <- ph_with(doc, GPC, location = ph_location_fullsize())
print(doc, target = "~/Downloads/nutrition_schematics.pptx")








