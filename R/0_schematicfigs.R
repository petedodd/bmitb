## schematic figures
library(here)
library(ggplot2)
library(data.table)
library(officer)
library(rvg)
library(patchwork)
library(MASS)

set.seed(1234)

## relative risk functions in common
source(here("R/riskfunctions.R"))

## ===== example dists: exaggerated
rrtxt <- mean(exp(BL(bmi1, t1, t2))) / mean(exp(BL(bmi0, t1, t2)))
(rrtxt <- round(rrtxt, digits = 2))


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
  annotate(
    geom = "text",
    label = paste0("RR = ", rrtxt),
    x = 30, y = 0.08, col = 2, size = 6
  ) +
  xlab("BMI (kg/m^2)") +
  ylab("Density") +
  theme_classic() +
  ggpubr::grids()

GP

ggsave(GP, file = here("output/eg_dist.png"), w = 6, h = 5)

DP <- dml(ggobj = GP) # convert


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

DPC <- dml(ggobj = GPC) # convert

## ======= lopoffs

## truncated renormalized
lopoff17 <- function(x) {
  w <- pgamma(17, shape = bmirefpop$k, scale = bmirefpop$theta)
  ifelse(x < 17, 0,
    dgamma(x, shape = bmirefpop$k, scale = bmirefpop$theta) /
      (1 - w)
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
    n = 1e3, lty = 3
  ) +
  geom_vline(xintercept = 17, col = 2, lty = 3) +
  xlab("BMI (kg/m^2)") +
  ylab("Density") +
  theme_classic() +
  ggpubr::grids()
GP2

ggsave(GP2, file = here("output/eg_lopoff.png"), w = 6, h = 5)

DP2 <- dml(ggobj = GP2) # convert


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

DP3 <- dml(ggobj = GP3) #convert

## save out relevant plots as PPT
doc <- read_pptx()
doc <- add_slide(doc, layout = "Blank")
doc <- ph_with(doc, DP, location = ph_location_fullsize())
doc <- add_slide(doc, layout = "Blank")
doc <- ph_with(doc, DP2, location = ph_location_fullsize())
doc <- add_slide(doc, layout = "Blank")
doc <- ph_with(doc, DP3, location = ph_location_fullsize())
doc <- add_slide(doc, layout = "Blank")
doc <- ph_with(doc, DPC, location = ph_location_fullsize())
print(doc, target = "~/Downloads/nutrition_schematics.pptx")

## ======= additional counterfactuals

## truncation -> U(17,h)
flat17 <- function(x, h = 25, fac = 1) {
  w <- pgamma(17, shape = bmirefpop$k, scale = bmirefpop$theta)
  ifelse(x < 17, 0,
    fac * dgamma(x, shape = bmirefpop$k, scale = bmirefpop$theta) +
      w * ifelse(x < h, 1 / (h - 17), 0)
  )
}


GP4 <- ggplot() +
  xlim(10, 45) +
  geom_function(
    fun = flat17,
    n = 1e3, col = 2
  ) +
  geom_function(
    fun = function(x) dgamma(x, shape = bmirefpop$k, scale = bmirefpop$theta),
    n = 1e3, lty = 3
  ) +
  geom_vline(xintercept = 17, col = 2, lty = 3) +
  xlab("BMI (kg/m^2)") +
  ylab("Density") +
  theme_classic() +
  ggpubr::grids()
GP4

ggsave(GP4, file = here("output/eg_flat.png"), w = 6, h = 5)


## truncation -> d(x-(h-17))
shift17 <- function(x, h = 25, fac = 1) {
  ifelse(x < 17, 0,
    fac * dgamma(x, shape = bmirefpop$k, scale = bmirefpop$theta) +
      ifelse(x < h,
        dgamma(x - (h - 17), shape = bmirefpop$k, scale = bmirefpop$theta),
        0
      )
  )
}

GP5 <- ggplot() +
  xlim(10, 45) +
  geom_function(
    fun = shift17,
    n = 1e3, col = 2
  ) +
  geom_function(
    fun = function(x) dgamma(x, shape = bmirefpop$k, scale = bmirefpop$theta),
    n = 1e3, lty = 3
  ) +
  geom_vline(xintercept = 17, col = 2, lty = 3) +
  xlab("BMI (kg/m^2)") +
  ylab("Density") +
  theme_classic() +
  ggpubr::grids()
GP5

ggsave(GP5, file = here("output/eg_shift.png"), w = 6, h = 5)

zero17 <- function(x, h = 25) {
  ifelse(x < 17, 0,
    dgamma(x, shape = bmirefpop$k, scale = bmirefpop$theta)
  )
}

extra17 <- function(x) {
  w <- pgamma(17, shape = bmirefpop$k, scale = bmirefpop$theta)
  ifelse(x < 17, 0,
    w * dgamma(x, shape = bmirefpop$k, scale = bmirefpop$theta) /
      (1 - w)
  )
}

GP00 <- ggplot() +
  xlim(10, 45) +
  geom_function(
    fun = zero17,
    n = 1e3, col = 2
  ) +
  geom_function(
    fun = function(x) dgamma(x, shape = bmirefpop$k, scale = bmirefpop$theta),
    n = 1e3, col = 1, lty = 3
  ) +
  geom_vline(xintercept = 17, col = 2, lty = 3) +
  xlab("BMI (kg/m^2)") +
  ylab("Density") +
  theme_classic() +
  ggpubr::grids()
GP00

GP01 <- ggplot() +
  xlim(10, 45) +
  geom_function(
    fun = extra17,
    n = 1e3, col = 2
  ) +
  geom_vline(xintercept = 17, col = 2, lty = 3) +
  xlab("BMI (kg/m^2)") +
  ylab("Density") +
  theme_classic() +
  ggpubr::grids()
GP01

GP02 <- ggplot() +
  xlim(10, 45) +
  geom_function(
    fun = flat17,
    args = list(fac = 0),
    n = 1e3, col = 2
  ) +
  geom_vline(xintercept = 17, col = 2, lty = 3) +
  geom_vline(xintercept = 25, col = 2, lty = 3) +
  xlab("BMI (kg/m^2)") +
  ylab("Density") +
  theme_classic() +
  ggpubr::grids()
GP02

GP03 <- ggplot() +
  xlim(10, 45) +
  geom_function(
    fun = shift17,
    args = list(fac = 0),
    n = 1e3, col = 2
  ) +
  geom_vline(xintercept = 17, col = 2, lty = 3) +
  geom_vline(xintercept = 25, col = 2, lty = 3) +
  xlab("BMI (kg/m^2)") +
  ylab("Density") +
  theme_classic() +
  ggpubr::grids()
GP03


## standardize y axes
ulim <- 0.105
GP00 <- GP00 + expand_limits(y = c(0, ulim))
GP01 <- GP01 + expand_limits(y = c(0, ulim))
GP01 <- GP01 + expand_limits(y = c(0, ulim))
GP02 <- GP02 + expand_limits(y = c(0, ulim))
GP03 <- GP03 + expand_limits(y = c(0, ulim))
GP2 <- GP2 + expand_limits(y = c(0, ulim))
GP4 <- GP4 + expand_limits(y = c(0, ulim))
GP5 <- GP5 + expand_limits(y = c(0, ulim))


## combine
GPall <- ((GP00 | GP01 | GP2) + plot_layout(tag_level = "new")) /
  ((GP00 | GP02 | GP4) + plot_layout(tag_level = "new")) /
  ((GP00 | GP03 | GP5) + plot_layout(tag_level = "new")) +
  plot_annotation(tag_levels = c("A", "1"))

GPall <- ggplotify::as.ggplot(GPall)

GPall <- GPall +
  annotate(
    geom = "text", label = "+", size = unit(14, "pt"),
    x = 1.05 / 3, y = 1 / 6
  ) +
  annotate(
    geom = "text", label = "+", size = unit(14, "pt"),
    x = 1.05 / 3, y = 1 / 6 + 1 / 3
  ) +
  annotate(
    geom = "text", label = "+", size = unit(14, "pt"),
    x = 1.05 / 3, y = 1 / 6 + 2 / 3
  ) +
  annotate(
    geom = "text", label = "=", size = unit(14, "pt"),
    x = 1.0 / 3 + 1 / 3, y = 1 / 6
  ) +
  annotate(
    geom = "text", label = "=", size = unit(14, "pt"),
    x = 1.0 / 3 + 1 / 3, y = 1 / 6 + 1 / 3
  ) +
  annotate(
    geom = "text", label = "=", size = unit(14, "pt"),
    x = 1.0 / 3 + 1 / 3, y = 1 / 6 + 2 / 3
  )
GPall

ggsave(GPall, file = here("output/eg_all.png"), w = 7, h = 7)

