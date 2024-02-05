# loading libraries
library(mgcv)
library(here)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)

# loading in the data
here()
data_comb <- read.csv("GAM Data/updated cb ratios.csv")

# renaming species codes
data_comb$svspp[data_comb$svspp == 15] <- 'Spiny.Dogfish'
data_comb$svspp[data_comb$svspp == 23] <- 'Winter.Skate'
data_comb$svspp[data_comb$svspp == 75] <- 'Pollock'
data_comb$svspp[data_comb$svspp == 76] <- 'White.Hake'
data_comb$svspp[data_comb$svspp == 77] <- 'Red.Hake'
data_comb$svspp[data_comb$svspp == 197] <- 'Goosefish'
data_comb$svspp[data_comb$svspp == 28] <- 'Thorny.Skate'
data_comb$svspp[data_comb$svspp == 73] <- 'Atlantic.Cod'
data_comb$svspp[data_comb$svspp == 164] <- 'Sea.Raven'
data_comb$svspp[data_comb$svspp == 72] <- 'Silver.Hake'
data_comb$svspp[data_comb$svspp == 103] <- 'Summer.Flounder'
data_comb$svspp[data_comb$svspp == 104] <- 'Fourspot.Flounder'
data_comb$svspp[data_comb$svspp == 108] <- 'Windowpane'
data_comb$svspp[data_comb$svspp == 112] <- 'Buckler.Dory'
data_comb$svspp[data_comb$svspp == 135] <- 'Blue.Fish'
data_comb$svspp[data_comb$svspp == 163] <- 'Longhorn.Sculpin'
data_comb$svspp[data_comb$svspp == 172] <- 'Striped.Searobin'

# filtering to all prey types 
data_comb_f <- filter(data_comb, prey == 'all')

# omitting all NA's
data_comb_f <- data_comb_f %>% na.omit(cb)

# omitting all 0's from the dataset
rows_to_keep <- data_comb_f$cb > 0.000000000
data_comb_f <- data_comb_f[rows_to_keep, ]

# starting time series in 1981 for sensitivity analysis
data_comb_f <- data_comb_f %>% 
  filter(year > 1980)

# log transforming c:b ratios 
data_comb_f$cb <- log(data_comb_f$cb)

data_comb_f$cb <- data_comb_f$cb + 3

# looking at c/b ratios over time in spring 
dc_spring <- filter(data_comb_f, seacat == 'SPRING')

## Spiny Dogfish
### GLM
# Fit a GAM model
dc_spring_s.d <- filter(dc_spring, svspp == 'Spiny.Dogfish')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_s.d$cb)/sqrt(length(dc_spring_s.d$cb)))

dc_spring_s.d$lower <- dc_spring_s.d$cb - se # calculating lower bound
dc_spring_s.d$upper <- dc_spring_s.d$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_spring_s.d) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_s.d$year, y = dc_spring_s.d$cb, dc_spring_s.d$bt, dc_spring_s.d$lower, dc_spring_s.d$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "bt","lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod1<- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Spiny Dogfish")

ggsave("spiny dogfish GAM.png", plot = mod1, width = 6, height = 5, units = "in")

# plot data
s.d_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Spiny Dogfish")
s.d_plot

### GAM
dc_spring_s.d <- filter(dc_spring, svspp == 'Spiny.Dogfish')
dc_spring_s.d$z_bt <- (dc_spring_s.d$bt - mean(dc_spring_s.d$bt))/(dc_spring_s.d$bt_sd)
target_sd <- 0.5
mod1 <- gam(cb ~ s(year) + s(bt) , data=dc_spring_s.d, family = 'gaussian')
summary(mod1)
mod1$aic
plot(mod1, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod1)

## Winter Skate
dc_spring_w.s <- filter(dc_spring, svspp == 'Winter.Skate')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_w.s$cb)/sqrt(length(dc_spring_w.s$cb)))

dc_spring_w.s$lower <- dc_spring_w.s$cb - se # calculating lower bound
dc_spring_w.s$upper <- dc_spring_w.s$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_spring_w.s) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_w.s$year, y = dc_spring_w.s$cb, dc_spring_w.s$lower, dc_spring_w.s$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod2<- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Winter Skate")

ggsave("winter skate GAM.png", plot = mod2, width = 6, height = 5, units = "in")

# plot data
w.s_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Winter Skate")
w.s_plot

### GAM
dc_spring_w.s <- filter(dc_spring, svspp == 'Winter.Skate')
dc_spring_w.s$z_bt <- (dc_spring_w.s$bt - mean(dc_spring_w.s$bt))/(dc_spring_w.s$bt_sd)
mod2 <- gam(cb ~ s(year) + s(bt), data=dc_spring_w.s, family = 'gaussian')
summary(mod2)
mod2$aic
plot(mod2, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod2)

## Pollock
dc_spring_p <- filter(dc_spring, svspp == 'Pollock')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_p$cb)/sqrt(length(dc_spring_p$cb)))

dc_spring_p$lower <- dc_spring_p$cb - se # calculating lower bound
dc_spring_p$upper <- dc_spring_p$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_spring_p) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_p$year, y = dc_spring_p$cb, dc_spring_p$lower, dc_spring_p$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod3<- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Pollock")

ggsave("pollock GAM.png", plot = mod3, width = 6, height = 5, units = "in")

# plot data
p_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Pollock")
p_plot

### GAM
dc_spring_p <- filter(dc_spring, svspp == 'Pollock')
dc_spring_p$z_bt <- (dc_spring_p$bt - mean(dc_spring_p$bt))/(dc_spring_p$bt_sd)
mod3 <- gam(cb ~ s(year) + s(bt), data=dc_spring_p, family = 'gaussian')
summary(mod3)
mod3$aic
plot(mod3, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod3)

## White Hake
dc_spring_w.h <- filter(dc_spring, svspp == 'White.Hake')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_w.h$cb)/sqrt(length(dc_spring_w.h$cb)))

dc_spring_w.h$lower <- dc_spring_w.h$cb - se # calculating lower bound
dc_spring_w.h$upper <- dc_spring_w.h$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_spring_w.h) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_w.h$year, y = dc_spring_w.h$cb, dc_spring_w.h$lower, dc_spring_w.h$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod4<- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "White Hake")

ggsave("white hake GAM.png", plot = mod4, width = 6, height = 5, units = "in")

# plot data
w.h_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("White Hake")
w.h_plot

### GAM
dc_spring_w.h <- filter(dc_spring, svspp == 'White.Hake')
dc_spring_w.h$z_bt <- (dc_spring_w.h$bt - mean(dc_spring_w.h$bt))/(dc_spring_w.h$bt_sd)
mod4 <- gam(cb ~ s(year) + s(bt), data=dc_spring_w.h, family = 'gaussian')
summary(mod4)
mod4$aic
plot(mod4, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod4)

## Red Hake
dc_spring_r.h <- filter(dc_spring, svspp == 'Red.Hake')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_r.h$cb)/sqrt(length(dc_spring_r.h$cb)))

dc_spring_r.h$lower <- dc_spring_r.h$cb - se # calculating lower bound
dc_spring_r.h$upper <- dc_spring_r.h$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_spring_r.h) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_r.h$year, y = dc_spring_r.h$cb, dc_spring_r.h$lower, dc_spring_r.h$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod5<- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Red Hake")

ggsave("red hake GAM.png", plot = mod5, width = 6, height = 5, units = "in")

# plot data
r.h_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Red Hake")
r.h_plot

### GAM
dc_spring_r.h <- filter(dc_spring, svspp == 'Red.Hake')
dc_spring_r.h$z_bt <- (dc_spring_r.h$bt - mean(dc_spring_r.h$bt))/(dc_spring_r.h$bt_sd)
mod5 <- gam(cb ~ s(year) + s(bt), data=dc_spring_r.h, family = 'gaussian')
summary(mod5)
mod5$aic
plot(mod5, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod5)

## Goosefish
dc_spring_g <- filter(dc_spring, svspp == 'Goosefish')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_g$cb)/sqrt(length(dc_spring_g$cb)))

dc_spring_g$lower <- dc_spring_g$cb - se # calculating lower bound
dc_spring_g$upper <- dc_spring_g$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_spring_g) # make predictions with confidence intervals 
model$aic
summary(model)
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_g$year, y = dc_spring_g$cb, dc_spring_g$lower, dc_spring_g$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod6 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Goosefish")

ggsave("goosefish GAM.png", plot = mod6, width = 6, height = 5, units = "in")

# plot data
g_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Goosefish")
g_plot

### GAM
dc_spring_g <- filter(dc_spring, svspp == 'Goosefish')
dc_spring_g$z_bt <- (dc_spring_g$bt - mean(dc_spring_g$bt))/(dc_spring_g$bt_sd)
mod6 <- gam(cb ~ s(year) + s(bt), data=dc_spring_g, family = 'gaussian')
summary(mod6)
mod6$aic
plot(mod6, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod6)

## Thorny Skate
dc_spring_t.s <- filter(dc_spring, svspp == 'Thorny.Skate')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_t.s$cb)/sqrt(length(dc_spring_t.s$cb)))

dc_spring_t.s$lower <- dc_spring_t.s$cb - se # calculating lower bound
dc_spring_t.s$upper <- dc_spring_t.s$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_spring_t.s) # make predictions with confidence intervals 
model$aic
summary(model)
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_t.s$year, y = dc_spring_t.s$cb, dc_spring_t.s$lower, dc_spring_t.s$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod7 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Thorny Skate")

ggsave("thorny skate GAM.png", plot = mod7, width = 6, height = 5, units = "in")

# plot data
t.s_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Thorny Skate")
t.s_plot

### GAM
dc_spring_t.s <- filter(dc_spring, svspp == 'Thorny.Skate')
dc_spring_t.s$z_bt <- (dc_spring_t.s$bt - mean(dc_spring_t.s$bt))/(dc_spring_t.s$bt_sd)
mod7 <- gam(cb ~ s(year) + s(bt), data=dc_spring_t.s, family = 'gaussian')
summary(mod7)
mod7$aic
plot(mod7, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod7)

## Atlantic Cod
dc_spring_a.c <- filter(dc_spring, svspp == 'Atlantic.Cod')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_a.c$cb)/sqrt(length(dc_spring_a.c$cb)))

dc_spring_a.c$lower <- dc_spring_a.c$cb - se # calculating lower bound
dc_spring_a.c$upper <- dc_spring_a.c$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_spring_a.c) # make predictions with confidence intervals 
model$aic
summary(model)
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_a.c$year, y = dc_spring_a.c$cb, dc_spring_a.c$lower, dc_spring_a.c$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod8 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Atlantic Cod")

ggsave("atlantic cod GAM.png", plot = mod8, width = 6, height = 5, units = "in")

# plot data
a.c_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Atlantic Cod")
a.c_plot

### GAM
dc_spring_a.c <- filter(dc_spring, svspp == 'Atlantic.Cod')
dc_spring_a.c$z_bt <- (dc_spring_a.c$bt - mean(dc_spring_a.c$bt))/(dc_spring_a.c$bt_sd)
mod8 <- gam(cb ~ s(year) + s(bt), data=dc_spring_a.c, family = 'gaussian')
summary(mod8)
mod8$aic
plot(mod8, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod8)

## Sea Raven
dc_spring_s.r <- filter(dc_spring, svspp == 'Sea.Raven')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_s.r$cb)/sqrt(length(dc_spring_s.r$cb)))

dc_spring_s.r$lower <- dc_spring_s.r$cb - se # calculating lower bound
dc_spring_s.r$upper <- dc_spring_s.r$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_spring_s.r) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_s.r$year, y = dc_spring_s.r$cb, dc_spring_s.r$lower, dc_spring_s.r$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod9 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Sea Raven")

ggsave("sea raven GAM.png", plot = mod9, width = 6, height = 5, units = "in")

# plot data
s.r_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Sea Raven")
s.r_plot

### GAM
dc_spring_s.r <- filter(dc_spring, svspp == 'Sea.Raven')
dc_spring_s.r$z_bt <- (dc_spring_s.r$bt - mean(dc_spring_s.r$bt))/(dc_spring_s.r$bt_sd)
mod9 <- gam(cb ~ s(year) + s(bt), data=dc_spring_s.r, family = 'gaussian')
summary(mod9)
mod9$aic
plot(mod9, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod9)

## Silver Hake
dc_spring_s.h <- filter(dc_spring, svspp == 'Silver.Hake')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_s.h$cb)/sqrt(length(dc_spring_s.h$cb)))

dc_spring_s.h$lower <- dc_spring_s.h$cb - se # calculating lower bound
dc_spring_s.h$upper <- dc_spring_s.h$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_spring_s.h) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_s.h$year, y = dc_spring_s.h$cb, dc_spring_s.h$lower, dc_spring_s.h$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod10 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Silver Hake")

ggsave("silver hake GAM.png", plot = mod10, width = 6, height = 5, units = "in")

# plot data
s.h_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Silver Hake")
s.h_plot

### GAM
dc_spring_s.h <- filter(dc_spring, svspp == 'Silver.Hake')
dc_spring_s.h$z_bt <- (dc_spring_s.h$bt - mean(dc_spring_s.h$bt))/(dc_spring_s.h$bt_sd)
mod10 <- gam(cb ~ s(year) + s(bt), data=dc_spring_s.h, family = 'gaussian')
summary(mod10)
mod10$aic
plot(mod10, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod10)

## Summer Flounder
dc_spring_s.f <- filter(dc_spring, svspp == 'Summer.Flounder')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_s.f$cb)/sqrt(length(dc_spring_s.f$cb)))

dc_spring_s.f$lower <- dc_spring_s.f$cb - se # calculating lower bound
dc_spring_s.f$upper <- dc_spring_s.f$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_spring_s.f) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_s.f$year, y = dc_spring_s.f$cb, dc_spring_s.f$lower, dc_spring_s.f$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod11<- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Summer Flounder")

ggsave("summr flounder GAM.png", plot = mod11, width = 6, height = 5, units = "in")

# plot data
s.f_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Summer Flounder")
s.f_plot

### GAM
dc_spring_s.f <- filter(dc_spring, svspp == 'Summer.Flounder')
dc_spring_s.f$z_bt <- (dc_spring_s.f$bt - mean(dc_spring_s.f$bt))/(dc_spring_s.f$bt_sd)
mod11 <- gam(cb ~ s(year) + s(bt), data=dc_spring_s.f, family = 'gaussian')
summary(mod11)
mod11$aic
plot(mod11, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod11)

## Fourspot Flounder
dc_spring_f.f <- filter(dc_spring, svspp == 'Fourspot.Flounder')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_f.f$cb)/sqrt(length(dc_spring_f.f$cb)))

dc_spring_f.f$lower <- dc_spring_f.f$cb - se # calculating lower bound
dc_spring_f.f$upper <- dc_spring_f.f$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_spring_f.f) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_f.f$year, y = dc_spring_f.f$cb, dc_spring_f.f$lower, dc_spring_f.f$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod12 <-ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Fourspot Flounder")

ggsave("fourspot flounder GAM.png", plot = mod12, width = 6, height = 5, units = "in")

# plot data
f.f_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Fourspot Flounder")
f.f_plot

### GAM
dc_spring_f.f <- filter(dc_spring, svspp == 'Fourspot.Flounder')
dc_spring_f.f$z_bt <- (dc_spring_f.f$bt - mean(dc_spring_f.f$bt))/(dc_spring_f.f$bt_sd)
mod12 <- gam(cb ~ s(year) + s(bt), data=dc_spring_f.f, family = 'gaussian')
summary(mod12)
mod12$aic
plot(mod12, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod12)

## Windowpane
dc_spring_w <- filter(dc_spring, svspp == 'Windowpane')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_w$cb)/sqrt(length(dc_spring_w$cb)))

dc_spring_w$lower <- dc_spring_w$cb - se # calculating lower bound
dc_spring_w$upper <- dc_spring_w$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_spring_w) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_w$year, y = dc_spring_w$cb, dc_spring_w$lower, dc_spring_w$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod13 <-ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Windowpane")

ggsave("windowpane GAM.png", plot = mod13, width = 6, height = 5, units = "in")

# plot data
w_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Windowpane")
w_plot

### GAM
dc_spring_w <- filter(dc_spring, svspp == 'Windowpane')
dc_spring_w$z_bt <- (dc_spring_w$bt - mean(dc_spring_w$bt))/(dc_spring_w$bt_sd)
mod13 <- gam(cb ~ s(year) + s(bt), data=dc_spring_w, family = 'gaussian')
summary(mod13)
mod13$aic
plot(mod13, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod13)

## Buckler Dory
dc_spring_b.d <- filter(dc_spring, svspp == 'Buckler.Dory')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_b.d$cb)/sqrt(length(dc_spring_b.d$cb)))

dc_spring_b.d$lower <- dc_spring_b.d$cb - se # calculating lower bound
dc_spring_b.d$upper <- dc_spring_b.d$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year, k = 1), data = dc_spring_b.d) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_b.d$year, y = dc_spring_b.d$cb, dc_spring_b.d$lower, dc_spring_b.d$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod14 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Buckler Dory")

ggsave("buckler dory GAM.png", plot = mod14, width = 6, height = 5, units = "in")

# plot data
b.d_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Buckler Dory")
b.d_plot

### GAM
dc_spring_b.d <- filter(dc_spring, svspp == 'Buckler.Dory')
dc_spring_b.d$z_bt <- (dc_spring_b.d$bt - mean(dc_spring_b.d$bt))/(dc_spring_b.d$bt_sd)
mod14 <- gam(cb ~ s(year) , data=dc_spring_b.d, family = 'gaussian')
mod14.a <- gam(cb ~ s(bt) , data=dc_spring_b.d, family = 'gaussian')
summary(mod14)
summary(mod14.a)
mod14$aic
mod14.a$aic
plot(mod14, pages=1, scheme=1, unconditional=TRUE)
plot(mod14.a, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod14)
gam.check(mod14.a)

## Blue Fish
dc_spring_b.f <- filter(dc_spring, svspp == 'Blue.Fish')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_b.f$cb)/sqrt(length(dc_spring_b.f$cb)))

dc_spring_b.f$lower <- dc_spring_b.f$cb - se # calculating lower bound
dc_spring_b.f$upper <- dc_spring_b.f$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year, k = 2), data = dc_spring_b.f) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_b.f$year, y = dc_spring_b.f$cb, dc_spring_b.f$lower, dc_spring_b.f$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod15 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Blue Fish")

# plot data
b.f_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1.5) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Blue Fish")
b.f_plot

### GAM
dc_spring_b.f <- filter(dc_spring, svspp == 'Blue.Fish')
dc_spring_b.f$z_bt <- (dc_spring_b.f$bt - mean(dc_spring_b.f$bt))/(dc_spring_b.f$bt_sd)
mod15 <- gam(cb ~ s(year) + s(bt), data=dc_spring_b.f, family = 'gaussian')
summary(mod15)
mod15$aic
plot(mod15, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod15)

## Longhorn Sculpin
dc_spring_l.s <- filter(dc_spring, svspp == 'Longhorn.Sculpin')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_l.s$cb)/sqrt(length(dc_spring_l.s$cb)))

dc_spring_l.s$lower <- dc_spring_l.s$cb - se # calculating lower bound
dc_spring_l.s$upper <- dc_spring_l.s$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_spring_l.s) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_l.s$year, y = dc_spring_l.s$cb, dc_spring_l.s$lower, dc_spring_l.s$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod16 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Longhorn Sculpin")

ggsave("longhorn sculpin GAM.png", plot = mod16, width = 6, height = 5, units = "in")

# plot data
l.s_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Longhorn Sculpin")
l.s_plot

### GAM
dc_spring_l.s <- filter(dc_spring, svspp == 'Longhorn.Sculpin')
dc_spring_l.s$z_bt <- (dc_spring_l.s$bt - mean(dc_spring_l.s$bt))/(dc_spring_l.s$bt_sd)
mod16 <- gam(cb ~ s(year) + s(z_bt), data=dc_spring_l.s, family = 'gaussian')
summary(mod16)
mod16$aic
plot(mod16, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod16)

## Striped Searobin
dc_spring_s.s <- filter(dc_spring, svspp == 'Striped.Searobin')

# calculating individual confidence intervals for points
se <- (sd(dc_spring_s.s$cb)/sqrt(length(dc_spring_s.s$cb)))

dc_spring_s.s$lower <- dc_spring_s.s$cb - se # calculating lower bound
dc_spring_s.s$upper <- dc_spring_s.s$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year, k = 1), data = dc_spring_s.s) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_spring_s.s$year, y = dc_spring_s.s$cb, dc_spring_s.s$lower, dc_spring_s.s$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod17 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  scale_x_continuous(breaks = seq(2005, 2017, by = 2)) +
  labs(title = "Striped Searobin")

ggsave("striped searobin GAM.png", plot = mod17, width = 6, height = 5, units = "in")

# plot data
s.s_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  scale_x_continuous(breaks = seq(2005, 2017, by = 2)) +
  ggtitle("Striped Searobin")
s.s_plot

### GAM
dc_spring_s.s <- filter(dc_spring, svspp == 'Striped.Searobin')
dc_spring_s.s$z_bt <- (dc_spring_s.s$bt - mean(dc_spring_s.s$bt))/(dc_spring_s.s$bt_sd)
mod17 <- gam(cb ~ s(year), data=dc_spring_s.s, family = 'gaussian')
mod17.a <- gam(cb ~ s(z_bt), data=dc_spring_s.s, family = 'gaussian')
summary(mod17)
summary(mod17.a)
mod17$aic
mod17.a$aic
plot(mod17, pages=1, scheme=1, unconditional=TRUE)
plot(mod17.a, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod17)
gam.check(mod17.a)

# combining all of the linear plots into one figure
glm_cb <- ggarrange(s.d_plot, w.s_plot, p_plot, w.h_plot, r.h_plot,
                    g_plot, t.s_plot, a.c_plot, s.r_plot, s.h_plot,
                    s.f_plot, f.f_plot, w_plot, b.d_plot, b.f_plot,
                    l.s_plot, s.s_plot, ncol = 5, nrow = 4)
glm_cb

ggsave("glms.png", plot = glm_cb, width = 15, height = 7, units = "in")

# visually comparing GAMS across species
## GAM of year v. cb --> I save each of these figures (total of 2 figures) once the matrix
## of 3x3 is filled in order to have two large figures that I can compare all species with
par(mfrow=c(3,3))
mod1 <- gam(cb ~ s(year), data=dc_spring_s.d, family = 'gaussian')
plot(mod1, scheme=1, unconditional=TRUE, main = "Spiny Dogfish") 
mod2 <- gam(cb ~ s(year), data=dc_spring_w.s, family = 'gaussian')
plot(mod2, scheme=1, unconditional=TRUE, main = "Winter Skate") 
mod3 <- gam(cb ~ s(year), data=dc_spring_p, family = 'gaussian')
plot(mod3, scheme=1, unconditional=TRUE, main = "Pollock") 
mod4 <- gam(cb ~ s(year), data=dc_spring_w.h, family = 'gaussian')
plot(mod4, scheme=1, unconditional=TRUE, main = "White Hake")
mod5 <- gam(cb ~ s(year), data=dc_spring_r.h, family = 'gaussian')
plot(mod5, scheme=1, unconditional=TRUE, main = "Red Hake")
mod6 <- gam(cb ~ s(year), data=dc_spring_g, family = 'gaussian')
plot(mod6, scheme=1, unconditional=TRUE, main = "Goosefish")
mod7 <- gam(cb ~ s(year), data=dc_spring_t.s, family = 'gaussian')
plot(mod7, scheme=1, unconditional=TRUE, main = "Thorny Skate")
mod8 <- gam(cb ~ s(year), data=dc_spring_a.c, family = 'gaussian')
plot(mod8, scheme=1, unconditional=TRUE, main = "Atlantic Cod")
mod9 <- gam(cb ~ s(year), data=dc_spring_s.r, family = 'gaussian')
plot(mod9, scheme=1, unconditional=TRUE, main = "Sea Raven")
mod10 <- gam(cb ~ s(year), data=dc_spring_s.h, family = 'gaussian')
plot(mod10, scheme=1, unconditional=TRUE, main = "Silver Hake")
mod11 <- gam(cb ~ s(year), data=dc_spring_s.f, family = 'gaussian')
plot(mod11, scheme=1, unconditional=TRUE, main = "Summer Flounder")
mod12 <- gam(cb ~ s(year), data=dc_spring_f.f, family = 'gaussian')
plot(mod12, scheme=1, unconditional=TRUE, main = "Fourspot Flounder")
mod13 <- gam(cb ~ s(year), data=dc_spring_w, family = 'gaussian')
plot(mod13, scheme=1, unconditional=TRUE, main = "Windowpane")
mod14 <- gam(cb ~ s(year), data=dc_spring_b.d, family = 'gaussian')
plot(mod14, scheme=1, unconditional=TRUE, main = "Buckler Dory")
mod15 <- gam(cb ~ s(year), data=dc_spring_b.f, family = 'gaussian')
plot(mod15, scheme=1, unconditional=TRUE, main = "Blue Fish")
mod16 <- gam(cb ~ s(year), data=dc_spring_l.s, family = 'gaussian')
plot(mod16, scheme=1, unconditional=TRUE, main = "Longhorn Sculpin")
mod17 <- gam(cb ~ s(year), data=dc_spring_s.s, family = 'gaussian')
plot(mod17, scheme = 1, unconditional=TRUE, main = "Striped Searobin")

## GAM of bt v. cb --> same note as above
par(mfrow = c(3,3))
mod1 <- gam(cb ~ s(bt), data=dc_spring_s.d, family = 'gaussian')
plot(mod1, scheme=1, unconditional=TRUE, main = "Spiny Dogfish")
mod2 <- gam(cb ~ s(bt), data=dc_spring_w.s, family = 'gaussian')
plot(mod2, scheme=1, unconditional=TRUE, main = "Winter Skate")
mod3 <- gam(cb ~ s(bt), data=dc_spring_p, family = 'gaussian')
plot(mod3, scheme=1, unconditional=TRUE, main = "Pollock")
mod4 <- gam(cb ~ s(bt), data=dc_spring_w.h, family = 'gaussian')
plot(mod4, scheme=1, unconditional=TRUE, main = "White Hake")
mod5 <- gam(cb ~ s(bt), data=dc_spring_r.h, family = 'gaussian')
plot(mod5, scheme=1, unconditional=TRUE, main = " Red Hake")
mod6 <- gam(cb ~ s(bt), data=dc_spring_g, family = 'gaussian')
plot(mod6, scheme=1, unconditional=TRUE, main = "Goosefish")
mod7 <- gam(cb ~ s(bt), data=dc_spring_t.s, family = 'gaussian')
plot(mod7, scheme=1, unconditional=TRUE, main = "Thorny Skate")
mod8 <- gam(cb ~ s(bt), data=dc_spring_a.c, family = 'gaussian')
plot(mod8, scheme=1, unconditional=TRUE, main = "Atlantic Cod")
mod9 <- gam(cb ~ s(bt), data=dc_spring_s.r, family = 'gaussian')
plot(mod9, scheme=1, unconditional=TRUE, main = "Sea Raven")
mod10 <- gam(cb ~ s(bt), data=dc_spring_s.h, family = 'gaussian')
plot(mod10, scheme=1, unconditional=TRUE, main = "Silver Hake")
mod11 <- gam(cb ~ s(bt), data=dc_spring_s.f, family = 'gaussian')
plot(mod11, scheme=1, unconditional=TRUE, main = "Summer Flounder")
mod12 <- gam(cb ~ s(bt), data=dc_spring_f.f, family = 'gaussian')
plot(mod12, scheme=1, unconditional=TRUE, main = "Fourspot Flounder")
mod13 <- gam(cb ~ s(bt), data=dc_spring_w, family = 'gaussian')
plot(mod13, scheme=1, unconditional=TRUE, main = "Windowpane")
mod14 <- gam(cb ~ s(bt) , data=dc_spring_b.d, family = 'gaussian')
plot(mod14, scheme=1, unconditional=TRUE, main = "Buckler Dory")
mod15 <- gam(cb ~ s(bt), data=dc_spring_b.f, family = 'gaussian')
plot(mod15, scheme=1, unconditional=TRUE, main = "Blue Fish")
mod16 <- gam(cb ~ s(z_bt), data=dc_spring_l.s, family = 'gaussian')
plot(mod16, scheme=1, unconditional=TRUE, main = "Longhorn Sculpin")
mod17 <- gam(cb ~ s(z_bt), data=dc_spring_s.s, family = 'gaussian')
plot(mod17, scheme = 1, unconditional=TRUE, main = "Striped Searobin")

# looking at c/b ratios over time in fall 
dc_fall <- filter(data_comb_f, seacat == 'FALL')

## Spiny Dogfish
dc_fall_s.d <- filter(dc_fall, svspp == 'Spiny.Dogfish')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_s.d$cb)/sqrt(length(dc_fall_s.d$cb)))

dc_fall_s.d$lower <- dc_fall_s.d$cb - se # calculating lower bound
dc_fall_s.d$upper <- dc_fall_s.d$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_s.d) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_s.d$year, y = dc_fall_s.d$cb, dc_fall_s.d$lower, dc_fall_s.d$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod1 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Spiny Dogfish")

ggsave("spiny dogfish fall GAM.png", plot = mod1, width = 6, height = 5, units = "in")

# plot data
s.d_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Spiny Dogfish")
s.d_plot

### GAM
dc_fall_s.d <- filter(dc_fall, svspp == 'Spiny.Dogfish')
dc_fall_s.d$z_bt <- (dc_fall_s.d$bt - mean(dc_fall_s.d$bt))/(dc_fall_s.d$bt_sd)
mod18 <- gam(cb ~ s(year) + s(bt), data=dc_fall_s.d, family = 'gaussian')
summary(mod18)
mod18$aic
plot(mod18, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod18)

## Winter Skate
dc_fall_w.s <- filter(dc_fall, svspp == 'Winter.Skate')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_w.s$cb)/sqrt(length(dc_fall_w.s$cb)))

dc_fall_w.s$lower <- dc_fall_w.s$cb - se # calculating lower bound
dc_fall_w.s$upper <- dc_fall_w.s$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_w.s) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_w.s$year, y = dc_fall_w.s$cb, dc_fall_w.s$lower, dc_fall_w.s$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod2 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Winter Skate")

ggsave("winter skate fall GAM.png", plot = mod2, width = 6, height = 5, units = "in")

# plot data
w.s_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Winter Skate")
w.s_plot

### GAM
dc_fall_w.s <- filter(dc_fall, svspp == 'Winter.Skate')
dc_fall_w.s$z_bt <- (dc_fall_w.s$bt - mean(dc_fall_w.s$bt))/(dc_fall_w.s$bt_sd)
mod19 <- gam(cb ~ s(year) + s(z_bt), data=dc_fall_w.s, family = 'gaussian')
summary(mod19)
mod19$aic
plot(mod19, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod19)

## Pollock
dc_fall_p <- filter(dc_fall, svspp == 'Pollock')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_p$cb)/sqrt(length(dc_fall_p$cb)))

dc_fall_p$lower <- dc_fall_p$cb - se # calculating lower bound
dc_fall_p$upper <- dc_fall_p$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_p) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_p$year, y = dc_fall_p$cb, dc_fall_p$lower, dc_fall_p$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod3 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Pollock")

ggsave("pollock fall GAM.png", plot = mod3, width = 6, height = 5, units = "in")

# plot data
p_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Pollock")
p_plot

### GAM
dc_fall_p <- filter(dc_fall, svspp == 'Pollock')
dc_fall_p$z_bt <- (dc_fall_p$bt - mean(dc_fall_p$bt))/(dc_fall_p$bt_sd)
mod20 <- gam(cb ~ s(year) + s(bt), data=dc_fall_p, family = 'gaussian')
summary(mod20)
mod20$aic
plot(mod20, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod20)

## White Hake
dc_fall_w.h <- filter(dc_fall, svspp == 'White.Hake')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_w.h$cb)/sqrt(length(dc_fall_w.h$cb)))

dc_fall_w.h$lower <- dc_fall_w.h$cb - se # calculating lower bound
dc_fall_w.h$upper <- dc_fall_w.h$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_w.h) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_w.h$year, y = dc_fall_w.h$cb, dc_fall_w.h$lower, dc_fall_w.h$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod4 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "White Hake")

ggsave("white hake fall GAM.png", plot = mod4, width = 6, height = 5, units = "in")

# plot data
w.h_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("White Hake")
w.h_plot

### GAM
dc_fall_w.h <- filter(dc_fall, svspp == 'White.Hake')
dc_fall_w.h$z_bt <- (dc_fall_w.h$bt - mean(dc_fall_w.h$bt))/(dc_fall_w.h$bt_sd)
mod21 <- gam(cb ~ s(year) + s(z_bt), data=dc_fall_w.h, family = 'gaussian')
summary(mod21)
mod21$aic
plot(mod21, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod21)

## Red Hake
dc_fall_r.h <- filter(dc_fall, svspp == 'Red.Hake')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_r.h$cb)/sqrt(length(dc_fall_r.h$cb)))

dc_fall_r.h$lower <- dc_fall_r.h$cb - se # calculating lower bound
dc_fall_r.h$upper <- dc_fall_r.h$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_r.h) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_r.h$year, y = dc_fall_r.h$cb, dc_fall_r.h$lower, dc_fall_r.h$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod5 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Red Hake")

ggsave("red hake fall GAM.png", plot = mod5, width = 6, height = 5, units = "in")

# plot data
r.h_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Red Hake")
r.h_plot

### GAM
dc_fall_r.h <- filter(dc_fall, svspp == 'Red.Hake')
dc_fall_r.h$z_bt <- (dc_fall_r.h$bt - mean(dc_fall_r.h$bt))/(dc_fall_r.h$bt_sd)
mod22 <- gam(cb ~ s(year) + s(bt), data=dc_fall_r.h, family = 'gaussian')
summary(mod22)
mod22$aic
plot(mod22, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod22)

## Goosefish
dc_fall_g <- filter(dc_fall, svspp == 'Goosefish')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_g$cb)/sqrt(length(dc_fall_g$cb)))

dc_fall_g$lower <- dc_fall_g$cb - se # calculating lower bound
dc_fall_g$upper <- dc_fall_g$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_g) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_g$year, y = dc_fall_g$cb, dc_fall_g$lower, dc_fall_g$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod6 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Goosefish")

ggsave("goosefish fall GAM.png", plot = mod6, width = 6, height = 5, units = "in")

# plot data
g_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Goosefish")
g_plot

### GAM
dc_fall_g <- filter(dc_fall, svspp == 'Goosefish')
dc_fall_g$z_bt <- (dc_fall_g$bt - mean(dc_fall_g$bt))/(dc_fall_g$bt_sd)
mod23 <- gam(cb ~ s(year) + s(z_bt), data=dc_fall_g, family = 'gaussian')
summary(mod23)
mod23$aic
plot(mod23, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod23)

## Thorny Skate
dc_fall_t.s <- filter(dc_fall, svspp == 'Thorny.Skate')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_t.s$cb)/sqrt(length(dc_fall_t.s$cb)))

dc_fall_t.s$lower <- dc_fall_t.s$cb - se # calculating lower bound
dc_fall_t.s$upper <- dc_fall_t.s$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_t.s) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_t.s$year, y = dc_fall_t.s$cb, dc_fall_t.s$lower, dc_fall_t.s$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod7 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Thorny Skate")

ggsave("thorny skate fall GAM.png", plot = mod7, width = 6, height = 5, units = "in")

# plot data
t.s_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Thorny Skate")
t.s_plot

### GAM
dc_fall_t.s <- filter(dc_fall, svspp == 'Thorny.Skate')
dc_fall_t.s$z_bt <- (dc_fall_t.s$bt - mean(dc_fall_t.s$bt))/(dc_fall_t.s$bt_sd)
mod24 <- gam(cb ~ s(year) + s(bt), data=dc_fall_t.s, family = 'gaussian')
summary(mod24)
mod24$aic
plot(mod24, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod24)

## Atlantic Cod
dc_fall_a.c <- filter(dc_fall, svspp == 'Atlantic.Cod')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_a.c$cb)/sqrt(length(dc_fall_a.c$cb)))

dc_fall_a.c$lower <- dc_fall_a.c$cb - se # calculating lower bound
dc_fall_a.c$upper <- dc_fall_a.c$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_a.c) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_a.c$year, y = dc_fall_a.c$cb, dc_fall_a.c$lower, dc_fall_a.c$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod8 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Atlantic Cod")

ggsave("atlantic cod fall GAM.png", plot = mod8, width = 6, height = 5, units = "in")

# plot data
a.c_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Atlantic Cod")
a.c_plot

### GAM
dc_fall_a.c <- filter(dc_fall, svspp == 'Atlantic.Cod')
dc_fall_a.c$z_bt <- (dc_fall_a.c$bt - mean(dc_fall_a.c$bt))/(dc_fall_a.c$bt_sd)
mod25 <- gam(cb ~ s(year) + s(bt), data=dc_fall_a.c, family = 'gaussian')
summary(mod25)
mod25$aic
plot(mod25, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod25)

## Sea Raven
dc_fall_s.r <- filter(dc_fall, svspp == 'Sea.Raven')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_s.r$cb)/sqrt(length(dc_fall_s.r$cb)))

dc_fall_s.r$lower <- dc_fall_s.r$cb - se # calculating lower bound
dc_fall_s.r$upper <- dc_fall_s.r$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_s.r) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_s.r$year, y = dc_fall_s.r$cb, dc_fall_s.r$lower, dc_fall_s.r$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod9 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Sea Raven")

ggsave("sea raven fall GAM.png", plot = mod9, width = 6, height = 5, units = "in")

# plot data
s.r_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Sea Robin")
s.r_plot

### GAM
dc_fall_s.r <- filter(dc_fall, svspp == 'Sea.Raven')
dc_fall_s.r$z_bt <- (dc_fall_s.r$bt - mean(dc_fall_s.r$bt))/(dc_fall_s.r$bt_sd)
mod26 <- gam(cb ~ s(year) + s(z_bt), data=dc_fall_s.r, family = 'gaussian')
summary(mod26)
mod26$aic
plot(mod26, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod26)

## Silver Hake
dc_fall_s.h <- filter(dc_fall, svspp == 'Silver.Hake')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_s.h$cb)/sqrt(length(dc_fall_s.h$cb)))

dc_fall_s.h$lower <- dc_fall_s.h$cb - se # calculating lower bound
dc_fall_s.h$upper <- dc_fall_s.h$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_s.h) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_s.h$year, y = dc_fall_s.h$cb, dc_fall_s.h$lower, dc_fall_s.h$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod10 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Silver Hake")

ggsave("silver hake fall GAM.png", plot = mod10, width = 6, height = 5, units = "in")

# plot data
s.h_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Silver Hake")
s.h_plot

### GAM
dc_fall_s.h <- filter(dc_fall, svspp == 'Silver.Hake')
dc_fall_s.h$z_bt <- (dc_fall_s.h$bt - mean(dc_fall_s.h$bt))/(dc_fall_s.h$bt_sd)
mod27 <- gam(cb ~ s(year) + s(z_bt), data=dc_fall_s.h, family = 'gaussian')
summary(mod27)
mod27$aic
plot(mod27, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod27)

## Summer Flounder
dc_fall_s.f <- filter(dc_fall, svspp == 'Summer.Flounder')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_s.f$cb)/sqrt(length(dc_fall_s.f$cb)))

dc_fall_s.f$lower <- dc_fall_s.f$cb - se # calculating lower bound
dc_fall_s.f$upper <- dc_fall_s.f$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_s.f) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_s.f$year, y = dc_fall_s.f$cb, dc_fall_s.f$lower, dc_fall_s.f$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod11 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Summer Flounder")

ggsave("summer flounder fall GAM.png", plot = mod11, width = 6, height = 5, units = "in")

# plot data
s.f_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Summer Flounder")
s.f_plot

### GAM
dc_fall_s.f <- filter(dc_fall, svspp == 'Summer.Flounder')
dc_fall_s.f$z_bt <- (dc_fall_s.f$bt - mean(dc_fall_s.f$bt))/(dc_fall_s.f$bt_sd)
mod28 <- gam(cb ~ s(year) + s(bt), data=dc_fall_s.f, family = 'gaussian')
summary(mod28)
mod28$aic
plot(mod28, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod28)

## Fourspot Flounder
dc_fall_f.f <- filter(dc_fall, svspp == 'Fourspot.Flounder')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_f.f$cb)/sqrt(length(dc_fall_f.f$cb)))

dc_fall_f.f$lower <- dc_fall_f.f$cb - se # calculating lower bound
dc_fall_f.f$upper <- dc_fall_f.f$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_f.f) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_f.f$year, y = dc_fall_f.f$cb, dc_fall_f.f$lower, dc_fall_f.f$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod12 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Fourspot Flounder")

ggsave("fourspot flounder fall GAM.png", plot = mod12, width = 6, height = 5, units = "in")

# plot data
f.f_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Fourspot Flounder")
f.f_plot

### GAM
dc_fall_f.f <- filter(dc_fall, svspp == 'Fourspot.Flounder')
dc_fall_f.f$z_bt <- (dc_fall_f.f$bt - mean(dc_fall_f.f$bt))/(dc_fall_f.f$bt_sd)
mod29 <- gam(cb ~ s(year) + s(z_bt), data=dc_fall_f.f, family = 'gaussian')
summary(mod29)
mod29$aic
plot(mod29, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod29)

## Windowpane
dc_fall_w <- filter(dc_fall, svspp == 'Windowpane')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_w$cb)/sqrt(length(dc_fall_w$cb)))

dc_fall_w$lower <- dc_fall_w$cb - se # calculating lower bound
dc_fall_w$upper <- dc_fall_w$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_w) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_w$year, y = dc_fall_w$cb, dc_fall_w$lower, dc_fall_w$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod13 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Windowpane")

ggsave("windowpane fall GAM.png", plot = mod13, width = 6, height = 5, units = "in")

# plot data
w_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Windowpane")
w_plot

### GAM
dc_fall_w <- filter(dc_fall, svspp == 'Windowpane')
dc_fall_w$z_bt <- (dc_fall_w$bt - mean(dc_fall_w$bt))/(dc_fall_w$bt_sd)
mod30 <- gam(cb ~ s(year) + s(z_bt), data=dc_fall_w, family = 'gaussian')
summary(mod30)
mod30$aic
plot(mod30, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod30)

## Buckler Dory
dc_fall_b.d <- filter(dc_fall, svspp == 'Buckler.Dory')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_b.d$cb)/sqrt(length(dc_fall_b.d$cb)))

dc_fall_b.d$lower <- dc_fall_b.d$cb - se # calculating lower bound
dc_fall_b.d$upper <- dc_fall_b.d$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year, k = 1), data = dc_fall_b.d) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_b.d$year, y = dc_fall_b.d$cb, dc_fall_b.d$lower, dc_fall_b.d$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod14 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Buckler Dory")

ggsave("buckler dory fall GAM.png", plot = mod14, width = 6, height = 5, units = "in")

# plot data
b.d_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Buckler Dory")
b.d_plot

### GAM
dc_fall_b.d <- filter(dc_fall, svspp == 'Buckler.Dory')
dc_fall_b.d$z_bt <- (dc_fall_b.d$bt - mean(dc_fall_b.d$bt))/(dc_fall_b.d$bt_sd)
mod31 <- gam(cb ~ s(year), data=dc_fall_b.d, family = 'gaussian')
mod31.a <- gam(cb ~ s(z_bt), data=dc_fall_b.d, family = 'gaussian')
summary(mod31)
summary(mod31.a)
mod31$aic
mod31.a$aic
plot(mod31, pages=1, scheme=1, unconditional=TRUE)
plot(mod31.a, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod14)
gam.check(mod14.a)

## Blue Fish
dc_fall_b.f <- filter(dc_fall, svspp == 'Blue.Fish')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_b.f$cb)/sqrt(length(dc_fall_b.f$cb)))

dc_fall_b.f$lower <- dc_fall_b.f$cb - se # calculating lower bound
dc_fall_b.f$upper <- dc_fall_b.f$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_b.f) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_b.f$year, y = dc_fall_b.f$cb, dc_fall_b.f$lower, dc_fall_b.f$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod15 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Blue Fish")

ggsave("blue fish fall GAM.png", plot = mod15, width = 6, height = 5, units = "in")

# plot data
b.f_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Blue Fish")
b.f_plot

### GAM
dc_fall_b.f <- filter(dc_fall, svspp == 'Blue.Fish')
dc_fall_b.f$z_bt <- (dc_fall_b.f$bt - mean(dc_fall_b.f$bt))/(dc_fall_b.f$bt_sd)
mod32 <- gam(cb ~ s(year), data=dc_fall_b.f, family = 'gaussian')
mod32.a <- gam(cb ~ s(z_bt), data=dc_fall_b.f, family = 'gaussian')
summary(mod32)
summary(mod32.a)
mod32$aic
mod32.a$aic
plot(mod32, pages=1, scheme=1, unconditional=TRUE)
plot(mod32.a, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod32)
gam.check(mod32.a)

## Longhorn Sculpin
dc_fall_l.s <- filter(dc_fall, svspp == 'Longhorn.Sculpin')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_l.s$cb)/sqrt(length(dc_fall_l.s$cb)))

dc_fall_l.s$lower <- dc_fall_l.s$cb - se # calculating lower bound
dc_fall_l.s$upper <- dc_fall_l.s$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_l.s) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_l.s$year, y = dc_fall_l.s$cb, dc_fall_l.s$lower, dc_fall_l.s$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod16 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Longhorn Sculpin")

ggsave("longhorn sculpin fall GAM.png", plot = mod16, width = 6, height = 5, units = "in")

# plot data
l.s_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Longhorn Sculpin")
l.s_plot

### GAM
dc_fall_l.s <- filter(dc_fall, svspp == 'Longhorn.Sculpin')
dc_fall_l.s$z_bt <- (dc_fall_l.s$bt - mean(dc_fall_l.s$bt))/(dc_fall_l.s$bt_sd)
mod33 <- gam(cb ~ s(year) + s(z_bt), data=dc_fall_l.s, family = 'gaussian')
summary(mod33)
mod33$aic
plot(mod33, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod33)

## Striped Searobin
dc_fall_s.s <- filter(dc_fall, svspp == 'Striped.Searobin')

# calculating individual confidence intervals for points
se <- (sd(dc_fall_s.s$cb)/sqrt(length(dc_fall_s.s$cb)))

dc_fall_s.s$lower <- dc_fall_s.s$cb - se # calculating lower bound
dc_fall_s.s$upper <- dc_fall_s.s$cb + se # calculating upper bound

# calculating confidence interval for the gam fit line itself
model <- gam(cb ~ s(year), data = dc_fall_s.s) # make predictions with confidence intervals 
model$aic
predictions <- predict(model, type = "response", se.fit = TRUE)

predicted_values <- predictions$fit   # Extract predicted values and confidence intervals
predicted_values_d <- as.data.frame(predicted_values)
lower_ci <- predicted_values - 1.96 * predictions$se.fit
lower_ci_d <- as.data.frame(lower_ci)
upper_ci <- predicted_values + 1.96 * predictions$se.fit
upper_ci_d <- as.data.frame(upper_ci)

# combine everything into a data frame
result_df <- data.frame(x = dc_fall_s.s$year, y = dc_fall_s.s$cb, dc_fall_s.s$lower, dc_fall_s.s$upper,predicted_values_d, lower_ci_d, upper_ci_d)
new_col_names <- c("year", "cb", "lower","upper","predicted_values_d", "lower_ci_d", "upper_ci_d")
colnames(result_df) <- new_col_names

mod17 <- ggplot(result_df, aes(x = year, y = cb)) +
  geom_line(aes(y = predicted_values), color = "darkslategray4") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "darkslategray4", alpha = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point() +
  labs(title = "Striped Searobin")

ggsave("striped searobin fall GAM.png", plot = mod17, width = 6, height = 5, units = "in")

# plot data
s.s_plot <- ggplot(data = result_df, aes(x = year, y = cb)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "darkorange3") +
  geom_point(color = "black", size = 1) +
  ecodata::geom_gls() +
  theme_bw() +
  theme(plot.title = element_text(size = 13)) +
  xlab("Year") +
  ylab("Consumption:Biomass") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 10),
        axis.title.x=element_text(size = 10)) +
  ggtitle("Striped Searobin")
s.s_plot

### GAM
dc_fall_s.s <- filter(dc_fall, svspp == 'Striped.Searobin')
dc_fall_s.s$z_bt <- (dc_fall_s.s$bt - mean(dc_fall_s.s$bt))/(dc_fall_s.s$bt_sd)
mod34 <- gam(cb ~ s(year), data=dc_fall_s.s, family = 'gaussian')
mod34.a <- gam(cb ~ s(bt), data=dc_fall_s.s, family = 'gaussian')
summary(mod34)
summary(mod34.a)
mod34$aic
mod34.a$aic
plot(mod34, pages=1, scheme=1, unconditional=TRUE)
plot(mod34.a, pages=1, scheme=1, unconditional=TRUE)
gam.check(mod34)
gam.check(mod34.a)

# combining all of the linear plots into one page
glm_cb_f <- ggarrange(s.d_plot, w.s_plot, p_plot, w.h_plot, r.h_plot,
                    g_plot, t.s_plot, a.c_plot, s.r_plot, s.h_plot,
                    s.f_plot, f.f_plot, w_plot, b.d_plot, b.f_plot,
                    l.s_plot, s.s_plot, ncol = 5, nrow = 4)
glm_cb_f

ggsave("glms fall.png", plot = glm_cb_f, width = 15, height = 7, units = "in")

# visually comparing GAMS across species
## GAM of year v. cb --> I save each of these figures (total of 2 figures) once 
## the matrix (3x3) is filled in order to have two large figures that I can compare all species with
par(mfrow=c(3,3))
mod18 <- gam(cb ~ s(year), data=dc_fall_s.d, family = 'gaussian')
plot(mod18, scheme=1, unconditional=TRUE, main = "Spiny Dogfish") 
mod19 <- gam(cb ~ s(year), data=dc_fall_w.s, family = 'gaussian')
plot(mod19, scheme=1, unconditional=TRUE, main = "Winter Skate") 
mod20 <- gam(cb ~ s(year), data=dc_fall_p, family = 'gaussian')
plot(mod20, scheme=1, unconditional=TRUE, main = "Pollock") 
mod21 <- gam(cb ~ s(year), data=dc_fall_w.h, family = 'gaussian')
plot(mod21, scheme=1, unconditional=TRUE, main = "White Hake")
mod22 <- gam(cb ~ s(year), data=dc_fall_r.h, family = 'gaussian')
plot(mod22, scheme=1, unconditional=TRUE, main = "Red Hake")
mod23 <- gam(cb ~ s(year), data=dc_fall_g, family = 'gaussian')
plot(mod23, scheme=1, unconditional=TRUE, main = "Goosefish")
mod24 <- gam(cb ~ s(year), data=dc_fall_t.s, family = 'gaussian')
plot(mod24, scheme=1, unconditional=TRUE, main = "Thorny Skate")
mod25 <- gam(cb ~ s(year), data=dc_fall_a.c, family = 'gaussian')
plot(mod25, scheme=1, unconditional=TRUE, main = "Atlantic Cod")
mod26 <- gam(cb ~ s(year), data=dc_fall_s.r, family = 'gaussian')
plot(mod26, scheme=1, unconditional=TRUE, main = "Sea Raven")
mod27 <- gam(cb ~ s(year), data=dc_fall_s.h, family = 'gaussian')
plot(mod27, scheme=1, unconditional=TRUE, main = "Silver Hake")
mod28 <- gam(cb ~ s(year), data=dc_fall_s.f, family = 'gaussian')
plot(mod28, scheme=1, unconditional=TRUE, main = "Summer Flounder")
mod29 <- gam(cb ~ s(year), data=dc_fall_f.f, family = 'gaussian')
plot(mod29, scheme=1, unconditional=TRUE, main = "Fourspot Flounder")
mod30 <- gam(cb ~ s(year), data=dc_fall_w, family = 'gaussian')
plot(mod30, scheme=1, unconditional=TRUE, main = "Windowpane")
mod31 <- gam(cb ~ s(year), data=dc_fall_b.d, family = 'gaussian')
plot(mod31, scheme=1, unconditional=TRUE, main = "Buckler Dory")
mod32 <- gam(cb ~ s(year), data=dc_fall_b.f, family = 'gaussian')
plot(mod32, scheme=1, unconditional=TRUE, main = "Blue Fish")
mod33 <- gam(cb ~ s(year), data=dc_fall_l.s, family = 'gaussian')
plot(mod33, scheme=1, unconditional=TRUE, main = "Longhorn Sculpin")
mod34 <- gam(cb ~ s(year), data=dc_fall_s.s, family = 'gaussian')
plot(mod34, scheme = 1, unconditional=TRUE, main = "Striped Searobin")

## GAM of bt v. cb --> same note as above
par(mfrow = c(3,3))
mod18 <- gam(cb ~ s(bt), data=dc_spring_s.d, family = 'gaussian')
plot(mod18, scheme=1, unconditional=TRUE, main = "Spiny Dogfish")
mod19 <- gam(cb ~ s(bt), data=dc_spring_w.s, family = 'gaussian')
plot(mod19, scheme=1, unconditional=TRUE, main = "Winter Skate")
mod20 <- gam(cb ~ s(bt), data=dc_spring_p, family = 'gaussian')
plot(mod20, scheme=1, unconditional=TRUE, main = "Pollock")
mod21 <- gam(cb ~ s(bt), data=dc_spring_w.h, family = 'gaussian')
plot(mod21, scheme=1, unconditional=TRUE, main = "White Hake")
mod22 <- gam(cb ~ s(bt), data=dc_spring_r.h, family = 'gaussian')
plot(mod22, scheme=1, unconditional=TRUE, main = " Red Hake")
mod23 <- gam(cb ~ s(bt), data=dc_spring_g, family = 'gaussian')
plot(mod23, scheme=1, unconditional=TRUE, main = "Goosefish")
mod24 <- gam(cb ~ s(bt), data=dc_spring_t.s, family = 'gaussian')
plot(mod24, scheme=1, unconditional=TRUE, main = "Thorny Skate")
mod25 <- gam(cb ~ s(bt), data=dc_spring_a.c, family = 'gaussian')
plot(mod25, scheme=1, unconditional=TRUE, main = "Atlantic Cod")
mod26 <- gam(cb ~ s(bt), data=dc_spring_s.r, family = 'gaussian')
plot(mod26, scheme=1, unconditional=TRUE, main = "Sea Raven")
mod27 <- gam(cb ~ s(bt), data=dc_spring_s.h, family = 'gaussian')
plot(mod27, scheme=1, unconditional=TRUE, main = "Silver Hake")
mod28 <- gam(cb ~ s(bt), data=dc_spring_s.f, family = 'gaussian')
plot(mod28, scheme=1, unconditional=TRUE, main = "Summer Flounder")
mod29 <- gam(cb ~ s(bt), data=dc_spring_f.f, family = 'gaussian')
plot(mod29, scheme=1, unconditional=TRUE, main = "Fourspot Flounder")
mod30 <- gam(cb ~ s(bt), data=dc_spring_w, family = 'gaussian')
plot(mod30, scheme=1, unconditional=TRUE, main = "Windowpane")
mod31 <- gam(cb ~ s(bt) , data=dc_spring_b.d, family = 'gaussian')
plot(mod31, scheme=1, unconditional=TRUE, main = "Buckler Dory")
mod32 <- gam(cb ~ s(bt), data=dc_spring_b.f, family = 'gaussian')
plot(mod32, scheme=1, unconditional=TRUE, main = "Blue Fish")
mod33 <- gam(cb ~ s(z_bt), data=dc_spring_l.s, family = 'gaussian')
plot(mod33, scheme=1, unconditional=TRUE, main = "Longhorn Sculpin")
mod34 <- gam(cb ~ s(z_bt), data=dc_spring_s.s, family = 'gaussian')
plot(mod34, scheme = 1, unconditional=TRUE, main = "Striped Searobin")

# filtering data to examine bottom temperature trend
dc_temp <- data_comb %>% 
  select('year','bt') %>% 
  na.omit('year','bt')

## bottom temperature trend across spring and fall
bt_plot <- ggplot(data = dc_temp, 
       aes(x = year, y = bt)) +
  geom_point(color = "black", size = 1) +
  geom_smooth(method = lm, color = "purple2", linewidth = 2) +
  theme_bw() +
  theme(plot.title = element_text(size = 15)) +
  xlab("Year") +
  ylab("Bottom Temperature") +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 13),
        axis.title.x=element_text(size = 13)) +
  ggtitle("Bottom Temperature Over Time")
bt_plot

## bottom temperature trend across fall
dc_fall <- filter(data_comb_f, seacat == 'FALL')

bt_fall_plot <- ggplot(data = dc_fall, 
                  aes(x = year, y = bt)) +
  geom_point(color = "black", size = 1) +
  geom_smooth(method = lm, color = "purple2", linewidth = 2, se = FALSE) +
  theme_bw() +
  theme(plot.title = element_text(size = 15)) +
  xlab("Year") +
  ylab("Bottom Temperature") +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 13),
        axis.title.x=element_text(size = 13)) +
  ggtitle("Fall Bottom Temperature Over Time")
bt_fall_plot

## bottom temperature trend across spring 
dc_spring <- filter(data_comb_f, seacat == 'SPRING')

bt_spring_plot <- ggplot(data = dc_spring, aes(x = year, y = bt)) +
  geom_point(color = "black", size = 1) +
  geom_smooth(method = lm, color = "purple2", linewidth = 2, se = FALSE) +
  theme_bw() +
  theme(plot.title = element_text(size = 15)) +
  xlab("Year") +
  ylab("Bottom Temperature") +
  theme(axis.text.x=element_text(size=10), 
        axis.text.y = element_text(size=10),
        axis.title.y=element_text(size = 13),
        axis.title.x=element_text(size = 13)) +
  ggtitle("Spring Bottom Temperature Over Time")
bt_spring_plot
