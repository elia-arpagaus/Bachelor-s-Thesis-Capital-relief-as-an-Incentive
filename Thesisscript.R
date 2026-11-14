###############################################################################
## Capital Relief as an Incentive
##
## This code is used for the calculations, the tables and the figures of the
## bachelor's thesis. In order of the sections / chapters.
##
## Elia E. Arpagaus
## Fall 2026
###############################################################################

# Setup ----
library(ggplot2)

theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.title = element_blank()))
update_geom_defaults("line", list(linewidth = 1.2))
update_geom_defaults("step", list(linewidth = 1.2))
update_geom_defaults("function", list(linewidth = 1.2))
green <- "#00802e"
blue <- "#73a5ad"
red <- "#eb696a"
grey_light <- "#e2d8c4"
grey <- "#9e978a"
grey_dark <- "#59574c"


# Baseline values ----
mu <- 1
sigma <- 1
delta <- 0.06
D <- mu / sigma
beta <- 0.5
c_cost <- 0.01



# 2.2 Risk measures and solvency rule ----
k_var <- function(alpha) qnorm(1 - alpha)
k_es <- function(alpha) dnorm(qnorm(1 - alpha)) / alpha
k_sii <- k_var(0.005)
k_sst <- k_es(0.01)
C0 <- function(k = k_sii, dist = D, s = sigma) s * (k - dist)
print(c(k_sii = k_sii, k_sst = k_sst, C0_sii = C0()))



# 2.3 Stringency of the European regimes ----
alphas <- c(0.005, 0.01, 0.025, 0.05)
print(data.frame(alpha = 100 * alphas, k_ES = k_es(alphas), equivalent_VaR_level = 100 * pnorm(-k_es(alphas))))

a_sii_es <- uniroot(function(al) k_es(al) - k_sii, c(0.001, 0.1), tol = 1e-10)$root
a_sst_var <- pnorm(-k_sst)
print(round(c(k_sii = k_sii, k_sst = k_sst, es_level_of_sii = a_sii_es, var_level_of_sst = a_sst_var), 4))
fig_stringency <- ggplot() +
  geom_function(aes(colour = "VaR"), fun = k_var, n = 500) +
  geom_function(aes(colour = "ES"), fun = k_es, n = 500) +
  annotate("segment", x = c(0.005, a_sst_var), xend = c(a_sii_es, 0.01), y = c(k_sii, k_sst), yend = c(k_sii, k_sst),
           linetype = "dashed", linewidth = 0.4) +
  annotate("segment", x = c(a_sii_es, a_sst_var), xend = c(a_sii_es, a_sst_var), y = -Inf, yend = c(k_sii, k_sst),
           linetype = "dotted", linewidth = 0.4) +
  annotate("point", x = c(0.005, 0.01), y = c(k_sii, k_sst), colour = c(red, green), size = 1.8) +
  annotate("text", x = c(0.005, 0.01), y = c(k_sii, k_sst), label = c("Solvency II", "SST"), size = 3, hjust = c(1, -0.15), vjust = c(1.8, -0.8)) +
  annotate("text", x = c(a_sii_es, a_sst_var), y = 1.4, label = sprintf("%.2f %%", 100 * c(a_sii_es, a_sst_var)), size = 2.8, hjust = c(-0.1, 1.1)) +
  scale_colour_manual(values = c(VaR = red, ES = green), breaks = c("VaR", "ES")) +
  scale_x_log10(limits = c(0.001, 0.1), breaks = c(0.001, 0.005, 0.01, 0.025, 0.05, 0.1), labels = function(al) paste0(100 * al, " %")) +
  labs(x = expression(Tail~level~alpha~(log~scale)), y = expression(Stringency~k[R](alpha)))
print(fig_stringency)
cat("The probability of default without capital is", 100 * pnorm(-D), "%\n")



# 2.4 Signal ----
kappa <- function(b, cc = c_cost, s = sigma) cc / s^2 * b / (1 - b)



# 2.5 Updating and timeline ----
s_post <- sigma * sqrt(1 - beta)
fig_posterior <- ggplot() +
  geom_function(aes(colour = "prior"), fun = dnorm, args = list(mean = mu, sd = sigma), n = 500) +
  geom_function(aes(colour = "good"), fun = dnorm, args = list(mean = mu + beta * sigma, sd = s_post), n = 500) +
  geom_function(aes(colour = "bad"), fun = dnorm, args = list(mean = mu - 2 * beta * sigma, sd = s_post), n = 500) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_colour_manual(values = c(prior = "black", good = green, bad = red), breaks = c("prior", "good", "bad"),
                      labels = c("Prior", expression(Good~signal), expression(Bad~signal))) +
  xlim(-3, 5) +
  labs(x = "Wealth W", y = "Density")
print(fig_posterior)


# 3.1 Capital relief ----
a <- function(b, k = k_sii, dist = D) (dist - k * sqrt(1 - b)) / sqrt(b)
betas <- c(0.1, 0.3, 0.5, 0.7, 0.9)
Ds <- c(0.5, 1, 1.5, 2, 2.5)
print(round(100 * setNames(pnorm(k_sii * (1 - sqrt(1 - betas)) / sqrt(betas)), betas), 1)) 
no_capital <- 100 * outer(betas, Ds, function(b, dd) pnorm(a(b, k_sii, dd)))
dimnames(no_capital) <- list(beta = betas, D = Ds); print(round(no_capital, 1)) 


# 3.2 Value of information ----
EC <- function(b, k = k_sii, dist = D, s = sigma) s * sqrt(b) * (dnorm(a(b, k, dist)) - a(b, k, dist) * pnorm(-a(b, k, dist)))
EC_bar <- function(b, k = k_sii, dist = D, s = sigma) s * (k * sqrt(1 - b) - dist)
print(c(C0 = C0(), EC = EC(beta), EC_bar = EC_bar(beta), EC_0.9 = EC(0.9), EC_bar_0.9 = EC_bar(0.9))) 

V_S <- function(b, k = k_sii, dist = D, s = sigma, d = delta) d * s * ((k - dist) - sqrt(b) * (dnorm(a(b, k, dist)) - a(b, k, dist) * pnorm(-a(b, k, dist))))
V_bar <- function(b, k = k_sii, s = sigma, d = delta) d * k * s * (1 - sqrt(1 - b))

bg <- seq(0.001, 0.99, length.out = 1000)
kg <- seq(D, 4, length.out = 300)
Dg <- seq(0, k_sii, length.out = 300)
panel_names <- c("Informativeness~beta", "Stringency~k[R]", "Distance~to~default~D")
panels <- rbind(data.frame(panel = panel_names[1], x = bg, floor = V_S(bg), nofloor = V_bar(bg)),
                data.frame(panel = panel_names[2], x = kg, floor = V_S(beta, k = kg), nofloor = V_bar(beta, k = kg)),
                data.frame(panel = panel_names[3], x = Dg, floor = V_S(beta, dist = Dg), nofloor = V_bar(beta)))
panels$panel <- factor(panels$panel, panel_names)
baseline <- data.frame(panel = factor(panel_names, panel_names), x = c(beta, k_sii, D), V = V_S(beta))
fig_theorem1 <- ggplot(panels, aes(x)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(aes(y = nofloor, colour = "without the floor", linetype = "without the floor")) +
  geom_line(aes(y = floor, colour = "with the floor", linetype = "with the floor")) +
  geom_point(data = baseline, aes(y = V), size = 1.5) +
  facet_wrap(~ panel, scales = "free_x", labeller = label_parsed, nrow = 1) +
  scale_colour_manual(values = c("with the floor" = "black", "without the floor" = grey)) +
  scale_linetype_manual(values = c("with the floor" = "solid", "without the floor" = "dashed")) +
  labs(x = NULL, y = expression(value~of~the~signal~V[S]), colour = NULL, linetype = NULL) +
  theme(panel.spacing.x = unit(1, "lines"))
print(fig_theorem1)

print(round(c(C0 = C0(), sigma_over_k = sigma / k_sii), 4))


# 3.3 Cost of information and chosen informativeness ----
net <- function(b, k = k_sii, dist = D, cc = c_cost, floor = TRUE) if (floor) V_S(b, k, dist) - kappa(b, cc) else V_bar(b, k) - kappa(b, cc)
b_grid <- seq(0.001, 0.999, by = 0.001)
beta_hat <- function(k = k_sii, dist = D, cc = c_cost, floor = TRUE) {
  i <- max(which(net(b_grid, k, dist, cc, floor) >= 0))
  uniroot(net, b_grid[i + 0:1], k = k, dist = dist, cc = cc, floor = floor, tol = 1e-10)$root
}
beta_star <- function(k = k_sii, dist = D, cc = c_cost, floor = TRUE) {
  i <- which.max(net(b_grid, k, dist, cc, floor))
  optimize(net, b_grid[i + c(-1, 1)], k = k, dist = dist, cc = cc, floor = floor, maximum = TRUE, tol = 1e-10)$maximum
}
b_hat <- beta_hat()
b_star <- beta_star()

curves <- data.frame(b = bg, V = V_S(bg), K = kappa(bg))
fig_cost <- ggplot(curves, aes(b)) +
  geom_ribbon(data = subset(curves, V >= K), aes(ymin = K, ymax = V), fill = grey_light) +
  geom_line(aes(y = V, colour = "value")) +
  geom_line(aes(y = K, colour = "cost")) +
  annotate("segment", x = b_star, xend = b_star, y = kappa(b_star), yend = V_S(b_star), colour = green, linewidth = 0.9) +
  annotate("segment", x = c(b_star, b_hat), xend = c(b_star, b_hat), y = -Inf, yend = c(kappa(b_star), V_S(b_hat)), linetype = "dotted", linewidth = 0.8) +
  annotate("point", x = c(b_star, b_hat), y = c(V_S(b_star), V_S(b_hat)), colour = c(green, "black"), size = 1.6) +
  annotate("text", x = b_star + 0.01, y = (kappa(b_star) + V_S(b_star)) / 2, label = "net gain", colour = green, hjust = 0, size = 3) +
  annotate("text", x = c(b_star, b_hat)+0.07, y = -0.003, label = sprintf(c("beta^'*' == %.2f", "hat(beta) == %.2f"), c(b_star, b_hat)),
           parse = TRUE, size = 3) +
  scale_colour_manual(values = c(value = "black", cost = red), breaks = c("value", "cost"),
                      labels = c(expression(Value~V[S](beta)), expression(Cost~kappa(beta)))) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.005, 0.12)) +
  labs(x = expression(Informativeness~beta), y = "Value and cost of the signal")
print(fig_cost)

print(round(c(beta_hat = b_hat, beta_star = b_star, beta_star_nofloor = beta_star(floor = FALSE)), 4))


# 3.4 Social value of information ----
alpha_fig <- 0.05
k_fig <- k_var(alpha_fig)
lgd_fig <- (dnorm(k_fig) / pnorm(-k_fig) - k_fig) * c(sigma, s_post)
xs <- seq(-2, 4, length.out = 1201)
lv <- c("Prior", "Posterior")
dens <- rbind(data.frame(x = xs, f = dnorm(xs, k_fig * sigma, sigma), what = lv[1]),
              data.frame(x = xs, f = dnorm(xs, k_fig * s_post, s_post), what = lv[2]))
fig_deficit <- ggplot(dens, aes(x, f, colour = what)) +
  geom_line() +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_colour_manual(values = c("black", green), breaks = lv) +
  scale_fill_manual(values = c("black", green), breaks = lv) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  annotate("segment", x = -lgd_fig, xend = -lgd_fig, y = 0, yend = Inf, colour = c("black", green), linetype = "dashed", linewidth = 0.7) +
  labs(x = expression(Wealth~plus~capital~W + C), y = "Density")
print(fig_deficit)

pd0 <- function(k = k_sii) pnorm(-k)
lgd0 <- function(k = k_sii, s = sigma) s * (dnorm(k) / pnorm(-k) - k)
epd_floor <- function(b, k = k_sii, dist = D, s = sigma) {
  t <- function(m) pmax(m, k * sqrt(1 - b)) / sqrt(1 - b)
  s * sqrt(1 - b) * integrate(function(m) dnorm(m, dist, sqrt(b)) * (dnorm(t(m)) - t(m) * pnorm(-t(m))), -Inf, Inf, rel.tol = 1e-8)$value
}

ks <- c(k_var(c(0.001, 0.005, 0.01, 0.025, 0.05)))
cs <- c(0.005, 0.01, 0.02, 0.03, 0.04)
deficit_change <- Vectorize(function(k, cc, floor) {
  b <- beta_star(k, D, cc, floor)
  epd <- if (floor) epd_floor(b, k, D) else sqrt(1 - b) * pd0(k) * lgd0(k)
  100 * (epd / (pd0(k) * lgd0(k)) - 1)
})
tab_deficit <- NULL
for (k in ks) tab_deficit <- rbind(tab_deficit, deficit_change(k, cs, TRUE), deficit_change(k, cs, FALSE))
dimnames(tab_deficit) <- list(paste(rep(round(ks, 2), each = 2), c("floor", "no floor")), c = cs)
print(round(tab_deficit, 1))


# 4.2 Model extension and intuition ----

Theta <- 5 * sigma
alpha_bar <- 0.005
p0 <- 0.004
p_after <- c(0.002, 0.006)
C_cat <- function(p, S = mu, b = beta, es = FALSE, alpha = alpha_bar) {
  m <- mu + b * (S - mu)
  s <- sigma * sqrt(1 - b)
  alpha_p <- pmax(alpha - p, 0) / (1 - p)
  w <- pmin(p / alpha, 1)
  ordinary <- (1 - p) / alpha * s * dnorm(k_var(alpha_p)) - (1 - w) * m 
  C_bar <- if (es) w * Theta + ordinary else pmin(Theta, s * k_var(alpha_p) - m)
  pmax(0, C_bar)
}
print(pnorm(-(Theta + mu) / sigma))

finding <- function(...) {
  C <- C_cat(p_after, ...)
  c(good = C[1], bad = C[2], expected = mean(C))
}
tab_cat <- rbind(VaR = finding(), ES = finding(es = TRUE))


S_fig <- mu + c(-2, 0, 1) * sigma
lv_rule <- c("Value-at-Risk (VaR)", "Expected Shortfall (ES)")
lv_signal <- c("bad", "neutral", "good")
lab_signal <- c(expression(Bad~Gaussian~signal), expression(Neutral~Gaussian~signal), expression(Good~Gaussian~signal))
cat_curves <- expand.grid(p = seq(0, 0.01, length.out = 1001), S = S_fig, es = c(FALSE, TRUE))
cat_curves$C <- with(cat_curves, ifelse(es, C_cat(p, S, es = TRUE), C_cat(p, S)))
cat_curves$rule <- factor(cat_curves$es, labels = lv_rule)
cat_curves$signal <- factor(cat_curves$S, labels = lv_signal)
chord <- data.frame(rule = factor(lv_rule, lv_rule), tab_cat)
fig_fat_tail <- ggplot(cat_curves, aes(p, C)) +
  geom_line(aes(colour = signal, linetype = signal, linewidth = signal)) +
  geom_segment(data = chord, aes(x = p_after[1], xend = p_after[2], y = good, yend = bad), linewidth = 1, colour = "black") +
  geom_point(data = chord, aes(x = p_after[1], y = good), size = 1.8, colour = "black") +
  geom_point(data = chord, aes(x = p_after[2], y = bad), size = 1.8, colour = "black") +
  geom_point(data = chord, aes(x = p0, y = expected, fill = "Expected capital after learning (neutral G. signal)"), shape = 23, size = 3, colour = "black", stroke = 1) +
  facet_wrap(~ rule) +
  scale_colour_manual(values = c(bad = red, neutral = "black", good = green), breaks = lv_signal, labels = lab_signal) +
  scale_fill_manual(values = c("Expected capital after learning (neutral G. signal)" = "#FEF24A")) +
  scale_linetype_manual(values = c(bad = "dashed", neutral = "solid", good = "dashed"), breaks = lv_signal, labels = lab_signal) +
  scale_x_continuous(breaks = seq(0, 0.01, 0.002), labels = function(p) paste0(100 * p, " %")) +
  scale_linewidth_manual(values = c(bad = 1, neutral = 1.5, good = 1), breaks = lv_signal, labels = lab_signal) +
  theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"), legend.margin = margin(0, 0, 0, 0)) +
  labs(x = "Probability of the catastrophe p", y = expression(Capital~required~C[Theta](S, p)))
print(fig_fat_tail)