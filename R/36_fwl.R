# 36: what the year effects take out of the static estimate
#     (Section 5.2)
#
# Adding year effects to the static within regression changes the
# coefficient on the public share. By the Frisch-Waugh-Lovell
# theorem the change is
#
#   b(country) - b(country, year) = sum_t n_t rbar_t g_t / sum_i r_i^2
#
# where r_i is the public share net of the other regressors and of
# the country effects, rbar_t its mean in year t, n_t the number of
# countries observed in year t, and g_t the year effects of the
# two-way fit. The year effects would move nothing if the public
# share had no common annual component (rbar_t = 0 for every t).
#
# Because the n_t rbar_t sum to zero, the total does not depend on
# which year the year effects are measured from, but its split into
# periods does. The split below measures them from their weighted
# mean, so that no single year serves as the reference.
#
# Reads  data/ through build_panel()
# Writes output/36_fwl.csv, 36_fwl_periods.csv

source("00_common.R")

X   <- KEY_X
OTH <- setdiff(X, "PUB")

fwl <- function(d) {
  f1 <- as.formula(paste("LE ~", paste(X, collapse = " + "), "+ factor(iso3c)"))
  m1 <- lm(f1, data = d)
  m2 <- lm(update(f1, . ~ . + factor(year)), data = d)
  gn <- grep("^factor\\(year\\)", names(coef(m2)), value = TRUE)
  g  <- c(0, coef(m2)[gn])
  names(g) <- c(min(d$year), sub("factor\\(year\\)", "", gn))
  r  <- resid(lm(as.formula(paste("PUB ~", paste(OTH, collapse = " + "),
                                  "+ factor(iso3c)")), data = d))
  rl <- resid(lm(as.formula(paste("LE ~", paste(OTH, collapse = " + "),
                                  "+ factor(iso3c)")), data = d))
  list(b1 = unname(coef(m1)["PUB"]), b2 = unname(coef(m2)["PUB"]),
       g = g, r = r, rl = rl, ss = sum(r^2))
}

dat <- build_panel(year_max = 2024, zero_rule = "multiplicative")
F   <- fwl(dat)
cat(sprintf("Panel: %d obs, %d countries\n", nrow(dat), n_distinct(dat$iso3c)))
cat(sprintf("static FE, country effects          : %+.4f\n", F$b1))
cat(sprintf("static FE, country and year effects : %+.4f\n", F$b2))
cat(sprintf("move                                : %+.4f\n", F$b1 - F$b2))
cat(sprintf("sum_t n_t rbar_t g_t / sum r^2      : %+.4f\n",
            sum(F$r * F$g[as.character(dat$year)]) / F$ss))
stopifnot(abs((F$b1 - F$b2) -
              sum(F$r * F$g[as.character(dat$year)]) / F$ss) < 1e-10)

# ---- 1. year by year ---------------------------------------
yt <- data.frame(year = dat$year, r = F$r, rl = F$rl,
                 pub_w = dat$PUB - ave(dat$PUB, dat$iso3c),
                 le_w  = dat$LE  - ave(dat$LE,  dat$iso3c)) %>%
  group_by(year) %>%
  summarise(n = n(), pub_net = mean(r), le_net = mean(rl),
            pub_within = mean(pub_w), le_within = mean(le_w),
            .groups = "drop") %>%
  mutate(g = unname(F$g[as.character(year)]),
         g_centred = g - weighted.mean(g, n),
         contribution = n * pub_net * g_centred / F$ss,
         contribution_from_first_year = n * pub_net * g / F$ss)
stopifnot(abs(sum(yt$contribution) - (F$b1 - F$b2)) < 1e-10)

cat("\n=== Common annual component of the public share and year effects ===\n")
cat("pub_net: public share net of the other regressors and the country\n")
cat("effects, mean by year (pp); g_centred: year effect in life expectancy\n")
cat("measured from its weighted mean (years).\n\n")
print(yt %>% transmute(year, n, pub_net = round(pub_net, 2),
                       g_centred = round(g_centred, 2),
                       contribution = round(contribution, 4)) %>%
        as.data.frame(), row.names = FALSE)

wc <- function(a, b, w) cov.wt(cbind(a, b), wt = w, cor = TRUE)$cor[1, 2]
cat(sprintf("\nweighted correlation of pub_net and the year effects   : %+.2f\n",
            wc(yt$pub_net, yt$g, yt$n)))
cat(sprintf("correlation of the raw within-country annual means\n"))
cat(sprintf("  of the public share and of life expectancy            : %+.2f\n",
            cor(yt$pub_within, yt$le_within)))
cat(sprintf("share of the variance of r explained by year means     : %.1f%%\n",
            100 * sum(yt$n * yt$pub_net^2) / F$ss))

# ---- 2. by period ------------------------------------------
per <- function(y) cut(y, c(-Inf, 2009, 2019, Inf),
                       labels = c("2000-2009", "2010-2019", "2020-2024"))
pt <- yt %>% mutate(period = per(year)) %>%
  group_by(period) %>%
  summarise(obs = sum(n), contribution = sum(contribution),
            contribution_from_first_year = sum(contribution_from_first_year),
            .groups = "drop") %>%
  mutate(period = as.character(period))
pt <- bind_rows(pt,
  pt %>% filter(period != "2020-2024") %>%
    summarise(period = "2000-2019", obs = sum(obs),
              contribution = sum(contribution),
              contribution_from_first_year = sum(contribution_from_first_year)),
  pt %>% summarise(period = "total", obs = sum(obs),
                   contribution = sum(contribution),
                   contribution_from_first_year = sum(contribution_from_first_year)))

cat("\n=== The move split by period ===\n")
cat("contribution: year effects measured from their weighted mean;\n")
cat("contribution_from_first_year: measured from 2000, for comparison.\n")
cat("Only the totals coincide.\n\n")
print(pt %>% mutate(across(starts_with("contribution"), ~ round(.x, 4))) %>%
        as.data.frame(), row.names = FALSE)

# ---- 3. the pre-pandemic window -----------------------------
# a split that needs no reference year: re-estimate without
# 2020-2024 and see how much of the move remains
Fp <- fwl(dat %>% filter(year <= 2019))
cat("\n=== The same move on 2000-2019 alone ===\n")
cat(sprintf("country effects %+.4f, country and year effects %+.4f, move %+.4f\n",
            Fp$b1, Fp$b2, Fp$b1 - Fp$b2))
cat(sprintf("(on 2000-2024 the move is %+.4f)\n", F$b1 - F$b2))

pt <- bind_rows(pt, data.frame(
  period = "re-estimated on 2000-2019", obs = sum(dat$year <= 2019),
  contribution = Fp$b1 - Fp$b2, contribution_from_first_year = NA_real_))

write.csv(yt, file.path(OUT, "36_fwl.csv"), row.names = FALSE)
write.csv(pt, file.path(OUT, "36_fwl_periods.csv"), row.names = FALSE)
cat("\nSaved: 36_fwl.csv, 36_fwl_periods.csv\n")
