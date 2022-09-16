test_that("cvine mixture works", {
  G <- diag(1, 4)
  G[upper.tri(G)] <- 0.6
  G[lower.tri(G)] <- 0.6
  xmod <- mpower::MixtureModel(method = "cvine", G = G, m = 50,
      cvine_marginals = list(binary = "qbinom(size=1, prob=0.2)",
                             count = "qpois(lambda = 1)",
                             ordinal = "qmultinom(probs=c(0.6, 0.2, 0.2))",
                             categorical = "qmultinom(probs=c(0.6, 0.2, 0.2))"),
      cvine_dtypes = list(categorical = "factor"))
  X <- genx(xmod, 20)
  expect_equal(dim(X), c(20, 4))
  expect_equal(sapply(X, class), c("binary"="numeric", "count"="numeric",
  "ordinal"="numeric", "categorical"="factor"))
  expect_equal(any(class(mplot(xmod)$hist[[1]]) %in% c("gg", "ggplot")), TRUE)
})

test_that("resampling mixture works", {
  data("NHANES")
  nhanes_demo <- NHANES %>%
    select(Education, HHIncome, HomeOwn, Poverty, Age) %>%
    filter(complete.cases(.)) %>%
    head(100)
  xmod <- mpower::MixtureModel(method = "resampling", data = nhanes_demo)
  X <- genx(xmod, 20)
  expect_equal(dim(X), c(20, 5))
  expect_equal(any(class(mplot(xmod)$hist[[1]]) %in% c("gg", "ggplot")), TRUE)
})

test_that("Gaussian copula mixture works", {
  data("NHANES")
  nhanes_demo <- NHANES %>%
    select(BMI, Education, HHIncome, HomeOwn, Poverty, Age) %>%
    filter(complete.cases(.)) %>%
    mutate(Education = case_when(Education == "8th Grade" ~ 1,
                                 Education == "9 - 11th Grade" ~ 2,
                                 Education == "High School" ~ 3,
                                 Education == "Some College" ~ 4,
                                 Education == "College Grad" ~ 5),
           HHIncome = case_when(HHIncome == "more 99999" ~ 12,
                                HHIncome == "75000-99999" ~ 11,
                                HHIncome == "65000-74999" ~ 10,
                                HHIncome == "55000-64999" ~ 9,
                                HHIncome == "45000-54999" ~ 8,
                                HHIncome == "35000-44999" ~ 7,
                                HHIncome == "25000-34999" ~ 6,
                                HHIncome == "20000-24999" ~ 5,
                                HHIncome == "15000-19999" ~ 4,
                                HHIncome == "10000-14999" ~ 3,
                                HHIncome == " 5000-9999" ~ 2,
                                HHIncome == " 0-4999" ~ 1),
            HomeOwn = ifelse(HomeOwn == "Own", 1, 0)) %>%
    head(100)
  xmod <- mpower::MixtureModel(method = "estimation" , data = nhanes_demo,
        sbg_args = list(nsamp = 100))
  X <- genx(xmod, 20)
  expect_equal(dim(X), c(20, 6))
  expect_equal(all(sapply(X, class) == "numeric"), TRUE)
  expect_equal(any(class(mplot(xmod)$hist[[1]]) %in% c("gg", "ggplot")), TRUE)
  expect_equal(any(class(mplot(xmod)$corr) %in% c("gg", "ggplot")), TRUE)
})
