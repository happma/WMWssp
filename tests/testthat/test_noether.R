context("Noether Formula")

# data
x <- c(315,375,356,374,412,418,445,403,431,410,391,475,379)
y <- x - 20

# calculate sample size, true result
result_noetherN <- 114.74467

test_that("function Noether", {
  expect_equivalent(WMWssp::WMWssp_noether(alpha = 0.05, power = 0.8, t =1/2, p = 0.349)$result[4, ], result_noetherN, tolerance=1e-4)
})
