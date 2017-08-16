context("googlesheets")

test_that("construct connector URL", {
  testinput = structure(
    list(
      x = c(434086L, 439796L),
      y = c(235265L, 238586L),
      z = c(167930L, 153790L),
      partner_skid = c(39254L, 1316642L)
    ),
    .Names = c("x", "y", "z", "partner_skid"),
    row.names = c(NA,-2L),
    class = c("tbl_df", "tbl", "data.frame")
  )

  testouput=c("https://neuropil.janelia.org/tracing/fafb/v13?pid=1&zp=167930&yp=235265&xp=434086&tool=tracingtool&active_skeleton_id=39254&sid0=5&s0=0",
    "https://neuropil.janelia.org/tracing/fafb/v13?pid=1&zp=153790&yp=238586&xp=439796&tool=tracingtool&active_skeleton_id=1316642&sid0=5&s0=0")
  expect_equal(connector_URL(testinput), testouput)

  expect_error(apply(testinput, 1, connector_URL))
})
