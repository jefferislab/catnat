.onLoad <- function (libname, pkgname) {

  # Note that the Morpho affine registrations need to be inverted
  # see e.g. elmr::tpsreg for discussion of source/reference conventions
  fcwb_chiangm <-
    reglist(
      solve(readRDS(system.file(
        "extdata/CMTKreg/InitialAffine/initialiseCMTKreg_ChiangMaleTowardsFCWB.rds",
        package = 'catnat'
      ))),
      cmtkreg(system.file(
        "extdata/CMTKreg/Registration/warp/FCWB_typicalbrainmale_01_warp_m0g80c8e1e-1x26r4.list/",
        package = 'catnat'
      )),
      solve(readRDS(system.file(
        "extdata/CMTKreg/InitialAffine/finalaffine_ChiangMaleTowardsFCWB.rds",
        package = 'catnat'
      )))
    )

  fcwb_chiangf <-
    reglist(
      solve(readRDS(system.file(
        "extdata/CMTKreg/InitialAffine/initialiseCMTKreg_ChiangFemaleTowardsFCWB.rds",
        package = 'catnat'
      ))),
      cmtkreg(system.file(
        "extdata/CMTKreg/Registration/warp/FCWB_typicalbrainfemale_01_warp_m0g80c8e1e-1x26r4.list/",
        package = 'catnat'
      )),
      solve(readRDS(system.file(
        "extdata/CMTKreg/InitialAffine/finalaffine_ChiangFemaleTowardsFCWB.rds",
        package = 'catnat'
      )))
    )

  fcwb_chiangm2=fcwb_chiangm[-3]
  fcwb_chiangf2=fcwb_chiangf[-3]

  nat.templatebrains::add_reglist(fcwb_chiangf, reference = nat.flybrains::FCWB, sample = 'chiangf')
  nat.templatebrains::add_reglist(fcwb_chiangm, reference = nat.flybrains::FCWB, sample = 'chiangm')

  nat.templatebrains::add_reglist(fcwb_chiangf2, reference = nat.flybrains::FCWB, sample = 'chiangf2')
  nat.templatebrains::add_reglist(fcwb_chiangm2, reference = nat.flybrains::FCWB, sample = 'chiangm2')
}
