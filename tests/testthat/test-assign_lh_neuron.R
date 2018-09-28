context("Assigning LH cell types")
test_that("Let's assign some neurons or known types",{
  truelhns = subset(lhns::most.lhns, anatomy.group!="notLHproper")
  tneurons = truelhns[1:10]
  aneurons = assign_lh_neuron(tneurons)
  expect_is(aneurons,"neuronlist")
  expect_identical(tneurons[,"cell.type"],as.character(aneurons[,"cell.type.1"]))
})
