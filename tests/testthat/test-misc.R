# Test funcs. in misc.R

test_that("simulateDataset gives expected output", {
    expect_equal(nrow(simulateDataset(n_genes = 30)), 30)
    expect_equal(nlevels(simulateDataset(n_rings = 4)$cluster), 4)
})
