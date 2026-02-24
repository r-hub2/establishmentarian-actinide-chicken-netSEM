# Setup file for testthat
# This file is run before tests to set up shared resources

# Ensure datasets are available
data(acrylic, package = "netSEM")
data(PVmodule, package = "netSEM")

# Suppress warnings during tests for cleaner output
options(warn = -1)
