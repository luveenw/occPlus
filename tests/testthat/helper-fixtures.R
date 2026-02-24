# Auto-sourced by testthat before every test file.
# Loads both package datasets into the test environment so every test file
# can reference `sampledata` and `fitmodel` without repeating the load calls.

data("sampledata",   package = "occPlus")   # -> object named `data` (sampledata)
data("sampleresults", package = "occPlus")  # -> object named `fitmodel`

# sampledata is stored under the variable name `data` inside the .rda — alias it
# to avoid masking base::data() and to make test code self-documenting.
sampledata <- data   # nolint: object_name
