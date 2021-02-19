test_type_requirement <- function(values, error, fn, ...) {
  for (value in values) {
    expect_error(
      do.call(fn, list(value, ...)), error)
  }
}
