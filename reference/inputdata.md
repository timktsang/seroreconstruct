# Example of input data

This is an example of the input data used in the `seroreconstruct`
function. This data frame illustrates the format of the input data.

## Usage

``` r
data(inputdata)
```

## Format

A data frame with 9 variables, where each row represents an individual:

- age_group:

  0: children, 1: adults, 2: older adults

- start_time:

  start of follow-up

- end_time:

  end of follow-up

- time1:

  date of first serum collection

- time2:

  date of second serum collection

- time3:

  date of third serum collection

- HAI_titer_1:

  HAI titer for first serum collection

- HAI_titer_2:

  HAI titer for second serum collection

- HAI_titer3:

  HAI titer for third serum collection
