# TOSDIF

Convenient function for the test of small differences in fit as proposed by MacCallum et al. (2006).

This function depends on "lavaan" and "semTools", both of which you'll probably already have loaded at this point.

The function only needs two arguments to be specified: `fit1` and `fit2`, which are lavaan fit-objects.
`groups` is set to 1 by default and can be changed if your fits have more than one group.
`alpha`, `rmseaa` and `rmseab` also have default values. I recommend against changing these, but you can, if you wish to.
`robust` is also by default set to be `TRUE`, but will only work if the fitted lavaan objects contain the information for robust variants of chi-squared. This can be achieved by setting `estimator = 'MLR'` when fitting the object.



Feel free to contact me for questions or suggestions.