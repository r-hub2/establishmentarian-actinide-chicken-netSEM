# netSEM 0.7.0

-   Refactored. Changes to codebase, but with same functionality.

-   netSEMp1() and netSEMp2() function call updates

    -   User no longer needs to follow a specific column order in data frame, but can explicitly declare in function argument.

-   Updates to output format.

    -   netSEMp1() "best relationship"" outputs are now in a data frame format with equation forms included.
    -   pseudoR2 calculations included for comparison with nls() model evaluation.
    -   netSEMp1() now saves model object only for best relationship.

-   paths() function now return all degradation pathways with varying lengths.

    -   Rank-ordering is evaluated by looking at the relationship with lowest adjusted R-squared.

# netSEM 0.6.2

-   Modified package level documentation to be consistent with new method as per request.

-   Updated author list: Added new contributors.

# netSEM 0.6.1

## Changes

-   Changed acknowledgements part in Description. No other changes to function scripts.

# netSEM 0.6.0

## New Features

-   Added principle 2 which functions on multiple regression. User can choose model evaluation criteria between AIC or BIC.

-   Changed Stressor & Response color scheme in netSEM network diagrams
