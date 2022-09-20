#' @importFrom DelayedArray type
array_type <- function(x) {
    switch(type(x),
        integer="integer",
        double="number",
        logical="boolean",
        character="string",
        "other"
    )
}
