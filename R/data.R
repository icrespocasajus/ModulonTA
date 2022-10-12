#' network.TILs Data
#'
#' A dataset containing the transcriptional network inferred using SCENIC pipeline on murine TILs transcriptomics.

#'
#' @format A dataframe encoding the network in 3 columns:'Source','Interaction','Target'
#' @source \url{https://www.nature.com/articles/s41467-021-23324-4}
#' @examples
#' network.TILs
"network.TILs"


#' network.LCMV Data
#'
#' A dataset containing the transcriptional network inferred using SCENIC pipeline on murine LCMV transcriptomics.

#'
#' @format A dataframe encoding the network in 3 columns:'Source','Interaction','Target'
#' @source \url{https://www.nature.com/articles/s41467-021-23324-4}
#' @examples
#' network.LCMV
"network.LCMV"



#' modulons.TILs Data
#'
#' A dataset containing the modulons inferred using modulon analysis on murine TILs transcriptomics.

#'
#' @format A list of 7 elements/modulons with the modulon constituent elements
#' @source \url{https://www.nature.com/articles/s41467-021-23324-4}
#' @examples
#' modulons.TILs
"modulons.TILs"

#' modulons.LCMV Data
#'
#' A dataset containing the modulons inferred using modulon analysis on murine LCMV transcriptomics.

#'
#' @format A list of 12 elements/modulons with the modulon constituent elements
#' @source \url{https://www.nature.com/articles/s41467-021-23324-4}
#' @examples
#' modulons.LCMV
"modulons.LCMV"



#' cc.TILs Data
#'
#' A dataset containing the connected components of each modulon.

#'
#' @format A list with 7 elements (one per modulon) containing all the connected components of the modulon
#' @source \url{https://github.com/icrespocasajus/ModulonCore}
#' @examples
#' cc.TILs
"cc.TILs"

#' cc.LCMV Data
#'
#' A dataset containing the connected components of each modulon.

#'
#' @format A list with 12 elements (one per modulon) containing all the connected components of the modulon
#' @source \url{https://github.com/icrespocasajus/ModulonCore}
#' @examples
#' cc.LCMV
"cc.LCMV"


#' Modulon.Cores.TILs Data
#'
#' A dataset containing the names of the modulon regulatory cores.

#'
#' @format A character vector of 5 elements.
#' @source \url{https://github.com/icrespocasajus/ModulonCore}
#' @examples
#' Modulon.Cores.TILs
"Modulon.Cores.TILs"

#' Modulon.Cores.LCMV Data
#'
#' A dataset containing the names of the modulon regulatory cores.

#'
#' @format A character vector of 11 elements.
#' @source \url{https://github.com/icrespocasajus/ModulonCore}
#' @examples
#' Modulon.Cores.LCMV
"Modulon.Cores.LCMV"


#' DA.TILs Data
#'
#' A dataset containing the OPLS-DA analysis score (WeightStarMN) normalized.

#'
#' @format A list of 5 elements.
#' @examples
#' DA.TILs
"DA.TILs"

#' DA.LCMV Data
#'
#' A dataset containing the OPLS-DA analysis score (WeightStarMN) normalized.

#'
#' @format A list of 6 elements.
#' @examples
#' DA.LCMV
"DA.LCMV"

#' RegAUC.TILs Data
#'
#' Regulon AUC values

#'
#' @format A matrix with 11379 rows (cells) and 583 columns (regulons).
#' @examples
#' RegAUC.TILs
"RegAUC.TILs"

#' RegAUC.LCMV Data
#'
#' Regulon AUC values

#'
#' @format A matrix with 4000 rows (cells) and 646 columns (regulons).
#' @examples
#' RegAUC.LCMV
"RegAUC.LCMV"