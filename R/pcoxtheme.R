#' Set theme for pcoxtime plots
#' 
#' Sets a theme for pcoxtime and other ggplot objects
#'
#' @return No return value, called for side effects (setting pcotime plotting theme).
#'
#' @examples
#' library(ggplot2)
#' pcoxtheme()
#' data(heart, package="survival")
#' lam <- 0.02
#' alp <- 1
#' pfit <- pcoxtime(Surv(start, stop, event) ~ age + year + surgery + transplant
#' 	, data = heart
#' 	, lambda = lam
#'		, alpha = alp
#'	)
#'
#' # Plot survival curves
#' psurv <- pcoxsurvfit(pfit)
#' plot(psurv)
#' @import ggplot2
#' @export

pcoxtheme <- function(){
   theme_set(theme_bw() +
      theme(panel.spacing = grid::unit(0,"lines")
      	, plot.title = element_text(hjust = 0.5)
			, legend.position = "bottom"
			, axis.ticks.y = element_blank()
			, axis.text.x = element_text(size = 12)
			, axis.text.y = element_text(size = 12)
			, axis.title.x = element_text(size = 12)
			, axis.title.y = element_text(size = 12)
			, legend.title = element_text(size = 13, hjust = 0.5)
			, legend.text = element_text(size = 13)
			, panel.grid.major = element_blank()
			, legend.key.size = unit(0.8, "cm")
			, legend.key = element_rect(fill = "white")
			, panel.spacing.y = unit(0.3, "lines")
			, panel.spacing.x = unit(1, "lines")
			, strip.background = element_blank()
			, panel.border = element_rect(colour = "grey"
				, fill = NA
				, size = 0.8
			)
			, strip.text.x = element_text(size = 11
				, colour = "black"
				, face = "bold"
			)
      )
   )
}


