#' Favorite Colorschemes
#'
#' Provides more colorscheme options in addition to basic ones: 
#' heat.colors, rainbow and so on. Hot and jet colorschemes are 
#' analagous to the MATLAB colorschemes with the same names. 
#' Blue and dusk colorschemes correspond to RColorBrewer "Blues" 
#' and inverse "RdYlBu".
#'
#' @param n Number of different colors in the colorscheme
#' @return Vector of colors as hex RGB strings
#' @name colorschemes
#' 
#' @examples
#' 
#' ## Linear
#' image(as.matrix(1:100), col=blue.colors(100))
#' image(as.matrix(1:100), col=dusk.colors(100))
#' image(as.matrix(1:100), col=hot.colors(100))
#' image(as.matrix(1:100), col=hot2.colors(100))
#' image(as.matrix(1:100), col=jet.colors(100))
#' image(as.matrix(1:100), col=jet2.colors(100))
#'
#' ## Random
#' N <- 10
#' set.seed(0)
#' mat <- matrix(rnorm(N^2), ncol=N)
#' image(mat, col=blue.colors(N^2))
#' image(mat, col=dusk.colors(N^2))
#' image(mat, col=hot.colors(N^2))
#' image(mat, col=hot2.colors(N^2))
#' image(mat, col=jet.colors(N^2))
#' image(mat, col=jet2.colors(N^2))
#' 
#' ## Airy pattern
#' x <- seq(-1,+1,length=100)
#' z <- outer(x,x,function(x,y){sin(10*sqrt(x^2+y^2))})
#' image(z, col=blue.colors(100))
#' image(z, col=dusk.colors(100))
#' image(z, col=hot.colors(100))
#' image(z, col=hot2.colors(100))
#' image(z, col=jet.colors(100))
#' image(z, col=jet2.colors(100))
NULL




# \name{colorschemes}
# \alias{hot.colors}
# \alias{hot2.colors}
# \alias{jet.colors}
# \alias{jet2.colors}
# \alias{blue.colors}
# \alias{dusk.colors}


# \description{
#     Provides more colorscheme options in addition to basic ones: heat.colors,
#     rainbow and so on. Hot and jet colorschemes are analagous to the MATLAB 
#     colorschemes with the same names. Blue and dusk colorschemes correspond to
#     RColorBrewer "Blues" and inverse "RdYlBu".
# }

# \usage{
#     hot.colors(n)
#     hot2.colors(n)
#     jet.colors(n)
#     jet2.colors(n)
#     blue.colors(n)
#     dusk.colors(n)
# }
# 
# \arguments{
#     \item{n}{Number of different colors in the colorscheme}
# }
# 
# \details{
#     The intended use is exactly the same way as basic pallets (e.g. heat.colors).
# }
# 
# \value{
#     Vector of colors in hex RGB strings.
# }

# \references{
#     \url{http://www.nature.com/nmeth/journal/v7/n8/pdf/nmeth0810-573.pdf}
#     \cr
#     \url{http://www.igorexchange.com/node/3817}
# }
# 
# \author{Vlad Petyuk
#         \email{petyuk@gmail.com}
# }
# 
# \seealso{
#     \code{\link{heat.colors}}
#     \code{\link{RColorBrewer}}
#     \code{\link{colorpanel}}
# }
# 
# \examples{
#     ## The use is exactly the same as palettes from grDevices
#     N <- 10
#     set.seed(0)
#     mat <- matrix(rnorm(N^2), ncol=N)
#     image(mat, col=jet.colors(N^2))
#     image(mat, col=jet2.colors(N^2))
#     image(mat, col=hot.colors(N^2))
#     image(mat, col=hot2.colors(N^2))
#     image(mat, col=blue.colors(N^2))
#     image(mat, col=dusk.colors(N^2))
#     
#     ## Example 2. To be Airy pattern
#     x <- seq(-1,+1,length=100)
#     z <- outer(x,x,function(x,y){sin(10*sqrt(x^2+y^2))})
#     image(z, col=jet.colors(100))
#     image(z, col=jet2.colors(100))
#     image(z, col=hot.colors(100))
#     image(z, col=hot2.colors(100))
#     image(z, col=blue.colors(100))
#     image(z, col=dusk.colors(100))
#     
#     ## cycle through the color schemes
#     demo(colorschemes)
# }
# 
# \keyword{color}




# http://www.wavemetrics.com/products/igorpro/creatinggraphs/colortab.htm



#' @export blue.colors
#' @describeIn colorschemes based on RColorBrewer "Blues"
blue.colors <- function(n)
    # matching RColorBrewer "Blues"
{
    colorset <- c("#F7FBFF",
                  "#DEEBF7",
                  "#C6DBEF",
                  "#9ECAE1",
                  "#6BAED6",
                  "#4292C6",
                  "#2171B5",
                  "#08519C",
                  "#08306B")
    grad <- colorRampPalette(colorset, interpolate = "spline")
    return(grad(n))
}


#' @export
#' @describeIn colorschemes based on reverse of RColorBrewer "RdYlBu". 
#'                          Looks like a dusk.
dusk.colors <- function(n)
    # reverse RColorBrewer "RbYlBu"
{
    colorset <- c("#313695",
                  "#4575B4",
                  "#74ADD1",
                  "#ABD9E9",
                  "#E0F3F8",
                  "#FFFFBF",
                  "#FEE090",
                  "#FDAE61",
                  "#F46D43",
                  "#D73027",
                  "#A50026")
    grad <- colorRampPalette(colorset, interpolate = "spline")
    return(grad(n))
}


#' @export
#' @describeIn colorschemes based on Matlab's "hot"
hot.colors <- function(n)
{
    colormat <- matrix(c(041,000,000,
                         082,000,000,
                         132,000,000,
                         173,000,000,
                         214,000,000,
                         255,000,000,
                         255,041,000,
                         255,082,000,
                         255,132,000,
                         255,165,000,
                         255,214,000,
                         255,255,000,
                         255,255,066,
                         255,255,132,
                         255,255,189,
                         255,255,255), ncol=3, byrow=TRUE)
    #    matplot(colormat, type='l', col=c("red","green","blue"),lty=1)
    #    lines(rowMeans(colormat), col='grey', lwd=5)
    colorset <- rgb(colormat, maxColorValue=255)
    grad <- colorRampPalette(colorset, interpolate = "spline")
    return(grad(n))
}


#' @export
#' @describeIn colorschemes similar to Matlab's "hot"
hot2.colors <- function(n){
    colorset <- c("black",
                  "red",
                  "orange",
                  "yellow",
                  "lightyellow")
    grad <- colorRampPalette(colorset, interpolate = "spline")
    return(grad(n))
}


#' @export
#' @describeIn colorschemes based on Matlab's "jet"
jet.colors <- function(n)
{
    colormat <- matrix(c(000,000,189,
                         000,000,255,
                         000,066,255,
                         000,132,255,
                         000,189,255,
                         000,255,255,
                         066,255,189,
                         132,255,132,
                         189,255,066,
                         255,255,000,
                         255,189,000,
                         255,132,000,
                         255,066,000,
                         255,000,000,
                         189,000,000,
                         132,000,000), ncol=3, byrow=TRUE)
    #    matplot(colormat, type='l', col=c("red","green","blue"),lty=1)
    #    lines(rowMeans(colormat), col='grey', lwd=5)
    colorset <- rgb(colormat, maxColorValue=255)
    grad <- colorRampPalette(colorset, interpolate = "spline")
    return(grad(n))
}


#' @export
#' @describeIn colorschemes similar to Matlab's "jet"
jet2.colors <- function(n)
   # kind of Matlab colorscheme
{
   colorset <- c("blue",
                 "#007FFF",
                 "cyan",
                 "#7FFF7F",
                 "yellow",
                 "#FF7F00",
                 "red")
   grad <- colorRampPalette(colorset, interpolate = "spline")
   return(grad(n))
}


# 
# 
# 
# 
# colorscheme("jet") returns function
# colorscheme("jet")(n) returns n colors
# colorscheme("jet")() returns base colors
# # base.colors("jet") returns set of base colors
# # the reason I wanted S4 is to use generics for plot/show
# # otherwise I'll have to come up with special function
# 
# 
# get.colorscheme("jet") returns function
# get.colorscheme("jet")(n) returns n colors
# get.colorscheme("jet")(NULL) returns base colors
# show.colorsheme("jet") # 1) base colors 2) gradient 3) rgb lines
# show.color(colorname or hex) # shows single color as rectangle
# 
# # I want to get and set it as matrix
# get.colorscheme.matrix("jet")
# set.colorscheme.matrix("jet", matrix)
# 
# # I want to get and set it as vector of colors
# get.colorscheme.colors("jet")
# set.colorscheme.colors("jet", matrix) # perhaps "[" generic
# 
# 
# get.colorfunc(jet)
# 
# 
# # what about S4 object?
# obj <- ColorScheme("jet") # sure that that would be two lines
# # well, may be S4 object can be called as function.
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
