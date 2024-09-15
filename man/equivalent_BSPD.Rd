\name{equivalent_BSPD}
\alias{equivalent_BSPD}
\title{Construct Equivalent Estimation Balanced Split Plot Designs}
\description{
The equivalent_BSPD() function constructs equivalent estimation balanced split plot designs for different experimental setups. The function takes inputs related to the number and levels of Whole Plot Factor(s) and Sub Plot Factor(s), the total number of experimental runs, and the bounds for a trend factor. Factor levels can be represented as integers, for example, -1 for low, 0 for medium, and 1 for high. A single whole plot factor or sub plot factor with levels 1 and -1 or 1,0 and -1 should be entered as list(c(1, -1)),list(c(1,0,-1)). If there are two whole plot factors or sub plot factors each with two levels, the input should be entered as list(c(1, -1), c(1, -1)), similarly for more number of factors. The trend factor value lies between 0 to 1. Higher the values, lesser the influence of trend effects on the run order. The function then generates equivalent estimation balanced split plot designs and evaluates them based on D-optimality and Dt-optimality criteria. Additionally, it explores the impact of trend effects and identifies the designs that maximize the trend factor value.
}
\usage{
equivalent_BSPD(Number_of_whole_plot_factors, Whole_Plot_factors_levels,
                 Number_of_sub_plot_factors, Sub_Plot_factors_levels,
                 Total_number_of_runs_required, lower_bound, upper_bound)
}
\arguments{
\item{Number_of_whole_plot_factors}{Number of Whole Plot Factor(s)}
\item{Whole_Plot_factors_levels}{Levels of Whole Plot Factor(s)}
\item{Number_of_sub_plot_factors}{Number of Sub Plot Factor(s)}
\item{Sub_Plot_factors_levels}{Levels of Sub Plot Factor(s)}
\item{Total_number_of_runs_required}{Total Number of Runs Required}
\item{lower_bound}{Lower bound of Trend Factor Value}
\item{upper_bound}{Upper bound of Trend Factor Value}
}
\value{
Returns a list with the following components:
\item{Total_equivalent_estimation_designs}{Total number of equivalent estimation balanced split plot designs for the given inputs.}
\item{Equivalent_Estimation_Designs}{All the equivalent estimation balanced split plot designs for the given inputs.}
\item{All_equivalent_estimation_designs_with_D_Dt_Trend_Factor}{All the equivalent estimation balanced split plot designs with D, Dt, and Trend Factor value for the given inputs.}
\item{Max_D_value}{Maximum D-value within the generated equivalent estimation balanced split plot designs.}
\item{D_optimal_designs}{Designs with the Maximum D-value within the generated equivalent estimation balanced split plot designs.}
\item{Max_Dt_value}{Maximum Dt-value within the generated equivalent estimation balanced split plot designs.}
\item{Dt_optimal_designs}{Designs with the Maximum Dt-value within the generated equivalent estimation balanced split plot designs.}
\item{Max_Trend_factor_value}{Maximum Trend Factor Value for the generated equivalent estimation balanced split plot designs.}
\item{Number_of_Designs_Max_Trend_Factor}{Number of equivalent estimation balanced split plot designs with Maximum Trend Factor value.}
\item{Equivalent_Estimation_designs_in_trend_factor_range}{Equivalent estimation balanced split plot designs within the specified range of trend factor.}
}
\references{
Bijoy Chanda, Arpan Bhowmik, Cini Varghese, Seema Jaggi, Anindita Datta, Baidya Nath Mandal and Soumen Pal (2024). Incomplete equivalent estimation split plot designs. Journal of Community Mobilization and Sustainable Development, 19(2), 357-360.

Harrison Macharia and Peter Goos (2010). D-optimal and D-efficient equivalent-estimation second order split-plot designs. Journal of Quality Technology, 42(4), 358-372.

Lieven Tack and Martina Vandebroek (2001). (Dt,C)-optimal run orders. Journal of Statistical Planning and Inference, 98, 293-310.
}
\examples{
# One whole plot factor at two levels and one sub plot factor at three levels in four runs
Result1 <- equivalent_BSPD(1, list(c(1, -1)), 1, list(c(1,0,-1)), 4, 0.1, 0.9)

# One whole plot factor at three levels and one sub plot factor at three levels in six runs
Result2 <- equivalent_BSPD(1, list(c(1, 0, -1)), 1, list(c(1, 0, -1)), 6, 0.5, 0.95)

# One whole plot factor at three levels and two sub plot factors each at two levels in six runs
Result3 <- equivalent_BSPD(1, list(c(1, 0, -1)), 2, list(c(1, -1), c(1, -1)), 6, 0.6, 0.9)

# Accessing results
Result1$Total_equivalent_estimation_designs
Result1$Equivalent_Estimation_Designs
Result1$All_equivalent_estimation_designs_with_D_Dt_Trend_Factor
Result1$Max_D_value
Result1$D_optimal_designs
Result1$Max_Dt_value
Result1$Dt_optimal_designs
Result1$Max_Trend_factor_value
Result1$Number_of_Designs_Max_Trend_Factor
Result1$Equivalent_Estimation_designs_in_trend_factor_range
}

