
equivalent_BSPD<-function(Number_of_whole_plot_factors,Whole_Plot_factors_levels,

                          Number_of_sub_plot_factors,Sub_Plot_factors_levels,

                          Total_number_of_runs_required,lower_bound, upper_bound){

  generate_combinations <- function(vectors, current_combination = NULL, index = 1) {

    if (index > length(vectors)) {

      return(matrix(unlist(current_combination), ncol = length(vectors), byrow = TRUE))

    } else {

      current_vector <- vectors[[index]]

      new_combinations <- lapply(current_vector, function(x) {

        new_current_combination <- append(current_combination, list(x))

        generate_combinations(vectors, new_current_combination, index + 1)

      })
      return(do.call(rbind, new_combinations))
    }
  }


  Whole_plot_treatment_Combinations<-generate_combinations(Whole_Plot_factors_levels)

  Sub_plot_treatment_combinations<-generate_combinations(Sub_Plot_factors_levels)

  Number_of_sub_plots_per_whole_plot<-Total_number_of_runs_required/nrow(Whole_plot_treatment_Combinations)

  First_matrix<-matrix(rep(Whole_plot_treatment_Combinations, each = nrow(Sub_plot_treatment_combinations)),ncol = ncol(Whole_plot_treatment_Combinations))


  Second_matrix<-matrix(rep(t(Sub_plot_treatment_combinations), nrow(Whole_plot_treatment_Combinations)),ncol = ncol(Sub_plot_treatment_combinations) ,byrow = TRUE)

  Initial_design<-cbind(First_matrix,Second_matrix)



  generate_permutations <- function(numbers, subset_size) {

    perms <- vector("list")

    permute <- function(current_permutation, remaining_numbers) {

      if (length(current_permutation) == subset_size) {

        perms <<- append(perms, list(current_permutation))

      } else {

        for (i in seq_along(remaining_numbers)) {

          next_permutation <- c(current_permutation, remaining_numbers[i])

          next_remaining_numbers <- remaining_numbers[-i]

          permute(next_permutation, next_remaining_numbers)
        }
      }
    }

    permute(c(), numbers)

    return(perms)
  }

  permutation<-generate_permutations(1:nrow(Sub_plot_treatment_combinations),Number_of_sub_plots_per_whole_plot)

  First_permutation<-matrix(unlist(permutation),ncol=Number_of_sub_plots_per_whole_plot,byrow = TRUE)

  Second_permutation<-generate_permutations(1:nrow(Whole_plot_treatment_Combinations),nrow(Whole_plot_treatment_Combinations))


  All_set<-list()

  for (i in 1:nrow(Whole_plot_treatment_Combinations)) {

    set<-list()


    for (j in 1:nrow(First_permutation)) {


      Sub_plots_per_whole_plot<-Initial_design[((i-1)*nrow(Sub_plot_treatment_combinations)+1):(i*nrow(Sub_plot_treatment_combinations)),][c(First_permutation[j,]),]

      set<-append(set,list(Sub_plots_per_whole_plot))


    }
    All_set<-append(All_set,list(set) )
  }



  combinations<-expand.grid(All_set)

  all_combinations<-apply(combinations,1,function(row) as.list(row))


  Third_permutation<-matrix(unlist(Second_permutation),ncol=nrow(Whole_plot_treatment_Combinations),byrow = TRUE)


  Z<-matrix(0,nrow=Total_number_of_runs_required,ncol =nrow(Whole_plot_treatment_Combinations) )

  for (i in 1:nrow(Whole_plot_treatment_Combinations)) {

    for (j in 1:Number_of_sub_plots_per_whole_plot) {

      k=j+(i-1)*Number_of_sub_plots_per_whole_plot

      Z[k,i]=1

    }

  }


  All_equivalent_estimation_designs_with_D_Dt_Trend_Factor<- list()
  Equivalent_Estimation_designs_in_trend_factor_range <- list()
  D_values <- c()

  Dt_values <- c()

  Max_D_value <- -Inf

  Max_Dt_value <- -Inf

  D_optimal_designs <- list()

  Dt_optimal_designs <- list()

  Total_equivalent_estimation_designs <- 0

  Equivalent_Estimation_Designs<-list()

  Max_Trend_factor_value <- -Inf

  Number_of_Designs_Max_Trend_Factor <- 0



  for (i in 1:length(all_combinations)) {

    for (j in 1:nrow(Third_permutation)) {

      permutation_indices <- as.numeric(Third_permutation[j, ])

      permutation <- all_combinations[[i]][permutation_indices]


      X <- do.call(rbind, permutation)

      X_with_mean<-cbind(matrix(1,nrow = Total_number_of_runs_required,ncol = 1),X)


      X_t <- t(X_with_mean)

      Z_t <- t(Z)

      D_1<- Z %*% Z_t

      B_1<- D_1 %*% X_with_mean

      Q_1<- X_t %*% B_1

      P_1<- X_t %*% X_with_mean

      if (det(P_1) != 0) {

        P_inv_1 <-solve(P_1)

        K_1 <- P_inv_1 %*% Q_1

        S_1<- X_with_mean %*% K_1

        if (identical(S_1, B_1)) {

          Total_equivalent_estimation_designs <- Total_equivalent_estimation_designs + 1

          Equivalent_Estimation_Designs<-append(Equivalent_Estimation_Designs,list(X))


          Design_mat<-X_with_mean

          Design_mat_t <- t(Design_mat)

          I_mat<- diag(Total_number_of_runs_required)

          V<- I_mat +D_1

          V_inv<- solve(V)

          A1 <- Design_mat_t%*%V_inv

          A2<- A1 %*% Design_mat


          D_value <- round(det(A2),2)

          D_values <- c(D_values, D_value)

          G_matrix <- poly(1:Total_number_of_runs_required,degree = 1,simple = TRUE)

          A3<-rbind(Design_mat_t,t(G_matrix))

          A4 <- cbind(Design_mat, G_matrix)

          A5 <- A3%*% V_inv

          A6<-A5%*%A4

          A7 <- t(G_matrix) %*% V_inv

          A8<-A7%*%G_matrix

          Final_value <- det(A6)/det(A8)

          Dt_value <- round(Final_value,2)

          threshold <- 1e-4


          if (Dt_value < threshold) {

            Dt_value <- 0
          }

          Dt_values <- c(Dt_values, Dt_value)

          Trend_factor <- round((Dt_value / D_value)^(1/(1+ncol(Initial_design))),2)



          All_equivalent_estimation_designs_with_D_Dt_Trend_Factor <- append(All_equivalent_estimation_designs_with_D_Dt_Trend_Factor,
                                                                             list(list(X,D_value=D_value,Dt_value=Dt_value,Trend_factor=Trend_factor)))




          if (D_value > Max_D_value) {



            Max_D_value <- D_value

            D_optimal_designs <- list(list(X, D_value = D_value, Dt_value = Dt_value, Trend_factor = Trend_factor))

          } else if (D_value == Max_D_value) {



            D_optimal_designs <- append(D_optimal_designs, list(list(X, D_value = D_value, Dt_value = Dt_value, Trend_factor = Trend_factor)))
          }



          if (Dt_value > Max_Dt_value) {


            Max_Dt_value <- Dt_value



            Dt_optimal_designs <- list(list(X, D_value = D_value, Dt_value = Dt_value, Trend_factor = Trend_factor))

          } else if (Dt_value == Max_Dt_value) {



            Dt_optimal_designs <- append(Dt_optimal_designs, list(list(X, D_value = D_value, Dt_value = Dt_value, Trend_factor = Trend_factor)))
          }



          if (!is.nan(Trend_factor) && Trend_factor > Max_Trend_factor_value) {



            Max_Trend_factor_value <- Trend_factor

            Number_of_Designs_Max_Trend_Factor <- 1

          } else if (!is.nan(Trend_factor) && Trend_factor == Max_Trend_factor_value) {



            Number_of_Designs_Max_Trend_Factor <- Number_of_Designs_Max_Trend_Factor + 1
          }


          if (!is.nan(Trend_factor) && Trend_factor >= lower_bound && Trend_factor <= upper_bound) { # trend factor should not be NAN



            design_with_Trend_factor <- list(X,D_value=D_value,

                                             Dt_value=Dt_value, Trend_factor = Trend_factor)

            Equivalent_Estimation_designs_in_trend_factor_range <- append(Equivalent_Estimation_designs_in_trend_factor_range, list(design_with_Trend_factor))
          }

        }


      }


    }

  }


  if (length(Equivalent_Estimation_designs_in_trend_factor_range) == 0) {

    Equivalent_Estimation_designs_in_trend_factor_range <- "No designs found within the specified range of trend factor"
  }


  if (length(unique(D_values)) == 1) {
    D_optimal_designs <- "No D optimal designs found."
  }


  if (length(unique(Dt_values)) == 1) {
    Dt_optimal_designs <- "No Dt optimal designs found."
  }

  return(list(Max_D_value = Max_D_value,
              D_optimal_designs=D_optimal_designs,
              Max_Dt_value = Max_Dt_value,
              Dt_optimal_designs = Dt_optimal_designs,
              Max_Trend_factor_value=Max_Trend_factor_value,
              Number_of_Designs_Max_Trend_Factor=Number_of_Designs_Max_Trend_Factor,
              Equivalent_Estimation_Designs=Equivalent_Estimation_Designs,
              Total_equivalent_estimation_designs=Total_equivalent_estimation_designs,
              All_equivalent_estimation_designs_with_D_Dt_Trend_Factor=All_equivalent_estimation_designs_with_D_Dt_Trend_Factor,
              Equivalent_Estimation_designs_in_trend_factor_range = Equivalent_Estimation_designs_in_trend_factor_range))
}

