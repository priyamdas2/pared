if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("The Rcpp package is required.")
}

if (!exists("anti_transformation", mode = "function")) {
  Rcpp::cppFunction('
    Rcpp::NumericVector anti_transformation(Rcpp::NumericVector x,
                                             Rcpp::NumericVector lb,
                                             Rcpp::NumericVector ub) {
      int n = x.size();
      Rcpp::NumericVector transformed_x(n);
      for(int i = 0; i < n; ++i) {
        transformed_x[i] = (x[i] - lb[i]) / (ub[i] - lb[i]);
      }
      return transformed_x;
    }
  ')
}

if (!exists("transformation", mode = "function")) {
  Rcpp::cppFunction('
    Rcpp::NumericVector transformation(Rcpp::NumericVector x,
                                        Rcpp::NumericVector lb,
                                        Rcpp::NumericVector ub) {
      int n = x.size();
      Rcpp::NumericVector transformed_x(n);
      for(int i = 0; i < n; ++i) {
        transformed_x[i] = x[i] * (ub[i] - lb[i]) + lb[i];
      }
      return transformed_x;
    }
  ')
}

if (!exists("update_x", mode = "function")) {
  Rcpp::cppFunction('
    Rcpp::NumericVector update_x(Rcpp::NumericVector x,
                                  double epsilon,
                                  double rho,
                                  double phi) {
      int n = x.size();
      int j = 0;
      Rcpp::NumericVector x_updated(2 * n);
      double current_value = 0;
      double updated_value = 0;
      double epsilon_temp = 0;
      double ff = 0;
      int f = 0;
      
      for(int i = 0; i < (2 * n); ++i) {
        j = std::ceil((i + 1) / 2.0) - 1;
        epsilon_temp = std::pow(-1.0, i) * epsilon;
        current_value = x[j];
        updated_value = current_value + epsilon_temp;
        
        if(updated_value > 1 && current_value < 1 - phi) {
          ff = std::log(epsilon_temp / (1 - current_value)) / std::log(rho);
          f = std::ceil(ff);
          epsilon_temp = epsilon_temp / std::pow(rho, f);
          updated_value = current_value + epsilon_temp;
        } else {
          if(updated_value < 0 && current_value > phi) {
            ff = std::log(-epsilon_temp / current_value) / std::log(rho);
            f = std::ceil(ff);
            epsilon_temp = epsilon_temp / std::pow(rho, f);
            updated_value = current_value + epsilon_temp;
          }
        }
        
        if(updated_value > 1 || updated_value < 0) {
          updated_value = current_value;
        }
        
        x_updated[i] = updated_value;
      }
      
      return x_updated;
    }
  ')
}

RMPSH_opt <- function(x0,
                      func,
                      lb,
                      ub,
                      rho_1 = 2,
                      rho_2 = 2,
                      phi = 10^(-6),
                      no_runs = 1000,
                      max_iter = 10000,
                      s_init = 2,
                      tol_fun = 10^(-6),
                      tol_fun_2 = 10^(-20),
                      max_time = 36000,
                      print_output = 0) {
  
  M <- length(x0)
  
  if (!is.function(func)) {
    stop("func must be a function.")
  }
  
  if (length(lb) != M || length(ub) != M) {
    stop("x0, lb, and ub must have the same length.")
  }
  
  if (any(!is.finite(x0)) || any(!is.finite(lb)) || any(!is.finite(ub))) {
    stop("x0, lb, and ub must contain only finite values.")
  }
  
  if (any(ub <= lb)) {
    stop("Each upper bound must be greater than the corresponding lower bound.")
  }
  
  if (any(x0 < lb) || any(x0 > ub)) {
    stop("Starting point x0 is outside the domain.")
  }
  
  start_value <- func(x0)
  start_time <- Sys.time()
  
  theta_array <- matrix(0, no_runs, M)
  each_run_solution <- array(0, no_runs)
  
  ill_condition <- 0
  theta <- anti_transformation(x0, lb, ub)
  
  for (iii in seq_len(no_runs)) {
    
    epsilon <- s_init
    
    if (iii == 1) {
      rho <- rho_1
      theta <- anti_transformation(x0, lb, ub)
      
      if (any(theta < 0) || any(theta > 1)) {
        stop("Transformed starting point is outside [0, 1]^M.")
      }
      
    } else {
      theta <- as.numeric(theta_array[iii - 1, ])
      rho <- rho_2
    }
    
    if (ill_condition == 1) {
      break
    }
    
    array_of_values <- array(0, max_iter)
    
    for (i in seq_len(max_iter)) {
      
      current_lh <- func(transformation(theta, lb, ub))
      
      possible_x_coords <- update_x(theta, epsilon, rho, phi)
      
      total_lh <- rep(1, 2 * M)
      
      for (kk in seq_len(2 * M)) {
        
        candidate_theta <- theta
        coord_number <- ceiling(kk / 2)
        
        if (possible_x_coords[kk] == theta[coord_number]) {
          total_lh[kk] <- current_lh
        } else {
          candidate_theta[coord_number] <- possible_x_coords[kk]
          total_lh[kk] <- func(transformation(candidate_theta, lb, ub))
        }
      }
      
      new_min <- min(total_lh)
      
      if (new_min < current_lh) {
        position_new_min <- which.min(total_lh)
        pos_of_theta <- ceiling(position_new_min / 2)
        theta[pos_of_theta] <- possible_x_coords[position_new_min]
      }
      
      array_of_values[i] <- min(new_min, current_lh)
      
      if (print_output == 1) {
        cat(
          "\n",
          "Run no. ",
          iii,
          ", iteration no. ",
          i,
          ", current fun value = ",
          array_of_values[i]
        )
      }
      
      if (i > 1) {
        if (abs(array_of_values[i] - array_of_values[i - 1]) < tol_fun) {
          if (abs(epsilon) > phi) {
            epsilon <- epsilon / rho
          } else {
            break
          }
        }
      }
      
      now_time <- Sys.time()
      now_time_spent <- as.numeric(difftime(now_time, start_time, units = "secs"))
      
      if (now_time_spent > max_time) {
        if (print_output == 1) {
          cat("\n")
          print("Time's up")
          cat("\n", "Starting objective function value: ", start_value)
          cat("\n", "Final objective function value: ", current_lh)
        }
        ill_condition <- 1
        break
      }
    }
    
    theta_array[iii, ] <- theta
    each_run_solution[iii] <- func(transformation(theta, lb, ub))
    
    if (iii > 1) {
      old_soln <- theta_array[iii - 1, ]
      new_soln <- theta_array[iii, ]
      
      if (norm(as.matrix(new_soln - old_soln), type = "F") < tol_fun_2) {
        break
      }
    }
  }
  
  end_time <- Sys.time()
  time_spent <- as.numeric(difftime(end_time, start_time, units = "secs"))
  final_point <- transformation(theta, lb, ub)
  final_value <- func(final_point)
  
  if (print_output == 1) {
    cat("\n", "Starting objective function value: ", start_value)
    cat("\n", "Final objective function value: ", final_value)
    cat("\n", "Total time required (in secs): ", time_spent)
    cat("\n", "Obtained minima point is: ", final_point)
    cat("\n")
  }
  
  return(final_point)
}