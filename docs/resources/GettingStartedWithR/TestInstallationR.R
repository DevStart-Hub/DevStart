# Complete R Environment Test Script
# Tests all essential packages for DevStart tutorials

cat("🔍 Testing R environment setup for DevStart...\n\n")

# List of packages to test
packages_to_test <- c("tidyverse", "easystats", "lme4", "lmerTest", 
                      "patchwork", "PupillometryR", "fitdistrplus")

# Initialize results tracking
test_results <- list()
all_passed <- TRUE

cat("📦 Testing package loading...\n")

# Test each package individually
for (pkg in packages_to_test) {
  cat(sprintf("Testing %s... ", pkg))
  
  result <- tryCatch({
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    cat("✅\n")
    TRUE
  }, error = function(e) {
    cat("❌\n")
    cat(sprintf("   Error: %s\n", conditionMessage(e)))
    FALSE
  })
  
  test_results[[pkg]] <- result
  if (!result) all_passed <- FALSE
}

cat("\n", rep("=", 50), "\n")

# If all packages loaded, run functionality tests
if (all_passed) {
  cat("🎉 ALL PACKAGES LOADED SUCCESSFULLY!\n\n")
  cat("🧪 Running functionality tests...\n\n")
  
  tryCatch({
    
    # Test 1: tidyverse functionality
    cat("📊 Test 1: tidyverse data manipulation...\n")
    library(dplyr)
    library(ggplot2)
    
    iris_summary <- iris %>%
      group_by(Species) %>%
      summarise(mean_sepal = mean(Sepal.Length), .groups = 'drop')
    
    basic_plot <- ggplot(iris, aes(x = Species, y = Sepal.Length)) +
      geom_boxplot() +
      theme_minimal()
    
    cat("   ✅ tidyverse working perfectly!\n\n")
    
    # Test 2: Mixed-effects model with lme4 and lmerTest
    cat("📈 Test 2: Mixed-effects modeling (lme4 + lmerTest)...\n")
    
    # Use built-in sleepstudy data
    data(sleepstudy, package = "lme4")
    model <- lmer(Reaction ~ Days + (Days|Subject), data = sleepstudy)
    
    # Test lmerTest functionality
    library(lmerTest)
    model_test <- lmer(Reaction ~ Days + (Days|Subject), data = sleepstudy)
    model_summary <- summary(model_test)  # Should include p-values
    
    cat("   ✅ lme4 and lmerTest working perfectly!\n\n")
    
    # Test 3: easystats functionality
    cat("🎯 Test 3: easystats model analysis...\n")
    library(modelbased)
    library(parameters)
    
    estimates <- estimate_means(model, by = "Days")
    model_params <- model_parameters(model)
    
    cat("   ✅ easystats working perfectly!\n\n")
    
    # Test 4: patchwork for plot combination
    cat("🎨 Test 4: patchwork plot combination...\n")
    library(patchwork)
    
    plot1 <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width)) + geom_point()
    plot2 <- ggplot(iris, aes(x = Petal.Length, y = Petal.Width)) + geom_point()
    combined_plot <- plot1 + plot2
    
    cat("   ✅ patchwork working perfectly!\n\n")
    
    # Test 5: fitdistrplus
    cat("📊 Test 5: fitdistrplus distribution fitting...\n")
    library(fitdistrplus)
    
    # Fit normal distribution to some data
    sample_data <- rnorm(100)
    fit_normal <- fitdist(sample_data, "norm")
    
    cat("   ✅ fitdistrplus working perfectly!\n\n")
    
    # Test 6: PupillometryR (basic load test)
    cat("👁️  Test 6: PupillometryR loading...\n")
    library(PupillometryR)
    
    cat("✅ PupillometryR loaded successfully!\n")
    cat("   Note: Full PupillometryR testing requires eye-tracking data\n\n")
    
    # Display some results
    cat("📋 Sample Results:\n")
    cat("=" %>% rep(40) %>% paste(collapse = ""), "\n")
    print(iris_summary)
    cat("\n")
    print(estimates[1:3, ])  # Show first few estimates
    
    # Create final celebration plot
    cat("\n🎨 Creating celebration plot...\n")
    
    # Create a beautiful success plot
    success_plot <- ggplot() +
      # Background
      geom_rect(aes(xmin = -1, xmax = 1, ymin = -1, ymax = 1), 
                fill = "lightgreen", alpha = 0.3) +
      # Success message
      annotate("text", x = 0, y = 0, 
               label = "🎉 ALL WENT WELL! 🎉\n\nYour R Environment\nis Ready for DevStart!", 
               size = 8, hjust = 0.5, vjust = 0.5, 
               fontface = "bold", color = "darkgreen") +
      # Decorative elements
      annotate("text", x = -0.8, y = 0.8, label = "✅", size = 10) +
      annotate("text", x = 0.8, y = 0.8, label = "✅", size = 10) +
      annotate("text", x = -0.8, y = -0.8, label = "🚀", size = 10) +
      annotate("text", x = 0.8, y = -0.8, label = "📊", size = 10) +
      # Clean theme
      theme_void() +
      theme(
        plot.background = element_rect(fill = "white", color = "darkgreen", size = 3),
        plot.margin = margin(20, 20, 20, 20)
      ) +
      coord_fixed() +
      xlim(-1, 1) + ylim(-1, 1)
    
    # Display the celebration plot
    print(success_plot)
    
    # Final success message
    cat("\n🎉 COMPLETE SUCCESS! 🎉\n")
    cat("=" %>% rep(50) %>% paste(collapse = ""), "\n")
    cat("✅ All packages are installed and working correctly!\n")
    cat("✅ tidyverse: Data manipulation and plotting ✓\n")
    cat("✅ easystats: Statistical analysis and modeling ✓\n")
    cat("✅ lme4 + lmerTest: Mixed-effects models ✓\n")
    cat("✅ patchwork: Plot combination ✓\n")
    cat("✅ fitdistrplus: Distribution fitting ✓\n")
    cat("✅ PupillometryR: Eye-tracking analysis ✓\n")
    cat("\nYou're ready for all DevStart tutorials! 🚀\n")
    
  }, error = function(e) {
    cat("❌ FUNCTIONALITY TEST FAILED!\n\n")
    cat("Error during testing:", conditionMessage(e), "\n\n")
    cat("The packages loaded but there might be compatibility issues.\n")
    cat("Try restarting R and running individual tests.\n")
  })
  
} else {
  # Show which packages failed
  cat("❌ SOME PACKAGES FAILED TO LOAD\n\n")
  
  failed_packages <- names(test_results)[!unlist(test_results)]
  successful_packages <- names(test_results)[unlist(test_results)]
  
  if (length(successful_packages) > 0) {
    cat("✅ Successfully loaded packages:\n")
    for (pkg in successful_packages) {
      cat(sprintf("   - %s\n", pkg))
    }
    cat("\n")
  }
  
  if (length(failed_packages) > 0) {
    cat("❌ Failed to load packages:\n")
    for (pkg in failed_packages) {
      cat(sprintf("   - %s\n", pkg))
    }
    cat("\n")
  }
  
  cat("🔧 TROUBLESHOOTING STEPS:\n")
  cat("=" %>% rep(30) %>% paste(collapse = ""), "\n")
  cat("1. Install missing packages:\n")
  cat(sprintf("   install.packages(c('%s'))\n", paste(failed_packages, collapse = "', '")))
  cat("\n2. If installation fails, try:\n")
  cat("   - Update R to the latest version\n")
  cat("   - Update RStudio\n")
  cat("   - Try installing packages one by one\n")
  cat("   - Check your internet connection\n")
}

cat("\n--- Test completed ---\n")