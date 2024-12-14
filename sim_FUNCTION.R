library(igraph)
library(ggplot2)
library(dplyr)
library(gridExtra)

# # Define parameter ranges
# n_vertices_range <- c(50)
# p_erdos_range <- c(0.04, 0.1, 0.2, 0.5, 0.7, 0.9, 0.95, 0.99)
# steps_range <- c(1, 4, 6, 10)
# seed_range <- c(155, 240, 400, 600, 921, 3332)

seed_range <- c(921)
n_vertices_range <- c(50)
p_erdos_range <- c(0.4)
steps_range <- c(8)



# Initialize results dataframe
results <- data.frame()

# FUNCTIONS ----

create_color_palette <- function(max_value = 20) {
  colorRampPalette(c(
    "white",  # Very light 
    "#BBDEFB",
    "#90CAF9",
    "#64B5F6",
    "#42A5F5",
    "#2196F3",
    "#1E88E5",
    "#1976D2",
    "#1565C0",
    "#0D47A1"   # Dark blue
  ))(max_value + 1)
}


# Add this function to print node values
print_node_values <- function(g, iteration) {
  cat("\nNode values at iteration", iteration, ":\n")
  
  # Create sequential node names if they don't exist
  if (is.null(V(g)$name)) {
    V(g)$name <- as.character(1:vcount(g))
  }
  
  node_data <- data.frame(
    Node = 1:vcount(g),  # Use node indices
    Value = V(g)$value,
    Failures = V(g)$failures
  )
  print(node_data)
  cat("\n")
}

visit_nodes <- function(g, start_node, steps) {
  # Helper function to check if target can be reached from any visited node
  can_reach_from_visited <- function(g, target, desired_path) {
    visited_nodes <- which(V(g)$value >= 1)
    
    cat("\nChecking reachability for target node:", target, "\n")
    cat("Currently visited nodes:", visited_nodes, "\n")
    cat("Desired path nodes:", desired_path, "\n")
    
    if(length(visited_nodes) == 0) {
      cat("No visited nodes available\n")
      return(FALSE)
    }
    
    valid_nodes <- unique(c(visited_nodes, desired_path))
    cat("Valid nodes for pathfinding:", valid_nodes, "\n")
    
    # Always consider the start node as reachable
    if(target %in% visited_nodes) {
      cat("Target is already visited - reachable\n")
      return(TRUE)
    }
    
    for(source in visited_nodes) {
      cat("\nTrying path from source", source, "to target", target, "\n")
      
      tryCatch({
        paths <- shortest_paths(g, from = source, to = target, mode = "out")
        path_vertices <- unlist(paths$vpath)
        
        if(length(path_vertices) > 0) {
          cat("Found path:", paste(path_vertices, collapse = " -> "), "\n")
          
          if(all(path_vertices %in% valid_nodes)) {
            cat("Path uses only valid nodes - reachable\n")
            return(TRUE)
          } else {
            cat("Path contains invalid nodes:", 
                paste(path_vertices[!path_vertices %in% valid_nodes], collapse = ", "), "\n")
          }
        } else {
          cat("No path found from", source, "to", target, "\n")
        }
        
      }, error = function(e) {
        cat("Error in pathfinding:", e$message, "\n")
      })
    }
    
    cat("Target", target, "is not reachable from any visited node\n")
    return(FALSE)
  }
  
  current_node <- start_node
  # First mark the start node as visited
  if(V(g)$value[start_node] == 0) {
    V(g)$value[start_node] <- 1
    cat("Marking start node", start_node, "as visited\n")
  }
  
  # Generate the desired path
  desired_path <- c(current_node)
  cat("Desired path:", current_node)
  
  for(i in 1:steps) {
    neighbourhood <- neighbors(g, current_node, mode = "all")
    if(length(neighbourhood) == 0) break
    
    neighbourhood <- rep(neighbourhood, 2)
    next_node <- sample(neighbourhood, 1)
    desired_path <- c(desired_path, next_node)
    cat(" ->", next_node)
    current_node <- next_node
  }
  cat("\n")
  
  # Now check reachability and update values for each node in the desired path
  cat("Checking reachability and updating values:\n")
  for(node in desired_path) {
    cat("\nProcessing node", node, "current value:", V(g)$value[node], "\n")
    
    # Check if this node can be reached from any visited node through the expanded subgraph
    if(can_reach_from_visited(g, node, desired_path)) {
      # Increment the value counter for this node
      V(g)$value[node] <- V(g)$value[node] + 1
      cat("Node", node, "is reachable - new value:", V(g)$value[node], "\n")
    } else {
      V(g)$failures[node] <- V(g)$failures[node] + 1
      cat("Node", node, "is not reachable - failure count:", V(g)$failures[node], "\n")
    }
  }
  
  
  
  return(list(
    graph = g,
    desired_path = desired_path
  ))
}


# Function to get the largest SCC nodes
get_largest_scc <- function(g) {
  # Get strongly connected components
  scc <- components(g, mode = "strong")
  
  # Create a table of component sizes
  comp_sizes <- table(scc$membership)
  
  # Find the index of the largest component
  largest_comp_idx <- which.max(comp_sizes)
  
  # Get nodes that belong to the largest component
  largest_scc_nodes <- which(scc$membership == largest_comp_idx)
  
  return(largest_scc_nodes)
}

# Function to get initial node from largest SCC
get_initial_node <- function(g) {
  largest_scc_nodes <- get_largest_scc(g)
  
  # Sample one node from the largest SCC
  largest_scc_nodes <- rep(largest_scc_nodes, 2)
  start_node <- sample(largest_scc_nodes, 1)
  
  cat("Selected start node:", start_node, "\n")
  return(start_node)
}

# Modified all_nodes_visited function with explicit SCC checking
all_nodes_visited <- function(g) {
  largest_scc_nodes <- get_largest_scc(g)
  
  # Check values of nodes in largest SCC
  scc_values <- V(g)$value[largest_scc_nodes]
  
  
  # Create a data frame for better debugging
  scc_status <- data.frame(
    node = largest_scc_nodes,
    value = scc_values
  )
  print(scc_status)
  
  return(all(scc_values >= 1))
}

# Curious node ----
find_highest_degree_node <- function(g) {
  if(sum(V(g)$value >= 1) == 0) {
    return(sample(V(g), 1))
  }
  
  # Get visited nodes and their degrees
  visited_nodes <- sort(which(V(g)$value >= 1))
  subg <- induced_subgraph(g, v=visited_nodes)
  degrees <- degree(subg, mode="all")
  #degrees <- 1/(degrees+1)                         # NOVELTY VS FAMILIARITY
  
  # Create a data frame to handle duplicates
  node_probs <- data.frame(
    node = visited_nodes,
    prob = softmax(degrees)
  )
  
  # Aggregate probabilities for duplicate nodes
  node_probs_agg <- aggregate(prob ~ node, data=node_probs, sum)
  
  # Sample node based on aggregated probabilities
  
  if (length(degrees) == 1) {
    highest_degree_node <- node_probs_agg$node
  } else {
    highest_degree_node <- sample(
      x = node_probs_agg$node,
      size = 1,
      prob = node_probs_agg$prob
    )
  }
  
  cat("\nvisited_nodes:", visited_nodes, "\n") 
  cat("\nDegrees:", degrees, "\n")
  cat("Probabilities:", node_probs_agg$prob, "\n")
  cat("Selected node:", highest_degree_node, "\n")
  
  return(highest_degree_node)
}

# Softmax ----
softmax <- function(x) {
  exp_x <- exp(x - max(x))  # Subtract max for numerical stability
  exp_x / sum(exp_x)
}



# Function to plot the graph with consistent layout ----
plot_graph_state <- function(g, iteration, fixed_layout, color_palette) {
  # Update colors based on current values
  V(g)$color <- color_palette[V(g)$value + 1]
  
  # Create labels showing both value and failures (not used)
  vertex_labels <- paste0(V(g)$name, "\n(", V(g)$value, "/", V(g)$failures, ")")
  
  plot(g,
       layout = fixed_layout,
       vertex.label = V(g)$name,  # Modified labels
       vertex.size = 15,
       vertex.label.cex = 0.8,
       vertex.label.color = "black",
       vertex.label.dist = 0,
       edge.arrow.size = 0.2,
       edge.curved = 0.2,
       main = paste("Iteration", iteration, "\n(value/failures)"))
}

# Define common theme
common_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )


# MAIN LOOPS ----

# Main loop through parameters
for(n_v in n_vertices_range) {
  for(p_e in p_erdos_range) {
    for(steps in steps_range) {
      for(seed in seed_range) {
        cat("\nRunning with parameters:")
        cat("\nn_vertices:", n_v)
        cat("\np_erdos:", p_e)
        cat("\nsteps:", steps)
        cat("\nseed:", seed, "\n")
        
        # Set seed
        set.seed(seed)
        
        # Create network
        g <- erdos.renyi.game(n = n_v, p = p_e, directed = TRUE)
        
        # Initialize values
        V(g)$name <- as.character(1:vcount(g))
        V(g)$value <- 0
        V(g)$failures <- 0
        
        # Generate layout
        fixed_layout <- layout_with_fr(g)
        
        # Create color palette
        color_palette <- create_color_palette()
        
        # Initialize tracking variables
        visited_nodes_count <- c(1)
        total_failures_count <- c(0)
        value_failure_ratios <- c(1)
        Familiarity <- c(0)
        iteration <- 1
        
        # Initialize start node
        start_node <- get_initial_node(g)
        V(g)$value[start_node] <- V(g)$value[start_node] + 1
        
        # Main network exploration loop
        while (!all_nodes_visited(g)) {
          start_node <- find_highest_degree_node(g)
          result <- visit_nodes(g, start_node, steps)
          g <- result$graph
          
          visited_nodes_count <- c(visited_nodes_count, sum(V(g)$value > 0))
          total_failures_count <- c(total_failures_count, sum(V(g)$failures))
          Familiarity <- c(Familiarity, sum(V(g)$value))
          
          iteration <- iteration + 1
          
          # Add timeout condition
          if(iteration > 300) {
            cat("\nTimeout reached for this parameter set\n")
            break
          }
        }
        
        # Create plot_data for current run
        plot_data <- data.frame(
          Iteration = 0:(length(visited_nodes_count)-1),
          VisitedNodes = visited_nodes_count,
          TotalFailures = total_failures_count,
          Familiarity = Familiarity,
          ValueFailureRatio = (visited_nodes_count - total_failures_count- (1/6)*Familiarity )/iteration
        )
        
        # Create the three progress plots
        progression_plot <- ggplot(plot_data, aes(x = Iteration, y = VisitedNodes)) +
          geom_line(color = "#2196F3", size = 1) +
          geom_point(color = "#2196F3", size = 3) +
          labs(
            title = paste("Learning Progress (n=", n_v, ", p=", p_e, ", steps=", steps, ", seed=", seed, ")"),
            x = "Iteration",
            y = "New nodes"
          ) +
          scale_y_continuous(
            limits = c(0,  n_v),
            breaks = seq(0,  n_v, by = 2)
          ) +
          common_theme
        
        failures_plot <- ggplot(plot_data, aes(x = Iteration, y = TotalFailures)) +
          geom_line(color = "#F44336", size = 1) +
          geom_point(color = "#F44336", size = 3) +
          labs(
            title = paste("Failures (n=", n_v, ", p=", p_e, ", steps=", steps, ", seed=", seed, ")"),
            x = "Iteration",
            y = "Number of Failures"
          ) +
          common_theme
        
        familiarity_plot <- ggplot(plot_data, aes(x = Iteration, y = Familiarity)) +
          geom_line(color = "#2196F3", size = 1) +
          geom_point(color = "#2196F3", size = 3) +
          labs(
            title = paste("Familiarity (n=", n_v, ", p=", p_e, ", steps=", steps, ", seed=", seed, ")"),
            x = "Iteration",
            y = "Familiarity"
          ) +
          scale_y_continuous(
            limits = c(0,  max(Familiarity)),
            breaks = seq(0,  max(Familiarity), by = 50)
          ) +
          common_theme
        
        ratio_plot <- ggplot(plot_data, aes(x = Iteration, y = ValueFailureRatio)) +
          geom_line(color = "#4CAF50", size = 1) +
          geom_point(color = "#4CAF50", size = 3) +
          labs(
            title = paste("Interest (n=", n_v, ", p=", p_e, ", steps=", steps, ", seed=", seed, ")"),
            x = "Iteration",
            y = " Learning - Failures"
          ) +
          common_theme
        
        
        # Display progress plots for current run
        grid.arrange(progression_plot, failures_plot, familiarity_plot, ratio_plot, ncol = 1)
        
        # Store results
        final_result <- data.frame(
          n_vertices = n_v,
          p_erdos = p_e,
          steps = steps,
          seed = seed,
          final_iterations = iteration,
          final_visited = tail(visited_nodes_count, 1),
          final_failures = tail(total_failures_count, 1),
          Familiarity = Familiarity,
          timeout = iteration > 300
        )
        
        results <- rbind(results, final_result)
      }
    }
  }
}

# Save results
write.csv(results, "network_analysis_results.csv", row.names = FALSE)

# Create summary visualizations
summary_plots <- list()

# 1. Iterations by probability
summary_plots[[1]] <- ggplot(results, aes(x = p_erdos, y = final_iterations, color = factor(steps))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~steps) +
  labs(title = "Iterations Required by Network Density",
       x = "Erdos-Renyi Probability",
       y = "Number of Iterations",
       color = "Vertices") +
  theme_minimal()

# 2. Success ratio by probability
summary_plots[[2]] <- ggplot(results, aes(x = p_erdos, y = Familiarity, color = factor(steps))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~steps) +
  labs(title = "Familiarity",
       x = "Erdos-Renyi Probability",
       y = "Familiarity",
       color = "Vertices") +
  theme_minimal()

# Display summary plots
grid.arrange(grobs = summary_plots, ncol = 1)