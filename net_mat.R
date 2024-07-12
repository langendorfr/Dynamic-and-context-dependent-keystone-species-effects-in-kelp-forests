#### Author: Ryan E. Langendorf
#### License: GPLv3

net_mat <- function (direct_effects, cutoff = 3, rowVScol = "row", node_specific = FALSE, directed = TRUE, verbose = FALSE) {

    ## From: library("reticulate")
    nx <- import("networkx")
    np <- import("numpy")
    main <- import_main()
    builtins <- import_builtins()

    ## Maintain col cause row orientation for finding paths and building the graph later on
    if (rowVScol == "col") {
        direct_effects = t(direct_effects)
    }

    if (nrow(direct_effects) != ncol(direct_effects)) {
        stop("The input direct effects matrix is not square.")
    }

    direct_effects_np <- np$array(direct_effects, order = "C")

    if (directed == TRUE) {
        graph_nx <- nx$from_numpy_matrix(direct_effects_np, parallel_edges = FALSE, create_using = nx$DiGraph)
    } else {
        graph_nx <- nx$from_numpy_matrix(direct_effects_np, parallel_edges = FALSE, create_using = nx$Graph)
    }

    net_effects <- matrix(NA, nrow = nrow(direct_effects), ncol = ncol(direct_effects))
    rownames(net_effects) = rownames(direct_effects)
    colnames(net_effects) = colnames(direct_effects)

    ## Row -> Col orientation for nx$all_simple_paths()
    for (row in 1:nrow(net_effects)) {
        for (col in 1:ncol(net_effects)) {
            if (verbose == TRUE) {
                print(c(row, col))
            }

            paths_nx <- nx$all_simple_paths(graph_nx, source=row-1, target=col-1, cutoff=cutoff)
            # paths_nx <- nx$all_shortest_paths(graph_nx, source=row-1, target=col-1, cutoff=cutoff)

            ## Some nodes have no paths of length <= cutoff, so return a net effect of zero
            if (length(paths_nx) == 0) {
                net_effects[row, col] = 0
            } else {
                paths_iterated <- reticulate::iterate(paths_nx)

                if (length(paths_iterated) == 0) {
                    net_effects[row, col] = 0
                } else if (row == col) {
                    net_effects[row, col] = 0
                } else {
                    paths <- list()
                    for (p in seq_along(paths_iterated)) {
                        p_collapsed <- do.call(rbind, paths_iterated[[p]]) %>% as.vector()

                        ## Add one because networkx (python) starts counting at zero
                        paths[[p]] = p_collapsed + 1
                    }

                    net_element <- 0
                    path_combined <- NULL
                    for (p in seq_along(paths)) {
                        path_addition <- 1
                        for (i in 2:length(paths[[p]])) {
                            r <- paths[[p]][i-1]
                            c <- paths[[p]][i]
                            path_addition = path_addition * direct_effects[r, c]
                        }

                        net_element = net_element + path_addition
                        path_combined = c(path_combined, path_addition)
                    }

                    net_effects[row, col] = net_element
                }
            }
        }
    }

    if (node_specific == TRUE) {
        net_effects = sweep(net_effects, 1, apply(net_effects, 1, function(x){max(abs(x), na.rm = TRUE)}), FUN = "/")

        ## If a row was all zeros, dividing by the max will result in NaN
        net_effects[is.na(net_effects)] = 0
    }

    ## Maintain col cause row orientation for finding paths and building the graph later on
    if (rowVScol == "col") {
        net_effects = t(net_effects)
    }

    return(net_effects)
}
