map_ranges_to_runs <- S4Vectors:::map_ranges_to_runs
map_positions_to_runs <- S4Vectors:::map_positions_to_runs

test_map_ranges_to_runs <- function()
{
    test_all_methods <- function(target, run_lens, start, width)
    {
        for (method in 1:3) {
            current <- map_ranges_to_runs(run_lens, start, width, method)
            checkIdentical(target, current)
        }
        current <- map_ranges_to_runs(run_lens, start, width)
        checkIdentical(target, current)
    }

    ## 0 range to map
    target <- list(integer(0), integer(0), integer(0), integer(0))
    test_all_methods(target, integer(0), integer(0), integer(0))
    test_all_methods(target, 15:10, integer(0), integer(0))

    ## 1 range to map
    target <- list(0L, 6L, 0L, 0L)
    test_all_methods(target, 15:10, 1L, sum(15:10))

    target <- list(0L, 1L, 0L, 12L)
    test_all_methods(target, 15:10, 1L, 3L)
    target <- list(0L, 1L, 12L, 0L)
    test_all_methods(target, 15:10, 13L, 3L)
    target <- list(0L, 2L, 13L, 13L)
    test_all_methods(target, 15:10, 14L, 3L)
    target <- list(0L, 2L, 13L, 0L)
    test_all_methods(target, 15:10, 14L, 16L)
    target <- list(0L, 3L, 14L, 12L)
    test_all_methods(target, 15:10, 15L, 16L)
    target <- list(1L, 2L, 0L, 11L)
    test_all_methods(target, 15:10, 16L, 16L)
    target <- list(5L, 1L, 8L, 0L)
    test_all_methods(target, 15:10, 74L, 2L)
    target <- list(5L, 1L, 0L, 0L)
    test_all_methods(target, 15:10, 66L, 10L)

    target <- list(1L, 3L, 11L, 2L)
    test_all_methods(target, c(9L, 15L, 17L, 11L), 21L, 30L)

    ## more than 1 range to map
    start <- 74:1
    width <- rep.int(2L, length(start))
    current <- map_ranges_to_runs(15:10, start, width)
    for (i in seq_along(start)) {
        target_i <- map_ranges_to_runs(15:10, start[i], width[i])
        current_i <- lapply(current, `[[`, i)
        checkIdentical(target_i, current_i)
    }
}

test_map_positions_to_runs <- function()
{
    test_all_methods <- function(run_lens, pos)
    {
        run_breakpoints <- cumsum(run_lens)
        target <- findInterval(pos - 1L, run_breakpoints) + 1L

        width <- rep.int(1L, length(pos))
        for (method in 1:3) {
            current <- map_positions_to_runs(run_lens, pos, method)
            checkIdentical(target, current)
            current <- map_ranges_to_runs(run_lens, pos, width, method)
            checkIdentical(target, current[[1L]] + 1L)
        }
        current <- map_positions_to_runs(run_lens, pos)
        checkIdentical(target, current)
        current <- map_ranges_to_runs(run_lens, pos, width)
        checkIdentical(target, current[[1L]] + 1L)

    }

    test_all_methods(integer(0), integer(0))
    test_all_methods(15:10, integer(0))
    test_all_methods(15:10, seq_len(sum(15:10)))
    test_all_methods(15:10, rev(seq_len(sum(15:10))))
}

