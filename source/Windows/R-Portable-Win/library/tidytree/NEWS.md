# tidytree 0.2.7

+ allow calling `MRCA` with only one node and will return the node itself (2019-08-30, Fri)

# tidytree 0.2.6

+ `nodeid` and `nodelab` methods for converting from label to node number and vice versa (2019-08-09, Fri)
+ allow determine MRCA of a vector of tips (2019-08-08, Thu)

# tidytree 0.2.5

+ convert elements of roxygen documentation to markdown (2019-05-05, Thu)

# tidytree 0.2.4

+ call `child.tbl_tree` instead of `child` in `offspring`, (2019-02-26, Tue)
  so that it works more robust for `data.frame`.

# tidytree 0.2.3

+ more parameter for `offspring` (2019-01-28, Mon)

# tidytree 0.2.2

+ mv vignette to [treedata-book](https://yulab-smu.github.io/treedata-book/) (2019-01-10, Thu)

# tidytree 0.2.1

+ `mutate.tbl_tree` method (2018-12-19, Wed)
  - <https://github.com/GuangchuangYu/tidytree/issues/7>
+ bug fixed in `child` 
  - <https://github.com/GuangchuangYu/tidytree/pull/8>

# tidytree 0.2.0

+ compatible with `tibble` v = 2.0.0 (2018-11-29, Thu)
  - change `as_data_frame` method to `as_tibble` since `as_data_frame` was deprecated in `tibble` and not exported as generics
  
# tidytree 0.1.9

+ `as_data_frame.phylo` works with `phylo$root.edge` (2018-06-13, Wed)

# tidytree 0.1.8

+ force `get_tree_data(treedata)$node` to be integer (2018-04-19, Thu)

# tidytree 0.1.7

+ `get.data`, `[` and `[[` methods (2018-02-26, Mon)
