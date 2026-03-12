library(phyloseq)
data(enterotype)
ig <- make_network(enterotype, type = "taxa", max.dist=0.3)

ig |> write_graph("data/enterotype_graph.graphml", format = "graphml")

plot(ig)
