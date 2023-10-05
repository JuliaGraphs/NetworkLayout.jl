module NetworkLayoutGraphsExt

if isdefined(Base, :get_extension)
    import NetworkLayout
    import Graphs
else
    import ..NetworkLayout
    import ..Graphs
end

function NetworkLayout.layout(l::NetworkLayout.AbstractLayout, g::Graphs.AbstractGraph)
    NetworkLayout.layout(l, Graphs.adjacency_matrix(g))
end
function NetworkLayout.LayoutIterator(l::NetworkLayout.IterativeLayout, g::Graphs.AbstractGraph)
    NetworkLayout.LayoutIterator(l, Graphs.adjacency_matrix(g))
end

end # module
