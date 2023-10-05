module NetworkLayoutLightGraphsExt

if isdefined(Base, :get_extension)
    import NetworkLayout
    import LightGraphs
else
    import ..NetworkLayout
    import ..LightGraphs
end

function NetworkLayout.layout(l::NetworkLayout.AbstractLayout, g::LightGraphs.AbstractGraph)
    NetworkLayout.layout(l, LightGraphs.adjacency_matrix(g))
end
function NetworkLayout.LayoutIterator(l::NetworkLayout.IterativeLayout, g::LightGraphs.AbstractGraph)
    NetworkLayout.LayoutIterator(l, LightGraphs.adjacency_matrix(g))
end

end # module
