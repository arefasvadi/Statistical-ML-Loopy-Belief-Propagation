function Z = sum_prod(A,w,its)
    %first we need to create appropriate clusters and singletons 
    %from adjacency matrix
    color_size = length(w)
    [rv_clusters,rv_singletons] = create_clusters(A,color_size)
    
    
end

