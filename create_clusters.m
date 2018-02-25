function [rv_clusters,rv_singletons] = create_clusters(A,color_size)
    %inefficient way of handling
    %transforming graph by swapping edges and vertices
    [nrows,ncols] = size(A);
    total_edges = 0;
    edge_mapping = zeros(nrows,ncols);
    current_max_deg = 0;
    rv_singletons = [];
    all_row_zeros = [];
    for i = 1:nrows
        for j = i+1:ncols
            if A(i,j) == 1
                total_edges = total_edges + 1;
                rv_singletons = [rv_singletons;zeros(1,color_size)]
                edge_mapping(i,j) = total_edges;
                edge_mapping(j,i) = total_edges;
            end
        end
        row_max_deg = 0;
        all_zeros = 0;
        for j = 1:ncols
            if(i == j)
                continue
            end
            if A(i,j) == 1
                if (all_zeros == 0)
                    all_zeros = 1;
                end
                row_max_deg = row_max_deg + 1;
            end
        end
        if (all_zeros == 1)
            all_row_zeros = [all_row_zeros,1]
        else
            all_row_zeros = [all_row_zeros,0]
        end
        if (current_max_deg < row_max_deg)
            current_max_deg = row_max_deg;
        end
        
    end
    %disp(edge_mapping)
    new_rvs = zeros(total_edges,total_edges);
    
    for i = 1:nrows
        for j = 1:nrows
            rv_index = edge_mapping(i,j);
            if ( rv_index == 0 )
                continue
            end
            for k = j+1:nrows
                if ( edge_mapping(i,k) == 0 )
                    continue
                end
                new_rvs(rv_index,edge_mapping(i,k)) = 1;
                new_rvs(edge_mapping(i,k),rv_index) = 1;
            end
            if(all_row_zeros(i) == 1)
                continue
            end
            
        end
    end 
    %now new_rvs hold the new graph!
    %disp(new_rvs)
    
end