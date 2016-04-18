function res = splitMass(num_mass_pieces, x_list, f_list)
% SPLITMASS Split density function into points with the same partial mass
%
% Given discretized function f_list defined on points of x_list.
% Let M be the total mass of f_list. Cut M into num_mass_pieces
% and return the positions of the pieces. I.e., return
%                           F^{-1}(y)
% where F is the cumulative distribution of f_list, F^{-1} is its
% inverse, and y are equidistant points between 0 and M.
    
    M = sum(f_list);
    F_list = prs_cdf(f_list);
    
    % Make an interpolation of the cdf
    F = griddedInterpolant(F_list); 
    
    % CDF corresponding to the points on the x-Axis that we want to find
    y_list = linspace(M / num_mass_pieces, M, num_mass_pieces);
        
    
        # Find between which points on the x axis y_list crosses F_list  
        # Store results in il ~ "intersection locations"
        il = np.asarray([sum(F_list < y)-1 for y in y_list], dtype=int)
            res = np.asarray([solve_linear(x_list[il[i]], x_list[il[i]+1], 
            F_list[il[i]], F_list[il[i]+1], 
            y_list[i]) for i in range(len(y_list))])
                
                return res

            
            
            
            end
            
            
%==== Function used to calculate the cumulative distribution function ====
function fout = prs_cdf(fin)
    fout = fin
    for i = 1 : 1 : length(fin) - 1
        fin(i+1) = f(i) + f(i+1);
    end
end % quadratic
%==========================================================================

