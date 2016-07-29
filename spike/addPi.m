function nonDecreasing = addPi(decreasingList)
    nonDecreasing = decreasingList;
    for i=2:length(nonDecreasing)
        while nonDecreasing(i-1) > nonDecreasing(i)
            nonDecreasing(i) = nonDecreasing(i) + pi;
        end
    end
end
