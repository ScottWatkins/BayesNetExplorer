function TwoSampleZtest(obs1::Int64, obs2::Int64, total1::Int64, total2::Int64)

    p1 = obs1 / total1
    p2 = obs2 / total2
    p  = (obs1 + obs2) / (total1 + total2)
    q  = 1 - p
    
    Z = (p1 - p2) / (sqrt( p * q * (1/total1 + 1/total2 ) ))
    Z = round(Z, digits=4)

    if Z < 0
        pval = round(2 * ( (cdf(Normal(), Z))), digits=4)
    else
        pval = round(2 * (1 - (cdf(Normal(), Z))), digits=4)
    end

    return Z, pval

end
