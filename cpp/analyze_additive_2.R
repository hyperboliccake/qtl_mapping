library(ggplot2)
library(qtl)

h2s = c(.001, .005, .01, .05, .1)
progeny_sizes = c(1000, 5000, 10000)#, 20000)

LOD_thresholds = c(4.4, 5.64)

# compare various final population sizes across a range of heritabilities
m1 = data.frame(
    LOD_threshold = rep(LOD_thresholds, each = length(progeny_sizes) * length(h2s)),
    n = rep(rep(progeny_sizes, each = length(h2s)), length(LOD_thresholds)),
    h2 = rep(rep(h2s, length(progeny_sizes)), length(LOD_thresholds)),
    avg_interval_size = NA,
    interval_ci_upper = NA,
    interval_ci_lower = NA,
    power = NA,
    power_ci_lower = NA,
    power_ci_upper = NA)

num_bootstrap = 100

print(m1)

# get lod scores into format
# chr, ps, sim_score, sim_score...

bayes_cred_int = function(LOD_scores, chr) {

    # compute 10^LOD and scale to have area 1
    LOD10 = 10^(LOD_scores[LOD_scores$chr == chr,]$LOD)

    x = LOD10 / sum(LOD10)

    o = order(x, decreasing = TRUE)
    cs = cumsum(x[o])

    prob = .95
    wh = -1
    for (i in 1:length(x)) {
        if(cs[i] >= prob) {
            wh = i
            break
        }
    }
        
    interval = range(o[1:wh])

    int_start_ps = (LOD_scores[LOD_scores$chr == chr,])[interval[1],]$ps
    int_end_ps = (LOD_scores[LOD_scores$chr == chr,])[interval[2],]$ps

    return(c(int_start_ps, int_end_ps))
}

for (LOD_threshold in LOD_thresholds) {
for (h in h2s) {
    for (size in progeny_sizes) {
        print(h)
        print(size)
        a = read.table(paste('out/sim_output_summary_', size, '_', size, '_1_', h, '_all.txt', sep=''), header=T)
        a_LOD = read.table(paste('out/sim_output_LOD_', size, '_', size, '_1_', h, '_all.txt', sep=''), header=T)

        num_sims = nrow(a)

        found = rep(0, num_sims)
        int_sizes = rep(-1, num_sims)

        for (sim in 1:num_sims) {
            current = a[sim,]
            #current_LOD = a_LOD[,sim+2]
            if (current$LOD_score >= LOD_threshold) {
                LOD_subset = data.frame(chr=a_LOD$chr, ps=a_LOD$ps, LOD=a_LOD[,sim+2])
                interval = bayes_cred_int(LOD_subset, current$pred_qtl_chr1)
                int_sizes[sim] = interval[2] - interval[1]
                if ((current$pred_qtl_chr1 == current$qtl_chr1) && (current$qtl_ps1 >= interval[1]) && (current$qtl_ps1 <= interval[2])) {
                    found[sim] = 1
                }
            }
        }
        avg_int_size = mean(int_sizes[int_sizes != -1])
        m1[(m1$n == size)&(m1$h2 == h)&(m1$LOD_threshold == LOD_threshold),]$avg_interval_size = avg_int_size
        se = sd(int_sizes[int_sizes != -1]) / sqrt(num_sims)
        m1[(m1$n == size)&(m1$h2 == h)&(m1$LOD_threshold == LOD_threshold),]$interval_ci_upper = avg_int_size + se * 1.96 
        m1[(m1$n == size)&(m1$h2 == h)&(m1$LOD_threshold == LOD_threshold),]$interval_ci_lower = avg_int_size - se * 1.96

        ## bootstrap power
        avg_power = rep(0, num_bootstrap)
        for (i in 1:num_bootstrap) {
            avg_power[i] = mean(found[sample(1:num_sims, replace=T)])
        }
        m1[(m1$n == size)&(m1$h2 == h)&(m1$LOD_threshold == LOD_threshold),]$power = mean(avg_power)
        m1[(m1$n == size)&(m1$h2 == h)&(m1$LOD_threshold == LOD_threshold),]$power_ci_upper = quantile(avg_power, .975)
        m1[(m1$n == size)&(m1$h2 == h)&(m1$LOD_threshold == LOD_threshold),]$power_ci_lower = quantile(avg_power, .025)
    }
}
}

m1$n = as.factor(m1$n)

save(m1, file="additive_2.Rdata")
