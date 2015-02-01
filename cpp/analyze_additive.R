library(ggplot2)

h2s = c(.001, .005, .01, .05, .1)#, .2)
int_pop_sizes = c(1, 10, 100, 1000, 5000, 10000) # just for final population of 10000
fin_pop_sizes = c(1000, 5000, 10000, 20000)

LOD_threshold = 4

# compare various final population sizes across a range of heritabilities
m1 = data.frame(n = rep(fin_pop_sizes, each = length(h2s)), h2 = rep(h2s, length(fin_pop_sizes)), avg_err = -1, std_err = -1, power = -1, avg_err_bootstrap = -1, upper_err_bootstrap = -1, lower_err_bootstrap=-1, power_bootstrap=-1, upper_power_bootstrap=-1, lower_power_bootstrap=-1)

num_bootstrap = 100

print(m1)

for (h in h2s) {
    for (size in fin_pop_sizes) {
        a = read.table(paste('out/sim_output_summary_', size, '_', size, '_1_', h, '_all.txt', sep=''), header=T)
        num_sims = nrow(a)
        print(num_sims)
        num_found = 0
        err = rep(0, num_sims)
        for (sim in 1:num_sims) {
            current = a[sim,]
            if (current$pred_qtl_chr1 == current$qtl_chr1) {
                if (current$LOD_score >= LOD_threshold) {
                    num_found = num_found + 1
                    err[sim] = abs(current$pred_qtl_ps1 - current$qtl_ps1)
                }
            }
            else if (current$LOD_score >= LOD_threshold) {
                print("significant and wrong chromosome\n")
                print(current$LOD_score)
            }
        }
        m1[(m1$n == size)&(m1$h2 == h),]$avg_err = mean(err)
        m1[(m1$n == size)&(m1$h2 == h),]$std_err = sd(err) / sqrt(num_found)
        m1[(m1$n == size)&(m1$h2 == h),]$power = num_found / num_sims

        # bootstrapping
        all_found = rep(0, num_bootstrap)
        all_err = rep(0, num_bootstrap)
        for (i in 1:num_bootstrap) {
            ab = a[sample(1:num_sims, replace=T),]
            num_found = 0
            err = rep(0, num_sims)
            for (sim in 1:num_sims) {
                current = ab[sim,]
                if (current$pred_qtl_chr1 == current$qtl_chr1) {
                    if (current$LOD_score >= LOD_threshold) {
                        num_found = num_found + 1
                        err[sim] = abs(current$pred_qtl_ps1 - current$qtl_ps1)
                    }
                }
            }
            all_err[i] = mean(err)
            all_found[i] = num_found
        }
        m1[(m1$n == size)&(m1$h2 == h),]$avg_err_bootstrap = mean(all_err)
        m1[(m1$n == size)&(m1$h2 == h),]$upper_err_bootstrap = quantile(all_err, .975)
        m1[(m1$n == size)&(m1$h2 == h),]$lower_err_bootstrap = quantile(all_err, .025)
        m1[(m1$n == size)&(m1$h2 == h),]$power_bootstrap = mean(all_found) / num_sims
        m1[(m1$n == size)&(m1$h2 == h),]$upper_power_bootstrap = quantile(all_found, .975) / num_sims
        m1[(m1$n == size)&(m1$h2 == h),]$lower_power_bootstrap = quantile(all_found, .025) / num_sims
    }
}

m1$n = as.factor(m1$n)

# resolution
ggplot(m1, aes(x = h2*100, y = avg_err, colour = n)) + geom_point() +
        theme(panel.background=element_rect(fill="white"), panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
                        axis.line=element_line(), axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black")) +
        xlab("heritability") + ylab("average resolution (bp)")

ggsave("additive_resolution_by_pop_size.pdf")



# compare using different intermediate population sizes for final population size of 10,000
m2 = data.frame(n = rep(int_pop_sizes, each = length(h2s)), h2 = rep(h2s, length(int_pop_sizes)), avg_err = -1, std_err = -1, power = -1, avg_err_bootstrap = -1, upper_err_bootstrap = -1, lower_err_bootstrap=-1, power_bootstrap=-1, upper_power_bootstrap=-1, lower_power_bootstrap=-1)

for (h in h2s) {
    for (size in int_pop_sizes) {
        print("*")
        print(h)
        print(size)
        a = read.table(paste('out/sim_output_summary_10000_', size, '_1_', h, '_all.txt', sep=''), header=T)
        num_sims = nrow(a)
        num_found = 0
        err = rep(0, num_sims)
        for (sim in 1:num_sims) {
            current = a[sim,]
            if (current$pred_qtl_chr1 == current$qtl_chr1) {
                if (current$LOD_score >= LOD_threshold) {
                    num_found = num_found + 1
                    err[sim] = abs(current$pred_qtl_ps1 - current$qtl_ps1)
                }
            }
            else if (current$LOD_score >= LOD_threshold) {
                print("significant and wrong chromosome\n")
                print(current$LOD_score)
            }
        }
        m2[(m2$n == size)&(m2$h2 == h),]$avg_err = mean(err)
        m2[(m2$n == size)&(m2$h2 == h),]$std_err = sd(err) / sqrt(num_found)
        m2[(m2$n == size)&(m2$h2 == h),]$power = num_found / num_sims

        ## bootstrapping
        all_found = rep(0, num_bootstrap)
        all_err = rep(0, num_bootstrap)
        for (i in 1:num_bootstrap) {
            ab = a[sample(1:num_sims, replace=T),]
            num_found = 0
            err = rep(0, num_sims)
            for (sim in 1:num_sims) {
                current = ab[sim,]
                if (current$pred_qtl_chr1 == current$qtl_chr1) {
                    if (current$LOD_score >= LOD_threshold) {
                        num_found = num_found + 1
                        err[sim] = abs(current$pred_qtl_ps1 - current$qtl_ps1)
                    }
                }
            }
            all_err[i] = mean(err)
            all_found[i] = num_found
        }
        m2[(m2$n == size)&(m2$h2 == h),]$avg_err_bootstrap = mean(all_err)
        m2[(m2$n == size)&(m2$h2 == h),]$upper_err_bootstrap = quantile(all_err, .975)
        m2[(m2$n == size)&(m2$h2 == h),]$lower_err_bootstrap = quantile(all_err, .025)
        m2[(m2$n == size)&(m2$h2 == h),]$power_bootstrap = mean(all_found) / num_sims
        m2[(m2$n == size)&(m2$h2 == h),]$upper_power_bootstrap = quantile(all_found, .975) / num_sims
        m2[(m2$n == size)&(m2$h2 == h),]$lower_power_bootstrap = quantile(all_found, .025) / num_sims
    }
}

m2$n = as.factor(m2$n)


ggplot(m2, aes(x = h2*100, y = avg_err, colour = n)) + geom_point() +
        theme(panel.background=element_rect(fill="white"), panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
                        axis.line=element_line(), axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black")) +
        xlab("heritability") + ylab("average resolution (bp)")

ggsave("additive_resolution_by_cross_size.pdf")


save(m1, m2, file = "additive.Rdata") 
