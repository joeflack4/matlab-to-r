get_addresses <- function(groups, i, s, d, lim){

# Initiate any sets not so far covered by s
if (!is.na(s)){
    fnames <- s
} else {
    fnames <- ''
}
for (ig in 1:length(groups)){
    gp <- groups[[ig]]
    for (ig2 in 1:length(gp)){
        browser()
        if (!is.na(gp[[ig2]],fnames)){
            (gp[[ig2]]) <- ''
        }
    }
}

if (length[groups] == 1){
    gp1 <- groups[[1]]
    for (ig1 in 1:length(gp1)){
        lim <- lim+1
        (gp1[[ig1]]) <- lim
        (gp1[[ig1]]) <- list(c((gp1[[ig1]]), lim))
        d[[lim]] <- gp1[[ig1]]
    }
}

if (length[groups] == 2){
    gp1 <- groups[[1]]; gp2 <- groups[[2]]

    for (ig1 in 1:length(gp1)){
        for (ig2 in 1:length(gp2)){
            lim <- lim+1
            (gp1[[ig1]])(gp2[[ig2]]) <- lim
            (gp1[[ig1]]) <- list(c(gp1[[ig1]], lim))
            (gp2[[ig2]]) <- list(c(gp2[[ig2]], lim))
            d[[lim]] <- c(gp1[[ig1]], ' ', gp2[[ig2]])
        }
    }
}

if (length(groups) == 3){
    gp1 <- groups[[1]]; gp2 <- groups[[2]]; gp3 <- groups[[3]]
    for (ig1 in 1:length(gp1)){
        for (ig2 in 1:length(gp2)){
            for (ig3 in 1:length(gp3)){
                lim <- lim+1
                (gp1[[ig1]])(gp2[[ig2]])(gp3[[ig3]]) <- lim
                (gp1[[ig1]]) <- list(c((gp1[[ig1]]), lim))
                (gp2[[ig2]]) <- list(c((gp2[[ig2]]), lim))
                (gp3[[ig3]]) <- list(c((gp3[[ig3]]), lim))
                d[[lim]] <- c(gp1[[ig1]],  ' ', gp2[[ig2]],  ' ', gp3[[ig3]])
            }
        }
    }
}

if (length(groups) == 4){
    gp1 <- groups[[1]]; gp2 <- groups[[2]]; gp3 <- groups[[3]]; gp4 <- groups[[4]]
    for (ig1 in 1:length(gp1)){
        for (ig2 in 1:length(gp2)){
            for (ig3 in 1:length(gp3)){
                for (ig4 in 1:length(gp4)){
                    lim <- lim+1
                    (gp1[[ig1]])(gp2[[ig2]])(gp3[[ig3]])(gp4[[ig4]]) <- lim
                    (gp1[[ig1]]) <- list(c((gp1[[ig1]]), lim))
                    (gp2[[ig2]]) <- list(c((gp2[[ig2]]), lim))
                    (gp3[[ig3]]) <- list(c((gp3[[ig3]]), lim))
                    (gp4[[ig4]]) <- list(c((gp4[[ig4]]), lim))
                    d[[lim]] <- c(gp1[[ig1]],  ' ', gp2[[ig2]],  ' ', gp3[[ig3]],  ' ', gp4[[ig4]])
                }
            }
        }
    }
}

if (length(groups) == 5){
    gp1 <- groups[[1]]; gp2 <- groups[[2]]; gp3 <- groups[[3]]; gp4 <- groups[[4]]; gp5 <- groups[[5]]
    for (ig1 in 1:length(gp1)){
        for (ig2 in 1:length(gp2)){
            for (ig3 in 1:length(gp3)){
                for (ig4 in 1:length(gp4)){
                    for (ig5 in 1:length(gp5)){
                    lim <- lim+1
                    (gp1[[ig1]])(gp2[[ig2]])(gp3[[ig3]])(gp4[[ig4]])(gp5[[ig5]]) <- lim
                    (gp1[[ig1]]) <- list(c((gp1[[ig1]]), lim))
                    (gp2[[ig2]]) <- list(c((gp2[[ig2]]), lim))
                    (gp3[[ig3]]) <- list(c((gp3[[ig3]]), lim))
                    (gp4[[ig4]]) <- list(c((gp4[[ig4]]), lim))
                    (gp5[[ig5]]) <- list(c((gp5[[ig5]]), lim))
                    d[[lim]] <- c(gp1[[ig1]],  ' ', gp2[[ig2]],  ' ', gp3[[ig3]],  ' ', gp4[[ig4]],  ' ', gp5[[ig5]])
                    }
                }
            }
        }
    }
}

nstates <- lim
}
