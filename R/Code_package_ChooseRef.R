########### Functions of ChooseRef package ######################



#' Similarities between references and restoration
#' 
#' @description Calculate similarity between plots and a group of reference plot
#'
#' @param RELEVES Variables of the restoration sites data matrix
#' @param REF Variables surveys of reference sites data matrix
#' @param METHOD Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis". Cf \code{\link{vegan::vegdist}} function 
#' @param BINARY (FALSE by default) If TRUE, data is converted to binary data
#' @param DUPLICATES (TRUE by default) Are REF plots also included in RELEVES data?
#' 
#' @return \item{Diss_Mean}{average dissimilarity by RELEVES}
#' @return \item{Diss_Min}{minimal dissimilarity by RELEVES}
#' @return \item{RelRef_order}{list of REF names in ascending order of dissimilarity for each RELEVES}
#' @return \item{DistRef_order}{list of REF distance values in ascending order of dissimilarity for each RELEVES}
#' 
#' @export
#' 
#' @examples # ------------  Creating the data needed for the example --------------------
#' library(vegan)
#' data("dune") #downloading of dune data (cf vegan)
#' data("dune.env") #downloading of dune.env data (cf vegan)
#' # keeping only the numeric variables :
#' dune.env <- data.frame(A1 = dune.env$A1, 
#'                        Moisture =  as.numeric(as.vector(dune.env$Moisture)),
#'                        Manure = as.numeric(as.vector(dune.env$Manure)))
#'                                           
#' # Creating a vector indicating which plots are the potential references and which ones are the restored sites
#' sites <- factor(c(rep("Rest",5),rep("Ref",15)))
#' # Creating a vector with the plot names
#' sites_names <- paste(sites,c(1:5,1:15))
#'
#' #Poviding names to rows (useful for the outputs)
#' row.names(dune.env) <- sites_names
#' row.names(dune) <- sites_names
#' dune.envRest <- dune.env[sites=="Rest",]
#' dune.envRef <- dune.env[sites=="Ref",]
#' 
#' # --------------  Calculating reference dissimilarities ---------------------
#' Distances <- DissRef3(RELEVES = dune.envRest, REF = dune.envRef, METHOD = "euclidean", DUPLICATES = FALSE)
#' Distances
#' 
#' # ----------------  Plotting reference dissimilarities ----------------------
#' Diss_Ref_Plot(RELEVES = dune.envRest, REF = dune.envRef, DISTANCES = Distances, LINK_NUMBER = "N_REF", N_REF = 3)
DissRef3 <- function (RELEVES, REF, METHOD = "bray", BINARY = FALSE, DUPLICATES = TRUE)
{
  # Fusion of tables:
  if(all.equal(names(RELEVES),names(REF))!=TRUE){warning("The two tables do not have the same variables")}
  Tableaux <- rbind(RELEVES,REF)
  
  # Dissimilarity matrix:
  matrice <- as.matrix(vegdist(Tableaux, method = METHOD, binary = BINARY))
  
  # Creation of variables:
  Diss_Mean <- NULL # average dissimilarity of RELEVES
  Diss_Min <- NULL # minimal dissimilarity of RELEVES
  
  # Rel_order will be, for each RELEVES, the list of the names of REF in ascending order of dissimilarity
  Rel_order <- data.frame(matrix(data = 0, nrow = nrow(RELEVES), ncol = nrow(REF)))
  row.names(Rel_order) <- row.names(RELEVES) # line names = RELEVES names
  
  # ‘Dist_order’: for each RELEVES, the list of distance values to REF in ascending order of dissimilarity
  # The object is identical to ‘Rel_order’
  Dist_order <- Rel_order
  
  # loop to calculate for each relevé
  for(i in 1:nrow(RELEVES))
  {
    # Ref distances (from the distance matrix)
    diss_qn <- matrice[i,c((nrow(RELEVES)+1):(nrow(RELEVES)+nrow(REF)))]
    
    # ‘if’ & ‘else’: to delete a zero in the distance matrix when there is a duplicate of data in RELEVES and in REF (in this case set DUPLICATES = TRUE) and we have a value 0 between the RELEVES i and at  least one of REF
    
    if(DUPLICATES==TRUE & min(diss_qn)==0)
    {
      diss_qnw0_2 <- diss_qn[-c(1:length(diss_qn))[diss_qn==0][1]]
    } else {
      diss_qnw0_2 <- diss_qn
    }
    
    # Placing the REF names in ascending order according to the distance to the RELEVES i
    Rel_order[i,] <- factor(row.names(REF)[order(diss_qn)])
    
    # Placing the distance values in ascending order according to the distance to RELEVES i
    Dist_order[i,] <- diss_qn[order(diss_qn)]
    
    # Average distance value
    Diss_Mean[i] <- mean(diss_qnw0_2)
    
    # Minimum distance value
    Diss_Min[i] <- min(diss_qnw0_2)
    
  }
  
  Output <- list(Diss_Mean,Diss_Min,Rel_order,Dist_order)
  names(Output) <- c("Diss_Mean","Diss_Min","RelRef_order","DistRef_order")
  return(Output)
}



#' Plot of references and restoration site dissimilarities
#' @description Display similarity link (previously calculated with \code{\link{ChosseRef::DissRef3}} function) between plots and a group of reference plot
#' 
#' @param RELEVES  Variables of the restoration sites data matrix
#' @param REF Variables of reference sites data matrix
#' @param DISTANCES Object resulting from the analysis of the \code{\link{ChosseRef::DissRef3}} function.
#' @param METHOD Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis". Cf \code{\link{vegan::vegdist}} function 
#' @param COUL_RELEVES Color of RELEVES names
#' @param COUL_Rel_variable If TRUE: color change for each RELEVES
#' @param COUL_REF Color of REF names
#' @param COUL_Seg Color of links
#' @param COUL_Seg_variable If TRUE: colors change for each link
#' @param LINK_NUMBER Number of link drawn. If "N_REF" then for each RELEVES, the number of link drawn is the number given in \code{N_REF}. If "DIST_MIN"
#' only reference that have distance to restoration sites lower than the value given in \code{DIST_MIN} are drawn.
#' @param N_REF Number of REF to display a link between REF and RELEVES
#' @param DIST_MIN Minimum distance value to draw a link between REF and RELEVES
#' @param VAL_DIST Providing or not (TRUE by default) the distance value display
#' @param DECAL Distance between the printing of the distance and the link
#' 
#' @export
#' 
#' @examples # ------------  Creating the data needed for the example --------------------
#' library(vegan)
#' data("dune") #downloading of dune data (cf vegan)
#' data("dune.env") #downloading of dune.env data (cf vegan)
#' # keeping only the numeric variables :
#' dune.env <- data.frame(A1 = dune.env$A1, 
#'                        Moisture =  as.numeric(as.vector(dune.env$Moisture)),
#'                        Manure = as.numeric(as.vector(dune.env$Manure)))
#'                                           
#' # Creating a vector indicating which plots are the potential references and which ones are the restored sites
#' sites <- factor(c(rep("Rest",5),rep("Ref",15)))
#' # Creating a vector with the plot names
#' sites_names <- paste(sites,c(1:5,1:15))
#'
#' #Poviding names to rows (useful for the outputs)
#' row.names(dune.env) <- sites_names
#' row.names(dune) <- sites_names
#' dune.envRest <- dune.env[sites=="Rest",]
#' dune.envRef <- dune.env[sites=="Ref",]
#' 
#' # --------------  Calculating reference dissimilarities ---------------------
#' Distances <- DissRef3(RELEVES = dune.envRest, REF = dune.envRef, METHOD = "euclidean", DUPLICATES = FALSE)
#' Distances
#' 
#' # ----------------  Plotting reference dissimilarities ----------------------
#' Diss_Ref_Plot(RELEVES = dune.envRest, REF = dune.envRef, DISTANCES = Distances, LINK_NUMBER = "N_REF", N_REF = 3)
Diss_Ref_Plot <- function(RELEVES, REF, DISTANCES,
                          METHOD = "euclidean",
                          COUL_RELEVES = 2, COUL_Rel_variable = TRUE,
                          COUL_REF = 1,
                          COUL_Seg = "#9C8809", COUL_Seg_variable = TRUE,
                          LINK_NUMBER="absent", N_REF="absent", 				  	   DIST_MIN="absent",
                          VAL_DIST = TRUE, DECAL = 0 )
  
{
  
  # Check that the names of RELEVES and REF are the same:
  if(all.equal(names(RELEVES),
               names(REF))!=TRUE){warning("REF and RELEVES don't have the same variables")}
  
  # Verification that the information of CHOICE NUMBER, N_REF and DIST_MIN coincide well:
  if(LINK_NUMBER!="N_REF" & LINK_NUMBER!="DIST_MIN"){warning("LINK_NUMBER is not correctly assigned")}
  if(LINK_NUMBER=="N_REF" & N_REF=="absent"){warning("N_REF is not correctly assigned")}
  if(LINK_NUMBER=="DIST_MIN" & DIST_MIN=="absent"){warning("DIST_MIN is not correctly assigned")}
  
  # Verification that DIST_MIN is greater than the minimum distance for each RELEVES:
  if(max(DISTANCES$DistRef_order[,1])>DIST_MIN){warning("DIST_MIN is too low and there is no reference close enough for each RELEVES")}
  
  # Fusion of tables
  Tableaux <- rbind(RELEVES,REF)
  
  # NMDS analysis:
  NMDS <- metaMDS(Tableaux, distance = METHOD, trace = FALSE)
  
  # Display samples of RELEVES and REF
  plot(NMDS$points, type="n", main="Choices of reference
       with environmental conditions")
  text(NMDS$points, labels = row.names(Tableaux),
       col = c(rep(COUL_RELEVES, nrow(RELEVES)), rep(COUL_REF, nrow(REF))))
  # Creation of the color palette if COUL_Seg_variable = TRUE:
  if(COUL_Seg_variable==TRUE){
    COUL_Seg <- colorRampPalette(colors=c("#0789e0","#9C093F","#9C0972","#91099C","#57099C","#26099C","#09199C","#09459C","#09819C","#099C8C","#099C5C","#099C26","#199C09","#4A9C09","#7B9C09","#9C9309","#9C6209","#9C3109","#9C0909"))(nrow(DISTANCES$RelRef_order))}else{COUL_Seg <- rep(COUL_Seg,nrow(DISTANCES$RelRef_order))
    }
  
  # Loop for each RELEVES:
  for (j in 1:nrow(RELEVES))
  {
    # Coordinates of relevé:
    Coo_REL <- NMDS$points[j,c(1:2)]
    
    # If we choose to present a number of references:
    if(LINK_NUMBER=="N_REF"){
      for(i in 1:N_REF)
      {
        Num_ref_proche <- c(1:nrow(REF))[row.names(REF)==DISTANCES$RelRef_order[j,i]]
        Coo_REF <- NMDS$points[nrow(RELEVES)+Num_ref_proche,c(1:2)]
        segments(x0 = Coo_REL[1], y0 = Coo_REL[2], x1 = Coo_REF[1],
                 y1 = Coo_REF[2], col = COUL_Seg[j])
        Dist_REF <- round(DISTANCES$DistRef_order[j,i],2)
        if(VAL_DIST==TRUE){text(x = mean(c(Coo_REF[1],Coo_REL[1]))+DECAL,
                                y = mean(c(Coo_REF[2],Coo_REL[2]))+DECAL,
                                Dist_REF, col="Grey40", cex = 0.7)}
      }
    }
    
    # If you choose a distance value for which to display the links:
    if(LINK_NUMBER=="DIST_MIN"){
      
      # Number of REF with distance less than DIST_MIN:
      N_REF <- length(which(DISTANCES$DistRef_order[1,]<DIST_MIN))
      for(i in 1:N_REF)
      {
        Num_ref_proche <- c(1:nrow(REF))[row.names(REF)==DISTANCES$RelRef_order[j,i]]
        Coo_REF <- NMDS$points[nrow(RELEVES)+Num_ref_proche,c(1:2)]
        segments(x0 = Coo_REL[1], y0 = Coo_REL[2], x1 = Coo_REF[1],
                 y1 = Coo_REF[2], col = COUL_Seg[j])
        Dist_REF <- round(DISTANCES$DistRef_order[j,i],2)
        if(VAL_DIST==TRUE){text(x = mean(c(Coo_REF[1],Coo_REL[1]))+DECAL,
                                y = mean(c(Coo_REF[2],Coo_REL[2]))+DECAL,
                                Dist_REF, col="Grey40", cex = 0.7)}
      }
    }
    
  }
}


