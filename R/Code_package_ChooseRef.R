#

Diss_Ref_Plot <- function(RELEVES, REF, DISTANCES,
                          METHOD = "euclidean",
                          COUL_RELEVES = 2, COUL_Rel_variable = TRUE,
                          COUL_REF = 1,
                          COUL_Seg = "#9C8809", COUL_Seg_variable = TRUE,
                          CHOIX_NOMBRE="absent", N_REF="absent", 				  	   DIST_MIN="absent",
                          VAL_DIST = TRUE, DECAL = 0 )

{

  # Check that the names of RELEVES and REF are the same:
  if(all.equal(names(RELEVES),
               names(REF))!=TRUE){warning("Les deux tableaux n'ont pas les mêmes espèces...")}

  # Verification that the information of CHOICE NUMBER, N_REF and DIST_MIN coincide well:
  if(CHOIX_NOMBRE!="N_REF" & CHOIX_NOMBRE!="DIST_MIN"){warning("La valeur de CHOIX_NOMBRE n'est pas correcte")}
  if(CHOIX_NOMBRE=="N_REF" & N_REF=="absent"){warning("L'argument N_REF n'est pas renseigné")}
  if(CHOIX_NOMBRE=="DIST_MIN" & DIST_MIN=="absent"){warning("L'argument DIST_MIN n'est pas renseigné")

    # Verification that DIST_MIN is greater than the minimum distance for each RELEVES:
    if(max(DISTANCES$DistRef_order[,1])>DIST_MIN){warning("DIST_MIN est trop faible et il n'y a pas une référence suffisamment proche pour chaque RELEVES")}

    # Fusion of tables
    Tableaux <- rbind(RELEVES,REF)

    # NMDS analysis:
    NMDS <- metaMDS(Tableaux, distance = METHOD)

    # Display samples of RELEVES and REF
    plot(NMDS$points, type="n", main="Choix des références en fonction
         des conditions environnementales")
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
      if(CHOIX_NOMBRE=="N_REF"){
        for(i in 1:N_REF)
        {
          Num_ref_proche <- c(1:nrow(REF))[row.names(REF)==DISTANCES$RelRef_order[j,i]]
          Coo_REF <- NMDS$points[nrow(RELEVES)+Num_ref_proche,c(1:2)]
          segments(x0 = Coo_REL[1], y0 = Coo_REL[2], x1 = Coo_REF[1],
                   y1 = Coo_REF[2], col = COUL_Seg[j])
          Dist_REF <- round(DISTANCES$DistRef_order[j,i],2)
          if(VAL_DIST==TRUE){text(x = mean(c(Coo_REF[1],Coo_REL[1]))+Decalx,
                                  y = mean(c(Coo_REF[2],Coo_REL[2]))+Decaly,
                                  Dist_REF, col="Grey40", cex = 0.7)}
        }
      }

      # If you choose a distance value for which to display the links:
      if(CHOIX_NOMBRE=="DIST_MIN"){

        # Number of REF with distance less than DIST_MIN:
        N_REF <- length(which(DISTANCES$DistRef_order[1,]<DIST_MIN))
        for(i in 1:N_REF)
        {
          Num_ref_proche <- c(1:nrow(REF))[row.names(REF)==DISTANCES$RelRef_order[j,i]]
          Coo_REF <- NMDS$points[nrow(RELEVES)+Num_ref_proche,c(1:2)]
          segments(x0 = Coo_REL[1], y0 = Coo_REL[2], x1 = Coo_REF[1],
                   y1 = Coo_REF[2], col = COUL_Seg[j])
          Dist_REF <- round(DISTANCES$DistRef_order[j,i],2)
          if(VAL_DIST==TRUE){text(x = mean(c(Coo_REF[1],Coo_REL[1]))+Decalx,
                                  y = mean(c(Coo_REF[2],Coo_REL[2]))+Decaly,
                                  Dist_REF, col="Grey40", cex = 0.7)}
        }
      }

    }
  }
