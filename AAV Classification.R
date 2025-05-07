#########################################
#### CLASSIFICATION ANCA-ASSOC VASC  ####
####   Sela Grays - Apr 25 2025      ####
#########################################
rm(list = ls())

## SET WORKING DIRECTORY
#setwd("Desktop/Peds Rheum RA/pedvas/Classification Manuscript/")
getwd()
`%notin%` <- Negate(`%in%`)


library(Gmisc)
library(irr)
library(pROC)
library(caret)
library(MASS)

#### LOAD NEW CLASSIFICATION SHEET ####
data <- read.csv("New classification sheet-USE THIS.csv", header = TRUE)


data$gender <- as.numeric(data$gender)

data$date_enroll <- as.Date(data$date_enroll, "%m/%d/%Y")
date_enroll_no_NA <- data[-which(is.na(data$date_enroll)),]

describeFactors(data$pvas_chest___3)

data$pvas_chest___3[data$lungs_s_node == 1|data$lungs_cav == 1|data$lungs_m_node == 1] = 1

describeFactors(data$pvas_chest___3)

describeFactors(data$pvas_chest___5)

data$pvas_chest___5[data$lungs_mig_infil == 1|data$lungs_fix_infil == 1] = 1

describeFactors(data$pvas_chest___5)

data$pvas_chest___4[data$lungs_effusion == 1] = 1

data$pvas_chest___2[data$lungs_sten == 1] = 1


data_missing <- data[which(data$dcvas_result3 == "" | data$EMAclass == "" | data$diagnosis == ""),]

data <- data[which(data$dcvas_result3 != "" 
                       | data$EMAclass != ""
                       | data$diagnosis != ""),]

data$cns[data$pvas_cns___1 == 1 | 
           data$pvas_cns___2 == 1 | 
           data$pvas_cns___3 == 1 |
           data$pvas_cns___4 == 1 | 
           data$pvas_cns___5 == 1 | 
           data$pvas_cns___6 == 1 | 
           data$pvas_cns___7 == 1 |
           data$pvas_cns___8 == 1 |
           data$pvas_cns___9 == 1] = 1

data$cns[is.na(data$cns)] = 0

data$ren[data$pvas_renal___1 == 1 | 
           data$pvas_renal___2 == 1 | 
           data$pvas_renal___3 == 1 |
           data$pvas_renal___4 == 1 | 
           data$pvas_renal___5 == 1 | 
           data$pvas_renal___6 == 1 | 
           data$pvas_renal___8 == 1 ] = 1

data$ren[is.na(data$ren)] = 0

data$gi[data$pvas_gi___1 == 1 | 
          data$pvas_gi___2 == 1 | 
          data$pvas_gi___3 == 1 | 
          data$pvas_gi___4 == 1] = 1

data$gi[is.na(data$gi)] = 0

data$cardio[data$pvas_cv___1 == 1 | 
              data$pvas_cv___2 == 1 | 
              data$pvas_cv___3 == 1 |
              data$pvas_cv___4 == 1 | 
              data$pvas_cv___5 == 1 | 
              data$pvas_cv___6 == 1 | 
              data$pvas_cv___7 == 1 |
              data$pvas_cv___8 == 1 |
              data$pvas_cv___9 == 1] = 1

data$cardio[is.na(data$cardio)] = 0

data$chest[data$pvas_ent___1 == 1 | 
             data$pvas_ent___2 == 1 | 
             data$pvas_ent___3 == 1 | 
             data$pvas_ent___4 == 1 | 
             data$pvas_ent___5 == 1] = 1

data$chest[is.na(data$chest)] = 0

data$ent[data$pvas_ent___1 == 1 | 
           data$pvas_ent___2 == 1 | 
           data$pvas_ent___3 == 1 | 
           data$pvas_ent___4 == 1 | 
           data$pvas_ent___5 == 1] = 1

data$ent[is.na(data$ent)] = 0

data$muc[data$pvas_mucous2___1 == 1 | 
           data$pvas_mucous2___2 == 1 | 
           data$pvas_mucous2___3 == 1 | 
           data$pvas_mucous2___4 == 1 |
           data$pvas_mucous2___5 == 1 | 
           data$pvas_mucous2___6 == 1 | 
           data$pvas_mucous2___7 == 1 | 
           data$pvas_mucous2___8 == 1 |
           data$pvas_mucous2___9 == 1 |
           data$pvas_mucous2___10 == 1] = 1

data$muc[is.na(data$muc)] = 0

data$cuta[data$pvas_cut2___1 == 1 | 
            data$pvas_cut2___2 == 1 | 
            data$pvas_cut2___3 == 1 | 
            data$pvas_cut2___4 == 1 |
            data$pvas_cut2___5 == 1 | 
            data$pvas_cut2___6 == 1 | 
            data$pvas_cut2___7 == 1 | 
            data$pvas_cut2___8 == 1 |
            data$pvas_cut2___9 == 1] = 1

data$cuta[is.na(data$cuta)] = 0

data$gen[data$pvas_gen1___1 == 1 | 
           data$pvas_gen1___2 == 1 | 
           data$pvas_gen1___3 == 1 | 
           data$pvas_gen1___4 == 1] = 1

data$gen[is.na(data$gen)] = 0


#### SEPERATE DIAGNOSIS ####

not_ACR <- data[which(data$dcvas_result3 == "" 
                      | data$dcvas_result3 == "EGPA"
                      | data$dcvas_result3 == "uncl AAV"
                      | data$dcvas_result3 == "unclassified"),]

all_ACR <- data[which(data$dcvas_result3 != "" 
                       & data$dcvas_result3 != "EGPA"
                       & data$dcvas_result3 != "uncl AAV"
                       & data$dcvas_result3 != "unclassified"),]


ACR <- data[which(data$dcvas_result3 != "" 
                      & data$dcvas_result3 != "GPA/MPA" 
                      & data$dcvas_result3 != "EGPA"
                      & data$dcvas_result3 != "uncl AAV"
                      & data$dcvas_result3 != "unclassified"),]

GPA <- data[which(data$dcvas_result3 != "GPA/MPA" 
                  & (data$dcvas_result3 == "GPA" | data$EMAclass == "GPA"
                     | data$dcvas_result3 == "GPA"
                  | data$diagnosis == 1 | data$diagnosis == 2)),]

GPA_inc <- data[which(data$EMAclass == "GPA"|data$dcvas_result3 == "GPA"|data$dcvas_result3 == "GPA/MPA"
                      |data$diagnosis == 1|data$diagnosis == 2),]

GPA_ACR <- data[which(data$dcvas_result3 == "GPA"),]

ACR_GPA_all <- data[which(data$dcvas_result3 == "GPA" | data$dcvas_result3 == "GPA/MPA"),]

GPA_EMA <- data[which(data$EMAclass == "GPA"),]

MPA <- data[which(data$EMAclass == "MPA"|data$dcvas_result3 == "MPA"|
                    data$diagnosis == 3 | data$diagnosis == 4),]

MPA_inc <-data[which(data$EMAclass == "MPA"|data$dcvas_result3 == "MPA"|data$dcvas_result3 == "GPA/MPA"|
                       data$diagnosis == 3 | data$diagnosis == 4),]

MPA_ACR <- data[which(data$dcvas_result3 == "MPA"),]

MPA_ACR_all <- data[which(data$dcvas_result3 == "MPA" | data$dcvas_result3 == "GPA/MPA"),]

MPA_EMA <- data[which(data$EMAclass == "MPA"),]

both_ACR <- data[which(data$dcvas_result3 == "GPA/MPA" ),]

EMA_UNC <- data[which(data$EMAclass == "UCV AAV" | data$EMAclass == "UCV"),]

UNC_ACR <- data[which(data$dcvas_result3 == "unclassified" | 
                        data$dcvas_result3 == "uncl AAV"),]

EMA <- data[which(data$EMAclass != "" 
                  & data$EMAclass != "PAN" 
                  & data$EMAclass != "EGPA"
                  & data$EMAclass != "UCV"
                  & data$EMAclass != "UCV AAV"),]

#### REPORT FORM VERIFICATION #####

all_gpa01 <- all_ACR[which(all_ACR$pvas_ent___1 == 1 #| data$chr_sinus = 1 not in dataset
                         | all_ACR$pvas_ent___2 == 1),]
gpa_gpa01 <- GPA_ACR[which(GPA_ACR$pvas_ent___1 == 1 #| data$chr_sinus = 1 not in dataset
                        | GPA_ACR$pvas_ent___2 == 1),]
mpa_gpa01 <- MPA_ACR[which(MPA_ACR$pvas_ent___1 == 1 #| data$chr_sinus = 1 not in dataset
                       | MPA_ACR$pvas_ent___2 == 1),]
both_gpa01 <- both_ACR[which(both_ACR$pvas_ent___1 == 1 #| data$chr_sinus = 1 not in dataset
                       | both_ACR$pvas_ent___2 == 1),]
not_gpa01 <- not_ACR[which(not_ACR$pvas_ent___1 == 1 #| data$chr_sinus = 1 not in dataset
                             | not_ACR$pvas_ent___2 == 1),]
data_gpa01 <- data[which(data$pvas_ent___1 == 1 #| data$chr_sinus = 1 not in dataset
                             | data$pvas_ent___2 == 1),]


all_gpa02 <- all_ACR[which(all_ACR$ent_nonpvas___9 == 1 | all_ACR$pvas_ent___3 == 1 
                         | all_ACR$u_air_trach_sten == 1 
                         | all_ACR$pvas_chest___2 == 1),]
gpa_gpa02 <- GPA_ACR[which(GPA_ACR$ent_nonpvas___9 == 1 | GPA_ACR$pvas_ent___3 == 1 
                         | GPA_ACR$u_air_trach_sten == 1 
                         | GPA_ACR$pvas_chest___2 == 1),]
mpa_gpa02 <- MPA_ACR[which(MPA_ACR$ent_nonpvas___9 == 1 | MPA_ACR$pvas_ent___3 == 1 
                         | MPA_ACR$u_air_trach_sten == 1 
                         | MPA_ACR$pvas_chest___2 == 1),]
both_gpa02 <- both_ACR[which(both_ACR$ent_nonpvas___9 == 1 | both_ACR$pvas_ent___3 == 1 
                         | both_ACR$u_air_trach_sten == 1 
                         | both_ACR$pvas_chest___2 == 1),]
not_gpa02 <- not_ACR[which(not_ACR$ent_nonpvas___9 == 1 | not_ACR$pvas_ent___3 == 1 
                             | not_ACR$u_air_trach_sten == 1 
                             | not_ACR$pvas_chest___2 == 1),]
data_gpa02 <- data[which(data$ent_nonpvas___9 == 1 | data$pvas_ent___3 == 1 
                             | data$u_air_trach_sten == 1 
                             | data$pvas_chest___2 == 1),]

all_gpa03 <- all_ACR[which(all_ACR$pvas_ent___4 == 1 | all_ACR$pvas_ent___5 == 1),]
gpa_gpa03 <- GPA_ACR[which(GPA_ACR$pvas_ent___4 == 1 | GPA_ACR$pvas_ent___5 == 1),]
mpa_gpa03 <- MPA_ACR[which(MPA_ACR$pvas_ent___4 == 1 | MPA_ACR$pvas_ent___5 == 1),]
both_gpa03 <- both_ACR[which(both_ACR$pvas_ent___4 == 1 | both_ACR$pvas_ent___5 == 1),]
not_gpa03 <- not_ACR[which(not_ACR$pvas_ent___4 == 1 | not_ACR$pvas_ent___5 == 1),]
data_gpa03 <- data[which(data$pvas_ent___4 == 1 | data$pvas_ent___5 == 1),]


MPO_all <- all_ACR[which(all_ACR$mpo_result == 2 | all_ACR$panca_result == 1),]
MPO_GPA <- GPA_ACR[which(GPA_ACR$mpo_result == 2 | GPA_ACR$panca_result == 1),]
MPO_MPA <- MPA_ACR[which(MPA_ACR$mpo_result == 2 | MPA_ACR$panca_result == 1),]
MPO_both <- both_ACR[which(both_ACR$mpo_result == 2 | both_ACR$panca_result == 1),]
MPO_not <- not_ACR[which(not_ACR$mpo_result == 2 | not_ACR$panca_result == 1),]
MPO_data <- data[which(data$mpo_result == 2 | data$panca_result == 1),]


all_gpa05 <- all_ACR[which(all_ACR$lungs_s_node == 1 | all_ACR$lungs_m_node == 1 
                         | all_ACR$lungs_cav == 1 | all_ACR$lungs_mass == 1),]
gpa_gpa05 <- GPA_ACR[which(GPA_ACR$lungs_s_node == 1 | GPA_ACR$lungs_m_node == 1 
                         | GPA_ACR$lungs_cav == 1 | GPA_ACR$lungs_mass == 1),]
mpa_gpa05 <- MPA_ACR[which(MPA_ACR$lungs_s_node == 1 | MPA_ACR$lungs_m_node == 1 
                         | MPA_ACR$lungs_cav == 1 | MPA_ACR$lungs_mass == 1),]
both_gpa05 <- both_ACR[which(both_ACR$lungs_s_node == 1 | both_ACR$lungs_m_node == 1 
                         | both_ACR$lungs_cav == 1 | both_ACR$lungs_mass == 1),]
not_gpa05 <- not_ACR[which(not_ACR$lungs_s_node == 1 | not_ACR$lungs_m_node == 1 
                             | not_ACR$lungs_cav == 1 | not_ACR$lungs_mass == 1),]
data_gpa05 <- data[which(data$lungs_s_node == 1 | data$lungs_m_node == 1 
                             | data$lungs_cav == 1 | data$lungs_mass == 1),]


all_mpa03 <- all_ACR[which(all_ACR$lungs_fibr == 1 | all_ACR$lungs_ild == 1),]
gpa_mpa03 <- GPA_ACR[which(GPA_ACR$lungs_fibr == 1 | GPA_ACR$lungs_ild == 1),]
mpa_mpa03 <- MPA_ACR[which(MPA_ACR$lungs_fibr == 1 | MPA_ACR$lungs_ild == 1),]
both_mpa03 <- both_ACR[which(both_ACR$lungs_fibr == 1 | both_ACR$lungs_ild == 1),]
not_mpa03 <- not_ACR[which(not_ACR$lungs_fibr == 1 | not_ACR$lungs_ild == 1),]
data_mpa03 <- data[which(data$lungs_fibr == 1 | data$lungs_ild == 1),]

all_gpa08 <- all_ACR[which(all_ACR$pauci == 1),]
gpa_gpa08 <- GPA_ACR[which(GPA_ACR$pauci == 1),]
mpa_gpa08 <- MPA_ACR[which(MPA_ACR$pauci == 1),]
both_gpa08 <- both_ACR[which(both_ACR$pauci == 1),]
not_gpa08 <- not_ACR[which(not_ACR$pauci == 1),]
data_gpa08 <- data[which(data$pauci == 1),]


all_mpa04 <- all_ACR[which(all_ACR$pauci == 1),]
gpa_mpa04 <- GPA_ACR[which(GPA_ACR$pauci == 1),]
mpa_mpa04 <- MPA_ACR[which(MPA_ACR$pauci == 1),]
both_mpa04 <- both_ACR[which(both_ACR$pauci == 1),]
not_mpa04 <- not_ACR[which(not_ACR$pauci == 1),]
data_mpa04 <- data[which(data$pauci == 1),]


all_gpa07 <- all_ACR[which(all_ACR$sinus_fluid == 1| all_ACR$sinus_mass_eff == 1 
                         | all_ACR$sinus_bd == 1 | all_ACR$ent_nonpvas___8 == 1 
                         | all_ACR$pvas_ent___2 == 1),]
gpa_gpa07 <- GPA_ACR[which(GPA_ACR$sinus_fluid == 1| GPA_ACR$sinus_mass_eff == 1 
                         | GPA_ACR$sinus_bd == 1 | GPA_ACR$ent_nonpvas___8 == 1 
                         | GPA_ACR$pvas_ent___2 == 1),]
mpa_gpa07 <- MPA_ACR[which(MPA_ACR$sinus_fluid == 1| MPA_ACR$sinus_mass_eff == 1 
                         | MPA_ACR$sinus_bd == 1 | MPA_ACR$ent_nonpvas___8 == 1 
                         | MPA_ACR$pvas_ent___2 == 1),]
both_gpa07 <- both_ACR[which(both_ACR$sinus_fluid == 1| both_ACR$sinus_mass_eff == 1 
                         | both_ACR$sinus_bd == 1 | both_ACR$ent_nonpvas___8 == 1 
                         | both_ACR$pvas_ent___2 == 1),]
not_gpa07 <- not_ACR[which(not_ACR$sinus_fluid == 1| not_ACR$sinus_mass_eff == 1 
                             | not_ACR$sinus_bd == 1 | not_ACR$ent_nonpvas___8 == 1 
                             | not_ACR$pvas_ent___2 == 1),]
data_gpa07 <- data[which(data$sinus_fluid == 1| data$sinus_mass_eff == 1 
                             | data$sinus_bd == 1 | data$ent_nonpvas___8 == 1 
                             | data$pvas_ent___2 == 1),]

all_gpa06 <- all_ACR[which(all_ACR$granulo == 1| all_ACR$bx_inflamm_skin___4 == 1 
                         | all_ACR$bx_inflamm_sinus___4 == 1 | all_ACR$bx_inflamm_u_air___4 == 1 
                         | all_ACR$bx_inflamm_l_air___4 == 1 | all_ACR$bx_inflamm_gt___4 == 1 
                         | all_ACR$bx_inflamm_cns___4 == 1 | all_ACR$bx_inflamm_pns___4 == 1 
                         | all_ACR$bx_inflamm_oth___4 == 1 | all_ACR$bx_infil_skin___3 == 1 
                         | all_ACR$bx_infil_sinus___3 == 1 | all_ACR$bx_infil_u_air___3 == 1 
                         | all_ACR$bx_infil_l_air___3 == 1 | all_ACR$bx_infil_gt___3 == 1 
                         | all_ACR$bx_infil_cns___3 == 1 | all_ACR$bx_infil_pns___3 == 1 
                         | all_ACR$bx_infil_oth___3 == 1 ),]
gpa_gpa06 <- GPA_ACR[which(GPA_ACR$granulo == 1| GPA_ACR$bx_inflamm_skin___4 == 1 
                         | GPA_ACR$bx_inflamm_sinus___4 == 1 | GPA_ACR$bx_inflamm_u_air___4 == 1 
                         | GPA_ACR$bx_inflamm_l_air___4 == 1 | GPA_ACR$bx_inflamm_gt___4 == 1 
                         | GPA_ACR$bx_inflamm_cns___4 == 1 | GPA_ACR$bx_inflamm_pns___4 == 1 
                         | GPA_ACR$bx_inflamm_oth___4 == 1 | GPA_ACR$bx_infil_skin___3 == 1 
                         | GPA_ACR$bx_infil_sinus___3 == 1 | GPA_ACR$bx_infil_u_air___3 == 1 
                         | GPA_ACR$bx_infil_l_air___3 == 1 | GPA_ACR$bx_infil_gt___3 == 1 
                         | GPA_ACR$bx_infil_cns___3 == 1 | GPA_ACR$bx_infil_pns___3 == 1 
                         | GPA_ACR$bx_infil_oth___3 == 1 ),]
mpa_gpa06 <- MPA_ACR[which(MPA_ACR$granulo == 1| MPA_ACR$bx_inflamm_skin___4 == 1 
                         | MPA_ACR$bx_inflamm_sinus___4 == 1 | MPA_ACR$bx_inflamm_u_air___4 == 1 
                         | MPA_ACR$bx_inflamm_l_air___4 == 1 | MPA_ACR$bx_inflamm_gt___4 == 1 
                         | MPA_ACR$bx_inflamm_cns___4 == 1 | MPA_ACR$bx_inflamm_pns___4 == 1 
                         | MPA_ACR$bx_inflamm_oth___4 == 1 | MPA_ACR$bx_infil_skin___3 == 1 
                         | MPA_ACR$bx_infil_sinus___3 == 1 | MPA_ACR$bx_infil_u_air___3 == 1 
                         | MPA_ACR$bx_infil_l_air___3 == 1 | MPA_ACR$bx_infil_gt___3 == 1 
                         | MPA_ACR$bx_infil_cns___3 == 1 | MPA_ACR$bx_infil_pns___3 == 1 
                         | MPA_ACR$bx_infil_oth___3 == 1 ),]
both_gpa06 <- both_ACR[which(both_ACR$granulo == 1| both_ACR$bx_inflamm_skin___4 == 1 
                         | both_ACR$bx_inflamm_sinus___4 == 1 | both_ACR$bx_inflamm_u_air___4 == 1 
                         | both_ACR$bx_inflamm_l_air___4 == 1 | both_ACR$bx_inflamm_gt___4 == 1 
                         | both_ACR$bx_inflamm_cns___4 == 1 | both_ACR$bx_inflamm_pns___4 == 1 
                         | both_ACR$bx_inflamm_oth___4 == 1 | both_ACR$bx_infil_skin___3 == 1 
                         | both_ACR$bx_infil_sinus___3 == 1 | both_ACR$bx_infil_u_air___3 == 1 
                         | both_ACR$bx_infil_l_air___3 == 1 | both_ACR$bx_infil_gt___3 == 1 
                         | both_ACR$bx_infil_cns___3 == 1 | both_ACR$bx_infil_pns___3 == 1 
                         | both_ACR$bx_infil_oth___3 == 1 ),]
not_gpa06 <- not_ACR[which(not_ACR$granulo == 1| not_ACR$bx_inflamm_skin___4 == 1 
                             | not_ACR$bx_inflamm_sinus___4 == 1 | not_ACR$bx_inflamm_u_air___4 == 1 
                             | not_ACR$bx_inflamm_l_air___4 == 1 | not_ACR$bx_inflamm_gt___4 == 1 
                             | not_ACR$bx_inflamm_cns___4 == 1 | not_ACR$bx_inflamm_pns___4 == 1 
                             | not_ACR$bx_inflamm_oth___4 == 1 | not_ACR$bx_infil_skin___3 == 1 
                             | not_ACR$bx_infil_sinus___3 == 1 | not_ACR$bx_infil_u_air___3 == 1 
                             | not_ACR$bx_infil_l_air___3 == 1 | not_ACR$bx_infil_gt___3 == 1 
                             | not_ACR$bx_infil_cns___3 == 1 | not_ACR$bx_infil_pns___3 == 1 
                             | not_ACR$bx_infil_oth___3 == 1 ),]
data_gpa06 <- data[which(data$granulo == 1| data$bx_inflamm_skin___4 == 1 
                             | data$bx_inflamm_sinus___4 == 1 | data$bx_inflamm_u_air___4 == 1 
                             | data$bx_inflamm_l_air___4 == 1 | data$bx_inflamm_gt___4 == 1 
                             | data$bx_inflamm_cns___4 == 1 | data$bx_inflamm_pns___4 == 1 
                             | data$bx_inflamm_oth___4 == 1 | data$bx_infil_skin___3 == 1 
                             | data$bx_infil_sinus___3 == 1 | data$bx_infil_u_air___3 == 1 
                             | data$bx_infil_l_air___3 == 1 | data$bx_infil_gt___3 == 1 
                             | data$bx_infil_cns___3 == 1 | data$bx_infil_pns___3 == 1 
                             | data$bx_infil_oth___3 == 1 ),]


#### PULL TABLE 3 STATS ####
## Section 1

describeFactors(GPA_ACR$gender)
describeFactors(MPA_ACR$gender)
describeFactors(both_ACR$gender)
describeFactors(all_ACR$gender)
describeFactors(not_ACR$gender)
describeFactors(data$gender)

summary(GPA_ACR$Age..years.)
summary(MPA_ACR$Age..years.)
summary(both_ACR$Age..years.)
summary(all_ACR$Age..years.)
summary(not_ACR$Age..years.)
summary(data$Age..years.)

summary(GPA_ACR$pvas_score)
summary(MPA_ACR$pvas_score)
summary(both_ACR$pvas_score)
summary(all_ACR$pvas_score)
summary(not_ACR$pvas_score)
summary(data$pvas_score)

## Section 2

# Nasal bloody discharge, ulcers, crusting, 
# congestion/ blockage or nasal septal defect/perforation
nrow(gpa_gpa01)
nrow(mpa_gpa01)
nrow(both_gpa01)
nrow(all_gpa01)
nrow(not_gpa01)
nrow(data_gpa01)

# describeFactors(data$dcvas_gpa01)
# describeFactors(GPA$dcvas_gpa01)
# describeFactors(MPA_ACR$dcvas_gpa01)
# describeFactors(both_ACR$dcvas_gpa01)

# Cartilaginous involvement 
nrow(gpa_gpa02)
nrow(mpa_gpa02)
nrow(both_gpa02)
nrow(all_gpa02)
nrow(not_gpa02)
nrow(data_gpa02)

# describeFactors(GPA$dcvas_gpa02)
# describeFactors(MPA_ACR$dcvas_gpa02)
# describeFactors(both_ACR$dcvas_gpa02)

# Conductive or sensorineural hearing loss
GPA_h_loss <- GPA[which(GPA$pvas_ent___4 == 1 & GPA$pvas_ent___5 == 1),]
GPA_c_loss <- GPA[which(GPA$pvas_ent___4 == 1 & GPA$pvas_ent___5 == 0),]
GPA_s_loss <- GPA[which(GPA$pvas_ent___4 == 0 & GPA$pvas_ent___5 == 1),]

nrow(GPA_h_loss) + nrow(GPA_c_loss) + nrow(GPA_s_loss)

# describeFactors(GPA$dcvas_gpa03)
# describeFactors(MPA$dcvas_gpa03)
# describeFactors(both_ACR$dcvas_gpa03)

nrow(gpa_gpa03)
nrow(mpa_gpa03)
nrow(both_gpa03)
nrow(all_gpa03)
nrow(not_gpa03)
nrow(data_gpa03)

## Section 3

#MPO ANCA
describeFactors(GPA_ACR$mpo_result)
describeFactors(MPA_ACR$mpo_result)
describeFactors(both_ACR$mpo_result)
describeFactors(all_ACR$mpo_result)
describeFactors(not_ACR$mpo_result)
describeFactors(data$mpo_result)


# describeFactors(GPA$panca_result)
# describeFactors(GPA$dcvas_gpa09)
# describeFactors(GPA$dcvas_mpa02)


#PR3 ANCA
describeFactors(GPA_ACR$pr3_result)
describeFactors(MPA_ACR$pr3_result)
describeFactors(both_ACR$pr3_result)
describeFactors(all_ACR$pr3_result)
describeFactors(not_ACR$pr3_result)
describeFactors(data$pr3_result)


# describeFactors(GPA_ACR$dcvas_gpa04)
# describeFactors(GPA$dcvas_mpa05)
# describeFactors(MPA$dcvas_mpa05)
# describeFactors(both_ACR$dcvas_mpa05)

# Pulmonary nodules, mass, or cavitation
# describeFactors(GPA$dcvas_gpa05)
# describeFactors(MPA$dcvas_gpa05)
# describeFactors(both_ACR$dcvas_gpa05)


nrow(gpa_gpa05)
nrow(mpa_gpa05)
nrow(both_gpa05)
nrow(all_gpa05)
nrow(not_gpa05)
nrow(data_gpa05)

# Fibrosis or interstitial lung disease 
# describeFactors(GPA$dcvas_mpa03)
# describeFactors(MPA$dcvas_mpa03)
# describeFactors(both_ACR$dcvas_mpa03)

nrow(gpa_mpa03)
nrow(mpa_mpa03)
nrow(both_mpa03)
nrow(all_mpa03)
nrow(not_mpa03)
nrow(data_mpa03)

# Pauci-immune glomerular nephritis
# describeFactors(GPA$dcvas_gpa08)
# describeFactors(GPA$dcvas_mpa04)
# describeFactors(MPA$dcvas_mpa04)
# describeFactors(both_ACR$dcvas_mpa04)


nrow(gpa_mpa04)
nrow(mpa_mpa04)
nrow(both_mpa04)
nrow(all_mpa04)
nrow(not_mpa04)
nrow(data_mpa04)

# Inflammation, consolidation or 
# effusion of nasal/paranasal sinuses or mastoiditis
# describeFactors(GPA$dcvas_gpa07)

nrow(gpa_gpa07)
nrow(mpa_gpa07)
nrow(both_gpa07)
nrow(not_gpa07)
nrow(all_gpa07)
nrow(data_gpa07)

#Granuloma, extravascular granulomatous inflammation or giant cells
# describeFactors(GPA$dcvas_gpa06)

nrow(all_gpa06)
nrow(gpa_gpa06)
nrow(mpa_gpa06)
nrow(both_gpa06)
nrow(not_gpa06)
nrow(data_gpa06)

### FIG 3a & 3b ####

MD_GPA <- data[which((data$diagnosis == 1|data$diagnosis == 2)
                     & data$dcvas_result3 != "GPA/MPA"),]

MD_GPA_only <- data[which((data$diagnosis == 1|data$diagnosis == 2) 
                          & data$EMAclass != "GPA"
                          & data$dcvas_result3 != "GPA"
                          & data$dcvas_result3 != "GPA/MPA"),]

EMA_GPA <- data[which(data$EMAclass == "GPA" & data$dcvas_result3 != "GPA/MPA"),]

EMA_GPA_only <- data[which(data$diagnosis != 1 & data$diagnosis != 2 
                           & data$EMAclass == "GPA"
                           & data$dcvas_result3 != "GPA"
                           & data$dcvas_result3 != "GPA/MPA"),]

ACR_GPA <- data[which(data$dcvas_result3 == "GPA" & data$dcvas_result3 != "GPA/MPA"),]

ACR_GPA_only <- data[which(data$diagnosis != 1 & data$diagnosis != 2 
                           & data$EMAclass != "GPA"
                           & data$dcvas_result3 == "GPA"
                           & data$dcvas_result3 != "GPA/MPA"),]

MD_EMA_GPA <- data[which((data$diagnosis == 1 | data$diagnosis == 2) 
                         & data$EMAclass == "GPA"
                         & data$dcvas_result3 != "GPA"
                         & data$dcvas_result3 != "GPA/MPA"),]

EMA_ACR_GPA <- data[which(data$diagnosis != 1 & data$diagnosis != 2 
                          & data$EMAclass == "GPA"
                          & data$dcvas_result3 == "GPA"
                          & data$dcvas_result3 != "GPA/MPA"),]

ACR_MD_GPA <- data[which((data$diagnosis == 1 | data$diagnosis == 2) 
                         & data$EMAclass != "GPA"
                         & data$dcvas_result3 == "GPA"
                         & data$dcvas_result3 != "GPA/MPA"),]

GPA_all <- data[which(data$EMAclass == "GPA" 
                      & data$dcvas_result3 == "GPA" 
                      & (data$diagnosis == 1 | data$diagnosis == 2)
                      & data$dcvas_result3 != "GPA/MPA"),]

GPA_any <- data[which((data$EMAclass == "GPA" 
                      | data$dcvas_result3 == "GPA" 
                      | data$diagnosis == 1 | data$diagnosis == 2)
                      & data$dcvas_result3 != "GPA/MPA"),]

MD_MPA <- data[which(data$diagnosis == 3|data$diagnosis == 4),]

MD_MPA_only <- data[which((data$diagnosis == 3|data$diagnosis == 4) 
                          & data$EMAclass != "MPA"
                          & data$dcvas_result3 != "MPA"
                          & data$dcvas_result3 != "GPA/MPA"),]

EMA_MPA <- data[which(data$EMAclass == "MPA"),]

EMA_MPA_only <- data[which(data$diagnosis != 3 & data$diagnosis != 4 
                           & data$EMAclass == "MPA"
                           & data$dcvas_result3 != "MPA"
                           & data$dcvas_result3 != "GPA/MPA"),]

ACR_MPA <- data[which(data$dcvas_result3 == "MPA"),]

ACR_MPA_all <- data[which(data$dcvas_result3 == "MPA" | data$dcvas_result3 == "GPA/MPA"),]

ACR_MPA_only <- data[which(data$diagnosis != 3 & data$diagnosis != 4 
                           & data$EMAclass != "MPA"
                           & data$dcvas_result3 == "MPA"
                           & data$dcvas_result3 != "GPA/MPA"),]

MD_EMA_MPA <- data[which((data$diagnosis == 3 | data$diagnosis == 4) 
                         & data$EMAclass == "MPA"
                         & data$dcvas_result3 != "MPA"
                         & data$dcvas_result3 != "GPA/MPA"),]

EMA_ACR_MPA <- data[which(data$diagnosis != 3 & data$diagnosis != 4 
                          & data$EMAclass == "MPA"
                          & data$dcvas_result3 == "MPA"
                          & data$dcvas_result3 != "GPA/MPA"),]

ACR_MD_MPA <- data[which((data$diagnosis == 3 | data$diagnosis == 4) 
                         & data$EMAclass != "MPA"
                         & data$dcvas_result3 == "MPA"
                         & data$dcvas_result3 != "GPA/MPA"),]

MPA_all <- data[which(data$EMAclass == "MPA" & data$dcvas_result3 == "MPA"
                      & data$dcvas_result3 != "GPA/MPA"
                      & (data$diagnosis == 3 | data$diagnosis == 4)
                      & data$dcvas_result3 != "GPA/MPA"),]

MPA_any <- data[which(data$EMAclass == "MPA" | data$dcvas_result3 == "MPA"  
                      | (data$diagnosis == 3 | data$diagnosis == 4)
                      & data$dcvas_result3 != "GPA/MPA"),] 

#### FIG 3a & 3b Alternative 1 ####
MD_GPA_inc <- data[which(data$diagnosis == 1|data$diagnosis == 2),]

MD_GPA_only_inc <- data[which((data$diagnosis == 1|data$diagnosis == 2) 
                              & data$EMAclass != "GPA"
                              & data$dcvas_result3 != "GPA"),]

EMA_GPA_inc <- data[which(data$EMAclass == "GPA"),]

EMA_GPA_only_inc <- data[which(data$diagnosis != 1 & data$diagnosis != 2 
                               & data$EMAclass == "GPA"
                               & data$dcvas_result3 != "GPA"),]

ACR_GPA_inc <- data[which(data$dcvas_result3 == "GPA"),]

ACR_GPA_only_inc <- data[which(data$diagnosis != 1 & data$diagnosis != 2 
                               & data$EMAclass != "GPA"
                               & data$dcvas_result3 == "GPA"),]

MD_EMA_GPA_inc <- data[which((data$diagnosis == 1 | data$diagnosis == 2) 
                             & data$EMAclass == "GPA"
                             & data$dcvas_result3 != "GPA"),]

EMA_ACR_GPA_inc <- data[which(data$diagnosis != 1 & data$diagnosis != 2 
                              & data$EMAclass == "GPA"
                              & data$dcvas_result3 == "GPA"),]

ACR_MD_GPA_inc <- data[which((data$diagnosis == 1 | data$diagnosis == 2) 
                             & data$EMAclass != "GPA"
                             & data$dcvas_result3 == "GPA"),]

GPA_all_inc <- data[which(data$EMAclass == "GPA" 
                          & data$dcvas_result3 == "GPA" 
                          & (data$diagnosis == 1 | data$diagnosis == 2)),]

GPA_any_inc <- data[which((data$EMAclass == "GPA" 
                       | data$dcvas_result3 == "GPA" 
                       | data$diagnosis == 1 | data$diagnosis == 2)),]

MD_MPA_inc <- data[which(data$diagnosis == 3|data$diagnosis == 4),]

MD_MPA_only_inc <- data[which((data$diagnosis == 3|data$diagnosis == 4) 
                              & data$EMAclass != "MPA"
                              & data$dcvas_result3 != "MPA"),]

EMA_MPA_inc <- data[which(data$EMAclass == "MPA"),]

EMA_MPA_only_inc <- data[which(data$diagnosis != 3 & data$diagnosis != 4 
                               & data$EMAclass == "MPA"
                               & data$dcvas_result3 != "MPA"),]

ACR_MPA_inc <- data[which(data$dcvas_result3 == "MPA"),]

ACR_MPA_only_inc <- data[which(data$diagnosis != 3 & data$diagnosis != 4 
                               & data$EMAclass != "MPA"
                               & data$dcvas_result3 == "MPA"),]

MD_EMA_MPA_inc <- data[which((data$diagnosis == 3 | data$diagnosis == 4) 
                             & data$EMAclass == "MPA"
                             & data$dcvas_result3 != "MPA"),]

EMA_ACR_MPA_inc <- data[which(data$diagnosis != 3 & data$diagnosis != 4 
                              & data$EMAclass == "MPA"
                              & data$dcvas_result3 == "MPA"),]

ACR_MD_MPA_inc <- data[which((data$diagnosis == 3 | data$diagnosis == 4) 
                             & data$EMAclass != "MPA"
                             & data$dcvas_result3 == "MPA"),]

MPA_all_inc <- data[which(data$EMAclass == "MPA" & data$dcvas_result3 == "MPA" 
                          & (data$diagnosis == 3 | data$diagnosis == 4)),]

#### FIG 3a & 3b Alternative 2 ####

MD_GPA_inc2 <- data[which(data$diagnosis == 1|data$diagnosis == 2),]

MD_GPA_only_inc2 <- data[which((data$diagnosis == 1|data$diagnosis == 2) 
                          & data$EMAclass != "GPA"
                          & data$dcvas_result3 != "GPA" & data$dcvas_result3 != "GPA/MPA"),]

EMA_GPA_inc2 <- data[which(data$EMAclass == "GPA"),]

EMA_GPA_only_inc2 <- data[which(data$diagnosis != 1 & data$diagnosis != 2 
                           & data$EMAclass == "GPA"
                           & data$dcvas_result3 != "GPA" & data$dcvas_result3 != "GPA/MPA"),]

ACR_GPA_inc2 <- data[which(data$dcvas_result3 == "GPA" | data$dcvas_result3 == "GPA/MPA"),]

ACR_GPA_only_inc2 <- data[which(data$diagnosis != 1 & data$diagnosis != 2 
                           & data$EMAclass != "GPA"
                           & (data$dcvas_result3 == "GPA"|data$dcvas_result3 == "GPA/MPA")),]

MD_EMA_GPA_inc2 <- data[which((data$diagnosis == 1 | data$diagnosis == 2) 
                         & data$EMAclass == "GPA"
                         & data$dcvas_result3 != "GPA"
                         & data$dcvas_result3 != "GPA/MPA"),]

EMA_ACR_GPA_inc2 <- data[which(data$diagnosis != 1 & data$diagnosis != 2 
                          & data$EMAclass == "GPA"
                          & (data$dcvas_result3 == "GPA"|data$dcvas_result3 == "GPA/MPA")),]

ACR_MD_GPA_inc2 <- data[which((data$diagnosis == 1 | data$diagnosis == 2) 
                         & data$EMAclass != "GPA"
                         & (data$dcvas_result3 == "GPA"|data$dcvas_result3 == "GPA/MPA")),]

GPA_all_inc2 <- data[which(data$EMAclass == "GPA" 
                      & (data$dcvas_result3 == "GPA" | data$dcvas_result3 == "GPA/MPA") 
                      & (data$diagnosis == 1 | data$diagnosis == 2)),]

GPA_any_inc_2 <- data[which((data$EMAclass == "GPA" 
                           | data$dcvas_result3 == "GPA" | data$dcvas_result3 == "GPA/MPA"
                           | data$diagnosis == 1 | data$diagnosis == 2)),]

MD_MPA_inc2 <- data[which(data$diagnosis == 3|data$diagnosis == 4),]

MD_MPA_only_inc2 <- data[which((data$diagnosis == 3|data$diagnosis == 4) 
                          & data$EMAclass != "MPA"
                          & data$dcvas_result3 != "MPA" & data$dcvas_result3 != "GPA/MPA"),]

EMA_MPA_inc2 <- data[which(data$EMAclass == "MPA"),]

EMA_MPA_only_inc2 <- data[which(data$diagnosis != 3 & data$diagnosis != 4 
                           & data$EMAclass == "MPA"
                           & data$dcvas_result3 != "MPA" & data$dcvas_result3 != "GPA/MPA"),]

ACR_MPA_inc2 <- data[which(data$dcvas_result3 == "MPA"),]

ACR_MPA_only_inc2 <- data[which(data$diagnosis != 3 & data$diagnosis != 4 
                           & data$EMAclass != "MPA"
                           & (data$dcvas_result3 == "MPA"|data$dcvas_result3 == "GPA/MPA")),]

MD_EMA_MPA_inc2 <- data[which((data$diagnosis == 3 | data$diagnosis == 4) 
                         & data$EMAclass == "MPA"
                         & data$dcvas_result3 != "MPA" & data$dcvas_result3 != "GPA/MPA"),]

EMA_ACR_MPA_inc2 <- data[which(data$diagnosis != 3 & data$diagnosis != 4 
                          & data$EMAclass == "MPA"
                          & (data$dcvas_result3 == "MPA"|data$dcvas_result3 == "GPA/MPA")),]

ACR_MD_MPA_inc2 <- data[which((data$diagnosis == 3 | data$diagnosis == 4) 
                         & data$EMAclass != "MPA"
                         & (data$dcvas_result3 == "MPA"|data$dcvas_result3 == "GPA/MPA")),]

MPA_all_inc2 <- data[which(data$EMAclass == "MPA" 
                      & (data$dcvas_result3 == "MPA"|data$dcvas_result3 == "GPA/MPA")
                      & (data$diagnosis == 3 | data$diagnosis == 4)),]


MPA_any_inc2 <- data[which(data$EMAclass == "MPA" | data$dcvas_result3 == "MPA"  
                      | data$diagnosis == 3 | data$diagnosis == 4
                      | data$dcvas_result3 == "GPA/MPA"),] 


### FIG 3c ####

#B
ACR_GPA_ACR_MPA <- data[which(data$EMAclass != "GPA"
                              & data$dcvas_result3 == "GPA/MPA"),]

#E
EMA_GPA_ACR_MPA_ACR_GPA <- data[which(data$EMAclass == "GPA" 
                                      & data$dcvas_result3 == "GPA/MPA"),]

#D
EMA_GPA_ACR_GPA <- data[which(data$EMAclass == "GPA" 
                              & data$dcvas_result3 == "GPA"),]

#F
EMA_GPA_ACR_MPA <- data[which(data$EMAclass == "GPA" 
                              & data$dcvas_result3 == "MPA"),]
#C
EMA_nonGPA_ACR_MPA <- data[which(data$EMAclass != "GPA" 
                                 & data$dcvas_result3 == "MPA"),]
#A
EMA_nonGPA_ACR_GPA <- data[which(data$EMAclass != "GPA" 
                                 & data$dcvas_result3 == "GPA"),]
#G
EMA_GPA_ACR_nonMPA <- data[which(data$EMAclass == "GPA" 
                                 & data$dcvas_result3 != "MPA"
                                 & data$dcvas_result3 != "GPA"
                                 & data$dcvas_result3 != "GPA/MPA"),]



all_GPAMPA <- data[which(data$EMAclass == "GPA"
                         | data$dcvas_result3 == "GPA"
                         | data$dcvas_result3 == "MPA"
                         | data$dcvas_result3 == "GPA/MPA"),]


describeFactors(data$dcvas_result3)

#d
EMA_GPA_ACR_nonMPA <- data[which(data$EMAclass == "GPA" 
                                 & data$dcvas_result3 != "MPA"
                                 & data$dcvas_result3 != "GPA/MPA"),]

#f
EMA_nonGPA_ACR_MPA <- data[which(data$EMAclass != "GPA" 
                                 & (data$dcvas_result3 == "MPA"|data$dcvas_result3 == "GPA/MPA")),]
#e
EMA_GPA_ACR_MPA_inc <- data[which(data$EMAclass == "GPA" 
                                  & (data$dcvas_result3 == "MPA"|data$dcvas_result3 == "GPA/MPA")),]



##### SUPPLEMENT TABLE VARIABLE INVESTIGATION ####
describeFactors(ACR$dcvas_result3)

describeFactors(EMA$EMAclass)

#### SUPP Table 1 ####
describeFactors(GPA_ACR$gen)
describeFactors(MPA_ACR$gen)
describeFactors(UNC_ACR$gen)

mod_gen_ACR <- glm(as.factor(dcvas_result3) ~ gen + Age..years. + gender, data = ACR, family="binomial")
summary(mod_gen_ACR)
exp(coef(mod_gen_ACR))
exp(confint(mod_gen_ACR))


#Myalgia
describeFactors(GPA_ACR$pvas_gen1___1)
describeFactors(MPA_ACR$pvas_gen1___1)
describeFactors(UNC_ACR$pvas_gen1___1)

mod_myalgia_ACR <- glm(as.factor(dcvas_result3) ~ pvas_gen1___1 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_myalgia_ACR)
exp(coef(mod_myalgia_ACR))
exp(confint(mod_myalgia_ACR))


#Arthralgia or arthritis
describeFactors(GPA_ACR$pvas_gen1___2)
describeFactors(MPA_ACR$pvas_gen1___2)
describeFactors(UNC_ACR$pvas_gen1___2)

mod_arthritis_ACR <- glm(as.factor(dcvas_result3) ~ pvas_gen1___2 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_arthritis_ACR)
exp(coef(mod_arthritis_ACR))
exp(confint(mod_arthritis_ACR))

#Fever >= 38 deg C
describeFactors(GPA_ACR$pvas_gen1___3)
describeFactors(MPA_ACR$pvas_gen1___3)
describeFactors(UNC_ACR$pvas_gen1___3)

mod_fever_ACR <- glm(as.factor(dcvas_result3) ~ pvas_gen1___3 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_fever_ACR)
exp(coef(mod_fever_ACR))
exp(confint(mod_fever_ACR))

#Weight loss >= 5% body weight
describeFactors(GPA_ACR$pvas_gen1___4)
describeFactors(MPA_ACR$pvas_gen1___4)
describeFactors(UNC_ACR$pvas_gen1___4)

mod_weight_ACR <- glm(as.factor(dcvas_result3) ~ pvas_gen1___4 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_weight_ACR)
exp(coef(mod_weight_ACR))
exp(confint(mod_weight_ACR))


describeFactors(GPA_ACR$cuta)
describeFactors(MPA_ACR$cuta)
describeFactors(UNC_ACR$cuta)

mod_cuta_ACR <- glm(as.factor(dcvas_result3) ~ cuta + Age..years. + gender, data = ACR, family="binomial")
summary(mod_cuta_ACR)
exp(coef(mod_cuta_ACR))
exp(confint(mod_cuta_ACR))


#Polymorphous exanthema
describeFactors(GPA_ACR$pvas_cut2___1)
describeFactors(MPA_ACR$pvas_cut2___1)
describeFactors(UNC_ACR$pvas_cut2___1)

mod_polye_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cut2___1 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_polye_ACR)
exp(coef(mod_polye_ACR))
exp(confint(mod_polye_ACR))


#Livedo
describeFactors(GPA_ACR$pvas_cut2___2)
describeFactors(MPA_ACR$pvas_cut2___2)
describeFactors(UNC_ACR$pvas_cut2___2)

mod_liv_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cut2___2 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_liv_ACR)
exp(coef(mod_liv_ACR))
exp(confint(mod_liv_ACR))

#Panniculitis
describeFactors(GPA_ACR$pvas_cut2___3)
describeFactors(MPA_ACR$pvas_cut2___3)
describeFactors(UNC_ACR$pvas_cut2___3)

mod_panni_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cut2___3 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_panni_ACR)
exp(coef(mod_panni_ACR))
exp(confint(mod_panni_ACR))

#Purpura
describeFactors(GPA_ACR$pvas_cut2___4)
describeFactors(MPA_ACR$pvas_cut2___4)
describeFactors(UNC_ACR$pvas_cut2___4)

mod_purp_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cut2___4 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_purp_ACR)
exp(coef(mod_purp_ACR))
exp(confint(mod_purp_ACR))


#Skin nodules
describeFactors(GPA_ACR$pvas_cut2___5)
describeFactors(MPA_ACR$pvas_cut2___5)
describeFactors(UNC_ACR$pvas_cut2___5)

mod_sn_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cut2___5 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_sn_ACR)
exp(coef(mod_sn_ACR))
exp(confint(mod_sn_ACR))

#Infarct (nail edge lesion, splinter haemorrhage)
describeFactors(GPA_ACR$pvas_cut2___6)
describeFactors(MPA_ACR$pvas_cut2___6)
describeFactors(UNC_ACR$pvas_cut2___6)

mod_infa_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cut2___6 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_infa_ACR)
exp(coef(mod_infa_ACR))
exp(confint(mod_infa_ACR))

#Ulcer (full-thickness necrosis)
describeFactors(GPA_ACR$pvas_cut2___7)
describeFactors(MPA_ACR$pvas_cut2___7)
describeFactors(UNC_ACR$pvas_cut2___7)

mod_ulcer_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cut2___7 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_ulcer_ACR)
exp(coef(mod_ulcer_ACR))
exp(confint(mod_ulcer_ACR))

#Gangrene (extensive necrosis)
describeFactors(GPA_ACR$pvas_cut2___8)
describeFactors(MPA_ACR$pvas_cut2___8)
describeFactors(UNC_ACR$pvas_cut2___8)

mod_gang_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cut2___8 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_gang_ACR)
exp(coef(mod_gang_ACR))
exp(confint(mod_gang_ACR))


#Other skin vasculitis (specify below)
describeFactors(GPA_ACR$pvas_cut2___9)
describeFactors(MPA_ACR$pvas_cut2___9)
describeFactors(UNC_ACR$pvas_cut2___9)

mod_osvasc_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cut2___9 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_osvasc_ACR)
exp(coef(mod_osvasc_ACR))
exp(confint(mod_osvasc_ACR))


describeFactors(GPA_ACR$muc)
describeFactors(MPA_ACR$muc)
describeFactors(UNC_ACR$muc)

mod_muc_ACR <- glm(as.factor(dcvas_result3) ~ muc + Age..years. + gender, data = ACR, family="binomial")
summary(mod_muc_ACR)
exp(coef(mod_muc_ACR))
exp(confint(mod_muc_ACR))


#Mouth ulcers/granulomata
describeFactors(GPA_ACR$pvas_mucous2___1)
describeFactors(MPA_ACR$pvas_mucous2___1)
describeFactors(UNC_ACR$pvas_mucous2___1)

mod_mouth_ACR <- glm(as.factor(dcvas_result3) ~ pvas_mucous2___1 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_mouth_ACR)
exp(coef(mod_mouth_ACR))
exp(confint(mod_mouth_ACR))

#Genital ulcers
describeFactors(GPA_ACR$pvas_mucous2___2)
describeFactors(MPA_ACR$pvas_mucous2___2)
describeFactors(UNC_ACR$pvas_mucous2___2)

mod_genital_ACR <- glm(as.factor(dcvas_result3) ~ pvas_mucous2___2 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_genital_ACR)
exp(coef(mod_genital_ACR))
exp(confint(mod_genital_ACR))

#Adnexal inflammation (sialadenitis/ dacryocystitis)
describeFactors(GPA_ACR$pvas_mucous2___3)
describeFactors(MPA_ACR$pvas_mucous2___3)
describeFactors(UNC_ACR$pvas_mucous2___3)

mod_adn_ACR <- glm(as.factor(dcvas_result3) ~ pvas_mucous2___3 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_adn_ACR)
exp(coef(mod_adn_ACR))
exp(confint(mod_adn_ACR))

#Significant proptosis
describeFactors(GPA_ACR$pvas_mucous2___4)
describeFactors(MPA_ACR$pvas_mucous2___4)
describeFactors(UNC_ACR$pvas_mucous2___4)

mod_prop_ACR <- glm(as.factor(dcvas_result3) ~ pvas_mucous2___4 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_prop_ACR)
exp(coef(mod_prop_ACR))
exp(confint(mod_prop_ACR))

#Red eye (Epi)scleritis
describeFactors(GPA_ACR$pvas_mucous2___5)
describeFactors(MPA_ACR$pvas_mucous2___5)
describeFactors(UNC_ACR$pvas_mucous2___5)

mod_redeyeepi_ACR <- glm(as.factor(dcvas_result3) ~ pvas_mucous2___5 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_redeyeepi_ACR)
exp(coef(mod_redeyeepi_ACR))
exp(confint(mod_redeyeepi_ACR))

#Red eye conjunctivitis/ blepharitis/keratitis
describeFactors(GPA_ACR$pvas_mucous2___6)
describeFactors(MPA_ACR$pvas_mucous2___6)
describeFactors(UNC_ACR$pvas_mucous2___6)

mod_redeyeconj_ACR <- glm(as.factor(dcvas_result3) ~ pvas_mucous2___6 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_redeyeconj_ACR)
exp(coef(mod_redeyeconj_ACR))
exp(confint(mod_redeyeconj_ACR))

#Uveitis
describeFactors(GPA_ACR$pvas_mucous2___7)
describeFactors(MPA_ACR$pvas_mucous2___7)
describeFactors(UNC_ACR$pvas_mucous2___7)

mod_uve_ACR <- glm(as.factor(dcvas_result3) ~ pvas_mucous2___7 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_uve_ACR)
exp(coef(mod_uve_ACR))
exp(confint(mod_uve_ACR))

#Blurred vision
describeFactors(GPA_ACR$pvas_mucous2___8)
describeFactors(MPA_ACR$pvas_mucous2___8)
describeFactors(UNC_ACR$pvas_mucous2___8)

mod_bv_ACR <- glm(as.factor(dcvas_result3) ~ pvas_mucous2___8 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_bv_ACR)
exp(coef(mod_bv_ACR))
exp(confint(mod_bv_ACR))

#Sudden visual loss
describeFactors(GPA_ACR$pvas_mucous2___9)
describeFactors(MPA_ACR$pvas_mucous2___9)
describeFactors(UNC_ACR$pvas_mucous2___9)

mod_svl_ACR <- glm(as.factor(dcvas_result3) ~ pvas_mucous2___9 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_svl_ACR)
exp(coef(mod_svl_ACR))
exp(confint(mod_svl_ACR))

#Retinal vasculitis/retinal vessel thrombosis/retinal exudates/haemorrhages
describeFactors(GPA_ACR$pvas_mucous2___10)
describeFactors(MPA_ACR$pvas_mucous2___10)
describeFactors(UNC_ACR$pvas_mucous2___10)

mod_rv_ACR <- glm(as.factor(dcvas_result3) ~ pvas_mucous2___10 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_rv_ACR)
exp(coef(mod_rv_ACR))
exp(confint(mod_rv_ACR))

describeFactors(GPA_ACR$ent)
describeFactors(MPA_ACR$ent)
describeFactors(UNC_ACR$ent)

mod_ent_ACR <- glm(as.factor(dcvas_result3) ~ ent + Age..years. + gender, data = ACR, family="binomial")
summary(mod_ent_ACR)
exp(coef(mod_ent_ACR))
exp(confint(mod_ent_ACR))

#Bloody nasal discharge/crusts/ulcers/granuloma
describeFactors(GPA_ACR$pvas_ent___1)
describeFactors(MPA_ACR$pvas_ent___1)
describeFactors(UNC_ACR$pvas_ent___1)

mod_bnd_ACR <- glm(as.factor(dcvas_result3) ~ pvas_ent___1 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_bnd_ACR)
exp(coef(mod_bnd_ACR))
exp(confint(mod_bnd_ACR))

#Paranasal sinus involvement
describeFactors(GPA_ACR$pvas_ent___2)
describeFactors(MPA_ACR$pvas_ent___2)
describeFactors(UNC_ACR$pvas_ent___2)

mod_psi_ACR <- glm(as.factor(dcvas_result3) ~ pvas_ent___2 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_psi_ACR)
exp(coef(mod_psi_ACR))
exp(confint(mod_psi_ACR))

#Subglottic stenosis/ hoarseness /stridor
describeFactors(GPA_ACR$pvas_ent___3)
describeFactors(MPA_ACR$pvas_ent___3)
describeFactors(UNC_ACR$pvas_ent___3)

mod_sgs_ACR <- glm(as.factor(dcvas_result3) ~ pvas_ent___3 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_sgs_ACR)
exp(coef(mod_sgs_ACR))
exp(confint(mod_sgs_ACR))

#Conductive hearing loss
describeFactors(GPA_ACR$pvas_ent___4)
describeFactors(MPA_ACR$pvas_ent___4)
describeFactors(UNC_ACR$pvas_ent___4)

mod_chl_ACR <- glm(as.factor(dcvas_result3) ~ pvas_ent___4 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_chl_ACR)
exp(coef(mod_chl_ACR))
exp(confint(mod_chl_ACR))

#Sensorineural hearing loss
describeFactors(GPA_ACR$pvas_ent___5)
describeFactors(MPA_ACR$pvas_ent___5)
describeFactors(UNC_ACR$pvas_ent___5)

mod_shl_ACR <- glm(as.factor(dcvas_result3) ~ pvas_ent___5 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_shl_ACR)
exp(coef(mod_shl_ACR))
exp(confint(mod_shl_ACR))


describeFactors(GPA_ACR$chest)
describeFactors(MPA_ACR$chest)
describeFactors(UNC_ACR$chest)

mod_chest_ACR <- glm(as.factor(dcvas_result3) ~ chest + Age..years. + gender, data = ACR, family="binomial")
summary(mod_chest_ACR)
exp(coef(mod_chest_ACR))
exp(confint(mod_chest_ACR))

#Wheeze or expiratory dyspnea
describeFactors(GPA_ACR$pvas_chest___1)
describeFactors(MPA_ACR$pvas_chest___1)
describeFactors(UNC_ACR$pvas_chest___1)

mod_wheeze_ACR <- glm(as.factor(dcvas_result3) ~ pvas_chest___1 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_wheeze_ACR)
exp(coef(mod_wheeze_ACR))
exp(confint(mod_wheeze_ACR))

#Endobronchial/endotracheal involvement
describeFactors(GPA_ACR$pvas_chest___2)
describeFactors(MPA_ACR$pvas_chest___2)
describeFactors(UNC_ACR$pvas_chest___2)

mod_endo_ACR <- glm(as.factor(dcvas_result3) ~ pvas_chest___2 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_endo_ACR)
exp(coef(mod_endo_ACR))
exp(confint(mod_endo_ACR))

#Nodules or cavities 
describeFactors(GPA_ACR$pvas_chest___3)
describeFactors(MPA_ACR$pvas_chest___3)
describeFactors(UNC_ACR$pvas_chest___3)

mod_nodu_ACR <- glm(as.factor(dcvas_result3) ~ pvas_chest___3 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_nodu_ACR)
exp(coef(mod_nodu_ACR))
exp(confint(mod_nodu_ACR))

#Pleural effusion/pleurisy 
describeFactors(GPA_ACR$pvas_chest___4)
describeFactors(MPA_ACR$pvas_chest___4)
describeFactors(UNC_ACR$pvas_chest___4)

mod_pe_ACR <- glm(as.factor(dcvas_result3) ~ pvas_chest___4 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_pe_ACR)
exp(coef(mod_pe_ACR))
exp(confint(mod_pe_ACR))

#Infiltrate
describeFactors(GPA_ACR$pvas_chest___5)
describeFactors(MPA_ACR$pvas_chest___5)
describeFactors(UNC_ACR$pvas_chest___5)

mod_inf_ACR <- glm(as.factor(dcvas_result3) ~ pvas_chest___5 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_inf_ACR)
exp(coef(mod_inf_ACR))
exp(confint(mod_inf_ACR))

#Massive haemoptysis/Alveolar haemorrhage
describeFactors(GPA_ACR$pvas_chest___6)
describeFactors(MPA_ACR$pvas_chest___6)
describeFactors(UNC_ACR$pvas_chest___6)

mod_mahaem_ACR <- glm(as.factor(dcvas_result3) ~ pvas_chest___6 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_mahaem_ACR)
exp(coef(mod_mahaem_ACR))
exp(confint(mod_mahaem_ACR))

#Respiratory failure
describeFactors(GPA_ACR$pvas_chest___7)
describeFactors(MPA_ACR$pvas_chest___7)
describeFactors(UNC_ACR$pvas_chest___7)

mod_rf_ACR <- glm(as.factor(dcvas_result3) ~ pvas_chest___7 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_rf_ACR)
exp(coef(mod_rf_ACR))
exp(confint(mod_rf_ACR))

describeFactors(GPA_ACR$cardio)
describeFactors(MPA_ACR$cardio)
describeFactors(UNC_ACR$cardio)

mod_cardio_ACR <- glm(as.factor(dcvas_result3) ~ cardio + Age..years. + gender, data = ACR, family="binomial")
summary(mod_cardio_ACR)
exp(coef(mod_cardio_ACR))
exp(confint(mod_cardio_ACR))

#Loss of pulses
describeFactors(GPA_ACR$pvas_cv___1)
describeFactors(MPA_ACR$pvas_cv___1)
describeFactors(UNC_ACR$pvas_cv___1)

mod_lossp_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cv___1 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_lossp_ACR)
exp(coef(mod_lossp_ACR))
exp(confint(mod_lossp_ACR))

#Bruits over accessible arteries
describeFactors(GPA_ACR$pvas_cv___2)
describeFactors(MPA_ACR$pvas_cv___2)
describeFactors(UNC_ACR$pvas_cv___2)

mod_bruits_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cv___2 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_bruits_ACR)
exp(coef(mod_bruits_ACR))
exp(confint(mod_bruits_ACR))

#Blood pressure discrepancy
describeFactors(GPA_ACR$pvas_cv___3)
describeFactors(MPA_ACR$pvas_cv___3)
describeFactors(UNC_ACR$pvas_cv___3)

mod_bp_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cv___3 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_bp_ACR)
exp(coef(mod_bp_ACR))
exp(confint(mod_bp_ACR))

#Claudication of extremities
describeFactors(GPA_ACR$pvas_cv___4)
describeFactors(MPA_ACR$pvas_cv___4)
describeFactors(UNC_ACR$pvas_cv___4)

mod_claud_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cv___4 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_claud_ACR)
exp(coef(mod_claud_ACR))
exp(confint(mod_claud_ACR))

#Ischaemic cardiac pain
describeFactors(GPA_ACR$pvas_cv___5)
describeFactors(MPA_ACR$pvas_cv___5)
describeFactors(UNC_ACR$pvas_cv___5)

mod_icp_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cv___5 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_icp_ACR)
exp(coef(mod_icp_ACR))
exp(confint(mod_icp_ACR))

#Cardiomyopathy
describeFactors(GPA_ACR$pvas_cv___6)
describeFactors(MPA_ACR$pvas_cv___6)
describeFactors(UNC_ACR$pvas_cv___6)

mod_cmpy_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cv___6 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_cmpy_ACR)
exp(coef(mod_cmpy_ACR))
exp(confint(mod_cmpy_ACR))

#Congestive cardiac failure
describeFactors(GPA_ACR$pvas_cv___7)
describeFactors(MPA_ACR$pvas_cv___7)
describeFactors(UNC_ACR$pvas_cv___7)

mod_ccf_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cv___7 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_ccf_ACR)
exp(coef(mod_ccf_ACR))
exp(confint(mod_ccf_ACR))

#Valvular heart disease
describeFactors(GPA_ACR$pvas_cv___8)
describeFactors(MPA_ACR$pvas_cv___8)
describeFactors(UNC_ACR$pvas_cv___8)

mod_vhd_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cv___8 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_vhd_ACR)
exp(coef(mod_vhd_ACR))
exp(confint(mod_vhd_ACR))

#Pericarditis
describeFactors(GPA_ACR$pvas_cv___9)
describeFactors(MPA_ACR$pvas_cv___9)
describeFactors(UNC_ACR$pvas_cv___9)

mod_pericard_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cv___9 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_pericard_ACR)
exp(coef(mod_pericard_ACR))
exp(confint(mod_pericard_ACR))

describeFactors(GPA_ACR$gi)
describeFactors(MPA_ACR$gi)
describeFactors(UNC_ACR$gi)

mod_gi_ACR <- glm(as.factor(dcvas_result3) ~ gi + Age..years. + gender, data = ACR, family="binomial")
summary(mod_gi_ACR)
exp(coef(mod_gi_ACR))
exp(confint(mod_gi_ACR))

#Abdominal pain
describeFactors(GPA_ACR$pvas_gi___1)
describeFactors(MPA_ACR$pvas_gi___1)
describeFactors(UNC_ACR$pvas_gi___1)

mod_abdp_ACR <- glm(as.factor(dcvas_result3) ~ pvas_gi___1 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_abdp_ACR)
exp(coef(mod_abdp_ACR))
exp(confint(mod_abdp_ACR))

#Peritonitis
describeFactors(GPA_ACR$pvas_gi___2)
describeFactors(MPA_ACR$pvas_gi___2)
describeFactors(UNC_ACR$pvas_gi___2)

mod_periton_ACR <- glm(as.factor(dcvas_result3) ~ pvas_gi___2 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_periton_ACR)
exp(coef(mod_periton_ACR))
exp(confint(mod_periton_ACR))

#Blood in stools or bloody diarrhea
describeFactors(GPA_ACR$pvas_gi___3)
describeFactors(MPA_ACR$pvas_gi___3)
describeFactors(UNC_ACR$pvas_gi___3)

mod_bloodstool_ACR <- glm(as.factor(dcvas_result3) ~ pvas_gi___3 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_bloodstool_ACR)
exp(coef(mod_bloodstool_ACR))
exp(confint(mod_bloodstool_ACR))

#Bowel ischaemia
describeFactors(GPA_ACR$pvas_gi___4)
describeFactors(MPA_ACR$pvas_gi___4)
describeFactors(UNC_ACR$pvas_gi___4)

mod_bowelisch_ACR <- glm(as.factor(dcvas_result3) ~ pvas_gi___4 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_bowelisch_ACR)
exp(coef(mod_bowelisch_ACR))
exp(confint(mod_bowelisch_ACR))

describeFactors(GPA_ACR$ren)
describeFactors(MPA_ACR$ren)
describeFactors(UNC_ACR$ren)

mod_ren_ACR <- glm(as.factor(dcvas_result3) ~ ren + Age..years. + gender, data = ACR, family="binomial")
summary(mod_ren_ACR)
exp(coef(mod_ren_ACR))
exp(confint(mod_ren_ACR))

#Hypertension >95th centile (for height)
describeFactors(GPA_ACR$pvas_renal___1)
describeFactors(MPA_ACR$pvas_renal___1)
describeFactors(UNC_ACR$pvas_renal___1)

mod_hyt_ACR <- glm(as.factor(dcvas_result3) ~ pvas_renal___1 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_hyt_ACR)
exp(coef(mod_hyt_ACR))
exp(confint(mod_hyt_ACR))

#Proteinuria >0.3 g/24h,>20mg/mmol creatinine
describeFactors(GPA_ACR$pvas_renal___2)
describeFactors(MPA_ACR$pvas_renal___2)
describeFactors(UNC_ACR$pvas_renal___2)

mod_puria_ACR <- glm(as.factor(dcvas_result3) ~ pvas_renal___2 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_puria_ACR)
exp(coef(mod_puria_ACR))
exp(confint(mod_puria_ACR))

#Haematuria >= 2+ or 5 rbc/hpf or red cell casts
describeFactors(GPA_ACR$pvas_renal___3)
describeFactors(MPA_ACR$pvas_renal___3)
describeFactors(UNC_ACR$pvas_renal___3)

mod_hemuria_ACR <- glm(as.factor(dcvas_result3) ~ pvas_renal___3 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_hemuria_ACR)
exp(coef(mod_hemuria_ACR))
exp(confint(mod_hemuria_ACR))

#GFR 50-80ml/min/1.73 m2 *s
describeFactors(GPA_ACR$pvas_renal___4)
describeFactors(MPA_ACR$pvas_renal___4)
describeFactors(UNC_ACR$pvas_renal___4)

mod_GFR1_ACR <- glm(as.factor(dcvas_result3) ~ pvas_renal___4 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_GFR1_ACR)
exp(coef(mod_GFR1_ACR))
exp(confint(mod_GFR1_ACR))

#GFR 15-49 ml/min/1.73 m2 *
describeFactors(GPA_ACR$pvas_renal___5)
describeFactors(MPA_ACR$pvas_renal___5)
describeFactors(UNC_ACR$pvas_renal___5)

mod_GFR2_ACR <- glm(as.factor(dcvas_result3) ~ pvas_renal___5 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_GFR2_ACR)
exp(coef(mod_GFR2_ACR))
exp(confint(mod_GFR2_ACR))

#GFR < 15 ml/min/1.73m2 *
describeFactors(GPA_ACR$pvas_renal___6)
describeFactors(MPA_ACR$pvas_renal___6)
describeFactors(UNC_ACR$pvas_renal___6)

mod_GFR3_ACR <- glm(as.factor(dcvas_result3) ~ pvas_renal___6 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_GFR3_ACR)
exp(coef(mod_GFR3_ACR))
exp(confint(mod_GFR3_ACR))

#Rise in creatinine > 10% or Creatinine clearance (GFR) fall > 25%
describeFactors(GPA_ACR$pvas_renal___8)
describeFactors(MPA_ACR$pvas_renal___8)
describeFactors(UNC_ACR$pvas_renal___8)

mod_creatinine_ACR <- glm(as.factor(dcvas_result3) ~ pvas_renal___8 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_creatinine_ACR)
exp(coef(mod_creatinine_ACR))
exp(confint(mod_creatinine_ACR))


describeFactors(GPA_ACR$cns)
describeFactors(MPA_ACR$cns)
describeFactors(UNC_ACR$cns)

mod_cns_ACR <- glm(as.factor(dcvas_result3) ~ cns + Age..years. + gender, data = ACR, family="binomial")
summary(mod_cns_ACR)
exp(coef(mod_cns_ACR))
exp(confint(mod_cns_ACR))

#Headache
describeFactors(GPA_ACR$pvas_cns___1)
describeFactors(MPA_ACR$pvas_cns___1)
describeFactors(UNC_ACR$pvas_cns___1)

mod_headache_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cns___1 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_headache_ACR)
exp(coef(mod_headache_ACR))
exp(confint(mod_headache_ACR))

#Meningitis/encephalitis
describeFactors(GPA_ACR$pvas_cns___2)
describeFactors(MPA_ACR$pvas_cns___2)
describeFactors(UNC_ACR$pvas_cns___2)

mod_meningitis_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cns___2 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_meningitis_ACR)
exp(coef(mod_meningitis_ACR))
exp(confint(mod_meningitis_ACR))

#Organic confusion/cognitive dysfunction
describeFactors(GPA_ACR$pvas_cns___3)
describeFactors(MPA_ACR$pvas_cns___3)
describeFactors(UNC_ACR$pvas_cns___3)

mod_orgconf_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cns___3 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_orgconf_ACR)
exp(coef(mod_orgconf_ACR))
exp(confint(mod_orgconf_ACR))

#Seizures (not hypertensive)
describeFactors(GPA_ACR$pvas_cns___4)
describeFactors(MPA_ACR$pvas_cns___4)
describeFactors(UNC_ACR$pvas_cns___4)

mod_seizures_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cns___4 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_seizures_ACR)
exp(coef(mod_seizures_ACR))
exp(confint(mod_seizures_ACR))

#Stroke
describeFactors(GPA_ACR$pvas_cns___5)
describeFactors(MPA_ACR$pvas_cns___5)
describeFactors(UNC_ACR$pvas_cns___5)

mod_stroke_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cns___5 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_stroke_ACR)
exp(coef(mod_stroke_ACR))
exp(confint(mod_stroke_ACR))

#Cord lesion
describeFactors(GPA_ACR$pvas_cns___6)
describeFactors(MPA_ACR$pvas_cns___6)
describeFactors(UNC_ACR$pvas_cns___6)

mod_cl_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cns___6 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_cl_ACR)
exp(coef(mod_cl_ACR))
exp(confint(mod_cl_ACR))

#Cranial nerve palsy
describeFactors(GPA_ACR$pvas_cns___7)
describeFactors(MPA_ACR$pvas_cns___7)
describeFactors(UNC_ACR$pvas_cns___7)

mod_cnp_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cns___7 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_cnp_ACR)
exp(coef(mod_cnp_ACR))
exp(confint(mod_cnp_ACR))

#Sensory peripheral neuropathy
describeFactors(GPA_ACR$pvas_cns___8)
describeFactors(MPA_ACR$pvas_cns___8)
describeFactors(UNC_ACR$pvas_cns___8)

mod_spn_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cns___8 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_spn_ACR)
exp(coef(mod_spn_ACR))
exp(confint(mod_spn_ACR))

#Motor mononeuritis multiplex
describeFactors(GPA_ACR$pvas_cns___9)
describeFactors(MPA_ACR$pvas_cns___9)
describeFactors(UNC_ACR$pvas_cns___9)

mod_mmm_ACR <- glm(as.factor(dcvas_result3) ~ pvas_cns___9 + Age..years. + gender, data = ACR, family="binomial")
summary(mod_mmm_ACR)
exp(coef(mod_mmm_ACR))
exp(confint(mod_mmm_ACR))

#### SUPP Table 2 ####
describeFactors(EMA_GPA$gen)
describeFactors(EMA_MPA$gen)
describeFactors(EMA_UNC$gen)

mod_gen_EMA <- glm(as.factor(EMAclass) ~ gen + Age..years. + gender, data = EMA, family="binomial")
summary(mod_gen_EMA)
exp(coef(mod_gen_EMA))
exp(confint(mod_gen_EMA))


#Myalgia
describeFactors(EMA_GPA$pvas_gen1___1)
describeFactors(EMA_MPA$pvas_gen1___1)
describeFactors(EMA_UNC$pvas_gen1___1)

mod_myalgia_EMA <- glm(as.factor(EMAclass) ~ pvas_gen1___1 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_myalgia_EMA)
exp(coef(mod_myalgia_EMA))
exp(confint(mod_myalgia_EMA))


#Arthralgia or arthritis
describeFactors(EMA_GPA$pvas_gen1___2)
describeFactors(EMA_MPA$pvas_gen1___2)
describeFactors(EMA_UNC$pvas_gen1___2)

mod_arthritis_EMA <- glm(as.factor(EMAclass) ~ pvas_gen1___2 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_arthritis_EMA)
exp(coef(mod_arthritis_EMA))
exp(confint(mod_arthritis_EMA))

#Fever >= 38 deg C
describeFactors(EMA_GPA$pvas_gen1___3)
describeFactors(EMA_MPA$pvas_gen1___3)
describeFactors(EMA_UNC$pvas_gen1___3)

mod_fever_EMA <- glm(as.factor(EMAclass) ~ pvas_gen1___3 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_fever_EMA)
exp(coef(mod_fever_EMA))
exp(confint(mod_fever_EMA))

#Weight loss >= 5% body weight
describeFactors(EMA_GPA$pvas_gen1___4)
describeFactors(EMA_MPA$pvas_gen1___4)
describeFactors(EMA_UNC$pvas_gen1___4)

mod_weight_EMA <- glm(as.factor(EMAclass) ~ pvas_gen1___4 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_weight_EMA)
exp(coef(mod_weight_EMA))
exp(confint(mod_weight_EMA))


describeFactors(EMA_GPA$cuta)
describeFactors(EMA_MPA$cuta)
describeFactors(EMA_UNC$cuta)

mod_cuta_EMA <- glm(as.factor(EMAclass) ~ cuta + Age..years. + gender, data = EMA, family="binomial")
summary(mod_cuta_EMA)
exp(coef(mod_cuta_EMA))
exp(confint(mod_cuta_EMA))


#Polymorphous exanthema
describeFactors(EMA_GPA$pvas_cut2___1)
describeFactors(EMA_MPA$pvas_cut2___1)
describeFactors(EMA_UNC$pvas_cut2___1)

mod_polye_EMA <- glm(as.factor(EMAclass) ~ pvas_cut2___1 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_polye_EMA)
exp(coef(mod_polye_EMA))
exp(confint(mod_polye_EMA))


#Livedo
describeFactors(EMA_GPA$pvas_cut2___2)
describeFactors(EMA_MPA$pvas_cut2___2)
describeFactors(EMA_UNC$pvas_cut2___2)

mod_liv_EMA <- glm(as.factor(EMAclass) ~ pvas_cut2___2 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_liv_EMA)
exp(coef(mod_liv_EMA))
exp(confint(mod_liv_EMA))

#Panniculitis
describeFactors(EMA_GPA$pvas_cut2___3)
describeFactors(EMA_MPA$pvas_cut2___3)
describeFactors(EMA_UNC$pvas_cut2___3)

mod_panni_EMA <- glm(as.factor(EMAclass) ~ pvas_cut2___3 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_panni_EMA)
exp(coef(mod_panni_EMA))
exp(confint(mod_panni_EMA))

#Purpura
describeFactors(EMA_GPA$pvas_cut2___4)
describeFactors(EMA_MPA$pvas_cut2___4)
describeFactors(EMA_UNC$pvas_cut2___4)

mod_purp_EMA <- glm(as.factor(EMAclass) ~ pvas_cut2___4 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_purp_EMA)
exp(coef(mod_purp_EMA))
exp(confint(mod_purp_EMA))


#Skin nodules
describeFactors(EMA_GPA$pvas_cut2___5)
describeFactors(EMA_MPA$pvas_cut2___5)
describeFactors(EMA_UNC$pvas_cut2___5)

mod_sn_EMA <- glm(as.factor(EMAclass) ~ pvas_cut2___5 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_sn_EMA)
exp(coef(mod_sn_EMA))
exp(confint(mod_sn_EMA))

#Infarct (nail edge lesion, splinter haemorrhage)
describeFactors(EMA_GPA$pvas_cut2___6)
describeFactors(EMA_MPA$pvas_cut2___6)
describeFactors(EMA_UNC$pvas_cut2___6)

mod_infa_EMA <- glm(as.factor(EMAclass) ~ pvas_cut2___6 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_infa_EMA)
exp(coef(mod_infa_EMA))
exp(confint(mod_infa_EMA))

#Ulcer (full-thickness necrosis)
describeFactors(EMA_GPA$pvas_cut2___7)
describeFactors(EMA_MPA$pvas_cut2___7)
describeFactors(EMA_UNC$pvas_cut2___7)

mod_ulcer_EMA <- glm(as.factor(EMAclass) ~ pvas_cut2___7 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_ulcer_EMA)
exp(coef(mod_ulcer_EMA))
exp(confint(mod_ulcer_EMA))

#Gangrene (extensive necrosis)
describeFactors(EMA_GPA$pvas_cut2___8)
describeFactors(EMA_MPA$pvas_cut2___8)
describeFactors(EMA_UNC$pvas_cut2___8)

mod_gang_EMA <- glm(as.factor(EMAclass) ~ pvas_cut2___8 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_gang_EMA)
exp(coef(mod_gang_EMA))
exp(confint(mod_gang_EMA))


#Other skin vasculitis (specify below)
describeFactors(EMA_GPA$pvas_cut2___9)
describeFactors(EMA_MPA$pvas_cut2___9)
describeFactors(EMA_UNC$pvas_cut2___9)

mod_osvasc_EMA <- glm(as.factor(EMAclass) ~ pvas_cut2___9 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_osvasc_EMA)
exp(coef(mod_osvasc_EMA))
exp(confint(mod_osvasc_EMA))

#Mucous membranes/eyes involvement 
describeFactors(EMA_GPA$muc)
describeFactors(EMA_MPA$muc)
describeFactors(EMA_UNC$muc)

mod_muc_EMA <- glm(as.factor(EMAclass) ~ muc + Age..years. + gender, data = EMA, family="binomial")
summary(mod_muc_EMA)
exp(coef(mod_muc_EMA))
exp(confint(mod_muc_EMA))

#Mouth ulcers/granulomata
describeFactors(EMA_GPA$pvas_mucous2___1)
describeFactors(EMA_MPA$pvas_mucous2___1)
describeFactors(EMA_UNC$pvas_mucous2___1)

mod_mouth_EMA <- glm(as.factor(EMAclass) ~ pvas_mucous2___1 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_mouth_EMA)
exp(coef(mod_mouth_EMA))
exp(confint(mod_mouth_EMA))

#Genital ulcers
describeFactors(EMA_GPA$pvas_mucous2___2)
describeFactors(EMA_MPA$pvas_mucous2___2)
describeFactors(EMA_UNC$pvas_mucous2___2)

mod_genital_EMA <- glm(as.factor(EMAclass) ~ pvas_mucous2___2 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_genital_EMA)
exp(coef(mod_genital_EMA))
exp(confint(mod_genital_EMA))

#Adnexal inflammation (sialadenitis/ dEMAyocystitis)
describeFactors(EMA_GPA$pvas_mucous2___3)
describeFactors(EMA_MPA$pvas_mucous2___3)
describeFactors(EMA_UNC$pvas_mucous2___3)

mod_adn_EMA <- glm(as.factor(EMAclass) ~ pvas_mucous2___3 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_adn_EMA)
exp(coef(mod_adn_EMA))
exp(confint(mod_adn_EMA))

#Significant proptosis
describeFactors(EMA_GPA$pvas_mucous2___4)
describeFactors(EMA_MPA$pvas_mucous2___4)
describeFactors(EMA_UNC$pvas_mucous2___4)

mod_prop_EMA <- glm(as.factor(EMAclass) ~ pvas_mucous2___4 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_prop_EMA)
exp(coef(mod_prop_EMA))
exp(confint(mod_prop_EMA))

#Red eye (Epi)scleritis
describeFactors(EMA_GPA$pvas_mucous2___5)
describeFactors(EMA_MPA$pvas_mucous2___5)
describeFactors(EMA_UNC$pvas_mucous2___5)

mod_redeyeepi_EMA <- glm(as.factor(EMAclass) ~ pvas_mucous2___5 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_redeyeepi_EMA)
exp(coef(mod_redeyeepi_EMA))
exp(confint(mod_redeyeepi_EMA))

#Red eye conjunctivitis/ blepharitis/keratitis
describeFactors(EMA_GPA$pvas_mucous2___6)
describeFactors(EMA_MPA$pvas_mucous2___6)
describeFactors(EMA_UNC$pvas_mucous2___6)

mod_redeyeconj_EMA <- glm(as.factor(EMAclass) ~ pvas_mucous2___6 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_redeyeconj_EMA)
exp(coef(mod_redeyeconj_EMA))
exp(confint(mod_redeyeconj_EMA))

#Uveitis
describeFactors(EMA_GPA$pvas_mucous2___7)
describeFactors(EMA_MPA$pvas_mucous2___7)
describeFactors(EMA_UNC$pvas_mucous2___7)

mod_uve_EMA <- glm(as.factor(EMAclass) ~ pvas_mucous2___7 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_uve_EMA)
exp(coef(mod_uve_EMA))
exp(confint(mod_uve_EMA))

#Blurred vision
describeFactors(EMA_GPA$pvas_mucous2___8)
describeFactors(EMA_MPA$pvas_mucous2___8)
describeFactors(EMA_UNC$pvas_mucous2___8)

mod_bv_EMA <- glm(as.factor(EMAclass) ~ pvas_mucous2___8 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_bv_EMA)
exp(coef(mod_bv_EMA))
exp(confint(mod_bv_EMA))

#Sudden visual loss
describeFactors(EMA_GPA$pvas_mucous2___9)
describeFactors(EMA_MPA$pvas_mucous2___9)
describeFactors(EMA_UNC$pvas_mucous2___9)

mod_svl_EMA <- glm(as.factor(EMAclass) ~ pvas_mucous2___9 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_svl_EMA)
exp(coef(mod_svl_EMA))
exp(confint(mod_svl_EMA))

#Retinal vasculitis/retinal vessel thrombosis/retinal exudates/haemorrhages
describeFactors(EMA_GPA$pvas_mucous2___10)
describeFactors(EMA_MPA$pvas_mucous2___10)
describeFactors(EMA_UNC$pvas_mucous2___10)

mod_rv_EMA <- glm(as.factor(EMAclass) ~ pvas_mucous2___10 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_rv_EMA)
exp(coef(mod_rv_EMA))
exp(confint(mod_rv_EMA))

describeFactors(EMA_GPA$ent)
describeFactors(EMA_MPA$ent)
describeFactors(EMA_UNC$ent)

mod_ent_EMA <- glm(as.factor(EMAclass) ~ ent + Age..years. + gender, data = EMA, family="binomial")
summary(mod_ent_EMA)
exp(coef(mod_ent_EMA))
exp(confint(mod_ent_EMA))

#Bloody nasal discharge/crusts/ulcers/granuloma
describeFactors(EMA_GPA$pvas_ent___1)
describeFactors(EMA_MPA$pvas_ent___1)
describeFactors(EMA_UNC$pvas_ent___1)

mod_bnd_EMA <- glm(as.factor(EMAclass) ~ pvas_ent___1 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_bnd_EMA)
exp(coef(mod_bnd_EMA))
exp(confint(mod_bnd_EMA))

#Paranasal sinus involvement
describeFactors(EMA_GPA$pvas_ent___2)
describeFactors(EMA_MPA$pvas_ent___2)
describeFactors(EMA_UNC$pvas_ent___2)

mod_psi_EMA <- glm(as.factor(EMAclass) ~ pvas_ent___2 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_psi_EMA)
exp(coef(mod_psi_EMA))
exp(confint(mod_psi_EMA))

#Subglottic stenosis/ hoarseness /stridor
describeFactors(EMA_GPA$pvas_ent___3)
describeFactors(EMA_MPA$pvas_ent___3)
describeFactors(EMA_UNC$pvas_ent___3)

mod_sgs_EMA <- glm(as.factor(EMAclass) ~ pvas_ent___3 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_sgs_EMA)
exp(coef(mod_sgs_EMA))
exp(confint(mod_sgs_EMA))

#Conductive hearing loss
describeFactors(EMA_GPA$pvas_ent___4)
describeFactors(EMA_MPA$pvas_ent___4)
describeFactors(EMA_UNC$pvas_ent___4)

mod_chl_EMA <- glm(as.factor(EMAclass) ~ pvas_ent___4 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_chl_EMA)
exp(coef(mod_chl_EMA))
exp(confint(mod_chl_EMA))

#Sensorineural hearing loss
describeFactors(EMA_GPA$pvas_ent___5)
describeFactors(EMA_MPA$pvas_ent___5)
describeFactors(EMA_UNC$pvas_ent___5)

mod_shl_EMA <- glm(as.factor(EMAclass) ~ pvas_ent___5 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_shl_EMA)
exp(coef(mod_shl_EMA))
exp(confint(mod_shl_EMA))


describeFactors(EMA_GPA$chest)
describeFactors(EMA_MPA$chest)
describeFactors(EMA_UNC$chest)

mod_chest_EMA <- glm(as.factor(EMAclass) ~ chest + Age..years. + gender, data = EMA, family="binomial")
summary(mod_chest_EMA)
exp(coef(mod_chest_EMA))
exp(confint(mod_chest_EMA))

#Wheeze or expiratory dyspnea
describeFactors(EMA_GPA$pvas_chest___1)
describeFactors(EMA_MPA$pvas_chest___1)
describeFactors(EMA_UNC$pvas_chest___1)

mod_wheeze_EMA <- glm(as.factor(EMAclass) ~ pvas_chest___1 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_wheeze_EMA)
exp(coef(mod_wheeze_EMA))
exp(confint(mod_wheeze_EMA))

#Endobronchial/endotracheal involvement
describeFactors(EMA_GPA$pvas_chest___2)
describeFactors(EMA_MPA$pvas_chest___2)
describeFactors(EMA_UNC$pvas_chest___2)

mod_endo_EMA <- glm(as.factor(EMAclass) ~ pvas_chest___2 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_endo_EMA)
exp(coef(mod_endo_EMA))
exp(confint(mod_endo_EMA))

#Nodules or cavities 
describeFactors(EMA_GPA$pvas_chest___3)
describeFactors(EMA_MPA$pvas_chest___3)
describeFactors(EMA_UNC$pvas_chest___3)

mod_nodu_EMA <- glm(as.factor(EMAclass) ~ pvas_chest___3 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_nodu_EMA)
exp(coef(mod_nodu_EMA))
exp(confint(mod_nodu_EMA))

#Pleural effusion/pleurisy 
describeFactors(EMA_GPA$pvas_chest___4)
describeFactors(EMA_MPA$pvas_chest___4)
describeFactors(EMA_UNC$pvas_chest___4)

mod_pe_EMA <- glm(as.factor(EMAclass) ~ pvas_chest___4 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_pe_EMA)
exp(coef(mod_pe_EMA))
exp(confint(mod_pe_EMA))

#Infiltrate
describeFactors(EMA_GPA$pvas_chest___5)
describeFactors(EMA_MPA$pvas_chest___5)
describeFactors(EMA_UNC$pvas_chest___5)

mod_inf_EMA <- glm(as.factor(EMAclass) ~ pvas_chest___5 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_inf_EMA)
exp(coef(mod_inf_EMA))
exp(confint(mod_inf_EMA))

#Massive haemoptysis/Alveolar haemorrhage
describeFactors(EMA_GPA$pvas_chest___6)
describeFactors(EMA_MPA$pvas_chest___6)
describeFactors(EMA_UNC$pvas_chest___6)

mod_mahaem_EMA <- glm(as.factor(EMAclass) ~ pvas_chest___6 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_mahaem_EMA)
exp(coef(mod_mahaem_EMA))
exp(confint(mod_mahaem_EMA))

#Respiratory failure
describeFactors(EMA_GPA$pvas_chest___7)
describeFactors(EMA_MPA$pvas_chest___7)
describeFactors(EMA_UNC$pvas_chest___7)

mod_rf_EMA <- glm(as.factor(EMAclass) ~ pvas_chest___7 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_rf_EMA)
exp(coef(mod_rf_EMA))
exp(confint(mod_rf_EMA))

describeFactors(EMA_GPA$cardio)
describeFactors(EMA_MPA$cardio)
describeFactors(EMA_UNC$cardio)

mod_cardio_EMA <- glm(as.factor(EMAclass) ~ cardio + Age..years. + gender, data = EMA, family="binomial")
summary(mod_cardio_EMA)
exp(coef(mod_cardio_EMA))
exp(confint(mod_cardio_EMA))

#Loss of pulses
describeFactors(EMA_GPA$pvas_cv___1)
describeFactors(EMA_MPA$pvas_cv___1)
describeFactors(EMA_UNC$pvas_cv___1)

mod_lossp_EMA <- glm(as.factor(EMAclass) ~ pvas_cv___1 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_lossp_EMA)
exp(coef(mod_lossp_EMA))
exp(confint(mod_lossp_EMA))

#Bruits over accessible arteries
describeFactors(EMA_GPA$pvas_cv___2)
describeFactors(EMA_MPA$pvas_cv___2)
describeFactors(EMA_UNC$pvas_cv___2)

mod_bruits_EMA <- glm(as.factor(EMAclass) ~ pvas_cv___2 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_bruits_EMA)
exp(coef(mod_bruits_EMA))
exp(confint(mod_bruits_EMA))

#Blood pressure discrepancy
describeFactors(EMA_GPA$pvas_cv___3)
describeFactors(EMA_MPA$pvas_cv___3)
describeFactors(EMA_UNC$pvas_cv___3)

mod_bp_EMA <- glm(as.factor(EMAclass) ~ pvas_cv___3 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_bp_EMA)
exp(coef(mod_bp_EMA))
exp(confint(mod_bp_EMA))

#Claudication of extremities
describeFactors(EMA_GPA$pvas_cv___4)
describeFactors(EMA_MPA$pvas_cv___4)
describeFactors(EMA_UNC$pvas_cv___4)

mod_claud_EMA <- glm(as.factor(EMAclass) ~ pvas_cv___4 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_claud_EMA)
exp(coef(mod_claud_EMA))
exp(confint(mod_claud_EMA))

#Ischaemic cardiac pain
describeFactors(EMA_GPA$pvas_cv___5)
describeFactors(EMA_MPA$pvas_cv___5)
describeFactors(EMA_UNC$pvas_cv___5)

mod_icp_EMA <- glm(as.factor(EMAclass) ~ pvas_cv___5 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_icp_EMA)
exp(coef(mod_icp_EMA))
exp(confint(mod_icp_EMA))

#Cardiomyopathy
describeFactors(EMA_GPA$pvas_cv___6)
describeFactors(EMA_MPA$pvas_cv___6)
describeFactors(EMA_UNC$pvas_cv___6)

mod_cmpy_EMA <- glm(as.factor(EMAclass) ~ pvas_cv___6 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_cmpy_EMA)
exp(coef(mod_cmpy_EMA))
exp(confint(mod_cmpy_EMA))

#Congestive cardiac failure
describeFactors(EMA_GPA$pvas_cv___7)
describeFactors(EMA_MPA$pvas_cv___7)
describeFactors(EMA_UNC$pvas_cv___7)

mod_ccf_EMA <- glm(as.factor(EMAclass) ~ pvas_cv___7 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_ccf_EMA)
exp(coef(mod_ccf_EMA))
exp(confint(mod_ccf_EMA))

#Valvular heart disease
describeFactors(EMA_GPA$pvas_cv___8)
describeFactors(EMA_MPA$pvas_cv___8)
describeFactors(EMA_UNC$pvas_cv___8)

mod_vhd_EMA <- glm(as.factor(EMAclass) ~ pvas_cv___8 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_vhd_EMA)
exp(coef(mod_vhd_EMA))
exp(confint(mod_vhd_EMA))

#Pericarditis
describeFactors(EMA_GPA$pvas_cv___9)
describeFactors(EMA_MPA$pvas_cv___9)
describeFactors(EMA_UNC$pvas_cv___9)

mod_pericard_EMA <- glm(as.factor(EMAclass) ~ pvas_cv___9 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_pericard_EMA)
exp(coef(mod_pericard_EMA))
exp(confint(mod_pericard_EMA))

describeFactors(EMA_GPA$gi)
describeFactors(EMA_MPA$gi)
describeFactors(EMA_UNC$gi)

mod_gi_EMA <- glm(as.factor(EMAclass) ~ gi + Age..years. + gender, data = EMA, family="binomial")
summary(mod_gi_EMA)
exp(coef(mod_gi_EMA))
exp(confint(mod_gi_EMA))

#Abdominal pain
describeFactors(EMA_GPA$pvas_gi___1)
describeFactors(EMA_MPA$pvas_gi___1)
describeFactors(EMA_UNC$pvas_gi___1)

mod_abdp_EMA <- glm(as.factor(EMAclass) ~ pvas_gi___1 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_abdp_EMA)
exp(coef(mod_abdp_EMA))
exp(confint(mod_abdp_EMA))

#Peritonitis
describeFactors(EMA_GPA$pvas_gi___2)
describeFactors(EMA_MPA$pvas_gi___2)
describeFactors(EMA_UNC$pvas_gi___2)

mod_periton_EMA <- glm(as.factor(EMAclass) ~ pvas_gi___2 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_periton_EMA)
exp(coef(mod_periton_EMA))
exp(confint(mod_periton_EMA))

#Blood in stools or bloody diarrhea
describeFactors(EMA_GPA$pvas_gi___3)
describeFactors(EMA_MPA$pvas_gi___3)
describeFactors(EMA_UNC$pvas_gi___3)

mod_bloodstool_EMA <- glm(as.factor(EMAclass) ~ pvas_gi___3 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_bloodstool_EMA)
exp(coef(mod_bloodstool_EMA))
exp(confint(mod_bloodstool_EMA))

#Bowel ischaemia
describeFactors(EMA_GPA$pvas_gi___4)
describeFactors(EMA_MPA$pvas_gi___4)
describeFactors(EMA_UNC$pvas_gi___4)

mod_bowelisch_EMA <- glm(as.factor(EMAclass) ~ pvas_gi___4 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_bowelisch_EMA)
exp(coef(mod_bowelisch_EMA))
exp(confint(mod_bowelisch_EMA))

describeFactors(EMA_GPA$ren)
describeFactors(EMA_MPA$ren)
describeFactors(EMA_UNC$ren)

mod_ren_EMA <- glm(as.factor(EMAclass) ~ ren + Age..years. + gender, data = EMA, family="binomial")
summary(mod_ren_EMA)
exp(coef(mod_ren_EMA))
exp(confint(mod_ren_EMA))

#Hypertension >95th centile (for height)
describeFactors(EMA_GPA$pvas_renal___1)
describeFactors(EMA_MPA$pvas_renal___1)
describeFactors(EMA_UNC$pvas_renal___1)

mod_hyt_EMA <- glm(as.factor(EMAclass) ~ pvas_renal___1 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_hyt_EMA)
exp(coef(mod_hyt_EMA))
exp(confint(mod_hyt_EMA))

#Proteinuria >0.3 g/24h,>20mg/mmol creatinine
describeFactors(EMA_GPA$pvas_renal___2)
describeFactors(EMA_MPA$pvas_renal___2)
describeFactors(EMA_UNC$pvas_renal___2)

mod_puria_EMA <- glm(as.factor(EMAclass) ~ pvas_renal___2 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_puria_EMA)
exp(coef(mod_puria_EMA))
exp(confint(mod_puria_EMA))

#Haematuria >= 2+ or 5 rbc/hpf or red cell casts
describeFactors(EMA_GPA$pvas_renal___3)
describeFactors(EMA_MPA$pvas_renal___3)
describeFactors(EMA_UNC$pvas_renal___3)

mod_hemuria_EMA <- glm(as.factor(EMAclass) ~ pvas_renal___3 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_hemuria_EMA)
exp(coef(mod_hemuria_EMA))
exp(confint(mod_hemuria_EMA))

#GFR 50-80ml/min/1.73 m2 *s
describeFactors(EMA_GPA$pvas_renal___4)
describeFactors(EMA_MPA$pvas_renal___4)
describeFactors(EMA_UNC$pvas_renal___4)

mod_GFR1_EMA <- glm(as.factor(EMAclass) ~ pvas_renal___4 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_GFR1_EMA)
exp(coef(mod_GFR1_EMA))
exp(confint(mod_GFR1_EMA))

#GFR 15-49 ml/min/1.73 m2 *
describeFactors(EMA_GPA$pvas_renal___5)
describeFactors(EMA_MPA$pvas_renal___5)
describeFactors(EMA_UNC$pvas_renal___5)

mod_GFR2_EMA <- glm(as.factor(EMAclass) ~ pvas_renal___5 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_GFR2_EMA)
exp(coef(mod_GFR2_EMA))
exp(confint(mod_GFR2_EMA))

#GFR < 15 ml/min/1.73m2 *
describeFactors(EMA_GPA$pvas_renal___6)
describeFactors(EMA_MPA$pvas_renal___6)
describeFactors(EMA_UNC$pvas_renal___6)

mod_GFR3_EMA <- glm(as.factor(EMAclass) ~ pvas_renal___6 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_GFR3_EMA)
exp(coef(mod_GFR3_EMA))
exp(confint(mod_GFR3_EMA))

#Rise in creatinine > 10% or Creatinine clearance (GFR) fall > 25%
describeFactors(EMA_GPA$pvas_renal___8)
describeFactors(EMA_MPA$pvas_renal___8)
describeFactors(EMA_UNC$pvas_renal___8)

mod_creatinine_EMA <- glm(as.factor(EMAclass) ~ pvas_renal___8 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_creatinine_EMA)
exp(coef(mod_creatinine_EMA))
exp(confint(mod_creatinine_EMA))


describeFactors(EMA_GPA$cns)
describeFactors(EMA_MPA$cns)
describeFactors(EMA_UNC$cns)

mod_cns_EMA <- glm(as.factor(EMAclass) ~ cns + Age..years. + gender, data = EMA, family="binomial")
summary(mod_cns_EMA)
exp(coef(mod_cns_EMA))
exp(confint(mod_cns_EMA))

#Headache
describeFactors(EMA_GPA$pvas_cns___1)
describeFactors(EMA_MPA$pvas_cns___1)
describeFactors(EMA_UNC$pvas_cns___1)

mod_headache_EMA <- glm(as.factor(EMAclass) ~ pvas_cns___1 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_headache_EMA)
exp(coef(mod_headache_EMA))
exp(confint(mod_headache_EMA))

#Meningitis/encephalitis
describeFactors(EMA_GPA$pvas_cns___2)
describeFactors(EMA_MPA$pvas_cns___2)
describeFactors(EMA_UNC$pvas_cns___2)

mod_meningitis_EMA <- glm(as.factor(EMAclass) ~ pvas_cns___2 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_meningitis_EMA)
exp(coef(mod_meningitis_EMA))
exp(confint(mod_meningitis_EMA))

#Organic confusion/cognitive dysfunction
describeFactors(EMA_GPA$pvas_cns___3)
describeFactors(EMA_MPA$pvas_cns___3)
describeFactors(EMA_UNC$pvas_cns___3)

mod_orgconf_EMA <- glm(as.factor(EMAclass) ~ pvas_cns___3 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_orgconf_EMA)
exp(coef(mod_orgconf_EMA))
exp(confint(mod_orgconf_EMA))

#Seizures (not hypertensive)
describeFactors(EMA_GPA$pvas_cns___4)
describeFactors(EMA_MPA$pvas_cns___4)
describeFactors(EMA_UNC$pvas_cns___4)

mod_seizures_EMA <- glm(as.factor(EMAclass) ~ pvas_cns___4 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_seizures_EMA)
exp(coef(mod_seizures_EMA))
exp(confint(mod_seizures_EMA))

#Stroke
describeFactors(EMA_GPA$pvas_cns___5)
describeFactors(EMA_MPA$pvas_cns___5)
describeFactors(EMA_UNC$pvas_cns___5)

mod_stroke_EMA <- glm(as.factor(EMAclass) ~ pvas_cns___5 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_stroke_EMA)
exp(coef(mod_stroke_EMA))
exp(confint(mod_stroke_EMA))

#Cord lesion
describeFactors(EMA_GPA$pvas_cns___6)
describeFactors(EMA_MPA$pvas_cns___6)
describeFactors(EMA_UNC$pvas_cns___6)

mod_cl_EMA <- glm(as.factor(EMAclass) ~ pvas_cns___6 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_cl_EMA)
exp(coef(mod_cl_EMA))
exp(confint(mod_cl_EMA))

#Cranial nerve palsy
describeFactors(EMA_GPA$pvas_cns___7)
describeFactors(EMA_MPA$pvas_cns___7)
describeFactors(EMA_UNC$pvas_cns___7)

mod_cnp_EMA <- glm(as.factor(EMAclass) ~ pvas_cns___7 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_cnp_EMA)
exp(coef(mod_cnp_EMA))
exp(confint(mod_cnp_EMA))

#Sensory peripheral neuropathy
describeFactors(EMA_GPA$pvas_cns___8)
describeFactors(EMA_MPA$pvas_cns___8)
describeFactors(EMA_UNC$pvas_cns___8)

mod_spn_EMA <- glm(as.factor(EMAclass) ~ pvas_cns___8 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_spn_EMA)
exp(coef(mod_spn_EMA))
exp(confint(mod_spn_EMA))

#Motor mononeuritis multiplex
describeFactors(EMA_GPA$pvas_cns___9)
describeFactors(EMA_MPA$pvas_cns___9)
describeFactors(EMA_UNC$pvas_cns___9)

mod_mmm_EMA <- glm(as.factor(EMAclass) ~ pvas_cns___9 + Age..years. + gender, data = EMA, family="binomial")
summary(mod_mmm_EMA)
exp(coef(mod_mmm_EMA))
exp(confint(mod_mmm_EMA))

#### SUPP Table 1a ####
#GENERAL
gen_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$gen == 1),])
gen_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$gen == 0),])
gen_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$gen == 1),])
gen_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$gen == 0),])

df_gen <- data.frame(
  gpa = c(gen_gpa_pres_acr,gen_gpa_abs_acr),
  mpa = c(gen_mpa_pres_acr,gen_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_gen, conf.int = TRUE)

#Myalgia
mya_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gen1___1 == 1),])
mya_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gen1___1 == 0),])
mya_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gen1___1 == 1),])
mya_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gen1___1 == 0),])

df_mya <- data.frame(
  gpa = c(mya_gpa_pres_acr,mya_gpa_abs_acr),
  mpa = c(mya_mpa_pres_acr,mya_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_mya, conf.int = TRUE)

#Arthralgia or arthritis
arth_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gen1___2 == 1),])
arth_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gen1___2 == 0),])
arth_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gen1___2 == 1),])
arth_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gen1___2 == 0),])

df_arth <- data.frame(
  gpa = c(arth_gpa_pres_acr,arth_gpa_abs_acr),
  mpa = c(arth_mpa_pres_acr,arth_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_arth, conf.int = TRUE)

#Fever >= 38 deg C
fev_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gen1___3 == 1),])
fev_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gen1___3 == 0),])
fev_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gen1___3 == 1),])
fev_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gen1___3 == 0),])

df_fev <- data.frame(
  gpa = c(fev_gpa_pres_acr,fev_gpa_abs_acr),
  mpa = c(fev_mpa_pres_acr,fev_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_fev, conf.int = TRUE)

#Weight loss >= 5% body weight
wl_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gen1___4 == 1),])
wl_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gen1___4 == 0),])
wl_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gen1___4 == 1),])
wl_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gen1___4 == 0),])

df_wl <- data.frame(
  gpa = c(wl_gpa_pres_acr,wl_gpa_abs_acr),
  mpa = c(wl_mpa_pres_acr,wl_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_wl, conf.int = TRUE)

#Cutaneous involvement
cuta_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$cuta == 1),])
cuta_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$cuta == 0),])
cuta_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$cuta == 1),])
cuta_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$cuta == 0),])

df_cuta <- data.frame(
  gpa = c(cuta_gpa_pres_acr,cuta_gpa_abs_acr),
  mpa = c(cuta_mpa_pres_acr,cuta_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_cuta, conf.int = TRUE)

#Polymorphous exanthema
pe_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___1 == 1),])
pe_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___1 == 0),])
pe_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___1 == 1),])
pe_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___1 == 0),])

df_pe <- data.frame(
  gpa = c(pe_gpa_pres_acr,pe_gpa_abs_acr),
  mpa = c(pe_mpa_pres_acr,pe_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_pe, conf.int = TRUE)

#Livedo
liv_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___2 == 1),])
liv_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___2 == 0),])
liv_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___2 == 1),])
liv_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___2 == 0),])

df_liv <- data.frame(
  gpa = c(liv_gpa_pres_acr,liv_gpa_abs_acr),
  mpa = c(liv_mpa_pres_acr,liv_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_liv, conf.int = TRUE)

#Panniculitis
pan_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___3 == 1),])
pan_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___3 == 0),])
pan_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___3 == 1),])
pan_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___3 == 0),])

df_pan <- data.frame(
  gpa = c(pan_gpa_pres_acr,pan_gpa_abs_acr),
  mpa = c(pan_mpa_pres_acr,pan_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_pan, conf.int = TRUE)

#Purpura
purp_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___4 == 1),])
purp_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___4 == 0),])
purp_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___4 == 1),])
purp_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___4 == 0),])

df_purp <- data.frame(
  gpa = c(purp_gpa_pres_acr,purp_gpa_pres_acr),
  mpa = c(purp_mpa_pres_acr,purp_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_purp, conf.int = TRUE)

#Skin nodules
sn_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___5 == 1),])
sn_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___5 == 0),])
sn_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___5 == 1),])
sn_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___5 == 0),])

df_sn <- data.frame(
  gpa = c(sn_gpa_pres_acr,sn_gpa_abs_acr),
  mpa = c(sn_mpa_pres_acr,sn_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_sn, conf.int = TRUE)

#Infarct (nail edge lesion, splinter haemorrhage)
inf_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___6 == 1),])
inf_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___6 == 0),])
inf_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___6 == 1),])
inf_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___6 == 0),])

df_inf <- data.frame(
  gpa = c(inf_gpa_pres_acr,inf_gpa_abs_acr),
  mpa = c(inf_mpa_pres_acr,inf_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_inf, conf.int = TRUE)

#Ulcer (full-thickness necrosis)
ulc_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___7 == 1),])
ulc_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___7 == 0),])
ulc_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___7 == 1),])
ulc_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___7 == 0),])

df_ulc <- data.frame(
  gpa = c(ulc_gpa_pres_acr,ulc_gpa_abs_acr),
  mpa = c(ulc_mpa_pres_acr,ulc_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_ulc, conf.int = TRUE)

#Gangrene (extensive necrosis)
gang_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___8 == 1),])
gang_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___8 == 0),])
gang_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___8 == 1),])
gang_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___8 == 0),])

df_gang <- data.frame(
  gpa = c(gang_gpa_pres_acr,gang_gpa_abs_acr),
  mpa = c(gang_mpa_pres_acr,gang_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_gang, conf.int = TRUE)

#Other skin vasculitis (specify below)
osv_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___9 == 1),])
osv_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cut2___9 == 0),])
osv_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___9 == 1),])
osv_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cut2___9 == 0),])

df_osv <- data.frame(
  gpa = c(osv_gpa_pres_acr,osv_gpa_abs_acr),
  mpa = c(osv_mpa_pres_acr,osv_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_osv, conf.int = TRUE)

#Mucous membranes/eyes involvement
muc_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$muc == 1),])
muc_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$muc == 0),])
muc_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$muc == 1),])
muc_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$muc == 0),])

df_muc <- data.frame(
  gpa = c(muc_gpa_pres_acr,muc_gpa_abs_acr),
  mpa = c(muc_mpa_pres_acr,muc_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_muc, conf.int = TRUE)

#Mouth ulcers/granulomata
mug_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___1 == 1),])
mug_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___1 == 0),])
mug_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___1 == 1),])
mug_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___1 == 0),])

df_mug <- data.frame(
  gpa = c(mug_gpa_pres_acr,mug_gpa_abs_acr),
  mpa = c(mug_mpa_pres_acr,mug_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_mug, conf.int = TRUE)

#Genital ulcers
gulc_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___2 == 1),])
gulc_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___2 == 0),])
gulc_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___2 == 1),])
gulc_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___2 == 0),])

df_gulc <- data.frame(
  gpa = c(gulc_gpa_pres_acr,gulc_gpa_abs_acr),
  mpa = c(gulc_mpa_pres_acr,gulc_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_gulc, conf.int = TRUE)

#Adnexal inflammation (sialadenitis/ dEMAyocystitis)
adn_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___3 == 1),])
adn_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___3 == 0),])
adn_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___3 == 1),])
adn_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___3 == 0),])

df_adn <- data.frame(
  gpa = c(adn_gpa_pres_acr,adn_gpa_abs_acr),
  mpa = c(adn_mpa_pres_acr,adn_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_adn, conf.int = TRUE)

#Significant proptosis
prop_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___4 == 1),])
prop_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___4 == 0),])
prop_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___4 == 1),])
prop_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___4 == 0),])

df_prop <- data.frame(
  gpa = c(prop_gpa_pres_acr,prop_gpa_abs_acr),
  mpa = c(prop_mpa_pres_acr,prop_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_prop, conf.int = TRUE)

#Red eye (Epi)scleritis
red_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___5 == 1),])
red_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___5 == 0),])
red_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___5 == 1),])
red_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___5 == 0),])

df_red <- data.frame(
  gpa = c(red_gpa_pres_acr,red_gpa_abs_acr),
  mpa = c(red_mpa_pres_acr,red_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_red, conf.int = TRUE)

#Red eye conjunctivitis/ blepharitis/keratitis
redcon_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___6 == 1),])
redcon_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___6 == 0),])
redcon_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___6 == 1),])
redcon_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___6 == 0),])

df_redcon <- data.frame(
  gpa = c(redcon_gpa_pres_acr,redcon_gpa_abs_acr),
  mpa = c(redcon_mpa_pres_acr,redcon_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_redcon, conf.int = TRUE)

#Uveitis
uve_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___7 == 1),])
uve_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___7 == 0),])
uve_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___7 == 1),])
uve_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___7 == 0),])

df_uve <- data.frame(
  gpa = c(uve_gpa_pres_acr,uve_gpa_abs_acr),
  mpa = c(uve_mpa_pres_acr,uve_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_uve, conf.int = TRUE)

#Blurred vision
bv_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___8 == 1),])
bv_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___8 == 0),])
bv_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___8 == 1),])
bv_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___8 == 0),])

df_bv <- data.frame(
  gpa = c(bv_gpa_pres_acr,bv_gpa_abs_acr),
  mpa = c(bv_mpa_pres_acr,bv_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_bv, conf.int = TRUE)

#Sudden visual loss
svl_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___9 == 1),])
svl_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___9 == 0),])
svl_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___9 == 1),])
svl_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___9 == 0),])

df_svl <- data.frame(
  gpa = c(svl_gpa_pres_acr,svl_gpa_abs_acr),
  mpa = c(svl_mpa_pres_acr,svl_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_svl, conf.int = TRUE)

#Retinal vasculitis/retinal vessel thrombosis/retinal exudates/haemorrhages
rvas_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___10 == 1),])
rvas_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_mucous2___10 == 0),])
rvas_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___10 == 1),])
rvas_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_mucous2___10 == 0),])

df_rvas <- data.frame(
  gpa = c(rvas_gpa_pres_acr,rvas_gpa_abs_acr),
  mpa = c(rvas_mpa_pres_acr,rvas_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_rvas, conf.int = TRUE)

#ENT
ent_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$ent == 1),])
ent_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$ent == 0),])
ent_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$ent == 1),])
ent_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$ent == 0),])

df_ent <- data.frame(
  gpa = c(ent_gpa_pres_acr,ent_gpa_abs_acr),
  mpa = c(ent_mpa_pres_acr,ent_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_ent, conf.int = TRUE)

#Bloody nasal discharge/crusts/ulcers/granuloma
bnd_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_ent___1 == 1),])
bnd_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_ent___1 == 0),])
bnd_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_ent___1 == 1),])
bnd_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_ent___1 == 0),])

df_bnd <- data.frame(
  gpa = c(bnd_gpa_pres_acr,bnd_gpa_abs_acr),
  mpa = c(bnd_mpa_pres_acr,bnd_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_bnd, conf.int = TRUE)

#Paranasal sinus involvement
psi_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_ent___2 == 1),])
psi_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_ent___2 == 0),])
psi_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_ent___2 == 1),])
psi_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_ent___2 == 0),])

df_psi <- data.frame(
  gpa = c(psi_gpa_pres_acr,psi_gpa_abs_acr),
  mpa = c(psi_mpa_pres_acr,psi_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_psi, conf.int = TRUE)

#Subglottic stenosis/ hoarseness /stridor
sgs_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_ent___3 == 1),])
sgs_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_ent___3 == 0),])
sgs_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_ent___3 == 1),])
sgs_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_ent___3 == 0),])

df_sgs <- data.frame(
  gpa = c(sgs_gpa_pres_acr,sgs_gpa_abs_acr),
  mpa = c(sgs_mpa_pres_acr,sgs_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_sgs, conf.int = TRUE)

#Conductive hearing loss
chl_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_ent___4 == 1),])
chl_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_ent___4 == 0),])
chl_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_ent___4 == 1),])
chl_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_ent___4 == 0),])

df_chl <- data.frame(
  gpa = c(chl_gpa_pres_acr,chl_gpa_abs_acr),
  mpa = c(chl_mpa_pres_acr,chl_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_chl, conf.int = TRUE)

#Sensorineural hearing loss
shl_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_ent___5 == 1),])
shl_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_ent___5 == 0),])
shl_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_ent___5 == 1),])
shl_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_ent___5 == 0),])

df_shl <- data.frame(
  gpa = c(shl_gpa_pres_acr,shl_gpa_abs_acr),
  mpa = c(shl_mpa_pres_acr,shl_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_shl, conf.int = TRUE)

#Chest Involvement
chest_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$chest == 1),])
chest_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$chest == 0),])
chest_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$chest == 1),])
chest_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$chest == 0),])

df_chest <- data.frame(
  gpa = c(chest_gpa_pres_acr,chest_gpa_abs_acr),
  mpa = c(chest_mpa_pres_acr,chest_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_chest, conf.int = TRUE)

#Wheeze or expiratory dyspnea
wheeze_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_chest___1 == 1),])
wheeze_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_chest___1 == 0),])
wheeze_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_chest___1 == 1),])
wheeze_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_chest___1 == 0),])

df_wheeze <- data.frame(
  gpa = c(wheeze_gpa_pres_acr,wheeze_gpa_abs_acr),
  mpa = c(wheeze_mpa_pres_acr,wheeze_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_wheeze, conf.int = TRUE)

#Endobronchial/endotracheal involvement
ebi_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_chest___2 == 1),])
ebi_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_chest___2 == 0),])
ebi_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_chest___2 == 1),])
ebi_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_chest___2 == 0),])

df_ebi <- data.frame(
  gpa = c(ebi_gpa_pres_acr,ebi_gpa_abs_acr),
  mpa = c(ebi_mpa_pres_acr,ebi_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_ebi, conf.int = TRUE)

#Nodules or cavities 
nod_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_chest___3 == 1),])
nod_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_chest___3 == 0),])
nod_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_chest___3 == 1),])
nod_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_chest___3 == 0),])

df_nod <- data.frame(
  gpa = c(nod_gpa_pres_acr,nod_gpa_abs_acr),
  mpa = c(nod_mpa_pres_acr,nod_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_nod, conf.int = TRUE)

#Pleural effusion/pleurisy
peff_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_chest___4 == 1),])
peff_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_chest___4 == 0),])
peff_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_chest___4 == 1),])
peff_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_chest___4 == 0),])

df_peff <- data.frame(
  gpa = c(peff_gpa_pres_acr,peff_gpa_abs_acr),
  mpa = c(peff_mpa_pres_acr,peff_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_peff, conf.int = TRUE)

#Infiltrate
inf_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_chest___5 == 1),])
inf_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_chest___5 == 0),])
inf_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_chest___5 == 1),])
inf_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_chest___5 == 0),])

df_inf <- data.frame(
  gpa = c(inf_gpa_pres_acr,inf_gpa_abs_acr),
  mpa = c(inf_mpa_pres_acr,inf_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_inf, conf.int = TRUE)

#Massive haemoptysis/Alveolar haemorrhage
mhaem_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_chest___6 == 1),])
mhaem_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_chest___6 == 0),])
mhaem_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_chest___6 == 1),])
mhaem_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_chest___6 == 0),])

df_mhaem <- data.frame(
  gpa = c(mhaem_gpa_pres_acr,mhaem_gpa_abs_acr),
  mpa = c(mhaem_mpa_pres_acr,mhaem_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_mhaem, conf.int = TRUE)

#Respiratory failure
respf_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_chest___7 == 1),])
respf_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_chest___7 == 0),])
respf_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_chest___7 == 1),])
respf_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_chest___7 == 0),])

df_respf <- data.frame(
  gpa = c(respf_gpa_pres_acr,respf_gpa_abs_acr),
  mpa = c(respf_mpa_pres_acr,respf_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_respf, conf.int = TRUE)

#Cardiac Involvement
cardio_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$cardio == 1),])
cardio_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$cardio == 0),])
cardio_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$cardio == 1),])
cardio_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$cardio == 0),])

df_cardio <- data.frame(
  gpa = c(cardio_gpa_pres_acr,cardio_gpa_abs_acr),
  mpa = c(cardio_mpa_pres_acr,cardio_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_cardio, conf.int = TRUE)

#Loss of pulses
lop_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___1 == 1),])
lop_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___1 == 0),])
lop_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___1 == 1),])
lop_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___1 == 0),])

df_lop <- data.frame(
  gpa = c(lop_gpa_pres_acr,lop_gpa_abs_acr),
  mpa = c(lop_mpa_pres_acr,lop_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_lop, conf.int = TRUE)

#Bruits over accessible arteries
boaa_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___2 == 1),])
boaa_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___2 == 0),])
boaa_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___2 == 1),])
boaa_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___2 == 0),])

df_boaa <- data.frame(
  gpa = c(boaa_gpa_pres_acr,boaa_gpa_abs_acr),
  mpa = c(boaa_mpa_pres_acr,boaa_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_boaa, conf.int = TRUE)

#Blood pressure discrepancy
bpd_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___3 == 1),])
bpd_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___3 == 0),])
bpd_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___3 == 1),])
bpd_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___3 == 0),])

df_bpd <- data.frame(
  gpa = c(bpd_gpa_pres_acr,bpd_gpa_abs_acr),
  mpa = c(bpd_mpa_pres_acr,bpd_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_bpd, conf.int = TRUE)

#Claudication of extremities
coe_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___4 == 1),])
coe_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___4 == 0),])
coe_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___4 == 1),])
coe_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___4 == 0),])

df_coe <- data.frame(
  gpa = c(coe_gpa_pres_acr,coe_gpa_abs_acr),
  mpa = c(coe_mpa_pres_acr,coe_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_coe, conf.int = TRUE)

#Ischaemic cardiac pain
icp_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___5 == 1),])
icp_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___5 == 0),])
icp_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___5 == 1),])
icp_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___5 == 0),])

df_icp <- data.frame(
  gpa = c(icp_gpa_pres_acr,icp_gpa_abs_acr),
  mpa = c(icp_mpa_pres_acr,icp_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_icp, conf.int = TRUE)

#Cardiomyopathy
cmpy_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___6 == 1),])
cmpy_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___6 == 0),])
cmpy_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___6 == 1),])
cmpy_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___6 == 0),])

df_cmpy <- data.frame(
  gpa = c(cmpy_gpa_pres_acr,cmpy_gpa_abs_acr),
  mpa = c(cmpy_mpa_pres_acr,cmpy_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_cmpy, conf.int = TRUE)

#Congestive cardiac failure
ccf_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___7 == 1),])
ccf_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___7 == 0),])
ccf_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___7 == 1),])
ccf_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___7 == 0),])

df_ccf <- data.frame(
  gpa = c(ccf_gpa_pres_acr,ccf_gpa_abs_acr),
  mpa = c(ccf_mpa_pres_acr,ccf_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_ccf, conf.int = TRUE)

#Valvular heart disease
vhd_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___8 == 1),])
vhd_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___8 == 0),])
vhd_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___8 == 1),])
vhd_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___8 == 0),])

df_vhd <- data.frame(
  gpa = c(vhd_gpa_pres_acr,vhd_gpa_abs_acr),
  mpa = c(vhd_mpa_pres_acr,vhd_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_vhd, conf.int = TRUE)

#Pericarditis
pcd_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___9 == 1),])
pcd_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cv___9 == 0),])
pcd_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___9 == 1),])
pcd_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cv___9 == 0),])

df_pcd <- data.frame(
  gpa = c(pcd_gpa_pres_acr,pcd_gpa_abs_acr),
  mpa = c(pcd_mpa_pres_acr,pcd_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_pcd, conf.int = TRUE)

#Gastrointestinal Involvement
gi_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$gi == 1),])
gi_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$gi == 0),])
gi_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$gi == 1),])
gi_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$gi == 0),])

df_gi <- data.frame(
  gpa = c(gi_gpa_pres_acr,gi_gpa_abs_acr),
  mpa = c(gi_mpa_pres_acr,gi_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_gi, conf.int = TRUE)

#Abdominal pain
ap_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gi___1 == 1),])
ap_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gi___1 == 0),])
ap_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gi___1 == 1),])
ap_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$gpvas_gi___1 == 0),])

df_ap <- data.frame(
  gpa = c(ap_gpa_pres_acr,ap_gpa_abs_acr),
  mpa = c(ap_mpa_pres_acr,ap_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_ap, conf.int = TRUE)

#Peritonitis
pos_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gi___2 == 1),])
pos_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gi___2 == 0),])
pos_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gi___2 == 1),])
pos_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gi___2 == 0),])

df_pos <- data.frame(
  gpa = c(pos_gpa_pres_acr,pos_gpa_abs_acr),
  mpa = c(pos_mpa_pres_acr,pos_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_pos, conf.int = TRUE)

#Blood in stools or bloody diarrhea
bis_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gi___3 == 1),])
bis_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gi___3 == 0),])
bis_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gi___3 == 1),])
bis_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gi___3 == 0),])

df_bis <- data.frame(
  gpa = c(bis_gpa_pres_acr,bis_gpa_abs_acr),
  mpa = c(bis_mpa_pres_acr,bis_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_bis, conf.int = TRUE)

#Bowel ischaemia
bi_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gi___4 == 1),])
bi_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_gi___4 == 0),])
bi_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gi___4 == 1),])
bi_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_gi___4 == 0),])

df_bi <- data.frame(
  gpa = c(bi_gpa_pres_acr,bi_gpa_abs_acr),
  mpa = c(bi_mpa_pres_acr,bi_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_bi, conf.int = TRUE)

#Renal Involvement
ren_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$ren == 1),])
ren_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$ren == 0),])
ren_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$ren == 1),])
ren_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$ren == 0),])

df_ren <- data.frame(
  gpa = c(ren_gpa_pres_acr,ren_gpa_abs_acr),
  mpa = c(ren_mpa_pres_acr,ren_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_ren, conf.int = TRUE)

#Hypertension >95th centile (for height)
hyp_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_renal___1 == 1),])
hyp_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_renal___1 == 0),])
hyp_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_renal___1 == 1),])
hyp_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_renal___1 == 0),])

df_hyp <- data.frame(
  gpa = c(hyp_gpa_pres_acr,hyp_gpa_abs_acr),
  mpa = c(hyp_mpa_pres_acr,hyp_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_hyp, conf.int = TRUE)

#Proteinuria >0.3 g/24h,>20mg/mmol creatinine
pru_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_renal___2 == 1),])
pru_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_renal___2 == 0),])
pru_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_renal___2 == 1),])
pru_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_renal___2 == 0),])

df_pru <- data.frame(
  gpa = c(pru_gpa_pres_acr,pru_gpa_abs_acr),
  mpa = c(pru_mpa_pres_acr,pru_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_pru, conf.int = TRUE)

#Haematuria >= 2+ or 5 rbc/hpf or red cell casts
haem_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_renal___3 == 1),])
haem_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_renal___3 == 0),])
haem_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_renal___3 == 1),])
haem_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_renal___3 == 0),])

df_haem <- data.frame(
  gpa = c(haem_gpa_pres_acr,haem_gpa_abs_acr),
  mpa = c(haem_mpa_pres_acr,haem_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_haem, conf.int = TRUE)

#GFR 50-80ml/min/1.73 m2 *s
GFR1_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_renal___4 == 1),])
GFR1_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_renal___4 == 0),])
GFR1_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_renal___4 == 1),])
GFR1_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_renal___4 == 0),])

df_GFR1 <- data.frame(
  gpa = c(GFR1_gpa_pres_acr,GFR1_gpa_abs_acr),
  mpa = c(GFR1_mpa_pres_acr,GFR1_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_GFR1, conf.int = TRUE)

#GFR 15-49 ml/min/1.73 m2 *
GFR2_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_renal___5 == 1),])
GFR2_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_renal___5 == 0),])
GFR2_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_renal___5 == 1),])
GFR2_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_renal___5 == 0),])

df_GFR2 <- data.frame(
  gpa = c(GFR2_gpa_pres_acr,GFR2_gpa_abs_acr),
  mpa = c(GFR2_mpa_pres_acr,GFR2_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_GFR2, conf.int = TRUE)

#GFR < 15 ml/min/1.73m2 *
GFR3_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_renal___6 == 1),])
GFR3_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_renal___6 == 0),])
GFR3_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_renal___6 == 1),])
GFR3_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_renal___6 == 0),])

df_GFR3 <- data.frame(
  gpa = c(GFR3_gpa_pres_acr,GFR3_gpa_abs_acr),
  mpa = c(GFR3_mpa_pres_acr,GFR3_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_GFR3, conf.int = TRUE)

#Rise in creatinine > 10% or Creatinine clearance (GFR) fall > 25%
creat_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_renal___8 == 1),])
creat_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_renal___8 == 0),])
creat_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_renal___8 == 1),])
creat_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_renal___8 == 0),])

df_creat <- data.frame(
  gpa = c(creat_gpa_pres_acr,creat_gpa_abs_acr),
  mpa = c(creat_mpa_pres_acr,creat_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_creat, conf.int = TRUE)

#Nervous system involvement 
cns_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$cns == 1),])
cns_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$cns == 0),])
cns_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$cns == 1),])
cns_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$cns == 0),])

df_cns <- data.frame(
  gpa = c(cns_gpa_pres_acr,cns_gpa_abs_acr),
  mpa = c(cns_mpa_pres_acr,cns_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_cns, conf.int = TRUE)

#Headache
head_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___1 == 1),])
head_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___1 == 0),])
head_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___1 == 1),])
head_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___1 == 0),])

df_head <- data.frame(
  gpa = c(head_gpa_pres_acr,head_gpa_abs_acr),
  mpa = c(head_mpa_pres_acr,head_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_head, conf.int = TRUE)

#Meningitis/encephalitis
meni_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___2 == 1),])
meni_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___2 == 0),])
meni_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___2 == 1),])
meni_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___2 == 0),])

df_meni <- data.frame(
  gpa = c(meni_gpa_pres_acr,meni_gpa_abs_acr),
  mpa = c(meni_mpa_pres_acr,meni_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_meni, conf.int = TRUE)

#Organic confusion/cognitive dysfunction
cog_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___3 == 1),])
cog_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___3 == 0),])
cog_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___3 == 1),])
cog_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___3 == 0),])

df_cog <- data.frame(
  gpa = c(cog_gpa_pres_acr,cog_gpa_abs_acr),
  mpa = c(cog_mpa_pres_acr,cog_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_cog, conf.int = TRUE)

#Seizures (not hypertensive)
siez_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___4 == 1),])
siez_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___4 == 0),])
siez_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___4 == 1),])
siez_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___4 == 0),])

df_siez <- data.frame(
  gpa = c(siez_gpa_pres_acr,siez_gpa_abs_acr),
  mpa = c(siez_mpa_pres_acr,siez_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_siez, conf.int = TRUE)

#Stroke
stroke_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___5 == 1),])
stroke_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___5 == 0),])
stroke_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___5 == 1),])
stroke_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___5 == 0),])

df_stroke <- data.frame(
  gpa = c(stroke_gpa_pres_acr,stroke_gpa_abs_acr),
  mpa = c(stroke_mpa_pres_acr,stroke_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_stroke, conf.int = TRUE)

#Cord lesion
cl_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___6 == 1),])
cl_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___6 == 0),])
cl_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___6 == 1),])
cl_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___6 == 0),])

df_cl <- data.frame(
  gpa = c(cl_gpa_pres_acr,cl_gpa_abs_acr),
  mpa = c(cl_mpa_pres_acr,cl_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_cl, conf.int = TRUE)

#Cranial nerve palsy
cnp_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___7 == 1),])
cnp_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___7 == 0),])
cnp_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___7 == 1),])
cnp_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___7 == 0),])

df_cnp <- data.frame(
  gpa = c(cnp_gpa_pres_acr,cnp_gpa_abs_acr),
  mpa = c(cnp_mpa_pres_acr,cnp_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_cnp, conf.int = TRUE)

#Sensory peripheral neuropathy
spn_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___8 == 1),])
spn_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___8 == 0),])
spn_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___8 == 1),])
spn_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___8 == 0),])

df_spn <- data.frame(
  gpa = c(spn_gpa_pres_acr,spn_gpa_abs_acr),
  mpa = c(spn_mpa_pres_acr,spn_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_spn, conf.int = TRUE)

#Motor mononeuritis multiplex
mmm_gpa_pres_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___9 == 1),])
mmm_gpa_abs_acr <- nrow(GPA_ACR[which(GPA_ACR$pvas_cns___9 == 0),])
mmm_mpa_pres_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___9 == 1),])
mmm_mpa_abs_acr <- nrow(MPA_ACR[which(MPA_ACR$pvas_cns___9 == 0),])

df_mmm <- data.frame(
  gpa = c(mmm_gpa_pres_acr,mmm_gpa_abs_acr),
  mpa = c(mmm_mpa_pres_acr,mmm_mpa_abs_acr),
  row.names = c('Present','Absent'))

fisher.test(df_mmm, conf.int = TRUE)

#### SUPP Table 2a ####
#GENERAL
gen_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$gen == 1),])
gen_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$gen == 0),])
gen_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$gen == 1),])
gen_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$gen == 0),])

df2_gen <- data.frame(
  gpa = c(gen_gpa_pres_ema,gen_gpa_abs_ema),
  mpa = c(gen_mpa_pres_ema,gen_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_gen, conf.int = TRUE)

#Myalgia
mya_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gen1___1 == 1),])
mya_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gen1___1 == 0),])
mya_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gen1___1 == 1),])
mya_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gen1___1 == 0),])

df2_mya <- data.frame(
  gpa = c(mya_gpa_pres_ema,mya_gpa_abs_ema),
  mpa = c(mya_mpa_pres_ema,mya_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_mya, conf.int = TRUE)

#Arthralgia or arthritis
arth_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gen1___2 == 1),])
arth_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gen1___2 == 0),])
arth_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gen1___2 == 1),])
arth_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gen1___2 == 0),])

df2_arth <- data.frame(
  gpa = c(arth_gpa_pres_ema,arth_gpa_abs_ema),
  mpa = c(arth_mpa_pres_ema,arth_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_arth, conf.int = TRUE)

#Fever >= 38 deg C
fev_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gen1___3 == 1),])
fev_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gen1___3 == 0),])
fev_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gen1___3 == 1),])
fev_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gen1___3 == 0),])

df2_fev <- data.frame(
  gpa = c(fev_gpa_pres_ema,fev_gpa_abs_ema),
  mpa = c(fev_mpa_pres_ema,fev_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_fev, conf.int = TRUE)

#Weight loss >= 5% body weight
wl_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gen1___4 == 1),])
wl_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gen1___4 == 0),])
wl_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gen1___4 == 1),])
wl_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gen1___4 == 0),])

df2_wl <- data.frame(
  gpa = c(wl_gpa_pres_ema,wl_gpa_abs_ema),
  mpa = c(wl_mpa_pres_ema,wl_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_wl, conf.int = TRUE)

#Cutaneous involvement
cuta_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$cuta == 1),])
cuta_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$cuta == 0),])
cuta_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$cuta == 1),])
cuta_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$cuta == 0),])

df2_cuta <- data.frame(
  gpa = c(cuta_gpa_pres_ema,cuta_gpa_abs_ema),
  mpa = c(cuta_mpa_pres_ema,cuta_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_cuta, conf.int = TRUE)

#Polymorphous exanthema
pe_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___1 == 1),])
pe_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___1 == 0),])
pe_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___1 == 1),])
pe_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___1 == 0),])

df2_pe <- data.frame(
  gpa = c(pe_gpa_pres_ema,pe_gpa_abs_ema),
  mpa = c(pe_mpa_pres_ema,pe_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_pe, conf.int = TRUE)

#Livedo
liv_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___2 == 1),])
liv_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___2 == 0),])
liv_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___2 == 1),])
liv_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___2 == 0),])

df2_liv <- data.frame(
  gpa = c(liv_gpa_pres_ema,liv_gpa_abs_ema),
  mpa = c(liv_mpa_pres_ema,liv_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_liv, conf.int = TRUE)

#Panniculitis
pan_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___3 == 1),])
pan_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___3 == 0),])
pan_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___3 == 1),])
pan_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___3 == 0),])

df2_pan <- data.frame(
  gpa = c(pan_gpa_pres_ema,pan_gpa_abs_ema),
  mpa = c(pan_mpa_pres_ema,pan_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_pan, conf.int = TRUE)

#Purpura
purp_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___4 == 1),])
purp_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___4 == 0),])
purp_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___4 == 1),])
purp_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___4 == 0),])

df2_purp <- data.frame(
  gpa = c(purp_gpa_pres_ema,purp_gpa_pres_ema),
  mpa = c(purp_mpa_pres_ema,purp_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_purp, conf.int = TRUE)

#Skin nodules
sn_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___5 == 1),])
sn_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___5 == 0),])
sn_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___5 == 1),])
sn_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___5 == 0),])

df2_sn <- data.frame(
  gpa = c(sn_gpa_pres_ema,sn_gpa_abs_ema),
  mpa = c(sn_mpa_pres_ema,sn_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_sn, conf.int = TRUE)

#Infarct (nail edge lesion, splinter haemorrhage)
inf_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___6 == 1),])
inf_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___6 == 0),])
inf_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___6 == 1),])
inf_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___6 == 0),])

df2_inf <- data.frame(
  gpa = c(inf_gpa_pres_ema,inf_gpa_abs_ema),
  mpa = c(inf_mpa_pres_ema,inf_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_inf, conf.int = TRUE)

#Ulcer (full-thickness necrosis)
ulc_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___7 == 1),])
ulc_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___7 == 0),])
ulc_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___7 == 1),])
ulc_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___7 == 0),])

df2_ulc <- data.frame(
  gpa = c(ulc_gpa_pres_ema,ulc_gpa_abs_ema),
  mpa = c(ulc_mpa_pres_ema,ulc_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_ulc, conf.int = TRUE)

#Gangrene (extensive necrosis)
gang_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___8 == 1),])
gang_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___8 == 0),])
gang_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___8 == 1),])
gang_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___8 == 0),])

df2_gang <- data.frame(
  gpa = c(gang_gpa_pres_ema,gang_gpa_abs_ema),
  mpa = c(gang_mpa_pres_ema,gang_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_gang, conf.int = TRUE)

#Other skin vasculitis (specify below)
osv_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___9 == 1),])
osv_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cut2___9 == 0),])
osv_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___9 == 1),])
osv_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cut2___9 == 0),])

df2_osv <- data.frame(
  gpa = c(osv_gpa_pres_ema,osv_gpa_abs_ema),
  mpa = c(osv_mpa_pres_ema,osv_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_osv, conf.int = TRUE)

#Mucous membranes/eyes involvement
muc_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$muc == 1),])
muc_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$muc == 0),])
muc_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$muc == 1),])
muc_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$muc == 0),])

df2_muc <- data.frame(
  gpa = c(muc_gpa_pres_ema,muc_gpa_abs_ema),
  mpa = c(muc_mpa_pres_ema,muc_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_muc, conf.int = TRUE)

#Mouth ulcers/granulomata
mug_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___1 == 1),])
mug_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___1 == 0),])
mug_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___1 == 1),])
mug_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___1 == 0),])

df2_mug <- data.frame(
  gpa = c(mug_gpa_pres_ema,mug_gpa_abs_ema),
  mpa = c(mug_mpa_pres_ema,mug_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_mug, conf.int = TRUE)

#Genital ulcers
gulc_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___2 == 1),])
gulc_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___2 == 0),])
gulc_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___2 == 1),])
gulc_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___2 == 0),])

df2_gulc <- data.frame(
  gpa = c(gulc_gpa_pres_ema,gulc_gpa_abs_ema),
  mpa = c(gulc_mpa_pres_ema,gulc_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_gulc, conf.int = TRUE)

#Adnexal inflammation (sialadenitis/ dEMAyocystitis)
adn_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___3 == 1),])
adn_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___3 == 0),])
adn_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___3 == 1),])
adn_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___3 == 0),])

df2_adn <- data.frame(
  gpa = c(adn_gpa_pres_ema,adn_gpa_abs_ema),
  mpa = c(adn_mpa_pres_ema,adn_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_adn, conf.int = TRUE)

#Significant proptosis
prop_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___4 == 1),])
prop_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___4 == 0),])
prop_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___4 == 1),])
prop_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___4 == 0),])

df2_prop <- data.frame(
  gpa = c(prop_gpa_pres_ema,prop_gpa_abs_ema),
  mpa = c(prop_mpa_pres_ema,prop_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_prop, conf.int = TRUE)

#Red eye (Epi)scleritis
red_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___5 == 1),])
red_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___5 == 0),])
red_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___5 == 1),])
red_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___5 == 0),])

df2_red <- data.frame(
  gpa = c(red_gpa_pres_ema,red_gpa_abs_ema),
  mpa = c(red_mpa_pres_ema,red_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_red, conf.int = TRUE)

#Red eye conjunctivitis/ blepharitis/keratitis
redcon_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___6 == 1),])
redcon_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___6 == 0),])
redcon_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___6 == 1),])
redcon_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___6 == 0),])

df2_redcon <- data.frame(
  gpa = c(redcon_gpa_pres_ema,redcon_gpa_abs_ema),
  mpa = c(redcon_mpa_pres_ema,redcon_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_redcon, conf.int = TRUE)

#Uveitis
uve_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___7 == 1),])
uve_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___7 == 0),])
uve_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___7 == 1),])
uve_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___7 == 0),])

df2_uve <- data.frame(
  gpa = c(uve_gpa_pres_ema,uve_gpa_abs_ema),
  mpa = c(uve_mpa_pres_ema,uve_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_uve, conf.int = TRUE)

#Blurred vision
bv_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___8 == 1),])
bv_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___8 == 0),])
bv_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___8 == 1),])
bv_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___8 == 0),])

df2_bv <- data.frame(
  gpa = c(bv_gpa_pres_ema,bv_gpa_abs_ema),
  mpa = c(bv_mpa_pres_ema,bv_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_bv, conf.int = TRUE)

#Sudden visual loss
svl_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___9 == 1),])
svl_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___9 == 0),])
svl_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___9 == 1),])
svl_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___9 == 0),])

df2_svl <- data.frame(
  gpa = c(svl_gpa_pres_ema,svl_gpa_abs_ema),
  mpa = c(svl_mpa_pres_ema,svl_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_svl, conf.int = TRUE)

#Retinal vasculitis/retinal vessel thrombosis/retinal exudates/haemorrhages
rvas_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___10 == 1),])
rvas_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_mucous2___10 == 0),])
rvas_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___10 == 1),])
rvas_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_mucous2___10 == 0),])

df2_rvas <- data.frame(
  gpa = c(rvas_gpa_pres_ema,rvas_gpa_abs_ema),
  mpa = c(rvas_mpa_pres_ema,rvas_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_rvas, conf.int = TRUE)

#ENT
ent_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$ent == 1),])
ent_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$ent == 0),])
ent_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$ent == 1),])
ent_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$ent == 0),])

df2_ent <- data.frame(
  gpa = c(ent_gpa_pres_ema,ent_gpa_abs_ema),
  mpa = c(ent_mpa_pres_ema,ent_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_ent, conf.int = TRUE)

#Bloody nasal discharge/crusts/ulcers/granuloma
bnd_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_ent___1 == 1),])
bnd_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_ent___1 == 0),])
bnd_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_ent___1 == 1),])
bnd_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_ent___1 == 0),])

df2_bnd <- data.frame(
  gpa = c(bnd_gpa_pres_ema,bnd_gpa_abs_ema),
  mpa = c(bnd_mpa_pres_ema,bnd_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_bnd, conf.int = TRUE)

#Paranasal sinus involvement
psi_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_ent___2 == 1),])
psi_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_ent___2 == 0),])
psi_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_ent___2 == 1),])
psi_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_ent___2 == 0),])

df2_psi <- data.frame(
  gpa = c(psi_gpa_pres_ema,psi_gpa_abs_ema),
  mpa = c(psi_mpa_pres_ema,psi_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_psi, conf.int = TRUE)

#Subglottic stenosis/ hoarseness /stridor
sgs_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_ent___3 == 1),])
sgs_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_ent___3 == 0),])
sgs_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_ent___3 == 1),])
sgs_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_ent___3 == 0),])

df2_sgs <- data.frame(
  gpa = c(sgs_gpa_pres_ema,sgs_gpa_abs_ema),
  mpa = c(sgs_mpa_pres_ema,sgs_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_sgs, conf.int = TRUE)

#Conductive hearing loss
chl_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_ent___4 == 1),])
chl_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_ent___4 == 0),])
chl_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_ent___4 == 1),])
chl_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_ent___4 == 0),])

df2_chl <- data.frame(
  gpa = c(chl_gpa_pres_ema,chl_gpa_abs_ema),
  mpa = c(chl_mpa_pres_ema,chl_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_chl, conf.int = TRUE)

#Sensorineural hearing loss
shl_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_ent___5 == 1),])
shl_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_ent___5 == 0),])
shl_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_ent___5 == 1),])
shl_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_ent___5 == 0),])

df2_shl <- data.frame(
  gpa = c(shl_gpa_pres_ema,shl_gpa_abs_ema),
  mpa = c(shl_mpa_pres_ema,shl_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_shl, conf.int = TRUE)

#Chest Involvement
chest_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$chest == 1),])
chest_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$chest == 0),])
chest_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$chest == 1),])
chest_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$chest == 0),])

df2_chest <- data.frame(
  gpa = c(chest_gpa_pres_ema,chest_gpa_abs_ema),
  mpa = c(chest_mpa_pres_ema,chest_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_chest, conf.int = TRUE)

#Wheeze or expiratory dyspnea
wheeze_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_chest___1 == 1),])
wheeze_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_chest___1 == 0),])
wheeze_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_chest___1 == 1),])
wheeze_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_chest___1 == 0),])

df2_wheeze <- data.frame(
  gpa = c(wheeze_gpa_pres_ema,wheeze_gpa_abs_ema),
  mpa = c(wheeze_mpa_pres_ema,wheeze_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_wheeze, conf.int = TRUE)

#Endobronchial/endotracheal involvement
ebi_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_chest___2 == 1),])
ebi_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_chest___2 == 0),])
ebi_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_chest___2 == 1),])
ebi_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_chest___2 == 0),])

df2_ebi <- data.frame(
  gpa = c(ebi_gpa_pres_ema,ebi_gpa_abs_ema),
  mpa = c(ebi_mpa_pres_ema,ebi_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_ebi, conf.int = TRUE)

#Nodules or cavities 
nod_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_chest___3 == 1),])
nod_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_chest___3 == 0),])
nod_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_chest___3 == 1),])
nod_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_chest___3 == 0),])

df2_nod <- data.frame(
  gpa = c(nod_gpa_pres_ema,nod_gpa_abs_ema),
  mpa = c(nod_mpa_pres_ema,nod_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_nod, conf.int = TRUE)

#Pleural effusion/pleurisy
peff_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_chest___4 == 1),])
peff_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_chest___4 == 0),])
peff_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_chest___4 == 1),])
peff_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_chest___4 == 0),])

df2_peff <- data.frame(
  gpa = c(peff_gpa_pres_ema,peff_gpa_abs_ema),
  mpa = c(peff_mpa_pres_ema,peff_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_peff, conf.int = TRUE)

#Infiltrate
inf_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_chest___5 == 1),])
inf_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_chest___5 == 0),])
inf_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_chest___5 == 1),])
inf_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_chest___5 == 0),])

df2_inf <- data.frame(
  gpa = c(inf_gpa_pres_ema,inf_gpa_abs_ema),
  mpa = c(inf_mpa_pres_ema,inf_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_inf, conf.int = TRUE)

#Massive haemoptysis/Alveolar haemorrhage
mhaem_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_chest___6 == 1),])
mhaem_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_chest___6 == 0),])
mhaem_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_chest___6 == 1),])
mhaem_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_chest___6 == 0),])

df2_mhaem <- data.frame(
  gpa = c(mhaem_gpa_pres_ema,mhaem_gpa_abs_ema),
  mpa = c(mhaem_mpa_pres_ema,mhaem_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_mhaem, conf.int = TRUE)

#Respiratory failure
respf_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_chest___7 == 1),])
respf_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_chest___7 == 0),])
respf_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_chest___7 == 1),])
respf_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_chest___7 == 0),])

df2_respf <- data.frame(
  gpa = c(respf_gpa_pres_ema,respf_gpa_abs_ema),
  mpa = c(respf_mpa_pres_ema,respf_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_respf, conf.int = TRUE)

#Cardiac Involvement
cardio_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$cardio == 1),])
cardio_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$cardio == 0),])
cardio_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$cardio == 1),])
cardio_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$cardio == 0),])

df2_cardio <- data.frame(
  gpa = c(cardio_gpa_pres_ema,cardio_gpa_abs_ema),
  mpa = c(cardio_mpa_pres_ema,cardio_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_cardio, conf.int = TRUE)

#Loss of pulses
lop_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___1 == 1),])
lop_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___1 == 0),])
lop_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___1 == 1),])
lop_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___1 == 0),])

df2_lop <- data.frame(
  gpa = c(lop_gpa_pres_ema,lop_gpa_abs_ema),
  mpa = c(lop_mpa_pres_ema,lop_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_lop, conf.int = TRUE)

#Bruits over accessible arteries
boaa_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___2 == 1),])
boaa_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___2 == 0),])
boaa_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___2 == 1),])
boaa_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___2 == 0),])

df2_boaa <- data.frame(
  gpa = c(boaa_gpa_pres_ema,boaa_gpa_abs_ema),
  mpa = c(boaa_mpa_pres_ema,boaa_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_boaa, conf.int = TRUE)

#Blood pressure discrepancy
bpd_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___3 == 1),])
bpd_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___3 == 0),])
bpd_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___3 == 1),])
bpd_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___3 == 0),])

df2_bpd <- data.frame(
  gpa = c(bpd_gpa_pres_ema,bpd_gpa_abs_ema),
  mpa = c(bpd_mpa_pres_ema,bpd_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_bpd, conf.int = TRUE)

#Claudication of extremities
coe_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___4 == 1),])
coe_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___4 == 0),])
coe_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___4 == 1),])
coe_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___4 == 0),])

df2_coe <- data.frame(
  gpa = c(coe_gpa_pres_ema,coe_gpa_abs_ema),
  mpa = c(coe_mpa_pres_ema,coe_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_coe, conf.int = TRUE)

#Ischaemic cardiac pain
icp_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___5 == 1),])
icp_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___5 == 0),])
icp_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___5 == 1),])
icp_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___5 == 0),])

df2_icp <- data.frame(
  gpa = c(icp_gpa_pres_ema,icp_gpa_abs_ema),
  mpa = c(icp_mpa_pres_ema,icp_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_icp, conf.int = TRUE)

#Cardiomyopathy
cmpy_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___6 == 1),])
cmpy_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___6 == 0),])
cmpy_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___6 == 1),])
cmpy_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___6 == 0),])

df2_cmpy <- data.frame(
  gpa = c(cmpy_gpa_pres_ema,cmpy_gpa_abs_ema),
  mpa = c(cmpy_mpa_pres_ema,cmpy_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_cmpy, conf.int = TRUE)

#Congestive cardiac failure
ccf_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___7 == 1),])
ccf_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___7 == 0),])
ccf_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___7 == 1),])
ccf_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___7 == 0),])

df2_ccf <- data.frame(
  gpa = c(ccf_gpa_pres_ema,ccf_gpa_abs_ema),
  mpa = c(ccf_mpa_pres_ema,ccf_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_ccf, conf.int = TRUE)

#Valvular heart disease
vhd_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___8 == 1),])
vhd_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___8 == 0),])
vhd_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___8 == 1),])
vhd_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___8 == 0),])

df2_vhd <- data.frame(
  gpa = c(vhd_gpa_pres_ema,vhd_gpa_abs_ema),
  mpa = c(vhd_mpa_pres_ema,vhd_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_vhd, conf.int = TRUE)

#Pericarditis
pcd_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___9 == 1),])
pcd_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cv___9 == 0),])
pcd_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___9 == 1),])
pcd_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cv___9 == 0),])

df2_pcd <- data.frame(
  gpa = c(pcd_gpa_pres_ema,pcd_gpa_abs_ema),
  mpa = c(pcd_mpa_pres_ema,pcd_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_pcd, conf.int = TRUE)

#Gastrointestinal Involvement
gi_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$gi == 1),])
gi_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$gi == 0),])
gi_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$gi == 1),])
gi_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$gi == 0),])

df2_gi <- data.frame(
  gpa = c(gi_gpa_pres_ema,gi_gpa_abs_ema),
  mpa = c(gi_mpa_pres_ema,gi_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_gi, conf.int = TRUE)

#Abdominal pain
ap_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gi___1 == 1),])
ap_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gi___1 == 0),])
ap_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gi___1 == 1),])
ap_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$gpvas_gi___1 == 0),])

df2_ap <- data.frame(
  gpa = c(ap_gpa_pres_ema,ap_gpa_abs_ema),
  mpa = c(ap_mpa_pres_ema,ap_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_ap, conf.int = TRUE)

#Peritonitis
pos_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gi___2 == 1),])
pos_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gi___2 == 0),])
pos_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gi___2 == 1),])
pos_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gi___2 == 0),])

df2_pos <- data.frame(
  gpa = c(pos_gpa_pres_ema,pos_gpa_abs_ema),
  mpa = c(pos_mpa_pres_ema,pos_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_pos, conf.int = TRUE)

#Blood in stools or bloody diarrhea
bis_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gi___3 == 1),])
bis_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gi___3 == 0),])
bis_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gi___3 == 1),])
bis_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gi___3 == 0),])

df2_bis <- data.frame(
  gpa = c(bis_gpa_pres_ema,bis_gpa_abs_ema),
  mpa = c(bis_mpa_pres_ema,bis_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_bis, conf.int = TRUE)

#Bowel ischaemia
bi_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gi___4 == 1),])
bi_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_gi___4 == 0),])
bi_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gi___4 == 1),])
bi_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_gi___4 == 0),])

df2_bi <- data.frame(
  gpa = c(bi_gpa_pres_ema,bi_gpa_abs_ema),
  mpa = c(bi_mpa_pres_ema,bi_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_bi, conf.int = TRUE)

#Renal Involvement
ren_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$ren == 1),])
ren_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$ren == 0),])
ren_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$ren == 1),])
ren_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$ren == 0),])

df2_ren <- data.frame(
  gpa = c(ren_gpa_pres_ema,ren_gpa_abs_ema),
  mpa = c(ren_mpa_pres_ema,ren_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_ren, conf.int = TRUE)

#Hypertension >95th centile (for height)
hyp_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_renal___1 == 1),])
hyp_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_renal___1 == 0),])
hyp_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_renal___1 == 1),])
hyp_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_renal___1 == 0),])

df2_hyp <- data.frame(
  gpa = c(hyp_gpa_pres_ema,hyp_gpa_abs_ema),
  mpa = c(hyp_mpa_pres_ema,hyp_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_hyp, conf.int = TRUE)

#Proteinuria >0.3 g/24h,>20mg/mmol creatinine
pru_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_renal___2 == 1),])
pru_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_renal___2 == 0),])
pru_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_renal___2 == 1),])
pru_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_renal___2 == 0),])

df2_pru <- data.frame(
  gpa = c(pru_gpa_pres_ema,pru_gpa_abs_ema),
  mpa = c(pru_mpa_pres_ema,pru_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_pru, conf.int = TRUE)

#Haematuria >= 2+ or 5 rbc/hpf or red cell casts
haem_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_renal___3 == 1),])
haem_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_renal___3 == 0),])
haem_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_renal___3 == 1),])
haem_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_renal___3 == 0),])

df2_haem <- data.frame(
  gpa = c(haem_gpa_pres_ema,haem_gpa_abs_ema),
  mpa = c(haem_mpa_pres_ema,haem_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_haem, conf.int = TRUE)

#GFR 50-80ml/min/1.73 m2 *s
GFR1_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_renal___4 == 1),])
GFR1_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_renal___4 == 0),])
GFR1_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_renal___4 == 1),])
GFR1_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_renal___4 == 0),])

df2_GFR1 <- data.frame(
  gpa = c(GFR1_gpa_pres_ema,GFR1_gpa_abs_ema),
  mpa = c(GFR1_mpa_pres_ema,GFR1_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_GFR1, conf.int = TRUE)

#GFR 15-49 ml/min/1.73 m2 *
GFR2_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_renal___5 == 1),])
GFR2_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_renal___5 == 0),])
GFR2_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_renal___5 == 1),])
GFR2_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_renal___5 == 0),])

df2_GFR2 <- data.frame(
  gpa = c(GFR2_gpa_pres_ema,GFR2_gpa_abs_ema),
  mpa = c(GFR2_mpa_pres_ema,GFR2_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_GFR2, conf.int = TRUE)

#GFR < 15 ml/min/1.73m2 *
GFR3_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_renal___6 == 1),])
GFR3_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_renal___6 == 0),])
GFR3_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_renal___6 == 1),])
GFR3_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_renal___6 == 0),])

df2_GFR3 <- data.frame(
  gpa = c(GFR3_gpa_pres_ema,GFR3_gpa_abs_ema),
  mpa = c(GFR3_mpa_pres_ema,GFR3_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_GFR3, conf.int = TRUE)

#Rise in creatinine > 10% or Creatinine clearance (GFR) fall > 25%
creat_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_renal___8 == 1),])
creat_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_renal___8 == 0),])
creat_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_renal___8 == 1),])
creat_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_renal___8 == 0),])

df2_creat <- data.frame(
  gpa = c(creat_gpa_pres_ema,creat_gpa_abs_ema),
  mpa = c(creat_mpa_pres_ema,creat_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_creat, conf.int = TRUE)

#Nervous system involvement 
cns_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$cns == 1),])
cns_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$cns == 0),])
cns_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$cns == 1),])
cns_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$cns == 0),])

df2_cns <- data.frame(
  gpa = c(cns_gpa_pres_ema,cns_gpa_abs_ema),
  mpa = c(cns_mpa_pres_ema,cns_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_cns, conf.int = TRUE)

#Headache
head_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___1 == 1),])
head_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___1 == 0),])
head_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___1 == 1),])
head_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___1 == 0),])

df2_head <- data.frame(
  gpa = c(head_gpa_pres_ema,head_gpa_abs_ema),
  mpa = c(head_mpa_pres_ema,head_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_head, conf.int = TRUE)

#Meningitis/encephalitis
meni_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___2 == 1),])
meni_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___2 == 0),])
meni_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___2 == 1),])
meni_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___2 == 0),])

df2_meni <- data.frame(
  gpa = c(meni_gpa_pres_ema,meni_gpa_abs_ema),
  mpa = c(meni_mpa_pres_ema,meni_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_meni, conf.int = TRUE)

#Organic confusion/cognitive dysfunction
cog_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___3 == 1),])
cog_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___3 == 0),])
cog_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___3 == 1),])
cog_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___3 == 0),])

df2_cog <- data.frame(
  gpa = c(cog_gpa_pres_ema,cog_gpa_abs_ema),
  mpa = c(cog_mpa_pres_ema,cog_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_cog, conf.int = TRUE)

#Seizures (not hypertensive)
siez_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___4 == 1),])
siez_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___4 == 0),])
siez_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___4 == 1),])
siez_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___4 == 0),])

df2_siez <- data.frame(
  gpa = c(siez_gpa_pres_ema,siez_gpa_abs_ema),
  mpa = c(siez_mpa_pres_ema,siez_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_siez, conf.int = TRUE)

#Stroke
stroke_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___5 == 1),])
stroke_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___5 == 0),])
stroke_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___5 == 1),])
stroke_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___5 == 0),])

df2_stroke <- data.frame(
  gpa = c(stroke_gpa_pres_ema,stroke_gpa_abs_ema),
  mpa = c(stroke_mpa_pres_ema,stroke_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_stroke, conf.int = TRUE)

#Cord lesion
cl_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___6 == 1),])
cl_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___6 == 0),])
cl_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___6 == 1),])
cl_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___6 == 0),])

df2_cl <- data.frame(
  gpa = c(cl_gpa_pres_ema,cl_gpa_abs_ema),
  mpa = c(cl_mpa_pres_ema,cl_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_cl, conf.int = TRUE)

#Cranial nerve palsy
cnp_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___7 == 1),])
cnp_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___7 == 0),])
cnp_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___7 == 1),])
cnp_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___7 == 0),])

df2_cnp <- data.frame(
  gpa = c(cnp_gpa_pres_ema,cnp_gpa_abs_ema),
  mpa = c(cnp_mpa_pres_ema,cnp_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_cnp, conf.int = TRUE)

#Sensory peripheral neuropathy
spn_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___8 == 1),])
spn_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___8 == 0),])
spn_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___8 == 1),])
spn_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___8 == 0),])

df2_spn <- data.frame(
  gpa = c(spn_gpa_pres_ema,spn_gpa_abs_ema),
  mpa = c(spn_mpa_pres_ema,spn_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_spn, conf.int = TRUE)

#Motor mononeuritis multiplex
mmm_gpa_pres_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___9 == 1),])
mmm_gpa_abs_ema <- nrow(GPA_EMA[which(GPA_EMA$pvas_cns___9 == 0),])
mmm_mpa_pres_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___9 == 1),])
mmm_mpa_abs_ema <- nrow(MPA_EMA[which(MPA_EMA$pvas_cns___9 == 0),])

df2_mmm <- data.frame(
  gpa = c(mmm_gpa_pres_ema,mmm_gpa_abs_ema),
  mpa = c(mmm_mpa_pres_ema,mmm_mpa_abs_ema),
  row.names = c('Present','Absent'))

fisher.test(df2_mmm, conf.int = TRUE)

#### FIG 4a ####
gpa_gen_ema <- nrow(EMA_GPA[which(EMA_GPA$gen == 1),])
  
gpa_gen_acr <- nrow(ACR_GPA_all[which(ACR_GPA_all$gen == 1),])

gpa_cuta_ema <- nrow(EMA_GPA[which(EMA_GPA$cuta == 1),])

gpa_cuta_acr <- nrow(ACR_GPA_all[which(ACR_GPA_all$cuta == 1),])

gpa_muc_ema <- nrow(EMA_GPA[which(EMA_GPA$muc == 1),])

gpa_muc_acr <- nrow(ACR_GPA_all[which(ACR_GPA_all$muc == 1),])

gpa_cardio_ema <- nrow(EMA_GPA[which(EMA_GPA$cardio == 1),])

gpa_cardio_acr <- nrow(ACR_GPA_all[which(ACR_GPA_all$cardio == 1),])

gpa_ent_ema <- nrow(EMA_GPA[which(EMA_GPA$ent == 1),])

gpa_ent_acr <- nrow(ACR_GPA_all[which(ACR_GPA_all$ent == 1),])

gpa_chest_ema <- nrow(EMA_GPA[which(EMA_GPA$chest == 1),])

gpa_chest_acr <- nrow(ACR_GPA_all[which(ACR_GPA_all$chest == 1),])

gpa_gi_ema <- nrow(EMA_GPA[which(EMA_GPA$gi == 1),])

gpa_gi_acr <- nrow(ACR_GPA_all[which(ACR_GPA_all$gi == 1),])

gpa_ren_ema <- nrow(EMA_GPA[which(EMA_GPA$ren == 1),])

gpa_ren_acr <- nrow(ACR_GPA_all[which(ACR_GPA_all$ren == 1),])

gpa_cns_ema <- nrow(EMA_GPA[which(EMA_GPA$cns == 1),])

gpa_cns_acr <- nrow(ACR_GPA_all[which(ACR_GPA_all$cns == 1),])

gpa_pr3_ema <- nrow(EMA_GPA[which(EMA_GPA$pr3_result == 2),])

gpa_pr3_acr <- nrow(ACR_GPA_all[which(ACR_GPA_all$pr3_result == 2),])

gpa_mpo_ema <- nrow(EMA_GPA[which(EMA_GPA$mpo_result == 2),])

gpa_mpo_acr <- nrow(ACR_GPA_all[which(ACR_GPA_all$mpo_result == 2),])



mpa_gen_ema <- nrow(EMA_MPA[which(EMA_MPA$gen == 1),])

mpa_gen_acr <- nrow(ACR_MPA_all[which(ACR_MPA_all$gen == 1),])

mpa_cuta_ema <- nrow(EMA_MPA[which(EMA_MPA$cuta == 1),])

mpa_cuta_acr <- nrow(ACR_MPA_all[which(ACR_MPA_all$cuta == 1),])

mpa_muc_ema <- nrow(EMA_MPA[which(EMA_MPA$muc == 1),])

mpa_muc_acr <- nrow(ACR_MPA_all[which(ACR_MPA_all$muc == 1),])

mpa_cardio_ema <- nrow(EMA_MPA[which(EMA_MPA$cardio == 1),])

mpa_cardio_acr <- nrow(ACR_MPA_all[which(ACR_MPA_all$cardio == 1),])

mpa_ent_ema <- nrow(EMA_MPA[which(EMA_MPA$ent == 1),])

mpa_ent_acr <- nrow(ACR_MPA_all[which(ACR_MPA_all$ent == 1),])

mpa_chest_ema <- nrow(EMA_MPA[which(EMA_MPA$chest == 1),])

mpa_chest_acr <- nrow(ACR_MPA_all[which(ACR_MPA_all$chest == 1),])

mpa_gi_ema <- nrow(EMA_MPA[which(EMA_MPA$gi == 1),])

mpa_gi_acr <- nrow(ACR_MPA_all[which(ACR_MPA_all$gi == 1),])

mpa_ren_ema <- nrow(EMA_MPA[which(EMA_MPA$ren == 1),])

mpa_ren_acr <- nrow(ACR_MPA_all[which(ACR_MPA_all$ren == 1),])

mpa_cns_ema <- nrow(EMA_MPA[which(EMA_MPA$cns == 1),])

mpa_cns_acr <- nrow(ACR_MPA_all[which(ACR_MPA_all$cns == 1),])

mpa_pr3_ema <- nrow(EMA_MPA[which(EMA_MPA$pr3_result == 2),])

mpa_pr3_acr <- nrow(ACR_MPA_all[which(ACR_MPA_all$pr3_result == 2),])

mpa_mpo_ema <- nrow(EMA_MPA[which(EMA_MPA$mpo_result == 2),])

mpa_mpo_acr <- nrow(ACR_MPA_all[which(ACR_MPA_all$mpo_result == 2),])


unc_gen_ema <- nrow(EMA_UNC[which(EMA_UNC$gen == 1),])

unc_gen_acr <- nrow(UNC_ACR[which(UNC_ACR$gen == 1),])

unc_cuta_ema <- nrow(EMA_UNC[which(EMA_UNC$cuta == 1),])

unc_cuta_acr <- nrow(UNC_ACR[which(UNC_ACR$cuta == 1),])

unc_muc_ema <- nrow(EMA_UNC[which(EMA_UNC$muc == 1),])

unc_muc_acr <- nrow(UNC_ACR[which(UNC_ACR$muc == 1),])

unc_cardio_ema <- nrow(EMA_UNC[which(EMA_UNC$cardio == 1),])

unc_cardio_acr <- nrow(UNC_ACR[which(UNC_ACR$cardio == 1),])

unc_ent_ema <- nrow(EMA_UNC[which(EMA_UNC$ent == 1),])

unc_ent_acr <- nrow(UNC_ACR[which(UNC_ACR$ent == 1),])

unc_chest_ema <- nrow(EMA_UNC[which(EMA_UNC$chest == 1),])

unc_chest_acr <- nrow(UNC_ACR[which(UNC_ACR$chest == 1),])

unc_gi_ema <- nrow(EMA_UNC[which(EMA_UNC$gi == 1),])

unc_gi_acr <- nrow(UNC_ACR[which(UNC_ACR$gi == 1),])

unc_ren_ema <- nrow(EMA_UNC[which(EMA_UNC$ren == 1),])

unc_ren_acr <- nrow(UNC_ACR[which(UNC_ACR$ren == 1),])

unc_cns_ema <- nrow(EMA_UNC[which(EMA_UNC$cns == 1),])

unc_cns_acr <- nrow(UNC_ACR[which(UNC_ACR$cns == 1),])

unc_pr3_ema <- nrow(EMA_UNC[which(EMA_UNC$pr3_result == 2),])

unc_pr3_acr <- nrow(UNC_ACR[which(UNC_ACR$pr3_result == 2),])

unc_mpo_ema <- nrow(EMA_UNC[which(EMA_UNC$mpo_result == 2),])

unc_mpo_acr <- nrow(UNC_ACR[which(UNC_ACR$mpo_result == 2),])

df2 <- data.frame(supp=rep(c("Ankara (n = 265)", "ACR/Eular (n = 287)"), each=11),
                  system=rep(c("General", "Cutaenous", "Mucous", "Cardiovascular",
                  "ENT", "Chest", "GI", "Renal", "CNS", "PR3", "MPO"),2),
                  len = c(gpa_gen_ema, gpa_cuta_ema,
                  gpa_muc_ema, gpa_muc_acr, gpa_cardio_ema, gpa_ent_ema, gpa_chest_ema,
                  gpa_gi_ema, gpa_ren_ema, gpa_cns_ema, gpa_pr3_ema, gpa_mpo_ema,
                  gpa_gen_acr, gpa_cuta_acr, gpa_cardio_acr, gpa_ent_acr, 
                  gpa_chest_acr, gpa_gi_acr,gpa_ren_acr, gpa_cns_acr, gpa_pr3_acr,
                  gpa_mpo_acr))

ggplot(data=df2, aes(x=factor(system, , level=c("General", "Cutaenous", "Mucous", "Cardiovascular",
                                                "ENT", "Chest", "GI", "Renal", "CNS", "PR3", "MPO")), y=len, fill=supp)) + geom_bar(stat="identity", color="black", width = 0.5, position=position_dodge()) + guides(x =  guide_axis(angle = 45)) + scale_fill_manual(values=c("white", "black")) + theme_minimal() + theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +  labs(title = "GPA",
       subtitle = waiver(),
       caption = waiver(),
       tag = waiver(),
       alt = waiver(),
       alt_insight = waiver(),
       x = "", y = "", fill = "")

#### FIG 4b ####

df3 <- data.frame(supp=rep(c("Ped-EMA (n = 72)", "ACR/Eular  (n = 159)"), each=11),
                  system=rep(c("General", "Cutaenous", "Mucous", "Cardiovascular",
                               "ENT", "Chest", "GI", "Renal", "CNS", "PR3", "MPO"),2),
                  len = c(mpa_gen_ema, mpa_cuta_ema,
                        mpa_muc_ema, mpa_muc_acr, mpa_cardio_ema, mpa_ent_ema, mpa_chest_ema,
                        mpa_gi_ema, mpa_ren_ema, mpa_cns_ema, mpa_pr3_ema, mpa_mpo_ema,
                        mpa_gen_acr, mpa_cuta_acr, mpa_cardio_acr, mpa_ent_acr, 
                        mpa_chest_acr, mpa_gi_acr,mpa_ren_acr, mpa_cns_acr, mpa_pr3_acr,
                        mpa_mpo_acr))


ggplot(data=df3, aes(x=factor(system, , level=c("General", "Cutaenous", "Mucous", "Cardiovascular",
                               "ENT", "Chest", "GI", "Renal", "CNS", "PR3", "MPO")), y=len, fill=supp)) + geom_bar(stat="identity", color="black", width = 0.5, position=position_dodge()) + guides(x =  guide_axis(angle = 45)) + scale_fill_manual(values=c("white", "black")) +
  theme_minimal() + theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) + labs(title = "MPA",
                         subtitle = waiver(),
                         caption = waiver(),
                         tag = waiver(),
                         alt = waiver(),
                         alt_insight = waiver(),
                         x = "", y = "", fill = "")  

#### FIG 4c ####

df4 <- data.frame(supp=rep(c("Ped-EMA (n = 202)", "ACR/Eular  (n = 147)"), each=11),
                  system=rep(c("General", "Cutaenous", "Mucous", "Cardiovascular",
                               "ENT", "Chest", "GI", "Renal", "CNS", "PR3", "MPO"),2),
                  len = c(unc_gen_ema, unc_cuta_ema,
                          unc_muc_ema, unc_muc_acr, unc_cardio_ema, unc_ent_ema, unc_chest_ema,
                          unc_gi_ema, unc_ren_ema, unc_cns_ema, unc_pr3_ema, unc_mpo_ema,
                          unc_gen_acr, unc_cuta_acr, unc_cardio_acr, unc_ent_acr, 
                          unc_chest_acr, unc_gi_acr,unc_ren_acr, unc_cns_acr, unc_pr3_acr,
                          unc_mpo_acr))

ggplot(data=df4, aes(x=factor(system, , level=c("General", "Cutaenous", "Mucous", "Cardiovascular",
                                                "ENT", "Chest", "GI", "Renal", "CNS", "PR3", "MPO")), y=len, fill=supp)) + geom_bar(stat="identity", color="black", width = 0.5, position=position_dodge()) + guides(x =  guide_axis(angle = 45)) + scale_fill_manual(values=c("white", "black")) +
  theme_minimal() + theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +  labs(title = "Unclassified",
                                                                                                   subtitle = waiver(),
                                                                                                   caption = waiver(),
                                                                                                   tag = waiver(),
                                                                                                   alt = waiver(),
                                                                                                   alt_insight = waiver(),
                                                                                                   x = "", y = "", fill = "")

#### KAPPA STAT CALC ####

## Reassigning diagnosis

data$EMAclass[data$EMAclass == "UCV" | data$EMAclass == "UCV AAV"] = "unclassified"
data$EMAclass[data$EMAclass == "GPA/MPA"] = "unclassified"

data$diagnosis[data$diagnosis == 6 | data$diagnosis >= 9] = "unclassified"
data$diagnosis[data$diagnosis == 1 | data$diagnosis == 2] = "GPA"
data$diagnosis[data$diagnosis == 3 | data$diagnosis == 4] = "MPA"
data$diagnosis[data$diagnosis == 5] = "EGPA"
data$diagnosis[data$diagnosis == 7 | data$diagnosis == 8] = "PAN"

### Cohen’s Kappa: Agreement with MD Diagnosis ###
# Compute Cohen's Kappa for EMAclass vs MD diagnosis
kappa_ema_vs_md <- kappa2(data.frame(data$diagnosis, data$EMAclass))

data$dcvas_result3[data$dcvas_result3 == "uncl AAV" 
                   | data$dcvas_result3 == "GPA/MPA"] = "unclassified"

data$diagnosis[data$diagnosis == "PAN"] = "unclassified"

# Compute Cohen's Kappa for dcvas_result3 vs MD diagnosis
kappa_dcvas_vs_md <- kappa2(data.frame(data$diagnosis, data$dcvas_result3))

data$EMAclass[data$EMAclass == "PAN"] = "unclassified"

# Compute Cohen's Kappa for EMAclass vs dcvas_result3
kappa_ema_vs_dcvas <- kappa2(data.frame(data$EMAclass, data$dcvas_result3))
# Print Kappa results
print(kappa_ema_vs_md)     # EMAclass vs MD Diagnosis
print(kappa_dcvas_vs_md)   # dcvas_result3 vs MD Diagnosis
print(kappa_ema_vs_dcvas)  # EMAclass vs dcvas_result3

### Fleiss’ Kappa: Overall Agreement Among All 3 Classifications ###
# Create a matrix with 3 raters
ratings <- as.matrix(data[, c("diagnosis", "EMAclass", "dcvas_result3")])
# Compute Fleiss' Kappa
fleiss_kappa <- kappam.fleiss(ratings)
# Print result
print(fleiss_kappa)

## GPA ##
### For sensitivity, specificity, PPV, NPV, accuracy

data$new_GPA_y_n_from_diagnosis <- as.integer(data$diagnosis == "GPA")
data$new_GPA_y_n_from_EMAclass <- as.integer(data$EMAclass == "GPA")
data$new_GPA_y_n_from_dcvas <- as.integer(data$dcvas_result3 == "GPA")

# Confusion Matrix for each classification system vs Gold Standard (new_GPA_y_n_from_diagnosis)
# Ensure the columns are factors (binary classification: "GPA" or "Not GPA")
data$new_GPA_y_n_from_diagnosis <- as.factor(data$new_GPA_y_n_from_diagnosis)
data$new_GPA_y_n_from_EMAclass <- as.factor(data$new_GPA_y_n_from_EMAclass)
data$new_GPA_y_n_from_dcvas <- as.factor(data$new_GPA_y_n_from_dcvas)

# Confusion matrix for new_GPA_y_n_from_EMAclass vs new_GPA_y_n_from_diagnosis (gold standard)
conf_matrix_gpa_ema_vs_diagnosis <- confusionMatrix(data$new_GPA_y_n_from_EMAclass, data$new_GPA_y_n_from_diagnosis)
# Confusion matrix for new_GPA_y_n_from_dcvas vs new_GPA_y_n_from_diagnosis (gold standard)
conf_matrix_gpa_dcvas_vs_diagnosis <- confusionMatrix(data$new_GPA_y_n_from_dcvas, data$new_GPA_y_n_from_diagnosis)

### Print Sensitivity, Specificity, PPV, NPV, Accuracy for each classification system ###
# For EMA-based classification vs gold standard
print("EMA-based GPA Classification vs Gold Standard (new_GPA_y_n_from_diagnosis):")
print(paste("Sensitivity: ", conf_matrix_gpa_ema_vs_diagnosis$byClass["Sensitivity"]))
print(paste("Specificity: ", conf_matrix_gpa_ema_vs_diagnosis$byClass["Specificity"]))
print(paste("PPV (Positive Predictive Value): ", conf_matrix_gpa_ema_vs_diagnosis$byClass["Pos Pred Value"]))
print(paste("NPV (Negative Predictive Value): ", conf_matrix_gpa_ema_vs_diagnosis$byClass["Neg Pred Value"]))
print(paste("Accuracy: ", conf_matrix_gpa_ema_vs_diagnosis$overall["Accuracy"]))

# For DCVAS-based classification vs gold standard
print("DCVAS-based GPA Classification vs Gold Standard (new_GPA_y_n_from_diagnosis):")
print(paste("Sensitivity: ", conf_matrix_gpa_dcvas_vs_diagnosis$byClass["Sensitivity"]))
print(paste("Specificity: ", conf_matrix_gpa_dcvas_vs_diagnosis$byClass["Specificity"]))
print(paste("PPV (Positive Predictive Value): ", conf_matrix_gpa_dcvas_vs_diagnosis$byClass["Pos Pred Value"]))
print(paste("NPV (Negative Predictive Value): ", conf_matrix_gpa_dcvas_vs_diagnosis$byClass["Neg Pred Value"]))
print(paste("Accuracy: ", conf_matrix_gpa_dcvas_vs_diagnosis$overall["Accuracy"]))

# ROC for EMA-based classification vs Gold Standard (new_GPA_y_n_from_diagnosis)
# roc_ema <- roc(data$new_GPA_y_n_from_diagnosis, as.numeric(as.character(data$new_GPA_y_n_from_EMAclass)))
roc_ema_gpa <- roc(data$new_GPA_y_n_from_diagnosis, as.numeric(data$new_GPA_y_n_from_EMAclass))
# ROC for DCVAS-based classification vs Gold Standard (new_GPA_y_n_from_diagnosis)
roc_dcvas_gpa <- roc(data$new_GPA_y_n_from_diagnosis, as.numeric(data$new_GPA_y_n_from_dcvas))


roc_ema_gpa <- roc(data$new_GPA_y_n_from_diagnosis, as.numeric(data$new_GPA_y_n_from_EMAclass), n.interp = 500)
roc_dcvas_gpa <- roc(data$new_GPA_y_n_from_diagnosis, as.numeric(data$new_GPA_y_n_from_dcvas), n.interp = 500)


# Extract ROC data for ggplot
roc_data_ema_gpa <- data.frame(
  fpr = 1 - roc_ema_gpa$specificities,
  tpr = roc_ema_gpa$sensitivities,
  model = "EMA-based"
)

roc_data_dcvas_gpa <- data.frame(
  fpr = 1 - roc_dcvas_gpa$specificities,
  tpr = roc_dcvas_gpa$sensitivities,
  model = "DCVAS-based"
)

# Combine data for both models
roc_data <- rbind(roc_data_ema_gpa, roc_data_dcvas_gpa)

# Plot the ROC curve using ggplot2
roc_curve_gpa <- ggplot(roc_data, aes(x = fpr, y = tpr, color = model)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # Diagonal line (random classifier)
  labs(title = "ROC Curve for GPA Classification Systems",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)",
       color = "Model") +
  theme_minimal()

print(roc_curve_gpa)

# Print AUC values (Higher AUC indicates better overall performance)
print(paste("AUC for EMA-based classification: ", auc(roc_ema_gpa)))
print(paste("AUC for DCVAS-based classification: ", auc(roc_dcvas_gpa)))

## MPA ##
### For sensitivity, specificity, PPV, NPV, accuracy
data$new_MPA_y_n_from_diagnosis <- as.integer(data$diagnosis == "MPA")
data$new_MPA_y_n_from_EMAclass <- as.integer(data$EMAclass == "MPA")
data$new_MPA_y_n_from_dcvas <- as.integer(data$dcvas_result3 == "MPA")

# Confusion Matrix for each classification system vs Gold Standard (new_GPA_y_n_from_diagnosis)
# Ensure the columns are factors (binary classification: "GPA" or "Not GPA")
data$new_MPA_y_n_from_diagnosis <- as.factor(data$new_MPA_y_n_from_diagnosis)
data$new_MPA_y_n_from_EMAclass <- as.factor(data$new_MPA_y_n_from_EMAclass)
data$new_MPA_y_n_from_dcvas <- as.factor(data$new_MPA_y_n_from_dcvas)

# Confusion matrix for new_GPA_y_n_from_EMAclass vs new_GPA_y_n_from_diagnosis (gold standard)
conf_matrix_mpa_ema_vs_diagnosis <- confusionMatrix(data$new_MPA_y_n_from_EMAclass, data$new_MPA_y_n_from_diagnosis)
# Confusion matrix for new_GPA_y_n_from_dcvas vs new_GPA_y_n_from_diagnosis (gold standard)
conf_matrix_mpa_dcvas_vs_diagnosis <- confusionMatrix(data$new_MPA_y_n_from_dcvas, data$new_MPA_y_n_from_diagnosis)

### Print Sensitivity, Specificity, PPV, NPV, Accuracy for each classification system ###
# For EMA-based classification vs gold standard
print("EMA-based MPA Classification vs Gold Standard (new_MPA_y_n_from_diagnosis):")
print(paste("Sensitivity: ", conf_matrix_mpa_ema_vs_diagnosis$byClass["Sensitivity"]))
print(paste("Specificity: ", conf_matrix_mpa_ema_vs_diagnosis$byClass["Specificity"]))
print(paste("PPV (Positive Predictive Value): ", conf_matrix_mpa_ema_vs_diagnosis$byClass["Pos Pred Value"]))
print(paste("NPV (Negative Predictive Value): ", conf_matrix_mpa_ema_vs_diagnosis$byClass["Neg Pred Value"]))
print(paste("Accuracy: ", conf_matrix_mpa_ema_vs_diagnosis$overall["Accuracy"]))

# For DCVAS-based classification vs gold standard
print("DCVAS-based MPA Classification vs Gold Standard (new_MPA_y_n_from_diagnosis):")
print(paste("Sensitivity: ", conf_matrix_mpa_dcvas_vs_diagnosis$byClass["Sensitivity"]))
print(paste("Specificity: ", conf_matrix_mpa_dcvas_vs_diagnosis$byClass["Specificity"]))
print(paste("PPV (Positive Predictive Value): ", conf_matrix_mpa_dcvas_vs_diagnosis$byClass["Pos Pred Value"]))
print(paste("NPV (Negative Predictive Value): ", conf_matrix_mpa_dcvas_vs_diagnosis$byClass["Neg Pred Value"]))
print(paste("Accuracy: ", conf_matrix_mpa_dcvas_vs_diagnosis$overall["Accuracy"]))

# ROC for EMA-based classification vs Gold Standard (new_GPA_y_n_from_diagnosis)
# roc_ema <- roc(data$new_GPA_y_n_from_diagnosis, as.numeric(as.character(data$new_GPA_y_n_from_EMAclass)))
roc_ema_mpa <- roc(data$new_MPA_y_n_from_diagnosis, as.numeric(data$new_MPA_y_n_from_EMAclass))
# ROC for DCVAS-based classification vs Gold Standard (new_GPA_y_n_from_diagnosis)
roc_dcvas_mpa <- roc(data$new_MPA_y_n_from_diagnosis, as.numeric(data$new_MPA_y_n_from_dcvas))


roc_ema_mpa <- roc(data$new_MPA_y_n_from_diagnosis, as.numeric(data$new_MPA_y_n_from_EMAclass), n.interp = 500)
roc_dcvas_mpa <- roc(data$new_MPA_y_n_from_diagnosis, as.numeric(data$new_MPA_y_n_from_dcvas), n.interp = 500)


# Extract ROC data for ggplot
roc_data_ema_mpa <- data.frame(
  fpr = 1 - roc_ema_mpa$specificities,
  tpr = roc_ema_mpa$sensitivities,
  model = "EMA-based"
)

roc_data_dcvas_mpa <- data.frame(
  fpr = 1 - roc_dcvas_mpa$specificities,
  tpr = roc_dcvas_mpa$sensitivities,
  model = "DCVAS-based"
)

# Combine data for both models
roc_data <- rbind(roc_data_ema_mpa, roc_data_dcvas_mpa)

# Plot the ROC curve using ggplot2
roc_curve_MPA <- ggplot(roc_data, aes(x = fpr, y = tpr, color = model)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # Diagonal line (random classifier)
  labs(title = "ROC Curve for MPA Classification Systems",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)",
       color = "Model") +
  theme_minimal()

print(roc_curve_MPA)

# Print AUC values (Higher AUC indicates better overall performance)
print(paste("AUC for EMA-based classification (MPA): ", auc(roc_ema_mpa)))
print(paste("AUC for DCVAS-based classification (MPA): ", auc(roc_dcvas_mpa)))


#### Additional Manuscript Stats #####

data <- read.csv("New classification sheet-USE THIS.csv", header = TRUE)

A <- data[which(data$dcvas_result3 == "GPA"
                | data$dcvas_result3 == "GPA/MPA"),]

B <- data[which(data$EMAclass == "GPA"
                & data$dcvas_result3 != "GPA"
                & data$dcvas_result3 != "GPA/MPA"),]

C <- data[which(data$EMAclass == "GPA"
                & data$dcvas_result3 != "GPA"
                & data$dcvas_result3 == "MPA" 
                & data$dcvas_result3 != "GPA/MPA"),]


D <- data[which(data$EMAclass == "GPA"
                & data$dcvas_result3 != "GPA"
                & data$dcvas_result3 == "uncl AAV" 
                & data$dcvas_result3 != "GPA/MPA"),]

E <- data[which(data$EMAclass == "GPA"),]

Fnum <- data[which((data$dcvas_result3 == "GPA" | data$dcvas_result3 == "GPA/MPA") 
                   & data$EMAclass != "GPA"),]

G <- data[which(data$dcvas_result3 == "GPA/MPA" 
                & data$EMAclass != "GPA"),]

H <- data[which(data$EMAclass != "GPA"
                & data$dcvas_result3 != "GPA"
                & data$diagnosis == "9"
                & data$dcvas_result3 != "GPA/MPA"),]

nlevels(as.factor(data$redcap_data_access_group))

data$diagnosis[data$diagnosis == "1"] = "GPA"

data$diagnosis[data$diagnosis == "2"] = "Lim GPA"

data$diagnosis[data$diagnosis == "3"|data$diagnosis == "4"] = "MPA"

data$diagnosis[data$diagnosis == "5"] = "EGPA"

data$diagnosis[data$diagnosis == "6"] = "ANCA+ Pauci-Imm GN"

data$diagnosis[data$diagnosis == "7"] = "PAN"

data$diagnosis[data$diagnosis == "8"] = "Cutaneous PAN"

data$diagnosis[data$diagnosis == "9"] = "Unc AAV"

data$diagnosis[data$diagnosis == "12"] = "Other"

describeFactors(data$diagnosis)

describeFactors(data$EMAclass)

describeFactors(data$dcvas_result3)

describeFactors(data$mpo_result)

describeFactors(data$pr3_result)

ACR_GPA_noPR3_noMPO <- nrow(data[which(data$dcvas_result3 == "GPA"
                                  & data$pr3_result == 1
                                  & data$mpo_result == 1),])

ACR_MPA_noPR3_noMPO <- nrow(data[which(data$dcvas_result3 == "MPA"
                                  & data$pr3_result == "1"
                                  & data$mpo_result == "1"),])

EMA_GPA_noPR3_noMPO <- nrow(data[which(data$EMAclass == "GPA"
                                  & data$pr3_result == "1"
                                  & data$mpo_result == "1"),])

EMA_MPA_noPR3_noMPO <- nrow(data[which(data$EMAclass == "MPA"
                                  & data$pr3_result == "1"
                                  & data$mpo_result == "1"),])

ACR_GPA_miss_PR3_MPO <- nrow(data[which(data$dcvas_result3 == "GPA"
                                       & is.na(data$pr3_result)
                                       & is.na(data$mpo_result)),])

ACR_MPA_miss_PR3_MPO <- nrow(data[which(data$dcvas_result3 == "MPA"
                                       & is.na(data$pr3_result)
                                       & is.na(data$mpo_result)),])

EMA_GPA_miss_PR3_MPO <- nrow(data[which(data$EMAclass == "GPA"
                                       & is.na(data$pr3_result)
                                       & is.na(data$mpo_result)),])

EMA_MPA_miss_PR3_MPO <- nrow(data[which(data$EMAclass == "MPA"
                                       & is.na(data$pr3_result)
                                       & is.na(data$mpo_result)),])

