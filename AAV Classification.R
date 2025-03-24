#########################################
#### CLASSIFICATION ANCA-ASSOC VASC  ####
####   Sela Grays - Jan 27 2025      ####
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


#### LOAD NEW CLASSIFICATION SHEET ####
data <- read.csv("New classification sheet-USE THIS.csv", header = TRUE)

data$gender <- as.numeric(data$gender)

data$date_enroll <- as.Date(data$date_enroll, "%m/%d/%Y")
date_enroll_no_NA <- data[-which(is.na(data$date_enroll)),]

data$pvas_chest___3[data$lungs_s_node == 1|data$lungs_cav == 1|data$lungs_m_node == 1] = 1

data$pvas_chest___5[data$lungs_mig_infil == 1|data$lungs_fix_infil == 1] = 1

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


data$ren[data$pvas_renal___1 == 1 | 
           data$pvas_renal___2 == 1 | 
           data$pvas_renal___3 == 1 |
           data$pvas_renal___4 == 1 | 
           data$pvas_renal___5 == 1 | 
           data$pvas_renal___6 == 1 | 
           data$pvas_renal___8 == 1 ] = 1

data$gi[data$pvas_gi___1 == 1 | 
          data$pvas_gi___2 == 1 | 
          data$pvas_gi___3 == 1 | 
          data$pvas_gi___4 == 1] = 1

data$cardio[data$pvas_cv___1 == 1 | 
              data$pvas_cv___2 == 1 | 
              data$pvas_cv___3 == 1 |
              data$pvas_cv___4 == 1 | 
              data$pvas_cv___5 == 1 | 
              data$pvas_cv___6 == 1 | 
              data$pvas_cv___7 == 1 |
              data$pvas_cv___8 == 1 |
              data$pvas_cv___9 == 1] = 1

data$chest[data$pvas_ent___1 == 1 | 
             data$pvas_ent___2 == 1 | 
             data$pvas_ent___3 == 1 | 
             data$pvas_ent___4 == 1 | 
             data$pvas_ent___5 == 1] = 1

data$ent[data$pvas_ent___1 == 1 | 
           data$pvas_ent___2 == 1 | 
           data$pvas_ent___3 == 1 | 
           data$pvas_ent___4 == 1 | 
           data$pvas_ent___5 == 1] = 1

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

data$cuta[data$pvas_cut2___1 == 1 | 
            data$pvas_cut2___2 == 1 | 
            data$pvas_cut2___3 == 1 | 
            data$pvas_cut2___4 == 1 |
            data$pvas_cut2___5 == 1 | 
            data$pvas_cut2___6 == 1 | 
            data$pvas_cut2___7 == 1 | 
            data$pvas_cut2___8 == 1 |
            data$pvas_cut2___9 == 1] = 1

data$gen[data$pvas_gen1___1 == 1 | 
           data$pvas_gen1___2 == 1 | 
           data$pvas_gen1___3 == 1 | 
           data$pvas_gen1___4 == 1] = 1


#### SEPERATE DIAGNOSIS ####

not_ACR <- data[which(data$dcvas_result3 == "" 
                      | data$dcvas_result3 == "EGPA"
                      | data$dcvas_result3 == "uncl AAV"
                      | data$dcvas_result3 == "unclassified"),]

all_ACR <- data[which(data$dcvas_result3 != "" 
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

GPA_EMA <- data[which(data$EMAclass == "GPA"),]

MPA <- data[which(data$EMAclass == "MPA"|data$dcvas_result3 == "MPA"|
                    data$diagnosis == 3 | data$diagnosis == 4),]

MPA_inc <-data[which(data$EMAclass == "MPA"|data$dcvas_result3 == "MPA"|data$dcvas_result3 == "GPA/MPA"|
                       data$diagnosis == 3 | data$diagnosis == 4),]

MPA_ACR <- data[which(data$dcvas_result3 == "MPA"),]

MPA_EMA <- data[which(data$EMAclass == "GPA"),]

both_ACR <- data[which(data$dcvas_result3 == "GPA/MPA" ),]

UNC_EMA <- data[which(data$EMAclass == "GPA"),]

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
describeFactors(all_ACR$dcvas_result3)
#### SUPP Table 1 ####


describeFactors(GPA_ACR$gen)
describeFactors(MPA_ACR$gen)
describeFactors(UNC_ACR$gen)

#Myalgia
describeFactors(GPA_ACR$pvas_gen1___1)
describeFactors(MPA_ACR$pvas_gen1___1)
describeFactors(UNC_ACR$pvas_gen1___1)
# glm(pvas_gen1___1 ~ dcvas_result3, data = all_ACR, family="binomial")
# summary(glm(pvas_gen1___1 ~ dcvas_result3, data = all_ACR, family="binomial"))
# exp(coef(glm(pvas_gen1___1 ~ dcvas_result3, data = all_ACR, family="binomial")))

#Arthralgia or arthritis
describeFactors(GPA_ACR$pvas_gen1___2)
describeFactors(MPA_ACR$pvas_gen1___2)
describeFactors(UNC_ACR$pvas_gen1___2)

#Fever >= 38 deg C
describeFactors(GPA_ACR$pvas_gen1___3)
describeFactors(MPA_ACR$pvas_gen1___3)
describeFactors(UNC_ACR$pvas_gen1___3)

#Weight loss >= 5% body weight
describeFactors(GPA_ACR$pvas_gen1___4)
describeFactors(MPA_ACR$pvas_gen1___4)
describeFactors(UNC_ACR$pvas_gen1___4)


describeFactors(GPA_ACR$cuta)
describeFactors(MPA_ACR$cuta)
describeFactors(UNC_ACR$cuta)

#Polymorphous exanthema
describeFactors(GPA_ACR$pvas_cut2___1)
describeFactors(MPA_ACR$pvas_cut2___1)
describeFactors(UNC_ACR$pvas_cut2___1)

#Livedo
describeFactors(GPA_ACR$pvas_cut2___2)
describeFactors(MPA_ACR$pvas_cut2___2)
describeFactors(UNC_ACR$pvas_cut2___2)

#Panniculitis
describeFactors(GPA_ACR$pvas_cut2___3)
describeFactors(MPA_ACR$pvas_cut2___3)
describeFactors(UNC_ACR$pvas_cut2___3)

#Purpura
describeFactors(GPA_ACR$pvas_cut2___4)
describeFactors(MPA_ACR$pvas_cut2___4)
describeFactors(UNC_ACR$pvas_cut2___4)

#Skin nodules
describeFactors(GPA_ACR$pvas_cut2___5)
describeFactors(MPA_ACR$pvas_cut2___5)
describeFactors(UNC_ACR$pvas_cut2___5)

#Infarct (nail edge lesion, splinter haemorrhage)
describeFactors(GPA_ACR$pvas_cut2___6)
describeFactors(MPA_ACR$pvas_cut2___6)
describeFactors(UNC_ACR$pvas_cut2___6)

#Ulcer (full-thickness necrosis)
describeFactors(GPA_ACR$pvas_cut2___7)
describeFactors(MPA_ACR$pvas_cut2___7)
describeFactors(UNC_ACR$pvas_cut2___7)

#Gangrene (extensive necrosis)
describeFactors(GPA_ACR$pvas_cut2___8)
describeFactors(MPA_ACR$pvas_cut2___8)
describeFactors(UNC_ACR$pvas_cut2___8)

#Other skin vasculitis (specify below)
describeFactors(GPA_ACR$pvas_cut2___9)
describeFactors(MPA_ACR$pvas_cut2___9)
describeFactors(UNC_ACR$pvas_cut2___9)



describeFactors(GPA_ACR$muc)
describeFactors(MPA_ACR$muc)
describeFactors(UNC_ACR$muc)

#Mouth ulcers/granulomata
describeFactors(GPA_ACR$pvas_mucous2___1)
describeFactors(MPA_ACR$pvas_mucous2___1)
describeFactors(UNC_ACR$pvas_mucous2___1)

#Genital ulcers
describeFactors(GPA_ACR$pvas_mucous2___2)
describeFactors(MPA_ACR$pvas_mucous2___2)
describeFactors(UNC_ACR$pvas_mucous2___2)

#Adnexal inflammation (sialadenitis/ dacryocystitis)
describeFactors(GPA_ACR$pvas_mucous2___3)
describeFactors(MPA_ACR$pvas_mucous2___3)
describeFactors(UNC_ACR$pvas_mucous2___3)

#Significant proptosis
describeFactors(GPA_ACR$pvas_mucous2___4)
describeFactors(MPA_ACR$pvas_mucous2___4)
describeFactors(UNC_ACR$pvas_mucous2___4)

#Red eye (Epi)scleritis
describeFactors(GPA_ACR$pvas_mucous2___5)
describeFactors(MPA_ACR$pvas_mucous2___5)
describeFactors(UNC_ACR$pvas_mucous2___5)

#Red eye conjunctivitis/ blepharitis/keratitis
describeFactors(GPA_ACR$pvas_mucous2___6)
describeFactors(MPA_ACR$pvas_mucous2___6)
describeFactors(UNC_ACR$pvas_mucous2___6)

#Uveitis
describeFactors(GPA_ACR$pvas_mucous2___7)
describeFactors(MPA_ACR$pvas_mucous2___7)
describeFactors(UNC_ACR$pvas_mucous2___7)

#Blurred vision
describeFactors(GPA_ACR$pvas_mucous2___8)
describeFactors(MPA_ACR$pvas_mucous2___8)
describeFactors(UNC_ACR$pvas_mucous2___8)

#Sudden visual loss
describeFactors(GPA_ACR$pvas_mucous2___9)
describeFactors(MPA_ACR$pvas_mucous2___9)
describeFactors(UNC_ACR$pvas_mucous2___9)

#Retinal vasculitis/retinal vessel thrombosis/retinal exudates/haemorrhages
describeFactors(GPA_ACR$pvas_mucous2___10)
describeFactors(MPA_ACR$pvas_mucous2___10)
describeFactors(UNC_ACR$pvas_mucous2___10)




describeFactors(GPA_ACR$ent)
describeFactors(MPA_ACR$ent)
describeFactors(UNC_ACR$ent)

#Bloody nasal discharge/crusts/ulcers/granuloma
describeFactors(GPA_ACR$pvas_ent___1)
describeFactors(MPA_ACR$pvas_ent___1)
describeFactors(UNC_ACR$pvas_ent___1)

#Paranasal sinus involvement
describeFactors(GPA_ACR$pvas_ent___2)
describeFactors(MPA_ACR$pvas_ent___2)
describeFactors(UNC_ACR$pvas_ent___2)

#ubglottic stenosis/ hoarseness /stridor
describeFactors(GPA_ACR$pvas_ent___3)
describeFactors(MPA_ACR$pvas_ent___3)
describeFactors(UNC_ACR$pvas_ent___3)

#Conductive hearing loss
describeFactors(GPA_ACR$pvas_ent___4)
describeFactors(MPA_ACR$pvas_ent___4)
describeFactors(UNC_ACR$pvas_ent___4)

#Sensorineural hearing loss
describeFactors(GPA_ACR$pvas_ent___5)
describeFactors(MPA_ACR$pvas_ent___5)
describeFactors(UNC_ACR$pvas_ent___5)




describeFactors(GPA_ACR$chest)
describeFactors(MPA_ACR$chest)
describeFactors(UNC_ACR$chest)

#Wheeze or expiratory dyspnea
describeFactors(GPA_ACR$pvas_chest___1)
describeFactors(MPA_ACR$pvas_chest___1)
describeFactors(UNC_ACR$pvas_chest___1)

#Endobronchial/endotracheal involvement
describeFactors(GPA_ACR$pvas_chest___2)
describeFactors(MPA_ACR$pvas_chest___2)
describeFactors(UNC_ACR$pvas_chest___2)

#Nodules or cavities 
describeFactors(GPA_ACR$pvas_chest___3)
describeFactors(MPA_ACR$pvas_chest___3)
describeFactors(UNC_ACR$pvas_chest___3)

#Massive hemoptysis/alveolar hemorrhage 
describeFactors(GPA_ACR$pvas_chest___4)
describeFactors(MPA_ACR$pvas_chest___4)
describeFactors(UNC_ACR$pvas_chest___4)

#Respiratory failure
describeFactors(GPA_ACR$pvas_chest___5)
describeFactors(MPA_ACR$pvas_chest___5)
describeFactors(UNC_ACR$pvas_chest___5)



describeFactors(GPA_ACR$cardio)
describeFactors(MPA_ACR$cardio)
describeFactors(UNC_ACR$cardio)

#Loss of pulses
describeFactors(GPA_ACR$pvas_cv___1)
describeFactors(MPA_ACR$pvas_cv___1)
describeFactors(UNC_ACR$pvas_cv___1)

#Bruits over accessible arteries
describeFactors(GPA_ACR$pvas_cv___2)
describeFactors(MPA_ACR$pvas_cv___2)
describeFactors(UNC_ACR$pvas_cv___2)

#Blood pressure discrepancy
describeFactors(GPA_ACR$pvas_cv___3)
describeFactors(MPA_ACR$pvas_cv___3)
describeFactors(UNC_ACR$pvas_cv___3)

#Claudication of extremities
describeFactors(GPA_ACR$pvas_cv___4)
describeFactors(MPA_ACR$pvas_cv___4)
describeFactors(UNC_ACR$pvas_cv___4)

#Ischaemic cardiac pain
describeFactors(GPA_ACR$pvas_cv___5)
describeFactors(MPA_ACR$pvas_cv___5)
describeFactors(UNC_ACR$pvas_cv___5)

#Cardiomyopathy
describeFactors(GPA_ACR$pvas_cv___6)
describeFactors(MPA_ACR$pvas_cv___6)
describeFactors(UNC_ACR$pvas_cv___6)

#Congestive cardiac failure
describeFactors(GPA_ACR$pvas_cv___7)
describeFactors(MPA_ACR$pvas_cv___7)
describeFactors(UNC_ACR$pvas_cv___7)

#Valvular heart disease
describeFactors(GPA_ACR$pvas_cv___8)
describeFactors(MPA_ACR$pvas_cv___8)
describeFactors(UNC_ACR$pvas_cv___8)

#Pericarditis
describeFactors(GPA_ACR$pvas_cv___9)
describeFactors(MPA_ACR$pvas_cv___9)
describeFactors(UNC_ACR$pvas_cv___9)



describeFactors(GPA_ACR$gi)
describeFactors(MPA_ACR$gi)
describeFactors(UNC_ACR$gi)

#Abdominal pain
describeFactors(GPA_ACR$pvas_gi___1)
describeFactors(MPA_ACR$pvas_gi___1)
describeFactors(UNC_ACR$pvas_gi___1)

#Peritonitis
describeFactors(GPA_ACR$pvas_gi___2)
describeFactors(MPA_ACR$pvas_gi___2)
describeFactors(UNC_ACR$pvas_gi___2)

#Blood in stools or bloody diarrhea
describeFactors(GPA_ACR$pvas_gi___3)
describeFactors(MPA_ACR$pvas_gi___3)
describeFactors(UNC_ACR$pvas_gi___3)

#Bowel ischaemia
describeFactors(GPA_ACR$pvas_gi___4)
describeFactors(MPA_ACR$pvas_gi___4)
describeFactors(UNC_ACR$pvas_gi___4)


describeFactors(GPA_ACR$ren)
describeFactors(MPA_ACR$ren)
describeFactors(UNC_ACR$ren)

#Hypertension >95th centile (for height)
describeFactors(GPA_ACR$pvas_renal___1)
describeFactors(MPA_ACR$pvas_renal___1)
describeFactors(UNC_ACR$pvas_renal___1)

#Proteinuria >0.3 g/24h,>20mg/mmol creatinine
describeFactors(GPA_ACR$pvas_renal___2)
describeFactors(MPA_ACR$pvas_renal___2)
describeFactors(UNC_ACR$pvas_renal___2)

#Haematuria >= 2+ or 5 rbc/hpf or red cell casts
describeFactors(GPA_ACR$pvas_renal___3)
describeFactors(MPA_ACR$pvas_renal___3)
describeFactors(UNC_ACR$pvas_renal___3)

#GFR 50-80ml/min/1.73 m2 *s
describeFactors(GPA_ACR$pvas_renal___4)
describeFactors(MPA_ACR$pvas_renal___4)
describeFactors(UNC_ACR$pvas_renal___4)

#GFR 15-49 ml/min/1.73 m2 *
describeFactors(GPA_ACR$pvas_renal___5)
describeFactors(MPA_ACR$pvas_renal___5)
describeFactors(UNC_ACR$pvas_renal___5)

#GFR < 15 ml/min/1.73m2 *
describeFactors(GPA_ACR$pvas_renal___6)
describeFactors(MPA_ACR$pvas_renal___6)
describeFactors(UNC_ACR$pvas_renal___6)

#Rise in creatinine > 10% or Creatinine clearance (GFR) fall > 25%
describeFactors(GPA_ACR$pvas_renal___8)
describeFactors(MPA_ACR$pvas_renal___8)
describeFactors(UNC_ACR$pvas_renal___8)



describeFactors(GPA_ACR$cns)
describeFactors(MPA_ACR$cns)
describeFactors(UNC_ACR$cns)

#Headache
describeFactors(GPA_ACR$pvas_cv___1)
describeFactors(MPA_ACR$pvas_cv___1)
describeFactors(UNC_ACR$pvas_cv___1)

#Meningitis/encephalitis
describeFactors(GPA_ACR$pvas_cv___2)
describeFactors(MPA_ACR$pvas_cv___2)
describeFactors(UNC_ACR$pvas_cv___2)

#Organic confusion/cognitive dysfunction
describeFactors(GPA_ACR$pvas_cv___3)
describeFactors(MPA_ACR$pvas_cv___3)
describeFactors(UNC_ACR$pvas_cv___3)

#Seizures (not hypertensive)
describeFactors(GPA_ACR$pvas_cv___4)
describeFactors(MPA_ACR$pvas_cv___4)
describeFactors(UNC_ACR$pvas_cv___4)

#Stroke
describeFactors(GPA_ACR$pvas_cv___5)
describeFactors(MPA_ACR$pvas_cv___5)
describeFactors(UNC_ACR$pvas_cv___5)

#Cord lesion
describeFactors(GPA_ACR$pvas_cv___6)
describeFactors(MPA_ACR$pvas_cv___6)
describeFactors(UNC_ACR$pvas_cv___6)

#Cranial nerve palsy
describeFactors(GPA_ACR$pvas_cv___7)
describeFactors(MPA_ACR$pvas_cv___7)
describeFactors(UNC_ACR$pvas_cv___7)

#Sensory peripheral neuropathy
describeFactors(GPA_ACR$pvas_cv___8)
describeFactors(MPA_ACR$pvas_cv___8)
describeFactors(UNC_ACR$pvas_cv___8)

#Motor mononeuritis multiplex
describeFactors(GPA_ACR$pvas_cv___9)
describeFactors(MPA_ACR$pvas_cv___9)
describeFactors(UNC_ACR$pvas_cv___9)

#### SUPP Table 2 ####

describeFactors(GPA_EMA$pvas_gen1___1)
describeFactors(MPA_EMA$pvas_gen1___1)
describeFactors(UNC_EMA$pvas_gen1___1)

#### FIG 4a ####

# df2 <- data.frame(supp=rep(c("Ankara", "ACR/Eular"), each=3),
#                   dose=rep(c("General", "Cutaenous", "Mucous", "Cardiovascular", "ENT", "Chest", "GI", "Renal", "CNS", "PR3", "MPO"),2),
#                   len=c(6.8, 15, 33, 4.2, 10, 29.5))
# 
# ggplot(data=df2, aes(x=dose, y=len, fill=supp)) + geom_bar(stat="identity", position=position_dodge())
# 

#### FIG 4b ####
#### FIG 4c ####
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

UNC_ACR <- data[which(data$dcvas_result3 == "unclassified"),]

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

