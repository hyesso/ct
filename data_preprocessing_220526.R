#### scan landmark modeling

## load packages ---------------------------------------------------------
pkg_list<-c("readxl", "writexl", "dplyr", "scatterplot3d", "plotly", "rockchalk", "summarytools","tidyverse", "cowplot", "corrplot", "car","DT",
            "sjPlot", "sjmisc", "sjlabelled", "htmltools", "kableExtra",  # for tab_model
            "plyr", "devtools", "Table1", 
            "caret", "fastDummies", "glmnet") #for regression
pacman::p_load(pkg_list, install=TRUE, update=FALSE, character.only = TRUE)

## read_data -------------------------------------------------------------
file_list <- list.files("/Users/jinhyesu/my_project/rawData/220511_ct/ct_data_file")
df_fin <- data.frame()

for(i in 1:length(file_list)){
  f_path <- paste0("/Users/jinhyesu/my_project/rawData/220511_ct/ct_data_file/", file_list[[i]], "/landmark_result.csv")
  f_name <- substr(file_list[[i]], 1, 3)
  df<- read.csv(f_path, header=FALSE)
  
  df_pre <- df[c(2:26),]
  df_pre <- df_pre %>% mutate(state = rep("pre", 25), 
                    case = rep(f_name, 25))
  colnames(df_pre) <- c("p_num", "x", "y", "z", "state", "case")
  
  df_post <- df[c(28:52),]
  df_post <- df_post %>% mutate(state = rep("post", 25), 
                              case = rep(f_name, 25))
  colnames(df_post) <- c("p_num", "x", "y", "z", "state", "case")
  
  df_diff <- data.frame(p_num=df_post$p_num, x=df_pre$x-df_post$x, y=df_pre$y-df_post$y, z=df_pre$z-df_post$z, state=rep("diff", 25), case=rep(f_name, 25))

  df_fin <-bind_rows(df_fin, df_pre, df_post, df_diff)
}

sub_diff<- df_fin[which(df_fin$state=="diff"),]

scatterplot3d(sub_diff[,2:4])

sub_pre <- df_fin[which(df_fin$state=="pre" & df_fin$case =="263"),]

scatterplot3d(sub_pre[c(1,6,11,16,21),])
scatterplot3d(sub_pre[c(2,7,12,17,22),])
scatterplot3d(sub_pre[c(3,8,13,18,23),])
scatterplot3d(sub_pre[c(4,9,14,19,24),])
scatterplot3d(sub_pre[c(5,10,15,20,25),])

test_pre<- sub_pre[c(4,9,14,19,24),]

p <- plot_ly(test_pre, x = test_pre$x, y = test_pre$y, z = test_pre$z)
               +              #color = iris$Species, colors = c('#BF382A', '#0C4B8E')) %>%
  +   add_markers() %>%
  +   layout(scene = list(xaxis = list(title = 'Sepal.Length'), #x축 제목설정
                          +                       yaxis = list(title = 'Sepal.Width'),  #y축 제목설정
                          +                       zaxis = list(title = 'Petal.Length'))) #z축 제목설정

p <- plot_ly(sub_pre, x = sub_pre$x, y = sub_pre$y, z = sub_pre$z)
## pre model 
# 3D scatter plot
s3d <- scatterplot3d(sub_pre[,2:4], type = "h", color = "blue",
                     angle=55, pch = 16)
# Add regression plane ; line
my.lm <- lm(sub_pre$z ~ sub_pre$x + sub_pre$y)
s3d$plane3d(my.lm)
# Add supplementary points
s3d$points3d(seq(10, 20, 2), seq(85, 60, -5), seq(60, 10, -10),
             col = "red", type = "h", pch = 8)

# 2 dimension
my.lm <- lm(sub_pre$z ~ poly(sub_pre$x,2) + poly(sub_pre$y,2))

# Add supplementary points
plotPlane(my.lm, sub_pre$x, sub_pre$y, pch=16, drawArrows=TRUE, 
          acol="red")

sub_post <- df_fin[which(df_fin$state=="post" & df_fin$case =="263"),]
s3d <- scatterplot3d(sub_post[,2:4], type = "h", color = "blue",
                     angle=55, pch = 16)


test_pre<- sub_pre[c(1,2,3,4,5),]
my.lm <- lm(test_pre$z ~ poly(test_pre$x,2) + poly(test_pre$y,2))

x <- seq(-150, 100, by = 10)
y <- seq(140,160 , by = 2)
plane <- outer(x, y, function(a, b){my.lm$coef[1] + 
    my.lm$coef[2]*a +my.lm$coef[3]*a^2 + my.lm$coef[4]*b+ my.lm$coef[5]*b^2})

plot_ly(test_pre, x = test_pre$x, y = test_pre$y, z = test_pre$z, opacity = 0.5) %>%
  add_markers() %>%
  add_surface(x = ~x, y = ~y, z = ~plane, showscale = FALSE)

remotes::install_github("mannfred/curvr")

#' Calculate total curvature from smoothing or interpolating splines.
#'
#'
#' @param landmark_matrix is a \code{matrix} object with \code{[,1]}
#' containing the x landmark coordinates and \code{[,2]} containing
#' the y landmark coordinates.
#'
#' @param x_range the lower and upper x-value bounds to
#' calculate curvature. Concatenate the lower and upper bounds using
#' \code{c()}, E.g. for lower = 1 and upper = 10, type \code{c(1,10)}.
#'
#' @param type either 'ip' for an interpolating spline or 'smooth' for a
#' smoothing spline. Uses \code{stats::spline()} or \code{stats::smooth.spline()}, respectively,
#' for curve fitting and estimating the derivatives. Default is \code{type = 'smooth'}.
#' See: ?spline and ?smooth.spline for details.
#'
#' @return a `list` with two named elements. `$Ktot` is the total curvature in radians. `$Ki` is a numeric vector of local curvature values.
#'
#' @examples
#'
#' # a landmark matrix describing a segment of the unit circle#'
#' x <- seq(0, 1, by = 0.01)
#' y <- sqrt(1-x^2)
#' mdat <- matrix(c(x, y), nrow = 101, ncol = 2)
#'
#' # total curvature between x=0 and x=sqrt(2)/2 should be approximately pi/4
#' abs(curvature_spline(mdat, c(0, sqrt(2)/2), type='smooth')$Ktot)
#'
#' @importFrom dplyr %>%
#' @importFrom stats smooth.spline predict splinefun spline integrate
#'
#' @export

curvature_spline <- function(landmark_matrix, x_range, type = 'smooth') {
  
  # extract/separate x and y coords
  x_coords <- landmark_matrix[, 1]
  y_coords <- landmark_matrix[, 2]
  
  # check that for every x there is a y
  if (length(x_coords) != length(y_coords)) {
    stop("every x coordinate must have a corresponding y coordinate")
  }
  
  if (type == 'smooth'){
    # fit a spline to landmark coordinates
    s0 <- smooth.spline(landmark_matrix)
    
    # first deriv values
    s1 <- predict(s0, deriv = 1)
    
    # fit spline func to first deriv values
    s1func <- splinefun(x = s1$x, y = s1$y)
    
    # second deriv func (1st deriv of s1)
    s2func <- splinefun(x = s1$x, y = s1func(s1$x, deriv = 1))
    
  } else if (type == 'ip'){
    
    # compute coordinates from a cubic spline fit
    s0 <- spline(landmark_matrix)
    
    # create a spline function from coordinates
    s0func <- splinefun(s0$x, s0$y)
    
    # estimate y coords of first derivative
    s1 <- s0func(s0$x, deriv = 1)
    
    # create a function for first derivative
    s1func <- splinefun(s0$x, s1)
    
    # create a function for second derivative
    s2func <- splinefun(x = s0$x, y = s1func(s0$x, deriv = 1))
    
  } else {
    
    stop("spline type must be 'ip' or 'smooth'")
    
  }
  
  # define K * ds
  k_fun <- function(x) {
    f1 <- s1func
    f2 <- s2func
    ((f2(x)) / ((1 + (f1(x)^2))^1.5)) * (sqrt(1 + (f1(x))^2))
  }
  
  # compute integral of K*ds
  Ktot <- integrate(k_fun, lower = x_range[1], upper = x_range[2])$value
  
  # compute point-wise K from second deriv function
  Ki <- s2func(landmark_matrix[,1])
  
  curvature <- list(Ktot = Ktot, Ki = Ki)
  return(curvature)
  
}

# a landmark matrix describing a segment of the unit circle#'
x <- seq(0, 1, by = 0.01)
y <- sqrt(1-x^2)
mdat <- matrix(c(x, y), nrow = 101, ncol = 2)

# total curvature between x=0 and x=sqrt(2)/2 should be approximately pi/4
abs(curvature_spline(mdat, c(0, sqrt(2)/2), type='smooth')$Ktot)

## 220616 generate dataframe as asked
df<- read.csv("/Users/jinhyesu/Downloads/landmark_result.csv", header=FALSE)

pre_df <- df[c(2:26),]  
post_df <- df[c(28:52),]  

write.csv(pre_df, file = "/Users/jinhyesu/Desktop/df_pre.csv")
write.csv(post_df, file = "/Users/jinhyesu/Desktop/df_post.csv")
  
## 220620 b-spline surface
require(rms)  # Harrell's gift to the R world.

# Better to keep the original names and do so within a dataframe.
df <- df_fin[-which(df_fin$case %in% c("331","356", "367","380")),] #total 62 cases

pred_p <- list()
for(i in unique(df$case)){
  df_pre <- df[which(df$case == i &df$state =="pre"),]
  att <- df_pre[c('x','y','z')]
  add <- datadist(att) 
  options(datadist="add") 
  
  mdl <- ols(z ~ rcs(x,3)*rcs(y,3) ,data=att)
  pred <- Predict(mdl, 'x','y', np=30)
  
  pred_p[[i]] <- bplot(pred,lfun=wireframe)
}

df <- df_fin[which(df_fin$state == "pre"& df_fin$case == 263),]
att <- df[c('x','y','z')]
add <- datadist(att)  # records ranges and descriptive info on data
options(datadist="add")  # need these for the rms functions

mdl <- ols(z ~ rcs(x,3)*rcs(y,3) ,data=att)
pred <- Predict(mdl, 'x','y', np=30)

bplot(pred,lfun=wireframe)

att_pca <- prcomp(att, center = FALSE,scale. = FALSE)
att_pda_df<- as.data.frame(att_pca$x)

add_df <- datadist(att_pda_df)  # records ranges and descriptive info on data
options(datadist = "add_df")

mdx <- ols(PC3 ~ rcs(PC1, 3)*rcs(PC2, 3), data=att_pda_df)

predx <- Predict(mdx, 'PC1', 'PC2', np= 25)

bplot(predx,lfun=wireframe)

## using pca, extract axis
library(factoextra)

df <- df_fin[which(df_fin$state == "pre"& df_fin$case == 263),]
x_pca <- prcomp(df[c(20:24),][,c(2:4)],center = FALSE,scale. = FALSE)
str(x_pca)

fviz_pca_var(x_pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_eig(x_pca)

# Centering and scaling the supplementary individuals
ind.scaled <- scale(df[,c(2:4)], 
                    center = x_pca$center,
                    scale = x_pca$scale)
as.matrix(df[,c(2:4)]) %*% x_pca$rotation


str(x_pca$rotation[,1])

att <- as.data.frame(ind.scaled)
add <- datadist(att)  # records ranges and descriptive info on data
options(datadist="add")  # need these for the rms functions

mdl <- ols(z ~ rcs(x,3)*rcs(y,3) ,data=att)
pred <- Predict(mdl, 'x','y', np=30)

bplot(pred,lfun=wireframe)

##y axis certeralized
y_pca <- prcomp(df[c(3,8,13,18,23),][,c(2:4)],center = FALSE,scale. = FALSE)
summary(y_pca)
screeplot(y_pca)

fviz_pca_var(y_pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_eig(y_pca)

rot_mat <- cbind(x_pca$rotation[,1],y_pca$rotation[,1], c(1,1,1))

# Centering and scaling the supplementary individuals
ind.scaled <- scale(df[,c(2:4)], 
                    center = x_pca$center,
                    scale = x_pca$scale)

att <- as.data.frame(ind.scaled)
add <- datadist(att)  # records ranges and descriptive info on data
options(datadist="add")  # need these for the rms functions

mdl <- ols(z ~ rcs(x,3)*rcs(y,3) ,data=att)
pred <- Predict(mdl, 'x','y', np=30)

bplot(pred,lfun=wireframe)

### SVD 
p0<- svd(scale(att))
p0$d
p0$v # v = a matrix whose columns contain the right singular vectors of x, present if nv > 0. Dimension c(p, nv).

#recovers the original matrix
(p0$u %*% diag(p0$d) %*% t(p0$v)) %>% head
scale(att) %>%  head # it's equals with origin matrix

p0$d^2/(nrow(p0$u) -1)

eigenValues = p0$d^2/ (nrow(p0$u)-1)


svdSummary<-function(svdRes,sf=4){
  if(is(svdRes,"prcomp")){
    eigenvalue=svdRes$sdev^2
  }else{
    #d=signif(svdRes$d,sf)
    eigenvalue= svdRes$d^2/(nrow(svdRes$u) - 1)
  }
  data.frame(cbind(
    eigenvalues=signif(eigenvalue,sf),
    sd = signif(sqrt(eigenvalue),sf),
    variance.percent = paste0(signif((eigenvalue/sum(eigenvalue)),2)*100,"%"),
    cumulative.variance.percent = paste0(cumsum(signif((eigenvalue/sum(eigenvalue)),2))*100,"%")))
}


eigSum.svd <-svdSummary(p0)
eigSum.svd 

#### 220625 practice
# princomp ; 상관계수행렬 또는 공분산행렬에 대한 Eigenvalue Decomposition(고유근분해)을 이용하여 계산
# prcomp ; 공분산행렬의 eigenvalue을 이용하지 않고, 원 데이터에 대해 SVD(Singular Value Decomposition : 특이값분해)를 수행하여 계산
y_princomp <- princomp(df[,c(2:4)],cor=F)
summary(y_princomp)

str(y_princomp)
y_princomp$scores

att <- df_pre[c('x','y','z')]
att<- as.data.frame(y_princomp$scores)
add <- datadist(att) 
options(datadist="add") 
mdl <- ols(Comp.3 ~ rcs(Comp.1,3)*rcs(Comp.2,3) ,data=att)
pred <- Predict(mdl, 'Comp.1','Comp.2', np=30)

bplot(pred,lfun=wireframe)

########220627 add pca results
########220628 add post graphics
df <- df_fin[-which(df_fin$case %in% c("331","356", "367","380")),] #total 62 cases

pre_p <- list(); pre_r <- list()
pre_pca_p <- list(); pre_pca_r <- list()

post_p <- list(); post_r <- list()
post_pca_p <- list(); post_pca_r <- list()
for(i in unique(df$case)){
  # pre graph
  df_pre <- df[which(df$case == i &df$state =="pre"),]
  att <- df_pre[c('x','y','z')]
  add <- datadist(att) 
  options(datadist="add") 
  
  mdl <- ols(z ~ rcs(x,3)*rcs(y,3) ,data=att)
  pred <- Predict(mdl, 'x','y', np=30)
  
  pre_p[[i]] <- bplot(pred,lfun=wireframe)
  pre_r[[i]] <- mdl
  
  df_pca <- princomp(df_pre[,c(2:4)], cor= F) 
  att<- as.data.frame(df_pca$scores)
  add <- datadist(att) 
  options(datadist="add") 
  
  mdl <- ols(Comp.3 ~ rcs(Comp.1,3)*rcs(Comp.2,3) ,data=att)
  pred <- Predict(mdl, 'Comp.1','Comp.2', np=30)
  
  pre_pca_p[[i]] <-bplot(pred,lfun=wireframe)
  pre_pca_r[[i]] <- mdl
  
  # post graph
  df_post <- df[which(df$case == i &df$state =="post"),]
  att <- df_post[c('x','y','z')]
  add <- datadist(att) 
  options(datadist="add") 
  
  mdl <- ols(z ~ rcs(x,3)*rcs(y,3) ,data=att)
  
  pred <- Predict(mdl, 'x','y', np=30)
  
  post_p[[i]] <- bplot(pred,lfun=wireframe)
  post_r[[i]] <- mdl
  
  df_pca <- princomp(df_post[,c(2:4)], cor= F) 
  att<- as.data.frame(df_pca$scores)
  add <- datadist(att) 
  options(datadist="add") 
  
  mdl <- ols(Comp.3 ~ rcs(Comp.1,3)*rcs(Comp.2,3) ,data=att)
  pred <- Predict(mdl, 'Comp.1','Comp.2', np=30)
  
  post_pca_p[[i]] <-bplot(pred,lfun=wireframe)
  post_pca_r[[i]] <- mdl
}

## read additional scan data(18cases ; 220722) ------------------------------------
add_num <- c("291","307","310", "311", "328", "330", "332", "334", "336", "343", "349", "353", "355", "357", "359", "363", "364", "377")
home_dir <-"/Volumes/Hutom/Scan/익명화/01(신촌세브란스)_위장관외과/Surgeon_01_형우진/2.기복형성모델/211028_그래픽스팀_scan_2020_anony/"
f_path <- c()
case_num <- c(); col_n <- c(); add_n <- c()
for(i in 1:length(add_num)){
  f_path <- paste0(home_dir ,add_num[i], "/landmark_result.csv")
  r_df<- read.csv(f_path, header=FALSE)
  add_n<- r_df[grep("Landmark", r_df$V1),]$V1
  col_n <- append(col_n, add_n)
  case_num <- append(case_num, rep(add_num[i], length(add_n)))
}

a<- data.frame(case_num, col_n)
write_xlsx(a, "/Users/jinhyesu/my_project/rawData/dummy/check_colnames_220722.xlsx")

#get columns ; pre, post 
get_col <- c("PrePneumo Landmark (Scan Object)", "PostPneumo Landmark (Scan Object)")

state_col <- c("pre", "post")
add_df <- data.frame()
for(i in 1:length(add_num)){
  f_path <- paste0(home_dir ,add_num[i], "/landmark_result.csv")
  r_df<- read.csv(f_path, header=FALSE)

  for(j in 1:length(get_col)){
    start_r <- as.numeric(rownames(r_df[which(r_df$V1 %in% get_col[j]),]))    
    df_dumm <- r_df[c(start_r+1:25), ]
    df_dumm$state <- rep(state_col[j], 25)
    df_dumm$case <- rep(add_num[i], 25)
    
    add_df <-bind_rows(add_df, df_dumm)
  }
}
colnames(add_df) <-c("p_num", "x", "y", "z", "state", "case")

write_xlsx(add_df, "/Users/jinhyesu/my_project/rawData/dummy/check_df_add_18cases(scan)_220722.xlsx")

## read clinical data_ct(66cases ; 220722) -------------------------------
raw_df <- as.data.frame(read_xlsx("/Users/jinhyesu/my_project/rawData/220511_ct/Hutom Data(Final).xlsx", sheet="Hutom_Dicom(Final)", range = "B2:T382"))
colnames(raw_df)

get_col <- c("PSM_No","SEX", "AGE", "Weight(수술당시)", "Height", "BMI(수술당시)", "출산 여부(여자)")
cli_df <- raw_df[,which(colnames(raw_df) %in% get_col)] 
colnames(cli_df) <- c("psm_no", "sex", "age", "weight", "height", "bmi", "childbirth") 

cli_df$childbirth[cli_df$childbirth == "X" | cli_df$childbirth == "x"] <- "F"
cli_df$childbirth[cli_df$childbirth == "O"| cli_df$childbirth == "o"] <- "T"
cli_df$childbirth[cli_df$childbirth == "?"] <- NA

# 66cases
file_list <- list.files("/Users/jinhyesu/my_project/rawData/220511_ct/ct_data_file")
case_no <- c()
for(i in 1:length(file_list)){
  case_no <- append(case_no, strsplit(file_list[i], split = "_")[[1]][1])
}
case_no <- as.numeric(case_no)


## read distance data_ct(66cases ; 220726) -------------------------------
dist_df <- as.data.frame(read_xlsx("/Users/jinhyesu/my_project/rawData/220511_ct/dist_landmark_66cases.xlsx", range = "B1:L1651"))
dist_df$p_num <- rep(c(0:24), 66)

case <- c()
for(i in 1:length(dist_df$ann)){
    ca <- strsplit(dist_df$ann[i], "_")[[1]][1]
    case <- append(case, ca)
  }
dist_df$ann <- as.numeric(case)

# matched 2 data frame ; left_join
cli_matched <- cli_df[which(cli_df$psm_no %in% case_no),]
df_mid <- left_join(dist_df, cli_matched, by=c("ann"= "psm_no"))

str(df_mid)

view(dfSummary(cli_matched[, -1]))
view(cli_matched %>% group_by(sex)%>% dfSummary())

view(dfSummary(dist_df[,!(colnames(dist_df) %in% c("ann", "p_num"))]))

st_options(descr.stats = "common")

dist_df_by_pnum <- stby(data  = dist_df, 
                               INDICES   = dist_df$p_num, 
                               FUN       = descr)



## plotting distance in each of points -------------------------------------
hist_list <- list()
for( i in c(1:25)){
  df_sliced <- df_mid %>% filter(df_mid$p_num ==i-1)
  hist_g <- ggplot(df_sliced, aes(x=dist)) + geom_histogram(alpha=0.5, position="identity")+theme_minimal()+ggtitle(paste0(i, "th distance distribution"))+
    geom_vline(aes(xintercept = mean(dist, na.rm = TRUE)),col='pink',size=2) +  geom_text(aes(x = mean(dist, na.rm = TRUE)), y = 1, 
                                                                                          label =paste("mean = ", round(mean(df_sliced$dist, na.rm = TRUE), 2)))+
    geom_vline(aes(xintercept = median(dist, na.rm = TRUE)),col='blue',size=2) +  geom_text(aes(x = median(dist, na.rm = TRUE)), y = 2, 
                                                                                            label =paste("median = ", round(median(df_sliced$dist, na.rm = TRUE), 2)))
  hist_list[[i]]<- print(hist_g)
}

## read fat, muscle data( ; 220726) --------------------------------------
fat_df <- as.data.frame(read.csv("/Users/jinhyesu/my_project/rawData/220511_ct/measurement_pneumo.csv"))
str(fat_df)

fat_df_1 <- fat_df %>% mutate(ann = sub("^.*_0", "", fat_df$name)) %>% 
  select(contains(c("ann","area"))) 

oldnames <- colnames(fat_df_1)
newnames <- c("ann", paste0("muscle_area_", c(1:5)), paste0("fat_area_", c(1:5)))

fat_df_2 <- fat_df_1 %>% rename_at(vars(oldnames), ~ newnames) %>%  
  mutate_if(is.numeric, ~./100) %>% mutate(ann = as.numeric(ann))

df_fin <- left_join(df_mid, fat_df_2, by="ann")
write_xlsx(df_fin, "/Users/jinhyesu/my_project/rawData/220511_ct/df_fin_220816.xlsx")
setdiff(unique(fat_df_2$ann), unique(df_mid$ann)) 
#209 211 212 214 217 218 219 220 221 223
#224 225 228 229 232 233 234 235 236 237 
#238 239 240 241 246 247 249 250 254 255 
#256 257 258 259 291 307 310 321 328 332 
#334 343 346 349 355 357 359 364

na_df<- subset(df_fin, is.na(df_fin$fat_area_1))
unique(na_df$ann) #7cases 301 331 356 365 367 377 380

# histogram of fat, muscle
df_muscle <- fat_df_2 %>% drop_na() %>% select(contains("muscle"))
hist_muscle <- ggplot(gather(df_muscle), aes(x=key, y=value)) +
  geom_bar(stat='identity')+theme_minimal()+ggtitle("5 sliced muscle")+xlab("")+ylab("")

df_fat <- fat_df_2 %>% drop_na() %>% select(contains("fat"))
hist_fat <- ggplot(gather(df_fat), aes(x=key, y=value)) +
  geom_bar(stat='identity')+theme_minimal()+ggtitle("5 sliced fat")+xlab("")+ylab("")

df_fat_mucsle <- gather(fat_df_2) %>% filter(key !="ann") %>% drop_na() %>% mutate(group = str_sub(key, start= -1), 
                                                                                   subgroup = str_sub(key, start=1, end=1))

hist_fat_muscle <- ggplot(df_fat_mucsle, aes(x=group, y=value, fill = subgroup)) +
  geom_col()+theme_minimal()+ggtitle("5 sliced fat and muscle")+xlab("")+ylab("")+scale_fill_discrete(name= "subgroup", labels = c("fat", "muscle"))+
  scale_x_discrete(label = paste0(c(1:5), "th sliced")) +scale_y_continuous(labels = scales::comma)
  

plot_grid(hist_fat_muscle, NULL, hist_fat, hist_muscle,
          ncol = 2
) #cowplot

## correlation, boxplot in each points -------------------------------------------

#shapiro test for all variables ; shapiro.df
col <- c(); shapiro <- c(); qq.list <- list()
df_fin_check <- df_fin
for(i in c(2:length(df_fin_check))){
  if(is.numeric(df_fin_check[[i]])){
    df_fin_check[[i]] <- as.numeric(df_fin_check[[i]])
    col<-append(col,colnames(df_fin_check[i]))
    shapiro<-append(shapiro,shapiro.test(df_fin_check[[i]])$p.value)
    #qq<-qqnorm(df_fin_check[[i]])
    #qq.list[[i]]<- print(qq)
  } else{
    col <- append(col, colnames(df_fin_check[i]))
    shapiro <- append(shapiro, "could_not")
    #qq.list[[i]] <- "could_not"
  }
}
shapiro.df<- data.frame(col, shapiro) # lower than 0.05 / abnormal distribution

#correlation; continuous variables
df_fin_cor <- df_fin
#df_fin_cor <- df_fin %>% filter(ann != 294)
df_fin_cor<- df_fin_cor %>% mutate(sex = if_else(sex == "M", 1, 2), # man = 1
                                   childbirth = if_else(childbirth == "F", 1, 2)) # childbirth, F =  1


conti_var <- c("age", "weight", "height", "bmi", paste0("muscle_area_", c(1:5)), paste0("fat_area_", c(1:5)))
cate_var <- c("sex", "childbirth")
bgr_color <- c("#FFFFFF","#fff8f7")# p-value significant

p_val <- c(); level_list<-c(); xvar<-c(); point_list<-c(); tau_val <- c(); gender_var<-c()
for(j in 1:25){
  df_fin_cor_sub <- df_fin_cor %>% filter(p_num == j-1)
  df_fin_cor_sub_m <- df_fin_cor_sub[which(df_fin_cor_sub$sex == 1),]
  df_fin_cor_sub_f <- df_fin_cor_sub[which(df_fin_cor_sub$sex == 2),]
  
  if(j >= 1 & j< 6){
    conti_var <- c("age", "weight", "height", "bmi", "muscle_area_1", "fat_area_1")
  }else if(j >= 6 & j< 11){
    conti_var <- c("age", "weight", "height", "bmi", "muscle_area_2", "fat_area_2")
  }else if(j >= 11 & j< 16){
    conti_var <- c("age", "weight", "height", "bmi", "muscle_area_3", "fat_area_3")
  }else if(j >= 16 & j< 21){
    conti_var <- c("age", "weight", "height", "bmi", "muscle_area_4", "fat_area_4")
  }else if(j >= 21 & j< 26){
    conti_var <- c("age", "weight", "height", "bmi", "muscle_area_5", "fat_area_5")
  }
  
  for(i in 1:length(conti_var)){
    #correlation by gender
    corr_m <- cor.test(df_fin_cor_sub_m$dist, df_fin_cor_sub_m[, colnames(df_fin_cor_sub_m) %in% conti_var[i]], method="kendall") #male
    corr_f <- cor.test(df_fin_cor_sub_f$dist, df_fin_cor_sub_f[, colnames(df_fin_cor_sub_f) %in% conti_var[i]], method="kendall") #female 
    
    cor_plot <- ggplot(df_fin_cor_sub, aes(df_fin_cor_sub$dist, df_fin_cor_sub[, colnames(df_fin_cor_sub)%in% conti_var[i]], 
                                           color=factor(sex, levels=c(1,2), labels=c("M", "W"))))+
      geom_point()+ geom_smooth(method = "lm", se = FALSE, aes(group=sex))+
      annotate("text",x=max(df_fin_cor_sub$dist)*0.6, y= max(df_fin_cor_sub[, colnames(df_fin_cor_sub) %in% conti_var[i]],  na.rm = TRUE)*0.6, 
               label = paste0("tau(M) =", round(corr_m[['estimate']],4), "  ",
                              "tau(W) =", round(corr_f[['estimate']],4), "\n", 
                              "p(M) =", format(round(corr_m[['p.value']], 4), nsmall=4, scientific = FALSE), "  ",
                              "p(W) =", format(round(corr_f[['p.value']], 4), nsmall=4, scientific = FALSE))) + 
      theme_minimal()+ggtitle(paste(j, "th point :", "distance, ", conti_var[i])) + xlab("distansce")+ ylab(conti_var[i])+
      theme(panel.background = element_rect(fill = ifelse(corr_m[['p.value']] <=0.05 | corr_f[['p.value']] <=0.05, bgr_color[2],bgr_color[1])))+ 
      scale_color_manual("gender", values= c("#a4dadc", "#ebbab6"))
    print(cor_plot)
    
    #generate df
    point_list <- append(point_list, rep(j,2))
    xvar <- append(xvar, rep(conti_var[i],2))
    p_val <- append(p_val, c(round(corr_m[['p.value']], 4), round(corr_f[['p.value']], 4)))
    tau_val <- append(tau_val, c(round(corr_m[['estimate']], 4), round(corr_m[['estimate']], 4)))
    gender_var <- append(gender_var, c("M", "F"))
  }
  for(k in 1:length(cate_var)){
    corr <- cor.test(df_fin_cor_sub$dist, df_fin_cor_sub[, colnames(df_fin_cor_sub)%in% cate_var[k]], method="kendall") #성별을 구분하여 상관분석하는 것이 의미없음.
    cor_plot <- ggplot(df_fin_cor_sub, aes(df_fin_cor_sub$dist, as.factor(df_fin_cor_sub[, colnames(df_fin_cor_sub)%in% cate_var[k]]), 
                                           color=factor(sex, levels=c(1,2), labels=c("M", "W"))))+
      geom_boxplot(outlier.colour="red", outlier.shape=18, outlier.size=4)+
      geom_jitter(shape=20, position=position_jitter(0.1))+
      geom_smooth(method = "lm", se = FALSE,  aes(group=1), colour="#DDA28F")+
      annotate("text",x=max(df_fin_cor_sub$dist)*0.6, y= max(df_fin_cor_sub[, colnames(df_fin_cor_sub) %in% cate_var[k]],  na.rm = TRUE)*0.6, 
               label = paste0("tau =", round(corr[['estimate']],4), "\n", "p =", format(round(corr[['p.value']], 4), nsmall=4, scientific = FALSE)))+ 
      theme_minimal()+ggtitle(paste(j, "th point :", "distance, ", cate_var[k])) + xlab("distansce")+ ylab(cate_var[k])+
      theme(panel.background = element_rect(fill = case_when(corr[['p.value']] <=0.05 ~ bgr_color[2],
                                                             corr[['p.value']] >0.05 ~ bgr_color[1])))+
      scale_color_manual("gender", values= c("#a4dadc", "#ebbab6"))

    print(cor_plot)
    
    #generate df
    point_list <- append(point_list, rep(j,2))
    xvar <- append(xvar, rep(conti_var[i],2))
    p_val <- append(p_val, rep(round(corr[['p.value']], 4), 2))
    tau_val <- append(tau_val, rep(round(corr[['estimate']], 4),2))
    gender_var <- append(gender_var, c("M", "F"))
    }
}

df_test_fin<- data.frame(point_list, xvar, gender_var, p_val, tau_val) %>% 
                      mutate(level_list = case_when(p_val >0.1 ~ 'not significant',
                                                                       p_val <= 0.1 & p_val > 0.05 ~ '* significant', 
                                                                       p_val <=0.05 ~ '** significant'))
df_test_fin$level_list <- as.factor(df_test_fin$level_list)
                    
## heatmap of the results ------------------------------------------------
# 25 points / clinical variables 
#corr_label <- c('weak : r < 0.3','moderate : 0.3 =< r < 0.7','strong : r >=0.7')
pval_label <- c('not significant : p > 0.1',
                '* significant : 0.05 < p <= 0.1',
                '** significant : p <= 0.05')

#p_val_color_list <-  c("#FAF6F2","#FFE4C4", "#FFA07A") #red base
p_val_color_list <- c("#FAF6F2","#BAD7FF", "#5E9DF2") # blue base

heat_g <- ggplot(data = df_test_fin,mapping = aes(y=factor(xvar,levels =c( "age", "weight", "height", "bmi", "sex", "childbirth", "fat_area_1", "fat_area_2",
                                                                            "fat_area_3", "fat_area_4","fat_area_5","muscle_area_1","muscle_area_2",
                                                                             "muscle_area_3", "muscle_area_4","muscle_area_5")),
                                                           x = factor(point_list, levels = c(1:25)) , fill =level_list,label = paste0(p_val)))+
  geom_tile(color = 'black')+ 
  #geom_text() + 
  scale_fill_manual(values = p_val_color_list,label = pval_label)+
  ggtitle('Kendall\'s correlation (result: p-value)')+
  theme_minimal()+
  ylab('')+
  xlab('')+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text( vjust = 0.5, hjust=1), 
        legend.position = "none")# axis.text.x=element_text(angle=90, hjust=1)

# a) male
heat_g_m <- ggplot(data = df_test_fin[which(df_test_fin$gender_var =="M"),],mapping = aes(y=factor(xvar,levels =c( "age", "weight", "height", "bmi", "sex", "childbirth", "fat_area_1", "fat_area_2",
                                                                           "fat_area_3", "fat_area_4","fat_area_5","muscle_area_1","muscle_area_2",
                                                                           "muscle_area_3", "muscle_area_4","muscle_area_5")),
                                                  x = factor(point_list, levels = c(1:25)) , fill =level_list,label = paste0(p_val)))+
  geom_tile(color = 'black')+ 
  #geom_text() + 
  scale_fill_manual(values = p_val_color_list,label = pval_label)+
  ggtitle('Kendall\'s correlation (result: p-value)')+
  theme_minimal()+
  ylab('')+
  xlab('')+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text( vjust = 0.5, hjust=1), 
        legend.position = "none")# axis.text.x=element_text(angle=90, hjust=1)

# b) female
heat_g_f <- ggplot(data = df_test_fin[which(df_test_fin$gender_var =="F"),],mapping = aes(y=factor(xvar,levels =c( "age", "weight", "height", "bmi", "sex", "childbirth", "fat_area_1", "fat_area_2",
                                                                                                                   "fat_area_3", "fat_area_4","fat_area_5","muscle_area_1","muscle_area_2",
                                                                                                                   "muscle_area_3", "muscle_area_4","muscle_area_5")),
                                                                                          x = factor(point_list, levels = c(1:25)) , fill =level_list,label = paste0(p_val)))+
  geom_tile(color = 'black')+ 
  #geom_text() + 
  scale_fill_manual(values = p_val_color_list,label = pval_label)+
  ggtitle('Kendall\'s correlation (result: p-value)')+
  theme_minimal()+
  ylab('')+
  xlab('')+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text( vjust = 0.5, hjust=1), 
        legend.position = "none")# axis.text.x=element_text(angle=90, hjust=1)


# tau ; positive, negative 
# among significant results(*, **)

df_test_fin2 <- df_test_fin %>% mutate(tau_level = case_when(tau_val > 0 & level_list =='* significant'~ "* positive",
                                                             tau_val > 0 & level_list =='** significant'~ "** positive",
                                                             tau_val < 0 & level_list =='* significant'~ "* negative",
                                                             tau_val < 0 & level_list =='** significant'~ "** negative",
                                                             level_list =='not significant' ~ 'non'))


tau_label <- c("non","* positive", "** positive","* negative","** negative")
tau_color_list <- c("#fcfcfc","#bfe0cb","#74b88d", "#FFE4C4","#FFA07A") #green, red
heat_g_esti <- ggplot(data = df_test_fin2,mapping = aes(y=factor(xvar,levels =c( "age", "weight", "height", "bmi", "sex", "childbirth", "fat_area_1", "fat_area_2",
                                                                           "fat_area_3", "fat_area_4","fat_area_5","muscle_area_1","muscle_area_2",
                                                                           "muscle_area_3", "muscle_area_4","muscle_area_5")),
                                                  x = factor(point_list, levels = c(1:25)) , fill =tau_level, label=tau_level))+
  geom_tile(color = 'black')+ 
  #geom_text() + 
  scale_fill_manual(values = tau_color_list,breaks= tau_label, label = tau_label)+
  ggtitle('Kendall\'s correlation (result: tau)')+
  theme_minimal()+
  ylab('')+
  xlab('')+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text( vjust = 0.5, hjust=1))

# a) male
heat_g_esti_m <- ggplot(data = df_test_fin2[which(df_test_fin2$gender =="M"),],mapping = aes(y=factor(xvar,levels =c( "age", "weight", "height", "bmi", "sex", "childbirth", "fat_area_1", "fat_area_2",
                                                                                 "fat_area_3", "fat_area_4","fat_area_5","muscle_area_1","muscle_area_2",
                                                                                 "muscle_area_3", "muscle_area_4","muscle_area_5")),
                                                        x = factor(point_list, levels = c(1:25)) , fill =tau_level, label=tau_level))+
  geom_tile(color = 'black')+ 
  #geom_text() + 
  scale_fill_manual(values = tau_color_list,breaks= tau_label, label = tau_label)+
  ggtitle('Kendall\'s correlation (result: tau)')+
  theme_minimal()+
  ylab('')+
  xlab('')+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text( vjust = 0.5, hjust=1))

# a) female
heat_g_esti_f <- ggplot(data = df_test_fin2[which(df_test_fin2$gender =="F"),],mapping = aes(y=factor(xvar,levels =c( "age", "weight", "height", "bmi", "sex", "childbirth", "fat_area_1", "fat_area_2",
                                                                                                                      "fat_area_3", "fat_area_4","fat_area_5","muscle_area_1","muscle_area_2",
                                                                                                                      "muscle_area_3", "muscle_area_4","muscle_area_5")),
                                                                                             x = factor(point_list, levels = c(1:25)) , fill =tau_level, label=tau_level))+
  geom_tile(color = 'black')+ 
  #geom_text() + 
  scale_fill_manual(values = tau_color_list,breaks= tau_label, label = tau_label)+
  ggtitle('Kendall\'s correlation (result: tau)')+
  theme_minimal()+
  ylab('')+
  xlab('')+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text( vjust = 0.5, hjust=1))



##check correlation between x variables-----------------------------------------------------
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ..., method="kendall")
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
df_cor_m <- df_fin_cor_sub_m %>% select(c("age", "weight", "height", "bmi", paste0("muscle_area_", 1:5), paste0("fat_area_", 1:5)))
p.mat1 <- cor.mtest(df_cor_m)
M <- cor(df_cor_m, use = "na.or.complete", method = "kendall")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
cor_m_xvar<- corrplot::corrplot(M, method="color", col=col(200), title ="male", 
         type="upper", #order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat1, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,number.cex = 0.6
)
df_cor_f <- df_fin_cor_sub_f %>% select(c("age", "weight", "height", "bmi", paste0("muscle_area_", 1:5), paste0("fat_area_", 1:5)))
p.mat2 <- cor.mtest(df_cor_f)
M2 <- cor(df_cor_f, use = "na.or.complete", method = "kendall")

cor_f_xvar<- corrplot::corrplot(M2, method="color", col=col(200), title ="female",
                 type="upper", #order="hclust", 
                 addCoef.col = "black", # Add coefficient of correlation
                 tl.col="black", tl.srt=45, #Text label color and rotation
                 # Combine with significance
                 p.mat = p.mat2, sig.level = 0.05, insig = "blank", 
                 # hide correlation coefficient on the principal diagonal
                 diag=FALSE 
)

# all; gender
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ..., method="kendall")
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat1 <- cor.mtest(df_fin_cor[,c(13:28)])
M <- cor(df_fin_cor[,c(13:28)], use = "na.or.complete", method = "kendall")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
cat("male")
corrplot::corrplot(M, method="color", col=col(200),  
                   type="upper", order="hclust", 
                   addCoef.col = "black", # Add coefficient of correlation
                   tl.col="black", tl.srt=45, #Text label color and rotation
                   # Combine with significance
                   p.mat = p.mat1, sig.level = 0.05, insig = "blank", 
                   # hide correlation coefficient on the principal diagonal
                   diag=FALSE,number.cex = 0.6
)

p.mat2 <- cor.mtest(df_fin_cor[,c(13:18)])
M2 <- cor(df_fin_cor[,c(13:18)], use = "na.or.complete", method = "kendall")
cat("female")
corrplot::corrplot(M2, method="color", col=col(200),  
                   type="upper", order="hclust", 
                   addCoef.col = "black", # Add coefficient of correlation
                   tl.col="black", tl.srt=45, #Text label color and rotation
                   # Combine with significance
                   p.mat = p.mat2, sig.level = 0.05, insig = "blank", 
                   # hide correlation coefficient on the principal diagonal
                   diag=FALSE 
)





# generate dummy variables

df_fin_reg <- df_fin %>% select("dist", "p_num","sex","age","weight","height","bmi","childbirth",
                                "muscle_area_1", "muscle_area_2","muscle_area_3", "muscle_area_4", "muscle_area_5", 
                                "fat_area_1", "fat_area_2", "fat_area_3", "fat_area_4", "fat_area_5")%>% 
                        mutate(sex = factor(sex, levels = c("M", "F")),
                              childbirth = factor(childbirth, levels = c("F", "T", "NA")))

reg_list <- list()
var_list <-c()
for(j in 1:25){
  base <- c("dist","sex","age", "weight", "height", "bmi")
  if(j >= 1 & j< 6){
    #extract childbirth 
    var_list <- c(base, "muscle_area_1", "fat_area_1") 
  }else if(j >= 6 & j< 11){
    var_list <- c(base, "muscle_area_2", "fat_area_2") 
  }else if(j >= 11 & j< 16){
    var_list <- c(base, "muscle_area_3", "fat_area_3") 
  }else if(j >= 16 & j< 21){
    var_list <- c(base, "muscle_area_4", "fat_area_4") 
  }else if(j >= 21 & j< 26){
    var_list <- c(base, "muscle_area_5", "fat_area_5") 
  }
  reg_list[[j]] <- list()
  df_fin_reg_sub <- df_fin_reg %>% filter(p_num == j-1) %>% select(var_list)
  
  lmfit <- lm(dist ~ . , data = df_fin_reg_sub)
  reg_list[[j]][['lm']] <- lmfit
  vif_df<- car::vif(lmfit)
  print(j); print(lmfit); print(vif_df)
  reg_list[[j]][['vif']] <- as.data.frame(vif_df)
}

df_fin_reg <- df_fin %>% select("dist", "p_num","sex","age","weight","height","bmi","childbirth",
                                "muscle_area_1", "muscle_area_2","muscle_area_3", "muscle_area_4", "muscle_area_5", 
                                "fat_area_1", "fat_area_2", "fat_area_3", "fat_area_4", "fat_area_5") %>% 
                      mutate(sex = if_else(sex == 1, 0, 1), # man = 0
                             childbirth = if_else(childbirth == 1, 0, 1)) # childbirth, F =  0
df_fin_reg <- df_fin %>% select("dist", "p_num","sex","age","weight","height","bmi","childbirth",
                                "muscle_area_1", "muscle_area_2","muscle_area_3", "muscle_area_4", "muscle_area_5", 
                                "fat_area_1", "fat_area_2", "fat_area_3", "fat_area_4", "fat_area_5")%>% 
  mutate(sex = factor(sex, levels = c("M", "F")),
         childbirth = factor(childbirth, levels = c("F", "T", "NA")))


sp_reg <- function(xvar){
  df_fin_reg <- df_fin %>% select("dist", "p_num","sex","age","weight","height","bmi","childbirth",
                                  "muscle_area_1", "muscle_area_2","muscle_area_3", "muscle_area_4", "muscle_area_5", 
                                  "fat_area_1", "fat_area_2", "fat_area_3", "fat_area_4", "fat_area_5")%>% 
    mutate(sex = factor(sex, levels = c("M", "F")),
           childbirth = factor(childbirth, levels = c("F", "T", "NA")))
  simple_list <- list()
  var_list <-c()
  for(j in 1:25){
    base <- c("dist")
    var_list <- append(base, xvar)
    simple_list[[j]] <- list()
    df_fin_reg_sub <- df_fin_reg %>% filter(p_num == j-1) %>% select(var_list)
    
    lmfit <- lm(dist ~ . , data = df_fin_reg_sub)
    simple_list[[j]][['lm']] <- lmfit
    #simple_list[[j]][['vif']] <- as.data.frame(vif(lmfit))
  }
  return(simple_list)
}

sp_sex<- sp_reg("sex")
sp_age <- sp_reg("age")
sp_weight <- sp_reg("weight")
sp_height <- sp_reg("height")
sp_bmi <- sp_reg("bmi")
sp_childbirth <- sp_reg("childbirth")


p_model <- function(model_list){
  for( j in c(1, 6, 11, 16, 21)){
    cat("\n\n#####", j, "-", j+4, "th_point \n\n")
    cat(tab_model(model_list[[j]][['lm']], model_list[[j+1]][['lm']],model_list[[j+2]][['lm']],model_list[[j+3]][['lm']],model_list[[j+4]][['lm']],
                  dv.labels = c(j, j+1,j+2,j+3, j+4), auto.label = FALSE,  show.ci = FALSE, 
                  CSS = list(css.thead="border-top:double black; font-weight:bold; font-size:1.2em;",
                             css.firsttablecol="font-size:0.9em;",
                             css.modelcolumn2 = "background-color : #E7E6D2;",
                             css.modelcolumn4 = "background-color : #E7E6D2;",
                             css.summary="font-weight:bold; font-sime:0.7em;"))$knitr)
    if(length(model_list[[j]]) ==2){
      print(knitr::kables(
        list(
          model_list[[j]][['vif']]  %>% 
            knitr::kable(booktabs=T, digits = 2, format = "html", 
                         col.names = j, full_width = FALSE, position = "right") %>% 
            kable_styling(full_width = F) %>% 
            row_spec(which(model_list[[j]][['vif']] >10), 
                     bold =T, color="white", background = "red"),
          model_list[[j+1]][['vif']]  %>% 
            knitr::kable(booktabs=T, digits = 2, format = "html", 
                         col.names = j+1, full_width = FALSE, position = "right",row.names = NA) %>% 
            kable_styling(full_width = F) %>% 
            column_spec(1:2, background ="#E7E6D2") %>% 
            row_spec(which(model_list[[j]][['vif']] >10), 
                     bold =T, color="white", background = "red"),
          model_list[[j+2]][['vif']]  %>% 
            knitr::kable(booktabs=T, digits = 2, format = "html", 
                         col.names = j+2, full_width = FALSE, position = "right",row.names = NA) %>% 
            kable_styling(full_width = F) %>% 
            row_spec(which(model_list[[j]][['vif']] >10), 
                     bold =T, color="white", background = "red"),
          model_list[[j+3]][['vif']]  %>% 
            knitr::kable(booktabs=T, digits = 2, format = "html",
                         col.names = j+3, full_width = FALSE, position = "right",row.names = NA) %>% 
            kable_styling(full_width = F) %>% 
            column_spec(1:2, background ="#E7E6D2") %>% 
            row_spec(which(model_list[[j]][['vif']] >10), 
                     bold =T, color="white", background = "red"),
          model_list[[j+4]][['vif']]  %>% 
            knitr::kable(booktabs=T, digits = 2, format = "html",  
                         col.names = j+4, full_width = FALSE, position = "right",row.names = NA) %>% 
            kable_styling(full_width = F) %>% 
            row_spec(which(model_list[[j]][['vif']] >10), 
                     bold =T, color="white", background = "red")
        )))
    }
    cat('\n\n<!-- -->\n\n')
    
  }
}


## fat, muscle distribution(220811) ----------------------------------------------
prior_each_of_hist<- function(){
  hist_mf <- function(df){
    var_list <- c(c(paste0("muscle_area_", 1:5)), c(paste0("fat_area_", 1:5)))
    
    mf_hist_list <- list()
    for(i in c(1:length(var_list))){
      df_sliced <- df %>% select(var_list[i])
      hist_g <- ggplot(df_sliced, aes_string(x=var_list[i])) + geom_histogram(bins = 15)+theme_minimal()+ggtitle(var_list[i])+
        scale_x_continuous(breaks = seq(0, max(df_sliced), round(max(df_sliced)/8, 0)))
      
      mf_hist_list[[i]]<- print(hist_g)
    }
    return(mf_hist_list)
  }
  
  #sex : all / 59
  df_mf <- df_fin %>% select("ann", "sex", "age", "weight", "height", "bmi", "childbirth", c(paste0("muscle_area_", 1:5)), c(paste0("fat_area_", 1:5))) %>% 
    distinct() %>% drop_na("muscle_area_1")
  head(df_mf); str(df_mf)
  
  hist_mf_all <- hist_mf(df_mf)
  
  #table(df_mf$sex)
  # sex : M / 33
  df_mf_m <- df_mf %>% filter(sex == "M")
  hist_mf_m <- hist_mf(df_mf_m)
  # sex : F / 26
  df_mf_f <- df_mf %>% filter(sex == "F")
  hist_mf_f <- hist_mf(df_mf_f)
}

# in one plot
for(i in c(1:length(var_list))){
  hist_g <- ggplot(df_mf, aes_string(x= var_list[i], fill="sex", color="sex")) + geom_histogram(bins = 15, position="identity", alpha=0.5)+
    ggtitle(var_list[i])+scale_x_continuous(breaks = seq(0, max(df_sliced), round(max(df_sliced)/8, 0)))
  
  mf_hist_list[[i]]<- print(hist_g)
}

## update distance grouped(220811) -------------------------------------------
# sex, age(60), bmi(18.5, 23) / 이상, 미만
df_fin_g <- df_fin %>% mutate(age_2group = factor(ifelse(age >= 60, 2, 1), levels = c(2, 1), labels = c("age>=60", "60>age")), 
                              bmi_3group = factor(case_when(bmi >= 23 ~ 3,
                                                            bmi < 23 & bmi >=18.5 ~2, 
                                                            bmi <18.5 ~ 1), levels = c(3, 2, 1), labels = c("bmi>=23", "18.5<=bmi<23", "bmi<18.5")))

hist_by_group <- function(xvar){
    var_list <- c("ann", "dist", "p_num", xvar)
    df_dist <- df_fin_g %>% select(var_list)
    
    hist_d_list <- list()
    for(j in c(1:25)){
      #color by xvar
      if(xvar == "sex"){
        color_values <- c("#a4dadc", "#ebbab6")
      }else if(xvar =="age_2group"){
        color_values <- c("#00AFBB", "#E7B800")
      }else if(xvar=="bmi_3group"){
        color_values <- c("#755b3b", "#eab676","#f2d3ad")
      }
      df_dist_i <- df_dist %>% filter(p_num == j-1)
      mu <- df_dist_i %>% group_by(df_dist_i[,colnames(df_dist_i) %in% xvar]) %>% dplyr::summarise(mean = mean(dist)) %>% as.data.frame()
      
      hist_g <- ggplot(df_dist_i, aes(x= dist, fill=df_dist_i[,colnames(df_dist_i) %in% xvar], 
                                      color=df_dist_i[,colnames(df_dist_i) %in% xvar])) + 
        geom_histogram(bins = 15, position="identity", alpha=0.5)+
        ggtitle(paste0(j-1, "th point"))+
        geom_vline(data=mu, aes(xintercept=mean, color=unique(df_dist_i[,colnames(df_dist_i) %in% xvar])),linetype="dashed")+ 
        scale_x_continuous(breaks = seq(0, max(df_dist_i$dist),round(max(df_dist_i$dist)/8, 0)))+
        theme(legend.title = element_blank())+ 
        scale_color_manual(values = color_values) +
        scale_fill_manual(values = color_values)
      
      hist_d_list[[j]]<- print(hist_g)
    }
    return(hist_d_list)
}
hist_gender <- hist_by_group("sex")
hist_age_2group <- hist_by_group("age_2group")
hist_bmi_3group <- hist_by_group("bmi_3group")

## fat, muscle percentage by ann(220811) -----------------------------------------
each_mf<- function(df){
  mf_dist <- list()
  for(i in c(1:5)){
    mf_dist[[i]] <- list()
    index <- c("ann", i)
    a<- df[,grep(paste(index, collapse="|"), colnames(df))]
    a<- a %>% mutate(sum = a[,colnames(a) %in% paste0("muscle_area_", i)] + a[,colnames(a) %in% paste0("fat_area_", i)])
    b<- gather(a, paste0("muscle_area_", i), paste0("fat_area_", i), key = 'state', value = 'area') %>% mutate(ann = as.factor(ann), 
                                                                                                               prop = round(area/sum*100, 0),
                                                                                                               state = as.factor(state))
    # Stacked
    mf_dist[[i]][['stack']] <- print(ggplot(b, aes(fill=state, y=area, x=ann)) + 
                                       geom_bar(position="stack", stat="identity")+  geom_text(aes(label=paste0(prop, "%")), position = position_stack(.5), color="white", size=3.5) +
                                       ylab("area(mm^2)") + xlab("case_number(#)") + ggtitle(paste0(i, "slide")))
    
    # filled_100 percent
    mf_dist[[i]][['fill']] <- print(ggplot(aes(x = ann, y = prop,fill = factor(state)), data = b) + geom_bar(stat="identity")+  
                                      geom_text(aes(label=paste0(prop, "%")), position = position_stack(.5), color="white", size=3.5)+scale_y_continuous(breaks = seq(0, 100, 5)) +
                                      ylab("percert(%)") + xlab("case_number(#)")+ ggtitle(paste0(i, "slide")))
  }
  return(mf_dist)
}

df_ann <- df_fin %>% select("ann", "sex", c(paste0("muscle_area_", 1:5)), c(paste0("fat_area_", 1:5))) %>% distinct() %>% drop_na()
df_ann_m <- df_ann %>%  filter(sex == "M") #33
m_mf <- each_mf(df_ann_m)

df_ann_f <- df_ann %>%  filter(sex == "F") #26
f_mf <- each_mf(df_ann_f)

## [hist] clinical variable's count(220812)  -------------------------------------
#age, weight, height, bmi, childbirth
cli_var<- c("age", "weight", "height", "bmi")
mu <- df_fin %>% group_by(sex) %>% dplyr::summarise(age_mean = mean(age), 
                                                    weight_mean = mean(weight),
                                                    height_mean = mean(height),
                                                    bmi_mean = mean(bmi)) %>% as.data.frame()
cli_hist <- list()
for(i in 1:length(cli_var)){

  hist_g <- ggplot(df_fin, aes(x= df_fin[,colnames(df_fin) %in% cli_var[i]], fill=sex, color=sex)) + 
  geom_histogram(bins = 15, position="identity", alpha=0.5)+ 
  geom_vline(data=mu, aes(xintercept=mu[,i+1], color=sex),linetype="dashed")+ 
  theme(legend.title = element_blank())+xlab(cli_var[i])
  cli_hist[[i]] <- print(hist_g)
}

## data normalization; scale(220816) --------------------------------------------
df_fin<- as.data.frame(read_xlsx("/Users/jinhyesu/my_project/rawData/220511_ct/df_fin_220816.xlsx"))

df_fin_scale <- df_fin %>% mutate(age_2group = factor(ifelse(age >= 60, 2, 1), levels = c(2, 1), labels = c("age>=60", "60>age")), 
                                             bmi_3group = factor(case_when(bmi >= 23 ~ 3,
                                                                           bmi < 23 & bmi >=18.5 ~2, 
                                                                           bmi <18.5 ~ 1), levels = c(3, 2, 1), labels = c("bmi>=23", "18.5<=bmi<23", "bmi<18.5"))) %>%
  mutate_at(c("dist","age", "weight", "height", "bmi", 
                                       paste0("muscle_area_", 1:5), paste0("fat_area_", 1:5)),
                                     ~(scale(.) %>% as.vector))

write_xlsx(df_fin_scale, "/Users/jinhyesu/my_project/rawData/220511_ct/df_fin_scale_220816.xlsx")


## data normalization and regression(220817) ---------------------------------
df_fin<- as.data.frame(read_xlsx("/Users/jinhyesu/my_project/rawData/220511_ct/df_fin_220816.xlsx"))
df_mid = preProcess(x = df_fin[which(colnames(df_fin) %in% c("dist","age", "weight", "height", "bmi", 
                                                                paste0("muscle_area_", 1:5), paste0("fat_area_", 1:5)))], method = "range")
df_minmax <- predict(df_mid, df_fin)
df_reg <- dummy_cols(df_minmax, select_columns = c("sex", "childbirth"), remove_first_dummy = TRUE, ignore_na=TRUE)
                             
df_fin_reg <- df_reg %>% select("dist", "p_num","age","weight","height","bmi", 
                                paste0("muscle_area_", 1:5), paste0("fat_area_", 1:5), 
                                "sex_M", "childbirth_T") #dummy variables


eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  
  
  # Model performance metrics
  table <- data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
  return(table)
}

ridge_list <- list()
lasso_list <- list()
var_list <-c()
for(j in 1:25){
  base <- c("dist","sex_M","age", "weight", "height", "bmi", "childbirth_T")
  if(j >= 1 & j< 6){
    #extract childbirth 
    var_list <- c(base, "muscle_area_1", "fat_area_1") 
  }else if(j >= 6 & j< 11){
    var_list <- c(base, "muscle_area_2", "fat_area_2") 
  }else if(j >= 11 & j< 16){
    var_list <- c(base, "muscle_area_3", "fat_area_3") 
  }else if(j >= 16 & j< 21){
    var_list <- c(base, "muscle_area_4", "fat_area_4") 
  }else if(j >= 21 & j< 26){
    var_list <- c(base, "muscle_area_5", "fat_area_5") 
  }
  lasso_list[[j]] <- list(); ridge_list[[j]] <- list()
  df_fin_reg_sub <- df_fin_reg %>% filter(p_num == j-1) %>% select(var_list)
  
  set.seed(100)
  index = sample(1:nrow(df_fin_reg_sub), 0.7*nrow(df_fin_reg_sub))
  
  train  = df_fin_reg_sub[index,]
  test  = df_fin_reg_sub[-index,]
  
  #a) ridge
  fit_1 <- glmnet(train[,-1], train[,1], family="gaussian", alpha=0) # alpha=0, ridge
  cv_1 <- cv.glmnet(as.matrix(train[,-1]), train[,1], family="gaussian", alpha=0)
    #plot(fit_2, xvar="norm")
  
  # plot with ggplot ; plot(fit_2, xvar="lambda") + title(paste0(j, "th coefficient"))
  beta = coef(fit_1)
  tmp <- as.data.frame(as.matrix(beta))
  tmp$coef <- row.names(tmp)
  tmp <- reshape::melt(tmp, id = "coef")
  tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
  tmp$lambda <- log(fit_1$lambda[tmp$variable+1]) # 람다 추출 및 log10변환

  ridge_list[[j]][['g_coef']] <- print(ggplot(tmp[tmp$coef != "(Intercept)",], aes(lambda, value, color = coef, linetype = coef)) + 
                                       geom_line(size=1) +  
                                       xlab("Lambda (log scale)") + 
                                       ylab("Coefficients") +
                                       guides(color = guide_legend(title = ""), 
                                                      linetype = guide_legend(title = "")) +
                                       theme_bw() + 
                                       theme(legend.key.width = unit(3,"lines")))
                                  
  ridge_list[[j]][['g_mse']] <- plot(cv_1)
  ridge_list[[j]][['lambda']] <- cv_1$lambda.min
  ridge_list[[j]][['best_coef']] <- as.matrix(coef(fit_1, s= cv_1$lambda.min))

  #use fitted best model to make predictions
  pred_y <- predict(cv_1, s ="lambda.min", newx = as.matrix(test[,-1]))
  pred_y_df <- as.data.frame(pred_y)
  
  ridge_list[[j]][['g_pred']] <- print(ggplot() + geom_point(mapping=aes(x= rownames(test), y= test[,1]),data=test)+
                                         geom_point(mapping=aes(x=rownames(train), y=train[,1]), data=train, colour ="pink")+
                                         geom_line(mapping=aes(x=rownames(pred_y_df), y=pred_y, group=1), data=pred_y_df, col="red")+
       #geom_text(mapping=aes(x=rownames(pred_y_df), y=pred_y, group=1), data=pred_y_df, label = round(pred_y, 2))+
       theme_classic()+xlab("test_y_values"))
  
  ridge_list[[j]][['eval']] <- eval_results(test[,1],as.vector(pred_y), test)
  
  #b) lasso
  fit_2 <- glmnet(train[,-1], train[,1], family="gaussian", alpha=1) # lasso, alpha=1
  cv_2 <- cv.glmnet(as.matrix(train[,-1]), train[,1], family="gaussian", alpha=1)
  #plot(fit_2, xvar="norm")
  
  # plot with ggplot ; plot(fit_2, xvar="lambda") + title(paste0(j, "th coefficient"))
  beta = coef(fit_2)
  tmp <- as.data.frame(as.matrix(beta))
  tmp$coef <- row.names(tmp)
  tmp <- reshape::melt(tmp, id = "coef")
  tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
  tmp$lambda <- log(fit_2$lambda[tmp$variable+1]) # 람다 추출 및 log10변환
  
  lasso_list[[j]][['g_coef']] <- print(ggplot(tmp[tmp$coef != "(Intercept)",], aes(lambda, value, color = coef, linetype = coef)) + 
                                         geom_line(size=1) +  
                                         xlab("Lambda (log scale)") + 
                                         ylab("Coefficients") +
                                         guides(color = guide_legend(title = ""), 
                                                linetype = guide_legend(title = "")) +
                                         theme_bw() + 
                                         theme(legend.key.width = unit(3,"lines")))
  
  lasso_list[[j]][['g_mse']] <- print(plot(cv_2))
  lasso_list[[j]][['lambda']] <- cv_2$lambda.min
  lasso_list[[j]][['best_coef']] <- as.matrix(coef(fit_2, s= cv_2$lambda.min))
  
  #use fitted best model to make predictions
  pred_y <- predict(cv_2, s ="lambda.min", newx = as.matrix(test[,-1]))
  pred_y_df <- as.data.frame(pred_y)
  
  lasso_list[[j]][['g_pred']] <- print(ggplot() + geom_point(mapping=aes(x= rownames(test), y= test[,1]),data=test)+
                                         geom_point(mapping=aes(x=rownames(train), y=train[,1]), data=train, colour ="pink")+
                                         geom_line(mapping=aes(x=rownames(pred_y_df), y=pred_y, group=1), data=pred_y_df, col="red")+
                                         #geom_text(mapping=aes(x=rownames(pred_y_df), y=pred_y, group=1), data=pred_y_df, label = round(pred_y, 2))+
                                         theme_classic()+xlab("test_y_values"))
  
  lasso_list[[j]][['eval']] <- eval_results(test[,1],as.vector(pred_y), test)
}

#사용하는 인덱스 번호 전달(정은님)
head(df_fin)
index_list<- as.data.frame(unique(df_fin$ann))
write.csv(index_list, "/Users/jinhyesu/my_project/rawData/dummy/ct_data_number(66cases)_220818.csv")

## read fat, muscle_updated(220822) -------------------------------------------
fat_updated<- as.data.frame(read.csv("/Users/jinhyesu/my_project/rawData/220511_ct/measurement_pneumo_220822.csv",header=TRUE)) # 106 cases, additional(width, height)
str(fat_updated)

body_h <- fat_updated %>% mutate(ann = as.numeric(sub("^.*_0", "", fat_updated$name))) %>% 
                                 select(contains(c("body_height")))

a<- gather(body_h, key= "groups", value="height", -ann)


ggplot(a, aes(x= height, fill=groups, color=groups)) + 
  geom_histogram(bins = 15, position="identity", alpha=0.5) 
  theme(legend.title = element_blank())+xlab("height")

body_w <- fat_updated %>% mutate(ann = as.numeric(sub("^.*_0", "", fat_updated$name))) %>% 
  select(contains(c("body_width")))

muscle_h <- fat_updated %>% mutate(ann = as.numeric(sub("^.*_0", "", fat_updated$name))) %>% 
  select(contains(c("muscle_height")))

muscle_w <- fat_updated %>% mutate(ann = as.numeric(sub("^.*_0", "", fat_updated$name))) %>% 
  select(contains(c("muscle_width")))

fat_h <- fat_updated %>% mutate(ann = as.numeric(sub("^.*_0", "", fat_updated$name))) %>% 
  select(contains(c("fat_height")))

fat_w <- fat_updated %>% mutate(ann = as.numeric(sub("^.*_0", "", fat_updated$name))) %>% 
  select(contains(c("fat_width")))

max_body_h<-c();max_body_w<-c();
max_mus_h<-c();max_mus_w<-c();
max_fat_h<-c();max_fat_w<-c();

for(i in c(1:106)){
  print(i)
  max_body_h<- append(max_body_h, colnames(body_h[which(body_h[i,] == max(body_h[i,]))]))
  max_body_w<- append(max_body_w, colnames(body_w[which(body_w[i,] == max(body_w[i,]))]))
  
  max_mus_h<- append(max_mus_h, colnames(muscle_h[which(muscle_h[i,] == max(muscle_h[i,]))]))
  max_mus_w<- append(max_mus_w, colnames(muscle_w[which(muscle_w[i,] == max(muscle_w[i,]))]))
  
  max_fat_h<- append(max_fat_h, colnames(fat_h[which(fat_h[i,] == max(fat_h[i,]))]))
  max_fat_w<- append(max_fat_w, colnames(fat_w[which(fat_w[i,] == max(fat_w[i,]))]))
}  

max_df<- data.frame(max_body_h, max_body_w, max_mus_h, max_mus_w, max_fat_h, max_fat_w)
