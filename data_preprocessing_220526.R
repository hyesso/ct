#### scan landmark modeling

## load packages ---------------------------------------------------------
pkg_list<-c("readxl", "writexl", "dplyr", "scatterplot3d", "plotly", "rockchalk", "summarytools")
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
df_fin <- left_join(dist_df, cli_matched, by=c("ann"= "psm_no"))

str(df_fin)

view(dfSummary(cli_matched[, -1]))
view(cli_matched %>% group_by(sex)%>% dfSummary())

view(dfSummary(dist_df[,!(colnames(dist_df) %in% c("ann", "p_num"))]))


(dist_df_by_pnum <- stby(data      = select(dist_df, -c(ann,p_num)), 
                               INDICES   = dist_df$p_num, 
                               FUN       = descr, 
                               stats     = "common", 
                               transpose = TRUE))
attr(a, 'data_info') <- NULL
str(dist_df_by_pnum)
attr(a, "fn_call")
view(a)
attr(a, "format_info")$headings <- FALSE
attr(a, "format_info")$style

## read fat, muscle data( ; 220726) --------------------------------------
fat_df <- as.data.frame(read.csv("/Users/jinhyesu/my_project/rawData/220511_ct/measurement_pneumo.csv"))
str(fat_df)


