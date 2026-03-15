# Perspective trapezoidal correction maintaining real aspect ratio
# www.overfitting.net
# https://www.overfitting.net/2026/03/correccion-de-perspectiva-preservando.html

library(tiff)
library(Rcpp)


classify_aspect_ratio <- function(ratio, tolerance = 0.05) {
    # Orientation invariant
    ratio_norm <- if (ratio < 1) 1 / ratio else ratio
    
    # Canonical aspect ratios
    targets <- c(
        "1:1"  = 1/1,
        "3:2"  = 3/2,
        "4:3"  = 4/3,
        "16:9" = 16/9
    )
    
    # Relative deviation
    rel_dev <- abs(ratio_norm - targets) / targets
    idx <- which.min(rel_dev)
    
    if (rel_dev[idx] < tolerance) {
        return(names(targets)[idx])
    } else {
        return("")
    }
}


# MAGIC ASPECT RATIO/FOCAL LENGTH FUNCTION

# Main code based on Stack Overflow Java code found here:
# https://stackoverflow.com/questions/38285229/calculating-aspect-ratio-of-perspective-transform-destination-image
# https://stackoverflow.com/questions/1194352/proportions-of-a-perspective-deformed-rectangle/1222855#1222855

# Which is the implementation of the paper equations (section 4) by Zhengyou Zhang:
# "Whiteboard It! Convert Whiteboard Content into an Electronic Document"

# This function analyzes a given trapezoid (originally a rectangle in the real world) found on an image to calculate:
#  o The aspect ratio of the dimensions of the real world rectangle
#  o The focal length (in FF mm) that was used to capture the image
#  o The FOV in the x,y,diagonal axes
get_aspectratio_focallength_FOV <- function(width, height, x, y) {
    # width and height are the dimensions of the image in pixels
    # x and y are vectors containing the coordinates of the four vertices of the trapezoid
    # expressed in (top-left, bottom-left, bottom-right, top-right) order
    
    require(png)
    
    # Obtain sensor format aspect ratio
    image_aspect_ratio = width / height
    aspect_ratio_label = classify_aspect_ratio(image_aspect_ratio)
    
    # Plot of trapezoid for checking
    img=matrix(0, ncol=width, nrow=height)
    draw_line <- function(M, x0, y0, x1, y1) {
        n <- max(abs(x1 - x0), abs(y1 - y0)) + 1
        xs <- round(seq(x0, x1, length.out = n))
        ys <- round(seq(y0, y1, length.out = n))
        M[cbind(ys, xs)] <- 1
        M
    }
    # Plot 3-pixel width lines
    for (dx in c(-1,0,1)) {
        for (dy in c(-1,0,1)) {
            img=draw_line(img, x[1]+dx, y[1]+dy, x[2]+dx, y[2]+dy)
            img=draw_line(img, x[2]+dx, y[2]+dy, x[3]+dx, y[3]+dy)
            img=draw_line(img, x[3]+dx, y[3]+dy, x[4]+dx, y[4]+dy)
            img=draw_line(img, x[4]+dx, y[4]+dy, x[1]+dx, y[1]+dy)
        }
    }
    writePNG(img, "trapezoid.png")
    
    
    # Let m1x,m1y...m4x,m4y be the (x,y) pixel coordinates of the 4 corners of the detected quadrangle
    # i.e. (m1x, m1y) are the cordinates of the first corner, (m2x, m2y) of the second corner and so on
    # m1, m2, m3, m4 follow the order: top-left, top-right, bottom-left, bottom-right
    # Let u0, v0 be the pixel coordinates of the principal point of the image
    # which for a normal camera will be the centre of the image: i.e. u0=IMAGEWIDTH/2, v0 =IMAGEHEIGHT/2
    # This assumption does not hold if the image has been cropped asymmetrically
    
    # Transform the image so the principal point is at (0,0) which makes the following equations much easier
    # Image center (principal point assumption)
    u0 <- width / 2
    v0 <- height / 2
    
    # Shift coordinates
    # NOTE: m1, m2, m3, m4 follow the order: top-left, top-right, bottom-left, bottom-right
    # as in the original paper "Whiteboard scanning and image enhancement" by Zhengyou Zhang
    m1x <- x[1] - u0; m1y <- y[1] - v0  # top-left
    m3x <- x[2] - u0; m3y <- y[2] - v0  # bottom-left
    m4x <- x[3] - u0; m4y <- y[3] - v0  # bottom-right
    m2x <- x[4] - u0; m2y <- y[4] - v0  # top-right
    
    # Compute k2
    k2 <- ((m1y - m4y) * m3x - (m1x - m4x) * m3y + m1x * m4y - m1y * m4x) /
          ((m2y - m4y) * m3x - (m2x - m4x) * m3y + m2x * m4y - m2y * m4x)
    # Compute k3
    k3 <- ((m1y - m4y) * m2x - (m1x - m4x) * m2y + m1x * m4y - m1y * m4x) /
          ((m3y - m4y) * m2x - (m3x - m4x) * m2y + m3x * m4y - m3y * m4x)
    # Focal length in pixels
    focal_length_px_squared <- -((k3 * m3y - m1y) * (k2 * m2y - m1y) + (k3 * m3x - m1x) * (k2 * m2x - m1x)) /
                                ((k3 - 1) * (k2 - 1))
    focal_length_px <- sqrt(focal_length_px_squared)
    
    
    # Calculate aspect ratio W/H
    
    # If k2==1 AND k3==1 the focal length equation is not solvable 
    # but the focal length is not needed to calculate the aspect ratio in that case
    # k2 and k3 become 1 when the rectangle is not distorted by perspective
    # The projective formula becomes numerically unstable when (k2−1)*(k3−1) -> 0 (affine case)
    # Check if projection is affine
    
    # D <- (k2 - 1)*(k3 - 1)
    # if (abs(k2 - 1) < eps && abs(k3 - 1) < eps) {
    #     # Case 1: affine
    # } else if (abs(D) < eps) {
    #     # Case 2: single vanishing direction
    # } else {
    #     # Case 3: full projective
    # }
    # # Both k’s near 1 → only Euclidean aspect ratio valid.
    # # Exactly one near 1 → nothing metrically identifiable.
    # # Neither near 1 → full recovery possible.
    
    eps <- 1e-6    
    if (abs(k2 - 1) < eps && abs(k3 - 1) < eps) {
        #if (k2 == 1 && k3 == 1) {
        # -------------------------------
        # Affine / weak-perspective case
        # -------------------------------
        # Parallel lines are approximately parallel in image
        # Foreshortening negligible, no projective correction needed
        cat("WARNING: affine / weak-perspective so aspect ratio derived from Euclidean distances, unstable focal length/FOV\n")
        whRatio <- sqrt(
            ((m2y - m1y)^2 + (m2x - m1x)^2) /  # observed width  (squared Euclidean distance between points m1 and m2)
            ((m3y - m1y)^2 + (m3x - m1x)^2)    # observed height (squared Euclidean distance between points m1 and m3)
        )
        # Focal length and FOV values will turn to NaN's
    } else {
        # ---------------------------------------
        # Full projective perspective case
        # ---------------------------------------
        # Corrects trapezoidal distortion using k2, k3, and focal length
        # Computes metric length ratio under perspective projection
        whRatio <- sqrt(
            ((k2 - 1)^2 +
                 (k2 * m2y - m1y)^2 / focal_length_px_squared +
                 (k2 * m2x - m1x)^2 / focal_length_px_squared) /
            ((k3 - 1)^2 +
                 (k3 * m3y - m1y)^2 / focal_length_px_squared +
                 (k3 * m3x - m1x)^2 / focal_length_px_squared)
        )
    }
    
    # FOV calculations (radians)
    FOV_x  <- 2 * atan((width  / 2) / focal_length_px)
    FOV_y  <- 2 * atan((height / 2) / focal_length_px)
    FOV_d  <- 2 * atan((sqrt(width^2 + height^2) / 2) / focal_length_px)
    
    rad2deg <- function(r) r * 180 / pi
    
    FOV_x_deg  <- rad2deg(FOV_x)
    FOV_y_deg  <- rad2deg(FOV_y)
    FOV_d_deg  <- rad2deg(FOV_d)
    
    # ---- FF equivalent focal length (mm) ----
    FF_diag_mm <- sqrt(36^2 + 24^2)  # FF sensor dimensions in mm (nominal FF sensor)
    # FF_diag_mm <- sqrt(35.8^2 + 23.9^2)  # FF sensor dimensions in mm (Sony A7 II)
    focal_length_FF_mm <- (FF_diag_mm / 2) / tan(FOV_d / 2)
    
    # ---- PRINT RESULTS ----
    cat(sprintf("Image aspect ratio (W/H): %.2f %s\n", image_aspect_ratio,
                ifelse(aspect_ratio_label=="", "", paste0("[",aspect_ratio_label,"]"))))
    cat(sprintf("Rectangle aspect ratio (W/H): %.2f\n", whRatio))
    cat(sprintf("Focal length FF: %.2fmm\n", focal_length_FF_mm))
    cat(sprintf("FOV X/Y/Diag: %.1f° / %.1f° / %.1f°\n", FOV_x_deg, FOV_y_deg, FOV_d_deg))
    
    return(list(
        image_aspect_ratio = image_aspect_ratio,
        rectangle_aspect_ratio = whRatio,
        # focal_length_pixels = focal_length_px,
        focal_length_FF_mm = focal_length_FF_mm,
        FOV_x_deg = FOV_x_deg,
        FOV_y_deg = FOV_y_deg,
        FOV_diag_deg = FOV_d_deg
    ))
}


# KEYSTONE CORRECTION FUNCTIONS

# Function that models the (xd,yd) -> (xu,yu) transformation through k (8 coefficients)
# NOTE: when using it, we'll swap the distorted and undistorted trapezoids because
# we want to model the transformation FROM CORRECTED coords (DST) -> TO UNCORRECTED coords (ORG)
solve.keystone = function(xd, yd, xu, yu) {
    # Solve 8 equations linear system: A * k = b -> k = inv(A) * b
    A=matrix(nrow=8, ncol=8)
    A[1,]=c(xd[1], yd[1], 1, 0,     0,     0, -xd[1]*xu[1], -yd[1]*xu[1])
    A[2,]=c(0,     0,     0, xd[1], yd[1], 1, -xd[1]*yu[1], -yd[1]*yu[1])
    A[3,]=c(xd[2], yd[2], 1, 0,     0,     0, -xd[2]*xu[2], -yd[2]*xu[2])
    A[4,]=c(0,     0,     0, xd[2], yd[2], 1, -xd[2]*yu[2], -yd[2]*yu[2])
    A[5,]=c(xd[3], yd[3], 1, 0,     0,     0, -xd[3]*xu[3], -yd[3]*xu[3])
    A[6,]=c(0,     0,     0, xd[3], yd[3], 1, -xd[3]*yu[3], -yd[3]*yu[3])
    A[7,]=c(xd[4], yd[4], 1, 0,     0,     0, -xd[4]*xu[4], -yd[4]*xu[4])
    A[8,]=c(0,     0,     0, xd[4], yd[4], 1, -xd[4]*yu[4], -yd[4]*yu[4])
    
    b=as.matrix(c(xu[1], yu[1], xu[2], yu[2], xu[3], yu[3], xu[4], yu[4]))

    k=solve(A, b)  # equivalent to inv(A) * b = solve(A) %*% b
    
    return(k)
}

# Undo distortion function
undo.keystone = function(xd, yd, k) {
    xu=(k[1]*xd+k[2]*yd+k[3]) / (k[7]*xd+k[8]*yd+1)
    yu=(k[4]*xd+k[5]*yd+k[6]) / (k[7]*xd+k[8]*yd+1)
    return(c(xu, yu))  # return pair (xu, yu)
}

# Keystone correction improvement and optimization:
# 1. C++ compilation of (x,y) nested loops over the output image
# 2. Bilinear interpolation instead of nearest neighbour interpolation
sourceCpp("keystone.cpp")  # keystone_correct_cpp() function



###################################################
# CALCULATE REAL ASPECT RATIO + PERSPECTIVE CORRECTION USING THAT ASPECT RATIO

# The result will be 100% the same as having used a Tilt-shit lens both in perspective and aspect ratio


#################################
# 1. Fuji X-S10 con 11mm (16,5mm eq.)

width=6240; height=4160
xu=c(1516,  982, 6181, 5254)  # top-left, bottom-left, bottom-right, top-right
yu=c(2391, 3609, 3372, 1873)
AR=get_aspectratio_focallength_FOV(width, height, xu, yu)

aspect_ratio=AR$rectangle_aspect_ratio


# Distorted points (source)
imgd=readTIFF("building.tif")

# Undistorted points (destination)
# Ad-hoc scaling/shifting adjustments to improve the corrected area
posx1=(xu[1]+xu[2])/2/1.7+1000  # top-left corner
posy1=(yu[1]+yu[4])/2/1.4+1300

posx2=(xu[3]+xu[4])/2/1.7+1000  # bottom-right corner
posy2=posy1+(posx2-posx1)/aspect_ratio  # THIS SENTENCE forces the desired aspect ratio

xd=c(posx1, posx1, posx2, posx2)
yd=c(posy1, posy2, posy2, posy1)


# Calculate k (8 coefficients) of the keystone correction 
k=solve.keystone(xd, yd, xu, yu)  # models the (xd,yd) -> (xu,yu) transformation

# Check
for (i in 1:4) print(undo.keystone(xd[i], yd[i], k))

# Plot trapezoid correction
png(paste0("correctionbuilding.png"), width=512, height=400)
    plot(c(xd, xd[1]), c(yd, yd[1]), type='l', col='red', asp=1,
         xlab='X', ylab='Y', xlim=c(1, ncol(imgd)), ylim=c(nrow(imgd), 1))
    lines(c(xu, xu[1]), c(yu, yu[1]), type='l', col='blue')
    for (i in 1:4) {
        lines(c(xd[i], xu[i]), c(yd[i], yu[i]), type='l', lty=3, col='darkgray')
    }
    abline(h=c(1,nrow(imgd)), v=c(1,ncol(imgd)))
dev.off()


# Write on imgc_cpp the keystone corrected version of imgd
imgc_cpp=keystone_correct_cpp(imgd*0, imgd, as.numeric(k))
writeTIFF(imgc_cpp, "correctedbuilding_aspectratio.tif", bits.per.sample=16)



#################################
# 2. Sony A7 II 50mm

width=6000; height=4000
xu=c(4207, 2511, 4181, 5780)  # top-left, bottom-left, bottom-right, top-right
yu=c( 857, 2244, 3336, 1649)
AR=get_aspectratio_focallength_FOV(width, height, xu, yu)

aspect_ratio=AR$rectangle_aspect_ratio


# Distorted points (source)
imgd=readTIFF("laminas.tif")

# Undistorted points (destination)
# Ad-hoc scaling/shifting adjustments to improve the corrected area
posx1=(xu[1]+xu[2])/2-800  # top-left corner
posy1=(yu[1]+yu[4])/2-1000

posx2=(xu[3]+xu[4])/2-800  # bottom-right corner
posy2=posy1+(posx2-posx1)/aspect_ratio  # THIS SENTENCE forces the desired aspect ratio

xd=c(posx1, posx1, posx2, posx2)
yd=c(posy1, posy2, posy2, posy1)


# Calculate k (8 coefficients) of the keystone correction 
k=solve.keystone(xd, yd, xu, yu)  # models the (xd,yd) -> (xu,yu) transformation

# Check
for (i in 1:4) print(undo.keystone(xd[i], yd[i], k))

# Plot trapezoid correction
png(paste0("correctionlaminas.png"), width=512, height=400)
    plot(c(xd, xd[1]), c(yd, yd[1]), type='l', col='red', asp=1,
         xlab='X', ylab='Y', xlim=c(1, ncol(imgd)), ylim=c(nrow(imgd), 1))
    lines(c(xu, xu[1]), c(yu, yu[1]), type='l', col='blue')
    for (i in 1:4) {
        lines(c(xd[i], xu[i]), c(yd[i], yu[i]), type='l', lty=3, col='darkgray')
    }
    abline(h=c(1,nrow(imgd)), v=c(1,ncol(imgd)))
dev.off()


# Write on imgc_cpp the keystone corrected version of imgd
imgc_cpp=keystone_correct_cpp(imgd*0, imgd, as.numeric(k))
writeTIFF(imgc_cpp, "correctedlaminas_aspectratio.tif", bits.per.sample=16)
