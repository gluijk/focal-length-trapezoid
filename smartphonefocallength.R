# Calculating the focal length of a mobile phone's camera
# www.overfitting.net
# https://www.overfitting.net/2026/03/calculando-la-distancia-focal-de-la.html

library(Rcpp)


# Quick Bresenham algorithm to draw lines over a matrix
# img is modified directly (no matrix copying)
cppFunction('
    void draw_line_bresenham(NumericMatrix img, int x0, int y0, int x1, int y1)
    {
        x0--; y0--; x1--; y1--;  // convert R indexing -> C indexing
    
        int dx = std::abs(x1 - x0);
        int dy = std::abs(y1 - y0);
    
        int sx = (x0 < x1) ? 1 : -1;
        int sy = (y0 < y1) ? 1 : -1;
    
        int err = dx - dy;
    
        while(!(x0 == x1 && y0 == y1))  // never plot last pixel to prevent plotting corner pixels twice
        {
            if(x0 >= 0 && x0 < img.ncol() && y0 >= 0 && y0 < img.nrow()) img(y0, x0) += 1.0;
    
            int e2 = 2 * err;
            if(e2 > -dy)
            {
                err -= dy;
                x0 += sx;
            }
            if(e2 < dx)
            {
                err += dx;
                y0 += sy;
            }
        }
    }
')



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
    # expressed in top-left, bottom-left, bottom-right, top-right order
    
    require(png)
    
    # Obtain image format aspect ratio
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
    # 3-pixel width lines
    for (dx in c(-1*0,0,1*0)) {
        for (dy in c(-1*0,0,1*0)) {
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
    # Let u0, v0 be the pixel coordinates of the principal point of the image (optical axis)
    # which for a normal camera will be the centre of the image: i.e. u0=IMAGEWIDTH/2, v0 =IMAGEHEIGHT/2
    # This assumption does not hold if the image has been cropped asymmetrically or shifted
    
    # Transform the image so the principal point is at (0,0) which makes the following equations much easier
    # Image center (principal point assumption)
    u0 <- width / 2
    v0 <- height / 2
    
    # Shift coordinates
    # NOTE: m1, m2, m3, m4 follow the order: top-left, top-right, bottom-left, bottom-right
    # as in the original paper "Whiteboard scanning and image enhancement" by Zhengyou Zhang
    # but (x,y) follow top-left, bottom-left, bottom-right, top-right
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
    diag_px <- sqrt(width^2 + height^2) 
    FOV_x  <- 2 * atan((width  / 2) / focal_length_px)
    FOV_y  <- 2 * atan((height / 2) / focal_length_px)
    FOV_d  <- 2 * atan((diag_px / 2) / focal_length_px)
    
    rad2deg <- function(r) r * 180 / pi
    FOV_x_deg  <- rad2deg(FOV_x)
    FOV_y_deg  <- rad2deg(FOV_y)
    FOV_d_deg  <- rad2deg(FOV_d)
    
    # ---- FF equivalent focal length (mm) ----
    FF_diag_mm <- sqrt(36^2 + 24^2)  # FF sensor dimensions in mm (nominal FF sensor)
    # FF_diag_mm <- sqrt(35.8^2 + 23.9^2)  # FF sensor dimensions in mm (Sony A7 II)
    focal_length_FF_mm <- (FF_diag_mm / 2) / tan(FOV_d / 2)
    # focal_length_FF_mm = focal_length_px * FF_diag_mm / diag_px

    # ---- PRINT RESULTS ----
    cat(sprintf("Image aspect ratio (W/H): %.2f %s\n", image_aspect_ratio,
                ifelse(aspect_ratio_label=="", "", paste0("[",aspect_ratio_label,"]"))))
    cat(sprintf("Rectangle aspect ratio (W/H): %.2f\n", whRatio))
    cat(sprintf("Focal length FF: %.2fmm\n", focal_length_FF_mm))
    cat(sprintf("FOV X/Y/Diag: %.1f° / %.1f° / %.1f°\n", FOV_x_deg, FOV_y_deg, FOV_d_deg))

    return(list(
        image_aspect_ratio = image_aspect_ratio,
        rectangle_aspect_ratio = whRatio,
        # focal_length_pixels = focal_length_px,  # image dependant value, not relevant
        focal_length_FF_mm = focal_length_FF_mm,
        FOV_x_deg = FOV_x_deg,
        FOV_y_deg = FOV_y_deg,
        FOV_diag_deg = FOV_d_deg
    ))
}


# Montecarlo version of get_aspectratio_focallength_FOV()
# It jitters x and y positions to study aspect ratio and focal length sensitivity
get_aspectratio_focallength_FOV_Montecarlo <- function(width, height, x, y, N = 10000, sd = 3, jitter = 'normal') {
    # width and height are the dimensions of the image in pixels
    # x and y are vectors containing the coordinates of the four vertices of the trapezoid
    # expressed in top-left, bottom-left, bottom-right, top-right order
    
    require(tiff)
    
    # Precalculations
    FF_diag_mm <- sqrt(36^2 + 24^2)  # FF sensor dimensions in mm (nominal FF sensor)
    diag_px <- sqrt(width^2 + height^2) 
    rad2deg <- function(r) r * 180 / pi

    # Obtain image format aspect ratio
    image_aspect_ratio = width / height
    aspect_ratio_label = classify_aspect_ratio(image_aspect_ratio)
    
    # Transform the image so the principal point is at (0,0) which makes the following equations much easier
    # Image center (principal point assumption)
    u0 <- width / 2
    v0 <- height / 2
    
    # Plot of trapezoid for checking
    img=matrix(0, ncol=width, nrow=height)

    # Build jitter matrix (first column preserves the original input values
    # so N=1 means the single point case)
    if (jitter == 'normal') {  # normal distribution
        Mx <- x + cbind(0, matrix(rnorm(4*(N-1), sd=sd), nrow=4))
        My <- y + cbind(0, matrix(rnorm(4*(N-1), sd=sd), nrow=4))
    } else {  # uniform distribution
        Mx <- x + cbind(0, matrix(runif(4*(N-1), min=-sd, max=sd), nrow=4))
        My <- y + cbind(0, matrix(runif(4*(N-1), min=-sd, max=sd), nrow=4))
    }
    
    whRatioAC=c()
    focal_length_FF_mmAC=c()
    for (i in 1:N) {
        
        # 1-pixel width lines
        draw_line_bresenham(img, round(Mx[1,i]), round(My[1,i]), round(Mx[2,i]), round(My[2,i]))
        draw_line_bresenham(img, round(Mx[2,i]), round(My[2,i]), round(Mx[3,i]), round(My[3,i]))
        draw_line_bresenham(img, round(Mx[3,i]), round(My[3,i]), round(Mx[4,i]), round(My[4,i]))
        draw_line_bresenham(img, round(Mx[4,i]), round(My[4,i]), round(Mx[1,i]), round(My[1,i]))
        
        # Shift coordinates
        m1x <- Mx[1,i] - u0; m1y <- My[1,i] - v0  # top-left
        m3x <- Mx[2,i] - u0; m3y <- My[2,i] - v0  # bottom-left
        m4x <- Mx[3,i] - u0; m4y <- My[3,i] - v0  # bottom-right
        m2x <- Mx[4,i] - u0; m2y <- My[4,i] - v0  # top-right
        
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
        # Full projective perspective case
        whRatio <- sqrt(
            ((k2 - 1)^2 +
                 (k2 * m2y - m1y)^2 / focal_length_px_squared +
                 (k2 * m2x - m1x)^2 / focal_length_px_squared) /
            ((k3 - 1)^2 +
                 (k3 * m3y - m1y)^2 / focal_length_px_squared +
                 (k3 * m3x - m1x)^2 / focal_length_px_squared)
        )
        
        focal_length_FF_mm = focal_length_px * FF_diag_mm / diag_px
        
        # Store valid calculations
        if (!is.nan(whRatio)) whRatioAC=c(whRatioAC, whRatio)
        if (!is.nan(focal_length_FF_mm)) focal_length_FF_mmAC=c(focal_length_FF_mmAC, focal_length_FF_mm)

    }  # end of Montecarlo loop
    writeTIFF(img/max(img), paste0("trapezoids_montecarlo", N, ".tiff"), bits.per.sample = 16)
    
    # Plotting histograms
    par(mfrow = c(2,1))
    
    # Focal length
    median_fl=median(focal_length_FF_mmAC)
    values=focal_length_FF_mmAC[focal_length_FF_mmAC <= 5*median_fl]
    xlim=c(min(values), max(values))
    hist(values, breaks=800, xlim=xlim, xlab="mm", ylab="", cex.main=0.9,
         main=paste0("Focal length FF: USR=", round(focal_length_FF_mmAC[1], 2),
                     "mm, Med=", round(median_fl, 2), "mm"))
    abline(v=focal_length_FF_mmAC[1], col='red')
    abline(v=median_fl, col='green')
    
    # Aspect ratio
    median_ar=median(whRatioAC)
    values=whRatioAC[whRatioAC <= 5*median_ar]
    xlim=c(min(values), max(values))
    hist(values, breaks=800, xlim=xlim, xlab="Aspect ratio", ylab="", cex.main=0.9,
         main=paste0("Aspect ratio: USR=", round(whRatioAC[1], 2),
                     ", Med=", round(median_ar, 2)))
    abline(v=whRatioAC[1], col='red')
    abline(v=median_ar, col='green')
    
    cat(sprintf("Stable focal length calculated values: %d / %.0f (%.1f%%)\n",
                length(focal_length_FF_mmAC), N,
                length(focal_length_FF_mmAC)/N*100))
    cat(sprintf("Stable aspect ratio calculated values: %d / %.0f (%.1f%%)\n",
                length(whRatioAC), N,
                length(whRatioAC)/N*100))

    return(list(
        image_aspect_ratio = image_aspect_ratio,
        rectangle_aspect_ratio = whRatioAC,
        focal_length_FF_mm = focal_length_FF_mmAC
    ))
}


###################################################
# EXAMPLES

# (x,y) coordinates expressed as: top-left, bottom-left, bottom-right, top-right
# We set sd as 1/2000 of its diagonal


# Main camera Samsung S20 FE (26mm eq., AR subject=122.4/95.6)
width=2993; height=3991
diag=(width^2+height^2)^0.5; sd=diag/2000
x=c(327,   60, 1784, 2875)
y=c(480, 2340, 3957,  484)
AR=get_aspectratio_focallength_FOV(width, height, x, y)
AR2=get_aspectratio_focallength_FOV_Montecarlo(width, height, x, y, N=500000, sd=sd)


# Main camera Samsung S24 FE (24mm eq.)
width=4039; height=3029
diag=(width^2+height^2)^0.5; sd=diag/2000
x=c(723, 1000, 3566, 2894)
y=c(314, 2584, 2281, 1163)
AR=get_aspectratio_focallength_FOV(width, height, x, y)
AR2=get_aspectratio_focallength_FOV_Montecarlo(width, height, x, y, N=500000, sd=sd)


# Main camera Samsung S24 Ultra (23mm eq.)
width=3960; height=2970
diag=(width^2+height^2)^0.5; sd=diag/2000
x=c(299, 1430, 3720, 2574)
y=c(959, 2566,  860,  231)
AR=get_aspectratio_focallength_FOV(width, height, x, y)
AR2=get_aspectratio_focallength_FOV_Montecarlo(width, height, x, y, N=500000, sd=5, random = 'unif')
