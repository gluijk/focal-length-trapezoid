# Derive focal length, aspect ratio and FOV from an image with some trapezoid
# www.overfitting.net
# https://www.overfitting.net/2026/02/calculando-la-distancia-focal-con-que.html


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



###################################################
# EXAMPLES

# (x,y) coordinates expressed as: top-left, bottom-left, bottom-right, top-right


# SYNTHETIC EXAMPLES: SPECIAL CASES

# Affine projection: square -> aspect ratio is readily available, focal length cannot be calculated
width=1000; height=1000
x=c(100, 100, 900, 900)
y=c(100, 900, 900, 100)
AR=get_aspectratio_focallength_FOV(width, height, x, y)  # AR=1, rest of values are NaN's

# Affine projection: rectangle -> aspect ratio is readily available, focal length cannot be calculated
width=1000; height=1000
x=c(100, 100, 900, 900)
y=c(100+100, 900-100, 900-100, 100+100)
AR=get_aspectratio_focallength_FOV(width, height, x, y)  # AR=1.33, rest of values are NaN's

# Single vanishing point -> nothing can be calculated
width=1000; height=1000
x=c(100, 400, 600, 900)
y=c(100, 900, 900, 100)
AR=get_aspectratio_focallength_FOV(width, height, x, y)



# SYNTHETIC EXAMPLES: CHECK ALGORITHM

# 24mm scene
# Rhinoceros 1
width=1434; height=956
x=c(217, 220, 830, 887)
y=c(138, 448, 379, 193)
AR=get_aspectratio_focallength_FOV(width, height, x, y)  # 3,02 (3), 23.89mm (24)

# Rhinoceros 2
width=1434; height=956
x=c(217, 220, 524, 612)
y=c(138, 448, 908, 484)
AR=get_aspectratio_focallength_FOV(width, height, x, y)  # 1.99 (2), 23.89mm (24)

# Rhinoceros 3
width=1434; height=956
x=c(220, 524, 1170, 830)
y=c(448, 908,  586, 379)
AR=get_aspectratio_focallength_FOV(width, height, x, y)  # 1.51 (1.5), 23.97mm (24)


# 50mm scene
# Rhinoceros 1
width=1434; height=956
x=c(106, 221, 775, 713)
y=c(427, 743, 382, 141)
AR=get_aspectratio_focallength_FOV(width, height, x, y)  # 3.02 (3), 50.97mm (50)

# Rhinoceros 2
width=1434; height=956
x=c(106, 221, 770, 659)
y=c(427, 743, 802, 400)
AR=get_aspectratio_focallength_FOV(width, height, x, y)  # 2.02 (2), 50.80mm (50)

# Rhinoceros 3
width=1434; height=956
x=c(221, 770, 1246, 775)
y=c(743, 802, 357, 382)
AR=get_aspectratio_focallength_FOV(width, height, x, y)  # 1.51 (1.5), 49.89mm (50)



# REAL WORLD EXAMPLES

# Panasonic S1R con 14-28mm a 14mm → 14,60mm estimados (several rectangles)
width=1920; height=1281
x=c(1207,  1158, 1786, 1684)
y=c(585, 775, 811, 554)
AR=get_aspectratio_focallength_FOV(width, height, x, y)  # 14.79mm

width=1920; height=1281
x=c(862, 315, 748, 999)
y=c( 77, 267, 809, 438)
AR=get_aspectratio_focallength_FOV(width, height, x, y)  # 14.50mm

width=1920; height=1281
x=c(1199, 1136, 1437, 1434)
y=c( 866, 1153, 1256, 897)
AR=get_aspectratio_focallength_FOV(width, height, x, y)  # 14.48mm

width=1920; height=1281
x=c(986,  866, 1093, 1172)
y=c(837, 1124, 1210, 859)
AR=get_aspectratio_focallength_FOV(width, height, x, y)  # 14.34mm

width=1920; height=1281
x=c(1569, 1591, 1950, 1865)
y=c( 985, 1178, 1277, 1035)
AR=get_aspectratio_focallength_FOV(width, height, x, y)  # 14.31mm

width=1920; height=1281
x=c( 682,  554,  638,  745)
y=c(1024, 1231, 1277, 1079)
AR=get_aspectratio_focallength_FOV(width, height, x, y)  # 15.06mm


# Fuji X-S10 con 16,5mm (eq.) → 16,84mm estimados
width=6240; height=4160
x=c(1516,  982, 6181, 5254)
y=c(2391, 3609, 3372, 1873)
AR=get_aspectratio_focallength_FOV(width, height, x, y)


# Ricoh GR IV con 28mm (eq.) → 28,26mm estimados
width=1280; height=1920
x=c(557,  80, 1028, 819)
y=c(493, 1826, 1525, 322)
AR=get_aspectratio_focallength_FOV(width, height, x, y)


# Nikon Z8 con 45mm → 46,51mm estimados
width=1920; height=1280
x=c(766,  752, 1395, 1376)
y=c(209, 684, 978, 652)
AR=get_aspectratio_focallength_FOV(width, height, x, y)


# Canon R5 II con 85mm → 82,80mm estimados
width=1920; height=1281
x=c(190,  150, 1049, 1077)
y=c(550, 1173, 1148, 283)
AR=get_aspectratio_focallength_FOV(width, height, x, y)


# Sony A7R V con 100mm → 95,87mm estimados
width=1920; height=1280
x=c(959,  958, 1595, 1581)
y=c(734, 1001, 1273, 1030)
AR=get_aspectratio_focallength_FOV(width, height, x, y)


# Sony A7R V con 200mm → 169,19mm estimados
width=4000; height=2667
x=c(3468, 3475, 3926, 3914)
y=c( 569, 1255, 1470,  804)
AR=get_aspectratio_focallength_FOV(width, height, x, y)


# 'Lunch atop a Skyscraper' → 43,92mm estimados
width=1920; height=1481
x=c(55,  64, 1774, 1783)
y=c(635, 734, 912, 780)
AR=get_aspectratio_focallength_FOV(width, height, x, y)


# Quake III Arena LEFT → 19,04mm estimados
width=2560; height=1440
x=c(348,  247, 736, 777)
y=c(475, 1160, 1111, 595)
AR=get_aspectratio_focallength_FOV(width, height, x, y)

# Quake III Arena RIGHT → 19,06mm estimados
width=2560; height=1440
x=c(1785,  1815, 2399, 2323)
y=c(750, 1149, 1190, 703)
AR=get_aspectratio_focallength_FOV(width, height, x, y)


# Escher LEFT → 30,33mm estimados
width=1159; height=1096
x=c(202,  60, 608, 608)
y=c(386, 945, 1035, 576)
AR=get_aspectratio_focallength_FOV(width, height, x, y)

# Escher - RIGHT -> 31.34mm estimados
width=1159; height=1096
x=c(608,  608, 1155, 952)
y=c(271, 912, 773, 30)
AR=get_aspectratio_focallength_FOV(width, height, x, y)

