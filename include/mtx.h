//https://developer.apple.com/legacy/library/samplecode/glut/Listings/gle_vvector_h.html

/* ========================================================== */
/* determinant of matrix
 *
 * Computes determinant of matrix m, returning d
 */
 
#define DETERMINANT_3X3(d,m)                    \
{                               \
   d = m[0][0] * (m[1][1]*m[2][2] - m[1][2] * m[2][1]);     \
   d -= m[0][1] * (m[1][0]*m[2][2] - m[1][2] * m[2][0]);    \
   d += m[0][2] * (m[1][0]*m[2][1] - m[1][1] * m[2][0]);    \
}

/* ========================================================== */
/* compute adjoint of matrix and scale
 *
 * Computes adjoint of matrix m, scales it by s, returning a
 */
 
#define SCALE_ADJOINT_3X3(a,s,m)                \
{                               \
   a[0][0] = (s) * (m[1][1] * m[2][2] - m[1][2] * m[2][1]); \
   a[1][0] = (s) * (m[1][2] * m[2][0] - m[1][0] * m[2][2]); \
   a[2][0] = (s) * (m[1][0] * m[2][1] - m[1][1] * m[2][0]); \
                                \
   a[0][1] = (s) * (m[0][2] * m[2][1] - m[0][1] * m[2][2]); \
   a[1][1] = (s) * (m[0][0] * m[2][2] - m[0][2] * m[2][0]); \
   a[2][1] = (s) * (m[0][1] * m[2][0] - m[0][0] * m[2][1]); \
                                \
   a[0][2] = (s) * (m[0][1] * m[1][2] - m[0][2] * m[1][1]); \
   a[1][2] = (s) * (m[0][2] * m[1][0] - m[0][0] * m[1][2]); \
   a[2][2] = (s) * (m[0][0] * m[1][1] - m[0][1] * m[1][0]); \
}

/* ========================================================== */
/* inverse of matrix 
 *
 * Compute inverse of matrix a, returning determinant m and 
 * inverse b
 */
 
#define INVERT_3X3(b,det,a)         \
{                       \
   double tmp;                  \
   DETERMINANT_3X3 (det, a);            \
   tmp = 1.0 / (det);               \
   SCALE_ADJOINT_3X3 (b, tmp, a);       \
}
