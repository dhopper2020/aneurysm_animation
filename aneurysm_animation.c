#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define UNUSED __attribute__((unused))
#define ROW 1000
#define COL 1000
#define MAX 1000000
#define NORMALS

double C441(double f)
{
   return ceil(f-0.00001);
}

double F441(double f)
{
   return floor(f+0.00001);
}

typedef struct
{
   double         X[3];
   double         Y[3];
   double         Z[3];
   double         color[3][3]; // color[2][0] is for V2, red channel
#ifdef NORMALS
   double         normals[3][3]; // normals[2][0] is for V2, x-component
   double         shading[3];
#endif
} Triangle;

typedef struct
{
   int numTriangles;
   Triangle *triangles;
} TriangleList;

typedef struct camera
{
    double near, far;
    double angle;
    double position[3];
    double focus[3];
    double up[3];
} Camera;

char *Read3Numbers(char *tmp, double *v1, double *v2, double *v3)
{
    *v1 = atof(tmp);
    while (*tmp != ' ')
       tmp++;
    tmp++; /* space */
    *v2 = atof(tmp);
    while (*tmp != ' ')
       tmp++;
    tmp++; /* space */
    *v3 = atof(tmp);
    while (*tmp != ' ' && *tmp != '\n')
       tmp++;
    return tmp;
}

TriangleList *Get3DTriangles()
{
   FILE *f = fopen("ws_tris.txt", "r");
   if (f == NULL)
   {
       fprintf(stderr, "You must place the ws_tris.txt file in the current directory.\n");
       exit(EXIT_FAILURE);
   }
   fseek(f, 0, SEEK_END);
   int numBytes = ftell(f);
   fseek(f, 0, SEEK_SET);
   if (numBytes != 3892295)
   {
       fprintf(stderr, "Your ws_tris.txt file is corrupted.  It should be 3892295 bytes, but you have %d.\n", numBytes);
       exit(EXIT_FAILURE);
   }

   char *buffer = (char *) malloc(numBytes);
   if (buffer == NULL)
   {
       fprintf(stderr, "Unable to allocate enough memory to load file.\n");
       exit(EXIT_FAILURE);
   }
   
   fread(buffer, sizeof(char), numBytes, f);

   char *tmp = buffer;
   int numTriangles = atoi(tmp);
   while (*tmp != '\n')
       tmp++;
   tmp++;
 
   if (numTriangles != 14702)
   {
       fprintf(stderr, "Issue with reading file -- can't establish number of triangles.\n");
       exit(EXIT_FAILURE);
   }

   TriangleList *tl = (TriangleList *) malloc(sizeof(TriangleList));
   tl->numTriangles = numTriangles;
   tl->triangles = (Triangle *) malloc(sizeof(Triangle)*tl->numTriangles);

   for (int i = 0 ; i < tl->numTriangles ; i++)
   {
       for (int j = 0 ; j < 3 ; j++)
       {
           double x, y, z;
           double r, g, b;
           double normals[3];

           tmp = Read3Numbers(tmp, &x, &y, &z);
           tmp += 3; /* space+slash+space */
           tmp = Read3Numbers(tmp, &r, &g, &b);
           tmp += 3; /* space+slash+space */
           tmp = Read3Numbers(tmp, normals+0, normals+1, normals+2);
           tmp++;    /* newline */

           tl->triangles[i].X[j] = x;
           tl->triangles[i].Y[j] = y;
           tl->triangles[i].Z[j] = z;
           tl->triangles[i].color[j][0] = r;
           tl->triangles[i].color[j][1] = g;
           tl->triangles[i].color[j][2] = b;
#ifdef NORMALS
           tl->triangles[i].normals[j][0] = normals[0];
           tl->triangles[i].normals[j][1] = normals[1];
           tl->triangles[i].normals[j][2] = normals[2];
#endif
       }
   }

   free(buffer);
   return tl;
}

typedef struct pixel
{
   unsigned char one;
   unsigned char two;
   unsigned char three;
} Pixel;

typedef struct image
{
   Pixel *array;
   double *depth_array;
}Image;

void convert(Image *image, char* file_name)
{
    FILE *fp;
    fp = fopen(file_name, "w");
    
    fprintf(fp, "%s", "P6\n1000 1000\n255\n");

    for (int x = ROW - 1; x >= 0; x--)
    {
        for (int y = 0; y < COL; y++)
        {
            fprintf(fp, "%c", image->array[x + COL * y].one);
            fprintf(fp, "%c", image->array[x + COL * y].two);
            fprintf(fp, "%c", image->array[x + COL * y].three);
        }
    }

    fclose(fp);
}

typedef struct lightingparameters
{
   double lightDir[3]; // The direction of the light source
   double Ka;           // The coefficient for ambient lighting.
   double Kd;           // The coefficient for diffuse lighting.
   double Ks;           // The coefficient for specular lighting.
   double alpha;        // The exponent term for specular lighting.
} LightingParameters;

LightingParameters GetLighting(Camera c)
{
   LightingParameters lp;
   lp.Ka = 0.3;
   lp.Kd = 0.7;
   lp.Ks = 2.8;
   lp.alpha = 50.5;
   lp.lightDir[0] = c.position[0]-c.focus[0];
   lp.lightDir[1] = c.position[1]-c.focus[1];
   lp.lightDir[2] = c.position[2]-c.focus[2];
   double mag = sqrt(lp.lightDir[0]*lp.lightDir[0]
                 + lp.lightDir[1]*lp.lightDir[1]
                 + lp.lightDir[2]*lp.lightDir[2]);
   if (mag > 0)
   {
     lp.lightDir[0] /= mag;
     lp.lightDir[1] /= mag;
     lp.lightDir[2] /= mag;
   }

   return lp;
}

double CalculateShading(LightingParameters lp, Triangle tri, int k, Camera c)
{
   //Getting view direction and normalizing it
   double viewDir[3];
   viewDir[0] = c.position[0] - tri.X[k];
   viewDir[1] = c.position[1] - tri.Y[k];
   viewDir[2] = c.position[2] - tri.Z[k];
   double norm = sqrt( (pow(viewDir[0], 2) + pow(viewDir[1], 2) + pow(viewDir[2], 2)) );
   for (int i = 0; i < 3; i++)
      viewDir[i] /= norm;

   //Getting diffuse
   double dot = lp.lightDir[0] * tri.normals[k][0] + lp.lightDir[1] * tri.normals[k][1] +
      lp.lightDir[2] * tri.normals[k][2];
   double diffuse = dot;
   if (diffuse < 0)
   {
      diffuse = 0;
   }

   //Getting Spectural light
   //Getting R
   double R[3];
   double tmp = dot * 2;
   for (int i = 0; i < 3; i++)
   {
      R[i] = tri.normals[k][i];
      R[i] *= tmp;
   }
   for (int i = 0; i < 3; i++)
   {
      R[i] -= lp.lightDir[i];
   }

   //Normalizing R
   norm = sqrt( (pow(R[0], 2) + pow(R[1], 2) + pow(R[2], 2)) );
   for (int i = 0; i < 3; i++)
      R[i] /= norm;

   //Getting cos(alpha)
   double cos_alpha = R[0] * viewDir[0] + R[1] * viewDir[1] + R[2] * viewDir[2];
   if (cos_alpha < 0)
   {
      cos_alpha = 0;
   }

   double specular = fabs(pow(cos_alpha, lp.alpha));

   double shading_amount = lp.Ka + lp.Kd * diffuse + lp.Ks * specular;
   return shading_amount;   
}

double lerpColors(double i, double x1, double x2, double c1, double c2)
{
   double t = (i - x1) / (x2 - x1);
   double z = c1 + t * (c2 - c1);
   return z;
}

void rasterizeUp(Triangle tl, Image *image)
{
   double rowMax = -INFINITY;
   double rowMin = INFINITY;
   int left = 0;
   int right = 0;
   int top;
   //finding the top of triangle
   for (int i = 0; i < 3; i++)
   {
      if (tl.Y[i] > rowMax)
      {
         rowMax = tl.Y[i];
         top = i;
      }
   }
   rowMax = F441(rowMax);

   //finding the left and the right of the triangle and setting the rowMin value
   switch(top) {
      case 0:
      if (tl.X[1] > tl.X[2])
      {
         right = 1;
         left = 2;
         rowMin = tl.Y[1];
      }
      else
      {
         right = 2;
         left = 1;
         rowMin = tl.Y[1];
      }
      break;
      case 1:
      if (tl.X[0] > tl.X[2])
      {
         right = 0;
         left = 2;
         rowMin = tl.Y[0];
      }
      else
      {
         right = 2;
         left = 0;
         rowMin = tl.Y[0];
      }
      break;
      case 2:
      if (tl.X[0] > tl.X[1])
      {
         right = 0;
         left = 1;
         rowMin = tl.Y[1];
      }
      else
      {
         right = 1;
         left = 0;
         rowMin = tl.Y[1];
      }
      break;
   }
   rowMin = C441(rowMin);


   double leftB, rightB, tmpX, tmpY, leftSlope, rightSlope, leftX, rightX;
   int checkLeft = 0;
   int checkRight = 0;
   //calculating the slope and b value for the left side of the triangle
   tmpY = tl.Y[top] - tl.Y[left];
   tmpX = tl.X[top] - tl.X[left];

   if (tmpX != 0)
   {
      leftSlope = tmpY / tmpX;
      leftB = tl.Y[top] - (tl.X[top] * leftSlope);
   }
   else
   {
      leftX = tl.X[top];
      checkLeft = 1;
   }
  

   //calculating the slope and b value for the right side of the triangle
   tmpY = tl.Y[top] - tl.Y[right];
   tmpX = tl.X[top] - tl.X[right];
   if (tmpX != 0)
   {
      rightSlope = tmpY / tmpX;
      rightB = tl.Y[top] - (tl.X[top] * rightSlope);
   }
   else
   {
      rightX = tl.X[top];
      checkRight = 1;
   }

   double leftEnd, rightEnd;
   for (int i = rowMin; i <= rowMax; i++)
   {
      //checking if out of bounds top and bottom
      if (i >= ROW)
      {
         break;
      }
      else if (i < 0)
      {
         continue;
      }

      if (!checkLeft)
      {
         leftEnd = ((i - leftB) / leftSlope);
      }
      else
      {
         leftEnd = leftX;
      }
      double start = C441(leftEnd);

      if (!checkRight)
      {
         rightEnd = ((i - rightB) / rightSlope);
      }
      else
      {
         rightEnd = rightX;
      }
      double end = F441(rightEnd);

      //lerping to get the z value of two intercepts to then find the z value of 
      //specific point along a row
      double t = (i - tl.Y[left]) / (tl.Y[top] - tl.Y[left]);
      double zL = tl.Z[left] + t * (tl.Z[top] - tl.Z[left]);

      double q = (i - tl.Y[right]) / (tl.Y[top] - tl.Y[right]);
      double zR = tl.Z[right] + q * (tl.Z[top] - tl.Z[right]);

      //Lerping colors from row to row
      double cLR = lerpColors(i, tl.Y[top], tl.Y[left], tl.color[top][0], tl.color[left][0]);
      double cLG = lerpColors(i, tl.Y[top], tl.Y[left], tl.color[top][1], tl.color[left][1]);
      double cLB = lerpColors(i, tl.Y[top], tl.Y[left], tl.color[top][2], tl.color[left][2]);

      double cRR = lerpColors(i, tl.Y[top], tl.Y[right], tl.color[top][0], tl.color[right][0]);
      double cRG = lerpColors(i, tl.Y[top], tl.Y[right], tl.color[top][1], tl.color[right][1]);
      double cRB = lerpColors(i, tl.Y[top], tl.Y[right], tl.color[top][2], tl.color[right][2]);

      double sL = tl.shading[left] + t * (tl.shading[top] - tl.shading[left]);
      double sR = tl.shading[right] + q * (tl.shading[top] - tl.shading[right]);
      
      for (int k = start; k <= end; k++)
      {
         if (k >= COL)
         {
            break;
         }
         else if (k < 0)
         {
            continue;
         }
         
         //Lerping to get z value of the current col
         double p = (k - leftEnd) / (rightEnd - leftEnd);
         double z = zL + p * (zR - zL);
         double s = sL + p * (sR - sL);

         //Lerping colors from col to col
         double c1 = lerpColors(k, leftEnd, rightEnd, cLR, cRR);
         double c2 = lerpColors(k, leftEnd, rightEnd, cLG, cRG);
         double c3 = lerpColors(k, leftEnd, rightEnd, cLB, cRB);

         c1 *= s;
         c2 *= s;
         c3 *= s;
         if (c1 < 0)
         {
            c1 = 0;
         }
         else if (c1 > 1)
         {
            c1 = 1;
         }
         if (c2 < 0)
         {
            c2 = 0;
         }
         else if (c2 > 1)
         {
            c2 = 1;
         }
         if (c3 < 0)
         {
            c3 = 0;
         }
         else if (c3 > 1)
         {
            c3 = 1;
         }

         //sanity check
         if (z > 0 || z < -1)
            abort();

         //If the z is greater than what is already there, overwrite it
         if (z > image->depth_array[i + COL * k])
         {
            image->array[i + COL * k].one = C441(c1 * 255);
            image->array[i + COL * k].two = C441(c2 * 255);
            image->array[i + COL * k].three = C441(c3 * 255);
            //reset depth array
            image->depth_array[i + COL * k] = z;
         }
      }
   }
}

void rasterizeDown(Triangle tl, Image *image)
{
   double rowMax = -INFINITY;
   double rowMin = INFINITY;
   int left = 0;
   int right = 0;
   int bot;
   //finding the bottom of triangle
   for (int i = 0; i < 3; i++)
   {
      if (tl.Y[i] < rowMin)
      {
         rowMin = tl.Y[i];
         bot = i;
      }
   }
   rowMin = C441(rowMin);

   //finding the left and the right of the triangle and setting the rowMax value
   switch(bot) {
      case 0:
      if (tl.X[1] > tl.X[2])
      {
         right = 1;
         left = 2;
         rowMax = tl.Y[1];
      }
      else
      {
         right = 2;
         left = 1;
         rowMax = tl.Y[1];
      }
      break;
      case 1:
      if (tl.X[0] > tl.X[2])
      {
         right = 0;
         left = 2;
         rowMax = tl.Y[0];
      }
      else
      {
         right = 2;
         left = 0;
         rowMax = tl.Y[0];
      }
      break;
      case 2:
      if (tl.X[0] > tl.X[1])
      {
         right = 0;
         left = 1;
         rowMax = tl.Y[1];
      }
      else
      {
         right = 1;
         left = 0;
         rowMax = tl.Y[1];
      }
      break;
   }
   rowMax = F441(rowMax);


   double leftB, rightB, tmpX, tmpY, leftX, rightX, leftSlope, rightSlope;
   int checkLeft = 0;
   int checkRight = 0;
   //calculating the slope and b value for the left side of the triangle
   tmpY = tl.Y[bot] - tl.Y[left];
   tmpX = tl.X[bot] - tl.X[left];
   if (tmpX != 0)
   {
      leftSlope = tmpY / tmpX;
      leftB = tl.Y[bot] - (tl.X[bot] * leftSlope);
   }
   else
   {
      leftX = tl.X[bot];
      checkLeft = 1;
   }
     

   //calculating the slope and b value for the right side of the triangle
   tmpY = tl.Y[bot] - tl.Y[right];
   tmpX = tl.X[bot] - tl.X[right];
   if (tmpX != 0)
   {
      rightSlope = tmpY / tmpX;
      rightB = tl.Y[bot] - (tl.X[bot] * rightSlope);
   }
   else
   {
      rightX = tl.X[bot];
      checkRight = 1;
   }
   

   double leftEnd, rightEnd;
   for (int i = rowMin; i <= rowMax; i++)
   {
      //checking if out of bounds on the top and bottom
      if (i >= ROW)
      {
         break;
      }
      else if (i < 0)
      {
         continue;
      }

      //checking if it is a left right triangle
      if (!checkLeft)
      {
         leftEnd = ((i - leftB) / leftSlope);
      }
      else
      {
         leftEnd = leftX;
      }
      double start = C441(leftEnd);

      //checking if it is a right right triagnle 
      if (!checkRight)
      {
         rightEnd = ((i - rightB) / rightSlope);
      }
      else
      {
         rightEnd = rightX;
      }
      double end = F441(rightEnd);

      //Lerping z coordinate of two intercepts to then find any z value along a row
      double t = (i - tl.Y[left]) / (tl.Y[bot] - tl.Y[left]);
      double zL = tl.Z[left] + t * (tl.Z[bot] - tl.Z[left]);

      double q = (i - tl.Y[right]) / (tl.Y[bot] - tl.Y[right]);
      double zR = tl.Z[right] + q * (tl.Z[bot] - tl.Z[right]);

      //Lerping colors of two intercepts to then find the color of a any point along a col
      double cLR = lerpColors(i, tl.Y[bot], tl.Y[left], tl.color[bot][0], tl.color[left][0]);
      double cLG = lerpColors(i, tl.Y[bot], tl.Y[left], tl.color[bot][1], tl.color[left][1]);
      double cLB = lerpColors(i, tl.Y[bot], tl.Y[left], tl.color[bot][2], tl.color[left][2]);

      double cRR = lerpColors(i, tl.Y[bot], tl.Y[right], tl.color[bot][0], tl.color[right][0]);
      double cRG = lerpColors(i, tl.Y[bot], tl.Y[right], tl.color[bot][1], tl.color[right][1]);
      double cRB = lerpColors(i, tl.Y[bot], tl.Y[right], tl.color[bot][2], tl.color[right][2]);
      
      double sL = tl.shading[left] + t * (tl.shading[bot] - tl.shading[left]);
      double sR = tl.shading[right] + q * (tl.shading[bot] - tl.shading[right]);

      for (int k = start; k <= end; k++)
      {
         if (k >= COL)
         {
            break;
         }
         else if (k < 0)
         {
            continue;
         }

         //find z value of the given col
         double p = (k - leftEnd) / (rightEnd - leftEnd);
         double z = zL + p * (zR - zL);
         double s = sL + p * (sR - sL);

         //finding the color of the given point 
         double c1 = lerpColors(k, leftEnd, rightEnd, cLR, cRR);
         double c2 = lerpColors(k, leftEnd, rightEnd, cLG, cRG);
         double c3 = lerpColors(k, leftEnd, rightEnd, cLB, cRB);

         c1 *= s;
         c2 *= s;
         c3 *= s;

         if (c1 < 0)
         {
            c1 = 0;
         }
         else if (c1 > 1)
         {
            c1 = 1;
         }
         if (c2 < 0)
         {
            c2 = 0;
         }
         else if (c2 > 1)
         {
            c2 = 1;
         }
         if (c3 < 0)
         {
            c3 = 0;
         }
         else if (c3 > 1)
         {
            c3 = 1;
         }

         //sanity check
         if (z > 0 || z < -1)
            abort();
         

         if (z > image->depth_array[i + COL * k])
         {
            image->array[i + COL * k].one = C441(c1 * 255);
            image->array[i + COL * k].two = C441(c2 * 255);
            image->array[i + COL * k].three = C441(c3 * 255);
            //resetting depth array
            image->depth_array[i + COL * k] = z;
         }
      }
   }
}

void rasterizeTriangle(Triangle tl, Image *image)
{
   double rowMax = -INFINITY;
   double rowMin = INFINITY;
   double pointY;
   int pointI;
   int check = 0;
   int top = -1;
   int bot = -1;

   //checking to see if it is not abritrary 
   if (tl.Y[0] == tl.Y[1])
   {
      if (tl.Y[0] < tl.Y[2])
      {
         rowMax = tl.Y[2];
         top = 2;
      }
      else
      {
         rowMin = tl.Y[2];
         bot = 2;
      }
      check = 1;
   }
   else if (tl.Y[0] == tl.Y[2])
   {
      if (tl.Y[0] < tl.Y[1])
      {
         rowMax = tl.Y[1];
         top = 1;
      }
      else
      {
         rowMin = tl.Y[1];
         bot = 1;
      }
      check = 1;
   }
   else if (tl.Y[1] == tl.Y[2])
   {
      if (tl.Y[1] < tl.Y[0])
      {
         rowMax = tl.Y[0];
         top = 0;
      }
      else
      {
         rowMin = tl.Y[0];
         bot = 0;
      }
      check = 1;
   }

   if (check == 1)
   {
      if (top >= 0)
      {
         rasterizeUp(tl, image);
      }
      else
      {
         rasterizeDown(tl, image);
      }
   }
   //dealing with arbitrary triangles
   else
   {
      //finding top and bottom
      for (int i = 0; i < 3; i++)
      {
         if (tl.Y[i] > rowMax)
         {
            rowMax = tl.Y[i];
            top = i;
         }
         if (tl.Y[i] < rowMin)
         {
            rowMin = tl.Y[i];
            bot = i;
         }
      }
      //finding the other point other than top and bottom
      for (int i = 0; i < 3; i++)
      {
         if (i != top && i != bot)
         {
            pointY = tl.Y[i];
            pointI = i;
         }
      }
      //Getting new point to split triangle up
      double b, tmpX, tmpY, slope, newPointX;
      double newPointY = tl.Y[pointI];

      int checkZero = 0;
      tmpY = tl.Y[top] - tl.Y[bot];
      tmpX = tl.X[top] - tl.X[bot];
      if (tmpX == 0)
      {
         checkZero = 1;
      }
      
      if (!checkZero)
      {
         slope = tmpY / tmpX;
         b = tl.Y[top] - (tl.X[top] * slope);

         newPointX = (pointY - b) / slope;
      }
      else
      {
         newPointX = tl.X[top];
      }
      
      //Creating a temporary trianlge that will be either a regular up or down triangle, that will
      //be passed to my rastorize up and down functions  
      Triangle triTmp;

      //lerping for z value of new point
      double t = (newPointY - tl.Y[bot]) / (tl.Y[top] - tl.Y[bot]);
      double newPointZ = tl.Z[bot] + t * (tl.Z[top] - tl.Z[bot]);
      double newPointS = tl.shading[bot] + t * (tl.shading[top] - tl.shading[bot]);
      
      //lerping to get color value of new vertex
      double cR = lerpColors(newPointY, tl.Y[top], tl.Y[bot], tl.color[top][0], tl.color[bot][0]);
      double cG = lerpColors(newPointY, tl.Y[top], tl.Y[bot], tl.color[top][1], tl.color[bot][1]);
      double cB = lerpColors(newPointY, tl.Y[top], tl.Y[bot], tl.color[top][2], tl.color[bot][2]);
      
      /*
      *   Initializing temporary triangle with all the values given from the triangle list and from
      *   the new vertex
      */
      //initializing the new point vertex
      triTmp.X[0] = newPointX;
      triTmp.Y[0] = newPointY;
      triTmp.Z[0] = newPointZ;
      triTmp.color[0][0] = cR;
      triTmp.color[0][1] = cG;
      triTmp.color[0][2] = cB;
      triTmp.shading[0] = newPointS;

      
      //initializing the middle vertex
      triTmp.X[1] = tl.X[pointI];
      triTmp.Y[1] = tl.Y[pointI];
      triTmp.Z[1] = tl.Z[pointI];
      triTmp.shading[1] = tl.shading[pointI];

      for (int i = 0; i < 3; i++)
      {
         triTmp.color[1][i] = tl.color[pointI][i];
      }

      //initializing top point vertex
      triTmp.X[2] = tl.X[top];
      triTmp.Y[2] = tl.Y[top];
      triTmp.Z[2] = tl.Z[top];
      triTmp.shading[2] = tl.shading[top];
      for (int i = 0; i < 3; i++)
      {
         triTmp.color[2][i] = tl.color[top][i];
      }
      

      //now that the temporary triangle is an up sided regualr triangle, I give it to my rastorizeUp
      //function 
      rasterizeUp(triTmp, image);
      //need to get the bottom triangle, so I change the top vertex to the bottom vertex because
      //all of the other points are still in place
      triTmp.X[2] = tl.X[bot];
      triTmp.Y[2] = tl.Y[bot];
      triTmp.Z[2] = tl.Z[bot];
      triTmp.shading[2] = tl.shading[bot];
      for (int i = 0; i < 3; i++)
      {
         triTmp.color[2][i] = tl.color[bot][i];
      }
      //now that the temporary triangle is a regular down triangle, I hand it off to my rastorizeDown
      //function
      rasterizeDown(triTmp, image);
   }
}

typedef struct
{
   double A[4][4];     // A[i][j] means row i, column j
} Matrix;


void PrintMatrix(Matrix m)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        printf("(%.7f %.7f %.7f %.7f)\n", m.A[i][0], m.A[i][1], m.A[i][2], m.A[i][3]);
    }
}

Matrix ComposeMatrices(Matrix M1, Matrix M2)
{
    Matrix m_out;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            m_out.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                m_out.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }
    return m_out;
}

void TransformPoint(Matrix m, const double *ptIn, double *ptOut)
{  
    ptOut[0] = ptIn[0]*m.A[0][0]
             + ptIn[1]*m.A[1][0]
             + ptIn[2]*m.A[2][0]
             + ptIn[3]*m.A[3][0];
    ptOut[1] = ptIn[0]*m.A[0][1]
             + ptIn[1]*m.A[1][1]
             + ptIn[2]*m.A[2][1]
             + ptIn[3]*m.A[3][1];
    ptOut[2] = ptIn[0]*m.A[0][2]
             + ptIn[1]*m.A[1][2]
             + ptIn[2]*m.A[2][2]
             + ptIn[3]*m.A[3][2];
    ptOut[3] = ptIn[0]*m.A[0][3]
             + ptIn[1]*m.A[1][3]
             + ptIn[2]*m.A[2][3]
             + ptIn[3]*m.A[3][3];
}

double SineParameterize(int curFrame, int nFrames, int ramp)
{  
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {        
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }        
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
} 

Camera GetCamera(int frame, int nframes)
{            
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0; 
    c.focus[1] = 0; 
    c.focus[2] = 0;
    c.up[0] = 0;    
    c.up[1] = 1;    
    c.up[2] = 0;    
    return c;       
}

Matrix GetViewTransform(Camera c)
{
   Matrix m;

   for (int i = 0; i < 4; i++)
   {
      for (int j = 0; j < 4; j++)
      {
         m.A[i][j] = 0;
      }
   }
   
   m.A[0][0] = 1 / tan((c.angle / 2));
   m.A[1][1] = 1 / tan((c.angle / 2));
   m.A[2][2] = ((c.far + c.near) / (c.far - c.near));
   m.A[2][3] = -1;
   m.A[3][2] = 2 * c.far * c.near / (c.far - c.near);

   return m;
}

Matrix GetCameraTransform(Camera c)
{   
   Matrix rv;
   double o[3], w[3], u[3], v[3];

   for (int i = 0; i < 3; i++)
      o[i] = c.position[i];

   for (int i = 0; i < 3; i++)
      w[i] = o[i] - c.focus[i];

   //normalizing w
   double norm = sqrt( (pow(w[0], 2) + pow(w[1], 2) + pow(w[2], 2)) );
   for (int i = 0; i < 3; i++)
      w[i] /= norm;

   //u = Up X w'
   u[0] = ((c.up[1] * w[2]) - (c.up[2] * w[1]));
   u[1] = ((w[0] * c.up[2]) - (c.up[0] * w[2]));
   u[2] = ((c.up[0] * w[1]) - (c.up[1] * w[0]));
   //normalizing u
   norm = sqrt( (pow(u[0], 2) + pow(u[1], 2) + pow(u[2], 2)) );
   for (int i = 0; i < 3; i++)
      u[i] /= norm;

   //v = w' X u'
   v[0] = ((w[1] * u[2]) - (w[2] * u[1]));
   v[1] = ((u[0] * w[2]) - (w[0] * u[2]));
   v[2] = ((w[0] * u[1]) - (w[1] * u[0]));
   //normalizing v
   norm = sqrt( (pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2)) );
   for (int i = 0; i < 3; i++)
      v[i] /= norm;

   /* u */
   rv.A[0][0] = u[0];
   rv.A[1][0] = u[1];
   rv.A[2][0] = u[2];
   /* v */
   rv.A[0][1] = v[0];
   rv.A[1][1] = v[1];
   rv.A[2][1] = v[2];
   /* w */
   rv.A[0][2] = w[0];
   rv.A[1][2] = w[1];
   rv.A[2][2] = w[2];
   /* ints */
   rv.A[0][3] = 0;
   rv.A[1][3] = 0;
   rv.A[2][3] = 0;
   rv.A[3][3] = 1;

   double t[3];
   for (int i = 0; i < 3; i++)
      t[i] = 0 - o[i];

   rv.A[3][0] = u[0] * t[0] + u[1] * t[1] + u[2] * t[2];
   rv.A[3][1] = v[0] * t[0] + v[1] * t[1] + v[2] * t[2];
   rv.A[3][2] = w[0] * t[0] + w[1] * t[1] + w[2] * t[2];
   
   return rv;
}

Matrix GetDeviceTransform()
{   
   Matrix rv;

   for (int i = 0; i < 4; i++)
   {
      for (int j = 0; j < 4; j++)
      {
         rv.A[i][j] = 0;
      }
   }

   rv.A[0][0] = ROW / 2;
   rv.A[1][1] = COL / 2;
   rv.A[2][2] = 1;
   rv.A[3][0] = ROW / 2;
   rv.A[3][1] = COL / 2;
   rv.A[3][3] = 1;

   return rv;
}

void InitializeScreen(Image *image)
{
   //initializing depth array
   for (int x = 0; x < ROW; x++)
   {
      for (int y = 0; y < COL; y++)
      {
         image->depth_array[x + COL * y] = -1;
      }
   }
   
   //make background black
   for (int x = 0; x < ROW; x++)
   {
      for (int y = 0; y < COL; y++)
      {
         image->array[x + COL * y].one = 0;
         image->array[x + COL * y].two = 0;
         image->array[x + COL * y].three = 0;
      }
   }
}

void TransformAndRenderTriangles(Camera c, TriangleList *tl, Image *image)
{
   Matrix camera_matrix = GetCameraTransform(c);
   Matrix view_matrix = GetViewTransform(c);
   Matrix device_matrix = GetDeviceTransform(c);

   Matrix result = ComposeMatrices(camera_matrix, view_matrix);
   Matrix final_result = ComposeMatrices(result, device_matrix);

   double list_in[4];
   double list_out[4];
   
   for (int i = 0; i < tl->numTriangles; i++)
   {
      Triangle triTmp = tl->triangles[i];
      LightingParameters lp = GetLighting(c);
      triTmp.shading[0] = CalculateShading(lp, triTmp, 0, c);
      triTmp.shading[1] = CalculateShading(lp, triTmp, 1, c);
      triTmp.shading[2]= CalculateShading(lp, triTmp, 2, c);
      for (int j = 0; j < 3; j++)
      {  
         list_in[0] = tl->triangles[i].X[j];
         list_in[1] = tl->triangles[i].Y[j];
         list_in[2] = tl->triangles[i].Z[j];
         list_in[3] = 1;
         TransformPoint(final_result, list_in, list_out);
         triTmp.X[j] = list_out[0] / list_out[3];
         triTmp.Y[j] = list_out[1] / list_out[3];
         triTmp.Z[j] = list_out[2] / list_out[3];
      }
      rasterizeTriangle(triTmp, image);
   }
}


int main(UNUSED int argc, UNUSED char const *argv[])
{
   //creating an image that will hold a 1D array that will be index like a 2D array to make the image
   //as well as a depth array to hold z values to know what is in front and back 
   Image *image = (Image *)malloc(sizeof(Image));
   image->array = (Pixel *)malloc(sizeof(Pixel) * MAX);
   image->depth_array = (double *)malloc(sizeof(double) * MAX);

   //Get triangle list
   TriangleList *tl = Get3DTriangles();

   for (int i = 0; i < 1000; i ++)
   {
      InitializeScreen(image);
      Camera c = GetCamera(i, 1000);
      TransformAndRenderTriangles(c, tl, image);
      char file_name[50];
      sprintf(file_name, "proj1F_frame%04d.pnm", i);
      convert(image, file_name);
   }

   free(image->array);
   free(image->depth_array);
   free(image);

   return 0;
}
