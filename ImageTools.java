import java.awt.*;
import java.awt.Image.*;
import javax.swing.ImageIcon;
import java.awt.image.MemoryImageSource;
import java.awt.image.PixelGrabber;
import java.util.Random;
import java.io.*;
public class ImageTools{
  //**************************************************************************
      /**
     * Computes the <b>partition coefficient</b> as proposed by Bezdek in
     * <i>"Bezdek JC. Cluster validity with fuzzy sets. J Cybern 1974;3:58–73"</i>
     * It considers the average sum of the squared fuzzy membership values
     * Better clustering partitions induce higher values of the coefficient
     * @param U The float matrix containing the fuzzy memberships, result of the fuzzy
     * clustering algorithm
     * @return the computed partition coefficient as a double
     */
  public static double mse(int [][][] input,double [][] centroids, int current_centroid, int[][] index_list){
    double total = 0;
    int count = 0;
    //Iterate through every pixel
    for(int row = 0; row < input[0].length; row++){
      for(int col = 0; col < input[0][0].length; col++) {
        int v = 0;
        if (index_list[row][col] == current_centroid) {
          for (int component = 1; component < input.length; component++) {
            v += Math.pow((input[component][row][col] - centroids[current_centroid][component]), 2);
          }
          count++;
        }
        total += v;
      }
    }
    //System.out.println("Total pixels for centroid"+current_centroid+"is:\t"+count);
    return ((Math.sqrt(total/count)));
  }

  public static double dbIndex(int [][][] input,double[][] kCenters, int[][] index_list){
    double [][] rij = new double[kCenters.length][kCenters.length];       //store all the r vals to be summed
    double [] r_max = new double[kCenters.length];       //to store all r max vals
    double [] s = new double[kCenters.length];
    double var;
    double result = 0;
    double r_val;
    double db_index_val = 0;
    int size;

    // CALCULATING THE Si VALUES
    for (int i = 0; i<kCenters.length; i++){
      for (int j = 0; j<kCenters.length; j++){
        rij[i][j] = 0;
      }
      s[i] = mse(input,kCenters,i,index_list);
      //System.out.println("This is the S["+i+"] value:"+ s[i]);
    }

    // CALCULATING Rij VALUES
    for (int i = 0; i<kCenters.length;i++){
      for(int j = 0; j<kCenters.length;j++){
        double d = 0;
        if(i!=j){
          for(int component = 1; component<input.length; component++) {
            d += Math.pow(kCenters[i][component] - kCenters[j][component], 2);
          }
        }
        d = Math.sqrt(d);
        rij[i][j] = ((s[i]+s[j])/d);
        //System.out.println("This is the R["+i+"]["+j+"] values:"+rij[i][j]);
      }
    }

    // FIND THE MAX Ri VALUE
    for (int i = 0; i<kCenters.length;i++){
      double max = 0;
      for (int j = 0; j<kCenters.length;j++){
        if(i!=j && max<rij[i][j]){
          max = rij[i][j];
        }
      }
      r_max[i] = max;
    }

    for (int i = 0; i<kCenters.length; i++) {
      db_index_val += r_max[i];
    }
    db_index_val /= kCenters.length;

    return db_index_val;
  }

  public static double bhattacharyya(int [][][] input, double[][] kCenters, int[][] index_list){
    double b_distance = 0;
    int[][][] cluster;


    return b_distance;
  }

  public static double[][] covariance(int [][][] input,double [][] centroids, int center1, int center2, int[][] index_list){
    int elementc1 = 0, elementc2 = 0;
    int c1index = 0 , c2index = 0;

    for(int row = 0; row<input[0].length; row++){
      for(int col = 0; col<input[0][0].length; col++){
        if(index_list[row][col] == center1){
          elementc1++;
        }
        else if(index_list[row][col] == center2){
          elementc2++;
        }
      }
    }

    int[][] cluster1 = new int[centroids[0].length][elementc1];
    int[][] cluster2 = new int[centroids[0].length][elementc2];

    double cov_matrix[][] = new double [elementc1][elementc2];

    for(int row = 0; row<input[0].length; row++) {
      for (int col = 0; col < input[0][0].length; col++) {
        if(index_list[row][col] == center1){
          for(int i = 0; i<3; i++){
            cluster1[i][c1index] = input[i][row][col];
            c1index++;
          }
        }
        else if(index_list[row][col] == center2){
          for(int i = 3; i<3; i++){
            cluster2[i][c2index] = input[i][row][col];
            c2index++;
          }
        }
      }
    }

    for(int c1 = 0; c1)

    return cov_matrix;
  }


  public static double iIndex(int[][][] input, double[][][] membership, double[][] kCenters, double fuzziness){

    double eK = 0;
    double eOne = 0;
    double max=Double.MIN_VALUE;

    for(int c = 0; c < kCenters.length; c++){
         for(int row=0;row<membership[0].length;row++){
           for(int column = 0; column < membership.length; column++)
           {

             double dist = 0;
             for(int d=1;d<kCenters[0].length;d++){
               dist+=(Math.pow(input[d][column][row]-kCenters[c][d],2));
             }
             dist=Math.sqrt(dist);
             eK += membership[column][row][c] * dist;
             if(c==0){
              eOne+= membership[column][row][c] * dist;
             }
           }

         }
       }

    for(int clus = 0; clus < kCenters.length; clus++){
      for(int clus2 = 0; clus2 < kCenters.length; clus2++){
        if(clus!=clus2){
          double dist = 0;
          for(int d=1;d<kCenters[0].length;d++){
            dist+=(Math.pow(kCenters[clus][d]-kCenters[clus2][d],2));
          }
          dist=Math.sqrt(dist);
          if(dist>max){
            max=dist;
          }
        }
      }
    }
    double res = Math.pow((1.0/kCenters.length)*(eOne/eK)*max,fuzziness);
    return res;


  }


    public static double compactnessAndSeparationMetric(int[][][] input, double[][][] membership,
                                                        double[][] kCenters, double fuzziness){

       double numerator = 0, denominator = 0, min = Double.MAX_VALUE;

       for(int c = 0; c < kCenters.length; c++){
         for(int row=0;row<membership[0].length;row++){
           for(int column = 0; column < membership.length; column++)
           {

             double dist = 0;
             for(int d=1;d<kCenters[0].length;d++){
               dist+=(Math.pow(input[d][column][row]-kCenters[c][d],2));
             }
             numerator += Math.pow(membership[column][row][c], fuzziness) * dist;
           }

         }
       }
       for(int clus = 0; clus < kCenters.length; clus++){
         for(int clus2 = 0; clus2 < kCenters.length; clus2++){
           if(clus!=clus2){
             double dist = 0;
             for(int d=1;d<kCenters[0].length;d++){
               dist+=(Math.pow(kCenters[clus][d]-kCenters[clus2][d],2));
             }
             dist=Math.sqrt(dist);

             if(dist<min){
               min=dist;
             }
           }
         }
       }
       denominator = min*membership.length*membership[0].length;
       double res = 1.0 * numerator/denominator;

     return res;
    }

    public static double[][][] getNormalize(double[][][] memb){
      double[][][] poss=memb;
      int clus_index=-1;
      for(int i=0;i<memb.length;i++){
        for(int j=0;j<memb[i].length;j++){
          double max=0;
          for(int k=0;k<memb[i][j].length;k++){//get max value in cluster set or sup
            if(memb[i][j][k]>max){
              max=memb[i][j][k];
              clus_index=k;
            }
          }
          for(int k=0;k<memb[i][j].length;k++){//update max value in cluster set
            poss[i][j][k]=memb[i][j][k]/memb[i][j][clus_index];
          }
        }
      }

      return poss;



    }

//**************************************************************************
  protected static void writeToFile(double membership[][][],String file){

  try{ BufferedWriter writer = new BufferedWriter(new FileWriter(file+".csv"));
    for (int row = 0; row < membership.length; row++){
      for (int column = 0; column < membership[row].length; column++)
      {
        writer.write(row+","+column+",");
        for(int i=0;i<membership[row][column].length;i++){
         // System.out.println(membership[column].length);
          writer.write(membership[row][column][i]+",");
        }
        writer.write("\n");
      }

    }
    writer.close();
  }catch(Exception e){
    System.out.println(e);
  }
}

  protected static double getEuclideanDist(double[] a,int[] b){
    double sum=0;
    for(int i=0;i<a.length;i++){
      sum+=Math.pow(a[i]-b[i],2);

    }

   sum=Math.sqrt(sum);
   return sum;

  }
  protected static int[][] getClustersColor(int cluster){//get Cluster Color to differentiate
    Random rnd=new Random();
    int[][] arr=new int[cluster][3];
    for(int i=0;i<arr.length;i++){
      for(int j=0;j<arr[0].length;j++){
        arr[i][j]=rnd.nextInt(255);
      }

    }
    return arr;
  }

  protected static boolean compareArray(double[][][]arr1,double[][][]arr2,double term){
  for(int i=0;i<arr1.length;i++){
    for(int j=0;j<arr1[i].length;j++){
      for(int k=0;k<arr1[i][j].length;k++){
        if((arr1[i][j][k]-arr2[i][j][k])>term){
          return false;
        }
      }
    }
  }
  return true;



}
    //**************************************************************************
//**************************************************************************
  protected static Dimension getImageDimension(Image imageIn)
  {

    // Convert image to an Image Icon to get the dimensions
    ImageIcon imageIconIn = new ImageIcon(imageIn);
    return new Dimension(imageIconIn.getIconWidth(), imageIconIn.getIconHeight());

  }


  protected static Image pixelsArrayToImage(int pixels[], Dimension imageInDimension)
  {

    // Generate an image object from the pixels array
    return Toolkit.getDefaultToolkit().createImage(new MemoryImageSource((int) imageInDimension.getWidth(), (int) imageInDimension.getHeight(), pixels, 0, (int) imageInDimension.getWidth()));

  }

  protected static int[] imageToPixelsArray(Image imageIn)
  {

    // Declare local storage
    Dimension imageInDimension = getImageDimension(imageIn);
    int imagePixelLength = (int) (imageInDimension.getWidth() * imageInDimension.getHeight());
    int pixels[] = new int[imagePixelLength];

    // Convert image to array via PixelGrabber
    PixelGrabber pixelGrabber = new PixelGrabber(imageIn, 0, 0, (int) imageInDimension.getWidth(), (int) imageInDimension.getHeight(), pixels, 0, (int) imageInDimension.getWidth());
    try
    {
      pixelGrabber.grabPixels();
    }
    catch(InterruptedException ie)
    {
      System.exit(1);
    }

    // Return the pixels array
    return pixels;

  }

  protected static int[] TRGBArrayToPixelsArray(int[][][] TRGB, Dimension imageInDimension)
  {

    // Declare local storage
    int imagePixelLength = (int) (imageInDimension.getWidth() * imageInDimension.getHeight());
    int pixels[] = new int[imagePixelLength];

    // Loop through TRGB array, assembling the TRGB values
    for (int column = 0, row = 0, pixelIndex = 0; pixelIndex < imagePixelLength; pixelIndex++)
    {

      // Retreive TRGB
      pixels[pixelIndex] = getTRGB(TRGB[0][column][row], TRGB[1][column][row], TRGB[2][column][row], TRGB[3][column][row]);

      // Calculate column and row indexes
      if (++column == imageInDimension.getWidth())
      {
        // Reset column and increment row
        column = 0;
        row++;
      }
    }

    // Return the pixels array
    return pixels;

  }

  protected static int[][][] pixelsArrayToTRGBArray(int[] pixels, Dimension imageInDimension)
  {

    // Declare local storage
    int imagePixelLength = (int) (imageInDimension.getWidth() * imageInDimension.getHeight());
    int TRGB[][][] = new int[4][(int) imageInDimension.getWidth()][(int) imageInDimension.getHeight()];

    // Convert pixel array to TRGB array
    for (int column = 0, row = 0, pixelIndex = 0; pixelIndex < imagePixelLength; pixelIndex++)
    {

      // Store transparency
      TRGB[0][column][row] = getTransparencyComponent(pixels[pixelIndex]);

      // Store red
      TRGB[1][column][row] = getRedComponent(pixels[pixelIndex]);

      // Store green
      TRGB[2][column][row] = getGreenComponent(pixels[pixelIndex]);

      // Store blue
      TRGB[3][column][row] = getBlueComponent(pixels[pixelIndex]);

      // Calculate column and row indexes
      if (++column == imageInDimension.getWidth())
      {
        // Reset column and increment row
        column = 0;
        row++;
      }

    }
    // Return the newly generated TRGB array
    return TRGB;

  }

  // TRGB Component retrieval functions //////////////////////////////
  protected static final int getTransparencyComponent(int pixel)
  {
    return (pixel >> 24) & 0xff; // Return Transparency
  }

  protected static final int getRedComponent(int pixel)
  {
    return (pixel >> 16) & 0xff; // Return Red
  }

  protected static final int getGreenComponent(int pixel)
  {
    return (pixel >> 8) & 0xff; // Return Green
  }

  protected static final int getBlueComponent(int pixel)
  {
    return pixel & 0xff; // Return Blue
  }

  protected static final int getTRGB(int transparency, int red, int green, int blue)
  {
    return (transparency << 24) | (red << 16) | (green << 8) | (blue); // Return TRGB
  }

  // Public


  // Inner Class(es) /////////////////////////////////////////////////
  // Private

  // Protected

  // Public

}
