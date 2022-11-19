#include<iostream>
#include<cmath>
#include<vector>
#include<fenv.h>




//===========================================================
/// Spherical coordinate fun
//===========================================================
namespace SphericalCoordinates
{

 /// Method for calculation of distances
 std::string Method="Exact";

 /// Mean n-s angle, theta, for simplified calculation
 double Mean_theta_deg_for_simple_method=50.0;
 
 /// Radius (Earth in km)
 double Radius=6378.388;

 /// Pi
 double Pi=4.0*atan(1.0);
 

 
 //===========================================================
/// Coordinate class
 //===========================================================
 class Coordinate
 {

 public:
  
  /// Constructor:
  Coordinate() : Theta(0.0), Phi(0.0)
   {}
  
  
  /// Constructor: (theta, phi) = longitude, latitude
  Coordinate(const double& theta_deg, const double& phi_deg) :
   Theta(theta_deg*Pi/180.0), Phi(phi_deg*Pi/180.0)
   {}

  /// phi (around the equator; -180 180)
  double phi_deg() const
   {
    return Phi*180.0/Pi;
   }

  
  /// theta (south pole to north pole): north pole: 90; south
  /// pole -90
  double theta_deg() const
   {
    return Theta*180.0/Pi;;
   }

  
  /// phi (around the equator; (-pi,pi))
  double& phi() 
   {
    return Phi;
   }

  
  /// theta (south pole to north pole): north pole: pi/2; south
  /// pole -pi/2)
  double& theta()
   {
    return Theta;
   }
  
  /// phi (around the equator; (-pi,pi))
  double phi()  const
   {
    return Phi;
   }

  
  /// theta (south pole to north pole): north pole: pi/2; south
  /// pole -pi/2)
  double theta() const
   {
    return Theta;
   }
  
  
 private:
  
  /// phi (around the equator; (-pi,pi))
  double Phi;
  
  
  /// theta (south pole to north pole): north pole: pi/2; south
  /// pole -pi/2)
  double Theta;
  
 };


  
  
 /// Compute distance
 double distance(const Coordinate& coord1,
                 const Coordinate& coord2)
 {
   
  if (Method=="Exact")
   {
    // https://www.kompf.de/gps/distcalc.html
    double lat1=coord1.theta();
    double lat2=coord2.theta();
    double lon1=coord1.phi();
    double lon2=coord2.phi();
    double ds=Radius*acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1));
    return ds;
   }
  else if (Method=="Simplest")
   {
    double dtheta=coord2.theta()-coord1.theta();
    double dphi  =coord2.phi()  -coord1.phi();

    // scale factors
    double h_phi=Radius*cos(Pi/180.0*Mean_theta_deg_for_simple_method);
    double h_theta=Radius;

    double ds=sqrt(pow(h_phi*dphi,2)+pow(h_theta*dtheta,2));
    return ds;
   }
  else
   {
    std::cout << "Never get here; Method = " << Method << std::endl;
    abort();
     
    //dummy return to shut up compiler
    return 1.0;
   }
 }
   


//=======================================================
/// hierher
//=======================================================
 class PointFinder
 {

 
 public:

  /// Constructor: Point1, Point2, required total distance
  PointFinder(const Coordinate& coord1,
              const Coordinate& coord2,
              const double& percentage_excess_total_distance) :
   Coord1(coord1), Coord2(coord2)
   {
    // Average angle for simple method
    Mean_theta_deg_for_simple_method=
     0.5*(Coord1.theta_deg()+Coord2.theta_deg());

    // distance between original points
    double direct_distance=distance(coord1,coord2);

    // Desired new length
    Required_total_distance=(1.0+percentage_excess_total_distance/100.0)*
     direct_distance;
            
    if (Required_total_distance<direct_distance)
     {
      std::cout << "Cannot be done without a tunnel!\n";
      abort();
     }
    else
     {
      std::cout << "percentage extra distance: "
                << fabs(Required_total_distance-direct_distance)/
       fabs(direct_distance)*100.0 << std::endl;
     }
    find_new_point();
   }


 private:

 
/// Residuals
  void residuals(const Coordinate& coord, std::vector<double>& res)
   {
    double dist1=distance(Coord1,coord);
    double dist2=distance(Coord2,coord);
    res[0]=dist1+dist2-Required_total_distance;
    res[1]=dist1-dist2;

    // // hierher newton solver test
    // res[0]=coord.theta()-1.0;
    // res[1]=coord.phi()-2.0;
   
   }

  // Jaobian
  void jacobian(const Coordinate& coord, std::vector<std::vector<double> >& jac)
   {
   
    if (Method=="Exact")
     {
     
      double theta=coord.theta();
      double phi=coord.phi();
     
      double theta1=Coord1.theta();
      double phi1=Coord1.phi();
     
      double theta2=Coord2.theta();
      double phi2=Coord2.phi();
     
      double jac00,jac11,jac01,jac10;
     

      // Maple
      jac00 = -Radius * (sin(theta1) * cos(theta) - cos(theta1) * sin(theta) * cos(phi - phi1)) * pow(-pow(sin(theta1) * sin(theta) + cos(theta1) * cos(theta) * cos(phi - phi1), 0.2e1) + 0.1e1, -0.1e1 / 0.2e1) - Radius * (sin(theta2) * cos(theta) - cos(theta2) * sin(theta) * cos(phi - phi2)) * pow(-pow(sin(theta2) * sin(theta) + cos(theta2) * cos(theta) * cos(phi - phi2), 0.2e1) + 0.1e1, -0.1e1 / 0.2e1);
      jac01 = Radius * cos(theta1) * cos(theta) * sin(phi - phi1) * pow(-pow(sin(theta1) * sin(theta) + cos(theta1) * cos(theta) * cos(phi - phi1), 0.2e1) + 0.1e1, -0.1e1 / 0.2e1) + Radius * cos(theta2) * cos(theta) * sin(phi - phi2) * pow(-pow(sin(theta2) * sin(theta) + cos(theta2) * cos(theta) * cos(phi - phi2), 0.2e1) + 0.1e1, -0.1e1 / 0.2e1);
      jac10 = -Radius * (sin(theta1) * cos(theta) - cos(theta1) * sin(theta) * cos(phi - phi1)) * pow(-pow(sin(theta1) * sin(theta) + cos(theta1) * cos(theta) * cos(phi - phi1), 0.2e1) + 0.1e1, -0.1e1 / 0.2e1) + Radius * (sin(theta2) * cos(theta) - cos(theta2) * sin(theta) * cos(phi - phi2)) * pow(-pow(sin(theta2) * sin(theta) + cos(theta2) * cos(theta) * cos(phi - phi2), 0.2e1) + 0.1e1, -0.1e1 / 0.2e1);
      jac11 = Radius * cos(theta1) * cos(theta) * sin(phi - phi1) * pow(-pow(sin(theta1) * sin(theta) + cos(theta1) * cos(theta) * cos(phi - phi1), 0.2e1) + 0.1e1, -0.1e1 / 0.2e1) - Radius * cos(theta2) * cos(theta) * sin(phi - phi2) * pow(-pow(sin(theta2) * sin(theta) + cos(theta2) * cos(theta) * cos(phi - phi2), 0.2e1) + 0.1e1, -0.1e1 / 0.2e1);


      jac[0][0]=jac00;
      jac[1][1]=jac11;
      jac[1][0]=jac10;
      jac[0][1]=jac01;

     }
    else if (Method=="Simplest")
     {
      double theta=coord.theta();
      double phi=coord.phi();
     
      double theta1=Coord1.theta();
      double phi1=Coord1.phi();
     
      double theta2=Coord2.theta();
      double phi2=Coord2.phi();
     
      double jac00,jac11,jac01,jac10;
     
      double cg=Mean_theta_deg_for_simple_method;

      
      jac00 = pow(pow(phi - phi1, 0.2e1) * Radius * Radius * pow(cos(0.1745329252e-1 * cg), 0.2e1) + pow(theta - theta1, 0.2e1) * Radius * Radius, -0.1e1 / 0.2e1) * (theta - theta1) * Radius * Radius + pow(pow(phi - phi2, 0.2e1) * Radius * Radius * pow(cos(0.1745329252e-1 * cg), 0.2e1) + pow(theta - theta2, 0.2e1) * Radius * Radius, -0.1e1 / 0.2e1) * (theta - theta2) * Radius * Radius;


      jac01 = pow(pow(phi - phi1, 0.2e1) * Radius * Radius * pow(cos(0.1745329252e-1 * cg), 0.2e1) + pow(theta - theta1, 0.2e1) * Radius * Radius, -0.1e1 / 0.2e1) * (phi - phi1) * Radius * Radius * pow(cos(0.1745329252e-1 * cg), 0.2e1) + pow(pow(phi - phi2, 0.2e1) * Radius * Radius * pow(cos(0.1745329252e-1 * cg), 0.2e1) + pow(theta - theta2, 0.2e1) * Radius * Radius, -0.1e1 / 0.2e1) * (phi - phi2) * Radius * Radius * pow(cos(0.1745329252e-1 * cg), 0.2e1);


      jac10 = pow(pow(phi - phi1, 0.2e1) * Radius * Radius * pow(cos(0.1745329252e-1 * cg), 0.2e1) + pow(theta - theta1, 0.2e1) * Radius * Radius, -0.1e1 / 0.2e1) * (theta - theta1) * Radius * Radius - pow(pow(phi - phi2, 0.2e1) * Radius * Radius * pow(cos(0.1745329252e-1 * cg), 0.2e1) + pow(theta - theta2, 0.2e1) * Radius * Radius, -0.1e1 / 0.2e1) * (theta - theta2) * Radius * Radius;


      jac11 = pow(pow(phi - phi1, 0.2e1) * Radius * Radius * pow(cos(0.1745329252e-1 * cg), 0.2e1) + pow(theta - theta1, 0.2e1) * Radius * Radius, -0.1e1 / 0.2e1) * (phi - phi1) * Radius * Radius * pow(cos(0.1745329252e-1 * cg), 0.2e1) - pow(pow(phi - phi2, 0.2e1) * Radius * Radius * pow(cos(0.1745329252e-1 * cg), 0.2e1) + pow(theta - theta2, 0.2e1) * Radius * Radius, -0.1e1 / 0.2e1) * (phi - phi2) * Radius * Radius * pow(cos(0.1745329252e-1 * cg), 0.2e1);

      jac[0][0]=jac00;
      jac[1][1]=jac11;
      jac[1][0]=jac10;
      jac[0][1]=jac01;



     }
    else
     {
      std::cout << "Never get here; Method = " << Method << std::endl;
      abort();
     
     }
   
   }
   

  /// Newton solver:
  void find_new_point()
   {
    // Newton tolerance on max. residual relative to orig length
    double fractional_tol=1.0e-4;
   
    
    // Create initial guess based on constant metric
    double phi_new=0.0;
    double theta_new=0.0;
    
    // Scale factors
    double h_phi=Radius*cos(Pi/180.0*Mean_theta_deg_for_simple_method);
    double h_theta=Radius;
    
    // Geometry
    double orig_distance=distance(Coord2,Coord1);
    double alpha=acos(orig_distance/Required_total_distance);

    double d_theta=Coord2.theta()-Coord1.theta();
    double d_phi  =Coord2.phi()  -Coord1.phi();
    double beta=atan2(h_theta*d_theta,h_phi*d_phi);
    double gamma=beta-alpha;
    
    phi_new  =Coord1.phi()  +0.5*Required_total_distance/h_phi  *cos(gamma);
    theta_new=Coord1.theta()+0.5*Required_total_distance/h_theta*sin(gamma);


    // // Initial guess: Halfway in spherical polars
    // theta_new=0.5*(Coord1.theta_deg()+Coord2.theta_deg());
    // phi_new=0.5*(Coord1.phi_deg()+Coord2.phi_deg());

    
    // Initial guess: 
    Coordinate coord=Coordinate(theta_new*180.0/Pi,phi_new*180.0/Pi);
             
    // Actual tolerance
    double tol=fractional_tol*orig_distance;
    
    
    // Now do both methds: approx to get a good initial guess; then exact
    for (unsigned method=0;method<2;method++)
     {
      Method="Simplest";
      if (method==1) Method="Exact";

      std::cout << "Doing method = " << Method << std::endl;

      // Flag
      bool converged=false;
              
      // Get residual
      std::vector<double> res(2,0.0);
      residuals(coord,res);
      double max_res=std::max(fabs(res[0]),fabs(res[1]));

      // hierher
      // std::cout << "Resid: "
      //           << res[0] << " "
      //           << res[1] << " "
      //           << max_res << std::endl;
      
      // Converged?
      if (max_res>tol)
       {
        std::cout << "Need to newton iterate" << std::endl;
        unsigned max_iter=50;
        for (unsigned iter=0;iter<max_iter;iter++)
         {
          
          std::vector<std::vector<double> > jac(2);
          jac[0].resize(2);
          jac[1].resize(2);
          bool use_fd=false; // hierher
          // Get the Jacobian by FDing
          if (use_fd)
           {
            
            double fd_step=1.0e-8;
            std::vector<double> res_pls(2,0.0);
            Coordinate coord_pls(coord);
            for (unsigned j=0;j<2;j++)
             {
              if (j==0)
               {
                coord_pls.theta()+=fd_step;
               }
              else
               {
                coord_pls.phi()+=fd_step;
               }
              residuals(coord_pls,res_pls);
              for (unsigned i=0;i<2;i++)
               {
                jac[i][j]=(res_pls[i]-res[i])/fd_step;
                // hierher
                // std::cout << "jac " << i << " " << j << " "
                //           << jac[i][j] << " RES: "
                //           << res_pls[i] << " "  << res[i]
                //           << std::endl;
               }
             }
           }
          // Analytical Jacobian
          else
           {
            jacobian(coord,jac);
            // hierher
            // std::cout << "JAC: " << jac[0][0] << " " << jac[0][1] << std::endl;
            // std::cout << "JAC: " << jac[1][0] << " " << jac[1][1] << std::endl; 
           }
          
          
          // Hand calculated soln, so we don't need a generic lu solver
          double det=jac[0][0]*jac[1][1]-jac[1][0]*jac[0][1];
          double dtheta=( jac[1][1]*res[0]-jac[0][1]*res[1])/det;
          double dphi  =(-jac[1][0]*res[0]+jac[0][0]*res[1])/det;
          
          coord.theta()-=dtheta;
          coord.phi()-=dphi;
          
          // Recompute residual
          residuals(coord,res);
          double max_res=std::max(fabs(res[0]),fabs(res[1]));

          // hierher
          // std::cout << "Iter: " << iter
          //           << " Coord: " << " "
          //           << coord.theta() << " "
          //           << coord.phi() << " "
          //           << " dCoord: " << " "
          //           << dtheta << " "
          //           << dphi << " "
          //           << " Resid: "
          //           << res[0] << " "
          //           << res[1] << " "
          //           << max_res << std::endl;
          if (max_res<tol)
           {
            converged=true;
            break;
           }
          
         } // end Newton iterations
       }
      // No iteration was required
      else
       {
        std::cout << "No need to newton iterate" << std::endl;
        converged=true;
       }
      
      // Not converged
      if (!converged)
       {
        std::cout << "Newton iteration didn't converge; die!" << std::endl;
        //abort();
       }
      else
       {
        std::cout << "Newton method converged! "
                  << coord.theta_deg() << " "
                  << coord.phi_deg() << " "
                  << std::endl;
        
        // Test:
        double final_dist1=distance(Coord1,coord);
        double final_dist2=distance(Coord2,coord);
        double orig_dist  =distance(Coord1,Coord2);
        
        std::cout << "Orig   distance via two   points: "
                  << orig_dist               << std::endl;
        std::cout << "Actual distance via three points: "
                  << final_dist1+final_dist2 << std::endl;
        std::cout << "Target distance via three points: "
                  <<  Required_total_distance << std::endl;

       }
     } // end of loop over method
    
   }
 
  /// Coords of first point
  Coordinate Coord1;
 
  /// Coords of second point
  Coordinate Coord2;

  /// Required distance
  double Required_total_distance;

 
 };



} // end namespace



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


int main()
{

 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

 using namespace SphericalCoordinates;

 double percentage_increase=15.0;

 std::cout << "\n====================================================\n\n";
 
  
 Coordinate berlin(52.5164,13.3777);
 Coordinate lisbon(38.692668,-9.177944);
 Method="Exact";
 double dist_berlin_lisbon=distance(berlin,lisbon);
 std::cout << "Berlin-Lisbon: " <<  dist_berlin_lisbon << std::endl;

 PointFinder(berlin,lisbon,percentage_increase);
  
 std::cout << "\n\n====================================================\n\n";

//exit(0);

 Coordinate rhm_bahnhof(49.9917,8.41321);
 Coordinate rhm_bruecke(50.0049,8.42182);
 Method="Exact";
 double dist_rhm=distance(rhm_bahnhof,rhm_bruecke);
 std::cout << "Rhm: " << dist_rhm << std::endl;

 PointFinder(rhm_bahnhof,rhm_bruecke,percentage_increase);
 
 std::cout << "\n\n====================================================\n\n";
 
 
 Coordinate bernd1(50.478238579,     8.2900669);
 Coordinate bernd2(50.47829943178574,8.289963);
 Method="Exact";
 double dist_bernd=distance(bernd1,bernd2);
 std::cout << "Bernd: " <<  dist_bernd<< std::endl;

 PointFinder(bernd1,bernd2,percentage_increase);
 

 return 0;

}
