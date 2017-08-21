#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int lane = 1; //0-left lane, 1-middle lane, 2-right lane
double ref_velocity = 0.0; //MPH

typedef struct _SensorFusion {
  double id;
  double x;
  double y;
  double vx;
  double vy;
  double s;
  double d;

}SensorFusion;

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

            int prev_path_size  = previous_path_x.size();


            if(prev_path_size > 0)
            {
              car_s = end_path_s;
            }

            bool decelerate  = false;
            bool lane_change = false;
            double max_velocity = 49.5;

            vector<SensorFusion> cars_left_ahead;
            vector<SensorFusion> cars_center_ahead;
            vector<SensorFusion> cars_right_ahead;

            vector<SensorFusion> cars_left_behind;
            vector<SensorFusion> cars_center_behind;
            vector<SensorFusion> cars_right_behind;

            for(int i = 0; i < sensor_fusion.size(); i++)
            {
              SensorFusion car;

              double id = sensor_fusion[i][0];
              double  x = sensor_fusion[i][1];
              double  y = sensor_fusion[i][2];
              double vx = sensor_fusion[i][3];
              double vy = sensor_fusion[i][4];
              double  s = sensor_fusion[i][5];
              double  d = sensor_fusion[i][6];

              double obj_speed = sqrt(vx*vx + vy*vy);
              double sensor_range = 400.0;

              car.id = id;
              car.x = x;
              car.y = y;
              car.vx = vx;
              car.vy = vy;
              car.s = s;
              car.d = d;

              if (((s - car_s) >= 0.0) && ((s - car_s) < sensor_range))
              {
                if (d < 4.0)
                {
                  cars_left_ahead.push_back(car);
                }
                else if ((d > 4.0) && (d < 8.0))
                {
                  cars_center_ahead.push_back(car);
                }
                else if(d > 8.0)
                {
                  cars_right_ahead.push_back(car);
                }
              }
              else if (((s - car_s) < 0.0) && ((s - car_s) > -sensor_range))
              {
                if (d < 4.0)
                {
                  cars_left_behind.push_back(car);
                }
                else if ((d > 4.0) && (d < 8.0))
                {
                  cars_center_behind.push_back(car);
                }
                else if(d > 8.0)
                {
                  cars_right_behind.push_back(car);
                }
              }

              if((d < (2+4*lane+2)) && (d > (2+4*lane-2)))
              {
                double vx = sensor_fusion[i][3];
                double vy = sensor_fusion[i][4];
                double obj_speed = sqrt(vx*vx + vy*vy);
                double obj_distance = sensor_fusion[i][5];

                double obj_distance_future = obj_distance + (0.02 * obj_speed * prev_path_size);

                if((obj_distance_future > car_s) && ((obj_distance_future - car_s) < 30.0))
                {
                  lane_change = true;
                  max_velocity = obj_speed;
                }
              }
            }

            printf("lane  - %d ", lane);
            printf("speed - %f ", ref_velocity);
            printf("Ahead - ");
            printf("|%3d|%3d|%3d| ", cars_left_ahead.size(), cars_center_ahead.size(), cars_right_ahead.size());
            printf("Behind - ");
            printf("|%3d|%3d|%3d| ", cars_left_behind.size(), cars_center_behind.size(), cars_right_behind.size());
            printf("Max vel - %5.2f", max_velocity);
            printf("\n");

            if(ref_velocity < max_velocity)
            {
              ref_velocity += (0.224 * 2);
            }
            else
            {
              ref_velocity -= (0.224 * 2);
            }

            if(lane_change == true) {

              int new_lane;
              double min_dist_ahead;
              double max_dist_behind;
              int car_behind;
              int car_ahead;
              double safe_dist = 20.0;

              printf("Try Lane Change!\n");
              printf("safe_dist = %f\n", safe_dist);

              //Check if adjacent lane is empty
              if(lane != 1)
              {
                min_dist_ahead = 10000;
                max_dist_behind = 10000;
                car_behind = -1;
                car_ahead = -1;

                new_lane = 1;
                for(int j = 0; j < cars_center_ahead.size(); j++)
                {
                  SensorFusion *car = (SensorFusion *)&cars_center_ahead[j];
                  int d = car->d;

                  if((d < (2+4*new_lane+2)) && (d > (2+4*new_lane-2)))
                  {
                    double obj_dst = car->s;

                    printf("Car ahead in center at distance = %f\n", (obj_dst - car_s));
                    if((obj_dst - car_s) < min_dist_ahead)
                    {
                      min_dist_ahead = (obj_dst - car_s);
                      car_ahead = car->id;
                    }
                  }
                }

                for(int j = 0; j < cars_center_behind.size(); j++)
                {
                  SensorFusion *car = (SensorFusion *)&cars_center_behind[j];
                  int d = car->d;

                  if((d < (2+4*new_lane+2)) && (d > (2+4*new_lane-2)))
                  {
                    double obj_dst = car->s;

                    printf("Car behind in center at distance = %f\n", (car_s - obj_dst));
                    if((car_s - obj_dst) < max_dist_behind)
                    {
                      max_dist_behind = (car_s - obj_dst);
                      car_behind = car->id;
                    }
                  }
                }

                if((car_ahead == -1) && (car_behind == -1)) {
                  lane = 1;
                }

                if((car_ahead == -1) && (car_behind >=0)) {
                  if(max_dist_behind > safe_dist) {
                    lane = new_lane;
                    printf("No cars ahead! max_dist_behind = %f  ", max_dist_behind);
                  }
                }

                if((car_behind == -1) && (car_ahead >=0)) {
                  if(min_dist_ahead > safe_dist) {
                    lane = new_lane;
                    printf("No cars behind! min_dist_ahead = %f  ", min_dist_ahead);
                  }
                }

                if((car_behind >= 0) && (car_ahead >=0)) {
                  if((min_dist_ahead > safe_dist) && (max_dist_behind > safe_dist)) {
                    lane = new_lane;
                    printf("Gap found! min_dist_ahead = %f, max_dist_behind = %f, safe_dist = %f \n", min_dist_ahead, max_dist_behind, safe_dist);
                  }
                }

              }
              else
              {
                bool left_lane_free = false;
                bool right_lane_free = false;

                min_dist_ahead = 10000;
                max_dist_behind = 10000;
                car_behind = -1;
                car_ahead = -1;

                new_lane = 0;
                for(int j = 0; j < cars_left_ahead.size(); j++)
                {
                  SensorFusion *car = (SensorFusion *)&cars_left_ahead[j];
                  int d = car->d;

                  if((d < (2+4*new_lane+2)) && (d > (2+4*new_lane-2)))
                  {
                    double obj_dst = car->s;

                    printf("Car ahead in left lane at distance = %f\n", (obj_dst - car_s));
                    if((obj_dst - car_s) < min_dist_ahead)
                    {
                      min_dist_ahead = (obj_dst - car_s);
                      car_ahead = car->id;
                    }
                  }
                }

                for(int j = 0; j < cars_left_behind.size(); j++)
                {
                  SensorFusion *car = (SensorFusion *)&cars_left_behind[j];
                  int d = car->d;
                  if((d < (2+4*new_lane+2)) && (d > (2+4*new_lane-2)))
                  {
                    double obj_dst = car->s;

                    printf("Car behind in left lane at distance = %f\n", (car_s - obj_dst));
                    if((car_s - obj_dst) < max_dist_behind)
                    {
                      max_dist_behind = (car_s - obj_dst);
                      car_behind = car->id;
                    }
                  }
                }

                if((car_ahead == -1) && (car_behind == -1)) {
                  left_lane_free = true;
                }

                if((car_ahead == -1) && (car_behind >=0)) {
                  if(max_dist_behind > safe_dist) {
                    left_lane_free = 1;
                    printf("No cars ahead! max_dist_behind = %f  \n", max_dist_behind);
                  }
                }

                if((car_behind == -1) && (car_ahead >=0)) {
                  if(min_dist_ahead > safe_dist) {
                    left_lane_free = 1;
                    printf("No cars behind! min_dist_ahead = %f  \n", min_dist_ahead);
                  }
                }

                if((car_behind >= 0) && (car_ahead >=0)) {
                  if((min_dist_ahead > safe_dist) && (max_dist_behind > safe_dist)) {
                    left_lane_free = 1;
                    printf("Gap found! min_dist_ahead = %f, max_dist_behind = %f, safe_dist = %f \n", min_dist_ahead, max_dist_behind, safe_dist);
                  }
                }

                min_dist_ahead = 10000;
                max_dist_behind = 10000;
                car_behind = -1;
                car_ahead = -1;

                new_lane = 2;
                for(int j = 0; j < cars_right_ahead.size(); j++)
                {
                  SensorFusion *car = (SensorFusion *)&cars_right_ahead[j];
                  int d = car->d;
                  if((d < (2+4*new_lane+2)) && (d > (2+4*new_lane-2)))
                  {
                    double obj_dst = car->s;

                    printf("Car ahead in right lane at distance = %f\n", (obj_dst - car_s));
                    if((obj_dst - car_s) < min_dist_ahead)
                    {
                      min_dist_ahead = (obj_dst - car_s);
                      car_ahead = car->id;
                    }
                  }
                }

                for(int j = 0; j < cars_right_behind.size(); j++)
                {
                  SensorFusion *car = (SensorFusion *)&cars_right_behind[j];
                  int d = car->d;
                  if((d < (2+4*new_lane+2)) && (d > (2+4*new_lane-2)))
                  {
                    double obj_dst = car->s;

                    printf("Car behind in right lane at distance = %f\n", (car_s - obj_dst));
                    if((car_s - obj_dst) < max_dist_behind)
                    {
                      max_dist_behind = (car_s - obj_dst);
                      car_behind = car->id;
                    }
                  }
                }

                if((car_ahead == -1) && (car_behind == -1)) {
                  right_lane_free = true;
                }

                if((car_ahead == -1) && (car_behind >=0)) {
                  if(max_dist_behind > safe_dist) {
                    right_lane_free = 1;
                    printf("No cars ahead! max_dist_behind = %f  \n", max_dist_behind);
                  }
                }

                if((car_behind == -1) && (car_ahead >=0)) {
                  if(min_dist_ahead > safe_dist) {
                    right_lane_free = 1;
                    printf("No cars behind! min_dist_ahead = %f  \n", min_dist_ahead);
                  }
                }

                if((car_behind >= 0) && (car_ahead >=0)) {
                  if((min_dist_ahead > safe_dist) && (max_dist_behind > safe_dist)) {
                    right_lane_free = 1;
                    printf("Gap found! min_dist_ahead = %f, max_dist_behind = %f, safe_dist = %f \n", min_dist_ahead, max_dist_behind, safe_dist);
                  }
                }

                printf("left_lane_free = %d, right_lane_free = %d\n", left_lane_free, right_lane_free);
                if(left_lane_free && right_lane_free)
                {
                  printf("Cars Left = %d, Right = %d\n", cars_left_ahead.size(), cars_right_ahead.size());
                  if(cars_left_ahead.size() > cars_right_ahead.size())
                  {
                    lane = 2;
                  }
                  else
                  {
                    lane = 0;
                  }
                }

                if(left_lane_free && !right_lane_free)
                {
                  lane = 0;
                }

                if(!left_lane_free && right_lane_free)
                {
                  lane = 2;
                }

                if(!left_lane_free && !right_lane_free)
                {
                  if((cars_left_ahead.size()  == 0) &&
                     (cars_right_ahead.size() == 0) &&
                     (cars_left_behind.size() == 0) &&
                     (cars_right_behind.size() == 0))
                  {
                    lane = 0;
                  }
                }
              }

              printf("Choosen lane %d\n", lane);
            }

            vector<double> pts_x;
            vector<double> pts_y;

            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);

            if(prev_path_size < 2)
            {
              pts_x.push_back(car_x - cos(car_yaw));
              pts_y.push_back(car_y - sin(car_yaw));

              pts_x.push_back(car_x);
              pts_y.push_back(car_y);
            }
            else
            {
              ref_x = previous_path_x[prev_path_size - 1];
              ref_y = previous_path_y[prev_path_size - 1];

              double ref_x_prev = previous_path_x[prev_path_size - 2];
              double ref_y_prev = previous_path_y[prev_path_size - 2];
              ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

              pts_x.push_back(ref_x_prev);
              pts_y.push_back(ref_y_prev);

              pts_x.push_back(ref_x);
              pts_y.push_back(ref_y);
            }

            vector<double> next_wp0 = getXY(car_s + 30, (2 + 4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp1 = getXY(car_s + 60, (2 + 4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp2 = getXY(car_s + 90, (2 + 4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

            pts_x.push_back(next_wp0[0]);
            pts_y.push_back(next_wp0[1]);

            pts_x.push_back(next_wp1[0]);
            pts_y.push_back(next_wp1[1]);

            pts_x.push_back(next_wp2[0]);
            pts_y.push_back(next_wp2[1]);

            for(int i = 0; i < pts_x.size(); i++)
            {
              double new_x = pts_x[i] - ref_x;
              double new_y = pts_y[i] - ref_y;

              pts_x[i] = (new_x * cos(0-ref_yaw)) - (new_y * sin(0-ref_yaw));
              pts_y[i] = (new_x * sin(0-ref_yaw)) + (new_y * cos(0-ref_yaw));

            }

            tk::spline s;

            s.set_points(pts_x, pts_y);

            vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	for(int i = 0; i < prev_path_size; i++)
          	{
          	  next_x_vals.push_back(previous_path_x[i]);
          	  next_y_vals.push_back(previous_path_y[i]);
          	}

          	double target_x = 30.0;
          	double target_y = s(target_x);

          	double target_dist = sqrt((target_x * target_x) + (target_y * target_y));
            double N = target_dist/(0.02 * ref_velocity / 2.24);
            double delta_x = target_x / N;

            double new_x = 0.0;
            double new_y = 0.0;

            for(int i = 0; i < (50 - prev_path_size); i++)
            {
              new_x = new_x + delta_x;
              new_y = s(new_x);

              double out_x = ref_x + (new_x * cos(ref_yaw)) - (new_y * sin(ref_yaw));
              double out_y = ref_y + (new_x * sin(ref_yaw)) + (new_y * cos(ref_yaw));

              next_x_vals.push_back(out_x);
              next_y_vals.push_back(out_y);
            }

          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconxds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































