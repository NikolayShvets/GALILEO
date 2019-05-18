#include "consumer.h"

Consumer::Consumer(long double latitude, long double longtitude, long double measurement_period, long double betta_0, long double betta_1, long double m, long double d)
{
    this->measure_vector.resize(5);
    this->measure_vector[0] = latitude;
    this->measure_vector[1] = longtitude;
    this->measure_vector[2] = this->Re;
    this->measure_vector[3] = betta_0;
    this->measure_vector[4] = betta_1;
    this->measurement_period = measurement_period;
    this->m = m;
    this->d = d;
    consumer_vector.resize(3);
}
Consumer::~Consumer()
{
   navigation_log_vector.erase(navigation_log_vector.begin(), navigation_log_vector.end());
}

void Consumer::get_derivatives(vec X_navigation_spacecraft, long double t, long double curr_lon, int measure_number)
{
    //составляем матрицу H(кол-во_измер Х кол-во_оцен_парам), матрицу D (кол-во_измер Х кол-во_измер)
    // H - матрица частных производных от модели измерений в каждый момент времени
    derivatives.resize(measure_number, 5);
    derivatives(measure_number-1,0) = (measure_vector[2]*(sinl(measure_vector[0])*(X_navigation_spacecraft[0]*cosl(curr_lon) + X_navigation_spacecraft[1]*sinl(measure_vector[2])-X_navigation_spacecraft[2]*cosl(measure_vector[0]))))/
            (sqrtl(powl((X_navigation_spacecraft[0] - measure_vector[2]*cosl(curr_lon)*cosl(measure_vector[0])), 2.0L) +
            powl((X_navigation_spacecraft[1] - measure_vector[2]*sinl(curr_lon)*cosl(measure_vector[0])),2.0L) +
            powl((X_navigation_spacecraft[2] - measure_vector[2]*sinl(measure_vector[0])),2.0L)));
    derivatives(measure_number-1,1) = (-1)*(measure_vector[2]*cosl(measure_vector[0])*(X_navigation_spacecraft[1]*cosl(curr_lon)-X_navigation_spacecraft[0]*sinl(curr_lon)))/
            (sqrtl(powl((X_navigation_spacecraft[0] - measure_vector[2]*cosl(curr_lon)*cosl(measure_vector[0])), 2.0L) +
            powl((X_navigation_spacecraft[1] - measure_vector[2]*sinl(curr_lon)*cosl(measure_vector[0])),2.0L) +
            powl((X_navigation_spacecraft[2] - measure_vector[2]*sinl(measure_vector[0])),2.0L)));
    derivatives(measure_number-1,2) = (measure_vector[2]-X_navigation_spacecraft[0]*cosl(curr_lon)*cosl(measure_vector[0])-X_navigation_spacecraft[1]*sinl(curr_lon)*cosl(measure_vector[0])-X_navigation_spacecraft[2]*sinl(measure_vector[0]))/
            (sqrtl(powl((X_navigation_spacecraft[0] - measure_vector[2]*cosl(curr_lon)*cosl(measure_vector[0])), 2.0L) +
            powl((X_navigation_spacecraft[1] - measure_vector[2]*sinl(curr_lon)*cosl(measure_vector[0])),2.0L) +
            powl((X_navigation_spacecraft[2] - measure_vector[2]*sinl(measure_vector[0])),2.0L)));
    derivatives(measure_number-1,3) = light_speed;
    derivatives(measure_number-1,4) = light_speed * t;
}

void Consumer::navigation(std::vector<mat> finish_modeling, bool init_dist, bool w_o_err, bool d_and_der)
{
    for(int k = 0; k < finish_modeling.size(); ++k)
    {
        if(!init_dist)
            navigation_log_vector.push_back(std::ofstream(std::to_string(k) + "_distance.txt"));
        long double temp_time{0.0L};
        long double current_longtitude{0.0L};
        for(int i = 0; i < finish_modeling[k].n_rows; ++i)
        {
            for (int j = 0; j < finish_modeling[k].n_cols; ++j)
            {
                vec temp_v(3);
                temp_v[0] = finish_modeling[k](i,1);
                temp_v[1] = finish_modeling[k](i,2);
                temp_v[2] = finish_modeling[k](i,3);
           //     std::cout<<temp_v<<std::endl;
                //доворачиваем потребитель
                current_longtitude = measure_vector[1] + omega*finish_modeling[k](i,0);
                //перевод координат радиус вектора к часам на повехрности в декартову систему
                consumer_vector[0] = measure_vector[2]*cosl(measure_vector[0])*cosl(current_longtitude);
                consumer_vector[1] = measure_vector[2]*cosl(measure_vector[0])*sinl(current_longtitude);
                consumer_vector[2] = measure_vector[2]*sinl(measure_vector[0]);
                //условие видимости
                long double alfa = acosl(as_scalar(dot(temp_v,consumer_vector))/(norm(temp_v)*norm(consumer_vector)));

                if (alfa < M_PI/3.0L) visibility = true; else visibility = false;

                if( (init_dist == true)&&(visibility)&&(temp_time + measurement_period <= finish_modeling[k](i,0)) )
                {
                    ++measure_number;
                    temp_time = finish_modeling[k](i,0);
                    distance = sqrtl(powl((temp_v[0] - measure_vector[2]*cosl(measure_vector[0])*cosl(current_longtitude)), 2.0L) +
                            powl((temp_v[1] - measure_vector[2]*cosl(measure_vector[0])*sinl(current_longtitude)),2.0L) +
                            powl((temp_v[2] -  measure_vector[2]*sinl(measure_vector[0])),2.0L)) +
                                (betta_0+betta_1*temp_time)*light_speed+
                            _generator.white_noise_generator(m, d);

                    init_distances.resize(measure_number, 3);
                    for(int q = measure_number - 1; q < measure_number; ++q)
                    {
                        init_distances(q,0) = temp_time;
                        init_distances(q,1) = k;
                        init_distances(q,2) = distance;
                    }
                }

                if ( (w_o_err)&&(temp_time + measurement_period <= finish_modeling[k](i,0)) )
                {
                    ++measure_number;
                    temp_time = finish_modeling[k](i,0);
                    distance = sqrtl(powl((temp_v[0] - measure_vector[2]*cosl(measure_vector[0])*cosl(current_longtitude)), 2.0L) +
                            powl((temp_v[1] - measure_vector[2]*cosl(measure_vector[0])*sinl(current_longtitude)),2.0L) +
                            powl((temp_v[2] -  measure_vector[2]*sinl(measure_vector[0])),2.0L))+
                            (betta_0+betta_1*temp_time)*light_speed;

                    true_distances.resize(measure_number, 3);
                    for(int q = measure_number - 1; q < measure_number; ++q)
                    {
                        true_distances(q,0) = temp_time;
                        true_distances(q,1) = k;
                        true_distances(q,2) = distance;
                    }
                    if(d_and_der)
                    {
                        get_derivatives(temp_v, temp_time, current_longtitude, measure_number);
                        //std::cout<<derivatives<<std::endl;
                    }
                }

                else break;
                if(!init_dist)
                {
                    navigation_log_vector.back()<<std::fixed;
                    navigation_log_vector.back()<<finish_modeling[k](i,0)<<" "<<distance<<std::endl;
                }
                break;
            }
        }
    }
}
