#include "consumer.h"

Consumer::Consumer(long double latitude, long double longtitude, long double measurement_period, long double betta_0, long double betta_1)
{
    this->latitude = latitude;
    this->longtitude = longtitude;
    this->betta_0 = betta_0;
    this->betta_1 = betta_1;
    this->measurement_period = measurement_period;
    consumer_vector.resize(3);
}
Consumer::~Consumer()
{
   navigation_log_vector.erase(navigation_log_vector.begin(), navigation_log_vector.end());
}

void Consumer::get_derivatives(TVector X_navigation_spacecraft, long double t, int measure_number, int satellite_number)
{
    //инкременитровать первый индекс, если изменился номер спутника
    //derivatives.back().resize(measure_number,5);

    if(current_s_num != satellite_number)
    {
        s = 1;
        ++k;
        TMatrix temp;
        derivatives.push_back(temp);
    }
    //ресайзить не на количество измерений в общем, а на количество измерений именно с этим номером спутника
derivatives[k].resize(/*measure_number*/s, 5);
    derivatives[k][/*measure_number*/s-1][0] = (X_navigation_spacecraft[0] - consumer_vector[0])/
            (sqrt(pow((X_navigation_spacecraft[0] - consumer_vector[0]), 2.0L) +
            pow((X_navigation_spacecraft[1] - consumer_vector[1]),2.0L) +
            pow((X_navigation_spacecraft[2] - consumer_vector[2]),2.0L)));
    derivatives[k][/*measure_number*/s-1][1] = (X_navigation_spacecraft[1] - consumer_vector[1])/
            (sqrt(pow((X_navigation_spacecraft[0] - consumer_vector[0]), 2.0L) +
            pow((X_navigation_spacecraft[1] - consumer_vector[1]),2.0L) +
            pow((X_navigation_spacecraft[2] - consumer_vector[2]),2.0L)));
    derivatives[k][/*measure_number*/s-1][2] = (X_navigation_spacecraft[2] - consumer_vector[2])/
            (sqrt(pow((X_navigation_spacecraft[0] - consumer_vector[0]), 2.0L) +
            pow((X_navigation_spacecraft[1] - consumer_vector[1]),2.0L) +
            pow((X_navigation_spacecraft[2] - consumer_vector[2]),2.0L)));
    derivatives[k][/*measure_number*/s-1][3] = light_speed;
    derivatives[k][/*measure_number*/s-1][4] = light_speed * t;

    current_s_num = satellite_number;
    ++s;
}

void Consumer::navigation(std::vector<TMatrix> finish_modeling, bool flag)
{
    for(int k = 0; k < finish_modeling.size(); ++k)
    {
        navigation_log_vector.push_back(std::ofstream(std::to_string(k) + "_distance.txt"));
        long double temp_time{0.0L};
        long double current_longtitude{0.0L};
        for(int i = 0; i < finish_modeling[k].row_count(); ++i)
        {
            for (int j = 0; j < finish_modeling[k].col_count(); ++j)
            {
                TVector temp_v(3);
                temp_v[0] = finish_modeling[k][i][1];
                temp_v[1] = finish_modeling[k][i][2];
                temp_v[2] = finish_modeling[k][i][3];
                //доворачиваем потребитель
                current_longtitude = longtitude + omega*finish_modeling[k][i][0];
                //перевод координат радиус вектора к часам на повехрности в декартову систему
                consumer_vector[0] = Re*cos(latitude)*cos(current_longtitude);
                consumer_vector[1] = Re*cos(latitude)*sin(current_longtitude);
                consumer_vector[2] = Re*sin(latitude);
                long double alfa = acos(temp_v*consumer_vector/(temp_v.length()*consumer_vector.length()));

                if (alfa < M_PI/3.0L) visibility = true; else visibility = false;

                if ( ((visibility == true)&&(temp_time + measurement_period <= finish_modeling[k][i][0])) )/*||
                     (flag == true)&&(temp_time + measurement_period <= finish_modeling[k][i][0]) )*/
                {//рассчет псевдодальности
                    ++measure_number;
                    temp_time = finish_modeling[k][i][0];
                    distance = sqrt(pow((finish_modeling[k][i][1] - consumer_vector[0]), 2.0L) +
                        pow((finish_modeling[k][i][2] - consumer_vector[1]),2.0L) +
                        pow((finish_modeling[k][i][3] - consumer_vector[2]),2.0L)) +
                                (betta_0+betta_1*finish_modeling[k][i][0])*light_speed +
                                _generator.white_noise_generator(0.0L, 1.0L/3.0L);

                    true_distances.resize(measure_number, 3);

                    for(int q = measure_number - 1; q < measure_number; ++q)
                    {
                        true_distances[q][0] = temp_time;
                        true_distances[q][1] = k;
                        true_distances[q][2] = distance;
                    }
                    //матрица производных
                    if(flag)
                    {
                        get_derivatives(temp_v, temp_time, measure_number, k );
                        //матрица дисперсий
                        D.resize(get_measure_number(),get_measure_number());
                        for(int d_i = 0; d_i < D.row_count(); ++d_i)
                        {
                            for(int d_j = 0; d_j < D.col_count(); ++d_j)
                            {
                                if(d_i == d_j)
                                    D[d_i][d_j] = 1.0L;
                                else
                                    D[d_i][d_j] = 0.0L;
                            }
                        }
                    }
                }
                else break;
                navigation_log_vector.back()<<std::fixed;
                navigation_log_vector.back()<<finish_modeling[k][i][0]<<" "<<distance<<std::endl;
                break;
            }
        }
    }
}
