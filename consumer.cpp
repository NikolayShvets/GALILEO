#include "consumer.h"

Consumer::Consumer(long double latitude, long double longtitude, long double frequency, long double betta_0, long double betta_1)
{
    this->latitude = latitude;
    this->longtitude = longtitude;
    this->betta_0 = betta_0;
    this->betta_1 = betta_1;
    this->frequency = frequency;
    consumer_vector.resize(3);
}
Consumer::~Consumer()
{
   navigation_log_vector.erase(navigation_log_vector.begin(), navigation_log_vector.end());
}
void Consumer::navigation(TMatrix result,int model_number)
{
    long double temp_time{0.0L};
    navigation_log_vector.push_back(std::ofstream(std::to_string(model_number) + "_distance.txt"));
    for(int i = 0; i < result.row_count(); ++i)
    {
        for(int j = 0; j < result.col_count(); ++j)
        {
            TVector temp_v(3);
            temp_v[0] = result[i][1];
            temp_v[1] = result[i][2];
            temp_v[2] = result[i][3];
            //доворачиваем потребитель
            longtitude = omega*result[i][0];
            //перевод координат радиус вектора к часам на повехрности в декартову систему
            consumer_vector[0] = Re*cos(latitude)*cos(longtitude);
            consumer_vector[1] = Re*cos(latitude)*sin(longtitude);
            consumer_vector[2] = Re*sin(latitude);
            long double alfa = acos(temp_v*consumer_vector/(temp_v.length()*consumer_vector.length()));

            if (alfa < M_PI/3.0L) visibility = true; else visibility = false;

            if ((visibility == true)&&(temp_time + frequency <= result[i][0]))
            {//рассчет псевдодальности
                temp_time = result[i][0];
                distance = sqrt(pow((result[i][1] - consumer_vector[0]), 2.0L) +
                    pow((result[i][2] - consumer_vector[1]),2.0L) +
                    pow((result[i][3] - consumer_vector[2]),2.0L)) +
                            (betta_0+betta_1*result[i][0])*light_speed +
                            _generator.white_noise_generator(0.0L, 10.0L);
            }
            else break;

            navigation_log_vector.back()<<std::fixed;
            navigation_log_vector.back()<<result[i][0]<<" "<<distance<<std::endl;
            break;
        }
    }

    navigation_log_vector.back().close();
}
