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

void Consumer::get_derivatives(vec X_navigation_spacecraft, vec consumer_AES, long double t, long double curr_lon, int measure_number)
{
    //составляем матрицу H(кол-во_измер Х кол-во_оцен_парам), матрицу D (кол-во_измер Х кол-во_измер)
    // H - матрица частных производных от модели измерений в каждый момент времени

    derivatives(measure_number-1,0) = (measure_vector[2]*(sinl(measure_vector[0])*(X_navigation_spacecraft[0]*cosl(curr_lon) + X_navigation_spacecraft[1]*sinl(curr_lon))-X_navigation_spacecraft[2]*cosl(measure_vector[0])))/
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

void Consumer::navigation(std::vector<mat> finish_modeling, bool w_err)
{

    int count = 0.0;
    for(int k = 0; k < finish_modeling.size(); ++k)
    {
        long double temp_time{0.0L};
        long double current_longtitude{0.0L};
        if (!w_err)
        {
            true_distances.resize(measure_number, 3);
            derivatives.resize(measure_number, 5);
        }

        for(int i = 0; i < finish_modeling[k].n_rows; ++i)
        {
            vec temp_v(3);
            vec consumer_AES(3);
            temp_v[0] = finish_modeling[k](i,1);
            temp_v[1] = finish_modeling[k](i,2);
            temp_v[2] = finish_modeling[k](i,3);

            //доворачиваем потребитель
            current_longtitude = measure_vector[1] + omega*finish_modeling[k](i,0);

            //перевод координат радиус вектора к часам на повехрности в декартову систему, получаем вектор потребителя
            consumer_vector[0] = measure_vector[2]*cosl(measure_vector[0])*cosl(current_longtitude);
            consumer_vector[1] = measure_vector[2]*cosl(measure_vector[0])*sinl(current_longtitude);
            consumer_vector[2] = measure_vector[2]*sinl(measure_vector[0]);

            //вектор между потребителем и спутником
            consumer_AES = temp_v - consumer_vector;

            //условие видимости
            long double alfa = acosl(as_scalar(dot(consumer_vector,consumer_AES))/(norm(consumer_vector)*norm(consumer_AES)));

            if (alfa < M_PI/2.0L) visibility = true; else visibility = false;

            if( (w_err)&&(visibility)&&(temp_time + measurement_period <= finish_modeling[k](i,0)) )
            {
                ++measure_number;
                temp_time = finish_modeling[k](i,0);

                distance = norm(consumer_AES) + (measure_vector[3]+measure_vector[4]*temp_time) * light_speed +_generator.white_noise_generator(m, d);

                init_distances.resize(measure_number, 3);
                init_distances(measure_number - 1,0) = temp_time;
                init_distances(measure_number - 1,1) = k;
                init_distances(measure_number - 1,2) = distance;
            }

            if ( (!w_err)&&(count < measure_number)&&(init_distances(count, 0) == finish_modeling[k](i,0))&&(k == init_distances(count, 1)) )
            {
                ++count;
                temp_time = finish_modeling[k](i,0);

                distance = norm(consumer_AES) + (measure_vector[3]+measure_vector[4]*temp_time)*light_speed;

                true_distances(count - 1,0) = temp_time;
                true_distances(count - 1,1) = k;
                true_distances(count - 1,2) = distance;

                get_derivatives(temp_v, consumer_AES, temp_time, current_longtitude, count);
            }
        }
    }
}
