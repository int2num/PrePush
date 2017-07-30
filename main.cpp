#include <iostream>
#include <iostream>
#include"Graph.h"
#include<sys/time.h>
int main()
{
    struct timeval starttime, endtime;
    PreflowPush plo=PreflowPush();
    ERGraph graph(2000,15,plo);
    gettimeofday(&endtime,NULL);
    cout<<"time is:"<<endtime.tv_sec- starttime.tv_sec<<endl;


}
