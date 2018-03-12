#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "solver.h"
#include "QDebug"
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    //Рисуем график y=x*x

     //Сгенерируем данные
     //Для этого создадим два массива точек:
     //один для созранения x координат точек,
     //а второй для y соответственно

     const double a = 0; //Начало интервала, где рисуем график по оси Ox
     const double b =  100; //Конец интервала, где рисуем график по оси Ox
     const double h = 1; //Шаг, с которым будем пробегать по оси Ox

     const int N=100; //Вычисляем количество точек, которые будем отрисовывать
     double x[N],y[N]; //Массивы координат точек
     for(int i=0;i<N;i++){
         x[i]=a+h*i;
     }
     //Вычисляем наши данные
     solver Solver;
     Solver.iterate(x,y,N);

     ui->widget->clearGraphs();//Если нужно, но очищаем все графики
     //Добавляем один график в widget
     ui->widget->addGraph();
     //Говорим, что отрисовать нужно график по нашим двум массивам x и y
     QVector<qreal> qx(N),qy(N);
     for(int i=0;i<N;i++){
         qx[i]=x[i];
         qy[i]=y[i];
       //  qDebug()<<qy[i];
     }

     ui->widget->graph(0)->setData(qx, qy);

     //Подписываем оси Ox и Oy
     ui->widget->xAxis->setLabel("x");
     ui->widget->yAxis->setLabel("y");

     //Установим область, которая будет показываться на графике
     ui->widget->xAxis->setRange(a, b);//Для оси Ox

     //Для показа границ по оси Oy сложнее, так как надо по правильному
     //вычислить минимальное и максимальное значение в векторах
     double minY = y[0], maxY = y[0];
     for (int i=1; i<N; i++)
     {
         if (y[i]<minY) minY = y[i];
         if (y[i]>maxY) maxY = y[i];
     }
     ui->widget->yAxis->setRange(minY, maxY);//Для оси Oy

     //И перерисуем график на нашем widget
     ui->widget->replot();
}
