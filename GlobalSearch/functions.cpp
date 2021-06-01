#include "functions.h"

using namespace std;

Interval::Interval()
{
  int s = 10000;
  zip_x = new pair<double, double>[s-1];
  z = new pair<double, double>[s - 1];
  R = new double[s - 1];
  HeapSize = 0;
}

Interval::Interval(int maxTrials)
{
  //�� ����� ������� ������� ������ ��������, ������� �� ����� ���������� maxTrials ��� �������� �� ������� ����� ���������� �������� ������
  zip_x = new pair<double, double>[4 * maxTrials];
  z = new pair<double, double>[4* maxTrials];
  R = new double[4 * maxTrials];
  HeapSize = 0;
}
Interval::~Interval()
{
  delete[] zip_x;
  delete[] R;
  delete[] z;
}

void Interval::out()
{
  int i = 0;
  int k = 1;
  while (i < HeapSize)
  {
    while ((i < k) && (i < HeapSize))
    {
      cout << R[i] << " (" << zip_x[i].first << ", " << zip_x[i].second << ")\t";
      i++;
    }
    cout << endl;
    k = k * 2 + 1;
  }
}

void Interval::sort_out()
{
  int key = 0;
  for (int i = 0; i < HeapSize; i++)
  {
    key = i;
    for (int j = i + 1; j < HeapSize - 1; j++)
    {
      if (R[j] > R[i])
        key = j;
    }
    cout << R[i] << " (" << zip_x[i].first << ", " << zip_x[i].second << ")\t";
  }
  cout << endl;
}

void Interval::AddElem(double _R, double x_l, double x_r, double z_l, double z_r)
{
  int i = HeapSize;
  R[i] = _R;
  zip_x[i].first = x_l;
  zip_x[i].second = x_r;
  z[i].first = z_l;
  z[i].second = z_r;
  int parent = (i - 1) / 2;

  while (parent >= 0 && i > 0)
  {
    if (R[i] > R[parent])
    {
      double temp = R[i];
      R[i] = R[parent];
      R[parent] = temp;

      temp = zip_x[i].first;
      zip_x[i].first = zip_x[parent].first;
      zip_x[parent].first = temp;

      temp = zip_x[i].second;
      zip_x[i].second = zip_x[parent].second;
      zip_x[parent].second = temp;

      temp = z[i].first;
      z[i].first = z[parent].first;
      z[parent].first = temp;

      temp = z[i].second;
      z[i].second = z[parent].second;
      z[parent].second = temp;
    }
    i = parent;
    parent = (i - 1) / 2;
  }
  HeapSize++;
}

void Interval::Update(double M, double r, int demension, double min)
{
  double N = 1.0 / demension;
  double delta = 0;
  double second = 0;
  double third = 0;

  for (int i = 0; i < HeapSize; i++)
  {
    if (zip_x[i].first == 0)
    {
      delta = pow((zip_x[i].second - zip_x[i].first), N);
      R[i] = 2 * delta - 4 * (z[i].second - min) / (r*M);
    }
    else if (zip_x[i].second == 1)
    {
      delta = pow((zip_x[i].second - zip_x[i].first), N);
      R[i] = 2 * delta - 4 * (z[i].first - min) / (r*M);
    }
    else
    {
      delta = pow((zip_x[i].second - zip_x[i].first), N);
      second = (z[i].second - z[i].first)*(z[i].second - z[i].first) / (r*r*M*M*delta);
      third = 2 * (z[i].second + z[i].first - 2 * min) / (r*M);
      R[i] = delta + second - third;
    }


    
    int parent = (i - 1) / 2;
    int j = i;

    while (parent >= 0 && j > 0)
    {
      if (R[j] > R[parent])
      {
        double temp = R[j];
        R[j] = R[parent];
        R[parent] = temp;

        temp = zip_x[j].first;
        zip_x[j].first = zip_x[parent].first;
        zip_x[parent].first = temp;

        temp = zip_x[j].second;
        zip_x[j].second = zip_x[parent].second;
        zip_x[parent].second = temp;

        temp = z[j].first;
        z[j].first = z[parent].first;
        z[parent].first = temp;

        temp = z[j].second;
        z[j].second = z[parent].second;
        z[parent].second = temp;
      }
      j = parent;
      parent = (j - 1) / 2;
    }
  }
}

double Interval::GetMaxPointLeft()
{
  return zip_x[0].first;
}

double Interval::GetMaxPointRight()
{
  return zip_x[0].second;
}

double Interval::GetMaxZLeft()
{
  return z[0].first;
}

double Interval::GetMaxZRight()
{
  return z[0].second;
}

double Interval::GetMax()
{
  double x = R[0];
  R[0] = R[HeapSize - 1];
  zip_x[0].first = zip_x[HeapSize - 1].first;
  zip_x[0].second = zip_x[HeapSize - 1].second;
  z[0].first = z[HeapSize - 1].first;
  z[0].second = z[HeapSize - 1].second;
  HeapSize--;
  Heapify(0);
  return(x);
}

int Interval::GetSize()
{
  return HeapSize;
}

void Interval::Heapify(int i)
{
  int left = 2 * i + 1;
  int right = 2 * i + 2;
  double temp;
  if (left < HeapSize)
  {
    if (R[i] < R[left])
    {
      temp = R[i];
      R[i] = R[left];
      R[left] = temp;

      temp = zip_x[i].first;
      zip_x[i].first = zip_x[left].first;
      zip_x[left].first = temp;

      temp = zip_x[i].second;
      zip_x[i].second = zip_x[left].second;
      zip_x[left].second = temp;

      temp = z[i].first;
      z[i].first = z[left].first;
      z[left].first = temp;

      temp = z[i].second;
      z[i].second = z[left].second;
      z[left].second = temp;

      Heapify(left);
    }
  }

  if (right < HeapSize)
  {
    if (R[i] < R[right])
    {
      temp = R[i];
      R[i] = R[right];
      R[right] = temp;

      temp = zip_x[i].first;
      zip_x[i].first = zip_x[right].first;
      zip_x[right].first = temp;

      temp = zip_x[i].second;
      zip_x[i].second = zip_x[right].second;
      zip_x[right].second = temp;

      temp = z[i].first;
      z[i].first = z[right].first;
      z[right].first = temp;

      temp = z[i].second;
      z[i].second = z[right].second;
      z[right].second = temp;

      Heapify(right);
    }
  }
}
//----------------------����� �������� ������--------------------------------------------------

double Perebor(int maxTrial, int demension, IProblem* problem, double* *bestX)
{
  int size;
  int rank;
  int delta;
  int ost;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double min, sum;

  double* lower_bounds = new double[demension];
  double* upper_bounds = new double[demension];
  problem->GetBounds(lower_bounds, upper_bounds);
  double* shag = new double[demension];
  for (int i =0; i < demension; i++)
    shag[i] = lower_bounds[i];
  double** masX = new double*[maxTrial + 1];
  for (int i = 0; i < maxTrial + 1; i++)
    masX[i] = new double[demension];
  for (int i = 0; i < demension; i++)
    masX[0][i] = shag[i];

  //���������� ���������� X � �����, � ������� ����� ����� ������� �������� �������
  for (int i = 1; i < maxTrial + 1; i++)
  {
    for (int j = 0; j < demension; j++) {
      if (shag[j] >= upper_bounds[j] + 0.1)
        break;
      shag[j] += (upper_bounds[j] - lower_bounds[j]) / (maxTrial + 1);
      masX[i][j] = shag[j];
    }
  }

  int flag = 0;
  int tmpsize = 0;
  ost = (maxTrial + 1) % size;
  if ((maxTrial + 1) <= size)
  {
    tmpsize = size;
    size = 1;
    ost = 0;
    flag = -1;
  }

  delta = (maxTrial + 1) / size;

  double prom_min = 0;
  if (rank == 0)
  {
    double* tmp = new double[demension];

    for (int i = 0; i < delta + ost; i++)
    {
      if (demension == 2) {
        tmp[0] = masX[i][0];
        for (int j = 0; j < delta + ost; j++) {
          tmp[1] = masX[j][1];

          sum = problem->CalculateFunctionals(tmp, 0);

          if ((j == 0)&&(i==0)) {
            prom_min = sum;
            (*bestX)[0] = tmp[0];
            (*bestX)[1] = tmp[1];
          }
          else
            if (sum < prom_min) {
              prom_min = sum;
              (*bestX)[0] = tmp[0];
              (*bestX)[1] = tmp[1];
            }
        }
      }
      if (demension == 1) {
        for (int j = 0; j < demension; j++)
          tmp[j] = masX[i][j];
        sum = problem->CalculateFunctionals(tmp, 0);

        if (i == 0) {
          prom_min = sum;
          (*bestX)[0] = tmp[0];
        }
        else
          if (sum < prom_min) {
            prom_min = sum;
            (*bestX)[0] = tmp[0];
          }
      }
    }
    min = prom_min;
    //delete[] tmp;
  }
  //else
  //{
  //  if (flag == 0)
  //  {
  //    double* tmp = new double[demension];
  //    for (int i = delta*rank + ost; i < delta* rank + ost + delta; i++)
  //    {
  //      if (demension == 2) {
  //        tmp[0] = masX[i][0];
  //        for (int j = delta * rank + ost; j < delta* rank + ost + delta; j++) {
  //          tmp[1] = masX[j][1];

  //          sum = problem->CalculateFunctionals(tmp, 0);

  //          if (i == delta * rank + ost) {
  //            prom_min = sum;
  //          }
  //          else
  //            if (sum < prom_min)
  //              prom_min = sum;
  //        }
  //      }
  //      if (demension == 1) {
  //        for (int j = 0; j < demension; j++)
  //          tmp[j] = masX[i][j];
  //        sum = problem->CalculateFunctionals(tmp, 0);

  //        if (i == 0) {
  //          prom_min = sum;
  //        }
  //        else
  //          if (sum < prom_min)
  //            prom_min = sum;
  //      }
  //    }
  //    //delete[] tmp;
  //  }
  //}

  if (flag == 0)
  {
    MPI_Reduce(&prom_min, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  }
  else
  {
    if (rank == 0)
    {
      return prom_min;
    }
    else
    {
      return -1;
    }
  }
  
  delete[] masX;
  delete[] lower_bounds;
  delete[] upper_bounds;
  return min;
}

void QuickSort(double* mas, int first, int last)
{
  double count;
  int f = first, l = last;
  double mid = mas[(f + l) / 2]; //���������� �������� ��������

  do
  {
    while (mas[f] < mid)
      f++;
    while (mas[l] > mid)
      l--;

    if (f <= l) //������������ ���������
    {
      count = mas[f];
      mas[f] = mas[l];
      mas[l] = count;

      f++;
      l--;
    }
  } while (f < l);
  if (first < l)
    QuickSort(mas, first, l);
  if (f < last)
    QuickSort(mas, f, last);
}

void QuickSort_D(double* mas, int first, int last, double* z)
{
  double count;
  int f = first, l = last;
  double mid = mas[(f + l) / 2]; //���������� �������� ��������

  do
  {
    while (mas[f] < mid)
      f++;
    while (mas[l] > mid)
      l--;

    if (f <= l) //������������ ���������
    {
      count = mas[f];
      mas[f] = mas[l];
      mas[l] = count;

      count = z[f];
      z[f] = z[l];
      z[l] = count;

      f++;
      l--;
    }
  } while (f < l);
  if (first < l)
    QuickSort_D(mas, first, l, z);
  if (f < last)
    QuickSort_D(mas, f, last, z);
}

int sng(double x) {
  if (x > 0)
    return 1;
  else
    return -1;
}

double AGP(int maxTrial, double accuracy, int demension, IProblem* problem)
{
  double* lower_bounds = new double[demension];
  double* upper_bounds = new double[demension];
  problem->GetBounds(lower_bounds, upper_bounds);

  double* X = new double[maxTrial];
  for (int i = 0; i < maxTrial; i++)
    X[i] = 0;
  X[0] = lower_bounds[0];
  X[1] = upper_bounds[0];

  //��� �������� ��������
  double* tmpmin = new double[demension];
  tmpmin[0] = X[0];
  double min = problem->CalculateFunctionals(tmpmin, 0);
  double prom_min;
  int iteration = 0;

  int first = 0;
  int last = 1;
  int k = 2;

  double* tmp = new double[demension];
  double z0;
  double z1;
  double M;

  int t = 1;
  double r = 4;
  double delta;
  double second;
  double third;
  double R;

  for (int i = 0; i < maxTrial; i++)
  {
    //����������� ����� �� ����������
    QuickSort(X, first, last);

    t = 1;

    //��������� ������ M ��� ����������� ��������� �������
    for (int j = 1; j < k; j++)
    {
      double m;
      tmp[0] = X[j - 1];
      z0 = problem->CalculateFunctionals(tmp, 0);
      tmp[0] = X[j];
      z1 = problem->CalculateFunctionals(tmp, 0);
      m = fabs(z1 - z0) / (X[j] - X[j - 1]);
      if (j == 1)
        M = m;
      if (m > M)
        M = m;
    }

    //�������� �������������� R(i);
    for (int j = 1; j < k; j++)
    {
      double Rtmp;
      tmp[0] = X[j - 1];
      z0 = problem->CalculateFunctionals(tmp, 0);
      tmp[0] = X[j];
      z1 = problem->CalculateFunctionals(tmp, 0);
      delta = X[j] - X[j - 1];
      second = (z1 - z0)*(z1 - z0) / (r*r*M*M*delta);
      third = 2 * (z1 + z0) / (r*M);
      Rtmp = delta + second - third;
      if (j == 1)
        R = Rtmp;
      if (Rtmp > R)
      {
        R = Rtmp;
        t = j;
      }
    }

    k++;
    tmp[0] = X[t - 1];
    z0 = problem->CalculateFunctionals(tmp, 0);
    tmp[0] = X[t];
    z1 = problem->CalculateFunctionals(tmp, 0);
    second = (X[t] + X[t - 1]) / 2;
    third = (z1 - z0) / (2 * r*M);
    X[k - 1] = second - third;

    tmpmin[0] = X[k - 1];
    prom_min = problem->CalculateFunctionals(tmpmin, 0);
    if (prom_min < min)
      min = prom_min;
    iteration++;

    if (fabs(X[k - 1] - X[t]) < accuracy)
      break;
    else if (fabs(X[k - 1] - X[t - 1]) < accuracy)
      break;
    last++;
  }

  printf("Number of iteration = %d\n", iteration);
  //delete[] tmp;
  //delete[] tmpmin;
  delete[] lower_bounds;
  delete[] upper_bounds;

  return min;
}


// ��� ������������ ������------------------------------------------------------------
double AGP_Space(int maxTrial, double accuracy, int demension, IProblem* problem, double* *BestX, int *iteration, double r, double* Truemin, int flag, char* argv[])
{
  double N = 1.0 / demension;

  double* lower_bounds = new double[demension];
  double* upper_bounds = new double[demension];
  problem->GetBounds(lower_bounds, upper_bounds);

  //����������� �����
  double** X = new double*[maxTrial];
  for (int i = 0; i < maxTrial; i++)
  {
    X[i] = new double[demension];
    for (int j = 0; j < demension; j++)
      X[i][j] = 0;
  }

  double z_l, z_r, z;

  Interval inter(maxTrial);

  //��� ���� ��������� �������� �� ����� ������� � ����� � ������� [0; 1] � �� �� ����� � ������ ��������
  double zip_x_l = 0;
  double zip_x_r = 1;
  double zip_x = 0;

  //��� �������� ��������
  double min;
  GetImage(zip_x_l, X[0], lower_bounds, upper_bounds, demension, 10);
  min = problem->CalculateFunctionals(X[0], 0);
  z_l = min;

  for (int i = 0; i < demension; i++)
    (*BestX)[i] = X[0][i];

  GetImage(zip_x_r, X[0], lower_bounds, upper_bounds, demension, 10);
  z_r = problem->CalculateFunctionals(X[0], 0);

  if (z_r < min)
  {
    min = z_r;
    for (int i = 0; i < demension; i++)
      (*BestX)[i] = X[0][i];
  }

  double prom_min;

  double M = 1;

  // ���� �� ���������� M
  int f = 0;

  double delta;
  double second;
  double third;
  double R = DBL_MIN;
  double p;

  ////��������� ���� ��� ������-------------------------------------------------------------------
  //std::string pointLogName = argv[3];
  //FILE* pointLog = fopen(pointLogName.c_str(), "w");
  ////--------------------------------------------------------------------------------------------

  for (int i = 0; i < maxTrial; i++)
  {
    /*inter.sort_out();
    cout << "-------------------------------------\n";*/
    if (inter.GetSize() != 0)
    {
      zip_x_l = inter.GetMaxPointLeft();
      zip_x_r = inter.GetMaxPointRight();
      z_l = inter.GetMaxZLeft();
      z_r = inter.GetMaxZRight();
    }

    if (zip_x_l != 0 && zip_x_r != 1) {
      second = (zip_x_r + zip_x_l) / 2;
      p = pow((fabs(z_r - z_l) / M), demension);
      third = p / 2 / r;
      zip_x = second - sng(z_r - z_l) * third;
    }
    else
      zip_x = (zip_x_r + zip_x_l) / 2;


    GetImage(zip_x, X[i], lower_bounds, upper_bounds, demension, 10);
    prom_min = problem->CalculateFunctionals(X[i], 0);
    z = prom_min;
    if (prom_min < min) {
      min = prom_min;
      for (int kk = 0; kk < demension; kk++)
        (*BestX)[kk] = X[i][kk];
    }
    (*iteration)++;

    if ((zip_x >= zip_x_r) || (zip_x <= zip_x_l))
      std::cout << "ERROR!!!!!!!\n";

    if (flag == 0) {
      p = pow(fabs(zip_x - zip_x_r), N);
      if (p < accuracy)
        break;
      p = pow(fabs(zip_x - zip_x_l), N);
      if (p < accuracy)
        break;
    }
    else {
      /*double maxfabs = fabs(Truemin[0] - (*BestX)[0]);
      if (fabs(Truemin[1] - (*BestX)[1]) > maxfabs)
        maxfabs = fabs(Truemin[1] - (*BestX)[1]);
      cout << i << "; " << maxfabs << "\n";*/
      int br = 0;
      for (int kk = 0; kk < demension; kk++)
        if (fabs(Truemin[kk] - (*BestX)[kk]) < accuracy)
          br++;
      if (br == demension)
        break;
    }

    //��������� ������ M ��� ����������� ��������� �������
    double m = 1;
    p = pow(fabs(zip_x - zip_x_r), N);
    m = fabs(z - z_r) / p;
    if (m > M)
    {
      M = m;
      f = 1;
    }
    p = pow(fabs(zip_x - zip_x_l), N);
    m = fabs(z - z_l) / p;
    if (m > M)
    {
      M = m;
      f = 1;
    }

    if (i == 0) {
      M = 1;
      f = 1;
    }

    //�������� �������������� R(i)
    R = DBL_MIN;

    if (inter.GetSize() != 0)
    {
      inter.GetMax();
      //cout << "After:\n";
      //inter.out();
      //cout << "-------------------------------------\n\n";
    }

    if (zip_x_l == 0)
    {
      delta = pow((zip_x - zip_x_l), N);
      R = 2 * delta - 4 * (z - min) / (r*M);
    }
    else
    {
      delta = pow((zip_x - zip_x_l), N);
      second = (z - z_l)*(z - z_l) / (r*r*M*M*delta);
      third = 2 * (z + z_l - 2 * min) / (r*M);
      R = delta + second - third;
    }

    inter.AddElem(R, zip_x_l, zip_x, z_l, z);

    //-----------------------------------------------
    R = DBL_MIN;

    if (zip_x_r == 1)
    {
      delta = pow((zip_x_r - zip_x), N);
      R = 2 * delta - 4 * (z - min) / (r*M);
    }
    else
    {
      delta = pow((zip_x_r - zip_x), N);
      second = (z_r - z)*(z_r - z) / (r*r*M*M*delta);
      third = 2 * (z_r + z - 2 * min) / (r*M);
      R = delta + second - third;
    }

    inter.AddElem(R, zip_x, zip_x_r, z, z_r);

    if (f == 1)
    {
      inter.Update(M, r, demension, min);
      f = 0;
    }

  }

  ////------------------------------------------------------------------------------------------
  //fprintf(pointLog, "%d\n", (*iteration));
  ////���������� � ���� �����----------------------------------------------------------------------------------
  //for (int i = 1; i < (*iteration); i++)
  //{
  //  for (int j = 0; j < demension; j++)
  //    fprintf(pointLog, "%lf ", X[i][j]);
  //  fprintf(pointLog, "\n");
  //}
  ////---------------------------------------------------------------------------------------------------------------

  ///// �������� ��������� �����
  //for (int j = 0; j < demension; j++)
  //  fprintf(pointLog, "%lf ", (*BestX)[j]);
  //fprintf(pointLog, "\n");

  ///// ��������� ����� ����������� ��������
  //for (int j = 0; j < demension; j++)
  //  fprintf(pointLog, "%lf ", Truemin[j]);
  //fclose(pointLog);
  ////------------------------------------------------------------------------------------------

  delete[] lower_bounds;
  delete[] upper_bounds;

  return min;
}

// ��� ������������ ������������� ������------------------------------------------------------------
double AGP_Space_OMP(int maxTrial, double accuracy, int demension, IProblem* problem, double* *BestX, int *iteration, double r, double* Truemin, int NumThr, int flag, char* argv[])
{
  double N = 1.0 / demension;

  double* lower_bounds = new double[demension];
  double* upper_bounds = new double[demension];
  problem->GetBounds(lower_bounds, upper_bounds);

  //����������� �����
  double** X = new double*[NumThr];
  for (int i = 0; i < NumThr; i++)
  {
    X[i] = new double[demension];
    for (int j = 0; j < demension; j++)
      X[i][j] = 0;
  }

  ////=======================================================================================================================================================
  //����� ��� ���������
  //int keyDraw = 0;
  //double** pointDraw = new double*[NumThr * maxTrial];
  //for (int i = 0; i < NumThr * maxTrial; i++)
  //{
  //  pointDraw[i] = new double[demension];
  //  for (int j = 0; j < demension; j++)
  //    pointDraw[i][j] = 0;
  //}
  ////��������� ���� ��� ������-------------------------------------------------------------------
  //std::string pointLogName = argv[3];
  //FILE* pointLog = fopen(pointLogName.c_str(), "w");
  ////=======================================================================================================================================================

  double* z_l = new double[NumThr];
  double* z_r = new double[NumThr];
  double* z = new double[NumThr];
  for (int i = 0; i < NumThr; i++)
  {
    z_l[i] = 0;
    z_r[i] = 0;
    z[i] = 0;
  }

  Interval inter(maxTrial);

  //��� ���� ��������� �������� �� ����� ������� � ����� � ������� [0; 1] � �� �� ����� � ������ ��������
  double* zip_x_r = new double[NumThr];
  double* zip_x_l = new double[NumThr];
  double* zip_x = new double[NumThr];
  double* prom_min = new double[NumThr];
  for (int i = 0; i < NumThr; i++)
  {
    zip_x_l[i] = 0;
    zip_x_r[i] = 1;
    zip_x[i] = 0;
    prom_min[i] = 0;
  }

  //��� �������� ��������
  double min;
  GetImage(zip_x_l[0], X[0], lower_bounds, upper_bounds, demension, 10);
  min = problem->CalculateFunctionals(X[0], 0);
  z_l[0] = min;

  for (int i = 0; i < demension; i++)
  {
    (*BestX)[i] = X[0][i];
    //pointDraw[0][i] = X[0][i];
  }

  GetImage(zip_x_r[0], X[0], lower_bounds, upper_bounds, demension, 10);
  z_r[0] = problem->CalculateFunctionals(X[0], 0);

  //��������, ��� �� �����, � �� ������
  //if (z_r[0] < min)
  //{
  //  min = z_r[0];
  //  for (int i = 0; i < demension; i++)
  //  {
  //    (*BestX)[i] = X[0][i];
  //    //pointDraw[0][i] = X[0][i];
  //  }
  //}

  double M = 1;
  
  // ���� �� ���������� M
  int f = 0;

  double delta;
  double second;
  double third;
  vector<pair<double, double> > R(NumThr);
  double p;

  // �������� ��� ������������ ������
  int paral;

  //��� ����, ����� ���������� �������������� � ������������ ������
  bool test = false;

  //���������� ��������� NumThr-1 ���������� �����
  for (int i = 0; i < NumThr; i++)
  {
    double h = 1.0 / (NumThr + 1);
    zip_x[0] = (i + 1) * h;


    GetImage(zip_x[0], X[0], lower_bounds, upper_bounds, demension, 10);
    prom_min[0] = problem->CalculateFunctionals(X[0], 0);
    if (prom_min[0] < min) {
      min = prom_min[0];
      for (int kk = 0; kk < demension; kk++)
      {
        (*BestX)[kk] = X[0][kk];
        //pointDraw[i][kk] = X[0][kk];
      }
    }
  }

  for (int i = 0; i < NumThr; i++)
  {
    if (inter.GetSize() != 0)
    {
      zip_x_l[0] = inter.GetMaxPointLeft();
      zip_x_r[0] = inter.GetMaxPointRight();
      z_l[0] = inter.GetMaxZLeft();
      z_r[0] = inter.GetMaxZRight();

      inter.GetMax();
    }

    /*if (zip_x_l[0] != 0 && zip_x_r[0] != 1) {
      second = (zip_x_r[0] + zip_x_l[0]) / 2;
      p = pow((fabs(z_r[0] - z_l[0]) / M), demension);
      third = p / 2 / r;
      zip_x[0] = second - sng(z_r[0] - z_l[0]) * third;
    }
    else
      zip_x[0] = (zip_x_r[0] + zip_x_l[0]) / 2;*/

    ///����� ������ ����� ������ � ������������ �����
    double h = 1.0 / (NumThr + 1);
     zip_x[0] = (i + 1) * h;


    GetImage(zip_x[0], X[0], lower_bounds, upper_bounds, demension, 10);
    prom_min[0] = problem->CalculateFunctionals(X[0], 0);
    z[0] = prom_min[0];

    (*iteration)++;

    if ((zip_x[0] >= zip_x_r[0]) || (zip_x[0] <= zip_x_l[0]))
      std::cout << "ERROR!!!!!!!\n";

    if (flag == 0) {
      p = pow(fabs(zip_x[0] - zip_x_r[0]), N);
      if (p < accuracy)
        return min;
      p = pow(fabs(zip_x[0] - zip_x_l[0]), N);
      if (p < accuracy)
        return min;
    }
    else {
      //double maxfabs = fabs(Truemin[0] - (*BestX)[0]);
      //if (fabs(Truemin[1] - (*BestX)[1]) > maxfabs)
      //  maxfabs = fabs(Truemin[1] - (*BestX)[1]);
      //cout << i << "; " << maxfabs << "\n";

      int br = 0;
      for (int kk = 0; kk < demension; kk++)
        if (fabs(Truemin[kk] - (*BestX)[kk]) < accuracy)
          br++;
      if (br == demension)
        return min;
    }

    //��������� ������ M ��� ����������� ��������� �������
    double m = 1;
    p = pow(fabs(zip_x[0] - zip_x_r[0]), N);
    m = fabs(z[0] - z_r[0]) / p;

    if (m > M)
    {
      M = m;
      f = 1;
    }
    p = pow(fabs(zip_x[0] - zip_x_l[0]), N);
    m = fabs(z[0] - z_l[0]) / p;
    if (m > M)
    {
      M = m;
      f = 1;
    }

    if (i == 0) {
      M = 1;
      f = 1;
    }

    ////////////////////////////////////////cout << "M = " << M << "\n";

    //�������� �������������� R(i)
    R[0].first = DBL_MIN;


    if (zip_x_l[0] == 0)
    {
      delta = pow((zip_x[0] - zip_x_l[0]), N);
      R[0].first = 2 * delta - 4 * (z[0] - min) / (r*M);
    }
    else
    {
      delta = pow((zip_x[0] - zip_x_l[0]), N);
      second = (z[0] - z_l[0])*(z[0] - z_l[0]) / (r*r*M*M*delta);
      third = 2 * (z[0] + z_l[0] - 2 * min) / (r*M);
      R[0].first = delta + second - third;
    }

    inter.AddElem(R[0].first, zip_x_l[0], zip_x[0], z_l[0], z[0]);

    //-----------------------------------------------
    R[0].second = DBL_MIN;

    if (zip_x_r[0] == 1)
    {
      delta = pow((zip_x_r[0] - zip_x[0]), N);
      R[0].second = 2 * delta - 4 * (z[0] - min) / (r*M);
    }
    else
    {
      delta = pow((zip_x_r[0] - zip_x[0]), N);
      second = (z_r[0] - z[0])*(z_r[0] - z[0]) / (r*r*M*M*delta);
      third = 2 * (z_r[0] + z[0] - 2 * min) / (r*M);
      R[0].second = delta + second - third;
    }

    inter.AddElem(R[0].second, zip_x[0], zip_x_r[0], z[0], z_r[0]);

    if (f == 1)
    {
      inter.Update(M, r, demension, min);
      f = 0;
    }
  }   

  //keyDraw = NumThr - 1;

  //���� ������ �������� ��� ����� ����������� ��������
  for (int i = NumThr; i < maxTrial; i++)
  {
    if (test == true)
      break;

    for (int j = 0; j < NumThr; j++)
      if (inter.GetSize() != 0)
      {
        zip_x_l[j] = inter.GetMaxPointLeft();
        zip_x_r[j] = inter.GetMaxPointRight();
        z_l[j] = inter.GetMaxZLeft();
        z_r[j] = inter.GetMaxZRight();

        inter.GetMax();
      }

//#pragma omp parallel for private(paral, second, p, third) num_threads(NumThr)
    for (paral = 0; paral < NumThr; paral++)
    {
      if (zip_x_l[paral] != 0 && zip_x_r[paral] != 1) {
        second = (zip_x_r[paral] + zip_x_l[paral]) / 2;
        p = pow((fabs(z_r[paral] - z_l[paral]) / M), demension);
        third = p / 2 / r;
        zip_x[paral] = second - sng(z_r[paral] - z_l[paral]) * third;
      }
      else
        zip_x[paral] = (zip_x_r[paral] + zip_x_l[paral]) / 2;

      if (zip_x_r[paral] - zip_x_l[paral] == 0)
        cout << "?????????????????????????????????????????????????????????????\n";
      GetImage(zip_x[paral], X[paral], lower_bounds, upper_bounds, demension, 10);
      prom_min[paral] = problem->CalculateFunctionals(X[paral], 0);
      z[paral] = prom_min[paral];
    }

    for (int j = 0; j < NumThr; j++)
    {
     // for (int kk = 0; kk < demension; kk++)
        //pointDraw[keyDraw][kk] = X[j][kk];
      //keyDraw++;
      if (prom_min[j] < min) {
        min = prom_min[j];
        for (int kk = 0; kk < demension; kk++)
        {
          (*BestX)[kk] = X[j][kk];
        }
      }
    }
    /*double maxfabs = fabs(Truemin[0] - (*BestX)[0]);
    if (fabs(Truemin[1] - (*BestX)[1]) > maxfabs)
      maxfabs = fabs(Truemin[1] - (*BestX)[1]);
    cout << i << "; " << maxfabs << "\n";*/

    (*iteration)++;

//#pragma omp parallel for shared(test, M, f) private(paral, second, p, third, delta) num_threads(NumThr)
    for (paral = 0; paral < NumThr; paral++)
    {
      if ((zip_x[paral] >= zip_x_r[paral]) || (zip_x[paral] <= zip_x_l[paral]))
        std::cout << "ERROR!!!!!!!\n";

      if (flag == 0) {
        p = pow(fabs(zip_x[paral] - zip_x_r[paral]), N);
        if (p < accuracy)
          test = true;
        p = pow(fabs(zip_x[paral] - zip_x_l[paral]), N);
        if (p < accuracy)
          test = true;
      }
      else {
        int br = 0;
        for (int kk = 0; kk < demension; kk++)
          if (fabs(Truemin[kk] - (*BestX)[kk]) < accuracy)
            br++;
        if (br == demension)
          test = true;
      }

      if (test != true)
      {
        //��������� ������ M ��� ����������� ��������� �������
        double m = 1;
        p = pow(fabs(zip_x[paral] - zip_x_r[paral]), N);
        m = fabs(z[paral] - z_r[paral]) / p;
        if (m > M)
        {
          M = m;
          f = 1;
        }
        p = pow(fabs(zip_x[paral] - zip_x_l[paral]), N);
        m = fabs(z[paral] - z_l[paral]) / p;
        if (m > M)
        {
          M = m;
          f = 1;
        }

        //�������� �������������� R(i)
        R[paral].first = DBL_MIN;

        if (zip_x_l[paral] == 0)
        {
          delta = pow((zip_x[paral] - zip_x_l[paral]), N);
          R[paral].first = 2 * delta - 4 * (z[paral] - min) / (r*M);
        }
        else
        {
          delta = pow((zip_x[paral] - zip_x_l[paral]), N);
          second = (z[paral] - z_l[paral])*(z[paral] - z_l[paral]) / (r*r*M*M*delta);
          third = 2 * (z[paral] + z_l[paral] - 2*min) / (r*M);
          R[paral].first = delta + second - third;
        }

        //-----------------------------------------------
        R[paral].second = DBL_MIN;

        if (zip_x_r[paral] == 1)
        {
          delta = pow((zip_x_r[paral] - zip_x[paral]), N);
          R[paral].second = 2 * delta - 4 * (z[paral] - min) / (r*M);
        }
        else
        {
          delta = pow((zip_x_r[paral] - zip_x[paral]), N);
          second = (z_r[paral] - z[paral])*(z_r[paral] - z[paral]) / (r*r*M*M*delta);
          third = 2 * (z_r[paral] + z[paral] - 2 * min) / (r*M);
          R[paral].second = delta + second - third;
        }
      }
    }

    //�������� �����
    /*double* testR = new double[inter.GetSize()];
    double* test_zip_x_l = new double[inter.GetSize()];
    double* test_zip_x_r = new double[inter.GetSize()];
    double* test_z_l = new double[inter.GetSize()];
    double* test_z_r = new double[inter.GetSize()];
    int size = inter.GetSize();
    for (int ii = 0; ii < size; ii++)
    {
      test_zip_x_l[ii] = inter.GetMaxPointLeft();
      cout << "zip_x_l = " << test_zip_x_l[ii] << "\t";
      test_zip_x_r[ii] = inter.GetMaxPointRight();
      cout << "zip_x_r = " << test_zip_x_r[ii] << "\t";
      test_z_l[ii] = inter.GetMaxZLeft();
      cout << "z_l = " << test_z_l[ii] << "\t";
      test_z_r[ii] = inter.GetMaxZRight();
      cout << "z_r = " << test_z_r[ii] << "\t";
      testR[ii] = inter.GetMax();
      cout << "R = " << testR[ii] << "\t";

      p = pow(fabs(test_zip_x_r[ii] - test_zip_x_l[ii]), N);
      double m = fabs(test_z_r[ii] - test_z_l[ii]) / p;
      cout << "m = " << m << "\n";
    }
    for (int ii = 0; ii < size; ii++)
    {
      for (int ii2 = 0; ii2 < size; ii2++)
      {
        if ((test_zip_x_l[ii] > test_zip_x_l[ii2]) && (test_zip_x_l[ii] < test_zip_x_r[ii2]))
          cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        if ((test_zip_x_r[ii] > test_zip_x_l[ii2]) && (test_zip_x_r[ii] < test_zip_x_r[ii2]))
          cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
      }
      inter.AddElem(testR[ii], test_zip_x_l[ii], test_zip_x_r[ii], test_z_l[ii], test_z_r[ii]);
    }
    cout << "\n";*/

    /////////////////////////////////////////////////////////////////////////////////////cout << "M = " << M << "\n";

    //������� ��� ����� ��������� � ���� � ��� ���������� ���������
    for (int j = 0; j < NumThr; j++)
    {
      inter.AddElem(R[j].first, zip_x_l[j], zip_x[j], z_l[j], z[j]);
      inter.AddElem(R[j].second, zip_x[j], zip_x_r[j], z[j], z_r[j]);
    }
    if (f == 1)
    {
      inter.Update(M, r, demension, min);
      f = 0;
    }
  }


  ////=======================================================================================================================================================
  //fprintf(pointLog, "%d\n", keyDraw);
  //////���������� � ���� �����----------------------------------------------------------------------------------
  //for (int i = 1; i < keyDraw; i++)
  //{
  //  for (int j = 0; j < demension; j++)
  //    fprintf(pointLog, "%lf ", pointDraw[i][j]);
  //  fprintf(pointLog, "\n");
  //}
  ////---------------------------------------------------------------------------------------------------------------

  /// �������� ��������� �����
  /*for (int j = 0; j < demension; j++)
    fprintf(pointLog, "%lf ", (*BestX)[j]);
  fprintf(pointLog, "\n");

  /// ��������� ����� ����������� ��������
  for (int j = 0; j < demension; j++)
    fprintf(pointLog, "%lf ", Truemin[j]);
  fclose(pointLog);*/
  ////=======================================================================================================================================================

  delete[] lower_bounds;
  delete[] upper_bounds;

  return min;
}
