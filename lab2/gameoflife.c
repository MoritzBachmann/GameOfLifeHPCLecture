#include <omp.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
// OPTIONAL: comment this out for No output
#define NO_OUTOUT
// OPTIONAL: comment this out for console output
// #define CONSOLE_OUTPUT

#define calcIndex(width, x, y) ((y) * (width) + (x))
#define ALIVE 1
#define DEAD 0

#define X 0
#define Y 1

#define START_TIMEMEASUREMENT(name)                                                               \
  struct timeval __FILE__##__func__##name##actualtime;                                            \
  gettimeofday(&__FILE__##__func__##name##actualtime, NULL);                                      \
  double __FILE__##__func__##name##s_time = (double)__FILE__##__func__##name##actualtime.tv_sec + \
                                            ((double)__FILE__##__func__##name##actualtime.tv_usec / 1000000.0)

#define END_TIMEMEASUREMENT(name, res)                        \
  gettimeofday(&__FILE__##__func__##name##actualtime, NULL);  \
  res = (double)__FILE__##__func__##name##actualtime.tv_sec + \
        ((double)__FILE__##__func__##name##actualtime.tv_usec / 1000000.0) - __FILE__##__func__##name##s_time

typedef uint8_t number_type;
typedef uint64_t header_type;
#define NUMBER_TYPE_VTK_NAME "UInt8"
#define HEADER_TYPE_VTK_NAME "UInt64"

void myexit(const char *s, ...)
{
  va_list args;
  va_start(args, s);
  vprintf(s, args);
  printf("\n");
  va_end(args);
  abort();
}

int testLittleEndian()
{
  int32_t test = 1;
  char *testdata = (char *)&test;
  if (testdata[0] == 1)
  {
    return 1;
  }
  return 0;
}

const char *vtk_header_template =
    "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"%s\" header_type=\"" HEADER_TYPE_VTK_NAME "\">\n"
    "  <ImageData WholeExtent=\"0 %d 0 %d 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 0\">\n"
    "    <Piece Extent=\"0 %d 0 %d 0 0\">\n"
    "      <PointData>\n"
    "      </PointData>\n"
    "      <CellData Scalars=\"GameOfLife\">\n"
    "        <DataArray type=\"" NUMBER_TYPE_VTK_NAME
    "\" Name=\"GameOfLife\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" offset=\"0\"/>\n"
    "      </CellData>\n"
    "    </Piece>\n"
    "  </ImageData>\n"
    "  <AppendedData encoding=\"raw\">\n"
    "   _";

char vtk_header[10000];
void create_vtk_header(char *header, int width, int height)
{
  snprintf(header, 10000, vtk_header_template, testLittleEndian() ? "LittleEndian" : "BigEndian", width, height, width,
           height);
}

char *vtk_tail = "\n  </AppendedData>\n"
                 "</VTKFile>\n";

void write_vtk_data(FILE *f, char *data, int length)
{
  if (fwrite(data, 1, length, f) != length)
  {
    myexit("Could not write vtk-Data");
  }
}

void write_field(number_type *currentfield, int width, int height, int timestep)
{
#ifdef NO_OUTOUT
  //printf("finished timestep %d\n", timestep);
  return;
#endif
#ifdef CONSOLE_OUTPUT
  printf("\033[H");
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
      printf(ALIVE == currentfield[calcIndex(width, x, y)] ? "\033[07m  \033[m" : "  ");
    printf("\033[E");
  }
  fflush(stdout);
  printf("\ntimestep=%d", timestep);
  usleep(80000);
#else
  if (timestep == 0)
  {
    mkdir("./gol/", 0777);
    create_vtk_header(vtk_header, width, height);
  }
  printf("writing timestep %d\n", timestep);
  FILE *fp; // The current file handle.
  char filename[1024];
  snprintf(filename, 1024, "./gol/gol-%05d.vti", timestep);
  fp = fopen(filename, "w");
  // write the header of the vti file
  write_vtk_data(fp, vtk_header, strlen(vtk_header));
  // write the length in bytes before writing the data
  header_type length = width * height * sizeof(number_type);
  write_vtk_data(fp, (char *)&length, sizeof(header_type));
  write_vtk_data(fp, (char *)currentfield, width * height * sizeof(number_type));
  // write the tail of the vti file
  write_vtk_data(fp, vtk_tail, strlen(vtk_tail));
  fclose(fp);
  //printf("finished writing timestep %d\n", timestep);
#endif
}

int count_neighbours(number_type *currentfield, int pos_x, int pos_y, int width)
{
  int count = 0;
  for (int i = pos_x - 1; i <= pos_x + 1; i++)
  {
    for (int j = pos_y - 1; j <= pos_y + 1; j++)
    {
      int t = currentfield[calcIndex(width, i, j)];
      if (t == ALIVE)
      {
        count++;
      }
    }
  }
  if (currentfield[calcIndex(width, pos_x, pos_y)] == ALIVE)
    count--;
  return count;
}

void evolve(number_type *currentfield, number_type *newfield, int starts[2], int ends[2], int width)
{
  // #pragma omp parallel for // collapse(2)
  for (int i = starts[0]; i < ends[0]; i++)
  {
    for (int j = starts[1]; j < ends[1]; j++)
    {
      number_type *current_cell = currentfield + calcIndex(width, i, j);
      number_type *new_cell = newfield + calcIndex(width, i, j);
      int neighbours = count_neighbours(currentfield, i, j, width);
      if (*current_cell == ALIVE)
      {
        if (neighbours < 2 || neighbours > 3)
          *new_cell = DEAD;
        else
          *new_cell = ALIVE;
      }
      else if (*current_cell == DEAD) // DEAD
      {
        if (neighbours == 3)
          *new_cell = ALIVE;
        else
          *new_cell = DEAD;
      }
      else
      {
        printf("Warn unexpected cel value \n");
      }
    }
  }
  // void evolve(number_type* currentfield, number_type* newfield, int width, int height) {
  // TODO traverse through each voxel and implement game of live logic and
  // parallelize using OpenMP.
  // HINT: use 'starts' and 'ends'
}

void filling_random(number_type *currentfield, int width, int height)
{
  int i;
  for (int y = 1; y < height - 1; y++)
  {
    for (int x = 1; x < width - 1; x++)
    {
      i = calcIndex(width, x, y);
      currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
    }
  }
}

void filling_runner(number_type *currentfield, int width, int height)
{
  int offset_x = width / 3;
  int offset_y = height / 2;
  currentfield[calcIndex(width, offset_x + 0, offset_y + 1)] = ALIVE;
  currentfield[calcIndex(width, offset_x + 1, offset_y + 2)] = ALIVE;
  currentfield[calcIndex(width, offset_x + 2, offset_y + 0)] = ALIVE;
  currentfield[calcIndex(width, offset_x + 2, offset_y + 1)] = ALIVE;
  currentfield[calcIndex(width, offset_x + 2, offset_y + 2)] = ALIVE;
}

void apply_periodic_boundaries(number_type *field, int width, int height)
{
  for (size_t i = 1; i < width - 1; i++)
  {
    field[calcIndex(width, i, 0)] = field[calcIndex(width, i, height - 2)];
    field[calcIndex(width, i, height - 1)] = field[calcIndex(width, i, 1)];
  }
  for (size_t i = 0; i < height; i++)
  {
    field[calcIndex(width, 0, i)] = field[calcIndex(width, width - 2, i)];
    field[calcIndex(width, width - 1, i)] = field[calcIndex(width, 1, i)];
  }
}

void game(int width, int height, int num_timesteps, int *decomposition)
{
  (void)decomposition; // required for task b suppress warning in task a
  number_type *currentfield = calloc(width * height, sizeof(number_type));
  number_type *newfield = calloc(width * height, sizeof(number_type));

  // TODO 1: use your favorite filling
  filling_random (currentfield, width, height);
  filling_runner(currentfield, width, height);
  int delta_height = (height - 2) / decomposition[Y];
  int delta_width = (width - 2) / decomposition[X];
  int starts[2];
  int ends[2];
  int time = 0;
  int size = decomposition[X] *decomposition[Y];
  int startsX[size];
  int startsY[size];
  int count = 0;
  write_field(currentfield, width, height, time);
  // TODO 3: implement periodic boundary condition
  apply_periodic_boundaries(currentfield, width, height);

for (int i = 1; i < height-delta_height; i= i+delta_height){
  for (int j = 1; j < width-delta_width; j= j+delta_width)
  {
    //printf("count: %d startY %d startX %d \n",count, i,j);
    startsY[count] = i;
    startsX[count] = j;
    count++;
  }
}

  for (time = 1; time <= num_timesteps; time++)
  {
    #pragma omp parallel num_threads(size) private(starts, ends)
    {
      starts[X] = startsX[omp_get_thread_num()];
      starts[Y] = startsY[omp_get_thread_num()];
      //printf("start thread: %d  X: %d Y: %d \n",omp_get_thread_num(), starts[X], starts[Y]);
      ends[X] = starts[X]+delta_width;
      ends[Y] = starts[Y]+delta_height;
      //printf("end thread: %d  X: %d Y: %d \n",omp_get_thread_num(), ends[X], ends[Y]);
      evolve(currentfield, newfield, starts, ends, width);
    }
  // TODO 3: implement periodic boundary condition
  apply_periodic_boundaries(newfield, width, height);
  write_field(newfield, width, height, time);
  // TODO 4: implement SWAP of the fields
  number_type *temp = currentfield;
  currentfield = newfield;
  newfield = temp;
}
free(currentfield);
free(newfield);
}

int main(int c, char **v)
{
/*
#pragma omp parallel
  {
    if (omp_get_thread_num() == 0)
    {
      printf("Running with %d threads\n", omp_get_num_threads());
    }
  }
  */

  int width, height, num_timesteps;
  int decomposition[2] = {1, 1};
  if (c == 4 || c == 6)
  {
    width = atoi(v[1]) + 2;     ///< read width + 2 boundary cells (low x, high x)
    height = atoi(v[2]) + 2;    ///< read height + 2 boundary cells (low y, high y)
    num_timesteps = atoi(v[3]); ///< read timesteps

    if (width <= 0)
    {
      width = 32; ///< default width
    }
    if (height <= 0)
    {
      height = 32; ///< default height
    }
    if (c == 6)
    {
      decomposition[X] = atoi(v[4]); ///< read number of threads in x
      decomposition[Y] = atoi(v[5]); ///< read number of threads in y
      if ((width-2) % decomposition[X] != 0 || (height-2) % decomposition[Y] != 0){
        myexit("The field must be divideble by the decomposition");
      }
    }
  }
  else
  {
    myexit("Too less arguments");
  }
  double elapsed_time;
  START_TIMEMEASUREMENT(measure_game_time);

  game(width, height, num_timesteps, decomposition);

  END_TIMEMEASUREMENT(measure_game_time, elapsed_time);
  printf("time elapsed: %lf sec\n", elapsed_time);
}
