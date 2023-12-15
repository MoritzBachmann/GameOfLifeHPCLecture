#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdbool.h>
#include <stdint.h>
#include <mpi.h>
// OPTIONAL: comment this out for console output
// #define CONSOLE_OUTPUT

#define calcIndex(width, x, y) ((y) * (width) + (x))
#define ALIVE 1
#define DEAD 0

#define START_TIMEMEASUREMENT(name)                                                                                         \
  struct timeval __FILE__##__func__##name##actualtime;                                                                      \
  gettimeofday(&__FILE__##__func__##name##actualtime, NULL);                                                                \
  double __FILE__##__func__##name##s_time = (double)__FILE__##__func__##name##actualtime.tv_sec +                           \
                                            ((double)__FILE__##__func__##name##actualtime.tv_usec / 1000000.0)

#define END_TIMEMEASUREMENT(name, res)                                                                                      \
  gettimeofday(&__FILE__##__func__##name##actualtime, NULL);                                                                \
  res = (double)__FILE__##__func__##name##actualtime.tv_sec +                                                               \
        ((double)__FILE__##__func__##name##actualtime.tv_usec / 1000000.0) - __FILE__##__func__##name##s_time

typedef uint8_t number_type;
typedef uint64_t header_type;
#define NUMBER_TYPE_VTK_NAME "UInt8"
#define HEADER_TYPE_VTK_NAME "UInt64"

void myexit(const char *s, ...) {
  va_list args;
  va_start(args, s);
  vprintf(s, args);
  printf("\n");
  va_end(args);
  abort();
}

int testLittleEndian() {
  int32_t test = 1;
  char *testdata = (char *)&test;
  if (testdata[0] == 1) {
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
void create_vtk_header(char *header, int width, int height) {
  snprintf(header, 10000, vtk_header_template, testLittleEndian() ? "LittleEndian" : "BigEndian", width, height, width,
           height);
}

char *vtk_tail = "\n  </AppendedData>\n"
                 "</VTKFile>\n";

void write_vtk_data(FILE *f, char *data, int length) {
  if (fwrite(data, 1, length, f) != length) {
    myexit("Could not write vtk-Data");
  }
}

void write_field(number_type *currentfield, int width, int height, int timestep) {
#ifdef CONSOLE_OUTPUT
  printf("\033[H");
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++)
      printf(ALIVE == currentfield[calcIndex(width, x, y)] ? "\033[07m  \033[m" : "  ");
    printf("\033[E");
  }
  fflush(stdout);
  printf("\ntimestep=%d", timestep);
  usleep(80000);
#else
  if (timestep == 0) {
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
  printf("finished writing timestep %d\n", timestep);
#endif
}


void evolve(number_type *currentfield, number_type *newfield, int width, int height) {
  // #pragma omp parallel for // collapse(2)
  for (int i = start_indices[0]; i < start_indices[0]+width; i++)
  {
    for (int j = start_indices[1]; j < start_indices[1]+height; j++)
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
}

void filling_random(number_type *currentfield, int width, int height) {
  int i;
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      i = calcIndex(width, x, y);
      currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
    }
  }
}

void filling_runner(number_type *currentfield, int width, int height) {
  int offset_x = width / 3;
  int offset_y = height / 2;
  currentfield[calcIndex(width, offset_x + 0, offset_y + 1)] = ALIVE;
  currentfield[calcIndex(width, offset_x + 1, offset_y + 2)] = ALIVE;
  currentfield[calcIndex(width, offset_x + 2, offset_y + 0)] = ALIVE;
  currentfield[calcIndex(width, offset_x + 2, offset_y + 1)] = ALIVE;
  currentfield[calcIndex(width, offset_x + 2, offset_y + 2)] = ALIVE;
}

void apply_periodic_boundaries(number_type *field, int width, int height) {

 

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

void game(int width, int height, int num_timesteps, int rank) {
  
  number_type *local_new_field = calloc(width * height, sizeof(number_type));
  // TODO 1: use your favorite filling
  // filling_random (currentfield, width, height);
  filling_runner(local_field, width, height);

  int time = 0;
  write_field(local_field, width, height, time);
  // TODO 4: implement periodic boundary condition
  int left, bottom, right, top;

  MPI_Cart_shift(MPI_COMM_WORLD, 0, -1, *rank, *left );
  MPI_Cart_shift(MPI_COMM_WORLD, 0, 1, *rank, *right );
  MPI_Cart_shift(MPI_COMM_WORLD, 1, -1, *rank, *left );
  MPI_Cart_shift(MPI_COMM_WORLD, 1, 1, *rank, *right );

  int row[2],col[2];
  row[0] =[width];
  row[1] = 0;
  col[0] =[0];
  col[1] = height; 
  MPI_Type_create_subarray(1, lsizes,row , [0][0]],
                         MPI_ORDER_C, MPI_FLOAT, &left_border);
  MPI_Type_commit(&left_Border);
  MPI_Type_create_subarray(1, lsizes,row , [width-1][0]],
                         MPI_ORDER_C, MPI_FLOAT, &right_border);
  MPI_Type_commit(&left_Border);
  MPI_Type_create_subarray(1, lsizes,col , [0][0]],
                         MPI_ORDER_C, MPI_FLOAT, &bottom_border);
  MPI_Type_commit(&left_Border);
  MPI_Type_create_subarray(1, lsizes,col , [0][hight-1]],
                         MPI_ORDER_C, MPI_FLOAT, &top_border);
  MPI_Type_commit(&left_Border);
  
  apply_periodic_boundaries(currentfield, width, height);

  for (time = 1; time <= num_timesteps; time++) {
    // TODO 2: implement evolve function (see above)
    evolve(local_field, newfield, width, height);
    
    //write_field(newfield, width, height, time);
    // TODO 4: implement periodic boundary condition
    apply_periodic_boundaries(newfield, width, height);
    // TODO 3: implement SWAP of the fields
  number_type *temp = currentfield;
  currentfield = newfield;
  newfield = temp;
  }

  free(currentfield);
  free(newfield);
}

int main(int c, char **v) {
  // TODO 5: implement MPI 
  
  MPI_Init(&c, &v);
  MPI_Comm cart_comm;
  MPI_Datatype filetype, memtype;

  int width, height, num_timesteps;
  int process_numX;
  int process_numY;
  int rank, size, i;
  if (c == 6) {
    width = atoi (v[1]); ///< read width + 2 boundary cells (low x, high x)
    height = atoi (v[2]); ///< read height + 2 boundary cells (low y, high y)
    num_timesteps = atoi (v[3]); ///< read timesteps
    
    if (width <= 0) {
      width = 32; ///< default width
    }
    if (height <= 0) {
      height = 32; ///< default height
    }
    process_numX = atoi (v[4]); ///< read number of processes in X
    process_numY = atoi (v[5]); ///< read number of processes in Y
    
  }
  else {
   myexit("Too less arguments");
  }
  
  // TODO 5a: get the global rank of the process and save it to rank_global
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (size != process_numX * process_numY) {
  printf("Communicator size must be %d \n",(c * process_numY));
  MPI_Abort(MPI_COMM_WORLD, 1);
  }
  // TODO 5b: get the number of processes and save it to num_tasks variable
  
  // TODO 5c: create a new cartesian communicator of the worker communicator and get the information.
  
  int gsizes[2] = {width, height};  // global size of the domain without boundaries
  int lsizes[2];
  
  lsize[0] = (width + process_numX - 1) / process_numX;
  lsize[1] =  (height + process_numY - 1) / process_numY;

  int dims[0] = process_numX;
  int dims[1] = process_numY;
  int periods[0] = periods[1] = 1;
  int coords[2];
  int start_indices[2];
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
  MPI_Comm_rank(cart_comm, &rank);
  MPI_Cart_coords(cart_comm, rank, 2, coords);

  start_indices[0] = coords[0] * lsizes[0];
  start_indices[1] = coords[1] * lsizes[1];
    /* TODO 5d: create and commit a subarray as a new filetype to describe the local
   *      worker field as a part of the global field.
   *      Use the global variable 'filetype'.
   * HINT: use MPI_Type_create_subarray and MPI_Type_commit functions
   */

  
  MPI_Type_create_subarray(2, gsizes, lsizes, start_indices,
                         MPI_ORDER_C, MPI_FLOAT, &filetype);
  MPI_Type_commit(&filetype);


   
  /* TODO 5e: Create a derived datatype that describes the layout of the inner local field
   *      in the memory buffer that includes the ghost layer (local field).
   *      This is another subarray datatype!
   *      Use the global variable 'memtype'.
  */
 int memsizes[2];
  memsizes[0] = lsizes[0] + 2; /* no. of rows in allocated array */
  memsizes[1] = lsizes[1] + 2; /* no. of columns in allocated array */
  intlocal_start_indices[2];
  local_start_indices[0] = local_start_indices[1] = 1; // one ghost layer on each side

  number_type *local_field = calloc(memsizes[0] * memsizes[1], sizeof(number_type));
  MPI_Type_create_subarray(2, memsizes, lsizes, start_indices,
                           MPI_ORDER_C, MPI_FLOAT, &memtype);
  
  game(lsizes[X], lsizes[Y], num_timesteps, gsizes, rank);
  
  MPI_Finalize();
}
