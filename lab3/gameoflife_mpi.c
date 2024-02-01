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
  printf("finished writing timestep %d\n", timestep);
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

void send_mpi(number_type *source, number_type *destination, int local_starts, int rows, int colums, int num_procs, int rank)
{
  for (int i = 1; i < colums; i++)
  {
    // printf("\n");
    for (int j = 1; j < rows; j += 1)
    {
      int local_pos = i * colums + j;
      int global_pos = local_starts + local_pos;
      if (rank == 0)
      {
        for (int sender = 1; sender < num_procs; sender++)
        {
          //printf("iD: %d - reciving %d from rank %d \n", rank, global_pos, sender);
          MPI_Status rec;
          MPI_Recv((destination + global_pos), 1, MPI_INT, sender, 99, MPI_COMM_WORLD, &rec);
        }
      }
      else
      {
        //printf("iD: %d - sendig location %d from rank \n", rank, local_pos);
        MPI_Send((source + local_pos), 1, MPI_INT, 0, 99, MPI_COMM_WORLD);
      }
    }
  }
}

void evolve(number_type *currentfield, number_type *newfield, int width, int height, int *start_indices)
{
  // #pragma omp parallel for // collapse(2)
  for (int i = start_indices[0]; i < start_indices[0] + width; i++)
  {
    for (int j = start_indices[1]; j < start_indices[1] + height; j++)
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
        printf("Warn unexpected cel value: %d at %d %d\n", *current_cell, i,j);
      }
    }
  }
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

void apply_periodic_boundaries(number_type *field, int width, int height, int *start_indices)
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

void game(int c, int process_numX, int process_numY, int width, int height, int num_timesteps)
{
  MPI_Comm cart_comm;
  MPI_File fh;
  MPI_Status status;
  MPI_Datatype filetype, memtype, left_border, right_border, bottom_border, top_border;
  int rank, size, i;
  // TODO 5a: get the global rank of the process and save it to rank_global
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (size != process_numX * process_numY)
  {
    printf("Communicator size is %d but must be %d \n", size, (process_numX * process_numY));
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  // TODO 5b: get the number of processes and save it to num_tasks variable

  // TODO 5c: create a new cartesian communicator of the worker communicator and get the information.

  int gsizes[2] = {width, height}; // global size of the domain without boundaries
  int lsizes[2];

  lsizes[0] = (width + process_numX - 1) / process_numX;
  lsizes[1] = (height + process_numY - 1) / process_numY;

  int dims[2];
  dims[0] = process_numX;
  dims[1] = process_numY;
  int periods[2];
  periods[0] = periods[1] = 1;
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
  if (rank == 0)
  {
  }
  MPI_Type_create_subarray(2, gsizes, lsizes, start_indices,
                           MPI_ORDER_C, MPI_INT, &filetype);
  MPI_Type_commit(&filetype);

  /* TODO 5e: Create a derived datatype that describes the layout of the inner local field
   *      in the memory buffer that includes the ghost layer (local field).
   *      This is another subarray datatype!
   *      Use the global variable 'memtype'.
   */
  int memsizes[2];
  memsizes[0] = lsizes[0] + 2; /* no. of rows in allocated array */
  memsizes[1] = lsizes[1] + 2; /* no. of columns in allocated array */
  int local_start_indices[2];
  local_start_indices[0] = local_start_indices[1] = 1; // one ghost layer on each side

  number_type *local_field = calloc(memsizes[0] * memsizes[1], sizeof(number_type));
  number_type *newfield = calloc(memsizes[0] * memsizes[1], sizeof(number_type));
  // probably i do not need this
  /*
  MPI_Type_create_subarray(2, memsizes, lsizes, start_indices,
                           MPI_ORDER_C, MPI_INT, &memtype);
                           */

  // TODO  use your favorite filling
  // filling_random (currentfield, width, height);
  filling_runner(local_field, lsizes[0], lsizes[1]);

  int time = 0;
  // write_field(local_field, width, height, time);
  //  TODO 4: implement periodic boundary condition

  enum DIRECTIONS
  {
    DOWN,
    UP,
    LEFT,
    RIGHT
  };
  char *neighbours_names[4] = {"down", "up", "left", "right"};
  int neighbours_ranks[4];

  // Let consider dims[0] = X, so the shift tells us our left and right neighbours
  MPI_Cart_shift(cart_comm, 0, 1, &neighbours_ranks[LEFT], &neighbours_ranks[RIGHT]);

  // Let consider dims[1] = Y, so the shift tells us our up and down neighbours
  MPI_Cart_shift(cart_comm, 1, 1, &neighbours_ranks[DOWN], &neighbours_ranks[UP]);
  /*
  MPI_Cart_shift(cart_comm, 0, -1, *left, *corner_b_l);
  MPI_Cart_shift(cart_comm, 0, 1, *left, *corner_t_l);
  MPI_Cart_shift(cart_comm, 1, -1, *right, *corner_b_r);
  MPI_Cart_shift(cart_comm, 1, 1, *right, *corner_t_r );
  */
  int row[2], col[2];
  row[0] = memsizes[0];
  row[1] = 1;
  col[0] = 1;
  col[1] = memsizes[1];
  int starts[2] = {0, 1};
  // printf("0. iD: %d war hier\n", rank);
  // printf("Starts:  %d  %d row: %d %d \n", starts[0], starts[1], row[0], row[1]);
  MPI_Type_create_subarray(2, memsizes, row, starts,
                           MPI_ORDER_C, MPI_INT, &bottom_border);
  MPI_Type_commit(&bottom_border);

  // printf(" 1. iD: %d war hier\n", rank);
  starts[0] = 1;
  starts[1] = 0;
  // printf("Starts:  %d  %d col: %d %d \n", starts[0], starts[1], col[0], col[1]);
  MPI_Type_create_subarray(2, memsizes, col, starts,
                           MPI_ORDER_C, MPI_INT, &left_border);
  MPI_Type_commit(&left_border);

  // printf(" 2. iD: %d war hier\n", rank);
  starts[1] = memsizes[1] - 1;
  starts[0] = 0;
  // printf("Starts:  %d  %d row: %d %d \n", starts[0], starts[1], row[0], row[1]);
  MPI_Type_create_subarray(2, memsizes, row, starts,
                           MPI_ORDER_C, MPI_INT, &top_border);
  MPI_Type_commit(&top_border);

  // printf("3. iD: %d war hier\n", rank);
  starts[0] = memsizes[1] - 1;
  starts[1] = 0;
  // printf("Starts:  %d  %d col: %d %d \n", starts[0], starts[1], col[0], col[1]);
  MPI_Type_create_subarray(2, memsizes, col, starts,
                           MPI_ORDER_C, MPI_INT, &right_border);
  MPI_Type_commit(&right_border);
  // printf("4. iD: %d war hier\n", rank);
  //MPI_Sendrecv(memsizes, 1, left_border, neighbours_ranks[LEFT], 1, memsizes, 1, right_border, neighbours_ranks[RIGHT], 1, cart_comm, MPI_STATUS_IGNORE);
  //MPI_Sendrecv(memsizes, 1, right_border, neighbours_ranks[RIGHT], 1, memsizes, 1, left_border, neighbours_ranks[LEFT], 1, cart_comm, MPI_STATUS_IGNORE);
  //MPI_Sendrecv(memsizes, 1, bottom_border, neighbours_ranks[DOWN], 1, memsizes, 1, top_border, neighbours_ranks[UP], 1, cart_comm, MPI_STATUS_IGNORE);
  //MPI_Sendrecv(memsizes, 1, top_border, neighbours_ranks[UP], 1, memsizes, 1, bottom_border, neighbours_ranks[DOWN], 1, cart_comm, MPI_STATUS_IGNORE);

  // apply_periodic_boundaries(currentfield, width, height);

  number_type *global_field = NULL;
  if (rank == 0)
  {
    global_field = calloc(gsizes[0] * gsizes[1], sizeof(number_type));
  }

  int local_start = start_indices[0] + start_indices[1] * width;

  send_mpi(local_field, global_field, local_start, lsizes[0], lsizes[1], process_numX * process_numY, rank);
  MPI_Barrier(MPI_COMM_WORLD);
  // MPI_Gather(memsizes,1,filetype,global_field,1, filetype, 0,cart_comm);
  printf("iD: %d - gather done\n", rank);

  for (time = 0; time <= num_timesteps; time++)
  {
    if (rank == 0)
    {
      write_field(global_field, width, height, time);

      printf("iD: %d - write done\n", rank);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // TODO  implement evolve function (see above)
    evolve(local_field, newfield, lsizes[0], lsizes[1], start_indices);

    // write_field(newfield, width, height, time);
    //  TODO 4: implement periodic boundary condition
    // apply_periodic_boundaries(newfield, width, height);
    // TODO 3: implement SWAP of the fields
    number_type *temp = local_field;
    local_field = newfield;
    newfield = temp;
  }

  free(local_field);
  free(newfield);
}

void _write(int *newfield, int rank, int width, int height, int time)
{
  // instead use MPI_File_write_all this
  /*
  if (myrank != 0)
  MPI_Send(buf, BUFSIZE, MPI_INT, 0, 99, MPI_COMM_WORLD);
  else {
    FILE *myfile;
    myfile = fopen("testfile", "w");
    fwrite(buf, sizeof(int), BUFSIZE, myfile);
    for (i = 1 ; i < numprocs; i++) {
    MPI_Recv(buf, BUFSIZE, MPI_INT, i, 99, MPI_COMM_WORLD, &status);
    fwrite(buf, sizeof(int), BUFSIZE, myfile);
    }
    fclose(myfile);
  }
  */
}

int main(int c, char **v)
{
  // TODO 5: implement MPI

  MPI_Init(&c, &v);

  int width, height, num_timesteps;
  int process_numX;
  int process_numY;

  if (c == 6)
  {
    width = atoi(v[1]);         ///< read width + 2 boundary cells (low x, high x)
    height = atoi(v[2]);        ///< read height + 2 boundary cells (low y, high y)
    num_timesteps = atoi(v[3]); ///< read timesteps

    if (width <= 0)
    {
      width = 32; ///< default width
    }
    if (height <= 0)
    {
      height = 32; ///< default height
    }
    process_numX = atoi(v[4]); ///< read number of processes in X
    process_numY = atoi(v[5]); ///< read number of processes in Y
  }
  else
  {
    myexit("Too less arguments");
  }

  game(c, process_numX, process_numY, width, height, num_timesteps);

  MPI_Finalize();
}
