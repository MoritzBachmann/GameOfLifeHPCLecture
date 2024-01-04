#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <immintrin.h>
#include <emmintrin.h>
#include <avxintrin.h>
#include <smmintrin.h>
// OPTIONAL: comment this out for console output
// #define CONSOLE_OUTPUT

#define calcIndex(width, x, y) ((y) * (width) + (x))
#define ALIVE 1
#define DEAD 0
#define CELLS 4 	// NUMBER OF CELLS TO VECTORIZE

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
typedef __m256i SIMD_TYPE;		// USE __m128i IF YOUR CPU DOES NOT SUPPORT AVX
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

SIMD_TYPE load_simd_vector(number_type* field_pointer) {
	return  _mm256_loadu_epi8(field_pointer);
}

void store_simd_vector(SIMD_TYPE source, number_type* destination) {
	_mm256_storeu_epi8(destination, source);
}
SIMD_TYPE add_simd_vector(SIMD_TYPE a, SIMD_TYPE b) {
	 return _mm256_adds_epu8(a, b);
}

SIMD_TYPE count_neibours(number_type** neighbor_pointer) {
  SIMD_TYPE n0 = load_simd_vector(neighbor_pointer[0]);
  SIMD_TYPE n1 = load_simd_vector(neighbor_pointer[1]);
  SIMD_TYPE n2 = load_simd_vector(neighbor_pointer[2]);
  SIMD_TYPE n3 = load_simd_vector(neighbor_pointer[3]);
  SIMD_TYPE n4 = load_simd_vector(neighbor_pointer[4]);
  SIMD_TYPE n5 = load_simd_vector(neighbor_pointer[5]);
  SIMD_TYPE n6 = load_simd_vector(neighbor_pointer[6]);
  SIMD_TYPE n7 = load_simd_vector(neighbor_pointer[7]);
  SIMD_TYPE n = _mm256_adds_epu8(n0, n1);
  n =_mm256_adds_epu8(n, n2);
  n =_mm256_adds_epu8(n, n3);
  n =_mm256_adds_epu8(n, n4);
  n =_mm256_adds_epu8(n, n5);
  n =_mm256_adds_epu8(n, n6);
  n =_mm256_adds_epu8(n, n7);
  return n;
}

void evolve(number_type *currentfield, number_type *newfield, int width, int height, number_type* source_pointer, number_type* destination_pointer, number_type** neighbor_pointer) {
  for (size_t i = 0; i < height; i++) {
    for (size_t j = 0; j < width; j += sizeof(SIMD_TYPE)) {
  source_pointer = currentfield + i*width+j;
  destination_pointer = currentfield + i*width+j;
  neighbor_pointer[0] = source_pointer -1;
  neighbor_pointer[1] = source_pointer +1;
  neighbor_pointer[2] = source_pointer - width;
  neighbor_pointer[3] = source_pointer - width-1;
  neighbor_pointer[4] = source_pointer - width+1;
  neighbor_pointer[5] = source_pointer + width;
  neighbor_pointer[6] = source_pointer + width-1;
  neighbor_pointer[7] = source_pointer + width+1;

  SIMD_TYPE current = load_simd_vector(source_pointer);
  SIMD_TYPE n = count_neibours(neighbor_pointer);
  SIMD_TYPE three = _mm256_set1_epi8(3);
  SIMD_TYPE two = _mm256_set1_epi8(2);
  SIMD_TYPE alive_con1 = _mm256_cmpeq_epi8(n,three); // if neibours == 3 then alive
  SIMD_TYPE alive_con2 = _mm256_or_epi64(current, _mm256_cmpeq_epi8(n,two)); // if current == alive and neigbours == 2 then alive
  SIMD_TYPE new = _mm256_or_epi64(alive_con1, alive_con2);
  store_simd_vector(new, destination_pointer);
    }
  }


  
  // TODO traverse through a vector of voxels (length CELLS) and use SIMD intrinsics to implement game of live logic
  // HINT: avoid boundaries
  // HINT: create a SIMD vector for currentfield, newfield and each neighbor (eight in total
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

void game(int width, int height, int num_timesteps) {
  number_type *currentfield = calloc(width * height, sizeof(number_type));
  number_type *newfield = calloc(width * height, sizeof(number_type));
  
  number_type* src; 
  number_type* dst; 

  number_type* neighbor_pointer[8];
  
  // TODO use your favorite filling
  // filling_random (currentfield, width, height);
  filling_runner(currentfield, width, height);

  int time = 0;
  write_field(currentfield, width, height, time);
  // TODO implement periodic boundary condition
  apply_periodic_boundaries(currentfield, width, height);

  for (time = 1; time <= num_timesteps; time++) {
    // TODO assign pointers to correct addresses
    
    // TODO implement evolve function (see above)
    evolve(currentfield, newfield, width, height, src, dst, neighbor_pointer);
    write_field(newfield, width, height, time);
    // TODO implement periodic boundary condition
    apply_periodic_boundaries(newfield, width, height);
    // TODO implement SWAP of the fields
  }

  free(currentfield);
  free(newfield);
}

int main(int c, char **v) {
  int width = 0, height = 0, num_timesteps;
  if (c == 4) {
    width = atoi(v[1]) + 2;     ///< read width + 2 boundary cells (low x, high x)
    height = atoi(v[2]) + 2;    ///< read height + 2 boundary cells (low y, high y)
    num_timesteps = atoi(v[3]); ///< read timesteps
    if (width <= 0) {
      width = 32; ///< default width
    }
    if (height <= 0) {
      height = 32; ///< default height
    }

    double elapsed_time;
    START_TIMEMEASUREMENT(measure_game_time);

    game(width, height, num_timesteps);

    END_TIMEMEASUREMENT(measure_game_time, elapsed_time);
    printf("time elapsed: %lf sec\n", elapsed_time);
  } else {
    myexit("Too less arguments, example: ./gameoflife <x size> <y size> <number of timesteps>");
  }
}