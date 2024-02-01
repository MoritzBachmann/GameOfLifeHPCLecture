      /*
      if (rank == 0)
      {
        global_field = calloc(gsizes[0] * gsizes[1], sizeof(number_type));
        for (int i = 1; i < process_numX * process_numY; i++)
        {
          MPI_Cart_coords(cart_comm, i, 2, coords);
          // cart_comm
          MPI_Status rec;
          printf("iD: %d - reciving  from rank %d \n", rank, i);
          printf("iD: %d - cords: %d %d \n", rank, coords[0], coords[1]);
          MPI_Recv(global_field, 1, filetype, i, 99, cart_comm, &rec);
          for (int i = 1; i < height; i++)
          {
            printf("\n");
            for (int j = 1; j < width; j += 1)
            {
              int pos = i * width + j;
              printf("%d ", (int)*(global_field + pos));
            }
          }
        }
      }
      else
      {
        // local_start_indices[0]][local_start_indices[1]
        printf("iD: %d - sendig  from rank \n", rank);
        MPI_Send(local_field, 1, filetype, 0, 99, cart_comm);
      }

    */