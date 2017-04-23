NAME =
CC = gcc
CCFLAGS = -O3 -ansi -pedantic -Wall
OBJ =

all: $(NAME)

$(NAME): $(NAME).c $(OBJ)
	$(CC) -o $(NAME) $(NAME).c $(OBJ) $(CCFLAGS)

%.o: %.c %.h
	$(CC) -c -o $@ $< $(CCFLAGS)

.PHONY: clean

clean:
	rm -f $(OBJ)
