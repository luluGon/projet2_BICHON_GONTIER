# Compiler
CC = gcc
CFLAGS = -O2 
LFLAGS = -lm
LAPFLAGS = -llapacke

# Directories
SRC_DIR = prog
OBJ_DIR = obj
INCLUDE_DIR = headers

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))

# Exécutable
EXECUTABLE = exe

# Compilation
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -c $< -o $@ 

# Création de l'exécutable
$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) -o $@  $^ $(LFLAGS) $(LAPFLAGS)

# Cibles pour nettoyer les fichiers objets et l'exécutable
clean:
	rm -rf $(OBJ_DIR) $(EXECUTABLE)

# Définition de la cible par défaut
all: $(EXECUTABLE)

# Lancement de l'exécutable
run: $(EXECUTABLE)
	./$(EXECUTABLE)

.PHONY: all clean run
