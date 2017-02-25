# create code for each example ready to compile

# reduced gravity configurations
# f_plane.f90 - boring test with flat interface on f-plane
sed '37s/.*/    integer, parameter :: nx = 100/' ../MIM.f90 > f_plane.f90
sed '38s/.*/    integer, parameter :: ny = 100/' f_plane.f90 > f_plane2.f90
sed '39s/.*/    integer, parameter :: layers = 1/' f_plane2.f90 > ./reduced_gravity/f_plane/f_plane.f90

rm f_plane.f90
rm f_plane2.f90


# beta_plane_bump - test with a bump on a beta-plane
sed '37s/.*/    integer, parameter :: nx = 100/' ../MIM.f90 > beta_plane_bump.f90
sed '38s/.*/    integer, parameter :: ny = 100/' beta_plane_bump.f90 > beta_plane_bump2.f90
sed '39s/.*/    integer, parameter :: layers = 1/' beta_plane_bump2.f90 > ./reduced_gravity/beta_plane_bump/beta_plane_bump.f90

rm beta_plane_bump.f90
rm beta_plane_bump2.f90



# beta_plane_gyre - twin-gyre simulation on a beta-plane
sed '37s/.*/    integer, parameter :: nx = 100/' ../MIM.f90 > beta_plane_gyre.f90
sed '38s/.*/    integer, parameter :: ny = 200/' beta_plane_gyre.f90 > beta_plane_gyre2.f90
sed '39s/.*/    integer, parameter :: layers = 1/' beta_plane_gyre2.f90 > ./reduced_gravity/beta_plane_gyre/beta_plane_gyre.f90

rm beta_plane_gyre.f90
rm beta_plane_gyre2.f90


# n layer configurations
# f_plane.f90 - boring test with flat interface on f-plane
sed '37s/.*/    integer, parameter :: nx = 100/' ../MIM.f90 > f_plane.f90
sed '38s/.*/    integer, parameter :: ny = 100/' f_plane.f90 > f_plane2.f90
sed '39s/.*/    integer, parameter :: layers = 2/' f_plane2.f90 > ./n_layer/f_plane/f_plane.f90

rm f_plane.f90
rm f_plane2.f90


# beta_plane_bump - test with a bump on a beta-plane
sed '37s/.*/    integer, parameter :: nx = 100/' ../MIM.f90 > beta_plane_bump.f90
sed '38s/.*/    integer, parameter :: ny = 100/' beta_plane_bump.f90 > beta_plane_bump2.f90
sed '39s/.*/    integer, parameter :: layers = 2/' beta_plane_bump2.f90 > ./n_layer/beta_plane_bump/beta_plane_bump.f90

rm beta_plane_bump.f90
rm beta_plane_bump2.f90



# beta_plane_gyre - twin-gyre simulation on a beta-plane
sed '37s/.*/    integer, parameter :: nx = 100/' ../MIM.f90 > beta_plane_gyre.f90
sed '38s/.*/    integer, parameter :: ny = 200/' beta_plane_gyre.f90 > beta_plane_gyre2.f90
sed '39s/.*/    integer, parameter :: layers = 2/' beta_plane_gyre2.f90 > ./n_layer/beta_plane_gyre/beta_plane_gyre.f90

rm beta_plane_gyre.f90
rm beta_plane_gyre2.f90


# run the example simulations
echo 'running reduced gravity models'
echo '   f plane example'
# run f_plane example
cd reduced_gravity/f_plane
rm -r output
rm -r renders
rm run_finished.txt

gfortran f_plane.f90 -o f_plane.out -Ofast
if [ $? -ne 0 ]; then
  echo "Error while compiling model."
  exit
fi
mkdir output
mkdir renders
./f_plane.out
if [ $? -ne 0 ]; then
  echo "Error while running model."
  exit
fi
cd ../


echo '   beta plane bump'
# run beta_plane_bump example
cd beta_plane_bump
rm -r output
rm -r renders
rm run_finished.txt

gfortran beta_plane_bump.f90 -o beta_plane_bump.out -Ofast
if [ $? -ne 0 ]; then
  echo "Error while compiling model."
  exit
fi
mkdir output
mkdir renders
./beta_plane_bump.out
if [ $? -ne 0 ]; then
  echo "Error while running model."
  exit
fi
cd ../

echo '   beta plane gyres'
echo '      no slip'
# run beta_plane_gyre example
# no slip
cd beta_plane_gyre
gfortran beta_plane_gyre.f90 -o beta_plane_gyre.out -Ofast
if [ $? -ne 0 ]; then
  echo "Error while compiling model."
  exit
fi

cd no_slip
rm -r output
rm -r renders
rm run_finished.txt

mkdir output
mkdir renders
../beta_plane_gyre.out
if [ $? -ne 0 ]; then
  echo "Error while running model."
  exit
fi

echo '      free slip'
cd ../free_slip
rm -r output
rm -r renders
rm run_finished.txt

mkdir output
mkdir renders
../beta_plane_gyre.out
if [ $? -ne 0 ]; then
  echo "Error while running model."
  exit
fi

cd ../../../


# Run n layer examples
echo 'Run n layer examples'

echo '   f plane example'
# run f_plane example
cd n_layer/f_plane
rm -r output
rm -r renders
rm run_finished.txt

gfortran f_plane.f90 -o f_plane.out -Ofast
if [ $? -ne 0 ]; then
  echo "Error while compiling model."
  exit
fi
mkdir output
mkdir renders
./f_plane.out
if [ $? -ne 0 ]; then
  echo "Error while running model."
  exit
fi
cd ../


echo '   beta plane bump'
# run beta_plane_bump example
cd beta_plane_bump
rm -r output
rm -r renders
rm run_finished.txt

gfortran beta_plane_bump.f90 -o beta_plane_bump.out -Ofast
if [ $? -ne 0 ]; then
  echo "Error while compiling model."
  exit
fi
mkdir output
mkdir renders
./beta_plane_bump.out
if [ $? -ne 0 ]; then
  echo "Error while running model."
  exit
fi
cd ../

echo '   beta plane gyres'
echo '      no slip'
# run beta_plane_gyre example
# no slip
cd beta_plane_gyre
gfortran beta_plane_gyre.f90 -o beta_plane_gyre.out -O1
if [ $? -ne 0 ]; then
  echo "Error while compiling model."
  exit
fi

cd no_slip
rm -r output
rm -r renders
rm run_finished.txt

mkdir output
mkdir renders
../beta_plane_gyre.out
if [ $? -ne 0 ]; then
  echo "Error while running model."
  exit
fi

echo '      free slip'
cd ../free_slip
rm -r output
rm -r renders
rm run_finished.txt

mkdir output
mkdir renders
../beta_plane_gyre.out
if [ $? -ne 0 ]; then
  echo "Error while running model."
  exit
fi
cd ../../
echo 'Finished running examples'