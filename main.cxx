#include <vtkCellArray.h>
#include <vtkProperty.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkSmartPointer.h>
#include <vtkMath.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCleanPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkOBJReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkParametricFunctionSource.h>
#include <vtkParametricSuperEllipsoid.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyDataMapper.h>
#include<vtkParametricEllipsoid.h>
#include<vtkCylinderSource.h>
#include<Windows.h>
#include <vtkUnstructuredGrid.h>
#include <vtkOutlineFilter.h>
#include <vtkClipPolyData.h>
#include <vtkPlane.h>
void liberer(double **tab,int nLigne){
 
  for(int i=0 ; i < nLigne ; i++)
    delete[] tab[i];
  delete[] tab;
}

int main ( int argc, char *argv[] )
{
                  /////////////////////////////////////////////////////////////////////////////////////////
                 //*******affichage du femur, son centroide, matrice contenant ses sommets**************///
                //////////////////////////////////////////////////////////////////////////////////////////

// lire le femur
vtkOBJReader *femur = vtkOBJReader::New();
femur->SetFileName("femur.obj");


 // Define a clipping plane
  vtkSmartPointer<vtkPlane> clipPlane = 
    vtkSmartPointer<vtkPlane>::New();
  clipPlane->SetNormal(1, 1,1);
  clipPlane->SetOrigin(-0.0040,0.00118,0.00365);
  
// Clip the source with the plane
  vtkSmartPointer<vtkClipPolyData> fem = 
    vtkSmartPointer<vtkClipPolyData>::New();
   fem->SetInputConnection(femur->GetOutputPort());
   fem->SetClipFunction(clipPlane);



// mettre le femur dans un polydata: fpol
vtkPolyData *fpol =femur->GetOutput();

femur->Update();
/////////////////longueur de la diagonal du bounding box

  double bounds[6];
  fpol->GetBounds(bounds);
  double xmin=0,xmax=0,ymin=0,ymax=0,zmin=0,zmax=0;
  xmin=bounds[0];
  xmax=bounds[1];
  ymin=bounds[2];
  ymax=bounds[3];
  zmin=bounds[4];
  zmax=bounds[5];
double diagonal=0;
diagonal=sqrt(pow(2,(xmax-xmin))+pow(2,(ymax-ymin))+pow(2,(zmax-zmin)));
std::cout  << "diagonal:  " <<  diagonal << std::endl;

///////////////////////////////////////////////////////

///allocation dynamique pr le tab femur
double **F;
int l=fpol->GetNumberOfPoints();
F=new double*[l];
for (int i=0;i<l;i++)
{F[i]=new double[3];}

//affichage des sommets du femur et calcul de son centroide C(XC,YC,ZC)
double XC=0, YC=0, ZC=0;
  for(vtkIdType i = 0; i < fpol->GetNumberOfPoints(); i++)
    {
    double p[3];
    fpol->GetPoint(i,p);
	 F[i][0]=p[0];
	 F[i][1]=p[1];
	 F[i][2]=p[2];
	//std::cout << "F[" << i <<"]" <<" : (" << F[i][0] << " " << F[i][1] << " " << F[i][2] << ")" << std::endl;

     XC+= p[0];
     YC+= p[1];
     ZC+= p[2];
  }

  //calcul du centroide    
 
 XC/= fpol->GetNumberOfPoints();
 YC/= fpol->GetNumberOfPoints();
 ZC/= fpol->GetNumberOfPoints();
 std::cout << "XC " << XC << std::endl;
 std::cout << "YC " << YC << std::endl;
 std::cout << "ZC " << ZC << std::endl;


// mapper et actor pour l'affichage du femur 
vtkPolyDataMapper *map = vtkPolyDataMapper::New(); 
map->SetInput(fpol);
vtkActor *factor = vtkActor::New();
factor->SetMapper(map);
factor->GetProperty()->SetColor(1,1,0);  
factor->GetProperty()->SetOpacity(0.6);
factor->GetProperty()->SetColor(1,1,1);
factor->RotateX(-90);

//affichage du centroide
vtkPoints *points = vtkPoints::New();
  const float p[3] = {XC,YC, ZC};
  vtkCellArray *vertex = vtkCellArray::New();
  vtkIdType pid[1];
  pid[0] = points->InsertNextPoint(p);
  vertex->InsertNextCell(1,pid);
  vtkPolyData *point = vtkPolyData::New();
  point->SetPoints(points);
   point->SetVerts(vertex); 
 
  
  // mapper et actor pour le centroide
  vtkPolyDataMapper *cmapp = vtkPolyDataMapper::New();
cmapp->SetInput(point);
 
 vtkActor *cactor =vtkActor::New();
  cactor->SetMapper(cmapp);
  cactor->GetProperty()->SetPointSize(7);
  cactor->GetProperty()->SetColor(1,0,0);

  // Create the outline
  vtkSmartPointer<vtkOutlineFilter> outline = 
    vtkSmartPointer<vtkOutlineFilter>::New();
  outline->SetInput(fpol);
  vtkSmartPointer<vtkPolyDataMapper> outlineMapper = 
    vtkSmartPointer<vtkPolyDataMapper>::New();
  outlineMapper->SetInput(outline->GetOutput());
  vtkSmartPointer<vtkActor>outlineActor = 
    vtkSmartPointer<vtkActor>::New();
  outlineActor->SetMapper(outlineMapper);
  outlineActor->GetProperty()->SetColor(0,0,0);
  outlineActor->RotateX(-90);
                  ////////////////////////////////////////////////////////////////////////////////////////////
                 //******************* calcul de la surface externe du maillage du femur*******************//
                ////////////////////////////////////////////////////////////////////////////////////////////

//definir un filtre qui prend comme Input un polydata et génére aussi un polydata comme Output  
vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
cleaner->SetInputConnection (fpol->GetProducerPort());
 
//générer un objet vtkDelaunay3D pour le femur 

 vtkDelaunay3D *delaunay3D = vtkDelaunay3D::New();
  delaunay3D->SetInputConnection (cleaner->GetOutputPort());
  delaunay3D->SetAlpha(0.0125);
  //delaunay3D->SetOffset(0.5);
 delaunay3D->SetTolerance(0);
 
  // mapper and actor for delaunay fem

  vtkDataSetMapper *delaunayMapper = vtkDataSetMapper::New();
  delaunayMapper->SetInputConnection(delaunay3D->GetOutputPort());
  vtkActor *delaunayActor = vtkActor::New();
  delaunayActor->SetMapper(delaunayMapper);
  delaunayActor->GetProperty()->SetColor(1,1,1);
  delaunayActor->GetProperty()->SetOpacity(0.2);
 

                  /////////////////////////////////////////////////////////////////////////////////////////
                 //*****************************modèle retréci du femur***********************************///
                //////////////////////////////////////////////////////////////////////////////////////////

vtkTransform *scale = vtkTransform::New();
scale->Scale(0.9,0.9,0.9);
scale->Translate(-0.005,0,0);
 vtkTransformPolyDataFilter *scaleFilter = vtkTransformPolyDataFilter::New();
 scaleFilter->SetInputConnection(fem->GetOutputPort());
 scaleFilter->SetTransform(scale);
 scaleFilter->Update();

 vtkPolyData *scalepol = vtkPolyData::New();
scalepol=scaleFilter->GetOutput();
scaleFilter->Update();
//delaunay pr femur retreci
  vtkDelaunay3D *scaledelaunay = vtkDelaunay3D::New();   
  scaledelaunay->SetInputConnection(scaleFilter->GetOutputPort());
  scaledelaunay->SetAlpha(0.009);
  scaledelaunay->SetTolerance(0.1);


   // mapper and actor for scale

  vtkDataSetMapper *scaleMapper = vtkDataSetMapper::New();
  scaleMapper->SetInputConnection(scaledelaunay->GetOutputPort());
  vtkActor *scaleActor = vtkActor::New();
  scaleActor->SetMapper(scaleMapper);
  scaleActor->GetProperty()->SetColor(1,1,0);
  scaleActor->GetProperty()->SetOpacity(0.1);
  scaleActor->RotateX(-90);
  scaleActor->SetPosition(0.005,0,0);
  scaleActor->SetOrigin(factor->GetOrigin());
  
/////////////////////////////
   // Creer a renderer, render window, and interactor
  vtkRenderer *Renderer = vtkRenderer::New();
 vtkRenderer *renderer = vtkRenderer::New();
  vtkRenderWindow *renderWindow = vtkRenderWindow::New();
  renderWindow->SetSize(600,600);
 
  renderWindow->AddRenderer(Renderer);
  
  vtkRenderWindowInteractor *renderWindowInteractor = vtkRenderWindowInteractor::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
 // Renderer->SetBackground(.3, .6, .3);
 
 // Render et interact

                  /////////////////////////////////////////////////////////////////////////////////////////
                 //*****************************générer un elipsoide***********************************///
                //////////////////////////////////////////////////////////////////////////////////////////
vtkParametricEllipsoid *ellips = vtkParametricEllipsoid::New();
  vtkParametricFunctionSource *pfs = vtkParametricFunctionSource::New();
  pfs->SetParametricFunction(ellips);
//  //modifier les paramètre de l'ellipsoide
//pfs->SetVResolution(100);
//pfs->SetUResolution(100);
//pfs->SetWResolution(100);
// mettre l'ellipsoide dans un polydata epol
  vtkPolyData *epol = vtkPolyData::New();
  epol=pfs->GetOutput();
  epol->Update();
//std::cout << epol->GetNumberOfPoints() << std::endl; 
vtkCylinderSource *implant = vtkCylinderSource::New();

  implant->SetRadius(0.001);
  implant->SetHeight(0.03);
  implant->SetResolution(100);
 
  // Create a mapper and actor
 vtkPolyDataMapper *imapper =vtkPolyDataMapper::New();
  imapper->SetInputConnection(implant->GetOutputPort());
  vtkActor *iactor =vtkActor::New();
  iactor->SetMapper(imapper);

implant->SetCenter(-0.004,0.00118,0.00365);
implant->Update();
iactor->GetProperty()->SetColor(0,0,1); 
////////////////////////////////////////////////
  
  // mapper and actor for ellips

  vtkDataSetMapper *eMapper = vtkDataSetMapper::New();
  eMapper->SetInputConnection(pfs->GetOutputPort());
  vtkActor *eActor = vtkActor::New();
  eActor->SetMapper(eMapper);
  //eActor->SetPosition(-0.009,-0.0008,-0.01);
  eActor->SetPosition(0.004,0.015,-0.001);
  //eActor->GetProperty()->SetOpacity(0.4);
   eActor->RotateZ(45);
eActor->SetOrigin(XC,YC,ZC);
 Renderer->SetBackground(.2, .3, .4);

for(double i=0.01; i<0.05; i+=0.01)
{	
	renderWindow->FullScreenOn();
	
ellips->SetXRadius(i); 
ellips->SetYRadius(0.01);
ellips->SetZRadius(0.01);
eActor->SetMapper(eMapper);
//factor->SetOrientation(1,1,1);
 eActor->GetProperty()->SetColor(0,1,1);
 //factor->RotateX(90);
// pfsActor->AddPosition(0.009,0.0008,0.01);
  std::cout << i <<std::endl;
  Renderer->AddActor(eActor);
 Renderer->AddActor(factor);
 Renderer->AddActor(outlineActor);
 Renderer->AddActor(scaleActor);

 renderWindow->Render();
Sleep(1000);
pfs->Update();

}




                  /////////////////////////////////////////////////////////////////////////////////////////
                 //*****************************enveloppe interne***********************************///
                //////////////////////////////////////////////////////////////////////////////////////////


epol->Update();
vtkPoints *vertices=vtkPoints::New();
for(vtkIdType i = 0; i <epol->GetNumberOfPoints(); i++)
    {   
         double S[3];
       epol->GetPoint(i,S);
  double pcoords[3], weights[3];
 
  vtkIdType cellId;
  
  int subId;
  cellId =scaledelaunay->GetOutput()->FindCell(S,NULL, 0, .1,subId, pcoords, weights);
 
pfs->Update();

  if (cellId>=0) ////les points sont à l'interieur
    {	
 
	 vertices->InsertNextPoint(S[0],S[1],S[2]);

    }
  
 }

vtkPolyData *pointpol=vtkPolyData::New();
pointpol->SetPoints(vertices);
pointpol->Update();

  //delaunay pr les points inside
vtkDelaunay3D *indelaunay = vtkDelaunay3D::New();
indelaunay->SetInputConnection (pointpol->GetProducerPort());
indelaunay->SetAlpha(0.5);


///mapper and actor for points inside
vtkDataSetMapper *inMapper = vtkDataSetMapper::New();
 inMapper->SetInputConnection(indelaunay->GetOutputPort());
 vtkActor *inActor = vtkActor::New();
 inActor->SetMapper(inMapper);
 inActor->RotateZ(35);
inActor->SetPosition(0.004,0.015,-0.01);
inActor->GetProperty()->SetColor(1,0,1); 
//inActor->SetOrigin(XC,YC,ZC);

//inActor->GetProperty()->SetRepresentationToPoints();

eActor->VisibilityOff();
 outlineActor->VisibilityOff();
 renderer->AddActor(iactor);
  
  Renderer->AddActor(inActor);

                  /////////////////////////////////////////////////////////////////////////////////////////
                 //*****************************generer un cylindre immplant***********************************///
                //////////////////////////////////////////////////////////////////////////////////////////


Sleep(1000);

	iactor->RotateZ(80);

for(double i=0.001; i<0.09; i+=0.01)
{	
renderWindow->FullScreenOn();
iactor->SetPosition(-i,i,0.03);
Renderer->AddActor(iactor);
renderWindow->Render();
Sleep(1000);
implant->Update();

}

delaunayActor->VisibilityOff();


                       ///////***visualisation***//////////
  
 //Renderer->AddActor(delaunayActor);//surface externe du femur
 //Renderer->AddActor(factor); //femur initial
 //Renderer->AddActor(cactor);//centroide du femur
   //Renderer->AddActor(tActor); //ellipsoide transformé
   //Renderer->AddActor(oActor); // ellipsoide original
   //Renderer->AddActor(scaleActor);//femur retreci
 //enveloppe interne
  renderWindow->Render();

 renderWindowInteractor->Start();

 
  liberer(F,l);
//  liberer(E,g);
  return EXIT_SUCCESS;
} 