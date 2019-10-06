#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "ExG4DetectorConstruction01.hh"
#include "ExG4PhysicsList00.hh"
#include "ExG4ActionInitialization01.hh"
int main()
{
// construct the default run manager
G4RunManager* runManager = new G4RunManager;
// set mandatory initialization classes
runManager->SetUserInitialization(new ExG4DetectorConstruction01);
runManager->SetUserInitialization(new ExG4PhysicsList00);
runManager->SetUserInitialization(new ExG4ActionInitialization01);
// initialize G4 kernel
runManager->Initialize();
// get the pointer to the UI manager and set verbosities
G4UImanager* UI = G4UImanager::GetUIpointer();
UI->ApplyCommand("/run/verbose 1");
UI->ApplyCommand("/event/verbose 1");
UI->ApplyCommand("/tracking/verbose 1");
// start a run
int numberOfEvent = 3;
runManager->BeamOn(numberOfEvent);
// job termination
delete runManager;
return 0;
}