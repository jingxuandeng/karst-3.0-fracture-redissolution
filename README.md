
 # Network model of dissolving porous material

 ## DESCRIPTION

 This project is dedicated to simulate porous materials subjected to dissolution processes.
 Porous material is modeled as a network of tubes (pore throats) and is represented by the class Network.
 Pressure drop applied to the network edges results in fluid (water plus reactants) flow through the system.
 Pore sizes can change due to chemical reactions with reagents dissolved in the fluid.
 The main reaction here is a dissolution one. Apart form that he second type of reaction, precipitation one,
  can be additionally tracked.
 The pore size can change in two ways - either grow due to dissolution or shrinker due to precipitation.

 ### REPRESENTATION OF POROUS MATERIAL

 Porous material is modeled here as a network of interconnected tubes called pores.
 Pores are represented by the class Pore.
 Each pore is spanned between two nodes of the network. Nodes are represented by the class Node.
 Additionally, apart of pores  and nodes, the simulation can also track the pieces of material
  that are subjected to the dissolution, so called grains, represented by the class Grain.


 ### EVOLUTION OF THE SYSTEM

 The main purpose of this project is to simulate the dynamics of the system.
 The function Network::evolution will perform T time steps consisted of

 - calculating pressure field in the nodes

 - calculating flow field in the pores

 - calculating concentration filed in the nodes and pores

 - updating pores and grains new shape

 - (additionally) changing topology of the network if some grains has vanished due to the dissolution

 The simulation ends either after T time steps (mostly debugging mode)
 or if other condition connected to the network properties is fulfilled.
 We typically use a condition for the breakthrough - the simulation ends when the dissolution pattern
 (pattern consisted of broad, dissolved pores) has reached the outlet of the system.



## Other remarks:

1. Branches:

    - master: two reactions without transversal diffusion;

    - pure diffusion: transversal diffusion is added for two reactions (works but have to be tested more deeply)

    - pure diffusion transversal diffusion works for first reaction (dissolution). 

    - other branches are connected with work in progress