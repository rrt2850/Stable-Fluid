param(
    [int]$GridWidth = 128,
    [int]$GridHeight = 128,
    [int]$EmitterCount = 8,
    [int]$SimulationSteps = 180
)

$ErrorActionPreference = 'Stop'

javac src/*.java
java -cp src Main $GridWidth $GridHeight $EmitterCount $SimulationSteps
