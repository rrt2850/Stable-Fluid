param(
    [int]$GridWidth = 128,
    [int]$GridHeight = 128,
    [int]$EmitterCount = 8,
    [int]$SimulationSteps = 180,
    [bool]$ExportVideo = $false,
    [bool]$LogEveryStep = $false
)

$ErrorActionPreference = 'Stop'

$BuildDir = 'out/classes'
New-Item -ItemType Directory -Force -Path $BuildDir | Out-Null

javac -d $BuildDir src/*.java
java -cp $BuildDir Main $GridWidth $GridHeight $EmitterCount $SimulationSteps $ExportVideo $LogEveryStep
