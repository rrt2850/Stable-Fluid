param(
    [int]$EmitterCount = 8
)

$ErrorActionPreference = 'Stop'

javac src/*.java
java -cp src Main $EmitterCount
