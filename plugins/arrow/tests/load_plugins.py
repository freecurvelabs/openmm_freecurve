import openmm

platform = openmm.Platform.getPlatformByName("Reference")
try:
  platform.loadPluginsFromDirectory("/home/igor/PROG_SRC/OPENMM_PYTHON_LINUX/lib")
except Exception as e:
  print(f"Failed to load plugins: {e}")
  
print("Loaded plugins:")
for plugin in openmm.pluginLoadedLibNames:
  print(plugin)



