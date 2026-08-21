function parsave(filename, dataStruct)
% parsave - Helper function to save heavy data safely inside parfor
  save(filename, '-struct', 'dataStruct', '-v7.3');
end