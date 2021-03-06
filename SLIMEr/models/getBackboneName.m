%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% backName = getBackboneName(specName)
%
% Benjam�n J. S�nchez. Last update: 2018-03-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function backName = getBackboneName(specName)

backName = '';
compName = specName(strfind(specName,'[')+1:strfind(specName,']')-1);
specName = specName(1:strfind(specName,'[')-2);

%Group 1: keep generic name
# do lado esquerdo temos a chave e do lado direito temos o valor para a chave do dicionário
# group1 = é um dicionário
group1 = {'1-phosphatidyl-1D-myo-inositol'                      'cytoplasm'
          'sn-2-acyl-1-lysophosphatidylinositol'                'endoplasmic reticulum membrane'
          'phosphatidyl-L-serine'                               'endoplasmic reticulum membrane'
          'phosphatidylcholine'                                 'endoplasmic reticulum membrane'
          'phosphatidylethanolamine'                            'endoplasmic reticulum membrane'
          'phosphatidate'                                       'endoplasmic reticulum membrane'
          'diglyceride'                                         'endoplasmic reticulum membrane'
          'triglyceride'                                        'endoplasmic reticulum membrane'
          'phosphatidylglycerol'                                'mitochondrial membrane'
          'cardiolipin'                                         'mitochondrial membrane'
          'ceramide'                                            'Golgi'
          'inositol-P-ceramide'                                 'Golgi'
          'inositol phosphomannosylinositol phosphoceramide'    'Golgi'
          'mannosylinositol phosphorylceramide'                 'Golgi'};

# as funções seguintes iteram sobre o grupo (sobre as chaves do dicionário)
for i = 1:length(group1)
    if startsWith(specName,group1{i,1}) && ~strcmp(specName,group1{i,1}) && ...
                  strcmp(compName,group1{i,2}) && ~contains(specName,'phosphate')
        backName = [group1{i,1} ' [' compName ']'];
    end
end

%Group 2: replace specific name by generic
group2 = {'palmitate'                       'fatty acid'                'cytoplasm'
          'palmitoleate'                    'fatty acid'                'cytoplasm'
          'stearate'                        'fatty acid'                'cytoplasm'
          'oleate'                          'fatty acid'                'cytoplasm'
          'ergosterol'                      'ergosterol'                'cytoplasm'
          'ergosteryl palmitoleate'         'ergosterol ester'          'endoplasmic reticulum membrane'
          'ergosteryl oleate'               'ergosterol ester'          'endoplasmic reticulum membrane'
          'sphinganine'                     'long-chain base'           'endoplasmic reticulum'
          'phytosphingosine'                'long-chain base'           'endoplasmic reticulum'
          'sphinganine 1-phosphate'         'long-chain base phosphate' 'endoplasmic reticulum'
          'phytosphingosine 1-phosphate'	'long-chain base phosphate' 'endoplasmic reticulum'};
for i = 1:length(group2)
    if strcmp(specName,group2{i,1}) && strcmp(compName,group2{i,3}) 
        backName = [group2{i,2} ' [' compName ']'];
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%