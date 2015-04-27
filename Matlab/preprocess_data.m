if ~exist('HumanPPI700', 'var')
    tic
    % HumanPPI = import_file('../output/HumanPPI.tab');
    % HumanPPI700 = import_file('../output/HumanPPI700.tab');
    % HumanPPI900 = import_file('../output/HumanPPI900.tab');
    % load('../output/HumanPPI.mat');
    load('../output/HumanPPI700.mat');
    load('../output/HumanPPI900.mat');
    fprintf('Data loaded in: %.3f ms\n', toc);
end

if ~exist('protein_set700', 'var') && ~exist('protein_set900', 'var') 
    tic
    protein_set700 = union(HumanPPI700.protein1, HumanPPI700.protein2);
    protein_set900 = union(HumanPPI900.protein1, HumanPPI900.protein2);
    fprintf('Protein set created in: %.3f ms\n', toc);
end

if ~exist('go_knowledge', 'var')
    tic
    go_knowledge = import_go_full('../data/9606_go_knowledge_full.tsv');
    fprintf('Go knowledge imported in: %.3f ms\n', toc);
end

tic
fprintf('Filtering data...\n');
go_string_ids = strcat(num2str(go_knowledge.var1), '.', go_knowledge.string_id);
mask_GO700 = ismember(go_string_ids, protein_set700);
mask_GO900 = ismember(go_string_ids, protein_set900);

HumanPPI700_GO = go_knowledge(mask_GO700, :);
HumanPPI700_GO_conf = HumanPPI700_GO(HumanPPI700_GO.confidentiality >= 3, :);
HumanPPI900_GO = go_knowledge(mask_GO900, :);
HumanPPI900_GO_conf = HumanPPI900_GO(HumanPPI900_GO.confidentiality >= 3, :);
fprintf('Data filtered in: %.3f ms\n', toc);
