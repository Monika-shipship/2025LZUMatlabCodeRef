% =========================================================================
% convert_figs_to_pdf.m
%
% 功能：
%   递归搜索当前文件夹（FigOutput）及其所有子文件夹中的所有 .fig 文件，
%   并将它们一一转换为 .pdf 文件。
%   转换后的 .pdf 文件将保存在与原始 .fig 文件相同的位置。
%
% 如何使用：
%   1. 将此 .m 文件保存在 'FigOutput' 文件夹中。
%   2. 在MATLAB中打开此文件。
%   3. 点击编辑器顶部的“运行”（Run）按钮。
% =========================================================================

fprintf('开始批量转换 .fig 到 .pdf ...\n');

% 1. 递归查找当前目录及子目录下的所有 .fig 文件
%    因为此脚本放在 FigOutput 中，所以会从 FigOutput 开始搜索
figFiles = dir('**\*.fig');

numFiles = length(figFiles);

if numFiles == 0
    fprintf('在当前目录或子目录中未找到任何 .fig 文件。\n');
    return;
end

fprintf('找到了 %d 个 .fig 文件。开始处理...\n\n', numFiles);

successCount = 0;
failCount = 0;

% 2. 循环遍历所有找到的文件
parfor i = 1:numFiles
    % 获取 .fig 文件的完整路径
    figFolder = figFiles(i).folder;
    figName = figFiles(i).name;
    fullFigPath = fullfile(figFolder, figName);
    
    % 3. 创建目标 .pdf 文件的路径
    %    [~, baseName, ~] = fileparts(figName) 会得到不带扩展名的文件名
    [~, baseName, ~] = fileparts(figName); 
    pdfName = [baseName, '.pdf'];
    fullPdfPath = fullfile(figFolder, pdfName);
    
    fprintf('正在处理 (%d/%d): %s\n', i, numFiles, fullFigPath);

    % 4. 使用 try-catch 块来防止单个文件出错导致整个脚本停止
    try
        % 5. 打开 .fig 文件（设置为不可见）
        hFig = openfig(fullFigPath, 'invisible');
        
        % 6. 将图形导出为 PDF (这是R2020a及以后版本推荐的方法)
        %    它通常比 saveas 提供更好的输出质量
        exportgraphics(hFig, fullPdfPath, 'ContentType', 'vector');
        
        % 7. 关闭已打开的图形，释放内存
        close(hFig);
        
        fprintf('  -> 成功保存到: %s\n', fullPdfPath);
        successCount = successCount + 1;
        
    catch ex
        fprintf('  *** 转换失败: %s\n', figName);
        fprintf('  *** 错误信息: %s\n', ex.message);
        failCount = failCount + 1;
        
        % 如果打开出错，确保句柄被关闭
        if exist('hFig', 'var') && ishandle(hFig)
            close(hFig);
        end
    end
    
    fprintf('--------------------------------------------------\n');
end

% 8. 打印最终总结
fprintf('\n批量转换完成。\n');
fprintf('成功: %d 个\n', successCount);
fprintf('失败: %d 个\n', failCount);
