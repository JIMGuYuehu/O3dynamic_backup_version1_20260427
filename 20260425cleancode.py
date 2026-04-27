import os
import nbformat
import shutil
from pathlib import Path

# --- 配置 ---
source_root = "." 
target_root = "20260425code_clean"
# 直接拷贝的后缀
copy_extensions = {".py", ".ncl", ".sh"} 
# -----------

def clean_ipynb(file_path, save_path):
    """清除 ipynb 文件的输出并保存，处理空文件或坏文件"""
    try:
        # 检查文件大小，如果是 0 字节直接跳过
        if os.path.getsize(file_path) == 0:
            print(f"Skipping empty file: {file_path}")
            return False
            
        with open(file_path, 'r', encoding='utf-8') as f:
            nb = nbformat.read(f, as_version=4)
        
        for cell in nb.cells:
            if cell.cell_type == 'code':
                cell.outputs = []
                cell.execution_count = None
        
        with open(save_path, 'w', encoding='utf-8') as f:
            nbformat.write(nb, f)
        return True
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return False

# 1. 创建目标根目录
target_path_obj = Path(target_root)
if target_path_obj.exists():
    shutil.rmtree(target_path_obj) # 如果目录已存在，先删除以确保干净
target_path_obj.mkdir()

print(f"开始全面整理代码到: {target_root} ...\n")

# 2. 递归遍历
for root, dirs, files in os.walk(source_root):
    # 排除目标目录本身
    if target_root in root:
        continue
    
    for filename in files:
        name, ext = os.path.splitext(filename)
        ext = ext.lower()
        
        if ext == ".ipynb" or ext in copy_extensions:
            src_path = os.path.join(root, filename)
            
            # 建立对应层级的目录
            rel_path = os.path.relpath(root, source_root)
            dest_dir = target_path_obj / rel_path
            dest_dir.mkdir(parents=True, exist_ok=True)
            dest_file_path = dest_dir / filename

            if ext == ".ipynb":
                clean_ipynb(src_path, dest_file_path)
            else:
                shutil.copy2(src_path, dest_file_path)

print(f"\n--- 整理完成！ ---")
print(f"请检查目录: {os.path.abspath(target_root)}")