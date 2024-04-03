#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 20:51:28 2024

@author: Rui
"""

import os
import sys ; sys.setrecursionlimit(sys.getrecursionlimit() * 5)
import pandas as pd
from tqdm import tqdm

def main():
    sys.path.append('/Users/Rui/anaconda3/envs/Bio_analysis/PY_documents')
    import message_check
    import Bio_analysis_redo
    print('程序开始执行')
    #current_directory = os.getcwd()
    #os.chdir('/Users/Rui/Documents/Bio_analysis')
    current_directory = os.getcwd()
    print("当前运行文件所在目录:", current_directory)
    progress = input("新建项目名称, 或者输入您已有的项目名称:")
    Taxonomy = input("请输入您研究对象的种属的Taxonomy编号(若未知请输入'n'):")
    if Taxonomy == 'n':
        species_name = input("请输入您要查询的物种名(英文)")  # 你要查询的物种名
        taxonomy_id = Bio_analysis_redo.get_taxonomy_id(species_name)
        print(f"Species: {species_name}, Taxonomy ID: {taxonomy_id}")
        Taxonomy = taxonomy_id
    Taxonomy_list = [int(Taxonomy)]
    sys.path.append(current_directory)
    try:
        os.makedirs(current_directory+'/'+progress)#新建项目文件夹
    except:
        pass

    this_program_path = current_directory+'/'+progress

    operation = input("您想执行什么操作？\n1.数据库粗筛\n2.DEG分析\n3.GO富集分析\n4.KEGG富集分析\n5.基因亚细胞结构定位\n输入项目的序号(阿拉伯数字,支持多个输入)即可:")
    
    matrix_save   = current_directory +'/'+ progress +'/'+ 'matrix'   +'/'
    analysis_save = current_directory +'/'+ progress +'/'+ 'analysis' +'/'
    save_path_GPL = current_directory +'/'+ progress +'/'+ 'platform' +'/'
    try:
        os.makedirs(analysis_save)#新建项目子文件夹
        os.makedirs(save_path_GPL)
        os.makedirs(matrix_save)
    except:
        pass
    
    if '1' not in operation:
        try:
            xiangxi_list = message_check.output_file_find_and_open(this_program_path, current_directory, progress)
            if xiangxi_list == None:
                xiangxi_list = input('请输入您要分析的数据库(用空格分隔开):').split(' ')
        except:
            xiangxi_list = input('请输入您要分析的数据库(用空格分隔开):').split(' ')
    
    if '1' in operation or '一' in operation:
        geo_ids = message_check.gds_result_file_check(input("请输入您从GEO网站得出的文件,注意包含路径和后缀名:"))
        have_plt_inf  = []
        GSE_GPL_dict  = {}
        GSE_PMID_dict = {}
        try:
            print("平台信息检索中\n")
            for geo_id in tqdm(geo_ids):
                #print(f"Checking platform info for {geo_id}:")
                geo_id_temp, GSE_GPL_dict_temp, GSE_PMID_dict_temp = message_check.get_platform_info(geo_id)
                have_plt_inf.append(geo_id_temp)
                GSE_GPL_dict.update(GSE_GPL_dict_temp)
                GSE_PMID_dict.update(GSE_PMID_dict_temp)
        except ConnectionError:
            print("平台信息检索步骤出错，请检查您的网络连接状况\nError in the first step, please check your network connection")
        # 存储表达矩阵数据库编号的列表
        database_numbers = [number for number in have_plt_inf if len(number) > 0] 
        # 你的数据库编号列表gse_numbers = [number for number in database_numbers if 'GSE' in number] 
        # 构建FTP链接并下载文件
        message_check.matrix_download(database_numbers, current_directory, progress)#表达矩阵下载
        Both_platform_matrix_exist = message_check.matrix_check(database_numbers, current_directory, progress)
        message_output_name = message_check.summary_of_PMID_PID(Both_platform_matrix_exist, current_directory, progress, GSE_GPL_dict, GSE_PMID_dict)
        
        message_check_df = pd.read_excel(message_output_name)
        xiangxi_list = list(message_check_df['GSE'])
    
        print("执行完毕！")
        print("您的最终输出结果将保存为：日期+output.xlsx")
        
        
        
    if '2' in operation or '二' in operation:
        print ("DEG analysis")
        # 指定文件夹路径
        folder_path = this_program_path
        try:
            # 用于存储包含 "output" 的文件路径的列表
            output_files = []
            # 遍历文件夹中的所有文件和子文件夹
            for root, dirs, files in os.walk(folder_path):
                for file_name in files:
                    # 检查文件名是否包含 "output"
                    if 'output' in file_name:
                        # 构造文件的完整路径
                        file_path = os.path.join(root, file_name)
                        # 将路径添加到列表中
                        output_files.append(file_path)
            # 打印包含 "output" 的文件路径列表
            message_output_name = output_files[0]
        except:
            message_output_name = input("未找到第一步筛选所得文件，请输入路径及名称:")
        try:
            # 用于存储包含 "output" 的文件夹路径的列表
            output_folders = []
            # 遍历文件夹中的所有文件和子文件夹
            for root, dirs, files in os.walk(folder_path):
                for dir_name in dirs:
                    # 检查文件夹名是否包含 "output"
                    if 'matrix' in dir_name:
                        # 构造文件夹的完整路径
                        folder_path = os.path.join(root, dir_name)
                        # 将路径添加到列表中
                        output_folders.append(folder_path)
            matrix_list_location = folder_path
        except:
            matrix_list_location = input("未找到矩阵保存文件夹，请输入路径及名称:")
            
        message_check_df = pd.read_excel(message_output_name)

        analysis_save = current_directory+'/'+progress+'/'+'analysis'+'/'
        save_path_GPL = current_directory+'/'+progress+'/'+'platform'+'/'
        matrix_path = matrix_list_location + '/'
        
        for xiangxi in tqdm(xiangxi_list): 
            try:
                matrix_path_temp = matrix_path + f'{xiangxi}_series_matrix.txt.gz'
                my_dict, new_expset = Bio_analysis_redo.Basic_pretreat(xiangxi, matrix_path_temp, message_check_df, save_path_GPL, analysis_save)
               
                best_result, best_combination, new_expset_combination = Bio_analysis_redo.DEG_analysis(xiangxi, analysis_save, new_expset, my_dict)
                
                Bio_analysis_redo.Cluster_PCA_analysis(analysis_save, xiangxi, best_combination, new_expset_combination, new_expset)
                    
                best_result.to_excel(analysis_save + xiangxi + "/DEG/DEG_result_"+xiangxi+".xlsx",index=True)
        
            except Exception as e:
                print(e)
                print(f"缺少注释文件，将跳过本数据库:{xiangxi}")
                continue
            
    
    if '3' in operation or '三' in operation: 
        print ("GO enrichment")
        try:
            os.makedirs(analysis_save + xiangxi + '/' + "GO/")
        except:
            pass
        #sys.path.append('/Users/Rui/anaconda3/envs/Bio_analysis')
        import Bio_analysis_redo
        id_mapper, rev_id_mapper, GO_items, goeaobj = Bio_analysis_redo.GO_pre("/go-basic.obo",
                                                                               '/gene2go',
                                                                               Taxonomy_list)
        for xiangxi in tqdm(xiangxi_list):
            try:

                df_v      = pd.read_excel(analysis_save + xiangxi + "/DEG/DEG_result_"+xiangxi+".xlsx", index_col = 0)
                up_gene   = [index for index, row in df_v.iterrows() if row['color'] == 'red']
                down_gene = [index for index, row in df_v.iterrows() if row['color'] == 'blue']
                up        =Bio_analysis_redo.GO_analysis(up_gene, id_mapper, goeaobj, GO_items, rev_id_mapper)
                down      =Bio_analysis_redo.GO_analysis(down_gene, id_mapper, goeaobj, GO_items, rev_id_mapper)

            except:
                print(f'未找到差异分析文件:{xiangxi}')
                
            try:
                plt_up = Bio_analysis_redo.GO_ANALYSIS_DRAW(up)
                plt_up.savefig(analysis_save + xiangxi + "/GO/GO_up_" +xiangxi+ '.tiff',dpi=300,bbox_inches='tight')
                up.to_excel(analysis_save + xiangxi + "/GO/GO_Up_"+xiangxi+".xlsx",index=True)
            except:
                pass
            
            try:
                plt_down = Bio_analysis_redo.GO_ANALYSIS_DRAW(down)
                plt_up.savefig(analysis_save + xiangxi + "/GO/GO_down_" +xiangxi+ '.tiff',dpi=300,bbox_inches='tight')
                down.to_excel(analysis_save + xiangxi +"/GO/GO_Down_"+xiangxi+".xlsx",index=True)
            except:
                pass
            

    if '4' in operation or '四' in operation: 
        print("KEGG enrichment")
        try:
            os.makedirs(analysis_save + xiangxi + '/' + "KEGG/")
        except:
            pass
        
        sys.path.append('/Users/Rui/anaconda3/envs/Bio_analysis')
        import Bio_analysis_redo
        for xiangxi in tqdm(xiangxi_list):    
            try:

                df_v = pd.read_excel(analysis_save + xiangxi + "/DEG/DEG_result_"+xiangxi+".xlsx", index_col = 0)
                up_gene = [index for index, row in df_v.iterrows() if row['color'] == 'red']
                down_gene = [index for index, row in df_v.iterrows() if row['color'] == 'blue']
            except:
                print(f'未找到差异分析文件:{xiangxi}')
                
            try:
                plt_up, KEGG_UP_RESULT = Bio_analysis_redo.KEGG_ANALYSIS(up_gene, analysis_save, xiangxi)
                plt_up.savefig(analysis_save + xiangxi + "/KEGG/KEGG_up_" +xiangxi+ '.tiff',dpi=300,bbox_inches='tight')
                KEGG_UP_RESULT.to_excel(analysis_save + xiangxi + "/KEGG/KEGG_Up_"+xiangxi+".xlsx",index=True)
            except:
                pass
            try:
                plt_down, KEGG_DOWN_RESULT = Bio_analysis_redo.KEGG_ANALYSIS(up_gene, analysis_save, xiangxi)
                plt_up.savefig(analysis_save + xiangxi + "/KEGG/KEGG_down_" +xiangxi+ '.tiff',dpi=300,bbox_inches='tight')
                KEGG_DOWN_RESULT.to_excel(analysis_save + xiangxi +"/KEGG/KEGG_Down_"+xiangxi+".xlsx",index=True)
            except:
                pass


    if '5' in operation or '五' in operation:
        print("Subcellular localization")
        try:
            os.makedirs(analysis_save + xiangxi + '/Subcellular_Localization/')
        except:
            pass
        file_pathway = "/Users/Rui/Downloads/All RNA subcellular localization data.txt"
        subcellular_gene_dict_whole = Bio_analysis_redo.dict_subcellular_create(file_pathway)
        
        species_name = Bio_analysis_redo.get_species_name(Taxonomy_list[0])
        try:
            subcellular_gene_dict = subcellular_gene_dict_whole[species_name]
        except:
            print("物种未收录")
        for xiangxi in xiangxi_list:    
            try:
                df_v = pd.read_excel(analysis_save + xiangxi + "/DEG/DEG_result_"+xiangxi+".xlsx", index_col = 0)
                up_gene = [index for index, row in df_v.iterrows() if row['color'] == 'red']
                down_gene = [index for index, row in df_v.iterrows() if row['color'] == 'blue']
                
            except:
                print(f'未找到差异分析文件:{xiangxi}')
                continue

            try:
                df_up, plt_up = Bio_analysis_redo.analyze_gene_localization(up_gene, species_name, subcellular_gene_dict)
                df_up.to_excel(analysis_save + xiangxi +'/Subcellular_Localization/'+"/Subcellular_Up_"+xiangxi+".xlsx",index=True)
                plt_up.savefig(analysis_save + xiangxi +'/Subcellular_Localization/'+'/Subcellular_Up_Venn.tiff',dpi=300,bbox_inches='tight')
            except:
                pass
            
            try:
                df_down, plt_down = Bio_analysis_redo.analyze_gene_localization(down_gene, species_name, subcellular_gene_dict)
                df_down.to_excel(analysis_save + xiangxi +'/Subcellular_Localization/'+"/Subcellular_Down_"+xiangxi+".xlsx",index=True)
                plt_down.savefig(analysis_save + xiangxi +'/Subcellular_Localization/'+'/Subcellular_Down_Venn.tiff',dpi=300,bbox_inches='tight')        
            
            except:
                pass


    for xiangxi in xiangxi_list:
        input_file = [analysis_save + xiangxi + "/DEG/DEG_result_"  +xiangxi+ ".xlsx",
                      analysis_save + xiangxi + "/GO/GO_Down_"      +xiangxi+ ".xlsx",
                      analysis_save + xiangxi + "/GO/GO_Up_"        +xiangxi+ ".xlsx",
                      analysis_save + xiangxi + "/KEGG/KEGG_Down_"  +xiangxi+ ".xlsx",
                      analysis_save + xiangxi + "/KEGG/KEGG_Up_"    +xiangxi+ ".xlsx",
                      
                      analysis_save + xiangxi + "/Subcellular_Localization/Subcellular_Up_"   +xiangxi+ ".xlsx",
                      analysis_save + xiangxi + "/Subcellular_Localization/Subcellular_Down_" +xiangxi+ ".xlsx"]

        keywords_to_search = ['NFE2L2', 'keap1', 'MFN1', 'MFN2', 'FIS1', 'Cofilin', 'DRP']
        for i in input_file:
            try:
                Bio_analysis_redo.highlight_rows_with_keywords(i, i, keywords_to_search)
            except:
                pass
# 增加调用main()函数
if __name__ == '__main__':
    main()