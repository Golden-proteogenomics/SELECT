#!/usr/bin/perl

use strict;
use warnings;
my $file = "SELEX-13-LFJ14896_L3_2.fq.gz";

# 打开管道，同时执行gunzip和perl脚本
open(my $pipe, "gunzip -c $file |") or die "无法打开管道: $!";

# 初始化计数器
my $count = 0;
my $is_read = 0;

# 打开压缩输出文件
my $output_fastq = "$file.output.fq.gz";
open(my $output_fh, "| gzip > $output_fastq") or die "无法创建输出文件 '$output_fastq': $!";

# 逐行读取管道中的输入
while (my $line = <$pipe>) {
    # 判断当前行是否为读取的第一行
    if ($line =~ /^@/) {
        $is_read = 1;
    }
    
    # 如果当前行为读取的第一行，则进行处理
    if ($is_read) {
        # 输出当前行到输出文件
        print $output_fh $line;
        $count++;
    }
    
    # 如果当前读取已经包含四行，则重置 $is_read 并检查是否已提取 2000 个读取
    if ($is_read && $count % 4 == 0) {
        $is_read = 0;
        last if $count >= 100000 * 4; # 每个读取包含四行
    }
}

# 关闭管道和文件句柄
close($pipe);
close($output_fh);

print "前 10000 个 FASTQ 读取已提取并压缩到 '$output_fastq' 文件中。\n";
