#!/usr/bin/ruby

require "zlib"

# 500 0.005 out.bed

#vcfFile = Zlib::GzipReader.open(ARGV[0])

#vcfFile = ARGF

$width = ARGV[0].to_i
$af = ARGV[1].to_f

$chr = nil

$rare_count = Array.new
$all_count = Array.new
$rare_hash = Hash.new
$all_hash = Hash.new
$file = File.open(ARGV[2], 'w')

def remove(array, pos)
	
	if array.size > 0 && array.first == pos
		array.delete_at 0
	end
	
end

def reset()
	$rare_count = Array.new
	$all_count = Array.new
	$rare_hash = Hash.new
	$all_hash = Hash.new
end


def print()
	puts "print chr#{$chr}"

	$rare_count.sort!
	$all_count.sort!

	start = [$rare_count.first - $width,0].max
	stop = $rare_count.last + $width

	print_pos = [start - $width,0].max
	print_rare = 0
	print_all = 0
	
	$rare_hash.default= $rare_count.size
	$all_hash.default= $all_count.size

	(start..stop).each do |pos|
		start_rare = $rare_hash[$rare_count.bsearch {|x| x >= [pos-$width,0].max}]
		stop_rare = $rare_hash[$rare_count.bsearch {|x| x > pos+$width}]
		rare = stop_rare - start_rare

		start_all = $all_hash[$all_count.bsearch {|x| x >= [pos-$width,0].max}]
		stop_all = $all_hash[$all_count.bsearch {|x| x > pos+$width}]
		all = stop_all - start_all
	
		if print_rare != rare || print_all != all
			
			score = 0
			score = print_rare.to_f / print_all.to_f if print_rare != 0
			$file.write "chr#{$chr}\t#{print_pos-1}\t#{pos-1}\t#{score}\t#{print_rare}\t#{print_all-print_rare}\n" if print_all != 0

			print_pos = [pos,0].max
			print_rare = rare
			print_all = all

		end
		
	end
	score = 0
	score = print_rare.to_f / print_all.to_f if print_rare != 0
	$file.write "chr#{$chr}\t#{print_pos-1}\t#{stop}\t#{score}\t#{print_rare}\t#{print_all-print_rare}\n" if print_all != 0

end

i_rare = 0
i_all = 0

while STDIN.gets
	
	line = $_ 
	#puts line

	# skip headers
	next if /^#(.*)$/.match(line)

	# get values
	m = /^(chr)?([0-9XYM]+)\s+([0-9]+)\s+.+;AF=(0\.[0-9]+|0|1).+$/.match(line)
	v_chr = m[2]
	v_pos = m[3].to_i
	v_af = m[4].to_f
	
	
	# find out if new chr
	if v_chr != $chr
		print() if !$chr.nil?
		reset()
		$chr = v_chr
	end

	$all_count << v_pos
	$all_hash[v_pos] = i_all
	i_all += 1
	if v_af <= $af 
		$rare_count << v_pos 
		$rare_hash[v_pos] = i_rare
		i_rare += 1
	end

end
print()
$file.close
exit