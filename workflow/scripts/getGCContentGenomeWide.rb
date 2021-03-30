
def remove(array, pos)
	
	if array.size > 0 && array.first == pos
		array.delete_at 0
	end
	
end


fastaFile = File.open(ARGV[0], 'r')

width = ARGV[1].to_i

chr = nil

gc_count = Array.new
n_count = Array.new
chars = Array.new

position = -1

fastaFile.each_line do |line|

	if m = /^>(chr.+)/.match(line)
		
		if !chr.nil?
			if correct_chr = /^(chr[0-9XYM]+)$/.match(chr)
				new_chr=correct_chr[1]
				(position-width+1..position).each do |p|

					if chars.first != 'N' && chars.first != 'n'
						gc = gc_count.size
						size = [[p-width,0].max,width].min + (position+width-p)-1 - n_count.size
						puts "#{new_chr}\t#{p}\t#{p+1}\t#{(gc.to_f/size.to_f).round(6)}"
					end

					remove(gc_count,p-width)
					remove(n_count,p-width)
					chars.delete_at(0)
				end
			end
			
		end
		gc_count = Array.new
		chars = Array.new
		n_count = Array.new
		chr = m[1].split(/\s+/)[0]
		position = -1
		next
	end

	line.strip!
	line.each_char do |c|
		position += 1
		chars << c
		gc_count << position if (c == 'G' || c == 'C' || c == 'g' || c == 'c')
		n_count << position if (c == 'N' || c == 'n')
	
		next if position - width < 0
		
		
		if chars.first != 'N' && chars.first != 'n'
			if correct_chr = /^(chr[0-9XYM]+)$/.match(chr)
				new_chr=correct_chr[1]
				gc = gc_count.size
				size = [[position-width,0].max,width].min + width + 1 - n_count.size
				puts "#{new_chr}\t#{position-width}\t#{position-width+1}\t#{(gc.to_f/size.to_f).round(6)}"
			end
		end
		
		remove(gc_count,position-2*width)
		remove(n_count,position-2*width)
		chars.delete_at(0)

	end

end

(position-width+1..position).each do |p|
	if chars.first != 'N' && chars.first != 'n'
		if correct_chr = /^(chr[0-9XYM]+)$/.match(chr)
			new_chr=correct_chr[1]
			gc = gc_count.size
			size = [[p-width,0].max,width].min + (position+width-p)-1 - n_count.size
			puts "#{new_chr}\t#{p}\t#{p+1}\t#{(gc.to_f/size.to_f).round(6)}"
		end
	end
		
	remove(gc_count,p-width)
	remove(n_count,p-width)
	chars.delete_at(0)
end
