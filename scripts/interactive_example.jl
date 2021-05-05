### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 1833a9ea-9d7e-400c-85aa-64a37b04ab8a
begin
	using PlutoUI
	using Plots
	using DiffDynProg
	using BioAlignments

	md"""
	#### Pairwise differentiable sequence alignment: an interactive example

	This interactive Pluto notebook contains simple, didactical and interactive
	examples for Differentiable Sequence Alignment in Julia. 
	
	Below, you can choose preset words, sequences or input your own protein sequences. 
	The Coronavirus S2 spike glycoproteins were collected from 
	Pfam (find [here](https://pfam.xfam.org/family/PF01601)). Finally, set all of 
	the parameters and let's align!
	"""
end

# ╔═╡ 3d3a3efb-198c-402a-87af-e749d950fe4b
md"""
Choose two words or sequences to align: 
$(@bind sequence1 Select(["alignment" => "alignment",
						"american" => "american",
						"gradient" => "gradient",
						"differential" => "differential",
						"certified" => "certified",
						"genetic" => "genetic",
						"grenada" => "grenada",
						"intelligence" => "intelligence",
						"differentiate" => "differentiate",
						"relative" => "relative",
"MFLILLISLPTAFAVIGDLKCTTASINDVDTGTPPISTATVDVTNGLGTYYVLDRVYLNTTLLLNGYYPISGSTYRNMALKGTLLLSTLWFKPPFLSDFNNGIFAKVKNTKVIKNGVMYSEFPAITIGSTFVNTSYSVVVLPHTTISDNKVLGLLEISVCQYTMCEYPHTICHPNLGNKRIELWHFDTDVVSCLYKRNFTYDVNADFLYFHFYQEGGTFYAYFTDTGFVTKFLFNVYLGTVLSHYYVMPLTCDSTLSLEYWVTPLTPRQYLLAFNQDGIIFNAADCMSDFMSEIKCKTQSIAPSTGVYELNGYTVQPIADVYRRIPNLPDCNIEAWLNDKSVPSPLNWERKTFSNCNFNMSSLMSFIQADSFTCNNIDAAKIYGMCFSSITIDKFAIPNGRKVDLQLGNLGYLQSFNYKIDTSATSCQLYYNLPAANVSVSRLNPSTWNRRFGFTEQSVFKPQPAGFFTAHDVVYAQHCFKAPTTFCPCKLNGSLCVGSGSGVDAGFKHTGIGTCPAGTNYLTCYNSVQCNCQCTPDPILSKATGPYKCPQTKYLVGVGEHCSGLAIKSDYCGGNPCTCQPQAFLGWSVDSCLQGDRCNIFANLILHGVNSGTTCSTDLQKANTDIIVGVCVNYDLYGISGQGIFVEVNATYYNSWQNLLYDSNGNLYGFRDYLTNRTFMIRSCYSGRVSAAFHANSSEPALLFRNIKCNYVFNNTFARQLQPINYFDSYLGCVVNADNSTASTVQTCDLTVGSGYCVDYSTQLRSRRAITTGYRFTNFEPFTVNAVNDSLQPVGGLYEIQIPSEFTIGNMEEFIQTSAPKVTIDCAAFVCGDYAACKSQLVEYGSFCDNINAILTEVNELLDTTQLQVANSLMNGVTLSTKLKDGVNFNVDDINFSPVLGCLGSDCNKASSRSAIEDLLFDKVKLADVGFVEAYNNCTGGAEIRDLICVQSYNGIKVLPPLLSENQISGYTLAATSASLFPPWTAAAGVPFYLNVQFRINGLGVTMDVLSQNQKLIANAFNNALDAIQEGFDATNSALAKIQAVVNANAEALNNLLQQLSNRFGAISSSLQEILSRLDALEAAAQIDRLINGRLTALNAYVSQQLSDSTLVKFSAAQAMEKVNECVKSQSSRINFCGNGNHIISLVQNAPYGLYFIHFSYVPTKYVTAKVSPGLCIAGDRGIAPKSGYFVNVNNTWMFTGSGYYYPEPITENNVAVMSTCAVNYTKAPDVMLNISTPNLPDFKEELDQWFKNQTSVAPDLSLGYINVTFLDIQDEMNRLQEALKVLNQSYINLKDIGTYEYYVKWPWYVWLLIGLAGVAVLVLLFFICCCTGCGTSCFKKCGGCCDDYTGHQELVIKTSHDD"=>"CoV_S2spike_var1"]))
$(@bind sequence2 Select(["alignment" => "alignment",
						"american" => "american",
						"gradient" => "gradient",
						"differential" => "differential",
						"certified" => "certified",
						"genetic" => "genetic",
						"grenada" => "grenada",
						"intelligence" => "intelligence",
						"differentiate" => "differentiate",
						"relative" => "relative",
"MIVLVTCLLLLCSYHTVLSTTNNECIQVNVTQLAGNENLIRDFLFSNFKEEGSVVVGGYYPTEVWYNCSRTARTTAFQYFNNIHAFYFVMEAMENSTGNARGKPLLFHVHGEPVSVIISAYRDDVQQRPLLKHGLVCITKNRHINYEQFTSNQWNSTCTGADRKIPFSVIPTDNGTKIYGLEWNDDFVTAYISGRSYHLNINTNWFNNVTLLYSRSSTATWEYSAAYAYQGVSNFTYYKLNNTNGLKTYELCEDYEHCTGYATNVFAPTSGGYIPDGFSFNNWFLLTNSSTFVSGRFVTNQPLLINCLWPVPSFGVAAQEFCFEGAQFSQCNGVSLNNTVDVIRFNLNFTADVQSGMGATVFSLNTTGGVILEISCYSDTVSESSSYSYGEIPFGITDGPRYCYVLYNGTALKYLGTLPPSVKEIAISKWGHFYINGYNFFSTFPIGCISFNLTTGVSGAFWTIAYTSYTEALVQVENTAIKNVTYCNSHINNIKCSQLTANLNNGFYPVASSEVGFVNKSVVLLPSFFTYTAVNITIDLGMKLSGYGQPIASTLSNITLPMQDNNTDVYCIRSNQFSVYVHSTCKSSLWDNIFNQDCTDVLEATAVIKTGTCPFSFDKLNNYLTFNKFCLSLSPVGANCKFDVAARTRTNEQVVRSLYVIYEEGDNIVGVPSDNSGLHDLSVLHLDSCTDYNIYGRTGVGIIRRTNSTLLSGLYYTSLSGDLLGFKNVSDGVIYSVTPCDVSAQAAVIDGAIVGAMTSINSELLGLTHWTTTPNFYYYSIYNYTSERTRGTAIDSNDVDCEPVITYSNIGVCKNGALVFINVTHSDGDVQPISTGNVTIPTNFTISVQVEYMQVYTTPVSIDCARYVCNGNPRCNKLLTQYVSACQTIEQALAMGARLENMEVDSMLFVSENALKLASVEAFNSTENLDPIYKEWPSIGGSWLGGLKDILPSHNSKRKYGSAIEDLLFDKVVTSGLGTVDEDYKRCTGGYDIADLVCAQYYNGIMVLPGVANADKMTMYTASLAGGITLGALGGGAVAIPFAVAVQARLNYVALQTDVLNKNQQILANAFNQAIGNITQAFGKVNDAIHQTSQGLATVAKALAKVQDVVNTQGQALSHLTVQLQNNFQAISSSISDIYNRLDELSADAQVDRLITGRLTALNAFVSQTLTRQAEVRASRQLAKDKVNECVRSQSQRFGFCGNGTHLFSLANAAPNGMIFFHTVLLPTAYETVTAWSGICASDGDRTFGLVVKDVQLTLFRNLDDKFYLTPRTMYQPRVATSSDFVQIEGCDVLFVNATVIDLPSIIPDYIDINQTVQDILENYRPNWTVPEFTLDIFNATYLNLTGEIDDLEFRSEKLHNTTVELAILIDNINNTLVNLEWLNRIETYVKWPWYVWLLIGLVVVFCIPLLLFCCFSTGCCGCIGCLGSCCHSICSRRQFENYEPIEKVHVH"=>"CoV_S2spike_var2"]))

Alternatively, supply your own protein sequences:
$(@bind protein1 TextField())
$(@bind protein2 TextField())

What type of alignment do you want to do?
$(@bind alignment_type Select(["global alignment", "local alignment"]))

Which max operator do you want? 
$(@bind maxop Select(["EntropyMax", "SquaredMax"], default="EntropyMax"))

Set the value for the smoothing parameter:
$(@bind lambda Slider(0.1:15, default=1; show_value=true))

Set the gap costs (by default, they are equal):
$(@bind gaps Slider(1:10, default=1; show_value=true))
$(@bind gapt Slider(1:10, default=1; show_value=true))

What substitution matrix do you want?
$(@bind submat Select(["BLOSUM45", "BLOSUM62", "BLOSUM90", "0/1 scoring"], default="BLOSUM62"))

Which plot(s) do you want to see?
$(@bind plot_type Select(["Overview", "Cost matrix θ", "DP matrix", "Overall gradient", "Gradient of θ", "Gradients of gap costs"], default="Overall gradient"))
"""

# ╔═╡ a8c1d51e-ac4f-4dfd-850b-61330c8abcad
begin
	((length(protein1) > 0) && (length(protein2) > 0) ? 
	zoomlimx = min(length(protein1)) :
	zoomlimx = min(length(sequence1)));
		
	((length(protein1) > 0) && (length(protein2) > 0) ? 
	zoomlimy = min(length(protein2)) :
	zoomlimy = min(length(sequence2)));
	
	md"""
	Zoom range for the right plot if applicable (x-axis & y-axis, respectively):
	
	$(@bind xrange RangeSlider(1:zoomlimx; show_value=true))
	
	$(@bind yrange RangeSlider(1:zoomlimy; show_value=true))
	"""
end

# ╔═╡ c31ccdbc-5b02-4287-8abb-defcd3c8b256
begin
	# select correct sequences: preset or own sequences
	if (length(protein1) > 0) && (length(protein2) > 0)
		s1 = protein1
		s2 = protein2
	else
		s1 = sequence1
		s2 = sequence2
	end
		
	# set the correct substitution matrix and compute theta
	if submat == "0/1 scoring"
		θ = [sᵢ == tⱼ for sᵢ in s1, tⱼ in s2] * 20.0 .- 5  # 0-1
	else
		S = (submat == "BLOSUM45") ? BLOSUM45 : 
		((submat == "BLOSUM62") ? BLOSUM62 : BLOSUM90)
		θ = [float(S[sᵢ,tⱼ]) for sᵢ in s1, tⱼ in s2] # blossum
	end
	
	# compute the alignment
	n, m = length(s1), length(s2)
	mo = (maxop == "EntropyMax") ? EntropyMax(lambda) : SquaredMax(lambda)
	if alignment_type == "global alignment"
		dp = DP(θ)
		D, E, Eθ, Egs, Egt = ∂NW_all(mo, θ, (gaps, gapt), dp)
		D = getD(dp)
		M = DiffDynProg.softmax(D)
	else
		dp = DP(θ)
		D, E, Eθ, Egs, Egt = ∂SW_all(mo, θ, (gaps, gapt), dp)
		D = getD(dp)
		M = DiffDynProg.softmax(D)
	end
	
	# compute the plot
	yticks = (1:n, split(s1, ""));
	xticks = (1:m, split(s2, ""));
	xstart = xrange[1]; xstop = xrange[end]
	ystart = yrange[1]; ystop = yrange[end]
	
	if plot_type == "Overview"
		l = @layout [a b; c d]
		p1 = heatmap(θ; yticks, xticks, yflip=true, color=:vik, clims=(-6.0, 6.0),
            axes=false)
		p2 = plot(heatmap(D; yticks, xticks, yflip=true, color=:Blues, axes=false))
		p3 = plot(heatmap(E; yticks, xticks, yflip=true, 
				color=cgrad([:white, :green]), clims=(0,1), axes=false))
		p4 = plot(heatmap(M; yticks, xticks, yflip=true, 
				color=cgrad([:white, :purple]), clims=(0,1), axes=false))
		plot(p1, p2, p3, p4, title=["θ" "D" "E" "M"], layout = l, size=(1200, 1000))
		
	elseif plot_type == "Cost matrix θ"
		l = @layout [a b]
		p1 = heatmap(θ; 
			yticks, xticks, yflip=true, color=:vik, clims=(-6.0, 6.0), axes=false)
		p2 = heatmap(θ[xstart:xstop, ystart:ystop];
			yticks, xticks, yflip=true, color=:vik, clims=(-6.0, 6.0), axes=false)
		plot(p1, p2, title=["θ" "θ zoomed"], layout = l, size=(1200, 500))
		
	elseif plot_type == "DP matrix" 
		l = @layout [a b]
		p1 = plot(heatmap(D;
				yticks, xticks, yflip=true, color=:Blues, axes=false))
		p2 = plot(heatmap(D[xstart:xstop, ystart:ystop];
				yticks, xticks, yflip=true, color=:Blues, axes=false))
		plot(p1, p2, title=["D" "D zoomed"], layout = l, size=(1200, 500))
		
	elseif plot_type == "Overall gradient"
		l = @layout [a b]
		p1 = plot(heatmap(E; 
				yticks, xticks, yflip=true, color=cgrad([:white, :green]), 
				clims=(0,1), axes=false))
		p2 = plot(heatmap(E[xstart:xstop, ystart:ystop]; 
				yticks, xticks, yflip=true, color=cgrad([:white, :green]), 
				clims=(0,1), axes=false))
		plot(p1, p2, title=["E" "E zoomed"], layout = l, size=(1200, 500))
		
	elseif plot_type == "Gradient of θ"
		l = @layout [a b]
		p1 = plot(heatmap(Eθ; 
				yticks, xticks, yflip=true, color=cgrad([:white, :green]), 
				clims=(0,1), axes=false))
		p2 = plot(heatmap(Eθ[xstart:xstop, ystart:ystop]; 
				yticks, xticks, yflip=true, color=cgrad([:white, :green]), 
				clims=(0,1), axes=false))
		plot(p1, p2, title=["Eθ" "Eθ zoomed"], layout = l, size=(1200, 500))
		
	elseif plot_type == "Gradients of gap costs"
		l = @layout [a b]
		p1 = plot(heatmap(Egs; 
				yticks, xticks, yflip=true, color=cgrad([:white, :green]), 
				clims=(0,1), axes=false))
		p2 = plot(heatmap(Egt; 
				yticks, xticks, yflip=true, color=cgrad([:white, :green]), 
				clims=(0,1), axes=false))
		plot(p1, p2, title=["Egs" "Egt"], layout = l, size=(1200, 500))
	end
end

# ╔═╡ Cell order:
# ╟─1833a9ea-9d7e-400c-85aa-64a37b04ab8a
# ╟─3d3a3efb-198c-402a-87af-e749d950fe4b
# ╟─a8c1d51e-ac4f-4dfd-850b-61330c8abcad
# ╟─c31ccdbc-5b02-4287-8abb-defcd3c8b256
