### A Pluto.jl notebook ###
# v0.14.5

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

	This interactive notebook illustrates Differentiable Sequence Alignment in Julia through interactive examples. 
	
	Below, you can choose preset words, sequences or input your own protein sequences. 
	The Coronavirus S2 spike glycoproteins were collected from 
	Pfam (find [here](https://pfam.xfam.org/family/PF01601)). Finally, set all of 
	the parameters and let's align!
	"""
end

# ╔═╡ 56d7d039-ca74-49ff-ab99-a642e75bf2ef
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
$(@bind seq_free1 TextField())
$(@bind seq_free2 TextField())
"""

# ╔═╡ 66607a15-6aab-408a-b30a-4b038706b447
begin
	s = (length(seq_free1) > 0) ? seq_free1 : sequence1;
	t = (length(seq_free2) > 0) ? seq_free2 : sequence1;
end

# ╔═╡ 3bb2cfd3-e30e-4276-b70d-070b3de0a033
n, m = length(s), length(t);

# ╔═╡ b36ef14e-088c-483f-9dbe-edbbd951830a
md"""
What type of alignment do you want to do?
$(@bind alignment_type Select(["global alignment", "local alignment"]))
"""

# ╔═╡ 0a9bccaa-72b2-4394-b561-1bdcb543a11d
md"""
Which max operator do you want to use? 
$(@bind maxop Select(["EntropyMax", "SquaredMax", "Hard"], default="EntropyMax"))
"""

# ╔═╡ b8ad5efe-f723-4c51-b8a8-fff26d67dd24
maxop  # check name

# ╔═╡ 77546a56-3b11-466e-a5ad-8bef6115ac97
md"""
Set the value for the smoothing parameter $\lambda$:
$(@bind lambda Slider(0.1:15, default=1; show_value=true))
"""

# ╔═╡ 3cab7e29-4750-4ca9-b3c1-910e16d00eaf
if (maxop == "EntropyMax")
	mo = EntropyMax(lambda)
elseif (maxop == "SquaredMax")
	mo = SquaredMax(lambda)
else
	mo = Max()
end

# ╔═╡ 419eea6b-0ace-406f-8977-1ecfc2c18747
md"""
Set the gap costs (by default, they are equal):
$(@bind gs Slider(0:10, default=1; show_value=true))
$(@bind gt Slider(0:10, default=1; show_value=true))
"""

# ╔═╡ ebc96485-9c73-4307-9799-531a1abdd7d4
md"""
Which substitution matrix do you to use?
$(@bind submat Select(["BLOSUM45", "BLOSUM62", "BLOSUM90", "0/1 scoring"], default="BLOSUM62"))
"""

# ╔═╡ 4f5709d1-9a78-4fbf-9e68-e7c7c7e85f26
begin
	# set the correct substitution matrix and compute theta
	if submat == "0/1 scoring"
		θ = [sᵢ == tⱼ for sᵢ in s, tⱼ in t]# * 20.0 .- 5  # 0-1
	else
		S = (submat == "BLOSUM45") ? BLOSUM45 : 
		((submat == "BLOSUM62") ? BLOSUM62 : BLOSUM90)
		θ = [float(S[sᵢ,tⱼ]) for sᵢ in s, tⱼ in t] # blossum
	end
end;

# ╔═╡ df424f8f-c454-4c0a-a3b6-cad807df77cf
begin 
	heatmap(θ; yflip=true, color=:vik, clims=(-6.0, 6.0), axes=false)
end

# ╔═╡ 352c789a-cac5-410a-81af-e1eee39797fc
∂NW_all(mo, θ, (gs, gt))

# ╔═╡ ce115b70-e2dc-460e-b8ca-6c1e9c7ea2b1
if alignment_type == "global alignment"
	∂NW_all(mo, θ, (gs, t))

# ╔═╡ 3d3a3efb-198c-402a-87af-e749d950fe4b
md"""



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

# ╔═╡ feeac491-ad6d-48d3-95d2-71c5cbfc09ec
md"Show sequence on plots: $(@bind show_seq CheckBox(default=true))"

# ╔═╡ c31ccdbc-5b02-4287-8abb-defcd3c8b256
begin

	
	
	# compute the alignment
	
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
# ╠═1833a9ea-9d7e-400c-85aa-64a37b04ab8a
# ╠═56d7d039-ca74-49ff-ab99-a642e75bf2ef
# ╠═66607a15-6aab-408a-b30a-4b038706b447
# ╠═3bb2cfd3-e30e-4276-b70d-070b3de0a033
# ╠═b36ef14e-088c-483f-9dbe-edbbd951830a
# ╠═0a9bccaa-72b2-4394-b561-1bdcb543a11d
# ╠═b8ad5efe-f723-4c51-b8a8-fff26d67dd24
# ╠═77546a56-3b11-466e-a5ad-8bef6115ac97
# ╠═3cab7e29-4750-4ca9-b3c1-910e16d00eaf
# ╠═419eea6b-0ace-406f-8977-1ecfc2c18747
# ╠═ebc96485-9c73-4307-9799-531a1abdd7d4
# ╠═4f5709d1-9a78-4fbf-9e68-e7c7c7e85f26
# ╠═df424f8f-c454-4c0a-a3b6-cad807df77cf
# ╠═352c789a-cac5-410a-81af-e1eee39797fc
# ╠═ce115b70-e2dc-460e-b8ca-6c1e9c7ea2b1
# ╠═3d3a3efb-198c-402a-87af-e749d950fe4b
# ╠═a8c1d51e-ac4f-4dfd-850b-61330c8abcad
# ╠═feeac491-ad6d-48d3-95d2-71c5cbfc09ec
# ╠═c31ccdbc-5b02-4287-8abb-defcd3c8b256
