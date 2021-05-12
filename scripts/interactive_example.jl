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
	## Pairwise differentiable sequence alignment: an interactive example

	This interactive notebook illustrates Differentiable Sequence Alignment in Julia through interactive examples. 
	
	Below, you can choose preset words, sequences or input your own protein sequences. 
	We have included real-world examples of two Coronavirus S2 spike glycoproteins, two zinc-finger domains (C2H2 type) and two fibronectin type III domains (collected from Pfam [here](https://pfam.xfam.org/family/PF01601), [here](https://pfam.xfam.org/family/PF00096) and [here](https://pfam.xfam.org/family/PF00041)). Afterwards choosing two sequences, set all of the parameters and let's align!
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
"MFLILLISLPTAFAVIGDLKCTTASINDVDTGTPPISTATVDVTNGLGTYYVLDRVYLNTTLLLNGYYPISGSTYRNMALKGTLLLSTLWFKPPFLSDFNNGIFAKVKNTKVIKNGVMYSEFPAITIGSTFVNTSYSVVVLPHTTISDNKVLGLLEISVCQYTMCEYPHTICHPNLGNKRIELWHFDTDVVSCLYKRNFTYDVNADFLYFHFYQEGGTFYAYFTDTGFVTKFLFNVYLGTVLSHYYVMPLTCDSTLSLEYWVTPLTPRQYLLAFNQDGIIFNAADCMSDFMSEIKCKTQSIAPSTGVYELNGYTVQPIADVYRRIPNLPDCNIEAWLNDKSVPSPLNWERKTFSNCNFNMSSLMSFIQADSFTCNNIDAAKIYGMCFSSITIDKFAIPNGRKVDLQLGNLGYLQSFNYKIDTSATSCQLYYNLPAANVSVSRLNPSTWNRRFGFTEQSVFKPQPAGFFTAHDVVYAQHCFKAPTTFCPCKLNGSLCVGSGSGVDAGFKHTGIGTCPAGTNYLTCYNSVQCNCQCTPDPILSKATGPYKCPQTKYLVGVGEHCSGLAIKSDYCGGNPCTCQPQAFLGWSVDSCLQGDRCNIFANLILHGVNSGTTCSTDLQKANTDIIVGVCVNYDLYGISGQGIFVEVNATYYNSWQNLLYDSNGNLYGFRDYLTNRTFMIRSCYSGRVSAAFHANSSEPALLFRNIKCNYVFNNTFARQLQPINYFDSYLGCVVNADNSTASTVQTCDLTVGSGYCVDYSTQLRSRRAITTGYRFTNFEPFTVNAVNDSLQPVGGLYEIQIPSEFTIGNMEEFIQTSAPKVTIDCAAFVCGDYAACKSQLVEYGSFCDNINAILTEVNELLDTTQLQVANSLMNGVTLSTKLKDGVNFNVDDINFSPVLGCLGSDCNKASSRSAIEDLLFDKVKLADVGFVEAYNNCTGGAEIRDLICVQSYNGIKVLPPLLSENQISGYTLAATSASLFPPWTAAAGVPFYLNVQFRINGLGVTMDVLSQNQKLIANAFNNALDAIQEGFDATNSALAKIQAVVNANAEALNNLLQQLSNRFGAISSSLQEILSRLDALEAAAQIDRLINGRLTALNAYVSQQLSDSTLVKFSAAQAMEKVNECVKSQSSRINFCGNGNHIISLVQNAPYGLYFIHFSYVPTKYVTAKVSPGLCIAGDRGIAPKSGYFVNVNNTWMFTGSGYYYPEPITENNVAVMSTCAVNYTKAPDVMLNISTPNLPDFKEELDQWFKNQTSVAPDLSLGYINVTFLDIQDEMNRLQEALKVLNQSYINLKDIGTYEYYVKWPWYVWLLIGLAGVAVLVLLFFICCCTGCGTSCFKKCGGCCDDYTGHQELVIKTSHDD" => "CoV_S2spike_var1", "FRCSECSRSFTHNSDLTAHMRKH" => "ZnFingerXFIN_XENLA", "DAPSNLRFLATTPNSLLVSWQPPRARITGYIIKYEKPGSPPREVVPRPRPGVTEATITGLEPGTEYTIQVIALKNNQKS" => "fibronectinFINC_BOVIN"]))
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
"MIVLVTCLLLLCSYHTVLSTTNNECIQVNVTQLAGNENLIRDFLFSNFKEEGSVVVGGYYPTEVWYNCSRTARTTAFQYFNNIHAFYFVMEAMENSTGNARGKPLLFHVHGEPVSVIISAYRDDVQQRPLLKHGLVCITKNRHINYEQFTSNQWNSTCTGADRKIPFSVIPTDNGTKIYGLEWNDDFVTAYISGRSYHLNINTNWFNNVTLLYSRSSTATWEYSAAYAYQGVSNFTYYKLNNTNGLKTYELCEDYEHCTGYATNVFAPTSGGYIPDGFSFNNWFLLTNSSTFVSGRFVTNQPLLINCLWPVPSFGVAAQEFCFEGAQFSQCNGVSLNNTVDVIRFNLNFTADVQSGMGATVFSLNTTGGVILEISCYSDTVSESSSYSYGEIPFGITDGPRYCYVLYNGTALKYLGTLPPSVKEIAISKWGHFYINGYNFFSTFPIGCISFNLTTGVSGAFWTIAYTSYTEALVQVENTAIKNVTYCNSHINNIKCSQLTANLNNGFYPVASSEVGFVNKSVVLLPSFFTYTAVNITIDLGMKLSGYGQPIASTLSNITLPMQDNNTDVYCIRSNQFSVYVHSTCKSSLWDNIFNQDCTDVLEATAVIKTGTCPFSFDKLNNYLTFNKFCLSLSPVGANCKFDVAARTRTNEQVVRSLYVIYEEGDNIVGVPSDNSGLHDLSVLHLDSCTDYNIYGRTGVGIIRRTNSTLLSGLYYTSLSGDLLGFKNVSDGVIYSVTPCDVSAQAAVIDGAIVGAMTSINSELLGLTHWTTTPNFYYYSIYNYTSERTRGTAIDSNDVDCEPVITYSNIGVCKNGALVFINVTHSDGDVQPISTGNVTIPTNFTISVQVEYMQVYTTPVSIDCARYVCNGNPRCNKLLTQYVSACQTIEQALAMGARLENMEVDSMLFVSENALKLASVEAFNSTENLDPIYKEWPSIGGSWLGGLKDILPSHNSKRKYGSAIEDLLFDKVVTSGLGTVDEDYKRCTGGYDIADLVCAQYYNGIMVLPGVANADKMTMYTASLAGGITLGALGGGAVAIPFAVAVQARLNYVALQTDVLNKNQQILANAFNQAIGNITQAFGKVNDAIHQTSQGLATVAKALAKVQDVVNTQGQALSHLTVQLQNNFQAISSSISDIYNRLDELSADAQVDRLITGRLTALNAFVSQTLTRQAEVRASRQLAKDKVNECVRSQSQRFGFCGNGTHLFSLANAAPNGMIFFHTVLLPTAYETVTAWSGICASDGDRTFGLVVKDVQLTLFRNLDDKFYLTPRTMYQPRVATSSDFVQIEGCDVLFVNATVIDLPSIIPDYIDINQTVQDILENYRPNWTVPEFTLDIFNATYLNLTGEIDDLEFRSEKLHNTTVELAILIDNINNTLVNLEWLNRIETYVKWPWYVWLLIGLVVVFCIPLLLFCCFSTGCCGCIGCLGSCCHSICSRRQFENYEPIEKVHVH"=>"CoV_S2spike_var2", "YTCGYCTEDSPSFPRPSLLESHISLMH" => "ZnFingerZN592", "DKVQGVSVSNSARSDYLRVSWVHATGDFDHYEVTIKNKNNFIQTKSIPKSENECVFVQLVPGRLYSVTVTTKSGQYE" => "fibronectinPTPRB_HUMAN"], default="gradient"))
"""

# ╔═╡ 52d9e9ec-c86d-496c-b47c-f5fd7b4a6af3
md"""
Alternatively, supply your own protein sequences:
$(@bind seq_free1 TextField())
$(@bind seq_free2 TextField())
"""

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

# ╔═╡ 77546a56-3b11-466e-a5ad-8bef6115ac97
md"""
Set the value for the smoothing parameter $\lambda$:
$(@bind lambda Slider(0.1:15, default=1; show_value=true))
"""

# ╔═╡ 419eea6b-0ace-406f-8977-1ecfc2c18747
md"""
Set the gap costs (by default, they are equal):
$(@bind gs Slider(0:10, default=1; show_value=true))
$(@bind gt Slider(0:10, default=1; show_value=true))
"""

# ╔═╡ ebc96485-9c73-4307-9799-531a1abdd7d4
md"""
Which substitution matrix do you want to use?
$(@bind submat Select(["BLOSUM45", "BLOSUM62", "BLOSUM90", "0/1 scoring"], default="BLOSUM62"))
"""

# ╔═╡ feeac491-ad6d-48d3-95d2-71c5cbfc09ec
md"Show sequence on plots: $(@bind show_seq CheckBox(default=true))"

# ╔═╡ d4a8ec90-0653-47c3-aa86-30ef0ae298fa
md"""
### Cost matrix $\theta$

The matrix $\theta$ quantifies the (dis)similarity between the amino acid residues of each of the two sequences. The matrix can be constructed based on the widely used BLOSUM substitution matrices or simple 0/1 scoring. The current matrix plotted below is constructed using $submat (can be changed in the parameter setting above).

If a subregion was selected using the sliders above, a second plot of this subregion appears below the complete plot.
"""

# ╔═╡ 67fe4cd2-f3cc-48be-825a-4fe493e19687
md"""
Values of the matrix (or subregion):
"""

# ╔═╡ 39ca6f98-0dac-429b-954b-f04c1e703db0
md"""
## Dynamic programming matrix

Below is the dynamic programming matrix $D$, which is computed based on the cost matrix $\theta$.

Again, if a subregion was selected using the sliders above, a second plot of this subregion appears below the complete plot.
"""

# ╔═╡ f5aa16a6-ad76-48b2-af88-2b0439b56a21
md""" 
## Gradients of the dynamic programming matrix

The global gradient of matrix $D$ is given below. If a subregion was selected, this subregion is also plotted underneath.
"""

# ╔═╡ f0c5d254-8f0f-41eb-b3b5-56c49cfac65c
md"""
This gradient can be decomposed into three components:
- a component depending on extending the alignment (match or mismatch), which is influenced by the substitution matrix;
- a component depending on inserting gaps in $s$;
- a component depending on inserting gaps in $t$.
"""

# ╔═╡ aefbe8d1-528e-46e0-aa45-b3294f8ef54d
md"""
We can also look at the gradient components in alternative ways. The plots below illustrate how the different substitutions between residues influence the alignemt, as well how the gradients related to the gaps change along the sequence (i.e. at what position is more likely to be introduced for each sequence).
"""

# ╔═╡ b5252d6b-db49-4e38-9fb7-2c4fedabf2e9
# Setting the correct sequences for alignment
begin
	s = (length(seq_free1) > 0) ? seq_free1 : sequence1;
	t = (length(seq_free2) > 0) ? seq_free2 : sequence2;
	n, m = length(s), length(t);
end;

# ╔═╡ e7735020-55f9-4dd2-bf78-67d60451aab3
begin
	
	md"""
	Optionally select a subregion of the alignment (zoom) to plot in detail (x-axis & y-axis, respectively):
	
	$(@bind yrange RangeSlider(1:n; show_value=true))
	
	$(@bind xrange RangeSlider(1:m; show_value=true))
	"""
end

# ╔═╡ 4ed12b25-4aca-41aa-93c9-44496c831e5e
# defining the subrange if applicable
begin
	subrange = !(extrema(yrange) == (1, n) && extrema(xrange) == (1, m))
	subregion = ([minimum(xrange), maximum(xrange), maximum(xrange), minimum(xrange), minimum(xrange)],
			[minimum(yrange), minimum(yrange), maximum(yrange), maximum(yrange), minimum(yrange)])
	xticks_sub = 1:length(xrange), split(t, "")[xrange]
	yticks_sub = 1:length(yrange), split(s, "")[yrange]
end;

# ╔═╡ dd3f27f5-0513-422d-a9b9-1da41c88b38c
# set the correct substitution matrix and compute theta
begin
	if submat == "0/1 scoring"
		θ = [sᵢ == tⱼ for sᵢ in s, tⱼ in t] * 20.0 .- 5.0  # 0-1
	else
		S = (submat == "BLOSUM45") ? BLOSUM45 : 
		((submat == "BLOSUM62") ? BLOSUM62 : BLOSUM90)
		θ = [float(S[sᵢ,tⱼ]) for sᵢ in s, tⱼ in t] # blossum
	end
end;

# ╔═╡ df424f8f-c454-4c0a-a3b6-cad807df77cf
begin 
	pθ = heatmap(θ; yflip=true, color=:vik, clims=(-6.0, 6.0), axes=false, title="Full matrix",  titlefontsize=11)
	show_seq && yticks!((1:n, split(s, "")))
	show_seq && xticks!((1:m, split(t, "")))
	if subrange
		plot!(pθ, subregion..., color="orange", lw=2, label="")
	end
	pθ
end

# ╔═╡ c220f01e-ecd5-4b17-be4c-4ebeb2062277
if subrange
	pθs = heatmap(θ[yrange, xrange]; yflip=true, color=:vik, clims=(-6.0, 6.0), axes=false, title="Subregion of matrix",  titlefontsize=11)
	show_seq && yticks!(yticks_sub...)
	show_seq && xticks!(xticks_sub...)
	pθs
end

# ╔═╡ 0251ca0d-f3a4-4758-a451-3b158b187992
θ[yrange, xrange]

# ╔═╡ 5ec50f03-9c0c-4dd6-b8dc-cfb17d9b32b0
# setting the correct max operator
if (maxop == "EntropyMax")
	mo = EntropyMax(lambda)
elseif (maxop == "SquaredMax")
	mo = SquaredMax(lambda)
else
	mo = Max()
end;

# ╔═╡ ce115b70-e2dc-460e-b8ca-6c1e9c7ea2b1
# Computing the differentiable alignment
begin 
	if alignment_type == "global alignment"
		D, E, Eθ, Egs, Egt = ∂NW_all(mo, θ, (gs, gt))
		v = last(D)
	else
		D, E, Eθ, Egs, Egt = ∂SW_all(mo, θ, (gs, gt))
		M = DiffDynProg.softmax(D)
		v = DiffDynProg.logsumexp(D)
	end
end;

# ╔═╡ 0a73991d-1ae7-4908-94ab-c0cbbf25a9c5
begin 
	pD = heatmap(D; yflip=true, color=:deep, axes=false, title="Full matrix",  titlefontsize=11)
	show_seq && yticks!((1:n, split(s, "")))
	show_seq && xticks!((1:m, split(t, "")))
	if subrange
		plot!(pD, subregion..., color="orange", lw=2, label="")
	end
	pD
end

# ╔═╡ 748e9149-800d-464f-87ea-6b51438baac0
if subrange
	pDs = heatmap(D[yrange,xrange]; yflip=true, color=:deep, axes=false, title="Subregion of matrix", titlefontsize=11)
	show_seq && yticks!(yticks_sub...)
	show_seq && xticks!(xticks_sub...)
	pDs
end

# ╔═╡ 46923e27-adc2-48a8-9ca0-979c3d32a5f3
if alignment_type == "global alignment"
	md"For global alignment, the alignment score is the final value in this matrix (bottom right). Here this value is $v."
else
	md"In local alignment, the alignment score is given by log-sum-exp of the matrix $D$. Here, this value is $v."
end

# ╔═╡ 6e45cee4-466d-415b-85c9-c5a93be88d4e
begin
	if alignment_type != "global alignment"
		md"""
		The change of this alignment score w.r.t. the elements of $D$ is given below.
		If a subregion is selected, the subregion is also plotted underneath.
		"""

		pM = plot(heatmap(M; yflip=true, 
					color=cgrad([:white, :purple]), clims=(0,1), axes=false), title="Full matrix", titlefontsize=11)
		show_seq && yticks!((1:n, split(s, "")))
		show_seq && xticks!((1:m, split(t, "")))
		if subrange
			plot!(pM, subregion..., color="orange", lw=2, label="")
		end
		pM
	end
end

# ╔═╡ 1f48c0da-a872-4674-85ed-21de2d4fe6ab
if alignment_type != "global alignment"
	if subrange
		pMs = plot(heatmap(M[yrange,xrange]; yflip=true, 
					color=cgrad([:white, :purple]), clims=(0,1), axes=false), title="Subregion of matrix",  titlefontsize=11)
		show_seq && yticks!(yticks_sub...)
		show_seq && xticks!(xticks_sub...)
		pMs
	end
end

# ╔═╡ 29f3cf02-d7b3-470a-9a4d-8a988e41f366
begin
	pE = plot(heatmap(E, yflip=true, 
				color=cgrad([:white, :green]), clims=(0,1), axes=false), title="Full matrix",  titlefontsize=11)
	show_seq && yticks!((1:n, split(s, "")))
	show_seq && xticks!((1:m, split(t, "")))
	if subrange
		plot!(pE, subregion..., color="orange", lw=2, label="")
	end
	pE
end

# ╔═╡ 3bc02cf8-62c2-4cf2-977b-75650e0c6cf4
if subrange
	pEs = plot(heatmap(E[yrange, xrange], yflip=true, 
				color=cgrad([:white, :green]), clims=(0,1), axes=false), title="Subregion of matrix",  titlefontsize=11)
	show_seq && yticks!(yticks_sub...)
	show_seq && xticks!(xticks_sub...)
	pEs
end

# ╔═╡ 6601ac61-dd7e-474f-8c0b-b0e81cb1b5fa
begin
	pEcomb = plot(
		heatmap(Eθ, yflip=true, 
					color=cgrad([:white, :green]), clims=(0,1), axes=false,
			title="component driven by alignment extension", titlefontsize=11),
		heatmap(Egs, yflip=true, 
					color=cgrad([:white, :green]), clims=(0,1), axes=false,
			title="component driven by gap in s", titlefontsize=11),
		heatmap(Egt, yflip=true, 
					color=cgrad([:white, :green]), clims=(0,1), axes=false,
			title="component driven by gap in t", titlefontsize=11),
		layout=(3, 1), size=(800, 1200))
	show_seq && yticks!((1:n, split(s, "")))
	show_seq && xticks!((1:m, split(t, "")))
	pEcomb
end

# ╔═╡ 85abc158-ce84-4eff-9d01-3a75d5483e1b
begin
	# compute gradient substitution matrix
	chars_s = s |> unique |> sort!
	chars_t = t |> unique |> sort!
	dS = zeros(length(chars_s), length(chars_t))
	for (i, sᵢ) in enumerate(s)
		for (j, tⱼ) in enumerate(t)
			dS[findfirst(isequal(sᵢ), chars_s), findfirst(isequal(tⱼ), chars_t)] += Eθ[i,j]
		end
	end
	
	# compute sum of gap-related components
	dgs = sum(Egs, dims=2)[:]
	dgt = sum(Egt, dims=1)[:]
	
	# make plot
	pS = plot(
		heatmap(dS, yflip=true, yticks=(1:length(chars_s), chars_s),
				xticks=(1:length(chars_t), chars_t),
					color=cgrad([:white, :green]), axes=false,
			title="gradient driven by substitution matrix", titlefontsize=12),
		
		plot(dgs, title="gradient driven by gap in s", label="", color=:green,
			xlabel="position in s", titlefontsize=12, xguidefontsize=10),
		plot(dgt, title="gradient driven by gap in t", label="", color=:green,
			xlabel="position in t", titlefontsize=12, xguidefontsize=10),
		layout=@layout([a{0.6w} [b;c]]), size=(800, 500)
		)
			
	
end

# ╔═╡ Cell order:
# ╟─1833a9ea-9d7e-400c-85aa-64a37b04ab8a
# ╟─56d7d039-ca74-49ff-ab99-a642e75bf2ef
# ╟─52d9e9ec-c86d-496c-b47c-f5fd7b4a6af3
# ╟─b36ef14e-088c-483f-9dbe-edbbd951830a
# ╟─0a9bccaa-72b2-4394-b561-1bdcb543a11d
# ╟─77546a56-3b11-466e-a5ad-8bef6115ac97
# ╟─419eea6b-0ace-406f-8977-1ecfc2c18747
# ╟─ebc96485-9c73-4307-9799-531a1abdd7d4
# ╟─feeac491-ad6d-48d3-95d2-71c5cbfc09ec
# ╟─e7735020-55f9-4dd2-bf78-67d60451aab3
# ╟─d4a8ec90-0653-47c3-aa86-30ef0ae298fa
# ╟─df424f8f-c454-4c0a-a3b6-cad807df77cf
# ╟─c220f01e-ecd5-4b17-be4c-4ebeb2062277
# ╟─67fe4cd2-f3cc-48be-825a-4fe493e19687
# ╟─0251ca0d-f3a4-4758-a451-3b158b187992
# ╟─39ca6f98-0dac-429b-954b-f04c1e703db0
# ╟─0a73991d-1ae7-4908-94ab-c0cbbf25a9c5
# ╟─748e9149-800d-464f-87ea-6b51438baac0
# ╟─46923e27-adc2-48a8-9ca0-979c3d32a5f3
# ╟─6e45cee4-466d-415b-85c9-c5a93be88d4e
# ╟─1f48c0da-a872-4674-85ed-21de2d4fe6ab
# ╟─f5aa16a6-ad76-48b2-af88-2b0439b56a21
# ╟─29f3cf02-d7b3-470a-9a4d-8a988e41f366
# ╟─3bc02cf8-62c2-4cf2-977b-75650e0c6cf4
# ╟─f0c5d254-8f0f-41eb-b3b5-56c49cfac65c
# ╟─6601ac61-dd7e-474f-8c0b-b0e81cb1b5fa
# ╟─aefbe8d1-528e-46e0-aa45-b3294f8ef54d
# ╟─85abc158-ce84-4eff-9d01-3a75d5483e1b
# ╟─b5252d6b-db49-4e38-9fb7-2c4fedabf2e9
# ╟─4ed12b25-4aca-41aa-93c9-44496c831e5e
# ╟─dd3f27f5-0513-422d-a9b9-1da41c88b38c
# ╟─5ec50f03-9c0c-4dd6-b8dc-cfb17d9b32b0
# ╟─ce115b70-e2dc-460e-b8ca-6c1e9c7ea2b1
