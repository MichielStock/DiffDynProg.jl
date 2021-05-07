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
"MIVLVTCLLLLCSYHTVLSTTNNECIQVNVTQLAGNENLIRDFLFSNFKEEGSVVVGGYYPTEVWYNCSRTARTTAFQYFNNIHAFYFVMEAMENSTGNARGKPLLFHVHGEPVSVIISAYRDDVQQRPLLKHGLVCITKNRHINYEQFTSNQWNSTCTGADRKIPFSVIPTDNGTKIYGLEWNDDFVTAYISGRSYHLNINTNWFNNVTLLYSRSSTATWEYSAAYAYQGVSNFTYYKLNNTNGLKTYELCEDYEHCTGYATNVFAPTSGGYIPDGFSFNNWFLLTNSSTFVSGRFVTNQPLLINCLWPVPSFGVAAQEFCFEGAQFSQCNGVSLNNTVDVIRFNLNFTADVQSGMGATVFSLNTTGGVILEISCYSDTVSESSSYSYGEIPFGITDGPRYCYVLYNGTALKYLGTLPPSVKEIAISKWGHFYINGYNFFSTFPIGCISFNLTTGVSGAFWTIAYTSYTEALVQVENTAIKNVTYCNSHINNIKCSQLTANLNNGFYPVASSEVGFVNKSVVLLPSFFTYTAVNITIDLGMKLSGYGQPIASTLSNITLPMQDNNTDVYCIRSNQFSVYVHSTCKSSLWDNIFNQDCTDVLEATAVIKTGTCPFSFDKLNNYLTFNKFCLSLSPVGANCKFDVAARTRTNEQVVRSLYVIYEEGDNIVGVPSDNSGLHDLSVLHLDSCTDYNIYGRTGVGIIRRTNSTLLSGLYYTSLSGDLLGFKNVSDGVIYSVTPCDVSAQAAVIDGAIVGAMTSINSELLGLTHWTTTPNFYYYSIYNYTSERTRGTAIDSNDVDCEPVITYSNIGVCKNGALVFINVTHSDGDVQPISTGNVTIPTNFTISVQVEYMQVYTTPVSIDCARYVCNGNPRCNKLLTQYVSACQTIEQALAMGARLENMEVDSMLFVSENALKLASVEAFNSTENLDPIYKEWPSIGGSWLGGLKDILPSHNSKRKYGSAIEDLLFDKVVTSGLGTVDEDYKRCTGGYDIADLVCAQYYNGIMVLPGVANADKMTMYTASLAGGITLGALGGGAVAIPFAVAVQARLNYVALQTDVLNKNQQILANAFNQAIGNITQAFGKVNDAIHQTSQGLATVAKALAKVQDVVNTQGQALSHLTVQLQNNFQAISSSISDIYNRLDELSADAQVDRLITGRLTALNAFVSQTLTRQAEVRASRQLAKDKVNECVRSQSQRFGFCGNGTHLFSLANAAPNGMIFFHTVLLPTAYETVTAWSGICASDGDRTFGLVVKDVQLTLFRNLDDKFYLTPRTMYQPRVATSSDFVQIEGCDVLFVNATVIDLPSIIPDYIDINQTVQDILENYRPNWTVPEFTLDIFNATYLNLTGEIDDLEFRSEKLHNTTVELAILIDNINNTLVNLEWLNRIETYVKWPWYVWLLIGLVVVFCIPLLLFCCFSTGCCGCIGCLGSCCHSICSRRQFENYEPIEKVHVH"=>"CoV_S2spike_var2"], default="gradient"))
"""

# ╔═╡ 52d9e9ec-c86d-496c-b47c-f5fd7b4a6af3
md"""
Alternatively, supply your own protein sequences:
$(@bind seq_free1 TextField())
$(@bind seq_free2 TextField())
"""

# ╔═╡ 66607a15-6aab-408a-b30a-4b038706b447
begin
	s = (length(seq_free1) > 0) ? seq_free1 : sequence1;
	t = (length(seq_free2) > 0) ? seq_free2 : sequence2;
end;

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

# ╔═╡ feeac491-ad6d-48d3-95d2-71c5cbfc09ec
md"Show sequence on plots: $(@bind show_seq CheckBox(default=true))"

# ╔═╡ e7735020-55f9-4dd2-bf78-67d60451aab3
begin
	
	md"""
	Zoom range for the right plot if applicable (x-axis & y-axis, respectively):
	
	$(@bind yrange RangeSlider(1:n; show_value=true))
	
	$(@bind xrange RangeSlider(1:m; show_value=true))
	"""
end

# ╔═╡ f530b954-e26f-45af-8826-cb69fee2d8a1
subrange = !(extrema(yrange) == (1, n) && extrema(xrange) == (1, m))

# ╔═╡ 963e255f-de9a-4ac3-a378-abf01396dbcb
subregion = ([minimum(xrange), maximum(xrange), maximum(xrange), minimum(xrange), minimum(xrange)],
			[minimum(yrange), minimum(yrange), maximum(yrange), maximum(yrange), minimum(yrange)])

# ╔═╡ 9de7d808-5a41-428c-8ff3-745a242db7ac
begin
	xticks_sub = 1:length(xrange), split(t, "")[xrange]
	yticks_sub = 1:length(yrange), split(s, "")[yrange]
end;

# ╔═╡ d4a8ec90-0653-47c3-aa86-30ef0ae298fa
md"""
## Paiwise substitution matrix

The matrix $\theta$ quantifies the (dis)similarity between the residues. Currently showing $submat.
"""

# ╔═╡ 4f5709d1-9a78-4fbf-9e68-e7c7c7e85f26
begin
	# set the correct substitution matrix and compute theta
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
	pθ = heatmap(θ; yflip=true, color=:vik, clims=(-6.0, 6.0), axes=false)
	show_seq && yticks!((1:n, split(s, "")))
	show_seq && xticks!((1:m, split(t, "")))
	if subrange
		plot!(pθ, subregion..., color="orange", lw=2, label="")
	end
	pθ
end

# ╔═╡ c220f01e-ecd5-4b17-be4c-4ebeb2062277
if subrange
	pθs = heatmap(θ[yrange, xrange]; yflip=true, color=:vik, clims=(-6.0, 6.0), axes=false)
	show_seq && yticks!(yticks_sub...)
	show_seq && xticks!(xticks_sub...)
	pθs
end

# ╔═╡ 0251ca0d-f3a4-4758-a451-3b158b187992
θ[yrange, xrange]

# ╔═╡ 39ca6f98-0dac-429b-954b-f04c1e703db0
md"""
## Dynamic programming matrix

Below is the dynamic programming matrix $D$. 
"""

# ╔═╡ 6e45cee4-466d-415b-85c9-c5a93be88d4e
if alignment_type != "global alignment"
	md"The change of this alignment score w.r.t. the elements of $D$ is given below."
end

# ╔═╡ f5aa16a6-ad76-48b2-af88-2b0439b56a21
md""" 
## Gradients

The gradient of the dynamic programming matrix is given below:
"""

# ╔═╡ f0c5d254-8f0f-41eb-b3b5-56c49cfac65c
md"""
This gradient can be decomposed in three components:
- a component depending on inserting gaps in $s$;
- a component depending on extending the alignment;
- a component depending on inserting gaps in $t$.
"""

# ╔═╡ e636e4a7-aa5e-4089-adec-ae01038cbe1d
chars_s = s |> unique |> sort!

# ╔═╡ af170735-63eb-4aac-9a63-e0349d55a805
chars_s

# ╔═╡ f9e960b2-ca06-48c5-9e1c-e58d8cf9420e
chars_t = t |> unique |> sort!

# ╔═╡ ce115b70-e2dc-460e-b8ca-6c1e9c7ea2b1
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
	pD = heatmap(D; yflip=true, color=:deep, axes=false)
	show_seq && yticks!((1:n, split(s, "")))
	show_seq && xticks!((1:m, split(t, "")))
	if subrange
		plot!(pD, subregion..., color="orange", lw=2, label="")
	end
	pD
end

# ╔═╡ 748e9149-800d-464f-87ea-6b51438baac0
if subrange
	pDs = heatmap(D[yrange,xrange]; yflip=true, color=:deep, axes=false)
	show_seq && yticks!(yticks_sub...)
	show_seq && xticks!(xticks_sub...)
	pDs
end

# ╔═╡ 46923e27-adc2-48a8-9ca0-979c3d32a5f3
if alignment_type == "global alignment"
	md"For global alignment, the alignment score is the final value in this matrix, here: $v"
else
	md"In local alignment, we take the log-sum-exp of the matrix $D$, obtaining an alignment score of $v"
end

# ╔═╡ abde38be-3988-4765-9f95-607efb779e42
if alignment_type != "global alignment"
	pM = plot(heatmap(M; yflip=true, 
				color=cgrad([:white, :purple]), clims=(0,1), axes=false))
	show_seq && yticks!((1:n, split(s, "")))
	show_seq && xticks!((1:m, split(t, "")))
	if subrange
		plot!(pM, subregion..., color="orange", lw=2, label="")
	end
	pM
end

# ╔═╡ 1f48c0da-a872-4674-85ed-21de2d4fe6ab
if subrange
	pMs = plot(heatmap(M[yrange,xrange]; yflip=true, 
				color=cgrad([:white, :purple]), clims=(0,1), axes=false))
	show_seq && yticks!(yticks_sub...)
	show_seq && xticks!(xticks_sub...)
	pMs
end

# ╔═╡ 29f3cf02-d7b3-470a-9a4d-8a988e41f366
begin
	pE = plot(heatmap(E, yflip=true, 
				color=cgrad([:white, :green]), clims=(0,1), axes=false))
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
				color=cgrad([:white, :green]), clims=(0,1), axes=false))
	show_seq && yticks!(yticks_sub...)
	show_seq && xticks!(xticks_sub...)
	pEs
end

# ╔═╡ 6601ac61-dd7e-474f-8c0b-b0e81cb1b5fa
begin
	pEcomb = plot(
		heatmap(Egs, yflip=true, 
					color=cgrad([:white, :green]), clims=(0,1), axes=false,
			title="gap in s"),
		heatmap(Eθ, yflip=true, 
					color=cgrad([:white, :green]), clims=(0,1), axes=false,
			title="extending alignment"),
		heatmap(Egt, yflip=true, 
					color=cgrad([:white, :green]), clims=(0,1), axes=false,
			title="gap in t"),
		layout=(3, 1), size=(800, 1200))
	show_seq && yticks!((1:n, split(s, "")))
	show_seq && xticks!((1:m, split(t, "")))
	pEcomb
end

# ╔═╡ 84676622-2bd2-4c62-99cb-1dd98f9518f9
begin
	dgs = sum(Egs, dims=2)[:]
	dgt = sum(Egt, dims=1)[:]
end;

# ╔═╡ 85abc158-ce84-4eff-9d01-3a75d5483e1b
begin
	# compute gradient substitution matrix
	dS = zeros(length(chars_s), length(chars_t))
	for (i, sᵢ) in enumerate(s)
		for (j, tⱼ) in enumerate(t)
			dS[findfirst(isequal(sᵢ), chars_s), findfirst(isequal(tⱼ), chars_t)] += Eθ[i,j]
		end
	end
	pS = plot(
		heatmap(dS, yflip=true, yticks=(1:length(chars_s), chars_s),
				xticks=(1:length(chars_t), chars_t),
					color=cgrad([:white, :green]), axes=false,
			title="gradient substitution matrix"),
		
		plot(dgs, title="gradient gs", label="", color=:green,
			xlabel="position in s"),
		plot(dgt, title="gradient gt", label="", color=:green,
			xlabel="position in t"),
		layout=@layout([a{0.6w} [b;c]]), size=(800, 500)
		)
			
	
end

# ╔═╡ Cell order:
# ╠═1833a9ea-9d7e-400c-85aa-64a37b04ab8a
# ╠═56d7d039-ca74-49ff-ab99-a642e75bf2ef
# ╠═52d9e9ec-c86d-496c-b47c-f5fd7b4a6af3
# ╠═66607a15-6aab-408a-b30a-4b038706b447
# ╠═3bb2cfd3-e30e-4276-b70d-070b3de0a033
# ╠═b36ef14e-088c-483f-9dbe-edbbd951830a
# ╠═0a9bccaa-72b2-4394-b561-1bdcb543a11d
# ╠═77546a56-3b11-466e-a5ad-8bef6115ac97
# ╠═3cab7e29-4750-4ca9-b3c1-910e16d00eaf
# ╠═419eea6b-0ace-406f-8977-1ecfc2c18747
# ╠═ebc96485-9c73-4307-9799-531a1abdd7d4
# ╠═feeac491-ad6d-48d3-95d2-71c5cbfc09ec
# ╠═e7735020-55f9-4dd2-bf78-67d60451aab3
# ╠═f530b954-e26f-45af-8826-cb69fee2d8a1
# ╠═963e255f-de9a-4ac3-a378-abf01396dbcb
# ╠═9de7d808-5a41-428c-8ff3-745a242db7ac
# ╠═d4a8ec90-0653-47c3-aa86-30ef0ae298fa
# ╠═4f5709d1-9a78-4fbf-9e68-e7c7c7e85f26
# ╠═df424f8f-c454-4c0a-a3b6-cad807df77cf
# ╠═c220f01e-ecd5-4b17-be4c-4ebeb2062277
# ╠═0251ca0d-f3a4-4758-a451-3b158b187992
# ╟─39ca6f98-0dac-429b-954b-f04c1e703db0
# ╠═0a73991d-1ae7-4908-94ab-c0cbbf25a9c5
# ╠═748e9149-800d-464f-87ea-6b51438baac0
# ╠═46923e27-adc2-48a8-9ca0-979c3d32a5f3
# ╠═6e45cee4-466d-415b-85c9-c5a93be88d4e
# ╠═abde38be-3988-4765-9f95-607efb779e42
# ╠═1f48c0da-a872-4674-85ed-21de2d4fe6ab
# ╠═f5aa16a6-ad76-48b2-af88-2b0439b56a21
# ╠═29f3cf02-d7b3-470a-9a4d-8a988e41f366
# ╠═3bc02cf8-62c2-4cf2-977b-75650e0c6cf4
# ╠═f0c5d254-8f0f-41eb-b3b5-56c49cfac65c
# ╠═6601ac61-dd7e-474f-8c0b-b0e81cb1b5fa
# ╠═84676622-2bd2-4c62-99cb-1dd98f9518f9
# ╠═af170735-63eb-4aac-9a63-e0349d55a805
# ╠═e636e4a7-aa5e-4089-adec-ae01038cbe1d
# ╠═f9e960b2-ca06-48c5-9e1c-e58d8cf9420e
# ╠═85abc158-ce84-4eff-9d01-3a75d5483e1b
# ╠═ce115b70-e2dc-460e-b8ca-6c1e9c7ea2b1
