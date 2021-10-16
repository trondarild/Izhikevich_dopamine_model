### A Pluto.jl notebook ###
# v0.16.1

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

# ╔═╡ 22821822-1b9f-11ec-33d8-e589cf58d1ac
using Plots

# ╔═╡ bf662a28-373b-4e4b-b8e8-b19be85b671b
using PlutoUI

# ╔═╡ 8092bfb1-f042-43ed-86fc-a9815ccfb555
using Accessors

# ╔═╡ 74589c68-0428-48e2-82bd-6bde8ff3def2
md"""
# Izhikevich model with dopamine support

These equations model GABAergic *medium spiny interneurons* in the direct- and indirect pathways of the basal ganglia.

### Original Iz. equations:

$C\dot{v} = k(v-v_r)(v-v_t)-u+I$

$\dot{u} = a[b(v-v_r)-u]$

$\text{if } v \ge v_{peak} \text{ then } v \leftarrow c, u \leftarrow u + d$

### Due to Humphries et al. 2009 :

$v_r \leftarrow v_r(1 + K\phi_1)$


$d \leftarrow d(1-L\phi_1)$

$k \leftarrow k(1-\alpha\phi_2)$

$I = I_{ampa} + B(v)I_{nmda} + I_{gaba}$

$I_z = \bar{g_z}h_z(E_z-v)$ 

$\dot{h_z} = \frac{-h_z}{\tau_z}$

$h_z(t) \leftarrow \frac{h_z(t)+S_z(t)}{\tau_z}$

$z \in [nmda, ampa, gaba]$

$B(v) = \frac{1}{1 + \frac{[Mg^{2+}]_0}{3.57}exp(-0.062v)}$

$I^{D1}_{nmda} = I_{nmda}(1+\beta_1\phi_1)$

$I^{D2}_{ampa} = I_{ampa}(1-\beta_2\phi_2)$


* where $\phi_1$ is proportion of activated D1 receptors; 
* and $\alpha$ parameterises size of appropriate reduction in $k$ to get the required f-I curve
* and $K$ is K+ current constant, $L$ is Ca2+ current constant (L-current)
* and $\bar{g_z}$ is maximum conductance
* and $E_z$ is reversal potential (?)
* and $\tau_z$ is synaptic time constant
* and $S_z(t)$ 
"""

# ╔═╡ 553cb9c7-2cf7-4fac-8cce-e16231e37ae6


# ╔═╡ ece420a5-0dee-4f2b-b154-4c5a8a9c7a76
function f_vr(v_r, d1, K)
	return v_r * (1 + K * d1)
end

# ╔═╡ 990c7dff-88da-46a1-a899-ae8ee0cd1366
function f_d(d, d1, L)
	return d * (1-L*d1)
end

# ╔═╡ 43a9a3e2-81b2-40e5-8622-b27bf88e5675
function f_k(k, d2, alpha)
	return k*(1-alpha*d2)
end

# ╔═╡ 9be7e945-6d91-4fd7-8bbe-2250e550bab0


# ╔═╡ f073c8ca-9918-4a41-adde-3becab1302ae
# update nmda current based on D1 activation
function f_ID1nmda(Inmda, d1, beta1)
	return Inmda * (1 + beta1 * d1)
end

# ╔═╡ dde99a00-db75-47ce-a835-1b18c74e2999
# update ampa current based on D2 activation
function f_ID2ampa(Iampa, d2, beta2)
	return Iampa * (1 - beta2 * d2)
end

# ╔═╡ 85802bb3-0684-40bc-a46a-1347ac3cbabb
#
# h ?
#

# ╔═╡ 9cae35bf-f3e5-4975-92c8-dcd0fa7dc6b7
function f_h(h_prev, sum_spikes, tau)
	return (h_prev + sum_spikes)/tau
end

# ╔═╡ 5ae88d5b-ee7b-4a40-8440-0eb59e610b27
#
# currents
#

# ╔═╡ 80c336df-b6c0-4791-a961-e96f437d961a
md"""
## Part-equation test
"""

# ╔═╡ 6467cfd6-d748-410c-9183-7f455482e317
data = Float16[]

# ╔═╡ c110f3f4-1f4e-4ea5-a774-390da7305e89
vtst = -50

# ╔═╡ 66d3d066-4e8f-429d-8d10-b83e13185623


# ╔═╡ 590fd0e9-7586-4e2b-9017-857a8bde87e8
md"""
# Constants
"""

# ╔═╡ 854e80d2-5f71-4c24-aa32-b1240feee027
md"""
Tuning results from Humphries 2009:
1. c C = 15.2
2. c $v_t$ = -29.7
3. c d = 91
4. c K = 0.0289
5. c L = 0.331
6. c $\alpha$ = 0.032
7. c $g_{ampa}$ = 6.9
7. c $g_{nmda}$ = 6.9 * 0.5
7. c $g_{gaba}$ = 6.9 / 1.4
8. c $\beta_1 = 6.3$
9. c $\beta_2 = 0.215$

"""

# ╔═╡ f744d017-5f8a-486c-a204-a5f9ba1e4e90
"""
# values from Iz model
begin
	const c_a = 0.1
	const c_b = 0.2
	const c_c = -65.0
	# const c_d = 8.0
	
	const c_thrs = 30.0
	const i_v = -65.0
	const i_u = -14.0
	const dt = 0.1
end
"""

# ╔═╡ 77040227-ebb1-4688-a617-f31fdb12c80a
# values from Humphries 2009, medium spiny neuron
begin
	const c_msn_a = 0.01
	const c_msn_b = -20
	const c_msn_c = -55.0
	const c_msn_thrs = 30.0
	
	const i_msn_v = -65.0
	const i_msn_u = -14.0
	
	const c_msn_k = 1
	const c_msn_vr = -80
	const c_msn_vpeak = 40
	const c_msn_Eampa = 0
	const c_msn_Enmda = 0
	const c_msn_Egaba = -60
	const c_msn_tau_ampa = 6
	const c_msn_tau_nmda = 160
	const c_msn_tau_gaba = 4
	const c_msn_Mg = 1
	
	const c_dt = 0.1
end

# ╔═╡ 73c30789-adba-45cd-b032-0df3bcf2279d
# Tuning results from Humphries 2009, medium spiny neuron
begin
    const c_msn_C = 15.2
	const c_msn_vt = -29.7
	const c_msn_d = 91
	const c_msn_K = 0.0289
	const c_msn_L = 0.331
	const c_msn_alpha = 0.032
	const c_msn_g_ampa = 6.9
	const c_msn_g_nmda = 0.5 * c_msn_g_ampa
	const c_msn_g_gaba = c_msn_g_ampa / 1.4
	const c_msn_beta_1 = 6.3
	const c_msn_beta_2 = 0.215
end

# ╔═╡ a78a41c7-c7d7-4fa2-954a-fc632008ee49
struct Iz_NeuronParameter
	c_a::Float16
	c_b::Float16
	c_c::Float16
	c_d::Float16
	c_thrs::Float16
	i_v::Float16
	i_u::Float16
end

# ╔═╡ b69a9769-0178-440b-809b-ca49918fdd9a
struct Hum_NeuronParameter
	c_iz::Iz_NeuronParameter
	
    c_k::Float16
	c_vr::Float16
	c_vpeak::Float16
	c_Eampa::Float16
	c_Enmda::Float16
	c_Egaba::Float16
	c_tau_ampa::Float16
	c_tau_nmda::Float16
	c_tau_gaba::Float16
	c_Mg::Float16
	
	c_C::Float16
	c_vt::Float16
	c_d::Float16
	c_K::Float16
	c_L::Float16
	c_alpha::Float16
	c_g_ampa::Float16
	c_g_nmda::Float16
	c_g_gaba::Float16
	c_beta_1::Float16
	c_beta_2::Float16
end

# ╔═╡ bc37bffc-b3f0-4014-982e-83d86c2c0a58
md"""
## Container for neural state
"""

# ╔═╡ cbd1b4a3-3246-4fa7-a248-f1ddf55395c7
struct NeuronState 
	v::Float16 # voltage
	u::Float16 # feedback current
	a::Float16 # 
	b::Float16 #
	c::Float16 #
	d::Float16 #
	d1::Float16 # proportion of D1 activation
	d2::Float16 # proportion of D2 activation
	
	h_ampa::Float16 # ampa current via spike count
	h_nmda::Float16
	h_gaba::Float16
end

# ╔═╡ dac654bd-9f19-42ca-9e9b-b444b2fbaa2d
md"""
## Neuronal computations
"""

# ╔═╡ df6382d8-3550-4937-95fb-50fc8aef341b
md"""
### Humphries dopamine modifications
"""

# ╔═╡ 46f92f9e-6ed3-43cc-881d-ed6c311e0f90
function f_B(v)
	return 1/(1+1/3.57 * exp(-0.062*v))
end

# ╔═╡ 5b935650-40bf-4ecd-ad6c-04b04acf5e51
# calculate total input current
function f_Itot(Iampa, Inmda, Igaba, v)
	return Iampa + f_B(v)*Inmda + Igaba
end

# ╔═╡ 2ca7f1cd-6ecf-4b61-8c90-83de6f7cec91
begin
	#data = Float16[]
	empty!(data)
	inp = range(-50,stop=50,length=100)
	
	for i in inp 
		tmp = f_B(i)
		append!(data, tmp)
	end
	plot(inp, data, title="B(v)")
end 

# ╔═╡ dfed20ee-7190-4dcd-a673-91924f356c3e
md"""
### Iz original
"""

# ╔═╡ ef55479b-b347-4fc7-a085-b2d726459907
function f_v(st, i, dt)
	v = st.v
	u = st.u
	return (0.04 * v * v + 5 * v + 140 - u + i) * dt
end

# ╔═╡ 8522a885-5e2f-448d-a336-64c1955a9cf3


# ╔═╡ ee249c45-d6cb-4ce8-9a64-87b1fcdf188f
function step(st, I, dt, c_thrs)
	ret_v = 0
	ret_u = 0
	if st.v >= c_thrs
		ret_v = st.c
		ret_u = st.u + st.d
	else
		dv = f_v(st, I, dt)
		du = f_u(st, dt)
		ret_v = st.v + dv
		ret_u = st.u + du
	end
	
	return NeuronState(ret_v, ret_u, st.a, st.b, st.c, st.d, st.d1, st.d2)
	
end


# ╔═╡ 868ec52d-856d-4324-b854-bbaf244e0c57
# access lenses for parameters
begin
	ln_a = @optic _.c_iz.c_a
	ln_b = @optic _.c_iz.c_b
	ln_thr = @optic _.c_iz.c_thrs
	ln_c = @optic _.c_iz.c_c
	ln_d = @optic _.c_iz.c_d
	ln_C = @optic _.c_C
	ln_k = @optic _.c_k
	ln_vr = @optic _.c_vr
	ln_vt = @optic _.c_vt
	ln_beta1 = @optic _.c_beta_1
	ln_beta2 = @optic _.c_beta_2
	ln_L = @optic _.c_L
	
	# assembling currents from conductance, spike inputs
	ln_Eampa = @optic _.c_Eampa
	ln_Enmda = @optic _.c_Enmda
	ln_Egaba = @optic _.c_Egaba
	ln_tau_ampa = @optic _.c_tau_ampa
	ln_tau_nmda = @optic _.c_tau_nmda
	ln_tau_gaba = @optic _.c_tau_gaba
	ln_g_ampa = @optic _.c_g_ampa
	ln_g_nmda = @optic _.c_g_nmda
	ln_g_gaba = @optic _.c_g_gaba
end

# ╔═╡ cd20e2ae-25e3-473f-9c9b-b005229c3c95
function f_h_ampa(h_ampa_prev, sum_ampa, params)
	return f_h(h_ampa_prev, sum_ampa, ln_tau_ampa(params))
end

# ╔═╡ c4253b6d-b02b-4d8c-873e-62f19551c007
function f_h_nmda(h_nmda_prev, sum_nmda, params)
	return f_h(h_nmda_prev, sum_nmda, ln_tau_nmda(params))
end

# ╔═╡ 16793c73-f437-402f-b326-d75a5d58d2b4
function f_h_gaba(h_gaba_prev, sum_gaba, params)
	return f_h(h_gaba_prev, sum_gaba, ln_tau_gaba(params))
end

# ╔═╡ 87a28783-cb8a-4c18-82ee-1aa7ea6cca2c
function f_I_ampa(h_ampa, v, params)
	return ln_g_ampa(params) * h_ampa * (ln_Eampa(params) - v)
end
	

# ╔═╡ 5f19bd11-fdb5-4fb9-98cb-7bdb63d3369f
function f_I_nmda(h_nmda, v, params)
	return ln_g_nmda(params) * h_nmda * (ln_Enmda(params) - v)
end

# ╔═╡ 2fba3958-2c5c-4539-a527-f44d0766e630
function f_I_gaba(h_gaba, v, params)
	return ln_g_gaba(params) * h_gaba * (ln_Egaba(params) - v)
end

# ╔═╡ 3996ad3b-5fcb-481c-8945-a38a68210088
function f_v_dopa(st, i, dt, params)
	v = st.v
	u = st.u
	C = ln_C(params)
	k = ln_k(params)
	vr = ln_vr(params)
	vt = ln_vt(params)
	return (k*(v-vr)*(v-vt) - u + i)*dt/C
end
	

# ╔═╡ 24f338ea-f4cb-45e1-b83b-c902f1bb64f9
function f_u_dopa(st, dt, params)
	v = st.v
	u = st.u
	a = ln_a(params)
	b = ln_b(params)
	vr = ln_vr(params)
	return a * (b * (v-vr) - u) * dt
end

# ╔═╡ 8ddc2629-b4a3-4aa1-a72b-002365ce703d

function step_dopa(st, Iampa, Inmda, Igaba, params)
	ret_v = 0
	ret_u = 0
	ret_d = st.d
	ret_hampa = st.h_ampa
	ret_hnmda = st.h_nmda
	ret_hgaba = st.h_gaba
	if st.v >= ln_thr(params)
		ret_v = ln_c(params)
		ret_d = f_d(ret_d, st.d1, ln_L(params))  # update also d
		# println("ret_d ", ret_d)
		ret_u = st.u + ret_d
	else
		# TODO update currents with h-functions
		hampa = f_h_ampa(st.h_ampa, Iampa, params)
		hnmda = f_h_nmda(st.h_nmda, Inmda, params)
		# println("hnmda ", hnmda)
		hgaba = f_h_gaba(st.h_gaba, Igaba, params)
		
		ampa_pre = f_I_ampa(hampa, st.v, params)
		# println("ampapre ", ampa_pre)
		nmda_pre = f_I_nmda(hnmda, st.v, params)
		# println("nmdapre ", nmda_pre)
		gaba_pre = f_I_gaba(hgaba, st.v, params)
		
		# update currents dept on dopa
		i_ampa = f_ID2ampa(ampa_pre, st.d2, ln_beta2(params))
		# println("i_ampa ", i_ampa)
		i_nmda = f_ID1nmda(nmda_pre, st.d1, ln_beta1(params))	
		# println("i_nmda ", i_nmda)
		Itot = f_Itot(i_ampa, i_nmda, gaba_pre, st.v)
		
		# update h current values
		dhampa = (-hampa/ln_tau_ampa(params))*c_dt
		dhnmda = (-hnmda/ln_tau_nmda(params))*c_dt
		dhgaba = (-hgaba/ln_tau_gaba(params))*c_dt
		dv = f_v_dopa(st, Itot, c_dt, params)
		du = f_u_dopa(st, c_dt, params)
		
		ret_hampa = hampa + dhampa
		ret_hnmda = hnmda + dhnmda
		ret_hgaba = hgaba + dhgaba
		ret_v = st.v + dv
		ret_u = st.u + du
		# println()
	end
	# TODO update h_ values
	return NeuronState(ret_v, ret_u, st.a, st.b, st.c, ret_d, st.d1, st.d2, ret_hampa, ret_hnmda, ret_hgaba)
	
end	

# ╔═╡ 9a0d0774-5134-4643-a4a6-2c7fb1842907
function injcur(a, t, mn, mx)
	return a * (mn < t < mx)
end

# ╔═╡ 6458986b-4b15-4998-aa96-36384379798e
md"""
# Simulation
"""

# ╔═╡ 00aae887-d1ea-4edf-b5e7-0c755b78d972
begin
	#ampl = 100.0 # amplify input current
	i_v=-65 # initial voltage
	i_u=-14 # initial u
end

# ╔═╡ 18fa260b-4789-487e-b279-b0b8c38753db
begin
	vlt = Float16[]
	ampa_data = Float16[]
	nmda_data = Float16[]
	gaba_data = Float16[]
end

# ╔═╡ 51038af2-38bb-4aef-b8d9-a6a3861d056f
injcurdata = Float16[]

# ╔═╡ cae7f940-b184-4a4c-8ee0-9f3a77f94271
# TODO: set up state and msn params; test

# ╔═╡ 21449dfe-196b-42a7-98f6-ccb747521daf
msn_params = Hum_NeuronParameter(
	Iz_NeuronParameter(
		c_msn_a,
		c_msn_b,
		c_msn_c,
		c_msn_d,
		c_msn_thrs,
		i_msn_v,
		i_msn_u),
	c_msn_k,
	c_msn_vr,
	c_msn_vpeak,
	c_msn_Eampa ,
	c_msn_Enmda ,
	c_msn_Egaba,
	c_msn_tau_ampa ,
	c_msn_tau_nmda,
	c_msn_tau_gaba ,
	c_msn_Mg,
	c_msn_C,
	c_msn_vt,
	c_msn_d,
	c_msn_K,
	c_msn_L,
	c_msn_alpha,
	c_msn_g_ampa,
	c_msn_g_nmda,
	c_msn_g_gaba,
	c_msn_beta_1,
	c_msn_beta_2)

# ╔═╡ fd7d426a-a0e4-4ed9-a916-f09209e3042b
d1_ui = @bind d1_val Slider(0:100)

# ╔═╡ dc193a59-7f70-48b2-8a60-0fb872283625
d1_val

# ╔═╡ e696e625-943f-431f-9776-339420f3db91
d2_ui = @bind d2_val Slider(0:100)

# ╔═╡ 532843c0-4dd2-49ab-9fe6-a5865d1d42a0
ampl_ui = @bind ampl Slider(-30:1530)

# ╔═╡ 477d77f2-8626-439a-8087-569e88f0f5ae
ampl

# ╔═╡ fc02eb71-f9fc-4d39-83c4-defa3f8fe613
ampa_ui = @bind ampa_frac Slider(0:100)

# ╔═╡ 46ce7229-e3fa-469e-820d-157cf9ef9d4a
ampa_frac

# ╔═╡ 5188c442-4205-46cc-8396-ef664b11f9b8
nmda_ui = @bind nmda_frac Slider(0:100)

# ╔═╡ 02f8b97a-0f6f-4daa-8fd6-e94f2dcb8875
nmda_frac

# ╔═╡ 759ee2c9-d2cc-4f21-a140-a1e0308892d8
gaba_ui = @bind gaba_frac Slider(0:100)

# ╔═╡ 7931cec0-610e-4836-b2b8-5fcc20d31b34
gaba_frac

# ╔═╡ 7013902a-b421-4abb-a1f2-df3fe58a62a1
state = NeuronState(i_v, i_u, 0.0, 0.0, 0.0, ln_d(msn_params), d1_val/100.0, d2_val/100.0, 10.9, 0.9, 0.9)

# ╔═╡ f7c89ffc-8230-4853-bbb3-d4287682e32a
f_v_dopa(state, 30, c_dt, msn_params)

# ╔═╡ 1901353b-62bd-4328-a6f3-13460c7a770d
begin
	ampl_ampa = ampl * ampa_frac/100.0
	ampl_nmda = ampl * nmda_frac/100.0
	ampl_gaba = -ampl * gaba_frac/100.0
end

# ╔═╡ c95d2ef5-eee1-4286-9702-e40c3282b5ae
begin
	empty!(vlt)
	empty!(ampa_data)
	empty!(nmda_data)
	empty!(gaba_data)
	nst = state
	for i in 2:500
		# push!(steps, step(steps[i-1], injcur(i), dt))
		i_ampa = injcur(ampl_ampa, i, 100, 130)
		i_nmda = injcur(ampl_nmda, i, 0, 200)
		i_gaba = injcur(ampl_gaba, i, 120, 130)
		push!(ampa_data, i_ampa)
		push!(nmda_data, i_nmda)
		push!(gaba_data, i_gaba)
		
		tmp = step_dopa(
			nst, 
			i_ampa,
			i_nmda,
			i_gaba,
			msn_params)
		push!(vlt, tmp.v)
		nst = tmp
		#push!(injcurdata, injcur(ampl, i))
		# push!(plt, state.v)
		# push!(vlt, steps[i].v)
		# println(state.v)
	end 
end

# ╔═╡ be927f12-e945-4aed-965a-62c0691ad23f
md"""
# Plots
"""

# ╔═╡ 007a14b3-583f-4df6-8ef4-eab29b882641
begin
	plot(vlt, title="Voltage", label="neuron")
	plot!(ampa_data, label="AMPA")
	plot!(nmda_data, label="NMDA")
	plot!(gaba_data, label="GABA")
	#plot!(injcurdata, title="Injection")
	#display(p1)
	#display(p2)
end

# ╔═╡ 215cfb07-96bb-4267-938a-efce8b90153a
md"""
# Test
"""

# ╔═╡ 8b3614b1-17d5-462e-a3a0-eb14df11dccc
st = NeuronState(i_v, i_u, c_a, c_b, c_c, c_d, 0.8, 0.8)

# ╔═╡ 7b4b5f56-3316-4ae5-b53d-4d588bda7975
f_v(st, 30, 1)

# ╔═╡ 7886a63b-928a-42a2-a4f5-bfd920eeb2ce
f_u(st, 1)

# ╔═╡ e9282a6f-3239-4e9c-b11e-e6d68aefde2d
st_b = step(st, 30, 1)

# ╔═╡ 2497e93a-3a45-46ac-a3b4-076047bb8c80
md"""
## References

* Humphries et al 2009 - Capturing dopaminergic modulation and bimodal membrane behaviour of striatal medium spiny neurons in accurate, reduced models

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Accessors = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Accessors = "~0.1.4"
Plots = "~1.22.1"
PlutoUI = "~0.7.10"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Accessors]]
deps = ["Compat", "CompositionsBase", "ConstructionBase", "Future", "LinearAlgebra", "MacroTools", "Requires", "Test"]
git-tree-sha1 = "ba270110280297b36d566cb19c948e6c724432c0"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
git-tree-sha1 = "bdf73eec6a88885256f282d48eafcad25d7de494"
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[Artifacts]]
deps = ["Pkg"]
git-tree-sha1 = "c30985d8821e0cd73870b17b0ed0ce6dc44cb744"
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.3.0"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c3598e525718abcc440f69cc6d5f60dda0a1b61e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.6+5"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "e2f47f6d8337369411569fd45ae5753ca10394c6"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.0+6"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "9995eb3977fbf67b86d0a0a0508e83017ded03f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.14.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "1a90210acd935f222ea19657f143004d2c2a1117"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.38.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8e695f735fca77e9708e795eda62afdb869cbb70"
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.3.4+0"

[[CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
git-tree-sha1 = "135bf1896be424235eadb17474b2a78331567f08"
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.5.1"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "92d8f9f208637e8d2d28c664051a00569c01493d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.1.5+1"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1402e52fcda25064f51c77a9655ce8680b76acf0"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.7+6"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "LibVPX_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3cc57ad0a213808473eafef4845a74766242e05f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.3.1+4"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "35895cf184ceaab11fd778b4590144034a167a2f"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.1+14"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "cbd58c9deb1d304f5a245a0b7eb841a2560cfec6"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.1+5"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0d20aed5b14dd4c9a2453c1b601d08e1149679cc"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.5+6"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "a199aefead29c3c2638c3571a9993b564109d45a"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.4+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "c2178cfbc0a5a552e16d097fae508f2024de61a3"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.59.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d59e8320c2747553788e4fc42231489cc602fa50"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.58.1+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "8c14294a079216000a0bdca5ec5a447f073ddc9d"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.20.1+7"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "04690cc5008b38ecbdfede949220bc7d9ba26397"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.59.0+4"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "60ed5f1643927479f845b0135bb369b031b541fa"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.14"

[[HypertextLiteral]]
git-tree-sha1 = "72053798e1be56026b81d4e2682dbe58922e5ec9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.0"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9aff0587d9603ea0de2c6f6300d9f9492bbefbd3"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.0.1+3"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "df381151e871f41ee86cee4f5f6fd598b8a68826"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.0+3"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f128cd6cd05ffd6d3df0523ed99b90ff6f9b349a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.0+3"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
git-tree-sha1 = "cdbe7465ab7b52358804713a53c7fe1dac3f8a3f"
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[LibCURL_jll]]
deps = ["LibSSH2_jll", "Libdl", "MbedTLS_jll", "Pkg", "Zlib_jll", "nghttp2_jll"]
git-tree-sha1 = "897d962c20031e6012bba7b3dcb7a667170dad17"
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.70.0+2"

[[LibGit2]]
deps = ["Printf"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Libdl", "MbedTLS_jll", "Pkg"]
git-tree-sha1 = "717705533148132e5466f2924b9a3657b16158e8"
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.9.0+3"

[[LibVPX_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "85fcc80c3052be96619affa2fe2e6d2da3908e11"
uuid = "dd192d2f-8180-539f-9fb4-cc70b1dcf69a"
version = "1.9.0+1"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "a2cd088a88c0d37eef7d209fd3d8712febce0d90"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.1+4"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "b391a18ab1170a2e568f9fb8d83bc7c780cb9999"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.5+4"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ec7f2e8ad5c9fa99fc773376cdbc86d9a5a23cb7"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.36.0+3"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cba7b560fcc00f8cd770fa85a498cbc1d63ff618"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.0+8"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51ad0c01c94c1ce48d5cad629425035ad030bfd5"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.34.0+3"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "291dd857901f94d683973cdf679984cdf73b56d0"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.1.0+2"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f879ae9edbaa2c74c922e8b85bb83cc84ea1450b"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.34.0+7"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "5a5bc6bf062f0f95e62d0fe0a2d99699fed82dd9"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.8"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0eef589dd1c26a3ac9d753fe1a8bcad63f956fa6"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.16.8+1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f1662575f7bf53c73c2bbc763bace4b024de822c"
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2021.1.19+0"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
git-tree-sha1 = "ed3157f48a05543cce9b241e1f2815f7e843d96e"
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "a42c0f138b9ebe8b58eba2271c5053773bde52d0"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.4+2"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "71bbbc616a1d710879f5a1021bcba65ffba6ce58"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.1+6"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f9d57f4126c39565e05a2b0264df99f497fc6f37"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.1+3"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1b556ad51dceefdbf30e86ffa8f528b73c7df2bb"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.42.0+4"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6a20a83c1ae86416f0a5de605eaea08a552844a3"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.0+0"

[[Pkg]]
deps = ["Dates", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "UUIDs"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "2537ed3c0ed5e03896927187f5f2ee6a4ab342db"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.14"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "4c2637482176b1c2fb99af4d83cb2ff0328fc33c"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.22.1"

[[PlutoUI]]
deps = ["Base64", "Dates", "HypertextLiteral", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "26b4d16873562469a0a1e6ae41d90dec9e51286d"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.10"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "16626cfabbf7206d60d84f2bf4725af7b37d4a77"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.2+0"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3240808c6d463ac46f1c1cd7638375cd22abbccb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.12"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8cbbc098554648c84f79a463c9ff0fd277144b6c"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.10"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
git-tree-sha1 = "44aaac2d2aec4a850302f9aa69127c74f0c3787e"
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "1162ce4a6c4b7e31e0e6b14486a6986951c73be9"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.2"

[[Test]]
deps = ["Distributed", "InteractiveUtils", "Logging", "Random"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "dc643a9b774da1c2781413fd7b6dcd2c56bb8056"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.17.0+4"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "be0db24f70aae7e2b89f2f3092e93b8606d659a6"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.10+3"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "2b3eac39df218762d2d005702d601cd44c997497"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.33+4"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "320228915c8debb12cb434c59057290f0834dbf6"
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.11+18"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "2c1332c54931e83f8f94d310fa447fd743e8d600"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.4.8+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "acc685bcf777b2202a904cdcb49ad34c2fa1880c"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.14.0+4"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7a5780a0d9c6864184b3a2eeeb833a0c871f00ab"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "0.1.6+4"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "6abbc424248097d69c0c87ba50fcb0753f93e0ee"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.37+6"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "fa14ac25af7a4b8a7f61b287a124df7aab601bcd"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.6+6"

[[nghttp2_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "8e2c44ab4d49ad9518f359ed8b62f83ba8beede4"
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.40.0+2"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d713c1ce4deac133e3334ee12f4adff07f81778f"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2020.7.14+2"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "487da2f8f2f0c8ee0e83f39d13037d6bbf0a45ab"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.0.0+3"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═22821822-1b9f-11ec-33d8-e589cf58d1ac
# ╠═bf662a28-373b-4e4b-b8e8-b19be85b671b
# ╠═8092bfb1-f042-43ed-86fc-a9815ccfb555
# ╟─74589c68-0428-48e2-82bd-6bde8ff3def2
# ╠═553cb9c7-2cf7-4fac-8cce-e16231e37ae6
# ╠═ece420a5-0dee-4f2b-b154-4c5a8a9c7a76
# ╠═990c7dff-88da-46a1-a899-ae8ee0cd1366
# ╠═43a9a3e2-81b2-40e5-8622-b27bf88e5675
# ╠═9be7e945-6d91-4fd7-8bbe-2250e550bab0
# ╠═f073c8ca-9918-4a41-adde-3becab1302ae
# ╠═dde99a00-db75-47ce-a835-1b18c74e2999
# ╠═5b935650-40bf-4ecd-ad6c-04b04acf5e51
# ╠═85802bb3-0684-40bc-a46a-1347ac3cbabb
# ╠═9cae35bf-f3e5-4975-92c8-dcd0fa7dc6b7
# ╠═cd20e2ae-25e3-473f-9c9b-b005229c3c95
# ╠═c4253b6d-b02b-4d8c-873e-62f19551c007
# ╠═16793c73-f437-402f-b326-d75a5d58d2b4
# ╠═5ae88d5b-ee7b-4a40-8440-0eb59e610b27
# ╠═87a28783-cb8a-4c18-82ee-1aa7ea6cca2c
# ╠═5f19bd11-fdb5-4fb9-98cb-7bdb63d3369f
# ╠═2fba3958-2c5c-4539-a527-f44d0766e630
# ╟─80c336df-b6c0-4791-a961-e96f437d961a
# ╠═6467cfd6-d748-410c-9183-7f455482e317
# ╠═c110f3f4-1f4e-4ea5-a774-390da7305e89
# ╠═66d3d066-4e8f-429d-8d10-b83e13185623
# ╟─590fd0e9-7586-4e2b-9017-857a8bde87e8
# ╟─854e80d2-5f71-4c24-aa32-b1240feee027
# ╟─f744d017-5f8a-486c-a204-a5f9ba1e4e90
# ╠═77040227-ebb1-4688-a617-f31fdb12c80a
# ╠═73c30789-adba-45cd-b032-0df3bcf2279d
# ╠═a78a41c7-c7d7-4fa2-954a-fc632008ee49
# ╠═b69a9769-0178-440b-809b-ca49918fdd9a
# ╟─bc37bffc-b3f0-4014-982e-83d86c2c0a58
# ╠═cbd1b4a3-3246-4fa7-a248-f1ddf55395c7
# ╟─dac654bd-9f19-42ca-9e9b-b444b2fbaa2d
# ╟─df6382d8-3550-4937-95fb-50fc8aef341b
# ╠═46f92f9e-6ed3-43cc-881d-ed6c311e0f90
# ╠═2ca7f1cd-6ecf-4b61-8c90-83de6f7cec91
# ╟─dfed20ee-7190-4dcd-a673-91924f356c3e
# ╠═ef55479b-b347-4fc7-a085-b2d726459907
# ╠═3996ad3b-5fcb-481c-8945-a38a68210088
# ╠═f7c89ffc-8230-4853-bbb3-d4287682e32a
# ╠═24f338ea-f4cb-45e1-b83b-c902f1bb64f9
# ╠═8522a885-5e2f-448d-a336-64c1955a9cf3
# ╠═ee249c45-d6cb-4ce8-9a64-87b1fcdf188f
# ╠═868ec52d-856d-4324-b854-bbaf244e0c57
# ╠═8ddc2629-b4a3-4aa1-a72b-002365ce703d
# ╠═9a0d0774-5134-4643-a4a6-2c7fb1842907
# ╟─6458986b-4b15-4998-aa96-36384379798e
# ╠═00aae887-d1ea-4edf-b5e7-0c755b78d972
# ╠═18fa260b-4789-487e-b279-b0b8c38753db
# ╠═51038af2-38bb-4aef-b8d9-a6a3861d056f
# ╠═cae7f940-b184-4a4c-8ee0-9f3a77f94271
# ╠═21449dfe-196b-42a7-98f6-ccb747521daf
# ╠═fd7d426a-a0e4-4ed9-a916-f09209e3042b
# ╠═dc193a59-7f70-48b2-8a60-0fb872283625
# ╠═e696e625-943f-431f-9776-339420f3db91
# ╠═532843c0-4dd2-49ab-9fe6-a5865d1d42a0
# ╠═477d77f2-8626-439a-8087-569e88f0f5ae
# ╠═fc02eb71-f9fc-4d39-83c4-defa3f8fe613
# ╠═46ce7229-e3fa-469e-820d-157cf9ef9d4a
# ╠═5188c442-4205-46cc-8396-ef664b11f9b8
# ╠═02f8b97a-0f6f-4daa-8fd6-e94f2dcb8875
# ╠═759ee2c9-d2cc-4f21-a140-a1e0308892d8
# ╠═7931cec0-610e-4836-b2b8-5fcc20d31b34
# ╠═7013902a-b421-4abb-a1f2-df3fe58a62a1
# ╠═1901353b-62bd-4328-a6f3-13460c7a770d
# ╠═c95d2ef5-eee1-4286-9702-e40c3282b5ae
# ╟─be927f12-e945-4aed-965a-62c0691ad23f
# ╠═007a14b3-583f-4df6-8ef4-eab29b882641
# ╟─215cfb07-96bb-4267-938a-efce8b90153a
# ╠═8b3614b1-17d5-462e-a3a0-eb14df11dccc
# ╠═7b4b5f56-3316-4ae5-b53d-4d588bda7975
# ╠═7886a63b-928a-42a2-a4f5-bfd920eeb2ce
# ╠═e9282a6f-3239-4e9c-b11e-e6d68aefde2d
# ╠═2497e93a-3a45-46ac-a3b4-076047bb8c80
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
