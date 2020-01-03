% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPR7
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPPR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(6) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (50->16), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->14)
	t27 = qJ(2) + pkin(8);
	t25 = sin(t27);
	t26 = cos(t27);
	t42 = -r_i_i_C(1) * t26 + r_i_i_C(2) * t25 - cos(qJ(2)) * pkin(2);
	t40 = r_i_i_C(3) + qJ(3) + pkin(6);
	t32 = cos(qJ(1));
	t39 = qJD(1) * t32;
	t37 = -pkin(1) + t42;
	t36 = sin(qJ(2)) * pkin(2) + r_i_i_C(1) * t25 + r_i_i_C(2) * t26;
	t30 = sin(qJ(1));
	t35 = t36 * t30;
	t34 = qJD(2) * t42;
	t33 = t36 * qJD(2);
	t1 = [t32 * qJD(3) + qJD(2) * t35 + (-t40 * t30 + t37 * t32) * qJD(1), qJD(1) * t35 + t32 * t34, t39, 0, 0; t30 * qJD(3) - t32 * t33 + (t37 * t30 + t40 * t32) * qJD(1), t30 * t34 - t36 * t39, qJD(1) * t30, 0, 0; 0, -t33, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:16
	% EndTime: 2019-12-31 19:37:16
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (103->20), mult. (160->28), div. (0->0), fcn. (111->6), ass. (0->17)
	t149 = qJ(2) + pkin(8);
	t147 = sin(t149);
	t148 = cos(t149);
	t165 = r_i_i_C(3) + qJ(4);
	t169 = pkin(3) - r_i_i_C(2);
	t158 = t169 * t147 - t165 * t148 + sin(qJ(2)) * pkin(2);
	t155 = -t158 * qJD(2) + t147 * qJD(4);
	t172 = t155 + (r_i_i_C(1) + qJ(3) + pkin(6)) * qJD(1);
	t171 = -t165 * t147 - t169 * t148 - cos(qJ(2)) * pkin(2);
	t152 = sin(qJ(1));
	t164 = qJD(1) * t152;
	t154 = cos(qJ(1));
	t163 = qJD(1) * t154;
	t162 = qJD(2) * t148;
	t157 = qJD(3) + (-pkin(1) + t171) * qJD(1);
	t156 = t171 * qJD(2) + qJD(4) * t148;
	t1 = [-t172 * t152 + t157 * t154, t156 * t154 + t158 * t164, t163, -t147 * t164 + t154 * t162, 0; t157 * t152 + t172 * t154, t156 * t152 - t158 * t163, t164, t147 * t163 + t152 * t162, 0; 0, t155, 0, qJD(2) * t147, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:16
	% EndTime: 2019-12-31 19:37:16
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (204->44), mult. (344->68), div. (0->0), fcn. (267->8), ass. (0->35)
	t208 = qJ(2) + pkin(8);
	t206 = sin(t208);
	t207 = cos(t208);
	t230 = pkin(3) + pkin(7) + r_i_i_C(3);
	t222 = t230 * t206 + sin(qJ(2)) * pkin(2);
	t249 = (-qJ(4) * t207 + t222) * qJD(2) - (pkin(4) + qJ(3) + pkin(6)) * qJD(1) - t206 * qJD(4);
	t210 = sin(qJ(5));
	t213 = cos(qJ(5));
	t223 = r_i_i_C(1) * t210 + r_i_i_C(2) * t213 + qJ(4);
	t247 = -t223 * t207 + t222;
	t226 = qJD(5) * t206 + qJD(1);
	t246 = t213 * t226;
	t244 = t226 * t210;
	t243 = -t230 * t207 - cos(qJ(2)) * pkin(2);
	t212 = sin(qJ(1));
	t237 = qJD(1) * t212;
	t215 = cos(qJ(1));
	t236 = qJD(1) * t215;
	t235 = qJD(2) * t206;
	t234 = qJD(2) * t207;
	t233 = qJD(2) * t213;
	t232 = qJD(5) * t207;
	t229 = t212 * t234;
	t228 = t215 * t234;
	t225 = -qJD(1) * t206 - qJD(5);
	t224 = t225 * t215;
	t220 = qJD(4) + (r_i_i_C(1) * t213 - r_i_i_C(2) * t210) * qJD(5);
	t219 = t225 * t212 + t228;
	t218 = qJD(3) + (-qJ(4) * t206 - pkin(1) + t243) * qJD(1);
	t216 = t220 * t207 + (-t223 * t206 + t243) * qJD(2);
	t204 = t219 * t210 + t215 * t246;
	t203 = t219 * t213 - t215 * t244;
	t202 = -t212 * t246 + (t224 - t229) * t210;
	t201 = t213 * t224 + (-t207 * t233 + t244) * t212;
	t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t249 * t212 + t218 * t215, t216 * t215 + t247 * t237, t236, -t206 * t237 + t228, t203 * r_i_i_C(1) - t204 * r_i_i_C(2); t204 * r_i_i_C(1) + t203 * r_i_i_C(2) + t218 * t212 - t249 * t215, t216 * t212 - t236 * t247, t237, t206 * t236 + t229, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2); 0, -qJD(2) * t247 + t220 * t206, 0, t235, (-t210 * t235 + t213 * t232) * r_i_i_C(2) + (t206 * t233 + t210 * t232) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end