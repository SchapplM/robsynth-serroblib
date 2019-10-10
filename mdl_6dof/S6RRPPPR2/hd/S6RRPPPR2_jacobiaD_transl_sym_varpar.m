% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPPR2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:19
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:54
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(7) + r_i_i_C(3);
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
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0, 0; 0, -t21, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:54
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (50->16), mult. (90->25), div. (0->0), fcn. (59->6), ass. (0->14)
	t27 = qJ(2) + pkin(9);
	t25 = sin(t27);
	t26 = cos(t27);
	t42 = -r_i_i_C(1) * t26 + r_i_i_C(2) * t25 - cos(qJ(2)) * pkin(2);
	t40 = r_i_i_C(3) + qJ(3) + pkin(7);
	t32 = cos(qJ(1));
	t39 = qJD(1) * t32;
	t37 = -pkin(1) + t42;
	t36 = sin(qJ(2)) * pkin(2) + r_i_i_C(1) * t25 + r_i_i_C(2) * t26;
	t30 = sin(qJ(1));
	t35 = t36 * t30;
	t34 = qJD(2) * t42;
	t33 = t36 * qJD(2);
	t1 = [t32 * qJD(3) + qJD(2) * t35 + (-t30 * t40 + t32 * t37) * qJD(1), qJD(1) * t35 + t32 * t34, t39, 0, 0, 0; t30 * qJD(3) - t32 * t33 + (t30 * t37 + t32 * t40) * qJD(1), t30 * t34 - t36 * t39, qJD(1) * t30, 0, 0, 0; 0, -t33, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:54
	% EndTime: 2019-10-10 09:19:54
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (103->20), mult. (160->28), div. (0->0), fcn. (111->6), ass. (0->17)
	t149 = qJ(2) + pkin(9);
	t147 = sin(t149);
	t148 = cos(t149);
	t165 = r_i_i_C(3) + qJ(4);
	t169 = pkin(3) - r_i_i_C(2);
	t158 = t169 * t147 - t165 * t148 + sin(qJ(2)) * pkin(2);
	t155 = -t158 * qJD(2) + t147 * qJD(4);
	t172 = t155 + (r_i_i_C(1) + qJ(3) + pkin(7)) * qJD(1);
	t171 = -t165 * t147 - t169 * t148 - cos(qJ(2)) * pkin(2);
	t152 = sin(qJ(1));
	t164 = qJD(1) * t152;
	t154 = cos(qJ(1));
	t163 = qJD(1) * t154;
	t162 = qJD(2) * t148;
	t157 = qJD(3) + (-pkin(1) + t171) * qJD(1);
	t156 = t171 * qJD(2) + qJD(4) * t148;
	t1 = [-t172 * t152 + t157 * t154, t156 * t154 + t158 * t164, t163, -t147 * t164 + t154 * t162, 0, 0; t157 * t152 + t172 * t154, t156 * t152 - t158 * t163, t164, t147 * t163 + t152 * t162, 0, 0; 0, t155, 0, qJD(2) * t147, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:55
	% EndTime: 2019-10-10 09:19:55
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (160->27), mult. (252->40), div. (0->0), fcn. (189->8), ass. (0->20)
	t174 = sin(pkin(10));
	t175 = cos(pkin(10));
	t173 = qJ(2) + pkin(9);
	t171 = sin(t173);
	t172 = cos(t173);
	t188 = r_i_i_C(1) * t174 + r_i_i_C(2) * t175 + qJ(4);
	t190 = pkin(3) + r_i_i_C(3) + qJ(5);
	t184 = t190 * t171 - t188 * t172 + sin(qJ(2)) * pkin(2);
	t181 = -t184 * qJD(2) + t171 * qJD(4) + t172 * qJD(5);
	t199 = t181 + (t175 * r_i_i_C(1) - t174 * r_i_i_C(2) + pkin(4) + pkin(7) + qJ(3)) * qJD(1);
	t198 = -t188 * t171 - t190 * t172 - cos(qJ(2)) * pkin(2);
	t178 = sin(qJ(1));
	t194 = qJD(1) * t178;
	t180 = cos(qJ(1));
	t193 = qJD(1) * t180;
	t192 = qJD(2) * t178;
	t191 = qJD(2) * t180;
	t183 = qJD(3) + (-pkin(1) + t198) * qJD(1);
	t182 = t198 * qJD(2) + qJD(4) * t172 - qJD(5) * t171;
	t1 = [-t199 * t178 + t183 * t180, t182 * t180 + t184 * t194, t193, -t171 * t194 + t172 * t191, -t171 * t191 - t172 * t194, 0; t183 * t178 + t199 * t180, t182 * t178 - t184 * t193, t194, t171 * t193 + t172 * t192, -t171 * t192 + t172 * t193, 0; 0, t181, 0, qJD(2) * t171, qJD(2) * t172, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:55
	% EndTime: 2019-10-10 09:19:55
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (310->51), mult. (399->76), div. (0->0), fcn. (313->10), ass. (0->38)
	t225 = qJ(2) + pkin(9);
	t221 = sin(t225);
	t223 = cos(t225);
	t249 = pkin(3) + r_i_i_C(3) + pkin(8) + qJ(5);
	t240 = t249 * t221 + sin(qJ(2)) * pkin(2);
	t246 = pkin(5) * sin(pkin(10)) + qJ(4);
	t250 = t223 * qJD(5);
	t268 = (-t246 * t223 + t240) * qJD(2) - (cos(pkin(10)) * pkin(5) + pkin(4) + qJ(3) + pkin(7)) * qJD(1) - t221 * qJD(4) - t250;
	t224 = pkin(10) + qJ(6);
	t220 = sin(t224);
	t222 = cos(t224);
	t239 = r_i_i_C(1) * t220 + r_i_i_C(2) * t222 + t246;
	t266 = -t239 * t223 + t240;
	t230 = sin(qJ(1));
	t245 = qJD(6) * t221 + qJD(1);
	t265 = t230 * t245;
	t232 = cos(qJ(1));
	t264 = t232 * t245;
	t261 = -t249 * t223 - cos(qJ(2)) * pkin(2);
	t256 = qJD(1) * t230;
	t255 = qJD(1) * t232;
	t254 = qJD(2) * t221;
	t253 = qJD(2) * t230;
	t252 = qJD(2) * t232;
	t251 = qJD(6) * t223;
	t248 = t223 * t253;
	t247 = t223 * t252;
	t244 = -qJD(1) * t221 - qJD(6);
	t238 = qJD(4) + (r_i_i_C(1) * t222 - r_i_i_C(2) * t220) * qJD(6);
	t237 = t244 * t232 - t248;
	t236 = t244 * t230 + t247;
	t235 = qJD(3) + (-t246 * t221 - pkin(1) + t261) * qJD(1);
	t233 = -qJD(5) * t221 + t238 * t223 + (-t239 * t221 + t261) * qJD(2);
	t217 = t236 * t220 + t222 * t264;
	t216 = -t220 * t264 + t236 * t222;
	t215 = t237 * t220 - t222 * t265;
	t214 = t220 * t265 + t237 * t222;
	t1 = [t215 * r_i_i_C(1) + t214 * r_i_i_C(2) + t268 * t230 + t235 * t232, t233 * t232 + t266 * t256, t255, -t221 * t256 + t247, -t221 * t252 - t223 * t256, t216 * r_i_i_C(1) - t217 * r_i_i_C(2); t217 * r_i_i_C(1) + t216 * r_i_i_C(2) + t235 * t230 - t268 * t232, t233 * t230 - t255 * t266, t256, t221 * t255 + t248, -t221 * t253 + t223 * t255, -t214 * r_i_i_C(1) + t215 * r_i_i_C(2); 0, -qJD(2) * t266 + t238 * t221 + t250, 0, t254, qJD(2) * t223, (-t220 * t254 + t222 * t251) * r_i_i_C(2) + (t220 * t251 + t222 * t254) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end