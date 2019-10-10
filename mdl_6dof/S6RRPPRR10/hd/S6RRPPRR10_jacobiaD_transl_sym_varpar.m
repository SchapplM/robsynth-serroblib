% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:52
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:05
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
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:05
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 09:52:06
	% EndTime: 2019-10-10 09:52:06
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (44->20), mult. (134->34), div. (0->0), fcn. (92->4), ass. (0->15)
	t139 = sin(qJ(2));
	t141 = cos(qJ(2));
	t151 = r_i_i_C(3) + qJ(3);
	t153 = pkin(2) - r_i_i_C(2);
	t154 = t153 * t139 - t151 * t141;
	t155 = t154 * qJD(2) - t139 * qJD(3);
	t152 = pkin(7) + r_i_i_C(1);
	t140 = sin(qJ(1));
	t150 = qJD(1) * t140;
	t142 = cos(qJ(1));
	t149 = qJD(1) * t142;
	t148 = qJD(2) * t142;
	t146 = -t151 * t139 - t153 * t141;
	t144 = -pkin(1) + t146;
	t1 = [t155 * t140 + (-t152 * t140 + t144 * t142) * qJD(1), (-t151 * t148 + t153 * t150) * t139 + (-t151 * t150 + (-t153 * qJD(2) + qJD(3)) * t142) * t141, -t139 * t150 + t141 * t148, 0, 0, 0; -t155 * t142 + (t144 * t140 + t152 * t142) * qJD(1), -t154 * t149 + (t146 * qJD(2) + qJD(3) * t141) * t140, t140 * qJD(2) * t141 + t139 * t149, 0, 0, 0; 0, -t155, qJD(2) * t139, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:06
	% EndTime: 2019-10-10 09:52:07
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (73->22), mult. (226->37), div. (0->0), fcn. (170->6), ass. (0->20)
	t162 = sin(pkin(10));
	t163 = cos(pkin(10));
	t164 = sin(qJ(2));
	t166 = cos(qJ(2));
	t175 = r_i_i_C(1) * t162 + r_i_i_C(2) * t163 + qJ(3);
	t177 = pkin(2) + r_i_i_C(3) + qJ(4);
	t172 = t177 * t164 - t175 * t166;
	t168 = -t172 * qJD(2) + t164 * qJD(3) + t166 * qJD(4);
	t184 = (t163 * r_i_i_C(1) - t162 * r_i_i_C(2) + pkin(3) + pkin(7)) * qJD(1) + t168;
	t165 = sin(qJ(1));
	t182 = qJD(1) * t165;
	t167 = cos(qJ(1));
	t181 = qJD(1) * t167;
	t180 = qJD(2) * t164;
	t179 = qJD(2) * t166;
	t178 = qJD(2) * t167;
	t173 = -t175 * t164 - t177 * t166;
	t170 = qJD(1) * (-pkin(1) + t173);
	t169 = t173 * qJD(2) + qJD(3) * t166 - qJD(4) * t164;
	t1 = [-t184 * t165 + t167 * t170, t169 * t167 + t172 * t182, -t164 * t182 + t166 * t178, -t164 * t178 - t166 * t182, 0, 0; t165 * t170 + t184 * t167, t169 * t165 - t172 * t181, t164 * t181 + t165 * t179, -t165 * t180 + t166 * t181, 0, 0; 0, t168, t180, t179, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:06
	% EndTime: 2019-10-10 09:52:07
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (188->46), mult. (373->73), div. (0->0), fcn. (294->8), ass. (0->37)
	t208 = sin(qJ(2));
	t210 = cos(qJ(2));
	t225 = sin(pkin(10)) * pkin(4) + qJ(3);
	t229 = t210 * qJD(4);
	t228 = pkin(2) + r_i_i_C(3) + pkin(8) + qJ(4);
	t238 = t228 * t208;
	t245 = (-t225 * t210 + t238) * qJD(2) - (pkin(7) + cos(pkin(10)) * pkin(4) + pkin(3)) * qJD(1) - t208 * qJD(3) - t229;
	t205 = pkin(10) + qJ(5);
	t203 = sin(t205);
	t204 = cos(t205);
	t218 = r_i_i_C(1) * t203 + r_i_i_C(2) * t204 + t225;
	t243 = -t218 * t210 + t238;
	t209 = sin(qJ(1));
	t224 = qJD(5) * t208 + qJD(1);
	t242 = t209 * t224;
	t211 = cos(qJ(1));
	t241 = t211 * t224;
	t235 = qJD(1) * t209;
	t234 = qJD(1) * t211;
	t233 = qJD(2) * t208;
	t232 = qJD(2) * t210;
	t231 = qJD(2) * t211;
	t230 = qJD(5) * t210;
	t227 = t209 * t232;
	t226 = t210 * t231;
	t223 = -qJD(1) * t208 - qJD(5);
	t221 = t228 * t210;
	t217 = qJD(3) + (r_i_i_C(1) * t204 - r_i_i_C(2) * t203) * qJD(5);
	t216 = t223 * t211 - t227;
	t215 = t223 * t209 + t226;
	t214 = qJD(1) * (-t225 * t208 - pkin(1) - t221);
	t212 = -qJD(4) * t208 + t217 * t210 + (-t218 * t208 - t221) * qJD(2);
	t201 = t215 * t203 + t204 * t241;
	t200 = -t203 * t241 + t215 * t204;
	t199 = t216 * t203 - t204 * t242;
	t198 = t203 * t242 + t216 * t204;
	t1 = [t199 * r_i_i_C(1) + t198 * r_i_i_C(2) + t209 * t245 + t211 * t214, t212 * t211 + t243 * t235, -t208 * t235 + t226, -t208 * t231 - t210 * t235, t200 * r_i_i_C(1) - t201 * r_i_i_C(2), 0; t201 * r_i_i_C(1) + t200 * r_i_i_C(2) + t209 * t214 - t211 * t245, t212 * t209 - t234 * t243, t208 * t234 + t227, -t209 * t233 + t210 * t234, -t198 * r_i_i_C(1) + t199 * r_i_i_C(2), 0; 0, -qJD(2) * t243 + t217 * t208 + t229, t233, t232, (-t203 * t233 + t204 * t230) * r_i_i_C(2) + (t203 * t230 + t204 * t233) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:07
	% EndTime: 2019-10-10 09:52:07
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (420->60), mult. (519->87), div. (0->0), fcn. (413->10), ass. (0->51)
	t242 = sin(qJ(2));
	t244 = cos(qJ(2));
	t240 = pkin(10) + qJ(5);
	t237 = cos(t240);
	t282 = pkin(5) * t237;
	t255 = qJD(5) * t282 + qJD(3);
	t265 = t244 * qJD(4);
	t236 = sin(t240);
	t275 = qJ(3) + pkin(5) * t236 + sin(pkin(10)) * pkin(4);
	t264 = pkin(2) + r_i_i_C(3) + pkin(9) + pkin(8) + qJ(4);
	t285 = t264 * t242;
	t293 = (-t275 * t244 + t285) * qJD(2) - (pkin(7) + t282 + cos(pkin(10)) * pkin(4) + pkin(3)) * qJD(1) - t255 * t242 - t265;
	t238 = qJ(6) + t240;
	t234 = sin(t238);
	t235 = cos(t238);
	t292 = r_i_i_C(1) * t234 + r_i_i_C(2) * t235;
	t291 = r_i_i_C(1) * t235 - r_i_i_C(2) * t234;
	t253 = t275 + t292;
	t289 = -t253 * t244 + t285;
	t243 = sin(qJ(1));
	t241 = qJD(5) + qJD(6);
	t261 = t241 * t242 + qJD(1);
	t288 = t243 * t261;
	t245 = cos(qJ(1));
	t287 = t245 * t261;
	t272 = qJD(1) * t242;
	t260 = -t241 - t272;
	t268 = qJD(2) * t244;
	t263 = t243 * t268;
	t251 = t260 * t245 - t263;
	t225 = t234 * t288 + t251 * t235;
	t226 = t251 * t234 - t235 * t288;
	t274 = -t225 * r_i_i_C(1) + t226 * r_i_i_C(2);
	t267 = qJD(2) * t245;
	t262 = t244 * t267;
	t250 = t260 * t243 + t262;
	t227 = -t234 * t287 + t250 * t235;
	t228 = t250 * t234 + t235 * t287;
	t273 = t227 * r_i_i_C(1) - t228 * r_i_i_C(2);
	t271 = qJD(1) * t243;
	t270 = qJD(1) * t245;
	t269 = qJD(2) * t242;
	t266 = qJD(5) * t236;
	t259 = qJD(5) + t272;
	t257 = t264 * t244;
	t254 = t236 * (-qJD(5) * t242 - qJD(1));
	t252 = t292 * t241 * t244 + t291 * t269;
	t249 = t291 * t241 + t255;
	t247 = -pkin(5) * t266 + (-t275 * t242 - pkin(1) - t257) * qJD(1);
	t246 = -qJD(4) * t242 + t249 * t244 + (-t253 * t242 - t257) * qJD(2);
	t1 = [t226 * r_i_i_C(1) + t225 * r_i_i_C(2) + t293 * t243 + t247 * t245, t246 * t245 + t289 * t271, -t242 * t271 + t262, -t242 * t267 - t244 * t271, (t245 * t254 + (-t259 * t243 + t262) * t237) * pkin(5) + t273, t273; t228 * r_i_i_C(1) + t227 * r_i_i_C(2) + t247 * t243 - t293 * t245, t246 * t243 - t270 * t289, t242 * t270 + t263, -t243 * t269 + t244 * t270, (t243 * t254 + (t259 * t245 + t263) * t237) * pkin(5) + t274, t274; 0, -qJD(2) * t289 + t249 * t242 + t265, t269, t268, (t237 * t269 + t244 * t266) * pkin(5) + t252, t252;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end