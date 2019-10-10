% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:23
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
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
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
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
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:25
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
	% StartTime: 2019-10-10 09:23:25
	% EndTime: 2019-10-10 09:23:25
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (73->22), mult. (226->37), div. (0->0), fcn. (170->6), ass. (0->20)
	t162 = sin(pkin(9));
	t163 = cos(pkin(9));
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
	% StartTime: 2019-10-10 09:23:25
	% EndTime: 2019-10-10 09:23:25
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (107->41), mult. (334->63), div. (0->0), fcn. (266->6), ass. (0->30)
	t176 = sin(qJ(2));
	t178 = cos(qJ(2));
	t189 = pkin(2) + r_i_i_C(2) + qJ(4);
	t174 = sin(pkin(9));
	t175 = cos(pkin(9));
	t199 = r_i_i_C(3) + qJ(5);
	t200 = pkin(4) + r_i_i_C(1);
	t202 = t200 * t174 - t199 * t175 + qJ(3);
	t181 = -t189 * t176 + t202 * t178;
	t185 = -qJD(5) * t175 + qJD(3);
	t206 = (t189 * qJD(2) - t185) * t176 - (pkin(3) + pkin(7)) * qJD(1) - (qJ(3) * qJD(2) + qJD(4)) * t178;
	t177 = sin(qJ(1));
	t198 = t175 * t177;
	t179 = cos(qJ(1));
	t197 = t175 * t179;
	t196 = t177 * t174;
	t195 = t179 * t174;
	t194 = qJD(1) * t177;
	t193 = qJD(1) * t179;
	t192 = qJD(2) * t176;
	t191 = qJD(2) * t178;
	t190 = qJD(2) * t179;
	t188 = t177 * t191;
	t187 = t178 * t190;
	t184 = t189 * t178;
	t182 = t174 * qJD(5) + (-qJ(3) * t176 - pkin(1) - t184) * qJD(1);
	t180 = -qJD(4) * t176 + t185 * t178 + (-t176 * t202 - t184) * qJD(2);
	t172 = -t175 * t187 + (t176 * t198 + t195) * qJD(1);
	t170 = -t175 * t188 + (-t176 * t197 + t196) * qJD(1);
	t1 = [t200 * (-t174 * t188 + (-t176 * t195 - t198) * qJD(1)) - t199 * t170 + t182 * t179 + t206 * t177, t180 * t179 - t181 * t194, -t176 * t194 + t187, -t176 * t190 - t178 * t194, t172, 0; t200 * (t174 * t187 + (-t176 * t196 + t197) * qJD(1)) + t199 * t172 + t182 * t177 - t206 * t179, t180 * t177 + t181 * t193, t176 * t193 + t188, -t177 * t192 + t178 * t193, t170, 0; 0, t181 * qJD(2) + t178 * qJD(4) + t185 * t176, t192, t191, -t175 * t192, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:25
	% EndTime: 2019-10-10 09:23:26
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (222->69), mult. (707->110), div. (0->0), fcn. (637->8), ass. (0->48)
	t249 = cos(qJ(2));
	t243 = sin(pkin(9));
	t244 = cos(pkin(9));
	t245 = sin(qJ(6));
	t248 = cos(qJ(6));
	t281 = pkin(4) + pkin(5);
	t257 = t248 * r_i_i_C(1) - t245 * r_i_i_C(2) + t281;
	t258 = t245 * r_i_i_C(1) + t248 * r_i_i_C(2) + qJ(5);
	t283 = t257 * t243 - t258 * t244 + qJ(3);
	t246 = sin(qJ(2));
	t267 = r_i_i_C(3) + pkin(8) - qJ(4) - pkin(2);
	t286 = t267 * t246;
	t252 = t283 * t249 + t286;
	t271 = t249 * qJD(4);
	t288 = (t249 * qJ(3) + t286) * qJD(2) + (pkin(3) + pkin(7)) * qJD(1) + t246 * qJD(3) + t271;
	t247 = sin(qJ(1));
	t273 = qJD(2) * t249;
	t269 = t247 * t273;
	t250 = cos(qJ(1));
	t275 = qJD(1) * t250;
	t285 = t246 * t275 + t269;
	t280 = t246 * t250;
	t279 = t247 * t243;
	t278 = t247 * t244;
	t276 = qJD(1) * t247;
	t274 = qJD(2) * t246;
	t272 = qJD(2) * t250;
	t268 = t249 * t272;
	t264 = t267 * t249;
	t238 = t243 * t250 + t246 * t278;
	t239 = t244 * t250 - t246 * t279;
	t263 = t238 * t248 - t239 * t245;
	t262 = t238 * t245 + t239 * t248;
	t261 = t243 * t248 - t244 * t245;
	t260 = -t243 * t245 - t244 * t248;
	t237 = t243 * t280 + t278;
	t256 = qJD(1) * (-qJ(3) * t246 - pkin(1) + t264);
	t255 = t260 * r_i_i_C(1) - t261 * r_i_i_C(2);
	t253 = -t244 * qJD(5) + t255 * qJD(6) + qJD(3);
	t251 = -t246 * qJD(4) + t253 * t249 + (-t246 * t283 + t264) * qJD(2);
	t236 = -t244 * t280 + t279;
	t235 = t239 * qJD(1) + t243 * t268;
	t234 = t238 * qJD(1) - t244 * t268;
	t233 = t237 * qJD(1) + t243 * t269;
	t232 = t243 * t276 - t285 * t244;
	t231 = t234 * t245 + t235 * t248 + (t236 * t248 - t237 * t245) * qJD(6);
	t230 = t234 * t248 - t235 * t245 + (-t236 * t245 - t237 * t248) * qJD(6);
	t1 = [t238 * qJD(5) - t258 * t232 - t257 * t233 + (t263 * r_i_i_C(1) - t262 * r_i_i_C(2)) * qJD(6) + t250 * t256 - t288 * t247, t251 * t250 - t252 * t276, -t246 * t276 + t268, -t246 * t272 - t249 * t276, t234, r_i_i_C(1) * t230 - r_i_i_C(2) * t231; t231 * r_i_i_C(1) + t230 * r_i_i_C(2) + t234 * qJ(5) + t236 * qJD(5) + t281 * t235 + t247 * t256 + t288 * t250, t251 * t247 + t252 * t275, t285, -t247 * t274 + t249 * t275, t232, (t232 * t248 - t233 * t245) * r_i_i_C(1) + (-t232 * t245 - t233 * t248) * r_i_i_C(2) + (t262 * r_i_i_C(1) + t263 * r_i_i_C(2)) * qJD(6); 0, t252 * qJD(2) + t253 * t246 + t271, t274, t273, -t244 * t274, (t261 * r_i_i_C(1) + t260 * r_i_i_C(2)) * t249 * qJD(6) + t255 * t274;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end