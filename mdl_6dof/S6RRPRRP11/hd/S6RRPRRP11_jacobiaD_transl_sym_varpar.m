% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRP11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP11_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP11_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP11_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
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
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
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
	% StartTime: 2019-10-10 10:46:44
	% EndTime: 2019-10-10 10:46:44
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
	% StartTime: 2019-10-10 10:46:44
	% EndTime: 2019-10-10 10:46:44
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (101->43), mult. (318->73), div. (0->0), fcn. (248->6), ass. (0->34)
	t199 = sin(qJ(2));
	t202 = cos(qJ(2));
	t216 = pkin(2) + pkin(8) + r_i_i_C(3);
	t213 = t216 * t199;
	t229 = (-qJ(3) * t202 + t213) * qJD(2) - t199 * qJD(3);
	t201 = cos(qJ(4));
	t211 = qJD(4) * t199 + qJD(1);
	t227 = t201 * t211;
	t198 = sin(qJ(4));
	t226 = t211 * t198;
	t225 = pkin(3) + pkin(7);
	t200 = sin(qJ(1));
	t223 = qJD(1) * t200;
	t203 = cos(qJ(1));
	t222 = qJD(1) * t203;
	t221 = qJD(2) * t199;
	t220 = qJD(2) * t202;
	t219 = qJD(2) * t203;
	t218 = qJD(4) * t202;
	t215 = t200 * t220;
	t214 = t202 * t219;
	t212 = t216 * t202;
	t210 = -qJD(1) * t199 - qJD(4);
	t209 = t210 * t203;
	t208 = r_i_i_C(1) * t198 + r_i_i_C(2) * t201 + qJ(3);
	t207 = -qJ(3) * t199 - pkin(1) - t212;
	t206 = qJD(3) + (r_i_i_C(1) * t201 - r_i_i_C(2) * t198) * qJD(4);
	t205 = t210 * t200 + t214;
	t204 = t208 * t202 - t213;
	t197 = t205 * t198 + t203 * t227;
	t196 = t205 * t201 - t203 * t226;
	t195 = -t200 * t227 + (t209 - t215) * t198;
	t194 = t201 * t209 + (-t201 * t220 + t226) * t200;
	t1 = [t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t229 * t200 + (-t225 * t200 + t207 * t203) * qJD(1), (-t208 * t219 + t216 * t223) * t199 + (-t208 * t223 + (-t216 * qJD(2) + t206) * t203) * t202, -t199 * t223 + t214, t196 * r_i_i_C(1) - t197 * r_i_i_C(2), 0, 0; t197 * r_i_i_C(1) + t196 * r_i_i_C(2) - t229 * t203 + (t207 * t200 + t225 * t203) * qJD(1), t204 * t222 + (t206 * t202 + (-t208 * t199 - t212) * qJD(2)) * t200, t199 * t222 + t215, -t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0, 0; 0, t204 * qJD(2) + t206 * t199, t221, (-t198 * t221 + t201 * t218) * r_i_i_C(2) + (t198 * t218 + t201 * t221) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:44
	% EndTime: 2019-10-10 10:46:44
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (279->55), mult. (490->85), div. (0->0), fcn. (384->8), ass. (0->48)
	t236 = sin(qJ(2));
	t239 = cos(qJ(2));
	t238 = cos(qJ(4));
	t271 = t238 * pkin(4);
	t250 = qJD(4) * t271 + qJD(3);
	t259 = pkin(2) + r_i_i_C(3) + pkin(9) + pkin(8);
	t252 = t259 * t236;
	t235 = sin(qJ(4));
	t256 = pkin(4) * t235 + qJ(3);
	t283 = (-t256 * t239 + t252) * qJD(2) - (pkin(7) + pkin(3) + t271) * qJD(1) - t250 * t236;
	t234 = qJ(4) + qJ(5);
	t231 = sin(t234);
	t232 = cos(t234);
	t281 = r_i_i_C(1) * t231 + r_i_i_C(2) * t232;
	t280 = r_i_i_C(1) * t232 - r_i_i_C(2) * t231;
	t237 = sin(qJ(1));
	t233 = qJD(4) + qJD(5);
	t255 = t233 * t236 + qJD(1);
	t279 = t237 * t255;
	t240 = cos(qJ(1));
	t278 = t240 * t255;
	t266 = qJD(1) * t236;
	t254 = -t233 - t266;
	t262 = qJD(2) * t239;
	t258 = t237 * t262;
	t246 = t254 * t240 - t258;
	t223 = t231 * t279 + t246 * t232;
	t224 = t246 * t231 - t232 * t279;
	t268 = -t223 * r_i_i_C(1) + t224 * r_i_i_C(2);
	t261 = qJD(2) * t240;
	t257 = t239 * t261;
	t245 = t254 * t237 + t257;
	t225 = -t231 * t278 + t245 * t232;
	t226 = t245 * t231 + t232 * t278;
	t267 = t225 * r_i_i_C(1) - t226 * r_i_i_C(2);
	t265 = qJD(1) * t237;
	t264 = qJD(1) * t240;
	t263 = qJD(2) * t236;
	t260 = qJD(4) * t235;
	t253 = qJD(4) + t266;
	t251 = t259 * t239;
	t249 = (-qJD(4) * t236 - qJD(1)) * t235;
	t248 = t281 * t233 * t239 + t280 * t263;
	t247 = t256 + t281;
	t244 = t280 * t233 + t250;
	t243 = t247 * t239 - t252;
	t242 = -pkin(4) * t260 + (-t256 * t236 - pkin(1) - t251) * qJD(1);
	t1 = [t224 * r_i_i_C(1) + t223 * r_i_i_C(2) + t283 * t237 + t242 * t240, (-t247 * t261 + t259 * t265) * t236 + (-t247 * t265 + (-t259 * qJD(2) + t244) * t240) * t239, -t236 * t265 + t257, (t240 * t249 + (-t253 * t237 + t257) * t238) * pkin(4) + t267, t267, 0; t226 * r_i_i_C(1) + t225 * r_i_i_C(2) + t242 * t237 - t283 * t240, t243 * t264 + (t244 * t239 + (-t247 * t236 - t251) * qJD(2)) * t237, t236 * t264 + t258, (t253 * t240 * t238 + (t238 * t262 + t249) * t237) * pkin(4) + t268, t268, 0; 0, t243 * qJD(2) + t244 * t236, t263, (t238 * t263 + t239 * t260) * pkin(4) + t248, t248, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:44
	% EndTime: 2019-10-10 10:46:44
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (396->66), mult. (592->93), div. (0->0), fcn. (462->8), ass. (0->56)
	t241 = qJ(4) + qJ(5);
	t236 = sin(t241);
	t237 = cos(t241);
	t244 = sin(qJ(1));
	t240 = qJD(4) + qJD(5);
	t243 = sin(qJ(2));
	t262 = t240 * t243 + qJD(1);
	t256 = t262 * t244;
	t247 = cos(qJ(1));
	t276 = qJD(1) * t243;
	t261 = t240 + t276;
	t246 = cos(qJ(2));
	t272 = qJD(2) * t246;
	t264 = t244 * t272;
	t287 = t261 * t247 + t264;
	t297 = -t236 * t256 + t287 * t237;
	t245 = cos(qJ(4));
	t233 = t245 * pkin(4) + pkin(5) * t237;
	t281 = pkin(4) * qJD(4);
	t285 = pkin(5) * t240;
	t228 = t237 * t285 + t245 * t281;
	t269 = qJD(3) + t228;
	t270 = t246 * qJD(6);
	t242 = sin(qJ(4));
	t232 = t242 * pkin(4) + pkin(5) * t236;
	t279 = qJ(3) + t232;
	t268 = pkin(2) + r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
	t289 = t268 * t243;
	t296 = (-t279 * t246 + t289) * qJD(2) - (pkin(7) + pkin(3) + t233) * qJD(1) - t269 * t243 - t270;
	t273 = qJD(2) * t243;
	t280 = t240 * t246;
	t295 = t236 * t280 + t237 * t273;
	t283 = r_i_i_C(2) * t237;
	t255 = r_i_i_C(1) * t236 + t279 + t283;
	t293 = -t255 * t246 + t289;
	t291 = t247 * t262;
	t284 = r_i_i_C(2) * t236;
	t224 = -t236 * t287 - t237 * t256;
	t278 = r_i_i_C(1) * t297 + t224 * r_i_i_C(2);
	t271 = qJD(2) * t247;
	t263 = t246 * t271;
	t251 = -t261 * t244 + t263;
	t225 = -t236 * t291 + t251 * t237;
	t226 = t251 * t236 + t237 * t291;
	t277 = t225 * r_i_i_C(1) - t226 * r_i_i_C(2);
	t275 = qJD(1) * t244;
	t274 = qJD(1) * t247;
	t266 = t295 * r_i_i_C(1) + t280 * t283;
	t259 = t268 * t246;
	t257 = t233 * t276 + t228;
	t254 = (r_i_i_C(1) * t237 - t284) * t240 + t269;
	t227 = -t236 * t285 - t242 * t281;
	t253 = -qJD(1) * t232 + t227 * t243 + t233 * t272;
	t250 = t227 + (-t279 * t243 - pkin(1) - t259) * qJD(1);
	t248 = -qJD(6) * t243 + t254 * t246 + (-t255 * t243 - t259) * qJD(2);
	t1 = [t224 * r_i_i_C(1) - r_i_i_C(2) * t297 + t296 * t244 + t250 * t247, t248 * t247 + t293 * t275, -t243 * t275 + t263, -t257 * t244 + t253 * t247 + t277, t225 * pkin(5) + t277, -t243 * t271 - t246 * t275; t226 * r_i_i_C(1) + t225 * r_i_i_C(2) + t250 * t244 - t296 * t247, t248 * t244 - t274 * t293, t243 * t274 + t264, t253 * t244 + t257 * t247 + t278, t297 * pkin(5) + t278, -t244 * t273 + t246 * t274; 0, -qJD(2) * t293 + t254 * t243 + t270, t273, -t246 * t227 + (t233 - t284) * t273 + t266, t295 * pkin(5) - t273 * t284 + t266, t272;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end