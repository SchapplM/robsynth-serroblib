% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:03
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:55
	% EndTime: 2019-10-10 13:03:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:55
	% EndTime: 2019-10-10 13:03:55
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
	% StartTime: 2019-10-10 13:03:55
	% EndTime: 2019-10-10 13:03:56
	% DurationCPUTime: 0.11s
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
	% StartTime: 2019-10-10 13:03:57
	% EndTime: 2019-10-10 13:03:57
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (83->38), mult. (270->73), div. (0->0), fcn. (211->6), ass. (0->31)
	t199 = sin(qJ(2));
	t202 = cos(qJ(2));
	t223 = pkin(8) + r_i_i_C(3);
	t213 = t223 * t202;
	t224 = -pkin(2) * t199 + t213;
	t201 = cos(qJ(3));
	t203 = cos(qJ(1));
	t221 = t201 * t203;
	t200 = sin(qJ(1));
	t220 = qJD(1) * t200;
	t219 = qJD(1) * t203;
	t218 = qJD(2) * t200;
	t217 = qJD(2) * t202;
	t216 = qJD(2) * t203;
	t215 = qJD(3) * t199;
	t214 = qJD(3) * t202;
	t212 = -qJD(1) + t214;
	t211 = qJD(1) * t202 - qJD(3);
	t198 = sin(qJ(3));
	t210 = r_i_i_C(1) * t198 + r_i_i_C(2) * t201;
	t209 = r_i_i_C(1) * t201 - r_i_i_C(2) * t198 + pkin(2);
	t208 = t212 * t198;
	t207 = -pkin(2) * t202 - t223 * t199 - pkin(1);
	t206 = qJD(2) * t209;
	t205 = t199 * t216 + t211 * t200;
	t204 = -t223 * qJD(2) + t210 * qJD(3);
	t197 = -t211 * t221 + (qJD(2) * t199 * t201 + t208) * t200;
	t196 = t212 * t201 * t200 + (-t199 * t218 + t211 * t203) * t198;
	t195 = t205 * t201 + t203 * t208;
	t194 = t205 * t198 - t212 * t221;
	t1 = [t197 * r_i_i_C(1) + t196 * r_i_i_C(2) - t224 * t218 + (-pkin(7) * t200 + t207 * t203) * qJD(1), (-t203 * t206 - t223 * t220) * t202 + (t204 * t203 + t209 * t220) * t199, t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0, 0, 0; -t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t224 * t216 + (pkin(7) * t203 + t207 * t200) * qJD(1), (-t200 * t206 + t223 * t219) * t202 + (t204 * t200 - t209 * t219) * t199, -t196 * r_i_i_C(1) + t197 * r_i_i_C(2), 0, 0, 0; 0, -t210 * t214 + (-t209 * t199 + t213) * qJD(2), (t198 * t215 - t201 * t217) * r_i_i_C(2) + (-t198 * t217 - t201 * t215) * r_i_i_C(1), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:57
	% EndTime: 2019-10-10 13:03:57
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (261->56), mult. (420->90), div. (0->0), fcn. (332->8), ass. (0->52)
	t235 = cos(qJ(3));
	t227 = t235 * pkin(3) + pkin(2);
	t233 = sin(qJ(2));
	t236 = cos(qJ(2));
	t267 = r_i_i_C(3) + pkin(9) + pkin(8);
	t251 = t267 * t236;
	t232 = sin(qJ(3));
	t266 = pkin(3) * qJD(3);
	t256 = t232 * t266;
	t276 = (-t227 * t233 + t251) * qJD(2) - t236 * t256;
	t237 = cos(qJ(1));
	t230 = qJD(3) + qJD(4);
	t250 = t230 * t236 - qJD(1);
	t274 = t237 * t250;
	t231 = qJ(3) + qJ(4);
	t228 = sin(t231);
	t229 = cos(t231);
	t268 = r_i_i_C(2) * t229;
	t246 = r_i_i_C(1) * t228 + t268;
	t273 = -t246 * t230 - t256;
	t260 = qJD(1) * t236;
	t249 = -t230 + t260;
	t234 = sin(qJ(1));
	t258 = qJD(2) * t233;
	t253 = t234 * t258;
	t272 = t249 * t237 - t253;
	t271 = pkin(3) * t232;
	t270 = r_i_i_C(1) * t229;
	t269 = r_i_i_C(2) * t228;
	t264 = t230 * t233;
	t252 = t237 * t258;
	t240 = t249 * t234 + t252;
	t222 = t240 * t228 - t229 * t274;
	t223 = t228 * t274 + t240 * t229;
	t263 = t222 * r_i_i_C(1) + t223 * r_i_i_C(2);
	t245 = t250 * t234;
	t224 = t272 * t228 + t229 * t245;
	t225 = t228 * t245 - t272 * t229;
	t262 = -t224 * r_i_i_C(1) + t225 * r_i_i_C(2);
	t261 = qJD(1) * t234;
	t259 = qJD(1) * t237;
	t257 = qJD(2) * t236;
	t255 = t235 * t266;
	t254 = pkin(7) + t271;
	t247 = -qJD(3) + t260;
	t244 = (-qJD(3) * t236 + qJD(1)) * t235;
	t243 = t227 - t269 + t270;
	t242 = -t227 * t236 - t267 * t233 - pkin(1);
	t241 = qJD(2) * t243;
	t239 = -t267 * qJD(2) - t273;
	t226 = t264 * t269;
	t1 = [t237 * t255 + t225 * r_i_i_C(1) + t224 * r_i_i_C(2) - t276 * t234 + (-t254 * t234 + t242 * t237) * qJD(1), (-t237 * t241 - t267 * t261) * t236 + (t239 * t237 + t243 * t261) * t233, (t237 * t244 + (t247 * t234 + t252) * t232) * pkin(3) + t263, t263, 0, 0; t234 * t255 - t223 * r_i_i_C(1) + t222 * r_i_i_C(2) + t276 * t237 + (t242 * t234 + t254 * t237) * qJD(1), (-t234 * t241 + t267 * t259) * t236 + (t239 * t234 - t243 * t259) * t233, (t234 * t244 + (-t247 * t237 + t253) * t232) * pkin(3) + t262, t262, 0, 0; 0, t273 * t236 + (-t243 * t233 + t251) * qJD(2), t226 + (-t230 * t270 - t255) * t233 + (-t246 - t271) * t257, -t257 * t268 + t226 + (-t228 * t257 - t229 * t264) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:57
	% EndTime: 2019-10-10 13:03:57
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (578->73), mult. (572->105), div. (0->0), fcn. (452->10), ass. (0->63)
	t247 = qJ(3) + qJ(4);
	t241 = sin(t247);
	t248 = sin(qJ(3));
	t282 = pkin(3) * qJD(3);
	t245 = qJD(3) + qJD(4);
	t288 = pkin(4) * t245;
	t231 = -t241 * t288 - t248 * t282;
	t242 = cos(t247);
	t251 = cos(qJ(3));
	t236 = t251 * pkin(3) + pkin(4) * t242;
	t234 = pkin(2) + t236;
	t249 = sin(qJ(2));
	t252 = cos(qJ(2));
	t283 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
	t267 = t283 * t252;
	t294 = t252 * t231 + (-t234 * t249 + t267) * qJD(2);
	t253 = cos(qJ(1));
	t240 = qJD(5) + t245;
	t266 = t240 * t252 - qJD(1);
	t292 = t253 * t266;
	t243 = qJ(5) + t247;
	t238 = sin(t243);
	t239 = cos(t243);
	t285 = r_i_i_C(2) * t239;
	t263 = r_i_i_C(1) * t238 + t285;
	t291 = -t240 * t263 + t231;
	t275 = qJD(1) * t252;
	t265 = -t240 + t275;
	t250 = sin(qJ(1));
	t273 = qJD(2) * t249;
	t269 = t250 * t273;
	t290 = t253 * t265 - t269;
	t289 = pkin(4) * t241;
	t287 = r_i_i_C(1) * t239;
	t286 = r_i_i_C(2) * t238;
	t235 = pkin(3) * t248 + t289;
	t284 = pkin(7) + t235;
	t280 = t240 * t249;
	t268 = t253 * t273;
	t255 = t250 * t265 + t268;
	t227 = t238 * t255 - t239 * t292;
	t228 = t238 * t292 + t239 * t255;
	t278 = t227 * r_i_i_C(1) + t228 * r_i_i_C(2);
	t260 = t266 * t250;
	t229 = t238 * t290 + t239 * t260;
	t230 = t238 * t260 - t239 * t290;
	t277 = -t229 * r_i_i_C(1) + t230 * r_i_i_C(2);
	t276 = qJD(1) * t250;
	t274 = qJD(1) * t253;
	t272 = qJD(2) * t252;
	t271 = t242 * t288;
	t270 = t240 * t287;
	t264 = -t245 + t275;
	t262 = t235 * t275 + t231;
	t261 = t242 * (-t245 * t252 + qJD(1));
	t259 = t234 - t286 + t287;
	t258 = -t234 * t252 - t249 * t283 - pkin(1);
	t257 = qJD(2) * t259;
	t232 = t251 * t282 + t271;
	t256 = qJD(1) * t236 - t232 * t252 + t235 * t273;
	t254 = -qJD(2) * t283 - t291;
	t233 = t280 * t286;
	t1 = [t230 * r_i_i_C(1) + t229 * r_i_i_C(2) + t253 * t232 - t294 * t250 + (-t250 * t284 + t253 * t258) * qJD(1), (-t253 * t257 - t276 * t283) * t252 + (t253 * t254 + t259 * t276) * t249, t250 * t262 + t253 * t256 + t278, (t253 * t261 + (t250 * t264 + t268) * t241) * pkin(4) + t278, t278, 0; -t228 * r_i_i_C(1) + t227 * r_i_i_C(2) + t250 * t232 + t294 * t253 + (t250 * t258 + t253 * t284) * qJD(1), (-t250 * t257 + t274 * t283) * t252 + (t250 * t254 - t259 * t274) * t249, t250 * t256 - t253 * t262 + t277, (t250 * t261 + (-t253 * t264 + t269) * t241) * pkin(4) + t277, t277, 0; 0, t291 * t252 + (-t249 * t259 + t267) * qJD(2), t233 + (-t232 - t270) * t249 + (-t235 - t263) * t272, t233 + (-t270 - t271) * t249 + (-t263 - t289) * t272, -t272 * t285 + t233 + (-t238 * t272 - t239 * t280) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:03:57
	% EndTime: 2019-10-10 13:03:57
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (800->81), mult. (691->108), div. (0->0), fcn. (542->10), ass. (0->68)
	t256 = qJ(3) + qJ(4);
	t253 = qJ(5) + t256;
	t247 = sin(t253);
	t259 = sin(qJ(1));
	t262 = cos(qJ(1));
	t255 = qJD(3) + qJD(4);
	t249 = qJD(5) + t255;
	t261 = cos(qJ(2));
	t289 = qJD(1) * t261;
	t277 = -t249 + t289;
	t258 = sin(qJ(2));
	t287 = qJD(2) * t258;
	t305 = -t259 * t287 + t277 * t262;
	t311 = t305 * t247;
	t250 = sin(t256);
	t301 = pkin(5) * t247;
	t302 = pkin(4) * t255;
	t237 = -t249 * t301 - t250 * t302;
	t257 = sin(qJ(3));
	t296 = pkin(3) * qJD(3);
	t235 = -t257 * t296 + t237;
	t248 = cos(t253);
	t251 = cos(t256);
	t244 = pkin(4) * t251 + pkin(5) * t248;
	t260 = cos(qJ(3));
	t241 = t260 * pkin(3) + t244;
	t239 = pkin(2) + t241;
	t243 = -pkin(4) * t250 - t301;
	t240 = t257 * pkin(3) - t243;
	t284 = t258 * qJD(6);
	t297 = r_i_i_C(3) + qJ(6) + pkin(10) + pkin(9) + pkin(8);
	t307 = t297 * t261;
	t310 = (-t239 * t258 + t307) * qJD(2) + (pkin(7) + t240) * qJD(1) + t261 * t235 + t284;
	t300 = r_i_i_C(2) * t247;
	t270 = r_i_i_C(1) * t248 + t239 - t300;
	t264 = -t270 * t258 + t307;
	t299 = r_i_i_C(2) * t248;
	t276 = r_i_i_C(1) * t247 + t299;
	t306 = t276 * t249 - t235;
	t303 = -pkin(5) - r_i_i_C(1);
	t294 = t248 * t249;
	t293 = t249 * t258;
	t285 = qJD(2) * t262;
	t266 = t258 * t285 + t277 * t259;
	t278 = t249 * t261 - qJD(1);
	t273 = t248 * t278;
	t231 = t266 * t247 - t262 * t273;
	t232 = t278 * t262 * t247 + t266 * t248;
	t292 = t231 * r_i_i_C(1) + t232 * r_i_i_C(2);
	t272 = t278 * t259;
	t233 = t248 * t272 + t311;
	t234 = t247 * t272 - t248 * t305;
	t291 = -t233 * r_i_i_C(1) + t234 * r_i_i_C(2);
	t290 = qJD(1) * t259;
	t288 = qJD(1) * t262;
	t286 = qJD(2) * t261;
	t283 = r_i_i_C(1) * t294;
	t281 = t297 * t258;
	t275 = t240 * t289 + t235;
	t274 = t243 * t289 - t237;
	t238 = -pkin(5) * t294 - t251 * t302;
	t236 = t260 * t296 - t238;
	t269 = qJD(1) * t241 - t236 * t261 + t240 * t287;
	t268 = qJD(1) * t244 + t238 * t261 - t243 * t287;
	t265 = t236 + (-t239 * t261 - pkin(1) - t281) * qJD(1);
	t263 = qJD(6) * t261 + t306 * t258 + (-t270 * t261 - t281) * qJD(2);
	t242 = t293 * t300;
	t1 = [t234 * r_i_i_C(1) + t233 * r_i_i_C(2) - t310 * t259 + t265 * t262, t263 * t262 - t264 * t290, t275 * t259 + t269 * t262 + t292, -t274 * t259 + t268 * t262 + t292, t231 * pkin(5) + t292, -t258 * t290 + t261 * t285; -t232 * r_i_i_C(1) + t231 * r_i_i_C(2) + t265 * t259 + t310 * t262, t263 * t259 + t264 * t288, t269 * t259 - t275 * t262 + t291, t268 * t259 + t274 * t262 + t291, (-t259 * t273 - t311) * pkin(5) + t291, t258 * t288 + t259 * t286; 0, t264 * qJD(2) - t306 * t261 + t284, t242 + (-t236 - t283) * t258 + (-t240 - t276) * t286, t242 + (t238 - t283) * t258 + (t243 - t276) * t286, t242 + t303 * t248 * t293 + (t303 * t247 - t299) * t286, t287;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end