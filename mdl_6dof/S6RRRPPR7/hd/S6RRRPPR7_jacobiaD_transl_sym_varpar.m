% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:27
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
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
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:46
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
	% StartTime: 2019-10-10 11:27:47
	% EndTime: 2019-10-10 11:27:47
	% DurationCPUTime: 0.20s
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
	% StartTime: 2019-10-10 11:27:47
	% EndTime: 2019-10-10 11:27:47
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (164->54), mult. (510->88), div. (0->0), fcn. (427->6), ass. (0->42)
	t248 = sin(qJ(3));
	t251 = cos(qJ(3));
	t281 = r_i_i_C(3) + qJ(4);
	t284 = pkin(3) + r_i_i_C(1);
	t285 = t284 * t248 - t281 * t251;
	t287 = -t285 * qJD(3) + qJD(4) * t248;
	t249 = sin(qJ(2));
	t252 = cos(qJ(2));
	t283 = pkin(8) + r_i_i_C(2);
	t268 = t283 * t252;
	t286 = -pkin(2) * t249 + t268;
	t259 = -t281 * t248 - t284 * t251;
	t256 = -pkin(2) + t259;
	t250 = sin(qJ(1));
	t280 = t250 * t248;
	t279 = t250 * t252;
	t253 = cos(qJ(1));
	t278 = t253 * t248;
	t277 = t253 * t251;
	t276 = qJD(1) * t250;
	t275 = qJD(1) * t253;
	t274 = qJD(2) * t250;
	t273 = qJD(2) * t252;
	t272 = qJD(2) * t253;
	t271 = qJD(3) * t251;
	t270 = qJD(3) * t253;
	t267 = t249 * t274;
	t266 = qJD(3) * t280;
	t265 = t249 * t272;
	t264 = t248 * t270;
	t263 = t251 * t270;
	t262 = t252 * t277 + t280;
	t261 = t248 * t279 + t277;
	t260 = -pkin(2) * t252 - t283 * t249 - pkin(1);
	t257 = t248 * t275 + t250 * t271;
	t255 = qJD(2) * t256;
	t254 = -t283 * qJD(2) - t287;
	t237 = t262 * qJD(1) - t251 * t267 - t252 * t266 - t263;
	t236 = -t248 * t267 - t251 * t276 + t257 * t252 - t264;
	t235 = t252 * t264 + (t252 * t276 + t265) * t251 - t257;
	t234 = t261 * qJD(1) + t248 * t265 - t252 * t263 - t266;
	t1 = [-t261 * qJD(4) - t284 * t237 - t281 * t236 - t286 * t274 + (-t250 * pkin(7) + t260 * t253) * qJD(1), (t253 * t255 - t283 * t276) * t252 + (t254 * t253 - t256 * t276) * t249, t262 * qJD(4) + t284 * t234 - t281 * t235, -t234, 0, 0; -(t250 * t251 - t252 * t278) * qJD(4) - t284 * t235 - t281 * t234 + t286 * t272 + (t253 * pkin(7) + t260 * t250) * qJD(1), (t250 * t255 + t283 * t275) * t252 + (t254 * t250 + t256 * t275) * t249, -(-t251 * t279 + t278) * qJD(4) + t281 * t237 - t284 * t236, t236, 0, 0; 0, t287 * t252 + (t256 * t249 + t268) * qJD(2), -t285 * t273 + (t259 * qJD(3) + t251 * qJD(4)) * t249, t248 * t273 + t249 * t271, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:47
	% EndTime: 2019-10-10 11:27:48
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (271->64), mult. (846->97), div. (0->0), fcn. (752->8), ass. (0->50)
	t262 = sin(qJ(2));
	t265 = cos(qJ(2));
	t290 = t262 * qJD(5);
	t261 = sin(qJ(3));
	t291 = t261 * qJD(4);
	t288 = -r_i_i_C(3) - qJ(5) + pkin(8);
	t304 = t288 * t265;
	t307 = (-t262 * pkin(2) + t304) * qJD(2) + t265 * t291 - t290;
	t264 = cos(qJ(3));
	t260 = cos(pkin(10));
	t259 = sin(pkin(10));
	t281 = t259 * r_i_i_C(2) - pkin(3) - pkin(4);
	t276 = t260 * r_i_i_C(1) - t281;
	t282 = -t259 * r_i_i_C(1) - qJ(4);
	t277 = t260 * r_i_i_C(2) - t282;
	t303 = t276 * t261 - t277 * t264;
	t306 = t303 * qJD(3) - t291;
	t271 = -t277 * t261 - t276 * t264;
	t269 = -pkin(2) + t271;
	t268 = t269 * t262 + t304;
	t263 = sin(qJ(1));
	t301 = t263 * t261;
	t300 = t263 * t265;
	t266 = cos(qJ(1));
	t299 = t264 * t266;
	t298 = qJD(1) * t263;
	t297 = qJD(1) * t266;
	t296 = qJD(2) * t262;
	t295 = qJD(2) * t265;
	t294 = qJD(2) * t266;
	t293 = qJD(3) * t264;
	t292 = qJD(3) * t266;
	t289 = t264 * qJD(4);
	t287 = t263 * t296;
	t286 = qJD(3) * t301;
	t285 = t262 * t294;
	t284 = t261 * t292;
	t283 = t264 * t292;
	t280 = t288 * t262;
	t275 = t265 * t299 + t301;
	t273 = t261 * t297 + t263 * t293;
	t272 = -pkin(2) * t265 - pkin(1) - t280;
	t267 = -t265 * qJD(5) + t306 * t262 + (t269 * t265 - t280) * qJD(2);
	t248 = t275 * qJD(1) - t264 * t287 - t265 * t286 - t283;
	t247 = -t261 * t287 - t264 * t298 + t273 * t265 - t284;
	t246 = t265 * t284 + (t265 * t298 + t285) * t264 - t273;
	t245 = t261 * t285 - t265 * t283 - t286 + (t261 * t300 + t299) * qJD(1);
	t244 = t247 * t260;
	t243 = t246 * t260;
	t1 = [-t266 * t289 - t244 * r_i_i_C(2) + t282 * t247 - t276 * t248 - t307 * t263 + (-pkin(7) * t263 + t272 * t266) * qJD(1), t267 * t266 - t268 * t298, -t243 * r_i_i_C(2) + t275 * qJD(4) + t276 * t245 + t282 * t246, -t245, t262 * t298 - t265 * t294, 0; -t263 * t289 - t243 * r_i_i_C(1) - t277 * t245 + t281 * t246 + t307 * t266 + (t266 * pkin(7) + t272 * t263) * qJD(1), t267 * t263 + t268 * t297, -t244 * r_i_i_C(1) - (t261 * t266 - t264 * t300) * qJD(4) + t277 * t248 + t281 * t247, t247, -t262 * t297 - t263 * t295, 0; 0, t268 * qJD(2) - t306 * t265 - t290, -t303 * t295 + (t271 * qJD(3) + t289) * t262, t261 * t295 + t262 * t293, -t296, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:47
	% EndTime: 2019-10-10 11:27:48
	% DurationCPUTime: 0.85s
	% Computational Cost: add. (553->96), mult. (1296->147), div. (0->0), fcn. (1206->10), ass. (0->62)
	t300 = sin(qJ(2));
	t299 = sin(qJ(3));
	t302 = cos(qJ(3));
	t296 = pkin(10) + qJ(6);
	t295 = cos(t296);
	t294 = sin(t296);
	t323 = sin(pkin(10)) * pkin(5) + qJ(4);
	t312 = -t294 * r_i_i_C(1) - t323;
	t308 = t295 * r_i_i_C(2) - t312;
	t343 = -pkin(3) - cos(pkin(10)) * pkin(5) - pkin(4);
	t320 = t294 * r_i_i_C(2) + t343;
	t311 = t295 * r_i_i_C(1) - t320;
	t345 = t308 * t299 + t311 * t302 + pkin(2);
	t303 = cos(qJ(2));
	t329 = -r_i_i_C(3) - pkin(9) - qJ(5) + pkin(8);
	t351 = t329 * t303;
	t357 = t345 * t300 - t351;
	t356 = (-t300 * pkin(2) + t351) * qJD(2) - t300 * qJD(5);
	t355 = (qJD(3) - qJD(6)) * t300;
	t304 = cos(qJ(1));
	t340 = t302 * t304;
	t301 = sin(qJ(1));
	t341 = t301 * t303;
	t277 = t299 * t341 + t340;
	t339 = t304 * t299;
	t278 = t302 * t341 - t339;
	t317 = t277 * t294 + t278 * t295;
	t318 = t277 * t295 - t278 * t294;
	t353 = (t317 * r_i_i_C(1) + t318 * r_i_i_C(2)) * qJD(6);
	t313 = t294 * t299 + t295 * t302;
	t314 = t294 * t302 - t295 * t299;
	t352 = -t299 * qJD(4) - t329 * qJD(2) + (t314 * r_i_i_C(1) + t313 * r_i_i_C(2)) * qJD(6) + (t311 * t299 - t308 * t302) * qJD(3);
	t333 = qJD(3) * t302;
	t337 = qJD(1) * t304;
	t310 = t299 * t337 + t301 * t333;
	t332 = qJD(3) * t304;
	t325 = t299 * t332;
	t336 = qJD(2) * t300;
	t328 = t301 * t336;
	t338 = qJD(1) * t301;
	t275 = -t299 * t328 - t302 * t338 + t310 * t303 - t325;
	t272 = t275 * t295;
	t342 = t301 * t299;
	t335 = qJD(2) * t303;
	t334 = qJD(2) * t304;
	t327 = qJD(3) * t342;
	t326 = t300 * t334;
	t324 = t302 * t332;
	t319 = -r_i_i_C(1) * (-t313 * t355 + t314 * t335) - r_i_i_C(2) * (t313 * t335 + t314 * t355);
	t279 = -t301 * t302 + t303 * t339;
	t280 = t303 * t340 + t342;
	t316 = t279 * t295 - t280 * t294;
	t315 = t279 * t294 + t280 * t295;
	t309 = -pkin(2) * t303 - t329 * t300 - pkin(1);
	t306 = -qJD(2) * t345 - qJD(5);
	t305 = t352 * t300 + t306 * t303;
	t276 = t280 * qJD(1) - t302 * t328 - t303 * t327 - t324;
	t274 = t303 * t325 + (t303 * t338 + t326) * t302 - t310;
	t273 = t277 * qJD(1) + t299 * t326 - t303 * t324 - t327;
	t269 = t316 * qJD(6) - t273 * t294 - t274 * t295;
	t268 = -t315 * qJD(6) - t273 * t295 + t274 * t294;
	t1 = [-t272 * r_i_i_C(2) - t277 * qJD(4) - t311 * t276 + t312 * t275 + (-t318 * r_i_i_C(1) + t317 * r_i_i_C(2)) * qJD(6) - t356 * t301 + (-t301 * pkin(7) + t309 * t304) * qJD(1), t305 * t304 + t357 * t338, t280 * qJD(4) - t308 * t274 + t311 * t273 + (t315 * r_i_i_C(1) + t316 * r_i_i_C(2)) * qJD(6), -t273, t300 * t338 - t303 * t334, r_i_i_C(1) * t268 - r_i_i_C(2) * t269; t269 * r_i_i_C(1) + t268 * r_i_i_C(2) + t279 * qJD(4) + t343 * t274 - t323 * t273 + t356 * t304 + (pkin(7) * t304 + t309 * t301) * qJD(1), t305 * t301 - t357 * t337, -t272 * r_i_i_C(1) + t278 * qJD(4) + t320 * t275 + t308 * t276 + t353, t275, -t300 * t337 - t301 * t335, (-t276 * t294 + t272) * r_i_i_C(1) + (-t275 * t294 - t276 * t295) * r_i_i_C(2) - t353; 0, t306 * t300 - t352 * t303, (t343 * t299 + t323 * t302) * t335 + (qJD(4) * t302 + (-t323 * t299 + t343 * t302) * qJD(3)) * t300 - t319, t299 * t335 + t300 * t333, -t336, t319;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end