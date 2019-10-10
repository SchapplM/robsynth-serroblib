% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:22
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:22:13
	% EndTime: 2019-10-10 11:22:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:22:13
	% EndTime: 2019-10-10 11:22:13
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
	% StartTime: 2019-10-10 11:22:13
	% EndTime: 2019-10-10 11:22:13
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
	% StartTime: 2019-10-10 11:22:14
	% EndTime: 2019-10-10 11:22:15
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
	% StartTime: 2019-10-10 11:22:14
	% EndTime: 2019-10-10 11:22:15
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (189->46), mult. (370->74), div. (0->0), fcn. (290->8), ass. (0->40)
	t213 = cos(qJ(3));
	t242 = t213 * pkin(3);
	t205 = pkin(2) + t242;
	t211 = sin(qJ(2));
	t214 = cos(qJ(2));
	t229 = qJD(3) * t214 - qJD(1);
	t234 = t211 * qJD(4);
	t210 = sin(qJ(3));
	t243 = pkin(3) * t210;
	t241 = r_i_i_C(3) + qJ(4) + pkin(8);
	t246 = t241 * t214;
	t249 = pkin(7) * qJD(1) + (-t205 * t211 + t246) * qJD(2) - t229 * t243 + t234;
	t208 = qJ(3) + pkin(10);
	t206 = sin(t208);
	t207 = cos(t208);
	t227 = r_i_i_C(1) * t207 - r_i_i_C(2) * t206;
	t223 = t205 + t227;
	t218 = -t223 * t211 + t246;
	t212 = sin(qJ(1));
	t224 = t229 * t212;
	t215 = cos(qJ(1));
	t225 = t229 * t215;
	t228 = qJD(1) * t214 - qJD(3);
	t237 = qJD(2) * t212;
	t245 = -t211 * t237 + t228 * t215;
	t239 = qJD(1) * t212;
	t238 = qJD(1) * t215;
	t236 = qJD(2) * t215;
	t235 = qJD(3) * t211;
	t232 = t241 * t211;
	t222 = r_i_i_C(1) * t206 + r_i_i_C(2) * t207 + t243;
	t221 = t222 * t214;
	t219 = t211 * t236 + t228 * t212;
	t217 = qJD(3) * t242 + (-t205 * t214 - pkin(1) - t232) * qJD(1);
	t216 = qJD(4) * t214 + t222 * t235 + (-t223 * t214 - t232) * qJD(2);
	t204 = t206 * t224 - t207 * t245;
	t203 = t245 * t206 + t207 * t224;
	t202 = t206 * t225 + t219 * t207;
	t201 = t219 * t206 - t207 * t225;
	t1 = [t204 * r_i_i_C(1) + t203 * r_i_i_C(2) - t249 * t212 + t217 * t215, t216 * t215 - t218 * t239, t201 * r_i_i_C(1) + t202 * r_i_i_C(2) + (t219 * t210 - t213 * t225) * pkin(3), -t211 * t239 + t214 * t236, 0, 0; -t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t217 * t212 + t249 * t215, t216 * t212 + t218 * t238, -t203 * r_i_i_C(1) + t204 * r_i_i_C(2) + (-t210 * t245 - t213 * t224) * pkin(3), t211 * t238 + t214 * t237, 0, 0; 0, t218 * qJD(2) - qJD(3) * t221 + t234, (-t227 - t242) * t235 - qJD(2) * t221, qJD(2) * t211, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:22:15
	% EndTime: 2019-10-10 11:22:16
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (354->67), mult. (610->99), div. (0->0), fcn. (506->8), ass. (0->51)
	t263 = cos(qJ(3));
	t301 = t263 * pkin(3);
	t255 = pkin(2) + t301;
	t260 = sin(qJ(3));
	t261 = sin(qJ(2));
	t264 = cos(qJ(2));
	t258 = qJ(3) + pkin(10);
	t256 = sin(t258);
	t300 = r_i_i_C(2) + qJ(4) + pkin(8);
	t273 = t300 * qJD(2) + t256 * qJD(5);
	t298 = pkin(3) * qJD(3);
	t302 = pkin(3) * t260;
	t311 = (-t260 * t298 + t273) * t264 + (pkin(7) + t302) * qJD(1) - (qJD(2) * t255 - qJD(4)) * t261;
	t257 = cos(t258);
	t299 = r_i_i_C(3) + qJ(5);
	t303 = r_i_i_C(1) + pkin(4);
	t304 = t303 * t256 - t299 * t257 + t302;
	t309 = t304 * qJD(3) - t273;
	t272 = -t299 * t256 - t303 * t257;
	t270 = -t255 + t272;
	t308 = t270 * t261 + t300 * t264;
	t262 = sin(qJ(1));
	t297 = t262 * t264;
	t265 = cos(qJ(1));
	t296 = t265 * t257;
	t295 = qJD(1) * t262;
	t294 = qJD(1) * t264;
	t293 = qJD(1) * t265;
	t292 = qJD(2) * t261;
	t291 = qJD(2) * t264;
	t290 = qJD(2) * t265;
	t289 = qJD(3) * t262;
	t288 = qJD(3) * t265;
	t286 = t257 * qJD(5);
	t283 = t262 * t292;
	t282 = t261 * t290;
	t281 = t256 * t289;
	t280 = t256 * t288;
	t279 = t257 * t288;
	t278 = -qJD(3) + t294;
	t276 = (-qJD(3) * t264 + qJD(1)) * t263;
	t275 = t262 * t256 + t264 * t296;
	t271 = t256 * t293 + t257 * t289;
	t268 = qJD(2) * t270 + qJD(4);
	t267 = t263 * t298 - t286 + (-t255 * t264 - t300 * t261 - pkin(1)) * qJD(1);
	t266 = t309 * t261 + t268 * t264;
	t244 = qJD(1) * t275 - t257 * t283 - t264 * t281 - t279;
	t243 = -t256 * t283 - t257 * t295 + t264 * t271 - t280;
	t242 = t264 * t280 + (t262 * t294 + t282) * t257 - t271;
	t241 = t256 * t282 - t264 * t279 - t281 + (t256 * t297 + t296) * qJD(1);
	t1 = [-t299 * t243 - t303 * t244 - t311 * t262 + t267 * t265, t266 * t265 - t308 * t295, t275 * qJD(5) - t299 * t242 + t303 * t241 + (t265 * t276 + (t262 * t278 + t282) * t260) * pkin(3), -t261 * t295 + t264 * t290, -t241, 0; -t299 * t241 - t303 * t242 + t267 * t262 + t311 * t265, t266 * t262 + t308 * t293, -(t265 * t256 - t257 * t297) * qJD(5) + t299 * t244 - t303 * t243 + (t262 * t276 + (-t265 * t278 + t283) * t260) * pkin(3), t261 * t293 + t262 * t291, t243, 0; 0, t268 * t261 - t309 * t264, -t304 * t291 + (t286 + (t272 - t301) * qJD(3)) * t261, t292, t261 * qJD(3) * t257 + t256 * t291, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:22:15
	% EndTime: 2019-10-10 11:22:16
	% DurationCPUTime: 0.86s
	% Computational Cost: add. (772->106), mult. (1258->159), div. (0->0), fcn. (1161->10), ass. (0->69)
	t304 = sin(qJ(2));
	t307 = cos(qJ(3));
	t352 = pkin(3) * t307;
	t297 = pkin(2) + t352;
	t300 = qJ(3) + pkin(10);
	t298 = sin(t300);
	t299 = cos(t300);
	t306 = cos(qJ(6));
	t302 = sin(qJ(6));
	t353 = -pkin(4) - pkin(5);
	t330 = t302 * r_i_i_C(2) + t353;
	t315 = t306 * r_i_i_C(1) - t330;
	t331 = -t302 * r_i_i_C(1) - qJ(5);
	t316 = t306 * r_i_i_C(2) - t331;
	t355 = t316 * t298 + t315 * t299 + t297;
	t308 = cos(qJ(2));
	t337 = -r_i_i_C(3) - pkin(9) + qJ(4) + pkin(8);
	t361 = t337 * t308;
	t366 = t355 * t304 - t361;
	t329 = qJD(3) * t308 - qJD(1);
	t303 = sin(qJ(3));
	t351 = t303 * pkin(3);
	t365 = qJD(1) * pkin(7) + (-t304 * t297 + t361) * qJD(2) + t304 * qJD(4) - t329 * t351;
	t305 = sin(qJ(1));
	t348 = t305 * t308;
	t309 = cos(qJ(1));
	t350 = t299 * t309;
	t281 = t298 * t348 + t350;
	t347 = t309 * t298;
	t282 = t299 * t348 - t347;
	t322 = t281 * t302 + t282 * t306;
	t323 = t281 * t306 - t282 * t302;
	t363 = (t322 * r_i_i_C(1) + t323 * r_i_i_C(2)) * qJD(6);
	t318 = t298 * t302 + t299 * t306;
	t319 = t298 * t306 - t299 * t302;
	t362 = (t315 * t298 - t316 * t299 + t351) * qJD(3) - (t319 * r_i_i_C(1) - t318 * r_i_i_C(2)) * qJD(6) - t337 * qJD(2) - t298 * qJD(5);
	t358 = (-qJD(3) + qJD(6)) * t304;
	t340 = qJD(3) * t305;
	t344 = qJD(1) * t309;
	t314 = t298 * t344 + t299 * t340;
	t339 = qJD(3) * t309;
	t333 = t298 * t339;
	t343 = qJD(2) * t304;
	t336 = t305 * t343;
	t346 = qJD(1) * t305;
	t279 = -t298 * t336 - t299 * t346 + t314 * t308 - t333;
	t276 = t279 * t306;
	t345 = qJD(1) * t308;
	t342 = qJD(2) * t308;
	t341 = qJD(2) * t309;
	t335 = t304 * t341;
	t334 = t298 * t340;
	t332 = t299 * t339;
	t328 = -qJD(3) + t345;
	t325 = -r_i_i_C(1) * (t318 * t358 - t319 * t342) - r_i_i_C(2) * (t318 * t342 + t319 * t358);
	t283 = -t305 * t299 + t308 * t347;
	t284 = t305 * t298 + t308 * t350;
	t321 = t283 * t306 - t284 * t302;
	t320 = t283 * t302 + t284 * t306;
	t317 = t329 * t307;
	t313 = qJD(3) * t352 + (-t297 * t308 - t337 * t304 - pkin(1)) * qJD(1);
	t311 = -qJD(2) * t355 + qJD(4);
	t310 = t362 * t304 + t311 * t308;
	t280 = t284 * qJD(1) - t299 * t336 - t308 * t334 - t332;
	t278 = t308 * t333 + (t305 * t345 + t335) * t299 - t314;
	t277 = t281 * qJD(1) + t298 * t335 - t308 * t332 - t334;
	t273 = t321 * qJD(6) - t277 * t302 - t278 * t306;
	t272 = -t320 * qJD(6) - t277 * t306 + t278 * t302;
	t1 = [-t276 * r_i_i_C(2) - t281 * qJD(5) + t331 * t279 - t315 * t280 + (-t323 * r_i_i_C(1) + t322 * r_i_i_C(2)) * qJD(6) + t313 * t309 - t365 * t305, t310 * t309 + t366 * t346, t284 * qJD(5) - t316 * t278 + t315 * t277 + (t320 * r_i_i_C(1) + t321 * r_i_i_C(2)) * qJD(6) + (-t309 * t317 + (t328 * t305 + t335) * t303) * pkin(3), -t304 * t346 + t308 * t341, -t277, r_i_i_C(1) * t272 - r_i_i_C(2) * t273; t273 * r_i_i_C(1) + t272 * r_i_i_C(2) - t277 * qJ(5) + t283 * qJD(5) + t353 * t278 + t313 * t305 + t365 * t309, t310 * t305 - t366 * t344, -t276 * r_i_i_C(1) + t282 * qJD(5) + t316 * t280 + t330 * t279 + t363 + (-t305 * t317 + (-t328 * t309 + t336) * t303) * pkin(3), t304 * t344 + t305 * t342, t279, (-t280 * t302 + t276) * r_i_i_C(1) + (-t279 * t302 - t280 * t306) * r_i_i_C(2) - t363; 0, t311 * t304 - t362 * t308, (t299 * qJ(5) + t353 * t298 - t351) * t342 + (qJD(5) * t299 + (-t298 * qJ(5) + t353 * t299 - t352) * qJD(3)) * t304 - t325, t343, t304 * qJD(3) * t299 + t298 * t342, t325;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end