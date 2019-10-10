% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRP3
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
% Datum: 2019-10-10 13:00
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:00:21
	% EndTime: 2019-10-10 13:00:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:00:21
	% EndTime: 2019-10-10 13:00:21
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
	% StartTime: 2019-10-10 13:00:21
	% EndTime: 2019-10-10 13:00:21
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
	% StartTime: 2019-10-10 13:00:21
	% EndTime: 2019-10-10 13:00:21
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (81->26), mult. (114->37), div. (0->0), fcn. (73->6), ass. (0->27)
	t38 = qJD(2) + qJD(3);
	t39 = qJ(2) + qJ(3);
	t37 = cos(t39);
	t59 = r_i_i_C(2) * t37;
	t36 = sin(t39);
	t61 = r_i_i_C(1) * t36;
	t49 = t59 + t61;
	t47 = t49 * t38;
	t40 = sin(qJ(2));
	t62 = pkin(2) * t40;
	t63 = qJD(2) * t62 + t47;
	t60 = r_i_i_C(2) * t36;
	t58 = r_i_i_C(3) + pkin(8) + pkin(7);
	t57 = t37 * t38;
	t41 = sin(qJ(1));
	t56 = qJD(1) * t41;
	t43 = cos(qJ(1));
	t55 = qJD(1) * t43;
	t42 = cos(qJ(2));
	t54 = qJD(2) * t42;
	t53 = r_i_i_C(1) * t57;
	t52 = t38 * t60;
	t51 = qJD(1) * t59;
	t48 = -t42 * pkin(2) - r_i_i_C(1) * t37 - pkin(1) + t60;
	t46 = t41 * t51 + t56 * t61 + (t52 - t53) * t43;
	t31 = t41 * t52;
	t1 = [t63 * t41 + (-t58 * t41 + t48 * t43) * qJD(1), (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0, 0, 0; -t63 * t43 + (t48 * t41 + t58 * t43) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0, 0, 0; 0, -t63, -t47, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:00:22
	% EndTime: 2019-10-10 13:00:23
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (270->58), mult. (394->95), div. (0->0), fcn. (297->8), ass. (0->54)
	t258 = qJ(2) + qJ(3);
	t255 = sin(t258);
	t262 = cos(qJ(4));
	t310 = r_i_i_C(1) * t262 + pkin(3);
	t313 = t255 * t310;
	t259 = sin(qJ(4));
	t293 = qJD(4) * t262;
	t256 = cos(t258);
	t257 = qJD(2) + qJD(3);
	t301 = t256 * t257;
	t312 = t255 * t293 + t259 * t301;
	t307 = pkin(9) + r_i_i_C(3);
	t288 = t307 * t256;
	t260 = sin(qJ(2));
	t302 = pkin(2) * qJD(2);
	t290 = t260 * t302;
	t305 = pkin(3) * t255;
	t311 = -t290 + (t288 - t305) * t257;
	t294 = qJD(4) * t259;
	t281 = t255 * t294;
	t308 = r_i_i_C(1) * t281 + t312 * r_i_i_C(2);
	t306 = pkin(2) * t260;
	t303 = r_i_i_C(2) * t259;
	t261 = sin(qJ(1));
	t300 = t257 * t261;
	t299 = t257 * t262;
	t264 = cos(qJ(1));
	t298 = t257 * t264;
	t297 = t262 * t264;
	t296 = qJD(1) * t261;
	t295 = qJD(1) * t264;
	t292 = t255 * t303;
	t291 = qJD(1) * t303;
	t289 = t307 * t255;
	t287 = t307 * t261;
	t286 = t255 * t299;
	t276 = qJD(4) * t256 - qJD(1);
	t275 = qJD(1) * t256 - qJD(4);
	t274 = t310 * t256;
	t273 = t310 * t264;
	t272 = t308 * t264 + t296 * t313;
	t271 = t276 * t259;
	t270 = t264 * t255 * t291 + t308 * t261 + t295 * t288;
	t263 = cos(qJ(2));
	t269 = -t263 * pkin(2) - pkin(3) * t256 - pkin(1) - t289;
	t268 = t255 * t298 + t275 * t261;
	t267 = -t263 * t302 + (-t274 - t289) * t257;
	t266 = -t256 * r_i_i_C(2) * t293 + (-t256 * t294 - t286) * r_i_i_C(1) + t307 * t301 + (-t305 + t292) * t257;
	t265 = -pkin(8) - pkin(7);
	t239 = -t275 * t297 + (t271 + t286) * t261;
	t238 = t276 * t262 * t261 + (-t255 * t300 + t275 * t264) * t259;
	t237 = t268 * t262 + t264 * t271;
	t236 = t268 * t259 - t276 * t297;
	t1 = [t239 * r_i_i_C(1) + t238 * r_i_i_C(2) - t311 * t261 + (t261 * t265 + t269 * t264) * qJD(1), (-t288 - t292 + t306) * t296 + t267 * t264 + t272, (-t261 * t291 - t307 * t298) * t255 + (-qJD(1) * t287 - t257 * t273) * t256 + t272, t236 * r_i_i_C(1) + t237 * r_i_i_C(2), 0, 0; -t237 * r_i_i_C(1) + t236 * r_i_i_C(2) + t311 * t264 + (t269 * t261 - t264 * t265) * qJD(1), (-t306 - t313) * t295 + t267 * t261 + t270, -t274 * t300 + (-qJD(1) * t273 - t257 * t287) * t255 + t270, -t238 * r_i_i_C(1) + t239 * r_i_i_C(2), 0, 0; 0, t266 - t290, t266, (-t256 * t299 + t281) * r_i_i_C(2) - t312 * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:00:22
	% EndTime: 2019-10-10 13:00:23
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (524->77), mult. (560->110), div. (0->0), fcn. (431->10), ass. (0->69)
	t291 = qJ(4) + qJ(5);
	t285 = sin(t291);
	t292 = qJ(2) + qJ(3);
	t286 = sin(t292);
	t288 = cos(t292);
	t290 = qJD(2) + qJD(3);
	t343 = t288 * t290;
	t287 = cos(t291);
	t289 = qJD(4) + qJD(5);
	t344 = t287 * t289;
	t364 = t285 * t343 + t286 * t344;
	t299 = -pkin(10) - pkin(9);
	t293 = sin(qJ(4));
	t347 = pkin(4) * qJD(4);
	t333 = t293 * t347;
	t363 = t290 * t299 + t333;
	t296 = cos(qJ(4));
	t283 = t296 * pkin(4) + pkin(3);
	t361 = r_i_i_C(1) * t287 + t283;
	t362 = t286 * t361 + t288 * t299;
	t294 = sin(qJ(2));
	t348 = pkin(2) * qJD(2);
	t330 = t294 * t348;
	t345 = t286 * t290;
	t349 = r_i_i_C(3) - t299;
	t353 = pkin(4) * t293;
	t360 = (t349 * t290 - t333) * t288 + (pkin(8) + pkin(7) + t353) * qJD(1) - t283 * t345 - t330;
	t346 = t285 * t286;
	t329 = t289 * t346;
	t359 = r_i_i_C(1) * t329 + t364 * r_i_i_C(2) + t363 * t286;
	t298 = cos(qJ(1));
	t317 = t288 * t289 - qJD(1);
	t358 = t298 * t317;
	t337 = qJD(1) * t288;
	t316 = -t289 + t337;
	t295 = sin(qJ(1));
	t326 = t295 * t345;
	t356 = t316 * t298 - t326;
	t354 = pkin(2) * t294;
	t351 = r_i_i_C(2) * t287;
	t350 = r_i_i_C(3) * t288;
	t341 = t290 * t298;
	t325 = t286 * t341;
	t304 = t316 * t295 + t325;
	t261 = t304 * t285 - t287 * t358;
	t262 = t285 * t358 + t304 * t287;
	t339 = t261 * r_i_i_C(1) + t262 * r_i_i_C(2);
	t311 = t317 * t295;
	t263 = t356 * t285 + t287 * t311;
	t264 = t285 * t311 - t356 * t287;
	t338 = -t263 * r_i_i_C(1) + t264 * r_i_i_C(2);
	t336 = qJD(1) * t295;
	t335 = qJD(1) * t298;
	t334 = r_i_i_C(2) * t346;
	t332 = t296 * t347;
	t331 = r_i_i_C(2) * qJD(1) * t285;
	t314 = -qJD(4) + t337;
	t313 = -r_i_i_C(1) * t285 - t351;
	t312 = t361 * t290;
	t310 = (-qJD(4) * t288 + qJD(1)) * t296;
	t309 = t298 * t286 * t331 + t359 * t295 + t335 * t350;
	t306 = t359 * t298 + t362 * t336;
	t305 = (-r_i_i_C(3) * t286 - t288 * t361) * t290;
	t297 = cos(qJ(2));
	t303 = t332 + (-t297 * pkin(2) - t283 * t288 - t349 * t286 - pkin(1)) * qJD(1);
	t302 = -t297 * t348 + t305;
	t301 = t290 * t334 + r_i_i_C(3) * t343 - t286 * t312 + (t313 * t289 - t363) * t288;
	t274 = r_i_i_C(2) * t329;
	t1 = [t264 * r_i_i_C(1) + t263 * r_i_i_C(2) - t360 * t295 + t303 * t298, (-t334 - t350 + t354) * t336 + t302 * t298 + t306, (-r_i_i_C(3) * t341 - t295 * t331) * t286 + (-r_i_i_C(3) * t336 - t298 * t312) * t288 + t306, (t298 * t310 + (t314 * t295 + t325) * t293) * pkin(4) + t339, t339, 0; -t262 * r_i_i_C(1) + t261 * r_i_i_C(2) + t303 * t295 + t360 * t298, t302 * t295 + (-t362 - t354) * t335 + t309, t295 * t305 - t335 * t362 + t309, (t295 * t310 + (-t314 * t298 + t326) * t293) * pkin(4) + t338, t338, 0; 0, t301 - t330, t301, t274 + (-r_i_i_C(1) * t344 - t332) * t286 + (t313 - t353) * t343, -t364 * r_i_i_C(1) - t343 * t351 + t274, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:00:22
	% EndTime: 2019-10-10 13:00:23
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (695->94), mult. (678->129), div. (0->0), fcn. (522->10), ass. (0->75)
	t302 = qJ(4) + qJ(5);
	t294 = sin(t302);
	t309 = cos(qJ(1));
	t300 = qJD(4) + qJD(5);
	t303 = qJ(2) + qJ(3);
	t297 = cos(t303);
	t344 = qJD(1) * t297;
	t327 = -t300 + t344;
	t295 = sin(t303);
	t301 = qJD(2) + qJD(3);
	t306 = sin(qJ(1));
	t349 = t301 * t306;
	t336 = t295 * t349;
	t365 = t327 * t309 - t336;
	t370 = t365 * t294;
	t296 = cos(t302);
	t307 = cos(qJ(4));
	t286 = t307 * pkin(4) + pkin(5) * t296;
	t284 = pkin(3) + t286;
	t369 = r_i_i_C(1) * t296 + t284;
	t304 = sin(qJ(4));
	t355 = pkin(4) * qJD(4);
	t361 = pkin(5) * t294;
	t275 = -t300 * t361 - t304 * t355;
	t285 = t304 * pkin(4) + t361;
	t305 = sin(qJ(2));
	t356 = pkin(2) * qJD(2);
	t339 = t305 * t356;
	t299 = -qJ(6) - pkin(10) - pkin(9);
	t357 = r_i_i_C(3) - t299;
	t368 = (t357 * t301 + t275) * t297 + (t285 + pkin(8) + pkin(7)) * qJD(1) - (t284 * t301 - qJD(6)) * t295 - t339;
	t351 = t297 * t301;
	t337 = t294 * t351;
	t354 = t295 * t300;
	t338 = t294 * t354;
	t353 = t296 * t300;
	t359 = r_i_i_C(2) * t295;
	t367 = r_i_i_C(1) * t338 + r_i_i_C(2) * t337 + qJD(6) * t297 + t353 * t359;
	t363 = -pkin(5) - r_i_i_C(1);
	t362 = pkin(2) * t305;
	t358 = r_i_i_C(3) * t297;
	t352 = t297 * t299;
	t350 = t301 * t295;
	t348 = t301 * t309;
	t335 = t295 * t348;
	t315 = t327 * t306 + t335;
	t328 = t297 * t300 - qJD(1);
	t322 = t296 * t328;
	t263 = t315 * t294 - t309 * t322;
	t264 = t328 * t309 * t294 + t315 * t296;
	t347 = t263 * r_i_i_C(1) + t264 * r_i_i_C(2);
	t321 = t328 * t306;
	t265 = t296 * t321 + t370;
	t266 = t294 * t321 - t296 * t365;
	t346 = -t265 * r_i_i_C(1) + t266 * r_i_i_C(2);
	t343 = qJD(1) * t306;
	t342 = qJD(1) * t309;
	t340 = t294 * t359;
	t334 = t295 * t343;
	t333 = t295 * t342;
	t325 = -r_i_i_C(1) * t294 - r_i_i_C(2) * t296;
	t324 = t285 * t344 + t275;
	t323 = t369 * t301;
	t320 = -t340 - t358;
	t319 = t294 * r_i_i_C(2) * t333 + t299 * t336 + t367 * t306 + t342 * t358;
	t318 = t299 * t335 + t367 * t309 + t369 * t334 + t343 * t352;
	t276 = pkin(5) * t353 + t307 * t355;
	t317 = qJD(1) * t286 - t276 * t297 + t285 * t350;
	t308 = cos(qJ(2));
	t314 = t276 + (-t308 * pkin(2) - t284 * t297 - t357 * t295 - pkin(1)) * qJD(1);
	t313 = -t275 * t295 + (-r_i_i_C(3) * t295 - t297 * t369) * t301;
	t312 = -t308 * t356 + t313;
	t311 = r_i_i_C(3) * t351 + t301 * t340 + (-t299 * t301 + t325 * t300 + t275) * t297 + (qJD(6) - t323) * t295;
	t281 = r_i_i_C(2) * t338;
	t1 = [t266 * r_i_i_C(1) + t265 * r_i_i_C(2) - t368 * t306 + t314 * t309, (t320 + t362) * t343 + t312 * t309 + t318, t313 * t309 + t320 * t343 + t318, t324 * t306 + t317 * t309 + t347, t263 * pkin(5) + t347, t297 * t348 - t334; -t264 * r_i_i_C(1) + t263 * r_i_i_C(2) + t314 * t306 + t368 * t309, (-t295 * t369 - t352 - t362) * t342 + t312 * t306 + t319, (-t299 * t342 - t306 * t323) * t297 + ((-r_i_i_C(3) * t301 - t275) * t306 - t369 * t342) * t295 + t319, t317 * t306 - t324 * t309 + t346, (-t306 * t322 - t370) * pkin(5) + t346, t297 * t349 + t333; 0, t311 - t339, t311, t281 + (-r_i_i_C(1) * t353 - t276) * t295 + (-t285 + t325) * t351, t281 + t363 * t337 + (-r_i_i_C(2) * t351 + t363 * t354) * t296, t350;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end