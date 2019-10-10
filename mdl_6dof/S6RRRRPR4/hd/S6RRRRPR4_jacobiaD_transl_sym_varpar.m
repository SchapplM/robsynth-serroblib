% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:36
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
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
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
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
	t1 = [t23 * t25 + (-t18 * t28 + t20 * t22) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0, 0; -t20 * t21 + (t18 * t22 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0, 0; 0, -t21, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
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
	% StartTime: 2019-10-10 12:36:52
	% EndTime: 2019-10-10 12:36:52
	% DurationCPUTime: 0.36s
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
	% StartTime: 2019-10-10 12:36:52
	% EndTime: 2019-10-10 12:36:53
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (443->75), mult. (515->106), div. (0->0), fcn. (394->10), ass. (0->62)
	t269 = qJ(4) + pkin(11);
	t264 = sin(t269);
	t265 = cos(t269);
	t334 = r_i_i_C(1) * t264 + r_i_i_C(2) * t265;
	t268 = qJD(2) + qJD(3);
	t322 = r_i_i_C(2) * t264;
	t333 = t268 * t322 + qJD(5);
	t275 = cos(qJ(4));
	t319 = t275 * pkin(4);
	t262 = pkin(3) + t319;
	t323 = r_i_i_C(1) * t265;
	t332 = t262 + t323;
	t270 = qJ(2) + qJ(3);
	t266 = sin(t270);
	t267 = cos(t270);
	t273 = sin(qJ(2));
	t317 = pkin(2) * qJD(2);
	t305 = t273 * t317;
	t272 = sin(qJ(4));
	t316 = pkin(4) * qJD(4);
	t306 = t272 * t316;
	t271 = -qJ(5) - pkin(9);
	t318 = r_i_i_C(3) - t271;
	t325 = pkin(4) * t272;
	t331 = (t318 * t268 - t306) * t267 + (pkin(8) + pkin(7) + t325) * qJD(1) - (t262 * t268 - qJD(5)) * t266 - t305;
	t309 = qJD(4) * t266;
	t330 = t266 * t306 + t333 * t267 + t334 * t309;
	t274 = sin(qJ(1));
	t293 = qJD(4) * t267 - qJD(1);
	t289 = t293 * t274;
	t277 = cos(qJ(1));
	t290 = t293 * t277;
	t292 = qJD(1) * t267 - qJD(4);
	t313 = t268 * t274;
	t304 = t266 * t313;
	t328 = t292 * t277 - t304;
	t326 = pkin(2) * t273;
	t320 = r_i_i_C(3) * t267;
	t315 = t267 * t268;
	t314 = t267 * t271;
	t312 = t268 * t277;
	t311 = qJD(1) * t274;
	t310 = qJD(1) * t277;
	t303 = t266 * t312;
	t302 = t266 * t311;
	t301 = t266 * t310;
	t291 = t332 * t268;
	t288 = -t325 - t334;
	t287 = -t266 * t332 - t314;
	t286 = t271 * t304 + t330 * t274 + t301 * t322 + t310 * t320;
	t285 = t271 * t303 + t330 * t277 + t332 * t302 + t311 * t314;
	t284 = (-r_i_i_C(3) * t266 - t267 * t332) * t268;
	t282 = t292 * t274 + t303;
	t276 = cos(qJ(2));
	t281 = t275 * t316 + (-t276 * pkin(2) - t262 * t267 - t318 * t266 - pkin(1)) * qJD(1);
	t280 = -t276 * t317 + t284;
	t279 = r_i_i_C(3) * t315 + (t288 * qJD(4) - t268 * t271) * t267 + (-t291 + t333) * t266;
	t241 = t264 * t289 - t265 * t328;
	t240 = t328 * t264 + t265 * t289;
	t239 = t264 * t290 + t282 * t265;
	t238 = t282 * t264 - t265 * t290;
	t1 = [t241 * r_i_i_C(1) + t240 * r_i_i_C(2) - t331 * t274 + t281 * t277, (-t266 * t322 - t320 + t326) * t311 + t280 * t277 + t285, (-r_i_i_C(3) * t312 - t311 * t322) * t266 + (-r_i_i_C(3) * t311 - t277 * t291) * t267 + t285, t238 * r_i_i_C(1) + t239 * r_i_i_C(2) + (t282 * t272 - t275 * t290) * pkin(4), t267 * t312 - t302, 0; -t239 * r_i_i_C(1) + t238 * r_i_i_C(2) + t281 * t274 + t331 * t277, t280 * t274 + (t287 - t326) * t310 + t286, t274 * t284 + t287 * t310 + t286, -t240 * r_i_i_C(1) + t241 * r_i_i_C(2) + (-t272 * t328 - t275 * t289) * pkin(4), t267 * t313 + t301, 0; 0, t279 - t305, t279, t288 * t315 + (-t319 + t322 - t323) * t309, t268 * t266, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:52
	% EndTime: 2019-10-10 12:36:53
	% DurationCPUTime: 0.64s
	% Computational Cost: add. (740->88), mult. (636->118), div. (0->0), fcn. (492->12), ass. (0->71)
	t304 = qJ(4) + pkin(11);
	t297 = qJ(6) + t304;
	t292 = sin(t297);
	t293 = cos(t297);
	t305 = qJ(2) + qJ(3);
	t299 = cos(t305);
	t303 = qJD(2) + qJD(3);
	t352 = t299 * t303;
	t298 = sin(t305);
	t302 = qJD(4) + qJD(6);
	t354 = t298 * t302;
	t369 = t292 * t352 + t293 * t354;
	t285 = pkin(5) * cos(t304) + cos(qJ(4)) * pkin(4);
	t283 = pkin(3) + t285;
	t360 = r_i_i_C(1) * t293;
	t368 = t283 + t360;
	t284 = sin(qJ(4)) * pkin(4) + pkin(5) * sin(t304);
	t280 = t284 * qJD(4);
	t307 = sin(qJ(2));
	t355 = pkin(2) * qJD(2);
	t340 = t307 * t355;
	t301 = -pkin(10) - qJ(5) - pkin(9);
	t356 = r_i_i_C(3) - t301;
	t367 = (t356 * t303 - t280) * t299 + (t284 + pkin(8) + pkin(7)) * qJD(1) - (t283 * t303 - qJD(5)) * t298 - t340;
	t339 = t292 * t354;
	t366 = r_i_i_C(1) * t339 + t369 * r_i_i_C(2) + qJD(5) * t299;
	t311 = cos(qJ(1));
	t328 = t299 * t302 - qJD(1);
	t365 = t311 * t328;
	t345 = qJD(1) * t299;
	t327 = -t302 + t345;
	t308 = sin(qJ(1));
	t350 = t303 * t308;
	t336 = t298 * t350;
	t363 = t327 * t311 - t336;
	t361 = pkin(2) * t307;
	t359 = r_i_i_C(2) * t292;
	t358 = r_i_i_C(2) * t293;
	t357 = r_i_i_C(3) * t299;
	t353 = t299 * t301;
	t351 = t303 * t298;
	t349 = t303 * t311;
	t335 = t298 * t349;
	t317 = t327 * t308 + t335;
	t262 = t317 * t292 - t293 * t365;
	t263 = t292 * t365 + t317 * t293;
	t348 = t262 * r_i_i_C(1) + t263 * r_i_i_C(2);
	t322 = t328 * t308;
	t264 = t363 * t292 + t293 * t322;
	t265 = t292 * t322 - t363 * t293;
	t347 = -t264 * r_i_i_C(1) + t265 * r_i_i_C(2);
	t344 = qJD(1) * t308;
	t343 = qJD(1) * t311;
	t341 = t298 * t359;
	t334 = t298 * t344;
	t333 = t298 * t343;
	t325 = -r_i_i_C(1) * t292 - t358;
	t324 = t284 * t345 - t280;
	t323 = t368 * t303;
	t321 = -t341 - t357;
	t320 = t301 * t336 + t366 * t308 + t333 * t359 + t343 * t357;
	t319 = t301 * t335 + t366 * t311 + t368 * t334 + t344 * t353;
	t281 = t285 * qJD(4);
	t318 = qJD(1) * t285 - t281 * t299 + t284 * t351;
	t310 = cos(qJ(2));
	t316 = t281 + (-t310 * pkin(2) - t283 * t299 - t356 * t298 - pkin(1)) * qJD(1);
	t315 = t280 * t298 + (-r_i_i_C(3) * t298 - t299 * t368) * t303;
	t314 = -t310 * t355 + t315;
	t313 = r_i_i_C(3) * t352 + t303 * t341 + (-t301 * t303 + t325 * t302 - t280) * t299 + (qJD(5) - t323) * t298;
	t276 = r_i_i_C(2) * t339;
	t1 = [t265 * r_i_i_C(1) + t264 * r_i_i_C(2) - t367 * t308 + t316 * t311, (t321 + t361) * t344 + t314 * t311 + t319, t315 * t311 + t321 * t344 + t319, t324 * t308 + t318 * t311 + t348, t299 * t349 - t334, t348; -t263 * r_i_i_C(1) + t262 * r_i_i_C(2) + t316 * t308 + t367 * t311, (-t298 * t368 - t353 - t361) * t343 + t314 * t308 + t320, (-t301 * t343 - t308 * t323) * t299 + ((-r_i_i_C(3) * t303 + t280) * t308 - t368 * t343) * t298 + t320, t318 * t308 - t324 * t311 + t347, t299 * t350 + t333, t347; 0, t313 - t340, t313, t276 + (-t302 * t360 - t281) * t298 + (-t284 + t325) * t352, t351, -t369 * r_i_i_C(1) - t352 * t358 + t276;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end