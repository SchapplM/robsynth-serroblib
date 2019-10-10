% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
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
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
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
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
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
	% StartTime: 2019-10-10 12:24:01
	% EndTime: 2019-10-10 12:24:02
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
	% StartTime: 2019-10-10 12:24:02
	% EndTime: 2019-10-10 12:24:02
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (474->75), mult. (702->108), div. (0->0), fcn. (565->8), ass. (0->61)
	t315 = sin(qJ(4));
	t366 = r_i_i_C(3) + qJ(5);
	t341 = t366 * t315;
	t377 = pkin(3) + t341;
	t314 = qJ(2) + qJ(3);
	t312 = cos(t314);
	t316 = sin(qJ(2));
	t365 = pkin(2) * qJD(2);
	t351 = t316 * t365;
	t311 = sin(t314);
	t313 = qJD(2) + qJD(3);
	t364 = t311 * t313;
	t355 = qJD(5) * t315;
	t370 = pkin(9) + r_i_i_C(1);
	t372 = t370 * t313 + t355;
	t376 = -pkin(3) * t364 + t372 * t312 - t351;
	t318 = cos(qJ(4));
	t367 = r_i_i_C(2) * t318;
	t375 = t311 * t367 + t370 * t312;
	t371 = pkin(4) - r_i_i_C(2);
	t352 = t313 * t367;
	t358 = qJD(4) * t315;
	t368 = pkin(4) * t311;
	t374 = t312 * t352 + t358 * t368;
	t373 = t366 * t318;
	t369 = pkin(2) * t316;
	t363 = t312 * t313;
	t317 = sin(qJ(1));
	t362 = t317 * t315;
	t320 = cos(qJ(1));
	t361 = t320 * t318;
	t360 = qJD(1) * t317;
	t359 = qJD(1) * t320;
	t357 = qJD(4) * t318;
	t356 = qJD(4) * t320;
	t354 = t318 * qJD(5);
	t350 = t317 * t364;
	t349 = t320 * t364;
	t347 = t318 * t360;
	t343 = t317 * t358;
	t342 = t318 * t356;
	t337 = t374 * t317 + t375 * t359;
	t336 = t377 * t311 * t360 + t374 * t320 + t347 * t368;
	t335 = t312 * t361 + t362;
	t319 = cos(qJ(2));
	t331 = -t319 * pkin(2) - pkin(3) * t312 - t370 * t311 - pkin(1);
	t330 = -pkin(4) * t318 - t377;
	t329 = t315 * t356 + t347;
	t328 = t315 * t359 + t317 * t357;
	t327 = t330 * t313;
	t326 = t312 * t327;
	t325 = (-r_i_i_C(2) * t315 - t373) * qJD(4) - t372;
	t324 = t370 * t363 + (t327 + t352) * t311 + (t366 * t357 - t371 * t358 + t355) * t312;
	t323 = t311 * t325 + t326;
	t322 = -t319 * t365 + t323;
	t321 = -pkin(8) - pkin(7);
	t281 = t335 * qJD(1) - t312 * t343 - t318 * t350 - t342;
	t280 = t312 * t328 - t315 * t350 - t329;
	t279 = t312 * t329 + t318 * t349 - t328;
	t278 = t315 * t349 - t312 * t342 - t343 + (t312 * t362 + t361) * qJD(1);
	t1 = [-t320 * t354 - t371 * t281 - t366 * t280 - t376 * t317 + (t317 * t321 + t320 * t331) * qJD(1), (-t375 + t369) * t360 + t322 * t320 + t336, t320 * t323 - t360 * t375 + t336, t335 * qJD(5) + t371 * t278 - t366 * t279, -t278, 0; -t317 * t354 - t371 * t279 - t366 * t278 + t376 * t320 + (t317 * t331 - t320 * t321) * qJD(1), (t311 * t330 - t369) * t359 + t322 * t317 + t337, t317 * t326 + (t317 * t325 + t330 * t359) * t311 + t337, -(-t317 * t312 * t318 + t320 * t315) * qJD(5) + t366 * t281 - t371 * t280, t280, 0; 0, t324 - t351, t324, (-t371 * t315 + t373) * t363 + (t354 + (-t371 * t318 - t341) * qJD(4)) * t311, t311 * t357 + t315 * t363, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:02
	% EndTime: 2019-10-10 12:24:02
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (635->81), mult. (928->113), div. (0->0), fcn. (759->8), ass. (0->61)
	t321 = sin(qJ(4));
	t374 = r_i_i_C(2) + qJ(5);
	t380 = t374 * t321;
	t324 = cos(qJ(4));
	t326 = cos(qJ(1));
	t361 = qJD(4) * t326;
	t323 = sin(qJ(1));
	t365 = qJD(1) * t323;
	t338 = t321 * t361 + t324 * t365;
	t320 = qJ(2) + qJ(3);
	t317 = sin(t320);
	t319 = qJD(2) + qJD(3);
	t318 = cos(t320);
	t359 = pkin(5) + pkin(9) + r_i_i_C(1);
	t346 = t359 * t318;
	t322 = sin(qJ(2));
	t373 = pkin(2) * qJD(2);
	t357 = t322 * t373;
	t379 = (-pkin(3) * t317 + t346) * t319 - t357;
	t358 = pkin(4) + r_i_i_C(3) + qJ(6);
	t377 = t317 * t358;
	t376 = pkin(2) * t322;
	t372 = t317 * t319;
	t371 = t318 * t319;
	t370 = t318 * t324;
	t369 = t323 * t321;
	t368 = t323 * t324;
	t367 = t326 * t321;
	t366 = t326 * t324;
	t364 = qJD(1) * t326;
	t363 = qJD(4) * t321;
	t362 = qJD(4) * t324;
	t360 = qJD(5) * t321;
	t356 = t323 * t372;
	t355 = t326 * t372;
	t350 = t323 * t363;
	t348 = t324 * t361;
	t347 = t359 * t317;
	t345 = t358 * t321;
	t340 = t364 * t346 + t350 * t377;
	t339 = t318 * t366 + t369;
	t280 = t318 * t369 + t366;
	t337 = t321 * t364 + t323 * t362;
	t325 = cos(qJ(2));
	t336 = -t325 * pkin(2) - pkin(3) * t318 - pkin(1) - t347;
	t335 = -t358 * t324 - t380;
	t334 = t338 * t377 + (pkin(3) + t380) * t317 * t365;
	t333 = -pkin(3) + t335;
	t332 = -t360 + (-t374 * qJD(4) - qJD(6)) * t324;
	t331 = t333 * t319;
	t330 = qJD(6) * t370 + t317 * t331 + t359 * t371 + (-qJD(4) * t345 + t374 * t362 + t360) * t318;
	t329 = t332 * t317 + (t333 * t318 - t347) * t319;
	t328 = -t325 * t373 + t329;
	t327 = -pkin(8) - pkin(7);
	t282 = t318 * t367 - t368;
	t281 = t318 * t368 - t367;
	t279 = t339 * qJD(1) - t318 * t350 - t324 * t356 - t348;
	t278 = t337 * t318 - t321 * t356 - t338;
	t277 = t338 * t318 + t324 * t355 - t337;
	t276 = t280 * qJD(1) - t318 * t348 + t321 * t355 - t350;
	t1 = [-t280 * qJD(5) - t281 * qJD(6) - t374 * t278 - t358 * t279 - t379 * t323 + (t323 * t327 + t336 * t326) * qJD(1), (-t346 + t376) * t365 + t328 * t326 + t334, t329 * t326 - t346 * t365 + t334, qJD(5) * t339 - t282 * qJD(6) + t358 * t276 - t374 * t277, -t276, -t277; t282 * qJD(5) + t339 * qJD(6) - t374 * t276 - t358 * t277 + t379 * t326 + (t336 * t323 - t326 * t327) * qJD(1), (t333 * t317 - t376) * t364 + t328 * t323 + t340, t323 * t318 * t331 + (t333 * t364 + (-t359 * t319 + t332) * t323) * t317 + t340, t281 * qJD(5) - t280 * qJD(6) - t358 * t278 + t374 * t279, t278, t279; 0, t330 - t357, t330, (t374 * t324 - t345) * t371 + (t335 * qJD(4) + qJD(5) * t324 - qJD(6) * t321) * t317, t317 * t362 + t321 * t371, -t317 * t363 + t319 * t370;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end