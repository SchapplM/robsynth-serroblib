% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
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
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
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
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
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
	% StartTime: 2019-10-10 12:20:29
	% EndTime: 2019-10-10 12:20:30
	% DurationCPUTime: 0.46s
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
	% StartTime: 2019-10-10 12:20:29
	% EndTime: 2019-10-10 12:20:30
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (443->75), mult. (515->106), div. (0->0), fcn. (394->10), ass. (0->62)
	t269 = qJ(4) + pkin(10);
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
	% StartTime: 2019-10-10 12:20:30
	% EndTime: 2019-10-10 12:20:31
	% DurationCPUTime: 0.62s
	% Computational Cost: add. (750->98), mult. (823->132), div. (0->0), fcn. (662->10), ass. (0->70)
	t388 = pkin(5) + r_i_i_C(1);
	t323 = qJ(4) + pkin(10);
	t318 = sin(t323);
	t382 = r_i_i_C(3) + qJ(6);
	t396 = t382 * t318;
	t319 = cos(t323);
	t326 = sin(qJ(4));
	t386 = pkin(4) * t326;
	t395 = -t388 * t318 + t382 * t319 - t386;
	t324 = qJ(2) + qJ(3);
	t320 = sin(t324);
	t321 = cos(t324);
	t380 = pkin(4) * qJD(4);
	t364 = t326 * t380;
	t393 = qJD(5) * t321 + t320 * t364;
	t328 = sin(qJ(1));
	t374 = qJD(1) * t328;
	t361 = t321 * t374;
	t322 = qJD(2) + qJD(3);
	t331 = cos(qJ(1));
	t377 = t322 * t331;
	t362 = t320 * t377;
	t392 = t361 + t362;
	t329 = cos(qJ(4));
	t384 = t329 * pkin(4);
	t316 = pkin(3) + t384;
	t327 = sin(qJ(2));
	t381 = pkin(2) * qJD(2);
	t365 = t327 * t381;
	t367 = t318 * qJD(6);
	t325 = -qJ(5) - pkin(9);
	t383 = r_i_i_C(2) - t325;
	t391 = (t383 * t322 - t364 + t367) * t321 + (pkin(8) + pkin(7) + t386) * qJD(1) - (t316 * t322 - qJD(5)) * t320 - t365;
	t387 = pkin(2) * t327;
	t385 = r_i_i_C(2) * t321;
	t379 = t321 * t322;
	t378 = t322 * t328;
	t376 = t328 * t318;
	t375 = t331 * t319;
	t373 = qJD(1) * t331;
	t372 = qJD(4) * t319;
	t371 = qJD(4) * t321;
	t370 = qJD(4) * t328;
	t369 = qJD(4) * t331;
	t366 = t319 * qJD(6);
	t363 = t320 * t378;
	t360 = t320 * t374;
	t358 = t318 * t370;
	t357 = t318 * t369;
	t356 = t319 * t369;
	t352 = qJD(1) * t321 - qJD(4);
	t347 = (qJD(1) - t371) * t329;
	t345 = t388 * t320 * t358 + t325 * t363 + t393 * t328 + t373 * t385;
	t343 = t321 * t375 + t376;
	t342 = -t388 * t319 - t396;
	t341 = t318 * t373 + t319 * t370;
	t340 = -t382 * t372 - t367;
	t339 = -t316 + t342;
	t338 = t392 * t325 + t393 * t331 + t388 * (t319 * t360 + t320 * t357) + (t316 + t396) * t360;
	t337 = t339 * t320 - t321 * t325;
	t330 = cos(qJ(2));
	t336 = t329 * t380 - t366 + (-t330 * pkin(2) - t316 * t321 - t383 * t320 - pkin(1)) * qJD(1);
	t335 = t340 * t320 + (-r_i_i_C(2) * t320 + t339 * t321) * t322;
	t334 = -t330 * t381 + t335;
	t333 = r_i_i_C(2) * t379 + t320 * qJD(5) + t321 * t367 + t337 * t322 + t395 * t371;
	t283 = t343 * qJD(1) - t319 * t363 - t321 * t358 - t356;
	t282 = -t318 * t363 - t319 * t374 + t341 * t321 - t357;
	t281 = t392 * t319 + t321 * t357 - t341;
	t280 = t318 * t362 - t321 * t356 - t358 + (t321 * t376 + t375) * qJD(1);
	t1 = [-t382 * t282 - t388 * t283 - t391 * t328 + t336 * t331, t338 + t334 * t331 + (-t385 + t387) * t374, -r_i_i_C(2) * t361 + t335 * t331 + t338, t343 * qJD(6) - t382 * t281 + t388 * t280 + (t331 * t347 + (t352 * t328 + t362) * t326) * pkin(4), t321 * t377 - t360, -t280; -t382 * t280 - t388 * t281 + t336 * t328 + t391 * t331, (t337 - t387) * t373 + t334 * t328 + t345, (-t325 * t373 + t339 * t378) * t321 + ((-r_i_i_C(2) * t322 + t340) * t328 + t339 * t373) * t320 + t345, -(-t328 * t321 * t319 + t331 * t318) * qJD(6) + t382 * t283 - t388 * t282 + (t328 * t347 + (-t352 * t331 + t363) * t326) * pkin(4), t320 * t373 + t321 * t378, t282; 0, t333 - t365, t333, t395 * t379 + (t366 + (t342 - t384) * qJD(4)) * t320, t322 * t320, t318 * t379 + t320 * t372;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end