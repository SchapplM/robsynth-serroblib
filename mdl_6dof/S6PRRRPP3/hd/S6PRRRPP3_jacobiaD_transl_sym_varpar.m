% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(10));
	t50 = sin(pkin(10));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:51
	% EndTime: 2019-10-09 22:44:51
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(10));
	t182 = cos(pkin(10));
	t185 = sin(qJ(2));
	t183 = cos(pkin(6));
	t187 = cos(qJ(2));
	t195 = t183 * t187;
	t201 = -t180 * t185 + t182 * t195;
	t200 = pkin(8) + r_i_i_C(3);
	t181 = sin(pkin(6));
	t184 = sin(qJ(3));
	t198 = t181 * t184;
	t186 = cos(qJ(3));
	t197 = t181 * t186;
	t196 = t183 * t185;
	t193 = t184 * r_i_i_C(1) + t186 * r_i_i_C(2);
	t192 = t186 * r_i_i_C(1) - t184 * r_i_i_C(2) + pkin(2);
	t176 = t180 * t187 + t182 * t196;
	t191 = t180 * t195 + t182 * t185;
	t190 = t180 * t196 - t182 * t187;
	t189 = qJD(3) * t193;
	t188 = qJD(2) * t192;
	t173 = t191 * qJD(2);
	t171 = t201 * qJD(2);
	t1 = [0, -t200 * t173 + t190 * t188 + t191 * t189, t193 * t173 + ((-t180 * t198 + t186 * t190) * r_i_i_C(1) + (-t180 * t197 - t184 * t190) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, t200 * t171 - t176 * t188 - t201 * t189, -t193 * t171 + ((-t176 * t186 + t182 * t198) * r_i_i_C(1) + (t176 * t184 + t182 * t197) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, (-t187 * t189 + (-t192 * t185 + t200 * t187) * qJD(2)) * t181, -t193 * t187 * t181 * qJD(2) + ((-t183 * t184 - t185 * t197) * r_i_i_C(1) + (-t183 * t186 + t185 * t198) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:52
	% EndTime: 2019-10-09 22:44:52
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (224->70), mult. (725->134), div. (0->0), fcn. (724->10), ass. (0->50)
	t291 = sin(qJ(3));
	t294 = cos(qJ(3));
	t290 = sin(qJ(4));
	t293 = cos(qJ(4));
	t305 = r_i_i_C(1) * t293 - r_i_i_C(2) * t290;
	t303 = pkin(3) + t305;
	t324 = pkin(9) + r_i_i_C(3);
	t326 = (t303 * t291 - t324 * t294) * qJD(3);
	t304 = r_i_i_C(1) * t290 + r_i_i_C(2) * t293;
	t321 = cos(pkin(6));
	t288 = sin(pkin(6));
	t320 = t288 * t291;
	t319 = t288 * t294;
	t295 = cos(qJ(2));
	t318 = t288 * t295;
	t292 = sin(qJ(2));
	t317 = qJD(2) * t292;
	t316 = qJD(2) * t295;
	t315 = qJD(4) * t290;
	t314 = qJD(4) * t293;
	t313 = qJD(4) * t294;
	t312 = t288 * t317;
	t311 = t288 * t316;
	t310 = t292 * t321;
	t309 = t295 * t321;
	t287 = sin(pkin(10));
	t307 = t287 * t310;
	t289 = cos(pkin(10));
	t306 = t289 * t309;
	t279 = t287 * t295 + t289 * t310;
	t302 = -t279 * t291 - t289 * t319;
	t301 = -t279 * t294 + t289 * t320;
	t281 = t289 * t295 - t307;
	t300 = -t281 * t291 + t287 * t319;
	t271 = t281 * t294 + t287 * t320;
	t299 = qJD(4) * t304;
	t280 = t287 * t309 + t289 * t292;
	t298 = -t292 * t320 + t321 * t294;
	t283 = t321 * t291 + t292 * t319;
	t297 = -t324 * t291 - t303 * t294 - pkin(2);
	t296 = t304 * t313 + t326;
	t278 = t287 * t292 - t306;
	t277 = -qJD(2) * t307 + t289 * t316;
	t276 = t280 * qJD(2);
	t275 = t279 * qJD(2);
	t274 = -qJD(2) * t306 + t287 * t317;
	t273 = t298 * qJD(3) + t294 * t311;
	t267 = t300 * qJD(3) - t276 * t294;
	t265 = t302 * qJD(3) - t274 * t294;
	t1 = [0, (-t276 * t290 + t281 * t314) * r_i_i_C(1) + (-t276 * t293 - t281 * t315) * r_i_i_C(2) - t276 * pkin(8) + t297 * t277 + t296 * t280, t324 * t267 - t300 * t299 + t303 * (-t271 * qJD(3) + t276 * t291), (-t267 * t290 + t277 * t293) * r_i_i_C(1) + (-t267 * t293 - t277 * t290) * r_i_i_C(2) + ((-t271 * t293 - t280 * t290) * r_i_i_C(1) + (t271 * t290 - t280 * t293) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (-t274 * t290 + t279 * t314) * r_i_i_C(1) + (-t274 * t293 - t279 * t315) * r_i_i_C(2) - t274 * pkin(8) + t297 * t275 + t296 * t278, t324 * t265 - t302 * t299 + t303 * (t301 * qJD(3) + t274 * t291), (-t265 * t290 + t275 * t293) * r_i_i_C(1) + (-t265 * t293 - t275 * t290) * r_i_i_C(2) + ((-t278 * t290 + t293 * t301) * r_i_i_C(1) + (-t278 * t293 - t290 * t301) * r_i_i_C(2)) * qJD(4), 0, 0; 0, ((t297 * qJD(2) + t305 * qJD(4)) * t292 + (qJD(2) * pkin(8) - t326 + t304 * (qJD(2) - t313)) * t295) * t288, t324 * t273 - t298 * t299 + t303 * (-t283 * qJD(3) - t291 * t311), (-t273 * t290 + t293 * t312) * r_i_i_C(1) + (-t273 * t293 - t290 * t312) * r_i_i_C(2) + ((-t283 * t293 + t290 * t318) * r_i_i_C(1) + (t283 * t290 + t293 * t318) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:53
	% EndTime: 2019-10-09 22:44:53
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (437->91), mult. (1367->160), div. (0->0), fcn. (1417->10), ass. (0->64)
	t353 = sin(qJ(3));
	t356 = cos(qJ(3));
	t395 = pkin(9) + r_i_i_C(1);
	t396 = -pkin(3) * t353 + t395 * t356;
	t394 = r_i_i_C(2) - pkin(4);
	t392 = r_i_i_C(3) + qJ(5);
	t391 = cos(pkin(6));
	t350 = sin(pkin(6));
	t390 = t350 * t353;
	t389 = t350 * t356;
	t352 = sin(qJ(4));
	t388 = t352 * t356;
	t357 = cos(qJ(2));
	t387 = t352 * t357;
	t355 = cos(qJ(4));
	t386 = t355 * t357;
	t354 = sin(qJ(2));
	t385 = qJD(2) * t354;
	t384 = qJD(2) * t357;
	t383 = qJD(3) * t353;
	t382 = qJD(3) * t357;
	t381 = qJD(4) * t356;
	t380 = t350 * t385;
	t379 = t350 * t384;
	t378 = t354 * t391;
	t377 = t357 * t391;
	t376 = -qJD(2) + t381;
	t349 = sin(pkin(10));
	t375 = t349 * t378;
	t351 = cos(pkin(10));
	t374 = t351 * t377;
	t336 = -qJD(2) * t374 + t349 * t385;
	t340 = t349 * t354 - t374;
	t373 = t340 * t381 - t336;
	t342 = t349 * t377 + t351 * t354;
	t338 = t342 * qJD(2);
	t372 = t342 * t381 - t338;
	t341 = t349 * t357 + t351 * t378;
	t368 = -t341 * t356 + t351 * t390;
	t371 = t340 * t352 - t355 * t368;
	t343 = t351 * t357 - t375;
	t333 = t343 * t356 + t349 * t390;
	t370 = t333 * t355 + t342 * t352;
	t369 = -t341 * t353 - t351 * t389;
	t367 = -t343 * t353 + t349 * t389;
	t345 = t391 * t353 + t354 * t389;
	t366 = -t345 * t355 + t350 * t387;
	t365 = -t356 * pkin(3) - t395 * t353 - pkin(2);
	t364 = qJD(3) * t396;
	t363 = -t354 * t390 + t391 * t356;
	t362 = t392 * t352 - t394 * t355 + pkin(3);
	t337 = t341 * qJD(2);
	t361 = qJD(4) * t341 - t337 * t356 + t340 * t383;
	t339 = -qJD(2) * t375 + t351 * t384;
	t360 = qJD(4) * t343 - t339 * t356 + t342 * t383;
	t359 = t353 * t382 + (qJD(2) * t356 - qJD(4)) * t354;
	t358 = qJD(5) * t352 + (t394 * t352 + t392 * t355) * qJD(4);
	t335 = t363 * qJD(3) + t356 * t379;
	t329 = t367 * qJD(3) - t338 * t356;
	t327 = t369 * qJD(3) - t336 * t356;
	t322 = -t366 * qJD(4) + t335 * t352 - t355 * t380;
	t316 = t370 * qJD(4) + t329 * t352 - t339 * t355;
	t314 = t371 * qJD(4) + t327 * t352 - t337 * t355;
	t1 = [0, -(t342 * t388 + t343 * t355) * qJD(5) - t338 * pkin(8) - t394 * (t372 * t352 + t360 * t355) + t392 * (t360 * t352 - t372 * t355) - t342 * t364 + t365 * t339, t395 * t329 + t358 * t367 + t362 * (-t333 * qJD(3) + t338 * t353), t370 * qJD(5) + t392 * (t329 * t355 + t339 * t352 + (-t333 * t352 + t342 * t355) * qJD(4)) + t394 * t316, t316, 0; 0, -(t340 * t388 + t341 * t355) * qJD(5) - t336 * pkin(8) - t394 * (t373 * t352 + t361 * t355) + t392 * (t361 * t352 - t373 * t355) - t340 * t364 + t365 * t337, t395 * t327 + t358 * t369 + t362 * (t368 * qJD(3) + t336 * t353), t371 * qJD(5) + t392 * (t327 * t355 + t337 * t352 + (t340 * t355 + t352 * t368) * qJD(4)) + t394 * t314, t314, 0; 0, (t394 * (t359 * t355 + t376 * t387) - t392 * (t359 * t352 - t376 * t386) - (t354 * t355 - t356 * t387) * qJD(5) + t396 * t382 + (t357 * pkin(8) + t365 * t354) * qJD(2)) * t350, t395 * t335 + t358 * t363 + t362 * (-t345 * qJD(3) - t353 * t379), -t366 * qJD(5) + t392 * (t352 * t380 + t335 * t355 + (-t345 * t352 - t350 * t386) * qJD(4)) + t394 * t322, t322, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:53
	% EndTime: 2019-10-09 22:44:54
	% DurationCPUTime: 0.66s
	% Computational Cost: add. (596->104), mult. (1844->182), div. (0->0), fcn. (1930->10), ass. (0->72)
	t352 = sin(qJ(2));
	t355 = cos(qJ(2));
	t348 = sin(pkin(10));
	t399 = cos(pkin(6));
	t378 = t348 * t399;
	t398 = cos(pkin(10));
	t338 = -t352 * t378 + t398 * t355;
	t351 = sin(qJ(3));
	t354 = cos(qJ(3));
	t384 = pkin(5) + pkin(9) + r_i_i_C(1);
	t402 = -pkin(3) * t351 + t384 * t354;
	t400 = r_i_i_C(2) + qJ(5);
	t349 = sin(pkin(6));
	t397 = t349 * t352;
	t396 = t349 * t354;
	t350 = sin(qJ(4));
	t395 = t350 * t354;
	t394 = t350 * t355;
	t353 = cos(qJ(4));
	t393 = t353 * t354;
	t392 = t353 * t355;
	t391 = t354 * t355;
	t390 = qJD(2) * t352;
	t389 = qJD(2) * t354;
	t388 = qJD(3) * t351;
	t387 = qJD(3) * t355;
	t386 = qJD(4) * t353;
	t385 = qJD(4) * t354;
	t383 = -r_i_i_C(3) - qJ(6) - pkin(4);
	t382 = t349 * t390;
	t381 = qJD(2) * t349 * t355;
	t380 = qJD(4) * t394;
	t379 = t351 * t387;
	t377 = t349 * t398;
	t372 = t399 * t398;
	t367 = t355 * t372;
	t331 = -qJD(2) * t367 + t348 * t390;
	t335 = t348 * t352 - t367;
	t374 = t335 * t385 - t331;
	t337 = t398 * t352 + t355 * t378;
	t333 = t337 * qJD(2);
	t373 = t337 * t385 - t333;
	t336 = t348 * t355 + t352 * t372;
	t361 = -t336 * t354 + t351 * t377;
	t371 = t335 * t350 - t353 * t361;
	t370 = -t335 * t353 - t350 * t361;
	t326 = t348 * t349 * t351 + t338 * t354;
	t369 = t326 * t353 + t337 * t350;
	t368 = t326 * t350 - t337 * t353;
	t366 = -t338 * t351 + t348 * t396;
	t340 = t399 * t351 + t352 * t396;
	t365 = t340 * t350 + t349 * t392;
	t364 = -t354 * pkin(3) - t384 * t351 - pkin(2);
	t363 = -t336 * t351 - t354 * t377;
	t362 = -t351 * t397 + t399 * t354;
	t360 = qJD(3) * t402;
	t332 = t336 * qJD(2);
	t359 = qJD(4) * t336 - t332 * t354 + t335 * t388;
	t334 = t338 * qJD(2);
	t358 = qJD(4) * t338 - t334 * t354 + t337 * t388;
	t357 = t400 * t350 - t383 * t353 + pkin(3);
	t356 = t350 * qJD(5) + t353 * qJD(6) + (t383 * t350 + t400 * t353) * qJD(4);
	t328 = t362 * qJD(3) + t354 * t381;
	t322 = t366 * qJD(3) - t333 * t354;
	t320 = t363 * qJD(3) - t331 * t354;
	t316 = -t365 * qJD(4) + t328 * t353 + t350 * t382;
	t315 = t328 * t350 + t340 * t386 - t349 * t380 - t353 * t382;
	t310 = -t368 * qJD(4) + t322 * t353 + t334 * t350;
	t309 = t369 * qJD(4) + t322 * t350 - t334 * t353;
	t308 = -t370 * qJD(4) + t320 * t353 + t332 * t350;
	t307 = t371 * qJD(4) + t320 * t350 - t332 * t353;
	t1 = [0, -(t337 * t393 - t338 * t350) * qJD(6) - (t337 * t395 + t338 * t353) * qJD(5) - t333 * pkin(8) + t400 * (t358 * t350 - t373 * t353) - t383 * (t373 * t350 + t358 * t353) - t337 * t360 + t364 * t334, t384 * t322 + t357 * (-t326 * qJD(3) + t333 * t351) + t356 * t366, t369 * qJD(5) - t368 * qJD(6) + t383 * t309 + t400 * t310, t309, t310; 0, -(t335 * t393 - t336 * t350) * qJD(6) - (t335 * t395 + t336 * t353) * qJD(5) - t331 * pkin(8) + t400 * (t359 * t350 - t374 * t353) - t383 * (t374 * t350 + t359 * t353) - t335 * t360 + t364 * t332, t384 * t320 + t357 * (t361 * qJD(3) + t331 * t351) + t356 * t363, t371 * qJD(5) - t370 * qJD(6) + t383 * t307 + t400 * t308, t307, t308; 0, t383 * (-t350 * t381 - t386 * t397) + (-t400 * ((qJD(2) - t385) * t392 + (t379 + (-qJD(4) + t389) * t352) * t350) + t383 * (t354 * t380 + (t352 * t389 + t379) * t353) - (-t350 * t352 - t353 * t391) * qJD(6) - (-t350 * t391 + t352 * t353) * qJD(5) + t402 * t387 + (t355 * pkin(8) + t364 * t352) * qJD(2)) * t349, t384 * t328 + t357 * (-t340 * qJD(3) - t351 * t381) + t356 * t362, -t365 * qJD(6) - (-t340 * t353 + t349 * t394) * qJD(5) + t400 * t316 + t383 * t315, t315, t316;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end