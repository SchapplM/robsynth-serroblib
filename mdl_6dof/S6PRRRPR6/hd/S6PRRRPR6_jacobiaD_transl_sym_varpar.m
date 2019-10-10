% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:56
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(11));
	t50 = sin(pkin(11));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:08
	% EndTime: 2019-10-09 22:56:08
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(11));
	t182 = cos(pkin(11));
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
	% StartTime: 2019-10-09 22:56:09
	% EndTime: 2019-10-09 22:56:09
	% DurationCPUTime: 0.50s
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
	t287 = sin(pkin(11));
	t307 = t287 * t310;
	t289 = cos(pkin(11));
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
	% StartTime: 2019-10-09 22:56:10
	% EndTime: 2019-10-09 22:56:10
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (437->91), mult. (1367->161), div. (0->0), fcn. (1417->10), ass. (0->63)
	t361 = sin(qJ(3));
	t362 = sin(qJ(2));
	t358 = sin(pkin(6));
	t364 = cos(qJ(3));
	t397 = t358 * t364;
	t400 = cos(pkin(6));
	t351 = t400 * t361 + t362 * t397;
	t363 = cos(qJ(4));
	t360 = sin(qJ(4));
	t365 = cos(qJ(2));
	t395 = t360 * t365;
	t407 = -t351 * t363 + t358 * t395;
	t391 = qJD(3) * t365;
	t406 = (qJD(2) * t364 - qJD(4)) * t362 + t361 * t391;
	t404 = pkin(9) + r_i_i_C(2);
	t405 = -pkin(3) * t361 + t404 * t364;
	t403 = r_i_i_C(1) + pkin(4);
	t401 = r_i_i_C(3) + qJ(5);
	t398 = t358 * t361;
	t396 = t360 * t364;
	t394 = qJD(2) * t362;
	t393 = qJD(2) * t365;
	t392 = qJD(3) * t361;
	t390 = qJD(4) * t364;
	t388 = t358 * t394;
	t387 = t358 * t393;
	t385 = t362 * t400;
	t384 = t365 * t400;
	t357 = sin(pkin(11));
	t382 = t357 * t385;
	t359 = cos(pkin(11));
	t381 = t359 * t384;
	t342 = -qJD(2) * t381 + t357 * t394;
	t346 = t357 * t362 - t381;
	t380 = t346 * t390 - t342;
	t348 = t357 * t384 + t359 * t362;
	t344 = t348 * qJD(2);
	t379 = t348 * t390 - t344;
	t347 = t357 * t365 + t359 * t385;
	t374 = -t347 * t364 + t359 * t398;
	t378 = t346 * t360 - t363 * t374;
	t349 = t359 * t365 - t382;
	t337 = t349 * t364 + t357 * t398;
	t377 = t337 * t363 + t348 * t360;
	t376 = (qJD(2) - t390) * t365;
	t375 = -t347 * t361 - t359 * t397;
	t373 = -t349 * t361 + t357 * t397;
	t372 = -t364 * pkin(3) - t404 * t361 - pkin(2);
	t371 = qJD(3) * t405;
	t370 = -t362 * t398 + t400 * t364;
	t369 = t401 * t360 + t403 * t363 + pkin(3);
	t343 = t347 * qJD(2);
	t368 = qJD(4) * t347 - t343 * t364 + t346 * t392;
	t345 = -qJD(2) * t382 + t359 * t393;
	t367 = qJD(4) * t349 - t345 * t364 + t348 * t392;
	t366 = qJD(5) * t360 + (-t403 * t360 + t401 * t363) * qJD(4);
	t339 = t370 * qJD(3) + t364 * t387;
	t333 = t373 * qJD(3) - t344 * t364;
	t331 = t375 * qJD(3) - t342 * t364;
	t326 = -t407 * qJD(4) + t339 * t360 - t363 * t388;
	t320 = t377 * qJD(4) + t333 * t360 - t345 * t363;
	t318 = t378 * qJD(4) + t331 * t360 - t343 * t363;
	t1 = [0, -(t348 * t396 + t349 * t363) * qJD(5) - t344 * pkin(8) + t403 * (t379 * t360 + t367 * t363) + t401 * (t367 * t360 - t379 * t363) - t348 * t371 + t372 * t345, t404 * t333 + t366 * t373 + t369 * (-t337 * qJD(3) + t344 * t361), t377 * qJD(5) + t401 * (t333 * t363 + t345 * t360 + (-t337 * t360 + t348 * t363) * qJD(4)) - t403 * t320, t320, 0; 0, -(t346 * t396 + t347 * t363) * qJD(5) - t342 * pkin(8) + t403 * (t380 * t360 + t368 * t363) + t401 * (t368 * t360 - t380 * t363) - t346 * t371 + t372 * t343, t404 * t331 + t366 * t375 + t369 * (t374 * qJD(3) + t342 * t361), t378 * qJD(5) + t401 * (t331 * t363 + t343 * t360 + (t346 * t363 + t360 * t374) * qJD(4)) - t403 * t318, t318, 0; 0, (t403 * (t360 * t376 - t406 * t363) - t401 * (t406 * t360 + t363 * t376) - (t362 * t363 - t364 * t395) * qJD(5) + t405 * t391 + (t365 * pkin(8) + t372 * t362) * qJD(2)) * t358, t404 * t339 + t366 * t370 + t369 * (-t351 * qJD(3) - t361 * t387), -t407 * qJD(5) + t401 * (t360 * t388 + t339 * t363 + (-t358 * t363 * t365 - t351 * t360) * qJD(4)) - t403 * t326, t326, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:10
	% EndTime: 2019-10-09 22:56:11
	% DurationCPUTime: 1.26s
	% Computational Cost: add. (988->143), mult. (3026->250), div. (0->0), fcn. (3289->12), ass. (0->85)
	t406 = sin(pkin(11));
	t411 = sin(qJ(2));
	t415 = cos(qJ(2));
	t461 = cos(pkin(11));
	t462 = cos(pkin(6));
	t432 = t462 * t461;
	t424 = t415 * t432;
	t391 = t406 * t411 - t424;
	t409 = sin(qJ(4));
	t413 = cos(qJ(4));
	t392 = t406 * t415 + t411 * t432;
	t410 = sin(qJ(3));
	t414 = cos(qJ(3));
	t407 = sin(pkin(6));
	t438 = t407 * t461;
	t421 = -t392 * t414 + t410 * t438;
	t363 = -t391 * t413 - t409 * t421;
	t364 = t391 * t409 - t413 * t421;
	t408 = sin(qJ(6));
	t412 = cos(qJ(6));
	t470 = ((t363 * t408 + t364 * t412) * r_i_i_C(1) + (t363 * t412 - t364 * t408) * r_i_i_C(2)) * qJD(6);
	t439 = t406 * t462;
	t394 = -t411 * t439 + t461 * t415;
	t458 = t407 * t410;
	t374 = t394 * t414 + t406 * t458;
	t393 = t461 * t411 + t415 * t439;
	t365 = t374 * t409 - t393 * t413;
	t366 = t374 * t413 + t393 * t409;
	t469 = ((t365 * t408 + t366 * t412) * r_i_i_C(1) + (t365 * t412 - t366 * t408) * r_i_i_C(2)) * qJD(6);
	t454 = t411 * t414;
	t396 = t407 * t454 + t462 * t410;
	t452 = t413 * t415;
	t377 = t396 * t409 + t407 * t452;
	t455 = t409 * t415;
	t460 = t396 * t413;
	t378 = -t407 * t455 + t460;
	t468 = ((t377 * t408 + t378 * t412) * r_i_i_C(1) + (t377 * t412 - t378 * t408) * r_i_i_C(2)) * qJD(6);
	t445 = -r_i_i_C(3) - pkin(10) + pkin(9);
	t467 = -pkin(3) * t410 + t445 * t414;
	t459 = t406 * t407;
	t457 = t409 * t411;
	t456 = t409 * t414;
	t453 = t413 * t414;
	t451 = t414 * t415;
	t450 = qJD(2) * t411;
	t449 = qJD(3) * t410;
	t448 = qJD(3) * t415;
	t447 = qJD(4) * t407;
	t446 = qJD(4) * t414;
	t444 = t411 * t458;
	t443 = t407 * t450;
	t442 = qJD(2) * t407 * t415;
	t441 = t410 * t448;
	t440 = t415 * t447;
	t435 = t414 * t438;
	t387 = -qJD(2) * t424 + t406 * t450;
	t434 = t391 * t446 - t387;
	t389 = t393 * qJD(2);
	t433 = t393 * t446 - t389;
	t425 = t408 * r_i_i_C(1) + t412 * r_i_i_C(2) + qJ(5);
	t423 = t412 * r_i_i_C(1) - t408 * r_i_i_C(2) + pkin(4) + pkin(5);
	t422 = -t414 * pkin(3) - t445 * t410 - pkin(2);
	t420 = qJD(3) * t467;
	t388 = t392 * qJD(2);
	t419 = qJD(4) * t392 - t388 * t414 + t391 * t449;
	t390 = t394 * qJD(2);
	t418 = qJD(4) * t394 - t390 * t414 + t393 * t449;
	t417 = t425 * t409 + t423 * t413 + pkin(3);
	t416 = t409 * qJD(5) + ((-t408 * t413 + t409 * t412) * r_i_i_C(1) + (-t408 * t409 - t412 * t413) * r_i_i_C(2)) * qJD(6) + (-t423 * t409 + t425 * t413) * qJD(4);
	t382 = (t413 * t451 + t457) * t407;
	t381 = (t409 * t451 - t411 * t413) * t407;
	t376 = -qJD(3) * t444 + (t462 * qJD(3) + t442) * t414;
	t370 = -t393 * t453 + t394 * t409;
	t369 = -t393 * t456 - t394 * t413;
	t368 = -t391 * t453 + t392 * t409;
	t367 = -t391 * t456 - t392 * t413;
	t362 = -t394 * t449 + (qJD(3) * t459 - t389) * t414;
	t360 = -qJD(3) * t435 - t387 * t414 - t392 * t449;
	t356 = -t377 * qJD(4) + t376 * t413 + t409 * t443;
	t355 = qJD(4) * t460 - t413 * t443 + (t376 - t440) * t409;
	t350 = -t365 * qJD(4) + t362 * t413 + t390 * t409;
	t349 = t366 * qJD(4) + t362 * t409 - t390 * t413;
	t348 = -t363 * qJD(4) + t360 * t413 + t388 * t409;
	t347 = t364 * qJD(4) + t360 * t409 - t388 * t413;
	t1 = [0, -t389 * pkin(8) + t369 * qJD(5) + t425 * (t418 * t409 - t433 * t413) + t423 * (t433 * t409 + t418 * t413) + ((t369 * t412 - t370 * t408) * r_i_i_C(1) + (-t369 * t408 - t370 * t412) * r_i_i_C(2)) * qJD(6) - t393 * t420 + t422 * t390, t445 * t362 + t417 * (-t374 * qJD(3) + t389 * t410) + t416 * (-t394 * t410 + t414 * t459), t366 * qJD(5) - t423 * t349 + t425 * t350 + t469, t349, (t349 * t412 - t350 * t408) * r_i_i_C(1) + (-t349 * t408 - t350 * t412) * r_i_i_C(2) - t469; 0, -t387 * pkin(8) + t367 * qJD(5) + t425 * (t419 * t409 - t434 * t413) + t423 * (t434 * t409 + t419 * t413) + ((t367 * t412 - t368 * t408) * r_i_i_C(1) + (-t367 * t408 - t368 * t412) * r_i_i_C(2)) * qJD(6) - t391 * t420 + t422 * t388, t445 * t360 + t417 * (t421 * qJD(3) + t387 * t410) + t416 * (-t392 * t410 - t435), t364 * qJD(5) - t423 * t347 + t425 * t348 + t470, t347, (t347 * t412 - t348 * t408) * r_i_i_C(1) + (-t347 * t408 - t348 * t412) * r_i_i_C(2) - t470; 0, t381 * qJD(5) - t425 * (-t440 * t453 - t447 * t457) + ((t381 * t412 - t382 * t408) * r_i_i_C(1) + (-t381 * t408 - t382 * t412) * r_i_i_C(2)) * qJD(6) + (-t425 * (t409 * t441 + (t409 * t454 + t452) * qJD(2)) + t423 * ((qJD(2) - t446) * t455 + (-t441 + (-qJD(2) * t414 + qJD(4)) * t411) * t413) + t467 * t448 + (pkin(8) * t415 + t422 * t411) * qJD(2)) * t407, t445 * t376 + t417 * (-t396 * qJD(3) - t410 * t442) + t416 * (t462 * t414 - t444), t378 * qJD(5) - t423 * t355 + t425 * t356 + t468, t355, (t355 * t412 - t356 * t408) * r_i_i_C(1) + (-t355 * t408 - t356 * t412) * r_i_i_C(2) - t468;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end