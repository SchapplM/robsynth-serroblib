% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:15
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
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
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(8) + r_i_i_C(3);
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
	% StartTime: 2019-10-10 11:15:13
	% EndTime: 2019-10-10 11:15:13
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (83->38), mult. (270->73), div. (0->0), fcn. (211->6), ass. (0->31)
	t199 = sin(qJ(2));
	t202 = cos(qJ(2));
	t223 = pkin(9) + r_i_i_C(3);
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
	t1 = [t197 * r_i_i_C(1) + t196 * r_i_i_C(2) - t224 * t218 + (-pkin(8) * t200 + t207 * t203) * qJD(1), (-t203 * t206 - t223 * t220) * t202 + (t204 * t203 + t209 * t220) * t199, t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0, 0, 0; -t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t224 * t216 + (pkin(8) * t203 + t207 * t200) * qJD(1), (-t200 * t206 + t223 * t219) * t202 + (t204 * t200 - t209 * t219) * t199, -t196 * r_i_i_C(1) + t197 * r_i_i_C(2), 0, 0, 0; 0, -t210 * t214 + (-t209 * t199 + t213) * qJD(2), (t198 * t215 - t201 * t217) * r_i_i_C(2) + (-t198 * t217 - t201 * t215) * r_i_i_C(1), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:15
	% EndTime: 2019-10-10 11:15:16
	% DurationCPUTime: 0.94s
	% Computational Cost: add. (291->88), mult. (998->154), div. (0->0), fcn. (918->10), ass. (0->56)
	t325 = sin(qJ(2));
	t320 = sin(pkin(10));
	t322 = cos(pkin(10));
	t324 = sin(qJ(3));
	t327 = cos(qJ(3));
	t323 = cos(pkin(6));
	t370 = t323 * t324;
	t321 = sin(pkin(6));
	t375 = r_i_i_C(3) + qJ(4);
	t382 = t321 * t375;
	t333 = -t324 * t382 + (t320 * t370 - t322 * t327) * r_i_i_C(1) + (t320 * t327 + t322 * t370) * r_i_i_C(2) - t327 * pkin(3);
	t332 = -pkin(2) + t333;
	t328 = cos(qJ(2));
	t340 = pkin(9) + (r_i_i_C(1) * t320 + r_i_i_C(2) * t322) * t321;
	t385 = t375 * t323 + t340;
	t335 = t385 * t328;
	t331 = t332 * t325 + t335;
	t386 = t323 * qJ(4) + t340;
	t368 = t323 * t327;
	t334 = t327 * t382 - (t320 * t368 + t322 * t324) * r_i_i_C(1) - (-t320 * t324 + t322 * t368) * r_i_i_C(2) - t324 * pkin(3);
	t356 = qJD(4) * t321;
	t384 = t334 * qJD(3) + t324 * t356;
	t377 = t325 * pkin(2);
	t326 = sin(qJ(1));
	t362 = qJD(1) * t328;
	t346 = -qJD(3) + t362;
	t347 = qJD(3) * t328 - qJD(1);
	t329 = cos(qJ(1));
	t358 = qJD(2) * t329;
	t351 = t325 * t358;
	t365 = t327 * t329;
	t311 = -t347 * t365 + (t346 * t326 + t351) * t324;
	t374 = t311 * t321;
	t373 = t320 * t323;
	t372 = t322 * t323;
	t369 = t323 * t325;
	t367 = t324 * t329;
	t366 = t326 * t328;
	t364 = t328 * t329;
	t363 = qJD(1) * t326;
	t361 = qJD(1) * t329;
	t360 = qJD(2) * t326;
	t359 = qJD(2) * t328;
	t357 = qJD(3) * t327;
	t355 = t323 * qJD(4);
	t354 = -pkin(2) * t328 - pkin(1);
	t352 = t325 * t355;
	t350 = qJD(3) * t367;
	t338 = t324 * t361 + t326 * t357;
	t336 = t385 * t325;
	t330 = t328 * t355 - t384 * t325 + (t332 * t328 - t336) * qJD(2);
	t314 = -t346 * t365 + (qJD(2) * t325 * t327 + t347 * t324) * t326;
	t313 = -t325 * t324 * t360 - t327 * t363 + t338 * t328 - t350;
	t312 = t328 * t350 + (t326 * t362 + t351) * t327 - t338;
	t310 = -t374 + (-t325 * t363 + t328 * t358) * t323;
	t1 = [(t313 * t373 + t314 * t322) * r_i_i_C(1) + (t313 * t372 - t314 * t320) * r_i_i_C(2) + t314 * pkin(3) - t326 * t352 + (-(t324 * t366 + t365) * qJD(4) - t375 * t313) * t321 + (-t335 + t377) * t360 + (-t326 * pkin(8) + (-t336 + t354) * t329) * qJD(1), t330 * t329 - t331 * t363, (t311 * t322 + t312 * t373) * r_i_i_C(1) + (-t311 * t320 + t312 * t372) * r_i_i_C(2) + t311 * pkin(3) + (-(-t326 * t324 - t327 * t364) * qJD(4) - t375 * t312) * t321, t310, 0, 0; (t311 * t373 - t312 * t322) * r_i_i_C(1) + (t311 * t372 + t312 * t320) * r_i_i_C(2) + t310 * r_i_i_C(3) - t312 * pkin(3) - qJ(4) * t374 + (-(-t324 * t364 + t326 * t327) * t321 + t329 * t369) * qJD(4) + (t386 * t328 - t377) * t358 + (t329 * pkin(8) + (-t386 * t325 + t354) * t326) * qJD(1), t330 * t326 + t331 * t361, (-t313 * t322 + t314 * t373) * r_i_i_C(1) + (t313 * t320 + t314 * t372) * r_i_i_C(2) - t313 * pkin(3) + (-(-t327 * t366 + t367) * qJD(4) - t375 * t314) * t321, t313 * t321 + (t325 * t361 + t326 * t359) * t323, 0, 0; 0, t331 * qJD(2) + t384 * t328 + t352, t334 * t359 + (t333 * qJD(3) + t327 * t356) * t325, t325 * t321 * t357 + (t321 * t324 * t328 + t369) * qJD(2), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:16
	% EndTime: 2019-10-10 11:15:17
	% DurationCPUTime: 1.45s
	% Computational Cost: add. (484->131), mult. (1653->224), div. (0->0), fcn. (1578->10), ass. (0->85)
	t408 = sin(qJ(2));
	t407 = sin(qJ(3));
	t410 = cos(qJ(3));
	t404 = sin(pkin(6));
	t476 = r_i_i_C(1) + qJ(4);
	t439 = t404 * t476;
	t428 = -t410 * pkin(3) - t407 * t439;
	t419 = -pkin(2) + t428;
	t411 = cos(qJ(2));
	t406 = cos(pkin(6));
	t485 = t476 * t406 + pkin(9);
	t434 = t485 * t411;
	t494 = t419 * t408 + t434;
	t412 = cos(qJ(1));
	t460 = t411 * t412;
	t409 = sin(qJ(1));
	t463 = t409 * t407;
	t394 = t410 * t460 + t463;
	t452 = qJD(3) * t412;
	t440 = t410 * t452;
	t444 = qJD(3) * t463;
	t448 = qJD(2) * t408 * t409;
	t390 = qJD(1) * t394 - t410 * t448 - t411 * t444 - t440;
	t403 = sin(pkin(10));
	t405 = cos(pkin(10));
	t453 = qJD(3) * t410;
	t443 = t409 * t453;
	t457 = qJD(1) * t412;
	t422 = t407 * t457 + t443;
	t441 = t407 * t452;
	t459 = qJD(1) * t409;
	t423 = t410 * t459 + t441;
	t389 = -t407 * t448 + t411 * t422 - t423;
	t456 = qJD(2) * t411;
	t447 = t409 * t456;
	t425 = t408 * t457 + t447;
	t483 = -t389 * t406 + t404 * t425;
	t492 = -t390 * t403 + t483 * t405;
	t491 = (t425 * t407 + t408 * t443) * t406 + (t411 * t457 - t448) * t404;
	t455 = qJD(2) * t412;
	t446 = t411 * t455;
	t424 = t408 * t459 - t446;
	t442 = t408 * t455;
	t458 = qJD(1) * t411;
	t426 = t409 * t458 + t442;
	t490 = (t407 * t424 - t408 * t440) * t406 + t426 * t404;
	t445 = t404 * t453;
	t450 = qJD(5) * t405;
	t451 = qJD(4) * t404;
	t470 = t403 * t410;
	t489 = t476 * (qJD(2) * t406 + t445) + t407 * (-pkin(3) * qJD(3) + t406 * t450 + t451) + qJD(2) * pkin(9) + qJD(5) * t470;
	t429 = -t406 * qJD(4) + t404 * t450;
	t486 = t429 * t408;
	t461 = t410 * t412;
	t387 = (-qJD(3) * t411 + qJD(1)) * t461 + (t442 + (-qJD(3) + t458) * t409) * t407;
	t484 = -t387 * t406 + t424 * t404;
	t479 = r_i_i_C(2) - pkin(4);
	t478 = t408 * pkin(2);
	t475 = r_i_i_C(3) + qJ(5);
	t471 = t403 * t406;
	t469 = t405 * t406;
	t468 = t406 * t407;
	t467 = t406 * t408;
	t466 = t406 * t410;
	t465 = t407 * t411;
	t464 = t408 * t410;
	t462 = t410 * t411;
	t454 = qJD(3) * t408;
	t449 = -pkin(2) * t411 - pkin(1);
	t438 = qJ(4) * t406 + pkin(9);
	t432 = -t403 * t407 + t405 * t466;
	t431 = t404 * t411 + t407 * t467;
	t421 = t432 * t411;
	t420 = (t403 * t466 + t405 * t407) * t411;
	t416 = t410 * t447 + (t410 * t457 - t444) * t408;
	t415 = t408 * t423 - t410 * t446;
	t414 = qJD(2) * t419 - t429;
	t413 = -t489 * t408 + t414 * t411;
	t393 = -t407 * t460 + t409 * t410;
	t392 = t407 * t412 - t409 * t462;
	t391 = t411 * t463 + t461;
	t388 = t410 * t426 + t411 * t441 - t422;
	t382 = -t387 * t404 - t406 * t424;
	t370 = -t388 * t403 + t484 * t405;
	t1 = [-(t391 * t469 - t392 * t403) * qJD(5) - t390 * pkin(3) - t479 * (-t390 * t405 - t403 * t483) + t475 * t492 + (-t391 * qJD(4) - t476 * t389) * t404 + (t486 + (-t434 + t478) * qJD(2)) * t409 + (-t409 * pkin(8) + (-t408 * t485 + t449) * t412) * qJD(1), -t479 * (-t490 * t403 + t415 * t405) + t475 * (t415 * t403 + t490 * t405) - t494 * t459 + t413 * t412, -(-t393 * t403 - t394 * t469) * qJD(5) + t387 * pkin(3) - t479 * (t387 * t405 + t388 * t471) + t475 * (t387 * t403 - t388 * t469) + (t394 * qJD(4) - t476 * t388) * t404, t382, t370, 0; t382 * r_i_i_C(1) - (t393 * t469 - t394 * t403) * qJD(5) - t388 * pkin(3) + (-qJ(4) * t387 - qJD(4) * t393) * t404 - t479 * (-t388 * t405 - t484 * t403) + t475 * t370 + (-t486 + (t411 * t438 - t478) * qJD(2)) * t412 + (t412 * pkin(8) + (-t408 * t438 + t449) * t409) * qJD(1), t479 * (-t491 * t403 + t416 * t405) - t475 * (t416 * t403 + t491 * t405) + t494 * t457 + t413 * t409, -(t391 * t403 + t392 * t469) * qJD(5) - t389 * pkin(3) - t479 * (-t389 * t405 - t390 * t471) + t475 * (-t389 * t403 + t390 * t469) + (-t392 * qJD(4) + t476 * t390) * t404, t389 * t404 + t406 * t425, -t492, 0; 0, t479 * (qJD(3) * t420 + (-t431 * t403 + t405 * t464) * qJD(2)) - t475 * (-qJD(3) * t421 + (t403 * t464 + t431 * t405) * qJD(2)) + t414 * t408 + t489 * t411, t479 * ((-t403 * t468 + t405 * t410) * t454 + qJD(2) * t420) - t475 * ((t405 * t468 + t470) * t454 - qJD(2) * t421) + (-pkin(3) * t407 + t410 * t439) * t456 + (qJD(3) * t428 + qJD(5) * t432 + t410 * t451) * t408, t408 * t445 + (t404 * t465 + t467) * qJD(2), t432 * t454 + (t403 * t462 + (-t404 * t408 + t406 * t465) * t405) * qJD(2), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:16
	% EndTime: 2019-10-10 11:15:17
	% DurationCPUTime: 1.59s
	% Computational Cost: add. (647->159), mult. (2194->269), div. (0->0), fcn. (2123->10), ass. (0->104)
	t428 = sin(qJ(2));
	t424 = sin(pkin(6));
	t427 = sin(qJ(3));
	t507 = t424 * t427;
	t430 = cos(qJ(3));
	t518 = t430 * pkin(3);
	t448 = qJ(4) * t507 + pkin(2) + t518;
	t431 = cos(qJ(2));
	t426 = cos(pkin(6));
	t463 = qJ(4) * t426 + pkin(9);
	t528 = t463 * t431;
	t536 = -t448 * t428 + t528;
	t432 = cos(qJ(1));
	t487 = qJD(2) * t432;
	t470 = t431 * t487;
	t429 = sin(qJ(1));
	t493 = qJD(1) * t429;
	t444 = t428 * t493 - t470;
	t484 = qJD(3) * t432;
	t464 = t430 * t484;
	t435 = t444 * t427 - t428 * t464;
	t472 = t428 * t487;
	t492 = qJD(1) * t431;
	t446 = t429 * t492 + t472;
	t535 = t446 * t424 + t435 * t426;
	t485 = qJD(3) * t430;
	t466 = t429 * t485;
	t491 = qJD(1) * t432;
	t442 = t427 * t491 + t466;
	t465 = t427 * t484;
	t443 = t430 * t493 + t465;
	t489 = qJD(2) * t429;
	t474 = t428 * t489;
	t398 = -t427 * t474 + t442 * t431 - t443;
	t494 = t431 * t432;
	t498 = t429 * t427;
	t407 = t430 * t494 + t498;
	t467 = qJD(3) * t498;
	t499 = t428 * t430;
	t473 = qJD(2) * t499;
	t399 = t407 * qJD(1) - t429 * t473 - t431 * t467 - t464;
	t423 = sin(pkin(10));
	t425 = cos(pkin(10));
	t488 = qJD(2) * t431;
	t471 = t429 * t488;
	t477 = t428 * t491;
	t445 = t471 + t477;
	t534 = -t399 * t423 + (-t398 * t426 + t445 * t424) * t425;
	t454 = qJD(5) * t425 - qJD(6) * t423;
	t469 = t424 * t485;
	t483 = t424 * qJD(4);
	t490 = qJD(2) * t426;
	t533 = (-qJD(3) * pkin(3) + t454 * t426 + t483) * t427 + (t469 + t490) * qJ(4) + (qJD(5) * t423 + qJD(6) * t425) * t430 + qJD(2) * pkin(9);
	t437 = t426 * qJD(4) - t454 * t424;
	t532 = qJD(2) * t528 + (-qJD(2) * pkin(2) + t437) * t428;
	t526 = (t430 * t491 - t467) * t428;
	t495 = t430 * t432;
	t396 = (-qJD(3) * t431 + qJD(1)) * t495 + (t472 + (-qJD(3) + t492) * t429) * t427;
	t525 = -t396 * t426 + t444 * t424;
	t523 = t398 * t424 + t445 * t426;
	t459 = t428 * t466;
	t521 = t445 * t427 + t459;
	t519 = r_i_i_C(1) + pkin(5);
	t517 = r_i_i_C(2) + qJ(5);
	t510 = t423 * t424;
	t509 = t423 * t426;
	t508 = t423 * t427;
	t506 = t424 * t428;
	t505 = t425 * t426;
	t503 = t426 * t427;
	t502 = t426 * t430;
	t501 = t427 * t428;
	t500 = t427 * t431;
	t497 = t429 * t431;
	t496 = t430 * t431;
	t486 = qJD(3) * t428;
	t481 = qJ(4) + t519;
	t480 = r_i_i_C(3) + qJ(6) + pkin(4);
	t479 = t425 * t496;
	t478 = t426 * t501;
	t475 = t431 * t491;
	t468 = t428 * t485;
	t462 = t423 * t477;
	t461 = t490 * t508;
	t460 = t488 * t510;
	t452 = t425 * t502 - t508;
	t451 = t423 * t502 + t425 * t427;
	t450 = t426 * t500 - t506;
	t447 = -t474 + t475;
	t441 = t452 * t431;
	t440 = t451 * t431;
	t439 = t398 * t509 - t399 * t425 - t424 * t462 - t429 * t460;
	t438 = -pkin(2) * t431 - t463 * t428 - pkin(1);
	t436 = t443 * t428 - t430 * t470;
	t434 = -t448 * qJD(2) + t437;
	t433 = -t533 * t428 + t434 * t431;
	t406 = -t427 * t494 + t429 * t430;
	t405 = t427 * t432 - t429 * t496;
	t404 = t427 * t497 + t495;
	t397 = t446 * t430 + t431 * t465 - t442;
	t388 = -t396 * t424 - t444 * t426;
	t377 = -t397 * t425 - t525 * t423;
	t376 = -t397 * t423 + t525 * t425;
	t1 = [-(-t404 * t509 - t405 * t425) * qJD(6) - (t404 * t505 - t405 * t423) * qJD(5) - t399 * pkin(3) + (-qJ(4) * t398 - qJD(4) * t404) * t424 - t519 * t523 + t517 * t534 + t480 * t439 - t532 * t429 + (-t429 * pkin(8) + t438 * t432) * qJD(1), t519 * (t435 * t424 - t446 * t426) + t517 * (t436 * t423 + t535 * t425) + t480 * (-t535 * t423 + t436 * t425) - t536 * t493 + t433 * t432, -(-t406 * t425 + t407 * t509) * qJD(6) - (-t406 * t423 - t407 * t505) * qJD(5) + t396 * pkin(3) + t517 * (t396 * t423 - t397 * t505) + t480 * (t396 * t425 + t397 * t509) + (t407 * qJD(4) - t481 * t397) * t424, t388, t376, t377; -(-t406 * t509 - t407 * t425) * qJD(6) - (t406 * t505 - t407 * t423) * qJD(5) - t397 * pkin(3) + (-qJ(4) * t396 - qJD(4) * t406) * t424 + t519 * t388 + t517 * t376 + t480 * t377 + t532 * t432 + (t432 * pkin(8) + t438 * t429) * qJD(1), t519 * (-t521 * t424 + t447 * t426) - t517 * ((t430 * t471 + t526) * t423 + (t447 * t424 + t521 * t426) * t425) - t480 * (-t459 * t509 - t461 * t497 - t462 * t503 - t475 * t510 + (t423 * t506 + t479) * t489 + t425 * t526) + t536 * t491 + t433 * t429, -(t404 * t425 - t405 * t509) * qJD(6) - (t404 * t423 + t405 * t505) * qJD(5) - t398 * pkin(3) + t517 * (-t398 * t423 + t399 * t505) + t480 * (-t398 * t425 - t399 * t509) + (-t405 * qJD(4) + t481 * t399) * t424, t523, -t534, -t439; 0, t519 * (t431 * t469 + (-t424 * t501 + t426 * t431) * qJD(2)) - t517 * (-qJD(3) * t441 + (t423 * t499 + (t424 * t431 + t478) * t425) * qJD(2)) - t480 * (qJD(3) * t440 + t425 * t473 - t428 * t461 - t460) + t434 * t428 + t533 * t431, -t517 * ((t423 * t430 + t425 * t503) * t486 - qJD(2) * t441) - t480 * (-qJD(3) * t423 * t478 + qJD(2) * t440 + t425 * t468) + (t481 * t430 * t424 - pkin(3) * t427) * t488 + (-t451 * qJD(6) + t452 * qJD(5) + t430 * t483 + (-t481 * t507 - t518) * qJD(3)) * t428, t424 * t468 + (t424 * t500 + t426 * t428) * qJD(2), t452 * t486 + (t423 * t496 + t450 * t425) * qJD(2), -t451 * t486 + (-t450 * t423 + t479) * qJD(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end