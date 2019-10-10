% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR15
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR15_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR15_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR15_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
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
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (35->18), mult. (110->35), div. (0->0), fcn. (94->6), ass. (0->20)
	t136 = sin(pkin(6));
	t151 = t136 * (pkin(9) + r_i_i_C(3));
	t138 = sin(qJ(2));
	t139 = sin(qJ(1));
	t149 = t138 * t139;
	t141 = cos(qJ(1));
	t148 = t138 * t141;
	t140 = cos(qJ(2));
	t147 = t139 * t140;
	t146 = t140 * t141;
	t137 = cos(pkin(6));
	t145 = -t137 * t146 + t149;
	t144 = t137 * t147 + t148;
	t143 = t137 * t148 + t147;
	t142 = t137 * t149 - t146;
	t135 = t142 * qJD(1) + t145 * qJD(2);
	t134 = t144 * qJD(1) + t143 * qJD(2);
	t133 = t143 * qJD(1) + t144 * qJD(2);
	t132 = t145 * qJD(1) + t142 * qJD(2);
	t1 = [t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(1) * t141 - t139 * t151) * qJD(1), t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0, 0, 0, 0; -t133 * r_i_i_C(1) + t132 * r_i_i_C(2) + (-pkin(1) * t139 + t141 * t151) * qJD(1), -t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t138 - r_i_i_C(2) * t140) * t136 * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:25
	% EndTime: 2019-10-10 12:18:25
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (229->67), mult. (725->118), div. (0->0), fcn. (712->10), ass. (0->53)
	t302 = cos(pkin(7));
	t304 = sin(qJ(3));
	t307 = cos(qJ(3));
	t321 = r_i_i_C(1) * t304 + r_i_i_C(2) * t307;
	t300 = sin(pkin(7));
	t340 = pkin(10) + r_i_i_C(3);
	t342 = t340 * t300;
	t345 = -t321 * t302 + t342;
	t303 = cos(pkin(6));
	t305 = sin(qJ(2));
	t309 = cos(qJ(1));
	t329 = t309 * t305;
	t306 = sin(qJ(1));
	t308 = cos(qJ(2));
	t331 = t306 * t308;
	t292 = t303 * t329 + t331;
	t317 = t303 * t331 + t329;
	t289 = t317 * qJD(1) + t292 * qJD(2);
	t301 = sin(pkin(6));
	t323 = qJD(1) * t300 * t301;
	t313 = -qJD(3) * t292 - t289 * t302 + t306 * t323;
	t344 = t313 * r_i_i_C(1);
	t341 = t340 * t302 + pkin(9);
	t339 = t301 * t306;
	t338 = t301 * t309;
	t337 = t302 * t304;
	t336 = t302 * t307;
	t335 = t304 * t305;
	t334 = t304 * t308;
	t333 = t305 * t307;
	t332 = t306 * t305;
	t330 = t307 * t308;
	t328 = t309 * t308;
	t327 = qJD(3) * t300;
	t324 = t303 * t332;
	t322 = t327 * t338;
	t290 = -qJD(1) * t324 - qJD(2) * t332 + (qJD(2) * t303 + qJD(1)) * t328;
	t291 = -t303 * t328 + t332;
	t320 = qJD(3) * t291 * t302 - t290;
	t319 = t307 * r_i_i_C(1) - t304 * r_i_i_C(2) + pkin(2);
	t318 = t300 * t339 - t302 * t317;
	t316 = t324 - t328;
	t287 = t291 * qJD(1) + t316 * qJD(2);
	t315 = t287 * t302 + t309 * t323;
	t314 = t320 + t322;
	t312 = t313 * r_i_i_C(2);
	t311 = (-t302 * t333 - t334) * r_i_i_C(1) + (t302 * t335 - t330) * r_i_i_C(2);
	t310 = (-t302 * t334 - t333) * r_i_i_C(1) + (-t302 * t330 + t335) * r_i_i_C(2);
	t295 = t307 * t322;
	t288 = t292 * qJD(1) + t317 * qJD(2);
	t286 = -t288 * t307 + t315 * t304 + (t304 * t316 + t318 * t307) * qJD(3);
	t285 = t288 * t304 + t315 * t307 + (-t318 * t304 + t307 * t316) * qJD(3);
	t1 = [-t290 * pkin(2) + t295 * r_i_i_C(1) - t289 * t342 + (-t309 * pkin(1) - t341 * t339) * qJD(1) + (t320 * r_i_i_C(1) - t312) * t307 + (-t314 * r_i_i_C(2) - t344) * t304, t319 * t287 - t345 * t288 + ((t304 * t317 + t316 * t336) * r_i_i_C(1) + (t307 * t317 - t316 * t337) * r_i_i_C(2)) * qJD(3), t285 * r_i_i_C(1) - t286 * r_i_i_C(2), 0, 0, 0; -t288 * pkin(2) + t286 * r_i_i_C(1) + t285 * r_i_i_C(2) - t287 * t342 + (-t306 * pkin(1) + t341 * t338) * qJD(1), -t319 * t289 + t345 * t290 + ((t291 * t304 - t292 * t336) * r_i_i_C(1) + (t291 * t307 + t292 * t337) * r_i_i_C(2)) * qJD(3), t295 * r_i_i_C(2) + (t320 * r_i_i_C(2) + t344) * t307 + (t314 * r_i_i_C(1) - t312) * t304, 0, 0, 0; 0, (t311 * qJD(3) + (-t305 * pkin(2) + t308 * t342 + t310) * qJD(2)) * t301, -t321 * t303 * t327 + (t311 * qJD(2) + t310 * qJD(3)) * t301, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:26
	% EndTime: 2019-10-10 12:18:26
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (470->85), mult. (1462->141), div. (0->0), fcn. (1486->10), ass. (0->59)
	t399 = cos(pkin(6));
	t401 = sin(qJ(2));
	t402 = sin(qJ(1));
	t429 = t402 * t401;
	t422 = t399 * t429;
	t404 = cos(qJ(2));
	t405 = cos(qJ(1));
	t425 = t405 * t404;
	t385 = -qJD(1) * t422 - qJD(2) * t429 + (qJD(2) * t399 + qJD(1)) * t425;
	t400 = sin(qJ(3));
	t403 = cos(qJ(3));
	t426 = t405 * t401;
	t428 = t402 * t404;
	t388 = t399 * t426 + t428;
	t387 = -t399 * t425 + t429;
	t396 = sin(pkin(7));
	t398 = cos(pkin(7));
	t397 = sin(pkin(6));
	t435 = t397 * t405;
	t417 = t387 * t398 + t396 * t435;
	t406 = -t388 * t403 + t417 * t400;
	t410 = t399 * t428 + t426;
	t384 = t410 * qJD(1) + t388 * qJD(2);
	t421 = qJD(1) * t396 * t397;
	t408 = -t384 * t398 + t402 * t421;
	t369 = -t406 * qJD(3) + t385 * t400 - t408 * t403;
	t407 = t388 * t400 + t417 * t403;
	t454 = -t407 * qJD(3) + t385 * t403 + t408 * t400;
	t446 = pkin(10) + r_i_i_C(1);
	t453 = t446 * t396;
	t450 = t446 * t398 + pkin(9);
	t409 = t422 - t425;
	t436 = t397 * t402;
	t415 = t396 * t436 - t398 * t410;
	t448 = t415 * t400 - t403 * t409;
	t445 = r_i_i_C(2) - pkin(3);
	t444 = r_i_i_C(3) + qJ(4);
	t439 = t409 * t400;
	t437 = t396 * t399;
	t434 = t398 * t400;
	t433 = t398 * t403;
	t432 = t400 * t401;
	t431 = t400 * t404;
	t430 = t401 * t403;
	t427 = t403 * t404;
	t420 = qJD(3) * t437;
	t419 = t405 * t421;
	t418 = t387 * t400 - t388 * t433;
	t416 = t400 * t410 + t409 * t433;
	t414 = t398 * t427 - t432;
	t413 = t398 * t430 + t431;
	t412 = t398 * t431 + t430;
	t411 = t398 * t432 - t427;
	t383 = t388 * qJD(1) + t410 * qJD(2);
	t382 = t387 * qJD(1) + t409 * qJD(2);
	t377 = t400 * t420 + (t413 * qJD(2) + t412 * qJD(3)) * t397;
	t368 = qJD(3) * t439 + (t382 * t398 + t419) * t400 + (t415 * qJD(3) - t383) * t403;
	t367 = t448 * qJD(3) - t382 * t433 - t383 * t400 - t403 * t419;
	t1 = [-t407 * qJD(4) - t385 * pkin(2) - t384 * t453 + t445 * t454 - t444 * t369 + (-t405 * pkin(1) - t450 * t436) * qJD(1), -t416 * qJD(4) + t382 * pkin(2) - t383 * t453 - t445 * (t416 * qJD(3) + t382 * t403 + t383 * t434) + t444 * (-t383 * t433 + t382 * t400 + (-t403 * t410 + t409 * t434) * qJD(3)), t448 * qJD(4) + t445 * t367 + t444 * t368, t367, 0, 0; -(t415 * t403 + t439) * qJD(4) - t383 * pkin(2) - t382 * t453 - t445 * t368 + t444 * t367 + (-t402 * pkin(1) + t450 * t435) * qJD(1), -t418 * qJD(4) - t384 * pkin(2) + t385 * t453 - t445 * (t418 * qJD(3) - t384 * t403 - t385 * t434) + t444 * (t385 * t433 - t384 * t400 + (-t387 * t403 - t388 * t434) * qJD(3)), -t406 * qJD(4) + t445 * t369 + t444 * t454, t369, 0, 0; 0, (t445 * (t412 * qJD(2) + t413 * qJD(3)) - t444 * (-t414 * qJD(2) + t411 * qJD(3)) + t413 * qJD(4) + (-t401 * pkin(2) + t404 * t453) * qJD(2)) * t397, -(-t412 * t397 - t400 * t437) * qJD(4) + t444 * (t403 * t420 + (-t411 * qJD(2) + t414 * qJD(3)) * t397) + t445 * t377, t377, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:27
	% EndTime: 2019-10-10 12:18:28
	% DurationCPUTime: 1.30s
	% Computational Cost: add. (901->146), mult. (2800->248), div. (0->0), fcn. (2920->12), ass. (0->91)
	t483 = sin(qJ(2));
	t484 = sin(qJ(1));
	t487 = cos(qJ(2));
	t488 = cos(qJ(1));
	t540 = cos(pkin(6));
	t512 = t488 * t540;
	t465 = t483 * t512 + t484 * t487;
	t513 = t484 * t540;
	t492 = t488 * t483 + t487 * t513;
	t454 = t492 * qJD(1) + t465 * qJD(2);
	t507 = t483 * t513;
	t523 = qJD(2) * t483;
	t525 = t488 * t487;
	t455 = -qJD(1) * t507 - t484 * t523 + (qJD(2) * t540 + qJD(1)) * t525;
	t464 = t484 * t483 - t487 * t512;
	t482 = sin(qJ(3));
	t486 = cos(qJ(3));
	t478 = sin(pkin(7));
	t479 = sin(pkin(6));
	t524 = qJD(1) * t479;
	t515 = t484 * t524;
	t510 = t478 * t515;
	t532 = t479 * t488;
	t518 = t478 * t532;
	t522 = qJD(3) * t482;
	t480 = cos(pkin(7));
	t530 = t480 * t486;
	t531 = t480 * t482;
	t535 = t465 * t486;
	t430 = (t464 * t531 - t535) * qJD(3) - t454 * t530 - t455 * t482 + t486 * t510 + t518 * t522;
	t538 = t454 * t478;
	t445 = t480 * t515 + t538;
	t481 = sin(qJ(5));
	t485 = cos(qJ(5));
	t549 = t430 * t481 - t445 * t485;
	t548 = t430 * t485 + t445 * t481;
	t500 = t464 * t480 + t518;
	t536 = t465 * t482;
	t545 = (t500 * t486 + t536) * qJD(3) - (-t454 * t480 + t510) * t482 - t455 * t486;
	t438 = t464 * t530 + t486 * t518 + t536;
	t458 = -t464 * t478 + t480 * t532;
	t544 = -t438 * t485 - t458 * t481;
	t543 = t438 * t481 - t458 * t485;
	t493 = t507 - t525;
	t452 = t464 * qJD(1) + t493 * qJD(2);
	t539 = t452 * t478;
	t534 = t478 * t479;
	t533 = t479 * t484;
	t529 = t482 * t483;
	t528 = t482 * t487;
	t527 = t483 * t486;
	t526 = t486 * t487;
	t521 = qJD(5) * t481;
	t520 = qJD(5) * t485;
	t519 = r_i_i_C(3) + pkin(11) + pkin(3);
	t517 = t480 * t526;
	t516 = pkin(10) * t480 + pkin(9);
	t514 = t488 * t524;
	t511 = t540 * t478;
	t509 = t478 * t514;
	t508 = t523 * t534;
	t506 = qJD(3) * t511;
	t505 = t485 * r_i_i_C(1) - t481 * r_i_i_C(2);
	t504 = -t481 * r_i_i_C(1) - t485 * r_i_i_C(2);
	t503 = qJ(4) - t504;
	t502 = pkin(4) + pkin(10) + t505;
	t450 = -t464 * t482 + t465 * t530;
	t451 = -t482 * t492 - t493 * t530;
	t499 = t478 * t533 - t480 * t492;
	t498 = t517 - t529;
	t497 = t480 * t527 + t528;
	t496 = t480 * t528 + t527;
	t495 = qJD(5) * t504;
	t491 = t505 * qJD(5) + qJD(4);
	t489 = t499 * t482 - t486 * t493;
	t463 = t540 * t480 - t487 * t534;
	t461 = t497 * t479;
	t460 = t478 * t492 + t480 * t533;
	t456 = -t498 * t479 - t486 * t511;
	t453 = t465 * qJD(1) + t492 * qJD(2);
	t446 = (-qJD(2) * t517 - qJD(3) * t526 + (qJD(3) * t480 + qJD(2)) * t529) * t479;
	t443 = t480 * t514 - t539;
	t441 = -t482 * t493 - t499 * t486;
	t436 = t482 * t506 + (t497 * qJD(2) + t496 * qJD(3)) * t479;
	t434 = t455 * t530 - t454 * t482 + (-t464 * t486 - t465 * t531) * qJD(3);
	t432 = -t453 * t530 + t452 * t482 + (-t486 * t492 + t493 * t531) * qJD(3);
	t427 = t493 * t522 + (t452 * t480 + t509) * t482 + (t499 * qJD(3) - t453) * t486;
	t426 = t489 * qJD(3) - t452 * t530 - t453 * t482 - t486 * t509;
	t425 = t426 * t481 + t443 * t485 + (t441 * t485 - t460 * t481) * qJD(5);
	t424 = t426 * t485 - t443 * t481 + (-t441 * t481 - t460 * t485) * qJD(5);
	t1 = [t549 * r_i_i_C(1) + t548 * r_i_i_C(2) - t445 * pkin(4) + t430 * qJ(4) - t438 * qJD(4) - t455 * pkin(2) - pkin(10) * t538 + (t544 * r_i_i_C(1) + t543 * r_i_i_C(2)) * qJD(5) + t519 * t545 + (-t488 * pkin(1) - t516 * t533) * qJD(1), (t432 * t481 + t451 * t520) * r_i_i_C(1) + (t432 * t485 - t451 * t521) * r_i_i_C(2) + t432 * qJ(4) + t451 * qJD(4) + t452 * pkin(2) + t519 * (-t451 * qJD(3) + t452 * t486 + t453 * t531) + (-t502 * t453 - t493 * t495) * t478, -t519 * t426 + t503 * t427 + t491 * t489, t426, t424 * r_i_i_C(1) - t425 * r_i_i_C(2), 0; -pkin(10) * t539 - t453 * pkin(2) + t443 * pkin(4) + t425 * r_i_i_C(1) + t424 * r_i_i_C(2) + t426 * qJ(4) + t441 * qJD(4) + t519 * t427 + (-pkin(1) * t484 + t516 * t532) * qJD(1), (t434 * t481 + t450 * t520) * r_i_i_C(1) + (t434 * t485 - t450 * t521) * r_i_i_C(2) + t434 * qJ(4) + t450 * qJD(4) - t454 * pkin(2) + t519 * (-t450 * qJD(3) - t454 * t486 - t455 * t531) + (t502 * t455 + t465 * t495) * t478, t491 * (-t500 * t482 + t535) - t503 * t545 + t519 * t430, -t430, -t548 * r_i_i_C(1) + t549 * r_i_i_C(2) + (-t543 * r_i_i_C(1) + t544 * r_i_i_C(2)) * qJD(5), 0; 0, (-t446 * t481 + t461 * t520) * r_i_i_C(1) + (-t446 * t485 - t461 * t521) * r_i_i_C(2) - t446 * qJ(4) + t461 * qJD(4) + (-t519 * (t496 * qJD(2) + t497 * qJD(3)) - pkin(2) * t523 + (t502 * t487 * qJD(2) + t483 * t495) * t478) * t479, t491 * (t496 * t479 + t482 * t511) + t503 * (t486 * t506 + (t498 * qJD(3) + (-t480 * t529 + t526) * qJD(2)) * t479) - t519 * t436, t436, (t436 * t485 - t481 * t508) * r_i_i_C(1) + (-t436 * t481 - t485 * t508) * r_i_i_C(2) + ((-t456 * t481 - t463 * t485) * r_i_i_C(1) + (-t456 * t485 + t463 * t481) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:33
	% EndTime: 2019-10-10 12:18:35
	% DurationCPUTime: 2.42s
	% Computational Cost: add. (2098->210), mult. (6474->344), div. (0->0), fcn. (7012->14), ass. (0->124)
	t690 = cos(qJ(2));
	t691 = cos(qJ(1));
	t752 = cos(pkin(6));
	t724 = t691 * t752;
	t686 = sin(qJ(2));
	t753 = sin(qJ(1));
	t727 = t753 * t686;
	t663 = -t690 * t724 + t727;
	t664 = t686 * t724 + t690 * t753;
	t685 = sin(qJ(3));
	t689 = cos(qJ(3));
	t680 = sin(pkin(7));
	t681 = sin(pkin(6));
	t743 = t681 * t691;
	t732 = t680 * t743;
	t682 = cos(pkin(7));
	t742 = t682 * t685;
	t633 = t663 * t742 - t664 * t689 + t685 * t732;
	t716 = t752 * t753;
	t694 = t691 * t686 + t690 * t716;
	t650 = qJD(1) * t694 + qJD(2) * t664;
	t707 = t686 * t716;
	t736 = t691 * t690;
	t651 = -qJD(1) * t707 - qJD(2) * t727 + (qJD(2) * t752 + qJD(1)) * t736;
	t729 = t681 * t753;
	t717 = qJD(1) * t729;
	t709 = t680 * t717;
	t741 = t682 * t689;
	t609 = t633 * qJD(3) - t650 * t741 - t651 * t685 + t689 * t709;
	t750 = t650 * t680;
	t638 = t682 * t717 + t750;
	t684 = sin(qJ(5));
	t688 = cos(qJ(5));
	t767 = t609 * t684 - t638 * t688;
	t766 = -t609 * t688 - t638 * t684;
	t748 = t664 * t685;
	t757 = t689 * (t663 * t682 + t732) + t748;
	t765 = t757 * qJD(3) - (-t650 * t682 + t709) * t685 - t651 * t689;
	t683 = sin(qJ(6));
	t764 = t633 * t683;
	t687 = cos(qJ(6));
	t763 = t633 * t687;
	t654 = -t663 * t680 + t682 * t743;
	t760 = t654 * t684;
	t759 = t654 * t688;
	t734 = (pkin(4) + pkin(10)) * t680;
	t756 = pkin(3) + pkin(11);
	t754 = r_i_i_C(3) + pkin(12);
	t695 = t707 - t736;
	t648 = qJD(1) * t663 + qJD(2) * t695;
	t751 = t648 * t680;
	t747 = t695 * t685;
	t746 = t680 * t681;
	t745 = t680 * t684;
	t744 = t680 * t688;
	t740 = t685 * t686;
	t739 = t685 * t690;
	t738 = t686 * t689;
	t737 = t689 * t690;
	t735 = qJD(2) * t681;
	t733 = t686 * t746;
	t731 = t681 * t740;
	t730 = t682 * t737;
	t728 = t682 * t753;
	t726 = qJD(1) * t743;
	t725 = t680 * t735;
	t723 = t752 * t680;
	t722 = t681 * t730;
	t721 = t680 * t729;
	t720 = t680 * t726;
	t719 = t686 * t725;
	t718 = t690 * t725;
	t715 = qJD(3) * t723;
	t714 = -t687 * r_i_i_C(1) + t683 * r_i_i_C(2);
	t713 = -t683 * r_i_i_C(1) - t687 * r_i_i_C(2);
	t630 = t663 * t741 + t689 * t732 + t748;
	t712 = t630 * t688 + t760;
	t621 = t630 * t684 - t759;
	t619 = -t684 * t757 + t759;
	t634 = -t689 * t721 + t694 * t741 - t747;
	t656 = t680 * t694 + t681 * t728;
	t711 = t634 * t688 - t656 * t684;
	t623 = t634 * t684 + t656 * t688;
	t652 = -t689 * t723 - t722 + t731;
	t662 = t682 * t752 - t690 * t746;
	t710 = t652 * t688 - t662 * t684;
	t629 = t652 * t684 + t662 * t688;
	t708 = pkin(5) - t714;
	t706 = -t713 + t756;
	t643 = -t663 * t685 + t664 * t741;
	t624 = t643 * t684 + t664 * t744;
	t645 = -t685 * t694 - t695 * t741;
	t625 = t645 * t684 - t695 * t744;
	t644 = -t663 * t689 - t664 * t742;
	t646 = -t689 * t694 + t695 * t742;
	t703 = t682 * t738 + t739;
	t702 = t682 * t739 + t738;
	t701 = -t682 * t740 + t737;
	t700 = qJD(6) * t714;
	t699 = qJD(6) * t713;
	t659 = t703 * t681;
	t647 = t659 * t684 + t688 * t733;
	t698 = -t682 * t694 + t721;
	t635 = t685 * t698 - t689 * t695;
	t693 = t684 * t708 - t688 * t754 + qJ(4);
	t692 = qJD(4) + t684 * t699 + (t684 * t754 + t688 * t708) * qJD(5);
	t660 = t701 * t681;
	t653 = t681 * t702 + t685 * t723;
	t649 = qJD(1) * t664 + qJD(2) * t694;
	t639 = -qJD(2) * t722 - t681 * qJD(3) * t737 + (qJD(3) * t682 + qJD(2)) * t731;
	t636 = t682 * t726 - t751;
	t627 = t689 * t715 + ((t730 - t740) * qJD(3) + t701 * qJD(2)) * t681;
	t626 = t685 * t715 + (qJD(2) * t703 + qJD(3) * t702) * t681;
	t615 = qJD(3) * t644 - t650 * t685 + t651 * t741;
	t613 = qJD(3) * t646 + t648 * t685 - t649 * t741;
	t612 = qJD(5) * t710 + t626 * t684 + t688 * t719;
	t606 = qJD(3) * t747 + (t648 * t682 + t720) * t685 + (qJD(3) * t698 - t649) * t689;
	t605 = qJD(3) * t635 - t648 * t741 - t649 * t685 - t689 * t720;
	t599 = qJD(5) * t712 - t767;
	t595 = qJD(5) * t711 + t605 * t684 + t636 * t688;
	t594 = qJD(5) * t623 - t605 * t688 + t636 * t684;
	t593 = t595 * t687 + t606 * t683 + (-t623 * t683 + t635 * t687) * qJD(6);
	t592 = -t595 * t683 + t606 * t687 + (-t623 * t687 - t635 * t683) * qJD(6);
	t1 = [-pkin(10) * t750 - t651 * pkin(2) - t638 * pkin(4) + t609 * qJ(4) - t757 * qJD(4) + t708 * ((-t688 * t757 - t760) * qJD(5) + t767) + t706 * t765 + t754 * (qJD(5) * t619 + t766) + ((-t619 * t683 + t763) * r_i_i_C(1) + (-t619 * t687 - t764) * r_i_i_C(2)) * qJD(6) + (-t691 * pkin(1) + (-pkin(9) * t753 - pkin(10) * t728) * t681) * qJD(1), t648 * pkin(2) + t613 * qJ(4) + t645 * qJD(4) - t649 * t734 + t708 * (-t649 * t744 + t613 * t684 + (t645 * t688 + t695 * t745) * qJD(5)) + t706 * (-qJD(3) * t645 + t648 * t689 + t649 * t742) - t754 * (-qJD(5) * t625 + t613 * t688 + t649 * t745) + ((-t625 * t683 + t646 * t687) * r_i_i_C(1) + (-t625 * t687 - t646 * t683) * r_i_i_C(2)) * qJD(6), -t605 * t706 + t606 * t693 + t634 * t700 + t635 * t692, t605, -t594 * t708 + t595 * t754 + t699 * t711, t592 * r_i_i_C(1) - t593 * r_i_i_C(2); -pkin(10) * t751 - t649 * pkin(2) + t636 * pkin(4) + t595 * pkin(5) + t593 * r_i_i_C(1) + t592 * r_i_i_C(2) + t605 * qJ(4) + t634 * qJD(4) + t756 * t606 + t754 * t594 + (-t753 * pkin(1) + (pkin(10) * t682 + pkin(9)) * t743) * qJD(1), -t650 * pkin(2) + t615 * qJ(4) + t643 * qJD(4) + t651 * t734 + t708 * (t651 * t744 + t615 * t684 + (t643 * t688 - t664 * t745) * qJD(5)) + t706 * (-qJD(3) * t643 - t650 * t689 - t651 * t742) - t754 * (-qJD(5) * t624 + t615 * t688 - t651 * t745) + ((-t624 * t683 + t644 * t687) * r_i_i_C(1) + (-t624 * t687 - t644 * t683) * r_i_i_C(2)) * qJD(6), t609 * t706 + t630 * t700 - t633 * t692 - t693 * t765, -t609, t754 * t599 + t712 * t699 + t708 * (-qJD(5) * t621 + t766), (-t599 * t683 - t687 * t765) * r_i_i_C(1) + (-t599 * t687 + t683 * t765) * r_i_i_C(2) + ((-t621 * t687 + t764) * r_i_i_C(1) + (t621 * t683 + t763) * r_i_i_C(2)) * qJD(6); 0, -t639 * qJ(4) + t659 * qJD(4) + t708 * (t688 * t718 - t639 * t684 + (t659 * t688 - t684 * t733) * qJD(5)) - t706 * (qJD(2) * t702 + qJD(3) * t703) * t681 + t754 * (qJD(5) * t647 + t639 * t688 + t684 * t718) + ((-t647 * t683 + t660 * t687) * r_i_i_C(1) + (-t647 * t687 - t660 * t683) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t686 + t690 * t734) * t735, -t626 * t706 + t627 * t693 + t652 * t700 + t653 * t692, t626, t754 * t612 + t710 * t699 + t708 * (-qJD(5) * t629 + t626 * t688 - t684 * t719), (-t612 * t683 + t627 * t687) * r_i_i_C(1) + (-t612 * t687 - t627 * t683) * r_i_i_C(2) + ((-t629 * t687 - t653 * t683) * r_i_i_C(1) + (t629 * t683 - t653 * t687) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end