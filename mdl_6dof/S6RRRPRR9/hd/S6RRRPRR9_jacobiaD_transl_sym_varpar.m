% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:08
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR9_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR9_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR9_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
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
	% StartTime: 2019-10-10 12:08:08
	% EndTime: 2019-10-10 12:08:08
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
	% StartTime: 2019-10-10 12:08:09
	% EndTime: 2019-10-10 12:08:09
	% DurationCPUTime: 0.46s
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
	t1 = [-t290 * pkin(2) + t295 * r_i_i_C(1) - t289 * t342 + (-t309 * pkin(1) - t339 * t341) * qJD(1) + (t320 * r_i_i_C(1) - t312) * t307 + (-t314 * r_i_i_C(2) - t344) * t304, t319 * t287 - t345 * t288 + ((t304 * t317 + t316 * t336) * r_i_i_C(1) + (t307 * t317 - t316 * t337) * r_i_i_C(2)) * qJD(3), t285 * r_i_i_C(1) - t286 * r_i_i_C(2), 0, 0, 0; -t288 * pkin(2) + t286 * r_i_i_C(1) + t285 * r_i_i_C(2) - t287 * t342 + (-t306 * pkin(1) + t338 * t341) * qJD(1), -t319 * t289 + t345 * t290 + ((t291 * t304 - t292 * t336) * r_i_i_C(1) + (t291 * t307 + t292 * t337) * r_i_i_C(2)) * qJD(3), t295 * r_i_i_C(2) + (t320 * r_i_i_C(2) + t344) * t307 + (t314 * r_i_i_C(1) - t312) * t304, 0, 0, 0; 0, (t311 * qJD(3) + (-t305 * pkin(2) + t308 * t342 + t310) * qJD(2)) * t301, -t321 * t303 * t327 + (t311 * qJD(2) + t310 * qJD(3)) * t301, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:09
	% EndTime: 2019-10-10 12:08:10
	% DurationCPUTime: 0.70s
	% Computational Cost: add. (398->123), mult. (1234->207), div. (0->0), fcn. (1244->12), ass. (0->66)
	t349 = sin(pkin(13));
	t357 = cos(qJ(3));
	t390 = cos(pkin(13));
	t370 = qJD(3) * t390;
	t354 = sin(qJ(3));
	t379 = qJD(3) * t354;
	t393 = t349 * t379 - t357 * t370;
	t392 = t354 * pkin(3);
	t391 = pkin(10) + qJ(4);
	t350 = sin(pkin(7));
	t351 = sin(pkin(6));
	t389 = t350 * t351;
	t358 = cos(qJ(2));
	t388 = t354 * t358;
	t355 = sin(qJ(2));
	t387 = t355 * t357;
	t356 = sin(qJ(1));
	t386 = t356 * t355;
	t385 = t356 * t358;
	t359 = cos(qJ(1));
	t384 = t359 * t355;
	t383 = t359 * t358;
	t382 = qJD(1) * t356;
	t381 = qJD(1) * t359;
	t380 = qJD(2) * t355;
	t378 = qJD(3) * t357;
	t377 = pkin(3) * t379;
	t376 = pkin(3) * t378;
	t353 = cos(pkin(6));
	t375 = t353 * t386;
	t374 = t351 * t382;
	t373 = t351 * t381;
	t319 = t393 * t350;
	t352 = cos(pkin(7));
	t371 = -t319 * r_i_i_C(1) + t352 * qJD(4) + t350 * t376;
	t338 = -t357 * t349 - t354 * t390;
	t324 = t338 * t350;
	t369 = -r_i_i_C(1) * t324 + t350 * t392 + t391 * t352 + pkin(9);
	t331 = -t353 * t383 + t386;
	t367 = t353 * t385 + t384;
	t332 = t353 * t384 + t385;
	t366 = t375 - t383;
	t365 = -t354 * t349 + t357 * t390;
	t318 = -qJD(1) * t375 - t356 * t380 + (qJD(2) * t353 + qJD(1)) * t383;
	t321 = t393 * t352;
	t330 = -t349 * t378 - t354 * t370;
	t364 = -t318 * t365 - t331 * t321 - t332 * t330;
	t361 = qJD(3) * t338;
	t322 = t352 * t361;
	t329 = t365 * qJD(3);
	t363 = -t318 * t338 + t331 * t322 + t332 * t329;
	t325 = t365 * t352;
	t326 = t338 * t352;
	t328 = -t391 * t350 + t352 * t392;
	t362 = -t326 * r_i_i_C(1) + t325 * r_i_i_C(2) - t350 * r_i_i_C(3) + t328;
	t315 = t331 * qJD(1) + t366 * qJD(2);
	t316 = t332 * qJD(1) + t367 * qJD(2);
	t360 = -t315 * t326 - t316 * t365 + t321 * t367 - t330 * t366;
	t348 = t357 * pkin(3) + pkin(2);
	t336 = -t350 * qJD(4) + t352 * t376;
	t323 = t365 * t350;
	t320 = t350 * t361;
	t317 = t367 * qJD(1) + t332 * qJD(2);
	t314 = -t315 * t350 + t352 * t373;
	t313 = t315 * t325 - t316 * t338 - t367 * t322 + t366 * t329 + (t320 * t356 + t323 * t381) * t351;
	t1 = [t364 * r_i_i_C(1) + t363 * r_i_i_C(2) - t318 * t348 + t332 * t377 + t331 * t336 - pkin(1) * t381 + t362 * t317 + ((t320 * r_i_i_C(2) + t371) * t359 + (-r_i_i_C(2) * t323 - r_i_i_C(3) * t352 - t369) * t382) * t351, (t315 * t365 - t321 * t366 - t330 * t367) * r_i_i_C(1) + (t315 * t338 + t322 * t366 + t329 * t367) * r_i_i_C(2) + t315 * t348 + t367 * t377 + t366 * t336 + t362 * t316, t313 * r_i_i_C(1) + ((t319 * t356 + t324 * t381) * t351 - t360) * r_i_i_C(2) + (t316 * t354 + (t315 * t352 + t350 * t373) * t357 + (t366 * t357 + (t352 * t367 - t356 * t389) * t354) * qJD(3)) * pkin(3), t314, 0, 0; t360 * r_i_i_C(1) + t313 * r_i_i_C(2) + t314 * r_i_i_C(3) - t316 * t348 + t366 * t377 + t315 * t328 - t367 * t336 - pkin(1) * t382 + (t371 * t356 + t369 * t381) * t351, (-t317 * t365 + t332 * t321 - t331 * t330) * r_i_i_C(1) + (-t317 * t338 - t332 * t322 + t331 * t329) * r_i_i_C(2) - t317 * t348 + t331 * t377 - t332 * t336 - t362 * t318, (-t317 * t325 - t363) * r_i_i_C(1) + (-t317 * t326 + t364) * r_i_i_C(2) + ((-t320 * t359 + t323 * t382) * r_i_i_C(1) + (-t319 * t359 + t324 * t382) * r_i_i_C(2)) * t351 + (-t318 * t354 + (-t317 * t352 + t350 * t374) * t357 + (-t332 * t357 + (t331 * t352 + t359 * t389) * t354) * qJD(3)) * pkin(3), t317 * t350 + t352 * t374, 0, 0; 0, ((t321 * t355 + t330 * t358) * r_i_i_C(1) + (-t322 * t355 - t329 * t358) * r_i_i_C(2) - t358 * t377 - t355 * t336 + ((-r_i_i_C(1) * t365 - t338 * r_i_i_C(2) - t348) * t355 - t362 * t358) * qJD(2)) * t351, (t320 * r_i_i_C(1) + t319 * r_i_i_C(2) - t350 * t377) * t353 + ((t322 * t358 - t329 * t355) * r_i_i_C(1) + (t321 * t358 - t330 * t355) * r_i_i_C(2) + ((-t325 * t355 + t338 * t358) * r_i_i_C(1) + (-t326 * t355 - t358 * t365) * r_i_i_C(2)) * qJD(2) + ((-t352 * t388 - t387) * qJD(3) + (-t352 * t387 - t388) * qJD(2)) * pkin(3)) * t351, t380 * t389, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:13
	% EndTime: 2019-10-10 12:08:14
	% DurationCPUTime: 1.96s
	% Computational Cost: add. (1116->193), mult. (3439->332), div. (0->0), fcn. (3690->14), ass. (0->103)
	t532 = sin(qJ(2));
	t533 = sin(qJ(1));
	t535 = cos(qJ(2));
	t536 = cos(qJ(1));
	t585 = cos(pkin(6));
	t560 = t536 * t585;
	t509 = t532 * t560 + t533 * t535;
	t561 = t533 * t585;
	t543 = t532 * t536 + t535 * t561;
	t489 = qJD(1) * t543 + qJD(2) * t509;
	t556 = t532 * t561;
	t574 = qJD(2) * t532;
	t577 = t536 * t535;
	t490 = -qJD(1) * t556 - t533 * t574 + (qJD(2) * t585 + qJD(1)) * t577;
	t527 = sin(pkin(7));
	t526 = sin(pkin(13));
	t584 = cos(pkin(13));
	t589 = cos(qJ(3));
	t555 = t589 * t584;
	t531 = sin(qJ(3));
	t572 = qJD(3) * t531;
	t592 = -qJD(3) * t555 + t526 * t572;
	t495 = t592 * t527;
	t529 = cos(pkin(7));
	t497 = t592 * t529;
	t562 = t531 * t584;
	t542 = t526 * t589 + t562;
	t500 = t542 * t527;
	t502 = t542 * t529;
	t563 = t589 * qJD(3);
	t507 = -qJD(3) * t562 - t526 * t563;
	t508 = t532 * t533 - t535 * t560;
	t528 = sin(pkin(6));
	t541 = -t526 * t531 + t555;
	t576 = qJD(1) * t533;
	t462 = -t489 * t502 + t490 * t541 + t508 * t497 + t509 * t507 + t528 * (t495 * t536 + t500 * t576);
	t565 = t528 * t576;
	t485 = t489 * t527 + t529 * t565;
	t530 = sin(qJ(5));
	t534 = cos(qJ(5));
	t599 = t462 * t530 - t485 * t534;
	t598 = -t462 * t534 - t485 * t530;
	t551 = t502 * t532 - t535 * t541;
	t595 = t585 * t495 + (qJD(2) * t551 + t497 * t535 - t507 * t532) * t528;
	t579 = t528 * t536;
	t476 = t500 * t579 + t502 * t508 - t509 * t541;
	t568 = t529 * t579;
	t491 = -t508 * t527 + t568;
	t594 = t476 * t534 + t491 * t530;
	t593 = -t476 * t530 + t491 * t534;
	t544 = t556 - t577;
	t487 = qJD(1) * t508 + qJD(2) * t544;
	t488 = qJD(1) * t509 + qJD(2) * t543;
	t575 = qJD(1) * t536;
	t459 = t487 * t502 - t488 * t541 + t543 * t497 - t544 * t507 + t528 * (-t495 * t533 + t500 * t575);
	t590 = r_i_i_C(3) + pkin(11);
	t588 = pkin(3) * t531;
	t586 = pkin(10) + qJ(4);
	t587 = t527 * t588 + t529 * t586 + pkin(9);
	t583 = t527 * t528;
	t582 = t527 * t530;
	t581 = t527 * t534;
	t580 = t528 * t533;
	t578 = t531 * t535;
	t573 = qJD(2) * t535;
	t571 = qJD(5) * t530;
	t570 = qJD(5) * t534;
	t569 = pkin(3) * t572;
	t567 = t529 * t589;
	t566 = t589 * t532;
	t558 = pkin(3) * t563;
	t557 = t574 * t583;
	t501 = t541 * t529;
	t553 = t501 * t535 - t532 * t542;
	t552 = t502 * t535 + t532 * t541;
	t550 = qJD(1) * t589 * t583;
	t549 = r_i_i_C(1) * t534 - r_i_i_C(2) * t530 + pkin(4);
	t547 = qJD(5) * (-r_i_i_C(1) * t530 - r_i_i_C(2) * t534);
	t540 = qJD(3) * t542;
	t496 = t527 * t540;
	t498 = t529 * t540;
	t499 = t541 * t527;
	t506 = t541 * qJD(3);
	t537 = t489 * t501 + t490 * t542 - t496 * t579 - t498 * t508 - t499 * t565 + t506 * t509;
	t525 = pkin(3) * t589 + pkin(2);
	t513 = -qJD(4) * t527 + t529 * t558;
	t512 = qJD(4) * t529 + t527 * t558;
	t505 = t529 * t585 - t535 * t583;
	t504 = -t527 * t586 + t529 * t588;
	t493 = t527 * t543 + t529 * t580;
	t486 = t551 * t528;
	t483 = qJD(1) * t568 - t487 * t527;
	t482 = t502 * t544 - t541 * t543;
	t481 = -t502 * t509 - t508 * t541;
	t480 = t500 * t585 + t528 * t552;
	t478 = t500 * t580 - t502 * t543 - t541 * t544;
	t473 = (-qJD(2) * t552 + t497 * t532 + t507 * t535) * t528;
	t468 = -t489 * t541 - t490 * t502 + t497 * t509 - t507 * t508;
	t466 = t487 * t541 + t488 * t502 - t497 * t544 - t507 * t543;
	t458 = t487 * t501 + t488 * t542 + t543 * t498 + t544 * t506 + (-t496 * t533 + t499 * t575) * t528;
	t456 = t459 * t534 + t483 * t530 + (-t478 * t530 + t493 * t534) * qJD(5);
	t455 = -t459 * t530 + t483 * t534 + (-t478 * t534 - t493 * t530) * qJD(5);
	t1 = [t598 * r_i_i_C(1) + t599 * r_i_i_C(2) - t462 * pkin(4) - t490 * t525 + t509 * t569 + t489 * t504 + t508 * t513 + t512 * t579 - t590 * t537 + (r_i_i_C(1) * t593 - r_i_i_C(2) * t594) * qJD(5) + (-t536 * pkin(1) - t580 * t587) * qJD(1), (t466 * t534 - t488 * t582) * r_i_i_C(1) + (-t466 * t530 - t488 * t581) * r_i_i_C(2) + t466 * pkin(4) + t487 * t525 + t543 * t569 + t488 * t504 + t544 * t513 - t590 * (-t487 * t542 + t488 * t501 - t498 * t544 + t506 * t543) + ((-t482 * t530 - t544 * t581) * r_i_i_C(1) + (-t482 * t534 + t544 * t582) * r_i_i_C(2)) * qJD(5), t590 * t459 + (t499 * t580 - t501 * t543 + t542 * t544) * t547 + t549 * t458 + (t536 * t550 + t487 * t567 + t488 * t531 + (t589 * t544 + (-t527 * t580 + t529 * t543) * t531) * qJD(3)) * pkin(3), t483, r_i_i_C(1) * t455 - r_i_i_C(2) * t456, 0; t544 * t569 + t512 * t580 + t459 * pkin(4) + t456 * r_i_i_C(1) + t455 * r_i_i_C(2) + t487 * t504 - t488 * t525 - t543 * t513 - t590 * t458 + (-pkin(1) * t533 + t579 * t587) * qJD(1), (t468 * t534 + t490 * t582) * r_i_i_C(1) + (-t468 * t530 + t490 * t581) * r_i_i_C(2) + t468 * pkin(4) - t489 * t525 + t508 * t569 - t490 * t504 - t509 * t513 - t590 * (t489 * t542 - t490 * t501 + t498 * t509 + t506 * t508) + ((-t481 * t530 + t509 * t581) * r_i_i_C(1) + (-t481 * t534 - t509 * t582) * r_i_i_C(2)) * qJD(5), t590 * t462 + (-t499 * t579 - t501 * t508 - t509 * t542) * t547 - t549 * t537 + (t533 * t550 - t489 * t567 - t490 * t531 + (-t589 * t509 + (t508 * t529 + t527 * t579) * t531) * qJD(3)) * pkin(3), t485, -t599 * r_i_i_C(1) + t598 * r_i_i_C(2) + (r_i_i_C(1) * t594 + r_i_i_C(2) * t593) * qJD(5), 0; 0, (t473 * t534 + t486 * t571) * r_i_i_C(1) + (-t473 * t530 + t486 * t570) * r_i_i_C(2) + t473 * pkin(4) + (t590 * (qJD(2) * t553 - t498 * t532 + t506 * t535) - t535 * t569 - t532 * t513 + (-t504 * t535 - t525 * t532) * qJD(2) + ((t530 * t573 + t532 * t570) * r_i_i_C(1) + (-t532 * t571 + t534 * t573) * r_i_i_C(2)) * t527) * t528, -t590 * t595 + (t499 * t585 + t528 * t553) * t547 + t549 * (-t585 * t496 + (-t498 * t535 - t506 * t532 + (-t501 * t532 - t535 * t542) * qJD(2)) * t528) + (-t585 * t527 * t572 + ((-t529 * t578 - t566) * qJD(3) + (-t529 * t566 - t578) * qJD(2)) * t528) * pkin(3), t557, (t530 * t595 + t534 * t557) * r_i_i_C(1) + (-t530 * t557 + t534 * t595) * r_i_i_C(2) + ((-t480 * t534 - t505 * t530) * r_i_i_C(1) + (t480 * t530 - t505 * t534) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:18
	% EndTime: 2019-10-10 12:08:22
	% DurationCPUTime: 3.90s
	% Computational Cost: add. (2842->270), mult. (8652->451), div. (0->0), fcn. (9703->16), ass. (0->142)
	t712 = sin(pkin(7));
	t711 = sin(pkin(13));
	t783 = cos(pkin(13));
	t788 = cos(qJ(3));
	t750 = t788 * t783;
	t717 = sin(qJ(3));
	t771 = qJD(3) * t717;
	t790 = -qJD(3) * t750 + t711 * t771;
	t680 = t790 * t712;
	t758 = t717 * t783;
	t700 = -t788 * t711 - t758;
	t685 = t700 * t712;
	t713 = sin(pkin(6));
	t723 = cos(qJ(1));
	t718 = sin(qJ(2));
	t719 = sin(qJ(1));
	t722 = cos(qJ(2));
	t784 = cos(pkin(6));
	t756 = t723 * t784;
	t694 = t718 * t756 + t719 * t722;
	t757 = t719 * t784;
	t733 = t723 * t718 + t722 * t757;
	t671 = t733 * qJD(1) + t694 * qJD(2);
	t751 = t718 * t757;
	t772 = qJD(2) * t718;
	t775 = t723 * t722;
	t672 = -qJD(1) * t751 - t719 * t772 + (qJD(2) * t784 + qJD(1)) * t775;
	t714 = cos(pkin(7));
	t682 = t790 * t714;
	t687 = t700 * t714;
	t759 = qJD(3) * t788;
	t692 = -qJD(3) * t758 - t711 * t759;
	t693 = t719 * t718 - t722 * t756;
	t732 = -t717 * t711 + t750;
	t728 = -t671 * t687 - t672 * t732 - t693 * t682 - t694 * t692;
	t774 = qJD(1) * t719;
	t625 = (t680 * t723 - t685 * t774) * t713 - t728;
	t762 = t713 * t774;
	t666 = t671 * t712 + t714 * t762;
	t716 = sin(qJ(5));
	t721 = cos(qJ(5));
	t777 = t713 * t723;
	t652 = -t685 * t777 - t693 * t687 - t694 * t732;
	t766 = t714 * t777;
	t673 = -t693 * t712 + t766;
	t749 = t652 * t716 - t673 * t721;
	t611 = t749 * qJD(5) + t625 * t721 + t666 * t716;
	t730 = qJD(3) * t700;
	t681 = t712 * t730;
	t683 = t714 * t730;
	t684 = t732 * t712;
	t686 = t732 * t714;
	t691 = t732 * qJD(3);
	t624 = -t671 * t686 + t672 * t700 - t693 * t683 - t694 * t691 + (-t681 * t723 + t684 * t774) * t713;
	t715 = sin(qJ(6));
	t720 = cos(qJ(6));
	t803 = t611 * t715 + t624 * t720;
	t802 = -t611 * t720 + t624 * t715;
	t636 = t652 * t721 + t673 * t716;
	t740 = -t684 * t777 - t686 * t693 + t694 * t700;
	t801 = -t636 * t715 + t720 * t740;
	t800 = t636 * t720 + t715 * t740;
	t799 = t636 * qJD(5) - t625 * t716 + t666 * t721;
	t734 = t751 - t775;
	t669 = t693 * qJD(1) + t734 * qJD(2);
	t670 = t694 * qJD(1) + t733 * qJD(2);
	t773 = qJD(1) * t723;
	t622 = -t669 * t687 - t670 * t732 + t733 * t682 - t734 * t692 + (-t680 * t719 - t685 * t773) * t713;
	t789 = r_i_i_C(3) + pkin(12);
	t787 = pkin(3) * t717;
	t785 = pkin(10) + qJ(4);
	t786 = t712 * t787 + t785 * t714 + pkin(9);
	t782 = t692 * t718;
	t781 = t712 * t713;
	t780 = t712 * t716;
	t779 = t712 * t721;
	t778 = t713 * t719;
	t776 = t717 * t722;
	t770 = qJD(6) * t715;
	t769 = qJD(6) * t720;
	t768 = pkin(3) * t771;
	t767 = t718 * t781;
	t765 = t722 * t781;
	t764 = t714 * t788;
	t763 = t788 * t718;
	t761 = t713 * t772;
	t755 = t784 * t680;
	t754 = pkin(3) * t759;
	t753 = qJD(2) * t765;
	t752 = t712 * t761;
	t675 = t712 * t733 + t714 * t778;
	t731 = -t685 * t778 + t687 * t733 - t732 * t734;
	t638 = t675 * t716 + t721 * t731;
	t748 = t675 * t721 - t716 * t731;
	t690 = t784 * t714 - t765;
	t745 = -t687 * t722 + t718 * t732;
	t726 = -t784 * t685 + t745 * t713;
	t645 = t690 * t716 + t721 * t726;
	t747 = t690 * t721 - t716 * t726;
	t746 = t686 * t722 + t700 * t718;
	t744 = t687 * t718 + t722 * t732;
	t743 = qJD(1) * t788 * t781;
	t742 = r_i_i_C(1) * t720 - r_i_i_C(2) * t715 + pkin(5);
	t660 = t687 * t694 - t693 * t732;
	t646 = t660 * t721 + t694 * t780;
	t662 = -t687 * t734 - t732 * t733;
	t647 = t662 * t721 - t734 * t780;
	t739 = qJD(6) * (-t715 * r_i_i_C(1) - t720 * r_i_i_C(2));
	t668 = t744 * t713;
	t663 = t668 * t721 + t716 * t767;
	t735 = qJD(1) * t766 - t669 * t712;
	t725 = t789 * t716 + t742 * t721 + pkin(4);
	t724 = t721 * t739 + (-t742 * t716 + t789 * t721) * qJD(5);
	t710 = t788 * pkin(3) + pkin(2);
	t698 = -t712 * qJD(4) + t714 * t754;
	t697 = t714 * qJD(4) + t712 * t754;
	t689 = -t785 * t712 + t714 * t787;
	t667 = (t686 * t718 - t700 * t722) * t713;
	t661 = -t686 * t734 + t700 * t733;
	t659 = t686 * t694 + t693 * t700;
	t657 = t784 * t684 + t746 * t713;
	t654 = t684 * t778 - t686 * t733 - t700 * t734;
	t643 = (-t745 * qJD(2) + t682 * t718 + t692 * t722) * t713;
	t642 = (t746 * qJD(2) + t683 * t718 + t691 * t722) * t713;
	t641 = -t755 + (t744 * qJD(2) - t682 * t722 + t782) * t713;
	t640 = t784 * t681 - t686 * t761 + (-t691 * t718 + (qJD(2) * t700 + t683) * t722) * t713;
	t639 = t755 - t687 * t761 + (-t782 + (-qJD(2) * t732 + t682) * t722) * t713;
	t633 = -t671 * t732 + t672 * t687 + t682 * t694 - t692 * t693;
	t632 = -t671 * t700 - t672 * t686 - t683 * t694 + t691 * t693;
	t631 = t669 * t732 - t670 * t687 - t682 * t734 - t692 * t733;
	t630 = t669 * t700 + t670 * t686 + t683 * t734 + t691 * t733;
	t629 = t716 * t753 + t643 * t721 + (-t668 * t716 + t721 * t767) * qJD(5);
	t623 = -t680 * t777 + t685 * t762 + t728;
	t621 = t669 * t686 - t670 * t700 - t733 * t683 + t734 * t691 + (t681 * t719 + t684 * t773) * t713;
	t619 = t747 * qJD(5) + t641 * t721 + t716 * t752;
	t617 = t672 * t780 + t633 * t721 + (-t660 * t716 + t694 * t779) * qJD(5);
	t615 = -t670 * t780 + t631 * t721 + (-t662 * t716 - t734 * t779) * qJD(5);
	t609 = t748 * qJD(5) + t622 * t721 + t735 * t716;
	t608 = t638 * qJD(5) + t622 * t716 - t735 * t721;
	t607 = t609 * t720 - t621 * t715 + (-t638 * t715 - t654 * t720) * qJD(6);
	t606 = -t609 * t715 - t621 * t720 + (-t638 * t720 + t654 * t715) * qJD(6);
	t1 = [t802 * r_i_i_C(1) + t803 * r_i_i_C(2) - t611 * pkin(5) - t625 * pkin(4) + t624 * pkin(11) - t672 * t710 + t694 * t768 + t671 * t689 + t693 * t698 + t697 * t777 + t789 * t799 + (t801 * r_i_i_C(1) - t800 * r_i_i_C(2)) * qJD(6) + (-t723 * pkin(1) - t786 * t778) * qJD(1), (t615 * t720 - t630 * t715) * r_i_i_C(1) + (-t615 * t715 - t630 * t720) * r_i_i_C(2) + t615 * pkin(5) + t631 * pkin(4) - t630 * pkin(11) + t669 * t710 + t733 * t768 + t670 * t689 + t734 * t698 + t789 * (t647 * qJD(5) + t631 * t716 + t670 * t779) + ((-t647 * t715 + t661 * t720) * r_i_i_C(1) + (-t647 * t720 - t661 * t715) * r_i_i_C(2)) * qJD(6), (t622 * t715 + t731 * t769) * r_i_i_C(1) + (t622 * t720 - t731 * t770) * r_i_i_C(2) + t622 * pkin(11) + t725 * t621 + t724 * t654 + (t723 * t743 + t669 * t764 + t670 * t717 + (t788 * t734 + (-t712 * t778 + t714 * t733) * t717) * qJD(3)) * pkin(3), t735, -t742 * t608 + t789 * t609 + t748 * t739, r_i_i_C(1) * t606 - r_i_i_C(2) * t607; t734 * t768 + t697 * t778 + t622 * pkin(4) + t609 * pkin(5) - t621 * pkin(11) + t607 * r_i_i_C(1) + t606 * r_i_i_C(2) + t669 * t689 - t670 * t710 - t733 * t698 + t789 * t608 + (-t719 * pkin(1) + t786 * t777) * qJD(1), (t617 * t720 - t632 * t715) * r_i_i_C(1) + (-t617 * t715 - t632 * t720) * r_i_i_C(2) + t617 * pkin(5) + t633 * pkin(4) - t632 * pkin(11) - t671 * t710 + t693 * t768 - t672 * t689 - t694 * t698 + t789 * (t646 * qJD(5) + t633 * t716 - t672 * t779) + ((-t646 * t715 + t659 * t720) * r_i_i_C(1) + (-t646 * t720 - t659 * t715) * r_i_i_C(2)) * qJD(6), (-t623 * t715 - t652 * t769) * r_i_i_C(1) + (-t623 * t720 + t652 * t770) * r_i_i_C(2) - t623 * pkin(11) + t725 * t624 + t724 * t740 + (t719 * t743 - t671 * t764 - t672 * t717 + (-t788 * t694 + (t693 * t714 + t712 * t777) * t717) * qJD(3)) * pkin(3), t666, t789 * t611 + t749 * t739 + t742 * t799, -t803 * r_i_i_C(1) + t802 * r_i_i_C(2) + (t800 * r_i_i_C(1) + t801 * r_i_i_C(2)) * qJD(6); 0, (t629 * t720 + t642 * t715) * r_i_i_C(1) + (-t629 * t715 + t642 * t720) * r_i_i_C(2) + t629 * pkin(5) + t643 * pkin(4) + t642 * pkin(11) + t789 * (t663 * qJD(5) + t643 * t716 - t721 * t753) + ((-t663 * t715 + t667 * t720) * r_i_i_C(1) + (-t663 * t720 - t667 * t715) * r_i_i_C(2)) * qJD(6) + (-t722 * t768 - t718 * t698 + (-t689 * t722 - t710 * t718) * qJD(2)) * t713, (-t639 * t715 + t726 * t769) * r_i_i_C(1) + (-t639 * t720 - t726 * t770) * r_i_i_C(2) - t639 * pkin(11) + t725 * t640 + t724 * t657 + (-t784 * t712 * t771 + ((-t714 * t776 - t763) * qJD(3) + (-t714 * t763 - t776) * qJD(2)) * t713) * pkin(3), t752, t789 * t619 + t747 * t739 + t742 * (-t645 * qJD(5) - t641 * t716 + t721 * t752), (-t619 * t715 - t640 * t720) * r_i_i_C(1) + (-t619 * t720 + t640 * t715) * r_i_i_C(2) + ((-t645 * t720 + t657 * t715) * r_i_i_C(1) + (t645 * t715 + t657 * t720) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end