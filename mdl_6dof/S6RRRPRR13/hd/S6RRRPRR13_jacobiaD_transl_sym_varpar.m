% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR13
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
% Datum: 2019-10-10 12:14
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR13_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR13_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR13_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
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
	% StartTime: 2019-10-10 12:14:18
	% EndTime: 2019-10-10 12:14:18
	% DurationCPUTime: 0.07s
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
	% StartTime: 2019-10-10 12:14:19
	% EndTime: 2019-10-10 12:14:20
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
	t1 = [-t290 * pkin(2) + t295 * r_i_i_C(1) - t289 * t342 + (-t309 * pkin(1) - t341 * t339) * qJD(1) + (t320 * r_i_i_C(1) - t312) * t307 + (-t314 * r_i_i_C(2) - t344) * t304, t319 * t287 - t345 * t288 + ((t304 * t317 + t316 * t336) * r_i_i_C(1) + (t307 * t317 - t316 * t337) * r_i_i_C(2)) * qJD(3), t285 * r_i_i_C(1) - t286 * r_i_i_C(2), 0, 0, 0; -t288 * pkin(2) + t286 * r_i_i_C(1) + t285 * r_i_i_C(2) - t287 * t342 + (-t306 * pkin(1) + t341 * t338) * qJD(1), -t319 * t289 + t345 * t290 + ((t291 * t304 - t292 * t336) * r_i_i_C(1) + (t291 * t307 + t292 * t337) * r_i_i_C(2)) * qJD(3), t295 * r_i_i_C(2) + (t320 * r_i_i_C(2) + t344) * t307 + (t314 * r_i_i_C(1) - t312) * t304, 0, 0, 0; 0, (t311 * qJD(3) + (-t305 * pkin(2) + t308 * t342 + t310) * qJD(2)) * t301, -t321 * t303 * t327 + (t311 * qJD(2) + t310 * qJD(3)) * t301, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:21
	% EndTime: 2019-10-10 12:14:22
	% DurationCPUTime: 0.81s
	% Computational Cost: add. (578->109), mult. (1811->183), div. (0->0), fcn. (1844->12), ass. (0->72)
	t445 = sin(qJ(2));
	t446 = sin(qJ(1));
	t496 = cos(pkin(6));
	t473 = t446 * t496;
	t467 = t445 * t473;
	t448 = cos(qJ(2));
	t449 = cos(qJ(1));
	t479 = t449 * t448;
	t481 = t446 * t445;
	t426 = -qJD(1) * t467 - qJD(2) * t481 + (qJD(2) * t496 + qJD(1)) * t479;
	t440 = sin(pkin(7));
	t441 = sin(pkin(6));
	t487 = t441 * t449;
	t477 = t440 * t487;
	t505 = -qJD(3) * t477 + t426;
	t497 = r_i_i_C(3) + qJ(4);
	t439 = sin(pkin(13));
	t442 = cos(pkin(13));
	t504 = -t442 * r_i_i_C(1) + t439 * r_i_i_C(2) - pkin(3);
	t444 = sin(qJ(3));
	t447 = cos(qJ(3));
	t454 = t467 - t479;
	t443 = cos(pkin(7));
	t453 = t449 * t445 + t448 * t473;
	t488 = t441 * t446;
	t459 = t440 * t488 - t443 * t453;
	t503 = t459 * t444 - t447 * t454;
	t472 = t449 * t496;
	t428 = t445 * t472 + t446 * t448;
	t425 = t453 * qJD(1) + t428 * qJD(2);
	t427 = -t448 * t472 + t481;
	t478 = qJD(1) * t441;
	t475 = t446 * t478;
	t470 = t440 * t475;
	t485 = t443 * t447;
	t486 = t443 * t444;
	t492 = t428 * t447;
	t502 = (t427 * t486 - t492) * qJD(3) - t425 * t485 + t447 * t470 - t505 * t444;
	t493 = t427 * t443;
	t501 = (-qJD(3) * t428 - t425 * t443 + t470) * t444 + (-qJD(3) * t493 + t505) * t447;
	t499 = t440 * pkin(10);
	t423 = t427 * qJD(1) + t454 * qJD(2);
	t495 = t423 * t440;
	t494 = t425 * t440;
	t490 = t439 * t440;
	t489 = t440 * t442;
	t484 = t444 * t445;
	t483 = t444 * t448;
	t482 = t445 * t447;
	t480 = t447 * t448;
	t476 = pkin(10) * t443 + pkin(9);
	t474 = t449 * t478;
	t471 = t496 * t440;
	t469 = t440 * t474;
	t466 = qJD(3) * t471;
	t463 = t427 * t444 - t428 * t485;
	t461 = t477 + t493;
	t460 = t444 * t453 + t454 * t485;
	t458 = t443 * t480 - t484;
	t457 = t443 * t482 + t483;
	t456 = t443 * t483 + t482;
	t455 = t443 * t484 - t480;
	t450 = t444 * t454 + t459 * t447;
	t424 = t428 * qJD(1) + t453 * qJD(2);
	t418 = -t443 * t475 - t494;
	t417 = t443 * t474 - t495;
	t415 = t444 * t466 + (t457 * qJD(2) + t456 * qJD(3)) * t441;
	t414 = t463 * qJD(3) - t425 * t447 - t426 * t486;
	t412 = t460 * qJD(3) + t423 * t447 + t424 * t486;
	t406 = -t424 * t447 + (t423 * t443 + t469) * t444 + t450 * qJD(3);
	t405 = t503 * qJD(3) - t423 * t485 - t424 * t444 - t447 * t469;
	t1 = [(t418 * t439 - t442 * t501) * r_i_i_C(1) + (t418 * t442 + t439 * t501) * r_i_i_C(2) - t501 * pkin(3) - (t428 * t444 + t461 * t447) * qJD(4) - t426 * pkin(2) - pkin(10) * t494 + t497 * t502 + (-t449 * pkin(1) - t476 * t488) * qJD(1), (t412 * t442 - t424 * t490) * r_i_i_C(1) + (-t412 * t439 - t424 * t489) * r_i_i_C(2) + t412 * pkin(3) - t460 * qJD(4) + t423 * pkin(2) - t424 * t499 + t497 * (-t424 * t485 + t423 * t444 + (-t447 * t453 + t454 * t486) * qJD(3)), t503 * qJD(4) + t405 * t504 + t497 * t406, t405, 0, 0; (t406 * t442 + t417 * t439) * r_i_i_C(1) + (-t406 * t439 + t417 * t442) * r_i_i_C(2) + t406 * pkin(3) - t450 * qJD(4) - t424 * pkin(2) - pkin(10) * t495 + t497 * t405 + (-t446 * pkin(1) + t476 * t487) * qJD(1), (t414 * t442 + t426 * t490) * r_i_i_C(1) + (-t414 * t439 + t426 * t489) * r_i_i_C(2) + t414 * pkin(3) - t463 * qJD(4) - t425 * pkin(2) + t426 * t499 + t497 * (t426 * t485 - t425 * t444 + (-t427 * t447 - t428 * t486) * qJD(3)), -(t461 * t444 - t492) * qJD(4) + t497 * t501 - t504 * t502, -t502, 0, 0; 0, (-t497 * (-t458 * qJD(2) + t455 * qJD(3)) - t504 * (-t456 * qJD(2) - t457 * qJD(3)) + t457 * qJD(4) + (-t445 * pkin(2) + (r_i_i_C(1) * t439 + r_i_i_C(2) * t442 + pkin(10)) * t448 * t440) * qJD(2)) * t441, -(-t456 * t441 - t444 * t471) * qJD(4) + t497 * (t447 * t466 + (-t455 * qJD(2) + t458 * qJD(3)) * t441) + t504 * t415, t415, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:21
	% EndTime: 2019-10-10 12:14:23
	% DurationCPUTime: 1.29s
	% Computational Cost: add. (904->146), mult. (2538->248), div. (0->0), fcn. (2654->14), ass. (0->97)
	t521 = sin(qJ(1));
	t520 = sin(qJ(2));
	t576 = cos(pkin(6));
	t553 = t521 * t576;
	t546 = t520 * t553;
	t560 = qJD(2) * t520;
	t523 = cos(qJ(2));
	t524 = cos(qJ(1));
	t562 = t524 * t523;
	t489 = -qJD(1) * t546 - t521 * t560 + (qJD(2) * t576 + qJD(1)) * t562;
	t515 = sin(pkin(7));
	t516 = sin(pkin(6));
	t569 = t516 * t524;
	t557 = t515 * t569;
	t587 = -qJD(3) * t557 + t489;
	t552 = t524 * t576;
	t498 = t520 * t552 + t521 * t523;
	t528 = t524 * t520 + t523 * t553;
	t488 = t528 * qJD(1) + t498 * qJD(2);
	t517 = cos(pkin(7));
	t519 = sin(qJ(3));
	t522 = cos(qJ(3));
	t561 = qJD(1) * t516;
	t555 = t521 * t561;
	t550 = t515 * t555;
	t497 = t521 * t520 - t523 * t552;
	t573 = t497 * t517;
	t463 = (-qJD(3) * t498 - t488 * t517 + t550) * t519 + (-qJD(3) * t573 + t587) * t522;
	t574 = t488 * t515;
	t479 = t517 * t555 + t574;
	t513 = pkin(13) + qJ(5);
	t511 = sin(t513);
	t512 = cos(t513);
	t586 = t463 * t511 - t479 * t512;
	t585 = -t463 * t512 - t479 * t511;
	t568 = t517 * t519;
	t539 = t497 * t568 - t498 * t522;
	t474 = t519 * t557 + t539;
	t492 = -t497 * t515 + t517 * t569;
	t582 = t474 * t512 + t492 * t511;
	t581 = -t474 * t511 + t492 * t512;
	t580 = (t557 + t573) * t522 + t498 * t519;
	t567 = t517 * t522;
	t579 = t539 * qJD(3) - t488 * t567 - t587 * t519 + t522 * t550;
	t578 = sin(pkin(13)) * pkin(4);
	t577 = r_i_i_C(3) + pkin(11) + qJ(4);
	t529 = t546 - t562;
	t486 = t497 * qJD(1) + t529 * qJD(2);
	t575 = t486 * t515;
	t571 = t515 * t516;
	t570 = t516 * t521;
	t566 = t519 * t520;
	t565 = t519 * t523;
	t564 = t520 * t522;
	t563 = t522 * t523;
	t559 = qJD(5) * t511;
	t558 = qJD(5) * t512;
	t556 = pkin(10) * t517 + pkin(9);
	t554 = t524 * t561;
	t551 = t576 * t515;
	t549 = t515 * t554;
	t548 = t560 * t571;
	t545 = qJD(3) * t551;
	t544 = t512 * r_i_i_C(1) - t511 * r_i_i_C(2);
	t543 = -t511 * r_i_i_C(1) - t512 * r_i_i_C(2);
	t510 = cos(pkin(13)) * pkin(4) + pkin(3);
	t541 = -t510 - t544;
	t540 = t497 * t519 - t498 * t567;
	t484 = -t497 * t522 - t498 * t568;
	t537 = t519 * t528 + t529 * t567;
	t485 = -t522 * t528 + t529 * t568;
	t536 = t515 * t570 - t517 * t528;
	t535 = t517 * t563 - t566;
	t534 = t517 * t564 + t565;
	t533 = t517 * t565 + t564;
	t532 = t517 * t566 - t563;
	t531 = qJD(5) * t544;
	t530 = qJD(5) * t543;
	t527 = pkin(10) - t543 + t578;
	t475 = t519 * t529 + t536 * t522;
	t476 = t536 * t519 - t522 * t529;
	t496 = t576 * t517 - t523 * t571;
	t495 = t532 * t516;
	t494 = t515 * t528 + t517 * t570;
	t491 = t533 * t516 + t519 * t551;
	t487 = t498 * qJD(1) + t528 * qJD(2);
	t481 = (-t533 * qJD(2) - t534 * qJD(3)) * t516;
	t477 = t517 * t554 - t575;
	t471 = t522 * t545 + (-t532 * qJD(2) + t535 * qJD(3)) * t516;
	t470 = t519 * t545 + (t534 * qJD(2) + t533 * qJD(3)) * t516;
	t469 = t540 * qJD(3) - t488 * t522 - t489 * t568;
	t467 = t537 * qJD(3) + t486 * t522 + t487 * t568;
	t461 = -t487 * t522 + (t486 * t517 + t549) * t519 + t475 * qJD(3);
	t460 = t476 * qJD(3) - t486 * t567 - t487 * t519 - t522 * t549;
	t459 = t461 * t512 + t477 * t511 + (-t476 * t511 + t494 * t512) * qJD(5);
	t458 = -t461 * t511 + t477 * t512 + (-t476 * t512 - t494 * t511) * qJD(5);
	t1 = [t585 * r_i_i_C(1) + t586 * r_i_i_C(2) - t463 * t510 - t580 * qJD(4) - t479 * t578 - t489 * pkin(2) - pkin(10) * t574 + t577 * t579 + (t581 * r_i_i_C(1) - t582 * r_i_i_C(2)) * qJD(5) + (-t524 * pkin(1) - t556 * t570) * qJD(1), (t467 * t512 - t485 * t559) * r_i_i_C(1) + (-t467 * t511 - t485 * t558) * r_i_i_C(2) + t467 * t510 - t537 * qJD(4) + t486 * pkin(2) + t577 * (t485 * qJD(3) + t486 * t519 - t487 * t567) + (-t527 * t487 - t529 * t531) * t515, t476 * qJD(4) + t541 * t460 + t577 * t461 + t475 * t530, t460, t458 * r_i_i_C(1) - t459 * r_i_i_C(2), 0; t477 * t578 - pkin(10) * t575 - t487 * pkin(2) + t459 * r_i_i_C(1) + t458 * r_i_i_C(2) - t475 * qJD(4) + t461 * t510 + t577 * t460 + (-pkin(1) * t521 + t556 * t569) * qJD(1), (t469 * t512 - t484 * t559) * r_i_i_C(1) + (-t469 * t511 - t484 * t558) * r_i_i_C(2) + t469 * t510 - t540 * qJD(4) - t488 * pkin(2) + t577 * (t484 * qJD(3) - t488 * t519 + t489 * t567) + (t527 * t489 + t498 * t531) * t515, -qJD(4) * t474 + t577 * t463 - t580 * t530 - t541 * t579, -t579, -t586 * r_i_i_C(1) + t585 * r_i_i_C(2) + (t582 * r_i_i_C(1) + t581 * r_i_i_C(2)) * qJD(5), 0; 0, (t481 * t512 + t495 * t559) * r_i_i_C(1) + (-t481 * t511 + t495 * t558) * r_i_i_C(2) + t481 * t510 + (-t577 * (-t535 * qJD(2) + t532 * qJD(3)) + t534 * qJD(4) - pkin(2) * t560 + (t527 * t523 * qJD(2) + t520 * t531) * t515) * t516, t491 * qJD(4) + t577 * t471 + (t535 * t516 + t522 * t551) * t530 + t541 * t470, t470, (-t471 * t511 + t512 * t548) * r_i_i_C(1) + (-t471 * t512 - t511 * t548) * r_i_i_C(2) + ((-t491 * t512 - t496 * t511) * r_i_i_C(1) + (t491 * t511 - t496 * t512) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:28
	% EndTime: 2019-10-10 12:14:32
	% DurationCPUTime: 3.39s
	% Computational Cost: add. (2289->229), mult. (6212->378), div. (0->0), fcn. (6746->16), ass. (0->133)
	t696 = sin(qJ(2));
	t698 = cos(qJ(2));
	t699 = cos(qJ(1));
	t768 = cos(pkin(6));
	t737 = t699 * t768;
	t770 = sin(qJ(1));
	t668 = t696 * t737 + t770 * t698;
	t728 = t768 * t770;
	t710 = t699 * t696 + t698 * t728;
	t655 = t710 * qJD(1) + t668 * qJD(2);
	t720 = t696 * t728;
	t743 = t770 * t696;
	t754 = t699 * t698;
	t656 = -qJD(1) * t720 - qJD(2) * t743 + (qJD(2) * t768 + qJD(1)) * t754;
	t667 = -t698 * t737 + t743;
	t690 = sin(pkin(7));
	t692 = cos(pkin(7));
	t695 = sin(qJ(3));
	t691 = sin(pkin(6));
	t755 = t699 * t691;
	t748 = t690 * t755;
	t771 = cos(qJ(3));
	t725 = t771 * t748;
	t740 = qJD(3) * t771;
	t731 = t692 * t740;
	t739 = t770 * qJD(1);
	t732 = t691 * t739;
	t614 = (-qJD(3) * t668 - t655 * t692 + t690 * t732) * t695 - qJD(3) * t725 + t656 * t771 - t667 * t731;
	t765 = t655 * t690;
	t644 = t692 * t732 + t765;
	t688 = pkin(13) + qJ(5);
	t686 = sin(t688);
	t687 = cos(t688);
	t659 = -t667 * t690 + t692 * t755;
	t678 = t695 * t748;
	t758 = t692 * t695;
	t783 = t667 * t758 - t668 * t771 + t678;
	t776 = t659 * t687 - t686 * t783;
	t606 = t776 * qJD(5) - t614 * t687 - t644 * t686;
	t760 = t690 * t691;
	t735 = t770 * t760;
	t718 = t771 * t735;
	t747 = t692 * t771;
	t615 = qJD(1) * t718 + t783 * qJD(3) - t655 * t747 - t656 * t695;
	t694 = sin(qJ(6));
	t697 = cos(qJ(6));
	t790 = t606 * t694 - t615 * t697;
	t789 = t606 * t697 + t615 * t694;
	t626 = -t659 * t686 - t687 * t783;
	t637 = t667 * t747 + t668 * t695 + t725;
	t788 = t626 * t694 - t637 * t697;
	t787 = t626 * t697 + t637 * t694;
	t786 = -t626 * qJD(5) - t614 * t686 + t644 * t687;
	t769 = pkin(4) * sin(pkin(13));
	t729 = (-t769 - pkin(10)) * t690;
	t772 = r_i_i_C(3) + pkin(12);
	t711 = t720 - t754;
	t774 = t695 * t711 + t718;
	t642 = -t711 * t771 + (-t692 * t710 + t735) * t695;
	t746 = t692 * t770;
	t661 = t690 * t710 + t691 * t746;
	t773 = -t642 * t686 + t661 * t687;
	t719 = qJD(6) * (t694 * r_i_i_C(1) + t697 * r_i_i_C(2));
	t762 = t686 * t690;
	t761 = t687 * t690;
	t759 = t691 * t698;
	t757 = t695 * t696;
	t756 = t695 * t698;
	t753 = qJD(2) * t691;
	t752 = qJD(6) * t694;
	t751 = qJD(6) * t697;
	t750 = t696 * t760;
	t749 = t691 * t757;
	t745 = t771 * t696;
	t744 = t771 * t698;
	t742 = qJD(1) * t755;
	t741 = t690 * t753;
	t738 = t690 * t768;
	t736 = t692 * t744;
	t734 = t696 * t741;
	t733 = t698 * t741;
	t730 = t695 * t738;
	t726 = t691 * t736;
	t629 = t642 * t687 + t661 * t686;
	t713 = t692 * t756 + t745;
	t658 = t713 * t691 + t730;
	t666 = -t690 * t759 + t768 * t692;
	t633 = t658 * t687 + t666 * t686;
	t723 = -t658 * t686 + t666 * t687;
	t722 = r_i_i_C(1) * t697 - r_i_i_C(2) * t694 + pkin(5);
	t721 = t771 * t738;
	t651 = -t667 * t771 - t668 * t758;
	t630 = t651 * t687 + t668 * t762;
	t653 = -t710 * t771 + t711 * t758;
	t631 = t653 * t687 - t711 * t762;
	t714 = -t692 * t757 + t744;
	t665 = t714 * t691;
	t649 = t665 * t687 + t686 * t750;
	t716 = t667 * t695 - t668 * t747;
	t715 = t695 * t710 + t711 * t747;
	t712 = t692 * t745 + t756;
	t708 = t711 * qJD(2);
	t685 = cos(pkin(13)) * pkin(4) + pkin(3);
	t705 = -t772 * t686 - t722 * t687 - t685;
	t704 = t687 * t719 + (t722 * t686 - t772 * t687) * qJD(5);
	t703 = t667 * qJD(1) + t708;
	t702 = t703 * t695;
	t701 = t703 * t771;
	t700 = t659 * qJD(1) - t690 * t708;
	t693 = -pkin(11) - qJ(4);
	t664 = t712 * t691;
	t657 = -t721 - t726 + t749;
	t654 = t668 * qJD(1) + t710 * qJD(2);
	t646 = (-t713 * qJD(2) - t712 * qJD(3)) * t691;
	t645 = -qJD(2) * t726 - t740 * t759 + (qJD(3) * t692 + qJD(2)) * t749;
	t641 = t710 * t747 - t774;
	t635 = qJD(3) * t721 + ((t736 - t757) * qJD(3) + t714 * qJD(2)) * t691;
	t634 = qJD(3) * t730 + (t712 * qJD(2) + t713 * qJD(3)) * t691;
	t624 = t686 * t733 + t646 * t687 + (-t665 * t686 + t687 * t750) * qJD(5);
	t622 = t716 * qJD(3) - t655 * t771 - t656 * t758;
	t621 = t651 * qJD(3) - t655 * t695 + t656 * t747;
	t620 = t715 * qJD(3) + t654 * t758 + t701;
	t619 = t653 * qJD(3) - t654 * t747 + t702;
	t618 = t723 * qJD(5) + t635 * t687 + t686 * t734;
	t612 = qJD(1) * t678 + t774 * qJD(3) - t654 * t771 + t692 * t702 - t710 * t731;
	t611 = -qJD(1) * t725 + t642 * qJD(3) - t654 * t695 - t692 * t701;
	t610 = t656 * t762 + t622 * t687 + (-t651 * t686 + t668 * t761) * qJD(5);
	t608 = -t654 * t762 + t620 * t687 + (-t653 * t686 - t711 * t761) * qJD(5);
	t602 = t773 * qJD(5) + t612 * t687 + t700 * t686;
	t601 = t629 * qJD(5) + t612 * t686 - t700 * t687;
	t600 = t602 * t697 + t611 * t694 + (-t629 * t694 + t641 * t697) * qJD(6);
	t599 = -t602 * t694 + t611 * t697 + (-t629 * t697 - t641 * t694) * qJD(6);
	t1 = [t789 * r_i_i_C(1) - t790 * r_i_i_C(2) + t606 * pkin(5) - t614 * t685 - t615 * t693 - t637 * qJD(4) - t644 * t769 - t656 * pkin(2) - pkin(10) * t765 + t772 * t786 + (t788 * r_i_i_C(1) + t787 * r_i_i_C(2)) * qJD(6) + (-t699 * pkin(1) + (-t770 * pkin(9) - pkin(10) * t746) * t691) * qJD(1), (t608 * t697 + t619 * t694 + (-t631 * t694 - t697 * t715) * qJD(6)) * r_i_i_C(1) + (-t608 * t694 + t619 * t697 + (-t631 * t697 + t694 * t715) * qJD(6)) * r_i_i_C(2) + t608 * pkin(5) + t620 * t685 - t619 * t693 - t715 * qJD(4) + t703 * pkin(2) + t654 * t729 + t772 * (t631 * qJD(5) + t620 * t686 + t654 * t761), (t612 * t694 + t642 * t751) * r_i_i_C(1) + (t612 * t697 - t642 * t752) * r_i_i_C(2) - t612 * t693 + t642 * qJD(4) + t705 * t611 + t704 * t641, t611, -t722 * t601 + t772 * t602 - t773 * t719, r_i_i_C(1) * t599 - r_i_i_C(2) * t600; -pkin(1) * t739 - t654 * pkin(2) + t602 * pkin(5) + pkin(9) * t742 + t600 * r_i_i_C(1) + t599 * r_i_i_C(2) + t641 * qJD(4) - t611 * t693 + t612 * t685 + t700 * t769 + t772 * t601 + (-t690 * t703 + t692 * t742) * pkin(10), (t610 * t697 + t621 * t694) * r_i_i_C(1) + (-t610 * t694 + t621 * t697) * r_i_i_C(2) + t610 * pkin(5) + t622 * t685 - t621 * t693 - t716 * qJD(4) - t655 * pkin(2) - t656 * t729 + t772 * (t630 * qJD(5) + t622 * t686 - t656 * t761) + ((-t630 * t694 - t697 * t716) * r_i_i_C(1) + (-t630 * t697 + t694 * t716) * r_i_i_C(2)) * qJD(6), (t614 * t694 - t751 * t783) * r_i_i_C(1) + (t614 * t697 + t752 * t783) * r_i_i_C(2) - t614 * t693 - t783 * qJD(4) - t705 * t615 + t704 * t637, -t615, -t606 * t772 + t776 * t719 + t722 * t786, t790 * r_i_i_C(1) + t789 * r_i_i_C(2) + (-t787 * r_i_i_C(1) + t788 * r_i_i_C(2)) * qJD(6); 0, (t624 * t697 - t645 * t694) * r_i_i_C(1) + (-t624 * t694 - t645 * t697) * r_i_i_C(2) + t624 * pkin(5) + t646 * t685 + t645 * t693 + t664 * qJD(4) + t772 * (t649 * qJD(5) + t646 * t686 - t687 * t733) + ((-t649 * t694 + t664 * t697) * r_i_i_C(1) + (-t649 * t697 - t664 * t694) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t696 - t698 * t729) * t753, (t635 * t694 + t658 * t751) * r_i_i_C(1) + (t635 * t697 - t658 * t752) * r_i_i_C(2) - t635 * t693 + t658 * qJD(4) + t705 * t634 + t704 * t657, t634, t772 * t618 - t723 * t719 + t722 * (-t633 * qJD(5) - t635 * t686 + t687 * t734), (-t618 * t694 + t634 * t697) * r_i_i_C(1) + (-t618 * t697 - t634 * t694) * r_i_i_C(2) + ((-t633 * t697 - t657 * t694) * r_i_i_C(1) + (t633 * t694 - t657 * t697) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end