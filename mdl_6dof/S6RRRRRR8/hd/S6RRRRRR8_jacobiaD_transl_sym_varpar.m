% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
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
	% StartTime: 2019-10-10 13:29:31
	% EndTime: 2019-10-10 13:29:31
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
	% StartTime: 2019-10-10 13:29:32
	% EndTime: 2019-10-10 13:29:33
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
	% StartTime: 2019-10-10 13:29:34
	% EndTime: 2019-10-10 13:29:35
	% DurationCPUTime: 1.23s
	% Computational Cost: add. (730->138), mult. (2280->247), div. (0->0), fcn. (2372->12), ass. (0->84)
	t502 = sin(qJ(2));
	t503 = sin(qJ(1));
	t506 = cos(qJ(2));
	t507 = cos(qJ(1));
	t553 = cos(pkin(6));
	t528 = t507 * t553;
	t487 = t503 * t502 - t506 * t528;
	t497 = sin(pkin(7));
	t499 = cos(pkin(7));
	t498 = sin(pkin(6));
	t546 = t498 * t507;
	t520 = t487 * t499 + t497 * t546;
	t529 = t503 * t553;
	t524 = t502 * t529;
	t537 = qJD(2) * t502;
	t539 = t507 * t506;
	t479 = -qJD(1) * t524 - t503 * t537 + (qJD(2) * t553 + qJD(1)) * t539;
	t488 = t502 * t528 + t503 * t506;
	t501 = sin(qJ(3));
	t505 = cos(qJ(3));
	t511 = t507 * t502 + t506 * t529;
	t478 = t511 * qJD(1) + t488 * qJD(2);
	t538 = qJD(1) * t498;
	t531 = t503 * t538;
	t513 = -t478 * t499 + t497 * t531;
	t454 = (-qJD(3) * t488 + t513) * t501 + (-t520 * qJD(3) + t479) * t505;
	t551 = t478 * t497;
	t470 = t499 * t531 + t551;
	t500 = sin(qJ(4));
	t504 = cos(qJ(4));
	t562 = t454 * t500 - t470 * t504;
	t561 = -t454 * t504 - t470 * t500;
	t465 = -t488 * t505 + t520 * t501;
	t482 = -t487 * t497 + t499 * t546;
	t560 = -t465 * t500 + t482 * t504;
	t559 = t465 * t504 + t482 * t500;
	t554 = r_i_i_C(3) + pkin(11);
	t512 = t524 - t539;
	t476 = t487 * qJD(1) + t512 * qJD(2);
	t552 = t476 * t497;
	t548 = t497 * t498;
	t547 = t498 * t503;
	t545 = t499 * t501;
	t544 = t499 * t505;
	t543 = t501 * t502;
	t542 = t501 * t506;
	t541 = t502 * t505;
	t540 = t505 * t506;
	t536 = qJD(2) * t506;
	t535 = qJD(4) * t500;
	t534 = qJD(4) * t504;
	t532 = pkin(10) * t499 + pkin(9);
	t530 = t507 * t538;
	t527 = t553 * t497;
	t526 = t497 * t530;
	t525 = t537 * t548;
	t523 = qJD(3) * t527;
	t521 = t504 * r_i_i_C(1) - t500 * r_i_i_C(2) + pkin(3);
	t474 = -t487 * t505 - t488 * t545;
	t475 = -t505 * t511 + t512 * t545;
	t519 = t497 * t547 - t499 * t511;
	t518 = t499 * t540 - t543;
	t517 = -t499 * t541 - t542;
	t516 = t499 * t542 + t541;
	t515 = t499 * t543 - t540;
	t514 = qJD(4) * (-t500 * r_i_i_C(1) - t504 * r_i_i_C(2));
	t509 = t501 * t512 + t519 * t505;
	t467 = t519 * t501 - t505 * t512;
	t508 = t465 * qJD(3) - t479 * t501 + t513 * t505;
	t486 = t553 * t499 - t506 * t548;
	t485 = t515 * t498;
	t484 = t497 * t511 + t499 * t547;
	t481 = t516 * t498 + t501 * t527;
	t477 = t488 * qJD(1) + t511 * qJD(2);
	t472 = (-t516 * qJD(2) + t517 * qJD(3)) * t498;
	t468 = t499 * t530 - t552;
	t462 = t505 * t523 + (-t515 * qJD(2) + t518 * qJD(3)) * t498;
	t460 = -t479 * t545 - t478 * t505 + (t487 * t501 - t488 * t544) * qJD(3);
	t458 = t477 * t545 + t476 * t505 + (t501 * t511 + t512 * t544) * qJD(3);
	t452 = -t477 * t505 + (t476 * t499 + t526) * t501 + t509 * qJD(3);
	t451 = t467 * qJD(3) - t476 * t544 - t477 * t501 - t505 * t526;
	t450 = t452 * t504 + t468 * t500 + (-t467 * t500 + t484 * t504) * qJD(4);
	t449 = -t452 * t500 + t468 * t504 + (-t467 * t504 - t484 * t500) * qJD(4);
	t1 = [t561 * r_i_i_C(1) + t562 * r_i_i_C(2) - t454 * pkin(3) - t479 * pkin(2) - pkin(10) * t551 + t554 * t508 + (t560 * r_i_i_C(1) - t559 * r_i_i_C(2)) * qJD(4) + (-t507 * pkin(1) - t532 * t547) * qJD(1), (t458 * t504 - t475 * t535) * r_i_i_C(1) + (-t458 * t500 - t475 * t534) * r_i_i_C(2) + t458 * pkin(3) + t476 * pkin(2) + t554 * (t475 * qJD(3) + t476 * t501 - t477 * t544) + ((-t477 * t500 - t512 * t534) * r_i_i_C(1) + (-t477 * t504 + t512 * t535) * r_i_i_C(2) - t477 * pkin(10)) * t497, -t521 * t451 + t554 * t452 + t509 * t514, t449 * r_i_i_C(1) - t450 * r_i_i_C(2), 0, 0; -pkin(10) * t552 - t477 * pkin(2) + t452 * pkin(3) + t450 * r_i_i_C(1) + t449 * r_i_i_C(2) + t554 * t451 + (-pkin(1) * t503 + t532 * t546) * qJD(1), (t460 * t504 - t474 * t535) * r_i_i_C(1) + (-t460 * t500 - t474 * t534) * r_i_i_C(2) + t460 * pkin(3) - t478 * pkin(2) + t554 * (t474 * qJD(3) - t478 * t501 + t479 * t544) + ((t479 * t500 + t488 * t534) * r_i_i_C(1) + (t479 * t504 - t488 * t535) * r_i_i_C(2) + t479 * pkin(10)) * t497, t554 * t454 + (-t488 * t501 - t520 * t505) * t514 + t521 * t508, -t562 * r_i_i_C(1) + t561 * r_i_i_C(2) + (t559 * r_i_i_C(1) + t560 * r_i_i_C(2)) * qJD(4), 0, 0; 0, (t472 * t504 + t485 * t535) * r_i_i_C(1) + (-t472 * t500 + t485 * t534) * r_i_i_C(2) + t472 * pkin(3) + (-t554 * (-t518 * qJD(2) + t515 * qJD(3)) - pkin(2) * t537 + ((t500 * t536 + t502 * t534) * r_i_i_C(1) + (-t502 * t535 + t504 * t536) * r_i_i_C(2) + pkin(10) * t536) * t497) * t498, t554 * t462 + (t518 * t498 + t505 * t527) * t514 + t521 * (-t501 * t523 + (t517 * qJD(2) - t516 * qJD(3)) * t498), (-t462 * t500 + t504 * t525) * r_i_i_C(1) + (-t462 * t504 - t500 * t525) * r_i_i_C(2) + ((-t481 * t504 - t486 * t500) * r_i_i_C(1) + (t481 * t500 - t486 * t504) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:34
	% EndTime: 2019-10-10 13:29:36
	% DurationCPUTime: 1.55s
	% Computational Cost: add. (1151->163), mult. (3096->278), div. (0->0), fcn. (3240->14), ass. (0->97)
	t544 = sin(qJ(1));
	t543 = sin(qJ(2));
	t609 = cos(pkin(6));
	t578 = t544 * t609;
	t567 = t543 * t578;
	t587 = qJD(2) * t543;
	t547 = cos(qJ(2));
	t548 = cos(qJ(1));
	t592 = t548 * t547;
	t515 = -qJD(1) * t567 - t544 * t587 + (qJD(2) * t609 + qJD(1)) * t592;
	t577 = t548 * t609;
	t524 = t543 * t577 + t544 * t547;
	t542 = sin(qJ(3));
	t546 = cos(qJ(3));
	t554 = t548 * t543 + t547 * t578;
	t514 = t554 * qJD(1) + t524 * qJD(2);
	t538 = sin(pkin(7));
	t540 = cos(pkin(7));
	t539 = sin(pkin(6));
	t588 = qJD(1) * t539;
	t581 = t544 * t588;
	t557 = -t514 * t540 + t538 * t581;
	t523 = t544 * t543 - t547 * t577;
	t599 = t539 * t548;
	t563 = t523 * t540 + t538 * t599;
	t488 = (-qJD(3) * t524 + t557) * t542 + (-t563 * qJD(3) + t515) * t546;
	t518 = -t523 * t538 + t540 * t599;
	t536 = qJD(4) + qJD(5);
	t613 = t518 * t536 - t488;
	t501 = -t524 * t546 + t563 * t542;
	t607 = t514 * t538;
	t506 = t540 * t581 + t607;
	t612 = -t501 * t536 - t506;
	t610 = r_i_i_C(3) + pkin(12) + pkin(11);
	t555 = t567 - t592;
	t512 = t523 * qJD(1) + t555 * qJD(2);
	t608 = t512 * t538;
	t537 = qJ(4) + qJ(5);
	t534 = sin(t537);
	t604 = t534 * t536;
	t535 = cos(t537);
	t603 = t535 * t536;
	t602 = t536 * t543;
	t601 = t538 * t539;
	t600 = t539 * t544;
	t598 = t540 * t542;
	t597 = t540 * t546;
	t596 = t542 * t543;
	t595 = t542 * t547;
	t594 = t543 * t546;
	t593 = t546 * t547;
	t562 = t538 * t600 - t540 * t554;
	t503 = t562 * t542 - t546 * t555;
	t580 = t548 * t588;
	t504 = t540 * t580 - t608;
	t570 = -t503 * t536 + t504;
	t513 = t524 * qJD(1) + t554 * qJD(2);
	t551 = t542 * t555 + t562 * t546;
	t569 = t538 * t580;
	t486 = -t513 * t546 + (t512 * t540 + t569) * t542 + t551 * qJD(3);
	t520 = t538 * t554 + t540 * t600;
	t576 = t520 * t536 + t486;
	t483 = -t576 * t534 + t570 * t535;
	t484 = t570 * t534 + t576 * t535;
	t591 = t483 * r_i_i_C(1) - t484 * r_i_i_C(2);
	t590 = (t534 * t613 - t535 * t612) * r_i_i_C(1) + (t534 * t612 + t535 * t613) * r_i_i_C(2);
	t559 = t540 * t595 + t594;
	t579 = t538 * t609;
	t517 = t559 * t539 + t542 * t579;
	t568 = t587 * t601;
	t556 = -t517 * t536 + t568;
	t558 = t540 * t596 - t593;
	t561 = t540 * t593 - t596;
	t566 = qJD(3) * t579;
	t498 = t546 * t566 + (-t558 * qJD(2) + t561 * qJD(3)) * t539;
	t522 = t609 * t540 - t547 * t601;
	t573 = -t522 * t536 - t498;
	t589 = (t573 * t534 + t556 * t535) * r_i_i_C(1) + (-t556 * t534 + t573 * t535) * r_i_i_C(2);
	t586 = qJD(2) * t547;
	t545 = cos(qJ(4));
	t585 = qJD(4) * t545;
	t541 = sin(qJ(4));
	t584 = qJD(4) * t541 * pkin(4);
	t582 = pkin(10) * t540 + pkin(9);
	t533 = pkin(4) * t545 + pkin(3);
	t564 = t535 * r_i_i_C(1) - t534 * r_i_i_C(2) + t533;
	t510 = -t523 * t546 - t524 * t598;
	t511 = -t546 * t554 + t555 * t598;
	t560 = -t540 * t594 - t595;
	t552 = -t584 + (-t534 * r_i_i_C(1) - t535 * r_i_i_C(2)) * t536;
	t550 = t501 * qJD(3) - t515 * t542 + t557 * t546;
	t521 = t558 * t539;
	t508 = (-t559 * qJD(2) + t560 * qJD(3)) * t539;
	t496 = -t515 * t598 - t514 * t546 + (t523 * t542 - t524 * t597) * qJD(3);
	t494 = t513 * t598 + t512 * t546 + (t542 * t554 + t555 * t597) * qJD(3);
	t485 = t503 * qJD(3) - t512 * t597 - t513 * t542 - t546 * t569;
	t1 = [-pkin(10) * t607 - t515 * pkin(2) - t488 * t533 + (r_i_i_C(1) * t613 + r_i_i_C(2) * t612) * t535 + (r_i_i_C(1) * t612 - r_i_i_C(2) * t613) * t534 + t610 * t550 + (-t548 * pkin(1) - t582 * t600) * qJD(1) + (-t506 * t541 + (-t501 * t541 + t518 * t545) * qJD(4)) * pkin(4), (t494 * t535 - t511 * t604) * r_i_i_C(1) + (-t494 * t534 - t511 * t603) * r_i_i_C(2) + t494 * t533 - t511 * t584 + t512 * pkin(2) + t610 * (t511 * qJD(3) + t512 * t542 - t513 * t597) + ((-t513 * t534 - t555 * t603) * r_i_i_C(1) + (-t513 * t535 + t555 * t604) * r_i_i_C(2) - t513 * pkin(10) + (-t513 * t541 - t555 * t585) * pkin(4)) * t538, -t564 * t485 + t610 * t486 + t552 * t551, (-t486 * t541 + t504 * t545 + (-t503 * t545 - t520 * t541) * qJD(4)) * pkin(4) + t591, t591, 0; -pkin(10) * t608 - t513 * pkin(2) + t484 * r_i_i_C(1) + t483 * r_i_i_C(2) + t486 * t533 + t610 * t485 + (-pkin(1) * t544 + t582 * t599) * qJD(1) + (t504 * t541 + (-t503 * t541 + t520 * t545) * qJD(4)) * pkin(4), (t496 * t535 - t510 * t604) * r_i_i_C(1) + (-t496 * t534 - t510 * t603) * r_i_i_C(2) + t496 * t533 - t510 * t584 - t514 * pkin(2) + t610 * (t510 * qJD(3) - t514 * t542 + t515 * t597) + ((t515 * t534 + t524 * t603) * r_i_i_C(1) + (t515 * t535 - t524 * t604) * r_i_i_C(2) + t515 * pkin(10) + (t515 * t541 + t524 * t585) * pkin(4)) * t538, t610 * t488 + t552 * (-t524 * t542 - t563 * t546) + t564 * t550, (-t488 * t541 + t506 * t545 + (t501 * t545 + t518 * t541) * qJD(4)) * pkin(4) + t590, t590, 0; 0, (t508 * t535 + t521 * t604) * r_i_i_C(1) + (-t508 * t534 + t521 * t603) * r_i_i_C(2) + t508 * t533 + t521 * t584 + (-t610 * (-t561 * qJD(2) + t558 * qJD(3)) - pkin(2) * t587 + ((t534 * t586 + t535 * t602) * r_i_i_C(1) + (-t534 * t602 + t535 * t586) * r_i_i_C(2) + pkin(10) * t586 + (t541 * t586 + t543 * t585) * pkin(4)) * t538) * t539, t610 * t498 + t552 * (t561 * t539 + t546 * t579) + t564 * (-t542 * t566 + (t560 * qJD(2) - t559 * qJD(3)) * t539), (t545 * t568 - t498 * t541 + (-t517 * t545 - t522 * t541) * qJD(4)) * pkin(4) + t589, t589, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:42
	% EndTime: 2019-10-10 13:29:46
	% DurationCPUTime: 3.91s
	% Computational Cost: add. (2908->254), mult. (7495->421), div. (0->0), fcn. (8135->16), ass. (0->159)
	t741 = qJ(4) + qJ(5);
	t739 = cos(t741);
	t748 = sin(qJ(2));
	t751 = cos(qJ(2));
	t752 = cos(qJ(1));
	t842 = cos(pkin(6));
	t803 = t752 * t842;
	t846 = sin(qJ(1));
	t721 = t748 * t803 + t846 * t751;
	t789 = t842 * t846;
	t768 = t752 * t748 + t751 * t789;
	t708 = t768 * qJD(1) + t721 * qJD(2);
	t780 = t748 * t789;
	t809 = t846 * t748;
	t824 = t752 * t751;
	t709 = -qJD(1) * t780 - qJD(2) * t809 + (qJD(2) * t842 + qJD(1)) * t824;
	t720 = -t751 * t803 + t809;
	t742 = sin(pkin(7));
	t744 = cos(pkin(7));
	t747 = sin(qJ(3));
	t743 = sin(pkin(6));
	t828 = t743 * t752;
	t816 = t742 * t828;
	t847 = cos(qJ(3));
	t785 = t847 * t816;
	t805 = t847 * qJD(3);
	t791 = t744 * t805;
	t804 = t846 * qJD(1);
	t792 = t743 * t804;
	t661 = (-qJD(3) * t721 - t708 * t744 + t742 * t792) * t747 - qJD(3) * t785 + t709 * t847 - t720 * t791;
	t712 = -t720 * t742 + t744 * t828;
	t740 = qJD(4) + qJD(5);
	t801 = -t712 * t740 + t661;
	t864 = t801 * t739;
	t814 = t721 * t847;
	t692 = (t720 * t744 + t816) * t747 - t814;
	t738 = sin(t741);
	t677 = t692 * t739 + t712 * t738;
	t813 = t744 * t847;
	t689 = t720 * t813 + t721 * t747 + t785;
	t745 = sin(qJ(6));
	t749 = cos(qJ(6));
	t863 = -t677 * t745 - t689 * t749;
	t862 = t677 * t749 - t689 * t745;
	t839 = t708 * t742;
	t696 = t744 * t792 + t839;
	t798 = t692 * t740 + t696;
	t861 = -t801 * t738 + t798 * t739;
	t832 = t742 * t743;
	t794 = t846 * t832;
	t777 = t847 * t794;
	t796 = t747 * t816;
	t827 = t744 * t747;
	t662 = qJD(1) * t777 - t708 * t813 - t709 * t747 + (t720 * t827 + t796 - t814) * qJD(3);
	t860 = t662 * t745;
	t859 = t662 * t749;
	t848 = r_i_i_C(3) + pkin(13);
	t857 = -t749 * r_i_i_C(1) - pkin(5);
	t823 = qJD(6) * t745;
	t820 = r_i_i_C(1) * t823;
	t822 = qJD(6) * t749;
	t852 = -t822 * r_i_i_C(2) - t820;
	t782 = -t745 * r_i_i_C(2) - t857;
	t769 = t780 - t824;
	t694 = -t769 * t847 + (-t744 * t768 + t794) * t747;
	t765 = t769 * qJD(2);
	t755 = t712 * qJD(1) - t742 * t765;
	t843 = pkin(4) * qJD(4);
	t851 = pkin(4) * t755 - t694 * t843;
	t850 = t747 * t769 + t777;
	t707 = t721 * qJD(1) + t768 * qJD(2);
	t759 = t720 * qJD(1) + t765;
	t758 = t759 * t747;
	t659 = qJD(1) * t796 + t850 * qJD(3) - t707 * t847 + t744 * t758 - t768 * t791;
	t812 = t744 * t846;
	t714 = t742 * t768 + t743 * t812;
	t834 = t739 * t740;
	t648 = t694 * t834 - t755 * t739 + (t714 * t740 + t659) * t738;
	t836 = t738 * t740;
	t649 = t659 * t739 - t694 * t836 + t714 * t834 + t755 * t738;
	t678 = -t694 * t738 + t714 * t739;
	t849 = (t648 * t745 - t678 * t822) * r_i_i_C(2) - t678 * t820 + t848 * t649 + t857 * t648;
	t746 = sin(qJ(4));
	t845 = pkin(4) * t746;
	t840 = t707 * t742;
	t838 = t709 * t742;
	t835 = t738 * t742;
	t833 = t740 * t742;
	t831 = t742 * t746;
	t750 = cos(qJ(4));
	t830 = t742 * t750;
	t829 = t743 * t751;
	t826 = t747 * t748;
	t825 = t747 * t751;
	t821 = t746 * t843;
	t818 = t714 * t843;
	t817 = t748 * t832;
	t815 = t743 * t826;
	t811 = t847 * t748;
	t810 = t847 * t751;
	t808 = qJD(1) * t828;
	t807 = qJD(2) * t832;
	t802 = t842 * t742;
	t773 = -t744 * t826 + t810;
	t781 = t847 * t802;
	t795 = t744 * t810;
	t685 = qJD(3) * t781 + ((t795 - t826) * qJD(3) + t773 * qJD(2)) * t743;
	t719 = -t742 * t829 + t842 * t744;
	t799 = t719 * t740 + t685;
	t797 = t830 * t843;
	t793 = t748 * t807;
	t790 = t747 * t802;
	t756 = t759 * t847;
	t774 = t747 * t768 + t769 * t813;
	t670 = t774 * qJD(3) + t707 * t827 + t756;
	t788 = -t769 * t833 + t670;
	t775 = t720 * t747 - t721 * t813;
	t672 = t775 * qJD(3) - t708 * t847 - t709 * t827;
	t787 = t721 * t833 + t672;
	t786 = t743 * t795;
	t702 = -t720 * t847 - t721 * t827;
	t784 = t702 * t740 - t838;
	t704 = -t768 * t847 + t769 * t827;
	t783 = t704 * t740 + t840;
	t771 = t744 * t811 + t825;
	t772 = t744 * t825 + t811;
	t698 = (-t772 * qJD(2) - t771 * qJD(3)) * t743;
	t779 = t740 * t817 + t698;
	t718 = t773 * t743;
	t770 = -t718 * t740 + t751 * t807;
	t737 = t750 * pkin(4) + pkin(3);
	t762 = -t848 * t738 - t782 * t739 - t737;
	t651 = t692 * t836 + t696 * t738 + t864;
	t761 = t852 * (t692 * t738 - t712 * t739) + t848 * t651 + t782 * t861;
	t711 = t772 * t743 + t790;
	t668 = -t711 * t836 + t738 * t793 + t799 * t739;
	t760 = t852 * (-t711 * t738 + t719 * t739) + t848 * t668 + t782 * ((-t711 * t740 + t793) * t739 - t799 * t738);
	t757 = t821 + (t745 * r_i_i_C(1) + t749 * r_i_i_C(2)) * t739 * qJD(6) + (t782 * t738 - t848 * t739) * t740;
	t753 = -pkin(12) - pkin(11);
	t717 = t771 * t743;
	t710 = -t781 - t786 + t815;
	t705 = t718 * t739 + t738 * t817;
	t697 = -qJD(2) * t786 - t805 * t829 + (qJD(3) * t744 + qJD(2)) * t815;
	t693 = t768 * t813 - t850;
	t687 = t711 * t739 + t719 * t738;
	t684 = qJD(3) * t790 + (t771 * qJD(2) + t772 * qJD(3)) * t743;
	t681 = t704 * t739 - t769 * t835;
	t680 = t702 * t739 + t721 * t835;
	t679 = t694 * t739 + t714 * t738;
	t674 = t770 * t738 + t779 * t739;
	t671 = t702 * qJD(3) - t708 * t747 + t709 * t813;
	t669 = t704 * qJD(3) - t707 * t813 + t758;
	t658 = -qJD(1) * t785 + t694 * qJD(3) - t707 * t747 - t744 * t756;
	t657 = -t784 * t738 + t787 * t739;
	t655 = -t783 * t738 + t788 * t739;
	t653 = -t798 * t738 - t864;
	t641 = t649 * t749 + t658 * t745 + (-t679 * t745 + t693 * t749) * qJD(6);
	t640 = -t649 * t745 + t658 * t749 + (-t679 * t749 - t693 * t745) * qJD(6);
	t1 = [(t653 * t749 + t860) * r_i_i_C(1) + (-t653 * t745 + t859) * r_i_i_C(2) + t653 * pkin(5) - t661 * t737 - t662 * t753 - t709 * pkin(2) - pkin(10) * t839 + t848 * t861 + (t863 * r_i_i_C(1) - t862 * r_i_i_C(2)) * qJD(6) + (-t752 * pkin(1) + (-t846 * pkin(9) - pkin(10) * t812) * t743) * qJD(1) + (-t696 * t746 + (-t692 * t746 + t712 * t750) * qJD(4)) * pkin(4), (t655 * t749 + t669 * t745 + (-t681 * t745 - t749 * t774) * qJD(6)) * r_i_i_C(1) + (-t655 * t745 + t669 * t749 + (-t681 * t749 + t745 * t774) * qJD(6)) * r_i_i_C(2) + t655 * pkin(5) + t670 * t737 - t704 * t821 - t669 * t753 - t707 * pkin(4) * t831 - t769 * t797 + t759 * pkin(2) - pkin(10) * t840 + t848 * (t788 * t738 + t783 * t739), (t659 * t745 + t694 * t822) * r_i_i_C(1) + (t659 * t749 - t694 * t823) * r_i_i_C(2) - t659 * t753 + t762 * t658 + t757 * t693, -t659 * t845 - t746 * t818 + t851 * t750 + t849, t849, t640 * r_i_i_C(1) - t641 * r_i_i_C(2); -pkin(1) * t804 - t707 * pkin(2) + t649 * pkin(5) + pkin(9) * t808 + t641 * r_i_i_C(1) + t640 * r_i_i_C(2) - t658 * t753 + t659 * t737 + t750 * t818 + t851 * t746 + t848 * t648 + (-t742 * t759 + t744 * t808) * pkin(10), (t657 * t749 + t671 * t745) * r_i_i_C(1) + (-t657 * t745 + t671 * t749) * r_i_i_C(2) + t657 * pkin(5) + t672 * t737 - t671 * t753 - t708 * pkin(2) + pkin(10) * t838 + t848 * (t787 * t738 + t784 * t739) + ((-t680 * t745 - t749 * t775) * r_i_i_C(1) + (-t680 * t749 + t745 * t775) * r_i_i_C(2)) * qJD(6) + (t709 * t831 + (-t702 * t746 + t721 * t830) * qJD(4)) * pkin(4), (t661 * t745 - t692 * t822) * r_i_i_C(1) + (t661 * t749 + t692 * t823) * r_i_i_C(2) - t661 * t753 - t762 * t662 + t757 * t689, (-t661 * t746 + t696 * t750 + (t692 * t750 + t712 * t746) * qJD(4)) * pkin(4) + t761, t761, (-t651 * t745 - t859) * r_i_i_C(1) + (-t651 * t749 + t860) * r_i_i_C(2) + (t862 * r_i_i_C(1) + t863 * r_i_i_C(2)) * qJD(6); 0, (t674 * t749 - t697 * t745) * r_i_i_C(1) + (-t674 * t745 - t697 * t749) * r_i_i_C(2) + t674 * pkin(5) + t698 * t737 - t718 * t821 + t697 * t753 + t848 * (t779 * t738 - t770 * t739) + ((-t705 * t745 + t717 * t749) * r_i_i_C(1) + (-t705 * t749 - t717 * t745) * r_i_i_C(2)) * qJD(6) + (t748 * t797 + (-pkin(2) * t748 + (pkin(10) + t845) * t751 * t742) * qJD(2)) * t743, (t685 * t745 + t711 * t822) * r_i_i_C(1) + (t685 * t749 - t711 * t823) * r_i_i_C(2) - t685 * t753 + t762 * t684 + t757 * t710, (t750 * t793 - t685 * t746 + (-t711 * t750 - t719 * t746) * qJD(4)) * pkin(4) + t760, t760, (-t668 * t745 + t684 * t749) * r_i_i_C(1) + (-t668 * t749 - t684 * t745) * r_i_i_C(2) + ((-t687 * t749 - t710 * t745) * r_i_i_C(1) + (t687 * t745 - t710 * t749) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end