% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:56
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR15_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR15_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR15_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
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
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
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
	% StartTime: 2019-10-10 12:56:27
	% EndTime: 2019-10-10 12:56:28
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
	% StartTime: 2019-10-10 12:56:30
	% EndTime: 2019-10-10 12:56:31
	% DurationCPUTime: 1.21s
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
	% StartTime: 2019-10-10 12:56:34
	% EndTime: 2019-10-10 12:56:36
	% DurationCPUTime: 1.72s
	% Computational Cost: add. (1291->149), mult. (3990->253), div. (0->0), fcn. (4251->12), ass. (0->97)
	t563 = sin(qJ(2));
	t564 = sin(qJ(1));
	t626 = cos(pkin(6));
	t598 = t564 * t626;
	t593 = t563 * t598;
	t567 = cos(qJ(2));
	t568 = cos(qJ(1));
	t608 = t568 * t567;
	t610 = t564 * t563;
	t540 = -qJD(1) * t593 - qJD(2) * t610 + (qJD(2) * t626 + qJD(1)) * t608;
	t597 = t568 * t626;
	t549 = t563 * t597 + t564 * t567;
	t562 = sin(qJ(3));
	t566 = cos(qJ(3));
	t574 = t568 * t563 + t567 * t598;
	t539 = t574 * qJD(1) + t549 * qJD(2);
	t558 = sin(pkin(7));
	t560 = cos(pkin(7));
	t559 = sin(pkin(6));
	t607 = qJD(1) * t559;
	t602 = t564 * t607;
	t576 = -t539 * t560 + t558 * t602;
	t548 = -t567 * t597 + t610;
	t616 = t559 * t568;
	t584 = t548 * t560 + t558 * t616;
	t510 = (-qJD(3) * t549 + t576) * t562 + (-t584 * qJD(3) + t540) * t566;
	t623 = t539 * t558;
	t531 = t560 * t602 + t623;
	t561 = sin(qJ(4));
	t565 = cos(qJ(4));
	t527 = -t549 * t566 + t584 * t562;
	t543 = -t548 * t558 + t560 * t616;
	t635 = t527 * t565 + t543 * t561;
	t640 = t635 * qJD(4) - t510 * t561 + t531 * t565;
	t636 = t527 * t561 - t543 * t565;
	t639 = t636 * qJD(4) + t510 * t565 + t531 * t561;
	t627 = r_i_i_C(3) + qJ(5);
	t629 = r_i_i_C(2) - pkin(4);
	t573 = t627 * t561 - t629 * t565 + pkin(3);
	t630 = r_i_i_C(1) + pkin(11);
	t628 = pkin(10) * t558;
	t575 = t593 - t608;
	t617 = t559 * t564;
	t583 = t558 * t617 - t560 * t574;
	t529 = t583 * t562 - t566 * t575;
	t625 = t529 * t561;
	t537 = t548 * qJD(1) + t575 * qJD(2);
	t624 = t537 * t558;
	t620 = t558 * t561;
	t619 = t558 * t565;
	t618 = t558 * t567;
	t615 = t560 * t562;
	t614 = t560 * t566;
	t613 = t562 * t563;
	t612 = t562 * t567;
	t611 = t563 * t566;
	t609 = t566 * t567;
	t606 = qJD(2) * t559;
	t605 = t558 * t559 * t563;
	t603 = pkin(10) * t560 + pkin(9);
	t601 = t568 * t607;
	t600 = t558 * t606;
	t599 = t558 * t626;
	t596 = t558 * t601;
	t595 = t563 * t600;
	t594 = t567 * t600;
	t592 = qJD(3) * t599;
	t545 = t558 * t574 + t560 * t617;
	t588 = t529 * t565 + t545 * t561;
	t580 = t560 * t612 + t611;
	t542 = t580 * t559 + t562 * t599;
	t547 = -t559 * t618 + t626 * t560;
	t587 = t542 * t565 + t547 * t561;
	t535 = -t548 * t566 - t549 * t615;
	t586 = -t535 * t561 + t549 * t619;
	t536 = -t566 * t574 + t575 * t615;
	t585 = -t536 * t561 - t575 * t619;
	t582 = t560 * t609 - t613;
	t581 = -t560 * t611 - t612;
	t579 = t560 * t613 - t609;
	t546 = t579 * t559;
	t578 = t546 * t561 + t565 * t605;
	t577 = t560 * t601 - t624;
	t571 = t562 * t575 + t583 * t566;
	t570 = qJD(5) * t561 + (t629 * t561 + t627 * t565) * qJD(4);
	t569 = t527 * qJD(3) - t540 * t562 + t576 * t566;
	t538 = t549 * qJD(1) + t574 * qJD(2);
	t533 = (-t580 * qJD(2) + t581 * qJD(3)) * t559;
	t522 = t566 * t592 + (-t579 * qJD(2) + t582 * qJD(3)) * t559;
	t518 = -t540 * t615 - t539 * t566 + (t548 * t562 - t549 * t614) * qJD(3);
	t516 = t538 * t615 + t537 * t566 + (t562 * t574 + t575 * t614) * qJD(3);
	t513 = t587 * qJD(4) + t522 * t561 - t565 * t595;
	t508 = -t538 * t566 + (t537 * t560 + t596) * t562 + t571 * qJD(3);
	t507 = t529 * qJD(3) - t537 * t614 - t538 * t562 - t566 * t596;
	t498 = -qJD(4) * t625 + t577 * t561 + (t545 * qJD(4) + t508) * t565;
	t497 = t588 * qJD(4) + t508 * t561 - t577 * t565;
	t1 = [t636 * qJD(5) - t510 * pkin(3) - t540 * pkin(2) - pkin(10) * t623 + t630 * t569 + t629 * t639 + t627 * t640 + (-t568 * pkin(1) - t603 * t617) * qJD(1), -t585 * qJD(5) + t516 * pkin(3) + t537 * pkin(2) - t538 * t628 + t630 * (t536 * qJD(3) + t537 * t562 - t538 * t614) - t629 * (t585 * qJD(4) + t516 * t565 - t538 * t620) + t627 * (t538 * t619 + t516 * t561 + (t536 * t565 - t575 * t620) * qJD(4)), -t573 * t507 + t630 * t508 + t570 * t571, t588 * qJD(5) + t629 * t497 + t627 * t498, t497, 0; -(t545 * t565 - t625) * qJD(5) + t508 * pkin(3) - t538 * pkin(2) - pkin(10) * t624 + t630 * t507 - t629 * t498 + t627 * t497 + (-t564 * pkin(1) + t603 * t616) * qJD(1), -t586 * qJD(5) + t518 * pkin(3) - t539 * pkin(2) + t540 * t628 + t630 * (t535 * qJD(3) - t539 * t562 + t540 * t614) - t629 * (t586 * qJD(4) + t518 * t565 + t540 * t620) + t627 * (-t540 * t619 + t518 * t561 + (t535 * t565 + t549 * t620) * qJD(4)), t630 * t510 + t570 * (-t549 * t562 - t584 * t566) + t573 * t569, -t635 * qJD(5) + t627 * t639 - t629 * t640, -t640, 0; 0, -t578 * qJD(5) + t533 * pkin(3) - t630 * (-t582 * qJD(2) + t579 * qJD(3)) * t559 - t629 * (t578 * qJD(4) + t533 * t565 + t561 * t594) + t627 * (-t565 * t594 + t533 * t561 + (-t546 * t565 + t561 * t605) * qJD(4)) + (-pkin(2) * t563 + pkin(10) * t618) * t606, t630 * t522 + t570 * (t582 * t559 + t566 * t599) + t573 * (-t562 * t592 + (t581 * qJD(2) - t580 * qJD(3)) * t559), t587 * qJD(5) + t627 * (t561 * t595 + t522 * t565 + (-t542 * t561 + t547 * t565) * qJD(4)) + t629 * t513, t513, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:35
	% EndTime: 2019-10-10 12:56:38
	% DurationCPUTime: 3.06s
	% Computational Cost: add. (2360->202), mult. (7270->336), div. (0->0), fcn. (7905->14), ass. (0->123)
	t670 = cos(qJ(2));
	t671 = cos(qJ(1));
	t747 = cos(pkin(6));
	t719 = t671 * t747;
	t667 = sin(qJ(2));
	t749 = sin(qJ(1));
	t724 = t749 * t667;
	t645 = -t670 * t719 + t724;
	t663 = cos(pkin(7));
	t666 = sin(qJ(3));
	t646 = t667 * t719 + t749 * t670;
	t750 = cos(qJ(3));
	t729 = t646 * t750;
	t661 = sin(pkin(7));
	t662 = sin(pkin(6));
	t739 = t662 * t671;
	t730 = t661 * t739;
	t619 = (t645 * t663 + t730) * t666 - t729;
	t637 = -t645 * t661 + t663 * t739;
	t665 = sin(qJ(4));
	t669 = cos(qJ(4));
	t604 = t619 * t665 - t637 * t669;
	t706 = t750 * t730;
	t728 = t663 * t750;
	t755 = -t645 * t728 - t706;
	t616 = t646 * t666 - t755;
	t664 = sin(qJ(6));
	t668 = cos(qJ(6));
	t768 = t604 * t664 - t616 * t668;
	t767 = t604 * t668 + t616 * t664;
	t710 = t747 * t749;
	t684 = t671 * t667 + t670 * t710;
	t633 = t684 * qJD(1) + t646 * qJD(2);
	t701 = t667 * t710;
	t735 = t671 * t670;
	t634 = -qJD(1) * t701 - qJD(2) * t724 + (qJD(2) * t747 + qJD(1)) * t735;
	t721 = t749 * qJD(1);
	t712 = t662 * t721;
	t591 = (-qJD(3) * t646 - t633 * t663 + t661 * t712) * t666 + t634 * t750 + t755 * qJD(3);
	t745 = t633 * t661;
	t686 = t663 * t712 + t745;
	t766 = t604 * qJD(4) + t591 * t669 + t686 * t665;
	t762 = t619 * t669 + t637 * t665;
	t765 = t762 * qJD(4) - t591 * t665 + t686 * t669;
	t754 = pkin(5) + pkin(11);
	t685 = t701 - t735;
	t740 = t662 * t661;
	t715 = t749 * t740;
	t621 = -t685 * t750 + (-t684 * t663 + t715) * t666;
	t727 = t663 * t749;
	t639 = t684 * t661 + t662 * t727;
	t753 = -t621 * t665 + t639 * t669;
	t733 = r_i_i_C(3) + pkin(12) + pkin(4);
	t709 = t668 * r_i_i_C(1) - t664 * r_i_i_C(2);
	t699 = t709 + t754;
	t679 = t684 * t750;
	t694 = t750 * t715;
	t752 = -t663 * t679 + t666 * t685 + t694;
	t718 = t666 * t730;
	t738 = t663 * t666;
	t751 = qJD(1) * t694 - t633 * t728 - t634 * t666 + (t645 * t738 + t718 - t729) * qJD(3);
	t687 = t709 * qJD(6) + qJD(5);
	t748 = pkin(10) * t661;
	t742 = t661 * t665;
	t741 = t661 * t669;
	t737 = t666 * t667;
	t736 = t666 * t670;
	t734 = qJD(2) * t662;
	t732 = t667 * t740;
	t731 = t662 * t737;
	t726 = t750 * t667;
	t725 = t750 * t670;
	t722 = t661 * t734;
	t720 = t661 * t747;
	t716 = t663 * t725;
	t714 = t667 * t722;
	t713 = t670 * t722;
	t711 = t666 * t720;
	t708 = -t664 * r_i_i_C(1) - t668 * r_i_i_C(2);
	t707 = t662 * t716;
	t704 = t621 * t669 + t639 * t665;
	t689 = t663 * t736 + t726;
	t636 = t689 * t662 + t711;
	t644 = t747 * t663 - t670 * t740;
	t703 = t636 * t669 + t644 * t665;
	t612 = t636 * t665 - t644 * t669;
	t702 = t750 * t720;
	t700 = qJ(5) - t708;
	t628 = -t645 * t750 - t646 * t738;
	t698 = -t628 * t665 + t646 * t741;
	t630 = t685 * t738 - t679;
	t697 = -t630 * t665 - t685 * t741;
	t695 = qJD(6) * t708;
	t690 = -t663 * t737 + t725;
	t643 = t690 * t662;
	t693 = -t643 * t665 + t669 * t732;
	t691 = t645 * t666 - t646 * t728;
	t688 = t663 * t726 + t736;
	t682 = t685 * qJD(2);
	t677 = -t700 * t665 - t733 * t669 - pkin(3);
	t629 = -t684 * t666 - t685 * t728;
	t676 = t645 * qJD(1) + t682;
	t675 = t676 * t666;
	t674 = -t687 * t665 + (t733 * t665 - t700 * t669) * qJD(4);
	t673 = t676 * t750;
	t672 = t637 * qJD(1) - t661 * t682;
	t642 = t688 * t662;
	t635 = -t702 - t707 + t731;
	t632 = t646 * qJD(1) + t684 * qJD(2);
	t624 = (-t689 * qJD(2) - t688 * qJD(3)) * t662;
	t610 = qJD(3) * t702 + ((t716 - t737) * qJD(3) + t690 * qJD(2)) * t662;
	t609 = qJD(3) * t711 + (t688 * qJD(2) + t689 * qJD(3)) * t662;
	t599 = t691 * qJD(3) - t633 * t750 - t634 * t738;
	t597 = -t629 * qJD(3) + t632 * t738 + t673;
	t594 = t703 * qJD(4) + t610 * t665 - t669 * t714;
	t589 = qJD(1) * t718 + t752 * qJD(3) - t632 * t750 + t663 * t675;
	t588 = -qJD(1) * t706 + t621 * qJD(3) - t632 * t666 - t663 * t673;
	t584 = t632 * t741 + t597 * t665 + (t630 * t669 - t685 * t742) * qJD(4);
	t579 = t753 * qJD(4) + t589 * t669 + t672 * t665;
	t578 = t704 * qJD(4) + t589 * t665 - t672 * t669;
	t577 = t578 * t664 + t588 * t668 + (t664 * t752 - t668 * t753) * qJD(6);
	t576 = t578 * t668 - t588 * t664 + (t664 * t753 + t668 * t752) * qJD(6);
	t1 = [-pkin(10) * t745 - t634 * pkin(2) - t591 * pkin(3) + t604 * qJD(5) + t700 * t765 + t699 * t751 + (t767 * r_i_i_C(1) - t768 * r_i_i_C(2)) * qJD(6) - t733 * t766 + (-t671 * pkin(1) + (-t749 * pkin(9) - pkin(10) * t727) * t662) * qJD(1), (t584 * t664 + (-t629 * t664 - t668 * t697) * qJD(6)) * r_i_i_C(1) + (t584 * t668 + (-t629 * t668 + t664 * t697) * qJD(6)) * r_i_i_C(2) + t584 * qJ(5) - t697 * qJD(5) + t597 * pkin(3) + t676 * pkin(2) - t632 * t748 + t699 * (t630 * qJD(3) - t632 * t728 + t675) + t733 * (t697 * qJD(4) + t597 * t669 - t632 * t742), t677 * t588 + t699 * t589 + t621 * t695 - t674 * t752, -t733 * t578 + t700 * t579 + t687 * t704, t578, r_i_i_C(1) * t576 - r_i_i_C(2) * t577; -pkin(1) * t721 - t632 * pkin(2) + t589 * pkin(3) + t577 * r_i_i_C(1) + t576 * r_i_i_C(2) + t578 * qJ(5) - t753 * qJD(5) - t676 * t748 + (pkin(10) * t663 + pkin(9)) * qJD(1) * t739 + t754 * t588 + t733 * t579, t634 * t748 - t633 * pkin(2) + t599 * pkin(3) - t698 * qJD(5) + t700 * (-t634 * t741 + t599 * t665 + (t628 * t669 + t646 * t742) * qJD(4)) + t699 * (t628 * qJD(3) - t633 * t666 + t634 * t728) + ((t664 * t691 - t668 * t698) * r_i_i_C(1) + (t664 * t698 + t668 * t691) * r_i_i_C(2)) * qJD(6) + t733 * (t698 * qJD(4) + t599 * t669 + t634 * t742), t699 * t591 + t674 * t616 - t619 * t695 - t677 * t751, -t687 * t762 + t700 * t766 + t733 * t765, -t765, (t664 * t751 - t668 * t765) * r_i_i_C(1) + (t664 * t765 + t668 * t751) * r_i_i_C(2) + (t768 * r_i_i_C(1) + t767 * r_i_i_C(2)) * qJD(6); 0, t624 * pkin(3) - t693 * qJD(5) + t700 * (-t669 * t713 + t624 * t665 + (t643 * t669 + t665 * t732) * qJD(4)) - t699 * (-qJD(2) * t707 - t662 * qJD(3) * t725 + (qJD(3) * t663 + qJD(2)) * t731) + ((-t642 * t664 - t668 * t693) * r_i_i_C(1) + (-t642 * t668 + t664 * t693) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t667 + t670 * t748) * t734 + t733 * (t693 * qJD(4) + t624 * t669 + t665 * t713), t677 * t609 + t699 * t610 + t674 * t635 + t636 * t695, t687 * t703 + t700 * (-t612 * qJD(4) + t610 * t669 + t665 * t714) - t733 * t594, t594, (t594 * t668 - t609 * t664) * r_i_i_C(1) + (-t594 * t664 - t609 * t668) * r_i_i_C(2) + ((-t612 * t664 - t635 * t668) * r_i_i_C(1) + (-t612 * t668 + t635 * t664) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end