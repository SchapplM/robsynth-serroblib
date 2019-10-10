% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR12_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR12_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR12_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
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
	% StartTime: 2019-10-10 12:49:59
	% EndTime: 2019-10-10 12:49:59
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
	% StartTime: 2019-10-10 12:50:00
	% EndTime: 2019-10-10 12:50:00
	% DurationCPUTime: 0.48s
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
	t1 = [-t290 * pkin(2) + t295 * r_i_i_C(1) - t289 * t342 + (-t309 * pkin(1) - t341 * t339) * qJD(1) + (t320 * r_i_i_C(1) - t312) * t307 + (-t314 * r_i_i_C(2) - t344) * t304, t319 * t287 - t345 * t288 + ((t304 * t317 + t316 * t336) * r_i_i_C(1) + (t307 * t317 - t316 * t337) * r_i_i_C(2)) * qJD(3), t285 * r_i_i_C(1) - t286 * r_i_i_C(2), 0, 0, 0; -t288 * pkin(2) + t286 * r_i_i_C(1) + t285 * r_i_i_C(2) - t287 * t342 + (-t306 * pkin(1) + t338 * t341) * qJD(1), -t319 * t289 + t345 * t290 + ((t291 * t304 - t292 * t336) * r_i_i_C(1) + (t291 * t307 + t292 * t337) * r_i_i_C(2)) * qJD(3), t295 * r_i_i_C(2) + (t320 * r_i_i_C(2) + t344) * t307 + (t314 * r_i_i_C(1) - t312) * t304, 0, 0, 0; 0, (t311 * qJD(3) + (-t305 * pkin(2) + t308 * t342 + t310) * qJD(2)) * t301, -t321 * t303 * t327 + (t311 * qJD(2) + t310 * qJD(3)) * t301, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:50:02
	% EndTime: 2019-10-10 12:50:03
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
	% StartTime: 2019-10-10 12:50:02
	% EndTime: 2019-10-10 12:50:04
	% DurationCPUTime: 1.63s
	% Computational Cost: add. (1013->180), mult. (2885->302), div. (0->0), fcn. (3026->14), ass. (0->99)
	t522 = sin(qJ(1));
	t521 = sin(qJ(2));
	t578 = cos(pkin(6));
	t551 = t522 * t578;
	t544 = t521 * t551;
	t562 = qJD(2) * t521;
	t525 = cos(qJ(2));
	t526 = cos(qJ(1));
	t564 = t526 * t525;
	t490 = -qJD(1) * t544 - t522 * t562 + (qJD(2) * t578 + qJD(1)) * t564;
	t515 = sin(pkin(7));
	t516 = sin(pkin(6));
	t571 = t516 * t526;
	t555 = t515 * t571;
	t589 = -qJD(3) * t555 + t490;
	t550 = t526 * t578;
	t499 = t521 * t550 + t522 * t525;
	t530 = t526 * t521 + t525 * t551;
	t489 = t530 * qJD(1) + t499 * qJD(2);
	t517 = cos(pkin(7));
	t520 = sin(qJ(3));
	t524 = cos(qJ(3));
	t563 = qJD(1) * t516;
	t553 = t522 * t563;
	t548 = t515 * t553;
	t498 = t522 * t521 - t525 * t550;
	t575 = t498 * t517;
	t464 = (-qJD(3) * t499 - t489 * t517 + t548) * t520 + (-qJD(3) * t575 + t589) * t524;
	t576 = t489 * t515;
	t480 = t517 * t553 + t576;
	t514 = qJ(4) + pkin(13);
	t512 = sin(t514);
	t513 = cos(t514);
	t588 = t464 * t512 - t480 * t513;
	t587 = -t464 * t513 - t480 * t512;
	t570 = t517 * t520;
	t539 = t498 * t570 - t499 * t524;
	t475 = t520 * t555 + t539;
	t493 = -t498 * t515 + t517 * t571;
	t584 = t475 * t513 + t493 * t512;
	t583 = -t475 * t512 + t493 * t513;
	t582 = (t555 + t575) * t524 + t499 * t520;
	t569 = t517 * t524;
	t581 = t539 * qJD(3) - t489 * t569 - t589 * t520 + t524 * t548;
	t519 = sin(qJ(4));
	t580 = t519 * pkin(4);
	t579 = r_i_i_C(3) + qJ(5) + pkin(11);
	t531 = t544 - t564;
	t487 = t498 * qJD(1) + t531 * qJD(2);
	t577 = t487 * t515;
	t573 = t515 * t516;
	t572 = t516 * t522;
	t568 = t520 * t521;
	t567 = t520 * t525;
	t566 = t521 * t524;
	t565 = t524 * t525;
	t561 = qJD(2) * t525;
	t560 = qJD(4) * t512;
	t559 = qJD(4) * t513;
	t558 = qJD(4) * t521;
	t523 = cos(qJ(4));
	t557 = qJD(4) * t523;
	t556 = qJD(4) * t580;
	t554 = pkin(10) * t517 + pkin(9);
	t552 = t526 * t563;
	t549 = t578 * t515;
	t547 = t515 * t552;
	t546 = t562 * t573;
	t543 = qJD(3) * t549;
	t511 = t523 * pkin(4) + pkin(3);
	t541 = -t513 * r_i_i_C(1) + t512 * r_i_i_C(2) - t511;
	t540 = t498 * t520 - t499 * t569;
	t485 = -t498 * t524 - t499 * t570;
	t537 = t520 * t530 + t531 * t569;
	t486 = -t524 * t530 + t531 * t570;
	t536 = t515 * t572 - t517 * t530;
	t535 = t517 * t565 - t568;
	t534 = t517 * t566 + t567;
	t533 = t517 * t567 + t566;
	t532 = t517 * t568 - t565;
	t528 = qJD(4) * (-t512 * r_i_i_C(1) - t513 * r_i_i_C(2) - t580);
	t476 = t520 * t531 + t536 * t524;
	t477 = t536 * t520 - t524 * t531;
	t497 = t578 * t517 - t525 * t573;
	t496 = t532 * t516;
	t495 = t515 * t530 + t517 * t572;
	t492 = t533 * t516 + t520 * t549;
	t488 = t499 * qJD(1) + t530 * qJD(2);
	t482 = (-t533 * qJD(2) - t534 * qJD(3)) * t516;
	t478 = t517 * t552 - t577;
	t472 = t524 * t543 + (-t532 * qJD(2) + t535 * qJD(3)) * t516;
	t471 = t520 * t543 + (t534 * qJD(2) + t533 * qJD(3)) * t516;
	t470 = t540 * qJD(3) - t489 * t524 - t490 * t570;
	t468 = t537 * qJD(3) + t487 * t524 + t488 * t570;
	t462 = -t488 * t524 + (t487 * t517 + t547) * t520 + t476 * qJD(3);
	t461 = t477 * qJD(3) - t487 * t569 - t488 * t520 - t524 * t547;
	t460 = t462 * t513 + t478 * t512 + (-t477 * t512 + t495 * t513) * qJD(4);
	t459 = -t462 * t512 + t478 * t513 + (-t477 * t513 - t495 * t512) * qJD(4);
	t1 = [t587 * r_i_i_C(1) + t588 * r_i_i_C(2) - t464 * t511 - t582 * qJD(5) - t480 * t580 - t490 * pkin(2) - pkin(10) * t576 + t579 * t581 + (-t526 * pkin(1) - t554 * t572) * qJD(1) + (t583 * r_i_i_C(1) - t584 * r_i_i_C(2) + (-t475 * t519 + t493 * t523) * pkin(4)) * qJD(4), (t468 * t513 - t486 * t560) * r_i_i_C(1) + (-t468 * t512 - t486 * t559) * r_i_i_C(2) + t468 * t511 - t486 * t556 - t537 * qJD(5) + t487 * pkin(2) + t579 * (t486 * qJD(3) + t487 * t520 - t488 * t569) + ((-t488 * t512 - t531 * t559) * r_i_i_C(1) + (-t488 * t513 + t531 * t560) * r_i_i_C(2) - t488 * pkin(10) + (-t488 * t519 - t531 * t557) * pkin(4)) * t515, t477 * qJD(5) + t541 * t461 + t579 * t462 + t476 * t528, t459 * r_i_i_C(1) - t460 * r_i_i_C(2) + (-t462 * t519 + t478 * t523 + (-t477 * t523 - t495 * t519) * qJD(4)) * pkin(4), t461, 0; -pkin(10) * t577 - t488 * pkin(2) + t460 * r_i_i_C(1) + t459 * r_i_i_C(2) - t476 * qJD(5) + t462 * t511 + t579 * t461 + (-pkin(1) * t522 + t554 * t571) * qJD(1) + (t478 * t519 + (-t477 * t519 + t495 * t523) * qJD(4)) * pkin(4), (t470 * t513 - t485 * t560) * r_i_i_C(1) + (-t470 * t512 - t485 * t559) * r_i_i_C(2) + t470 * t511 - t485 * t556 - t540 * qJD(5) - t489 * pkin(2) + t579 * (t485 * qJD(3) - t489 * t520 + t490 * t569) + ((t490 * t512 + t499 * t559) * r_i_i_C(1) + (t490 * t513 - t499 * t560) * r_i_i_C(2) + t490 * pkin(10) + (t490 * t519 + t499 * t557) * pkin(4)) * t515, -qJD(5) * t475 + t579 * t464 - t582 * t528 - t541 * t581, -t588 * r_i_i_C(1) + t587 * r_i_i_C(2) + (t584 * r_i_i_C(1) + t583 * r_i_i_C(2)) * qJD(4) + (-t464 * t519 + t480 * t523 + (t475 * t523 + t493 * t519) * qJD(4)) * pkin(4), -t581, 0; 0, (t482 * t513 + t496 * t560) * r_i_i_C(1) + (-t482 * t512 + t496 * t559) * r_i_i_C(2) + t482 * t511 + t496 * t556 + (-t579 * (-t535 * qJD(2) + t532 * qJD(3)) + t534 * qJD(5) - pkin(2) * t562 + ((t512 * t561 + t513 * t558) * r_i_i_C(1) + (-t512 * t558 + t513 * t561) * r_i_i_C(2) + pkin(10) * t561 + (t519 * t561 + t521 * t557) * pkin(4)) * t515) * t516, t492 * qJD(5) + t579 * t472 + t541 * t471 + (t535 * t516 + t524 * t549) * t528, (-t472 * t512 + t513 * t546) * r_i_i_C(1) + (-t472 * t513 - t512 * t546) * r_i_i_C(2) + ((-t492 * t513 - t497 * t512) * r_i_i_C(1) + (t492 * t512 - t497 * t513) * r_i_i_C(2)) * qJD(4) + (t523 * t546 - t472 * t519 + (-t492 * t523 - t497 * t519) * qJD(4)) * pkin(4), t471, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:50:09
	% EndTime: 2019-10-10 12:50:13
	% DurationCPUTime: 3.93s
	% Computational Cost: add. (2398->252), mult. (6559->414), div. (0->0), fcn. (7118->16), ass. (0->138)
	t696 = sin(qJ(2));
	t699 = cos(qJ(2));
	t700 = cos(qJ(1));
	t771 = cos(pkin(6));
	t740 = t700 * t771;
	t777 = sin(qJ(1));
	t669 = t696 * t740 + t777 * t699;
	t731 = t771 * t777;
	t714 = t700 * t696 + t699 * t731;
	t656 = t714 * qJD(1) + t669 * qJD(2);
	t723 = t696 * t731;
	t745 = t777 * t696;
	t759 = t700 * t699;
	t657 = -qJD(1) * t723 - qJD(2) * t745 + (qJD(2) * t771 + qJD(1)) * t759;
	t689 = sin(pkin(7));
	t691 = cos(pkin(7));
	t695 = sin(qJ(3));
	t690 = sin(pkin(6));
	t742 = t777 * qJD(1);
	t733 = t690 * t742;
	t778 = cos(qJ(3));
	t668 = -t699 * t740 + t745;
	t763 = t690 * t700;
	t750 = t689 * t763;
	t728 = t778 * t750;
	t749 = t691 * t778;
	t783 = -t668 * t749 - t728;
	t615 = (-qJD(3) * t669 - t656 * t691 + t689 * t733) * t695 + t657 * t778 + t783 * qJD(3);
	t770 = t656 * t689;
	t645 = t691 * t733 + t770;
	t688 = qJ(4) + pkin(13);
	t686 = sin(t688);
	t687 = cos(t688);
	t660 = -t668 * t689 + t691 * t763;
	t678 = t695 * t750;
	t762 = t691 * t695;
	t791 = t668 * t762 - t669 * t778 + t678;
	t782 = t660 * t687 - t686 * t791;
	t607 = qJD(4) * t782 - t615 * t687 - t645 * t686;
	t764 = t690 * t689;
	t736 = t777 * t764;
	t721 = t778 * t736;
	t616 = qJD(1) * t721 + t791 * qJD(3) - t656 * t749 - t657 * t695;
	t693 = sin(qJ(6));
	t697 = cos(qJ(6));
	t798 = t607 * t693 - t616 * t697;
	t797 = t607 * t697 + t616 * t693;
	t627 = -t660 * t686 - t687 * t791;
	t638 = t669 * t695 - t783;
	t796 = t627 * t693 - t638 * t697;
	t795 = t627 * t697 + t638 * t693;
	t794 = -qJD(4) * t627 - t615 * t686 + t645 * t687;
	t715 = t723 - t759;
	t643 = -t715 * t778 + (-t714 * t691 + t736) * t695;
	t712 = t715 * qJD(2);
	t702 = t660 * qJD(1) - t689 * t712;
	t790 = -t643 * qJD(4) + t702;
	t779 = r_i_i_C(3) + pkin(12);
	t781 = t790 * pkin(4);
	t722 = qJD(6) * (t693 * r_i_i_C(1) + t697 * r_i_i_C(2));
	t725 = t697 * r_i_i_C(1) - t693 * r_i_i_C(2) + pkin(5);
	t709 = t714 * t778;
	t780 = -t691 * t709 + t695 * t715 + t721;
	t776 = t689 * pkin(10);
	t694 = sin(qJ(4));
	t774 = t694 * pkin(4);
	t772 = pkin(4) * qJD(4);
	t748 = t691 * t777;
	t662 = t714 * t689 + t690 * t748;
	t769 = t662 * t687;
	t767 = t686 * t689;
	t766 = t687 * t689;
	t765 = t689 * t694;
	t761 = t695 * t696;
	t760 = t695 * t699;
	t758 = qJD(6) * t693;
	t757 = qJD(6) * t697;
	t698 = cos(qJ(4));
	t755 = t698 * t772;
	t754 = t694 * t772;
	t752 = t696 * t764;
	t751 = t690 * t761;
	t747 = t778 * t696;
	t746 = t778 * t699;
	t743 = qJD(2) * t764;
	t741 = t689 * t771;
	t739 = t689 * t755;
	t737 = t691 * t746;
	t735 = t696 * t743;
	t734 = t699 * t743;
	t732 = t695 * t741;
	t729 = t690 * t737;
	t630 = t643 * t687 + t662 * t686;
	t717 = t691 * t760 + t747;
	t659 = t717 * t690 + t732;
	t667 = t771 * t691 - t699 * t764;
	t634 = t659 * t687 + t667 * t686;
	t726 = -t659 * t686 + t667 * t687;
	t724 = t778 * t741;
	t652 = -t668 * t778 - t669 * t762;
	t631 = t652 * t687 + t669 * t767;
	t654 = t715 * t762 - t709;
	t632 = t654 * t687 - t715 * t767;
	t718 = -t691 * t761 + t746;
	t666 = t718 * t690;
	t650 = t666 * t687 + t686 * t752;
	t719 = t668 * t695 - t669 * t749;
	t716 = t691 * t747 + t760;
	t685 = pkin(4) * t698 + pkin(3);
	t708 = -t779 * t686 - t725 * t687 - t685;
	t653 = -t714 * t695 - t715 * t749;
	t706 = t668 * qJD(1) + t712;
	t705 = t687 * t722 + (t725 * t686 - t779 * t687 + t774) * qJD(4);
	t704 = t706 * t695;
	t703 = t706 * t778;
	t692 = -qJ(5) - pkin(11);
	t665 = t716 * t690;
	t658 = -t724 - t729 + t751;
	t655 = t669 * qJD(1) + t714 * qJD(2);
	t647 = (-t717 * qJD(2) - t716 * qJD(3)) * t690;
	t646 = -qJD(2) * t729 - t690 * qJD(3) * t746 + (qJD(3) * t691 + qJD(2)) * t751;
	t636 = qJD(3) * t724 + ((t737 - t761) * qJD(3) + t718 * qJD(2)) * t690;
	t635 = qJD(3) * t732 + (t716 * qJD(2) + t717 * qJD(3)) * t690;
	t625 = t686 * t734 + t647 * t687 + (-t666 * t686 + t687 * t752) * qJD(4);
	t623 = t719 * qJD(3) - t656 * t778 - t657 * t762;
	t622 = t652 * qJD(3) - t656 * t695 + t657 * t749;
	t621 = -t653 * qJD(3) + t655 * t762 + t703;
	t620 = t654 * qJD(3) - t655 * t749 + t704;
	t619 = t726 * qJD(4) + t636 * t687 + t686 * t735;
	t613 = qJD(1) * t678 + t780 * qJD(3) - t655 * t778 + t691 * t704;
	t612 = -qJD(1) * t728 + t643 * qJD(3) - t655 * t695 - t691 * t703;
	t611 = t657 * t767 + t623 * t687 + (-t652 * t686 + t669 * t766) * qJD(4);
	t609 = -t655 * t767 + t621 * t687 + (-t654 * t686 - t715 * t766) * qJD(4);
	t603 = qJD(4) * t769 + t613 * t687 + t790 * t686;
	t602 = t630 * qJD(4) + t613 * t686 - t702 * t687;
	t601 = t603 * t697 + t612 * t693 + (-t630 * t693 - t697 * t780) * qJD(6);
	t600 = -t603 * t693 + t612 * t697 + (-t630 * t697 + t693 * t780) * qJD(6);
	t1 = [t797 * r_i_i_C(1) - t798 * r_i_i_C(2) + t607 * pkin(5) - t615 * t685 - t616 * t692 - t638 * qJD(5) - t657 * pkin(2) - pkin(10) * t770 + t779 * t794 + (t796 * r_i_i_C(1) + t795 * r_i_i_C(2)) * qJD(6) + (-t700 * pkin(1) + (-t777 * pkin(9) - pkin(10) * t748) * t690) * qJD(1) + (-t645 * t694 + (t660 * t698 - t694 * t791) * qJD(4)) * pkin(4), (t609 * t697 + t620 * t693 + (-t632 * t693 + t653 * t697) * qJD(6)) * r_i_i_C(1) + (-t609 * t693 + t620 * t697 + (-t632 * t697 - t653 * t693) * qJD(6)) * r_i_i_C(2) + t609 * pkin(5) + t621 * t685 - t654 * t754 - t620 * t692 + t653 * qJD(5) - t715 * t739 + t706 * pkin(2) + (-t765 * pkin(4) - t776) * t655 + t779 * (t632 * qJD(4) + t621 * t686 + t655 * t766), (t613 * t693 + t643 * t757) * r_i_i_C(1) + (t613 * t697 - t643 * t758) * r_i_i_C(2) - t613 * t692 + t643 * qJD(5) + t708 * t612 - t705 * t780, -t613 * t774 - t662 * t754 + t781 * t698 + (-r_i_i_C(1) * t758 - r_i_i_C(2) * t757) * (-t643 * t686 + t769) + t779 * t603 - t725 * t602, t612, r_i_i_C(1) * t600 - r_i_i_C(2) * t601; -pkin(1) * t742 - t655 * pkin(2) + t603 * pkin(5) + t601 * r_i_i_C(1) + t600 * r_i_i_C(2) - t780 * qJD(5) - t612 * t692 + t613 * t685 + t662 * t755 - t706 * t776 + (pkin(10) * t691 + pkin(9)) * qJD(1) * t763 + t781 * t694 + t779 * t602, (t611 * t697 + t622 * t693) * r_i_i_C(1) + (-t611 * t693 + t622 * t697) * r_i_i_C(2) + t611 * pkin(5) + t623 * t685 - t622 * t692 - t719 * qJD(5) - t656 * pkin(2) + t657 * t776 + t779 * (t631 * qJD(4) + t623 * t686 - t657 * t766) + ((-t631 * t693 - t697 * t719) * r_i_i_C(1) + (-t631 * t697 + t693 * t719) * r_i_i_C(2)) * qJD(6) + (t657 * t765 + (t669 * t689 * t698 - t652 * t694) * qJD(4)) * pkin(4), (t615 * t693 - t757 * t791) * r_i_i_C(1) + (t615 * t697 + t758 * t791) * r_i_i_C(2) - t615 * t692 - t791 * qJD(5) - t708 * t616 + t705 * t638, -t779 * t607 + t782 * t722 + t725 * t794 + (-t615 * t694 + t645 * t698 + (t660 * t694 + t698 * t791) * qJD(4)) * pkin(4), -t616, t798 * r_i_i_C(1) + t797 * r_i_i_C(2) + (-t795 * r_i_i_C(1) + t796 * r_i_i_C(2)) * qJD(6); 0, (t625 * t697 - t646 * t693) * r_i_i_C(1) + (-t625 * t693 - t646 * t697) * r_i_i_C(2) + t625 * pkin(5) + t647 * t685 - t666 * t754 + t646 * t692 + t665 * qJD(5) + t779 * (qJD(4) * t650 + t647 * t686 - t687 * t734) + ((-t650 * t693 + t665 * t697) * r_i_i_C(1) + (-t650 * t697 - t665 * t693) * r_i_i_C(2)) * qJD(6) + (t696 * t739 + (-pkin(2) * t696 + (pkin(10) + t774) * t699 * t689) * qJD(2)) * t690, (t636 * t693 + t659 * t757) * r_i_i_C(1) + (t636 * t697 - t659 * t758) * r_i_i_C(2) - t636 * t692 + t659 * qJD(5) + t708 * t635 + t705 * t658, t779 * t619 - t726 * t722 + t725 * (-qJD(4) * t634 - t636 * t686 + t687 * t735) + (t698 * t735 - t636 * t694 + (-t659 * t698 - t667 * t694) * qJD(4)) * pkin(4), t635, (-t619 * t693 + t635 * t697) * r_i_i_C(1) + (-t619 * t697 - t635 * t693) * r_i_i_C(2) + ((-t634 * t697 - t658 * t693) * r_i_i_C(1) + (t634 * t693 - t658 * t697) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end