% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR14
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
% Datum: 2019-10-10 12:54
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR14_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR14_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR14_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
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
	% StartTime: 2019-10-10 12:54:13
	% EndTime: 2019-10-10 12:54:13
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
	% StartTime: 2019-10-10 12:54:14
	% EndTime: 2019-10-10 12:54:15
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
	% StartTime: 2019-10-10 12:54:16
	% EndTime: 2019-10-10 12:54:18
	% DurationCPUTime: 1.22s
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
	% StartTime: 2019-10-10 12:54:22
	% EndTime: 2019-10-10 12:54:24
	% DurationCPUTime: 2.18s
	% Computational Cost: add. (1597->175), mult. (4956->293), div. (0->0), fcn. (5291->14), ass. (0->112)
	t609 = sin(qJ(2));
	t610 = sin(qJ(1));
	t680 = cos(pkin(6));
	t650 = t610 * t680;
	t644 = t609 * t650;
	t613 = cos(qJ(2));
	t614 = cos(qJ(1));
	t661 = t614 * t613;
	t663 = t610 * t609;
	t581 = -qJD(1) * t644 - qJD(2) * t663 + (qJD(2) * t680 + qJD(1)) * t661;
	t649 = t614 * t680;
	t589 = -t613 * t649 + t663;
	t590 = t609 * t649 + t610 * t613;
	t608 = sin(qJ(3));
	t612 = cos(qJ(3));
	t623 = t614 * t609 + t613 * t650;
	t580 = t623 * qJD(1) + t590 * qJD(2);
	t603 = sin(pkin(7));
	t606 = cos(pkin(7));
	t604 = sin(pkin(6));
	t660 = qJD(1) * t610;
	t655 = t604 * t660;
	t625 = -t580 * t606 + t603 * t655;
	t657 = qJD(3) * t612;
	t652 = t604 * t657;
	t645 = t603 * t652;
	t658 = qJD(3) * t606;
	t553 = (-qJD(3) * t590 + t625) * t608 - (t589 * t658 - t581) * t612 - t614 * t645;
	t607 = sin(qJ(4));
	t611 = cos(qJ(4));
	t678 = t580 * t603;
	t626 = t606 * t655 + t678;
	t669 = t604 * t614;
	t633 = t589 * t606 + t603 * t669;
	t570 = -t590 * t612 + t633 * t608;
	t584 = -t589 * t603 + t606 * t669;
	t690 = t570 * t611 + t584 * t607;
	t695 = t690 * qJD(4) - t553 * t607 + t626 * t611;
	t691 = t570 * t607 - t584 * t611;
	t694 = t691 * qJD(4) + t553 * t611 + t626 * t607;
	t681 = r_i_i_C(3) + qJ(5);
	t685 = pkin(10) * t606 + pkin(9);
	t624 = t644 - t661;
	t670 = t604 * t610;
	t675 = t623 * t606;
	t632 = t603 * t670 - t675;
	t572 = t632 * t608 - t612 * t624;
	t586 = t603 * t623 + t606 * t670;
	t684 = -t572 * t607 + t586 * t611;
	t602 = sin(pkin(13));
	t605 = cos(pkin(13));
	t637 = t605 * r_i_i_C(1) - t602 * r_i_i_C(2) + pkin(4);
	t620 = t681 * t607 + t637 * t611 + pkin(3);
	t683 = pkin(10) * t603;
	t674 = t624 * t608;
	t673 = t603 * t607;
	t672 = t603 * t611;
	t671 = t604 * t603;
	t668 = t606 * t608;
	t667 = t606 * t612;
	t666 = t608 * t609;
	t665 = t608 * t613;
	t664 = t609 * t612;
	t662 = t612 * t613;
	t659 = qJD(2) * t604;
	t656 = t609 * t671;
	t654 = qJD(1) * t669;
	t653 = t613 * t659;
	t651 = t603 * t680;
	t648 = t603 * t654;
	t647 = qJD(2) * t656;
	t646 = t603 * t653;
	t643 = qJD(3) * t651;
	t639 = t572 * t611 + t586 * t607;
	t629 = t606 * t665 + t664;
	t583 = t629 * t604 + t608 * t651;
	t588 = t680 * t606 - t613 * t671;
	t638 = t583 * t611 + t588 * t607;
	t636 = t602 * r_i_i_C(1) + t605 * r_i_i_C(2) + pkin(11);
	t577 = -t589 * t612 - t590 * t668;
	t635 = -t577 * t607 + t590 * t672;
	t578 = -t612 * t623 + t624 * t668;
	t634 = -t578 * t607 - t624 * t672;
	t631 = t606 * t662 - t666;
	t630 = -t606 * t664 - t665;
	t628 = -t606 * t666 + t662;
	t587 = t628 * t604;
	t627 = -t587 * t607 + t611 * t656;
	t622 = t624 * qJD(2);
	t619 = t607 * qJD(5) + (-t637 * t607 + t681 * t611) * qJD(4);
	t618 = t589 * qJD(1) + t622;
	t617 = t618 * t608;
	t616 = t618 * t612;
	t554 = t570 * qJD(3) - t581 * t608 + t625 * t612;
	t615 = t584 * qJD(1) - t603 * t622;
	t579 = t590 * qJD(1) + t623 * qJD(2);
	t575 = (-t629 * qJD(2) + t630 * qJD(3)) * t604;
	t574 = -t653 * t667 - t613 * t652 + (qJD(2) + t658) * t604 * t666;
	t565 = t612 * t643 + (t628 * qJD(2) + t631 * qJD(3)) * t604;
	t563 = t627 * qJD(4) + t575 * t611 + t607 * t646;
	t561 = -t581 * t668 - t580 * t612 + (t589 * t608 - t590 * t667) * qJD(3);
	t560 = t577 * qJD(3) - t580 * t608 + t581 * t667;
	t559 = t616 + t579 * t668 + (t608 * t623 + t624 * t667) * qJD(3);
	t558 = t578 * qJD(3) - t579 * t667 + t617;
	t556 = t638 * qJD(4) + t565 * t607 - t611 * t647;
	t551 = qJD(3) * t674 - t579 * t612 + t606 * t617 + t608 * t648 + t610 * t645 - t657 * t675;
	t550 = t572 * qJD(3) - t579 * t608 - t606 * t616 - t612 * t648;
	t549 = t635 * qJD(4) + t561 * t611 + t581 * t673;
	t547 = t634 * qJD(4) + t559 * t611 - t579 * t673;
	t541 = t684 * qJD(4) + t551 * t611 + t615 * t607;
	t540 = t639 * qJD(4) + t551 * t607 - t615 * t611;
	t1 = [(t554 * t602 - t605 * t694) * r_i_i_C(1) + (t554 * t605 + t602 * t694) * r_i_i_C(2) - t694 * pkin(4) + t691 * qJD(5) - t553 * pkin(3) + t554 * pkin(11) - t581 * pkin(2) - pkin(10) * t678 + t681 * t695 + (-t614 * pkin(1) - t685 * t670) * qJD(1), (t547 * t605 + t558 * t602) * r_i_i_C(1) + (-t547 * t602 + t558 * t605) * r_i_i_C(2) + t547 * pkin(4) - t634 * qJD(5) + t559 * pkin(3) + t558 * pkin(11) + t618 * pkin(2) - t579 * t683 + t681 * (t579 * t672 + t559 * t607 + (t578 * t611 - t624 * t673) * qJD(4)), t636 * t551 + t619 * (t632 * t612 + t674) - t620 * t550, t639 * qJD(5) - t637 * t540 + t681 * t541, t540, 0; (t541 * t605 + t550 * t602) * r_i_i_C(1) + (-t541 * t602 + t550 * t605) * r_i_i_C(2) + t541 * pkin(4) - t684 * qJD(5) + t551 * pkin(3) + t550 * pkin(11) - t579 * pkin(2) - t618 * t683 - pkin(1) * t660 + t685 * t654 + t681 * t540, (t549 * t605 + t560 * t602) * r_i_i_C(1) + (-t549 * t602 + t560 * t605) * r_i_i_C(2) + t549 * pkin(4) - t635 * qJD(5) + t561 * pkin(3) + t560 * pkin(11) - t580 * pkin(2) + t581 * t683 + t681 * (-t581 * t672 + t561 * t607 + (t577 * t611 + t590 * t673) * qJD(4)), t636 * t553 + t619 * (-t590 * t608 - t633 * t612) + t620 * t554, -t690 * qJD(5) + t637 * t695 + t681 * t694, -t695, 0; 0, (t563 * t605 - t574 * t602) * r_i_i_C(1) + (-t563 * t602 - t574 * t605) * r_i_i_C(2) + t563 * pkin(4) - t627 * qJD(5) + t575 * pkin(3) - t574 * pkin(11) + t681 * (-t611 * t646 + t575 * t607 + (t587 * t611 + t607 * t656) * qJD(4)) + (-pkin(2) * t609 + t613 * t683) * t659, t636 * t565 + t619 * (t631 * t604 + t612 * t651) + t620 * (-t608 * t643 + (t630 * qJD(2) - t629 * qJD(3)) * t604), t638 * qJD(5) + t681 * (t607 * t647 + t565 * t611 + (-t583 * t607 + t588 * t611) * qJD(4)) - t637 * t556, t556, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:23
	% EndTime: 2019-10-10 12:54:26
	% DurationCPUTime: 3.15s
	% Computational Cost: add. (2277->209), mult. (6623->346), div. (0->0), fcn. (7210->16), ass. (0->131)
	t688 = cos(qJ(2));
	t689 = cos(qJ(1));
	t765 = cos(pkin(6));
	t737 = t689 * t765;
	t686 = sin(qJ(2));
	t769 = sin(qJ(1));
	t742 = t769 * t686;
	t657 = -t688 * t737 + t742;
	t682 = cos(pkin(7));
	t685 = sin(qJ(3));
	t658 = t686 * t737 + t769 * t688;
	t770 = cos(qJ(3));
	t747 = t658 * t770;
	t680 = sin(pkin(7));
	t681 = sin(pkin(6));
	t756 = t681 * t689;
	t749 = t680 * t756;
	t631 = (t657 * t682 + t749) * t685 - t747;
	t649 = -t657 * t680 + t682 * t756;
	t684 = sin(qJ(4));
	t687 = cos(qJ(4));
	t617 = t631 * t687 + t649 * t684;
	t723 = t770 * t749;
	t746 = t682 * t770;
	t628 = t657 * t746 + t658 * t685 + t723;
	t678 = pkin(13) + qJ(6);
	t676 = sin(t678);
	t677 = cos(t678);
	t784 = -t617 * t676 - t628 * t677;
	t783 = t617 * t677 - t628 * t676;
	t727 = t765 * t769;
	t700 = t689 * t686 + t688 * t727;
	t645 = t700 * qJD(1) + t658 * qJD(2);
	t718 = t686 * t727;
	t752 = t689 * t688;
	t646 = -qJD(1) * t718 - qJD(2) * t742 + (qJD(2) * t765 + qJD(1)) * t752;
	t739 = t770 * qJD(3);
	t729 = t682 * t739;
	t738 = t769 * qJD(1);
	t730 = t681 * t738;
	t604 = (-qJD(3) * t658 - t645 * t682 + t680 * t730) * t685 - qJD(3) * t723 + t646 * t770 - t657 * t729;
	t763 = t645 * t680;
	t702 = t682 * t730 + t763;
	t782 = t617 * qJD(4) - t604 * t684 + t702 * t687;
	t722 = t631 * t684 - t649 * t687;
	t594 = t722 * qJD(4) + t604 * t687 + t702 * t684;
	t773 = pkin(11) + sin(pkin(13)) * pkin(5);
	t766 = r_i_i_C(3) + pkin(12) + qJ(5);
	t701 = t718 - t752;
	t760 = t680 * t681;
	t733 = t769 * t760;
	t711 = t770 * t733;
	t772 = t685 * t701 + t711;
	t633 = -t701 * t770 + (-t682 * t700 + t733) * t685;
	t745 = t682 * t769;
	t651 = t680 * t700 + t681 * t745;
	t618 = -t633 * t684 + t651 * t687;
	t725 = t676 * r_i_i_C(1) + t677 * r_i_i_C(2);
	t712 = qJD(6) * t725;
	t703 = t725 + t773;
	t735 = t685 * t749;
	t755 = t682 * t685;
	t771 = qJD(1) * t711 - t645 * t746 - t646 * t685 + (t657 * t755 + t735 - t747) * qJD(3);
	t768 = pkin(10) * t680;
	t759 = t680 * t684;
	t758 = t680 * t687;
	t757 = t681 * t688;
	t754 = t685 * t686;
	t753 = t685 * t688;
	t751 = qJD(2) * t681;
	t750 = t686 * t760;
	t748 = t681 * t754;
	t744 = t770 * t686;
	t743 = t770 * t688;
	t740 = t680 * t751;
	t736 = t765 * t680;
	t734 = t682 * t743;
	t732 = t686 * t740;
	t731 = t688 * t740;
	t728 = t685 * t736;
	t726 = t677 * r_i_i_C(1) - t676 * r_i_i_C(2);
	t724 = t681 * t734;
	t619 = t633 * t687 + t651 * t684;
	t705 = t682 * t753 + t744;
	t648 = t705 * t681 + t728;
	t656 = -t680 * t757 + t765 * t682;
	t625 = t648 * t687 + t656 * t684;
	t720 = -t648 * t684 + t656 * t687;
	t719 = t770 * t736;
	t675 = cos(pkin(13)) * pkin(5) + pkin(4);
	t717 = t675 + t726;
	t640 = -t657 * t770 - t658 * t755;
	t716 = -t640 * t684 + t658 * t758;
	t620 = t640 * t687 + t658 * t759;
	t642 = -t700 * t770 + t701 * t755;
	t715 = -t642 * t684 - t701 * t758;
	t621 = t642 * t687 - t701 * t759;
	t713 = qJD(6) * t726;
	t706 = -t682 * t754 + t743;
	t655 = t706 * t681;
	t710 = -t655 * t684 + t687 * t750;
	t643 = t655 * t687 + t684 * t750;
	t708 = t657 * t685 - t658 * t746;
	t707 = t685 * t700 + t701 * t746;
	t704 = t682 * t744 + t753;
	t698 = t701 * qJD(2);
	t695 = -t766 * t684 - t717 * t687 - pkin(3);
	t694 = t657 * qJD(1) + t698;
	t693 = t694 * t685;
	t692 = -t684 * qJD(5) + t687 * t712 + (t717 * t684 - t766 * t687) * qJD(4);
	t691 = t694 * t770;
	t690 = t649 * qJD(1) - t680 * t698;
	t654 = t704 * t681;
	t647 = -t719 - t724 + t748;
	t644 = t658 * qJD(1) + t700 * qJD(2);
	t636 = (-t705 * qJD(2) - t704 * qJD(3)) * t681;
	t632 = t700 * t746 - t772;
	t623 = qJD(3) * t719 + ((t734 - t754) * qJD(3) + t706 * qJD(2)) * t681;
	t622 = qJD(3) * t728 + (t704 * qJD(2) + t705 * qJD(3)) * t681;
	t612 = t708 * qJD(3) - t645 * t770 - t646 * t755;
	t610 = t707 * qJD(3) + t644 * t755 + t691;
	t608 = t720 * qJD(4) + t623 * t687 + t684 * t732;
	t607 = t625 * qJD(4) + t623 * t684 - t687 * t732;
	t602 = qJD(1) * t735 + t772 * qJD(3) - t644 * t770 + t682 * t693 - t700 * t729;
	t601 = -qJD(1) * t723 + t633 * qJD(3) - t644 * t685 - t682 * t691;
	t598 = t715 * qJD(4) + t610 * t687 - t644 * t759;
	t592 = t618 * qJD(4) + t602 * t687 + t690 * t684;
	t591 = t619 * qJD(4) + t602 * t684 - t690 * t687;
	t590 = t592 * t677 + t601 * t676 + (-t619 * t676 + t632 * t677) * qJD(6);
	t589 = -t592 * t676 + t601 * t677 + (-t619 * t677 - t632 * t676) * qJD(6);
	t1 = [t722 * qJD(5) - t604 * pkin(3) - t646 * pkin(2) - pkin(10) * t763 - t717 * t594 + t703 * t771 + t766 * t782 + (t784 * r_i_i_C(1) - t783 * r_i_i_C(2)) * qJD(6) + (-t689 * pkin(1) + (-t769 * pkin(9) - pkin(10) * t745) * t681) * qJD(1), (t598 * t677 + (-t621 * t676 - t677 * t707) * qJD(6)) * r_i_i_C(1) + (-t598 * t676 + (-t621 * t677 + t676 * t707) * qJD(6)) * r_i_i_C(2) + t598 * t675 - t715 * qJD(5) + t610 * pkin(3) + t694 * pkin(2) - t644 * t768 + t703 * (t642 * qJD(3) - t644 * t746 + t693) + t766 * (t621 * qJD(4) + t610 * t684 + t644 * t758), t695 * t601 + t703 * t602 + t692 * t632 + t633 * t713, t619 * qJD(5) - t717 * t591 + t766 * t592 - t618 * t712, t591, t589 * r_i_i_C(1) - t590 * r_i_i_C(2); -pkin(1) * t738 - t644 * pkin(2) + t602 * pkin(3) + t590 * r_i_i_C(1) + t589 * r_i_i_C(2) - t618 * qJD(5) + t592 * t675 - t694 * t768 + (pkin(10) * t682 + pkin(9)) * qJD(1) * t756 + t773 * t601 + t766 * t591, -t716 * qJD(5) + t612 * pkin(3) - t645 * pkin(2) + t646 * t768 + t717 * (t716 * qJD(4) + t612 * t687 + t646 * t759) + t703 * (t640 * qJD(3) - t645 * t685 + t646 * t746) + t766 * (t620 * qJD(4) + t612 * t684 - t646 * t758) + ((-t620 * t676 - t677 * t708) * r_i_i_C(1) + (-t620 * t677 + t676 * t708) * r_i_i_C(2)) * qJD(6), t703 * t604 + t692 * t628 - t631 * t713 - t695 * t771, -qJD(5) * t617 + t766 * t594 - t722 * t712 + t717 * t782, -t782, (-t594 * t676 - t677 * t771) * r_i_i_C(1) + (-t594 * t677 + t676 * t771) * r_i_i_C(2) + (t783 * r_i_i_C(1) + t784 * r_i_i_C(2)) * qJD(6); 0, -t710 * qJD(5) + t636 * pkin(3) + t717 * (t710 * qJD(4) + t636 * t687 + t684 * t731) - t703 * (-qJD(2) * t724 - t739 * t757 + (qJD(3) * t682 + qJD(2)) * t748) + t766 * (t643 * qJD(4) + t636 * t684 - t687 * t731) + ((-t643 * t676 + t654 * t677) * r_i_i_C(1) + (-t643 * t677 - t654 * t676) * r_i_i_C(2)) * qJD(6) + (-pkin(2) * t686 + t688 * t768) * t751, t695 * t622 + t703 * t623 + t692 * t647 + t648 * t713, t625 * qJD(5) - t717 * t607 + t766 * t608 - t720 * t712, t607, (-t608 * t676 + t622 * t677) * r_i_i_C(1) + (-t608 * t677 - t622 * t676) * r_i_i_C(2) + ((-t625 * t677 - t647 * t676) * r_i_i_C(1) + (t625 * t676 - t647 * t677) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end