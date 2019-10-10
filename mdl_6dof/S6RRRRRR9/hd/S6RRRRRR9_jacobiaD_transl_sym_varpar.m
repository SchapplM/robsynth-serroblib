% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRR9
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
% Datum: 2019-10-10 13:32
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR9_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR9_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR9_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
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
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
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
	% StartTime: 2019-10-10 13:32:12
	% EndTime: 2019-10-10 13:32:13
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
	% StartTime: 2019-10-10 13:32:15
	% EndTime: 2019-10-10 13:32:16
	% DurationCPUTime: 1.33s
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
	% StartTime: 2019-10-10 13:32:21
	% EndTime: 2019-10-10 13:32:24
	% DurationCPUTime: 3.17s
	% Computational Cost: add. (1927->216), mult. (5954->367), div. (0->0), fcn. (6464->14), ass. (0->129)
	t678 = sin(qJ(2));
	t681 = cos(qJ(2));
	t682 = cos(qJ(1));
	t752 = cos(pkin(6));
	t721 = t682 * t752;
	t754 = sin(qJ(1));
	t656 = t678 * t721 + t754 * t681;
	t712 = t752 * t754;
	t693 = t682 * t678 + t681 * t712;
	t643 = t693 * qJD(1) + t656 * qJD(2);
	t704 = t678 * t712;
	t727 = t754 * t678;
	t739 = t682 * t681;
	t644 = -qJD(1) * t704 - qJD(2) * t727 + (qJD(2) * t752 + qJD(1)) * t739;
	t655 = -t681 * t721 + t727;
	t672 = sin(pkin(7));
	t674 = cos(pkin(7));
	t677 = sin(qJ(3));
	t673 = sin(pkin(6));
	t743 = t673 * t682;
	t734 = t672 * t743;
	t755 = cos(qJ(3));
	t709 = t755 * t734;
	t723 = qJD(3) * t755;
	t714 = t674 * t723;
	t724 = qJD(1) * t754;
	t715 = t673 * t724;
	t602 = (-qJD(3) * t656 - t643 * t674 + t672 * t715) * t677 - qJD(3) * t709 + t644 * t755 - t655 * t714;
	t750 = t643 * t672;
	t632 = t674 * t715 + t750;
	t676 = sin(qJ(4));
	t680 = cos(qJ(4));
	t732 = t656 * t755;
	t628 = (t655 * t674 + t734) * t677 - t732;
	t647 = -t655 * t672 + t674 * t743;
	t708 = t628 * t676 - t647 * t680;
	t594 = -t708 * qJD(4) - t602 * t680 - t632 * t676;
	t747 = t672 * t673;
	t718 = t754 * t747;
	t701 = t755 * t718;
	t720 = t677 * t734;
	t731 = t674 * t755;
	t742 = t674 * t677;
	t603 = qJD(1) * t701 - t643 * t731 - t644 * t677 + (t655 * t742 + t720 - t732) * qJD(3);
	t675 = sin(qJ(5));
	t679 = cos(qJ(5));
	t774 = t594 * t675 - t603 * t679;
	t773 = t594 * t679 + t603 * t675;
	t615 = t628 * t680 + t647 * t676;
	t625 = t655 * t731 + t656 * t677 + t709;
	t772 = -t615 * t675 - t625 * t679;
	t771 = t615 * t679 - t625 * t675;
	t770 = t615 * qJD(4) - t602 * t676 + t632 * t680;
	t756 = r_i_i_C(3) + pkin(12);
	t694 = t704 - t739;
	t758 = t677 * t694 + t701;
	t630 = -t694 * t755 + (-t674 * t693 + t718) * t677;
	t730 = t674 * t754;
	t649 = t672 * t693 + t673 * t730;
	t757 = -t630 * t676 + t649 * t680;
	t702 = qJD(5) * (r_i_i_C(1) * t675 + r_i_i_C(2) * t679);
	t753 = pkin(10) * t672;
	t746 = t672 * t676;
	t745 = t672 * t680;
	t744 = t673 * t681;
	t741 = t677 * t678;
	t740 = t677 * t681;
	t738 = qJD(2) * t673;
	t737 = qJD(5) * t675;
	t736 = qJD(5) * t679;
	t735 = t678 * t747;
	t733 = t673 * t741;
	t729 = t755 * t678;
	t728 = t755 * t681;
	t725 = t672 * t738;
	t722 = t672 * t752;
	t719 = t674 * t728;
	t717 = t678 * t725;
	t716 = t681 * t725;
	t713 = t677 * t722;
	t710 = t673 * t719;
	t617 = t630 * t680 + t649 * t676;
	t696 = t674 * t740 + t729;
	t646 = t696 * t673 + t713;
	t654 = -t672 * t744 + t752 * t674;
	t623 = t646 * t680 + t654 * t676;
	t707 = -t646 * t676 + t654 * t680;
	t706 = r_i_i_C(1) * t679 - r_i_i_C(2) * t675 + pkin(4);
	t705 = t755 * t722;
	t638 = -t655 * t755 - t656 * t742;
	t618 = t638 * t680 + t656 * t746;
	t640 = -t693 * t755 + t694 * t742;
	t619 = t640 * t680 - t694 * t746;
	t697 = -t674 * t741 + t728;
	t653 = t697 * t673;
	t641 = t653 * t680 + t676 * t735;
	t699 = t655 * t677 - t656 * t731;
	t698 = t677 * t693 + t694 * t731;
	t695 = t674 * t729 + t740;
	t691 = t694 * qJD(2);
	t688 = -t756 * t676 - t706 * t680 - pkin(3);
	t687 = t680 * t702 + (t706 * t676 - t756 * t680) * qJD(4);
	t686 = t655 * qJD(1) + t691;
	t685 = t686 * t677;
	t684 = t686 * t755;
	t683 = t647 * qJD(1) - t672 * t691;
	t652 = t695 * t673;
	t645 = -t705 - t710 + t733;
	t642 = t656 * qJD(1) + t693 * qJD(2);
	t634 = (-t696 * qJD(2) - t695 * qJD(3)) * t673;
	t633 = -qJD(2) * t710 - t723 * t744 + (qJD(3) * t674 + qJD(2)) * t733;
	t629 = t693 * t731 - t758;
	t621 = qJD(3) * t705 + ((t719 - t741) * qJD(3) + t697 * qJD(2)) * t673;
	t620 = qJD(3) * t713 + (t695 * qJD(2) + t696 * qJD(3)) * t673;
	t612 = t676 * t716 + t634 * t680 + (-t653 * t676 + t680 * t735) * qJD(4);
	t610 = t699 * qJD(3) - t643 * t755 - t644 * t742;
	t609 = t638 * qJD(3) - t643 * t677 + t644 * t731;
	t608 = t698 * qJD(3) + t642 * t742 + t684;
	t607 = t640 * qJD(3) - t642 * t731 + t685;
	t606 = t707 * qJD(4) + t621 * t680 + t676 * t717;
	t600 = qJD(1) * t720 + t758 * qJD(3) - t642 * t755 + t674 * t685 - t693 * t714;
	t599 = -qJD(1) * t709 + t630 * qJD(3) - t642 * t677 - t674 * t684;
	t598 = t644 * t746 + t610 * t680 + (-t638 * t676 + t656 * t745) * qJD(4);
	t596 = -t642 * t746 + t608 * t680 + (-t640 * t676 - t694 * t745) * qJD(4);
	t590 = t757 * qJD(4) + t600 * t680 + t683 * t676;
	t589 = t617 * qJD(4) + t600 * t676 - t683 * t680;
	t588 = t590 * t679 + t599 * t675 + (-t617 * t675 + t629 * t679) * qJD(5);
	t587 = -t590 * t675 + t599 * t679 + (-t617 * t679 - t629 * t675) * qJD(5);
	t1 = [t773 * r_i_i_C(1) - t774 * r_i_i_C(2) + t594 * pkin(4) - t602 * pkin(3) + t603 * pkin(11) - t644 * pkin(2) - pkin(10) * t750 + t756 * t770 + (t772 * r_i_i_C(1) - t771 * r_i_i_C(2)) * qJD(5) + (-t682 * pkin(1) + (-t754 * pkin(9) - pkin(10) * t730) * t673) * qJD(1), (t596 * t679 + t607 * t675 + (-t619 * t675 - t679 * t698) * qJD(5)) * r_i_i_C(1) + (-t596 * t675 + t607 * t679 + (-t619 * t679 + t675 * t698) * qJD(5)) * r_i_i_C(2) + t596 * pkin(4) + t608 * pkin(3) + t607 * pkin(11) + t686 * pkin(2) - t642 * t753 + t756 * (t619 * qJD(4) + t608 * t676 + t642 * t745), (t600 * t675 + t630 * t736) * r_i_i_C(1) + (t600 * t679 - t630 * t737) * r_i_i_C(2) + t600 * pkin(11) + t688 * t599 + t687 * t629, -t706 * t589 + t756 * t590 - t757 * t702, r_i_i_C(1) * t587 - r_i_i_C(2) * t588, 0; -pkin(1) * t724 - t642 * pkin(2) + t600 * pkin(3) + t590 * pkin(4) + t599 * pkin(11) + t588 * r_i_i_C(1) + t587 * r_i_i_C(2) - t686 * t753 + (pkin(10) * t674 + pkin(9)) * qJD(1) * t743 + t756 * t589, (t598 * t679 + t609 * t675) * r_i_i_C(1) + (-t598 * t675 + t609 * t679) * r_i_i_C(2) + t598 * pkin(4) + t610 * pkin(3) + t609 * pkin(11) - t643 * pkin(2) + t644 * t753 + t756 * (t618 * qJD(4) + t610 * t676 - t644 * t745) + ((-t618 * t675 - t679 * t699) * r_i_i_C(1) + (-t618 * t679 + t675 * t699) * r_i_i_C(2)) * qJD(5), (t602 * t675 - t628 * t736) * r_i_i_C(1) + (t602 * t679 + t628 * t737) * r_i_i_C(2) + t602 * pkin(11) - t688 * t603 + t687 * t625, -t594 * t756 - t708 * t702 + t706 * t770, t774 * r_i_i_C(1) + t773 * r_i_i_C(2) + (t771 * r_i_i_C(1) + t772 * r_i_i_C(2)) * qJD(5), 0; 0, (t612 * t679 - t633 * t675) * r_i_i_C(1) + (-t612 * t675 - t633 * t679) * r_i_i_C(2) + t612 * pkin(4) + t634 * pkin(3) - t633 * pkin(11) + t756 * (t641 * qJD(4) + t634 * t676 - t680 * t716) + ((-t641 * t675 + t652 * t679) * r_i_i_C(1) + (-t641 * t679 - t652 * t675) * r_i_i_C(2)) * qJD(5) + (-pkin(2) * t678 + t681 * t753) * t738, (t621 * t675 + t646 * t736) * r_i_i_C(1) + (t621 * t679 - t646 * t737) * r_i_i_C(2) + t621 * pkin(11) + t688 * t620 + t687 * t645, t756 * t606 - t707 * t702 + t706 * (-t623 * qJD(4) - t621 * t676 + t680 * t717), (-t606 * t675 + t620 * t679) * r_i_i_C(1) + (-t606 * t679 - t620 * t675) * r_i_i_C(2) + ((-t623 * t679 - t645 * t675) * r_i_i_C(1) + (t623 * t675 - t645 * t679) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:32:21
	% EndTime: 2019-10-10 13:32:25
	% DurationCPUTime: 3.62s
	% Computational Cost: add. (2727->240), mult. (7763->391), div. (0->0), fcn. (8472->16), ass. (0->150)
	t720 = sin(qJ(2));
	t723 = cos(qJ(2));
	t724 = cos(qJ(1));
	t815 = cos(pkin(6));
	t778 = t724 * t815;
	t819 = sin(qJ(1));
	t693 = t720 * t778 + t819 * t723;
	t755 = t815 * t819;
	t737 = t724 * t720 + t723 * t755;
	t680 = t737 * qJD(1) + t693 * qJD(2);
	t748 = t720 * t755;
	t784 = t819 * t720;
	t800 = t724 * t723;
	t681 = -qJD(1) * t748 - qJD(2) * t784 + (qJD(2) * t815 + qJD(1)) * t800;
	t692 = -t723 * t778 + t784;
	t714 = sin(pkin(7));
	t716 = cos(pkin(7));
	t719 = sin(qJ(3));
	t715 = sin(pkin(6));
	t804 = t715 * t724;
	t791 = t714 * t804;
	t820 = cos(qJ(3));
	t752 = t820 * t791;
	t781 = qJD(3) * t820;
	t757 = t716 * t781;
	t780 = t819 * qJD(1);
	t758 = t715 * t780;
	t639 = (-qJD(3) * t693 - t680 * t716 + t714 * t758) * t719 - qJD(3) * t752 + t681 * t820 - t692 * t757;
	t813 = t680 * t714;
	t669 = t716 * t758 + t813;
	t718 = sin(qJ(4));
	t722 = cos(qJ(4));
	t789 = t693 * t820;
	t665 = (t692 * t716 + t791) * t719 - t789;
	t684 = -t692 * t714 + t716 * t804;
	t751 = t665 * t718 - t684 * t722;
	t631 = -t751 * qJD(4) - t639 * t722 - t669 * t718;
	t788 = t716 * t820;
	t662 = t692 * t788 + t693 * t719 + t752;
	t712 = qJD(5) + qJD(6);
	t835 = -t662 * t712 + t631;
	t808 = t714 * t715;
	t761 = t819 * t808;
	t745 = t820 * t761;
	t763 = t719 * t791;
	t803 = t716 * t719;
	t640 = qJD(1) * t745 - t680 * t788 - t681 * t719 + (t692 * t803 + t763 - t789) * qJD(3);
	t652 = t665 * t722 + t684 * t718;
	t834 = -t652 * t712 + t640;
	t833 = t652 * qJD(4) - t639 * t718 + t669 * t722;
	t816 = r_i_i_C(3) + pkin(13) + pkin(12);
	t717 = sin(qJ(5));
	t817 = t717 * pkin(5);
	t824 = t817 + pkin(11);
	t738 = t748 - t800;
	t823 = t719 * t738 + t745;
	t667 = -t738 * t820 + (-t716 * t737 + t761) * t719;
	t787 = t716 * t819;
	t686 = t714 * t737 + t715 * t787;
	t822 = -t667 * t718 + t686 * t722;
	t713 = qJ(5) + qJ(6);
	t710 = sin(t713);
	t711 = cos(t713);
	t794 = qJD(5) * t817;
	t821 = (t710 * r_i_i_C(1) + t711 * r_i_i_C(2)) * t712 + t794;
	t818 = pkin(10) * t714;
	t810 = t710 * t712;
	t809 = t711 * t712;
	t807 = t714 * t718;
	t806 = t714 * t722;
	t805 = t715 * t723;
	t802 = t719 * t720;
	t801 = t719 * t723;
	t679 = t693 * qJD(1) + t737 * qJD(2);
	t734 = t738 * qJD(2);
	t730 = t692 * qJD(1) + t734;
	t728 = t730 * t820;
	t636 = -qJD(1) * t752 + t667 * qJD(3) - t679 * t719 - t716 * t728;
	t654 = t667 * t722 + t686 * t718;
	t772 = -t654 * t712 + t636;
	t729 = t730 * t719;
	t637 = qJD(1) * t763 + t823 * qJD(3) - t679 * t820 + t716 * t729 - t737 * t757;
	t726 = t684 * qJD(1) - t714 * t734;
	t627 = t822 * qJD(4) + t637 * t722 + t726 * t718;
	t666 = t737 * t788 - t823;
	t777 = t666 * t712 + t627;
	t622 = -t777 * t710 + t772 * t711;
	t623 = t772 * t710 + t777 * t711;
	t799 = t622 * r_i_i_C(1) - t623 * r_i_i_C(2);
	t798 = (t710 * t835 - t711 * t834) * r_i_i_C(1) + (t710 * t834 + t711 * t835) * r_i_i_C(2);
	t786 = t820 * t720;
	t739 = t716 * t786 + t801;
	t740 = t716 * t801 + t786;
	t779 = t714 * t815;
	t756 = t719 * t779;
	t657 = qJD(3) * t756 + (t739 * qJD(2) + t740 * qJD(3)) * t715;
	t683 = t740 * t715 + t756;
	t691 = -t714 * t805 + t815 * t716;
	t660 = t683 * t722 + t691 * t718;
	t765 = t660 * t712 - t657;
	t785 = t820 * t723;
	t741 = -t716 * t802 + t785;
	t749 = t820 * t779;
	t762 = t716 * t785;
	t658 = qJD(3) * t749 + ((t762 - t802) * qJD(3) + t741 * qJD(2)) * t715;
	t750 = -t683 * t718 + t691 * t722;
	t796 = qJD(2) * t715;
	t782 = t714 * t796;
	t760 = t720 * t782;
	t643 = t750 * qJD(4) + t658 * t722 + t718 * t760;
	t753 = t715 * t762;
	t790 = t715 * t802;
	t682 = -t749 - t753 + t790;
	t769 = -t682 * t712 - t643;
	t797 = (t769 * t710 - t765 * t711) * r_i_i_C(1) + (t765 * t710 + t769 * t711) * r_i_i_C(2);
	t721 = cos(qJ(5));
	t795 = qJD(5) * t721;
	t793 = pkin(5) * t795;
	t792 = t720 * t808;
	t742 = t719 * t737 + t738 * t788;
	t645 = t742 * qJD(3) + t679 * t803 + t728;
	t677 = -t737 * t820 + t738 * t803;
	t633 = -t679 * t807 + t645 * t722 + (-t677 * t718 - t738 * t806) * qJD(4);
	t774 = -t712 * t742 + t633;
	t743 = t692 * t719 - t693 * t788;
	t647 = t743 * qJD(3) - t680 * t820 - t681 * t803;
	t675 = -t692 * t820 - t693 * t803;
	t635 = t681 * t807 + t647 * t722 + (-t675 * t718 + t693 * t806) * qJD(4);
	t773 = -t712 * t743 + t635;
	t644 = t677 * qJD(3) - t679 * t788 + t729;
	t656 = t677 * t722 - t738 * t807;
	t768 = -t656 * t712 + t644;
	t646 = t675 * qJD(3) - t680 * t719 + t681 * t788;
	t655 = t675 * t722 + t693 * t807;
	t767 = -t655 * t712 + t646;
	t671 = (-t740 * qJD(2) - t739 * qJD(3)) * t715;
	t690 = t741 * t715;
	t759 = t723 * t782;
	t649 = t718 * t759 + t671 * t722 + (-t690 * t718 + t722 * t792) * qJD(4);
	t689 = t739 * t715;
	t766 = t689 * t712 + t649;
	t670 = -qJD(2) * t753 - t781 * t805 + (qJD(3) * t716 + qJD(2)) * t790;
	t678 = t690 * t722 + t718 * t792;
	t764 = -t678 * t712 - t670;
	t709 = pkin(5) * t721 + pkin(4);
	t747 = t711 * r_i_i_C(1) - t710 * r_i_i_C(2) + t709;
	t731 = -t816 * t718 - t747 * t722 - pkin(3);
	t727 = t821 * t722 + (t747 * t718 - t816 * t722) * qJD(4);
	t626 = t654 * qJD(4) + t637 * t718 - t726 * t722;
	t1 = [-pkin(10) * t813 - t681 * pkin(2) - t639 * pkin(3) + t640 * pkin(11) + t631 * t709 + (r_i_i_C(1) * t835 + r_i_i_C(2) * t834) * t711 + (r_i_i_C(1) * t834 - r_i_i_C(2) * t835) * t710 + t816 * t833 + (-t724 * pkin(1) + (-t819 * pkin(9) - pkin(10) * t787) * t715) * qJD(1) + (t640 * t717 + (-t652 * t717 - t662 * t721) * qJD(5)) * pkin(5), (t768 * t710 + t774 * t711) * r_i_i_C(1) + (-t774 * t710 + t768 * t711) * r_i_i_C(2) + t633 * t709 - t656 * t794 - t742 * t793 + t645 * pkin(3) + t730 * pkin(2) - t679 * t818 + t824 * t644 + t816 * (t656 * qJD(4) + t645 * t718 + t679 * t806), (t637 * t710 + t667 * t809) * r_i_i_C(1) + (t637 * t711 - t667 * t810) * r_i_i_C(2) + t637 * pkin(11) + (t637 * t717 + t667 * t795) * pkin(5) + t731 * t636 + t727 * t666, -t747 * t626 + t816 * t627 - t821 * t822, (-t627 * t717 + t636 * t721 + (-t654 * t721 - t666 * t717) * qJD(5)) * pkin(5) + t799, t799; -pkin(1) * t780 - t679 * pkin(2) + t637 * pkin(3) + t623 * r_i_i_C(1) + t622 * r_i_i_C(2) + t627 * t709 - t654 * t794 + t666 * t793 - t730 * t818 + (pkin(10) * t716 + pkin(9)) * qJD(1) * t804 + t824 * t636 + t816 * t626, t681 * t818 - t680 * pkin(2) + t647 * pkin(3) + t646 * pkin(11) + t635 * t709 + (t773 * r_i_i_C(1) + t767 * r_i_i_C(2)) * t711 + (t767 * r_i_i_C(1) - t773 * r_i_i_C(2)) * t710 + t816 * (t655 * qJD(4) + t647 * t718 - t681 * t806) + (t646 * t717 + (-t655 * t717 - t721 * t743) * qJD(5)) * pkin(5), (t639 * t710 - t665 * t809) * r_i_i_C(1) + (t639 * t711 + t665 * t810) * r_i_i_C(2) + t639 * pkin(11) + (t639 * t717 - t665 * t795) * pkin(5) - t731 * t640 + t727 * t662, -t631 * t816 + t747 * t833 - t821 * t751, (t631 * t717 - t640 * t721 + (t652 * t721 - t662 * t717) * qJD(5)) * pkin(5) + t798, t798; 0, t671 * pkin(3) - t670 * pkin(11) + t649 * t709 + (t766 * r_i_i_C(1) + t764 * r_i_i_C(2)) * t711 + (t764 * r_i_i_C(1) - t766 * r_i_i_C(2)) * t710 + t816 * (t678 * qJD(4) + t671 * t718 - t722 * t759) + (-pkin(2) * t720 + t723 * t818) * t796 + (-t670 * t717 + (-t678 * t717 + t689 * t721) * qJD(5)) * pkin(5), (t658 * t710 + t683 * t809) * r_i_i_C(1) + (t658 * t711 - t683 * t810) * r_i_i_C(2) + t658 * pkin(11) + (t658 * t717 + t683 * t795) * pkin(5) + t731 * t657 + t727 * t682, t816 * t643 - t821 * t750 + t747 * (-qJD(4) * t660 - t658 * t718 + t722 * t760), (-t643 * t717 + t657 * t721 + (-t660 * t721 - t682 * t717) * qJD(5)) * pkin(5) + t797, t797;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end