% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:16
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP12_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP12_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP12_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:47
	% EndTime: 2019-10-10 13:15:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:47
	% EndTime: 2019-10-10 13:15:47
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
	% StartTime: 2019-10-10 13:15:48
	% EndTime: 2019-10-10 13:15:48
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
	% StartTime: 2019-10-10 13:15:49
	% EndTime: 2019-10-10 13:15:49
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
	% StartTime: 2019-10-10 13:15:51
	% EndTime: 2019-10-10 13:15:52
	% DurationCPUTime: 1.26s
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
	% StartTime: 2019-10-10 13:15:58
	% EndTime: 2019-10-10 13:16:01
	% DurationCPUTime: 3.16s
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
	% StartTime: 2019-10-10 13:16:02
	% EndTime: 2019-10-10 13:16:06
	% DurationCPUTime: 4.04s
	% Computational Cost: add. (3260->241), mult. (9997->391), div. (0->0), fcn. (11050->14), ass. (0->145)
	t774 = sin(qJ(2));
	t777 = cos(qJ(2));
	t778 = cos(qJ(1));
	t861 = cos(pkin(6));
	t830 = t778 * t861;
	t864 = sin(qJ(1));
	t752 = t774 * t830 + t864 * t777;
	t821 = t861 * t864;
	t791 = t778 * t774 + t777 * t821;
	t739 = t791 * qJD(1) + t752 * qJD(2);
	t804 = t774 * t821;
	t836 = t864 * t774;
	t848 = t778 * t777;
	t740 = -qJD(1) * t804 - qJD(2) * t836 + (qJD(2) * t861 + qJD(1)) * t848;
	t751 = -t777 * t830 + t836;
	t768 = sin(pkin(7));
	t770 = cos(pkin(7));
	t773 = sin(qJ(3));
	t769 = sin(pkin(6));
	t853 = t769 * t778;
	t843 = t768 * t853;
	t865 = cos(qJ(3));
	t816 = t865 * t843;
	t832 = qJD(3) * t865;
	t823 = t770 * t832;
	t833 = qJD(1) * t864;
	t824 = t769 * t833;
	t697 = (-qJD(3) * t752 - t739 * t770 + t768 * t824) * t773 - qJD(3) * t816 + t740 * t865 - t751 * t823;
	t860 = t739 * t768;
	t728 = t770 * t824 + t860;
	t772 = sin(qJ(4));
	t776 = cos(qJ(4));
	t841 = t752 * t865;
	t724 = (t751 * t770 + t843) * t773 - t841;
	t743 = -t751 * t768 + t770 * t853;
	t808 = t724 * t772 - t743 * t776;
	t681 = -t808 * qJD(4) - t697 * t776 - t728 * t772;
	t857 = t768 * t769;
	t827 = t864 * t857;
	t802 = t865 * t827;
	t829 = t773 * t843;
	t840 = t770 * t865;
	t852 = t770 * t773;
	t698 = qJD(1) * t802 - t739 * t840 - t740 * t773 + (t751 * t852 + t829 - t841) * qJD(3);
	t771 = sin(qJ(5));
	t775 = cos(qJ(5));
	t710 = t724 * t776 + t743 * t772;
	t721 = t751 * t840 + t752 * t773 + t816;
	t881 = t710 * t775 - t721 * t771;
	t886 = t881 * qJD(5) + t681 * t771 - t698 * t775;
	t882 = t710 * t771 + t721 * t775;
	t885 = t882 * qJD(5) - t681 * t775 - t698 * t771;
	t880 = t710 * qJD(4) - t697 * t772 + t728 * t776;
	t867 = r_i_i_C(1) + pkin(5);
	t866 = r_i_i_C(2) + pkin(12);
	t862 = r_i_i_C(3) + qJ(6);
	t792 = t804 - t848;
	t868 = t773 * t792 + t802;
	t793 = t862 * t771 + t867 * t775 + pkin(4);
	t863 = pkin(10) * t768;
	t839 = t770 * t864;
	t745 = t768 * t791 + t769 * t839;
	t859 = t745 * t776;
	t856 = t768 * t772;
	t855 = t768 * t776;
	t854 = t769 * t777;
	t851 = t771 * t776;
	t850 = t773 * t774;
	t849 = t773 * t777;
	t847 = qJD(2) * t769;
	t846 = qJD(4) * t772;
	t845 = qJD(5) * t776;
	t844 = t774 * t857;
	t842 = t769 * t850;
	t838 = t865 * t774;
	t837 = t865 * t777;
	t834 = t768 * t847;
	t831 = t768 * t861;
	t828 = t770 * t837;
	t826 = t774 * t834;
	t825 = t777 * t834;
	t822 = t773 * t831;
	t738 = t752 * qJD(1) + t791 * qJD(2);
	t786 = t792 * qJD(2);
	t782 = t751 * qJD(1) + t786;
	t781 = t782 * t773;
	t695 = qJD(1) * t829 + t868 * qJD(3) - t738 * t865 + t770 * t781 - t791 * t823;
	t725 = t791 * t840 - t868;
	t820 = t725 * t845 + t695;
	t819 = t721 * t845 + t697;
	t797 = -t770 * t850 + t837;
	t805 = t865 * t831;
	t717 = qJD(3) * t805 + ((t828 - t850) * qJD(3) + t797 * qJD(2)) * t769;
	t817 = t769 * t828;
	t741 = -t805 - t817 + t842;
	t818 = t741 * t845 + t717;
	t726 = -t792 * t865 + (-t770 * t791 + t827) * t773;
	t712 = t726 * t776 + t745 * t772;
	t813 = t712 * t775 + t725 * t771;
	t812 = -t712 * t771 + t725 * t775;
	t734 = -t751 * t865 - t752 * t852;
	t713 = t734 * t776 + t752 * t856;
	t799 = t751 * t773 - t752 * t840;
	t811 = -t713 * t771 - t775 * t799;
	t736 = -t791 * t865 + t792 * t852;
	t714 = t736 * t776 - t792 * t856;
	t798 = t773 * t791 + t792 * t840;
	t810 = -t714 * t771 - t775 * t798;
	t796 = t770 * t849 + t838;
	t742 = t796 * t769 + t822;
	t750 = -t768 * t854 + t861 * t770;
	t719 = t742 * t776 + t750 * t772;
	t809 = t719 * t775 + t741 * t771;
	t749 = t797 * t769;
	t737 = t749 * t776 + t772 * t844;
	t795 = t770 * t838 + t849;
	t748 = t795 * t769;
	t807 = -t737 * t771 + t748 * t775;
	t806 = -t742 * t772 + t750 * t776;
	t801 = -t776 * pkin(4) - t866 * t772 - pkin(3);
	t794 = qJD(4) * (pkin(4) * t772 - t866 * t776);
	t780 = t782 * t865;
	t694 = -qJD(1) * t816 + t726 * qJD(3) - t738 * t773 - t770 * t780;
	t790 = qJD(5) * t726 - t694 * t776 + t725 * t846;
	t789 = -qJD(5) * t724 + t698 * t776 + t721 * t846;
	t716 = qJD(3) * t822 + (t795 * qJD(2) + t796 * qJD(3)) * t769;
	t788 = qJD(5) * t742 - t716 * t776 + t741 * t846;
	t783 = qJD(6) * t771 + (-t867 * t771 + t862 * t775) * qJD(5);
	t779 = t743 * qJD(1) - t768 * t786;
	t730 = (-t796 * qJD(2) - t795 * qJD(3)) * t769;
	t729 = -qJD(2) * t817 - t832 * t854 + (qJD(3) * t770 + qJD(2)) * t842;
	t707 = t772 * t825 + t730 * t776 + (-t749 * t772 + t776 * t844) * qJD(4);
	t705 = t799 * qJD(3) - t739 * t865 - t740 * t852;
	t704 = t734 * qJD(3) - t739 * t773 + t740 * t840;
	t703 = t798 * qJD(3) + t738 * t852 + t780;
	t702 = t736 * qJD(3) - t738 * t840 + t781;
	t701 = t806 * qJD(4) + t717 * t776 + t772 * t826;
	t687 = t740 * t856 + t705 * t776 + (-t734 * t772 + t752 * t855) * qJD(4);
	t685 = -t738 * t856 + t703 * t776 + (-t736 * t772 - t792 * t855) * qJD(4);
	t682 = t809 * qJD(5) + t701 * t771 - t716 * t775;
	t677 = qJD(4) * t859 + t695 * t776 - t726 * t846 + t779 * t772;
	t676 = t712 * qJD(4) + t695 * t772 - t779 * t776;
	t663 = t812 * qJD(5) + t677 * t775 + t694 * t771;
	t662 = t813 * qJD(5) + t677 * t771 - t694 * t775;
	t1 = [t882 * qJD(6) + t681 * pkin(4) - t697 * pkin(3) + t698 * pkin(11) - t740 * pkin(2) - pkin(10) * t860 + t866 * t880 - t867 * t885 + t862 * t886 + (-t778 * pkin(1) + (-t864 * pkin(9) - pkin(10) * t839) * t769) * qJD(1), t782 * pkin(2) + t703 * pkin(3) + t685 * pkin(4) + t702 * pkin(11) - t810 * qJD(6) - t738 * t863 + t866 * (t714 * qJD(4) + t703 * t772 + t738 * t855) + t867 * (t810 * qJD(5) + t685 * t775 + t702 * t771) + t862 * (t685 * t771 - t702 * t775 + (t714 * t775 - t771 * t798) * qJD(5)), -(t725 * t851 + t726 * t775) * qJD(6) + t695 * pkin(11) + t867 * (t820 * t771 + t790 * t775) + t862 * (t790 * t771 - t820 * t775) + t725 * t794 + t801 * t694, t866 * t677 + t783 * (-t726 * t772 + t859) - t793 * t676, t813 * qJD(6) - t867 * t662 + t862 * t663, t662; -pkin(1) * t833 - t738 * pkin(2) + t695 * pkin(3) + t677 * pkin(4) + t694 * pkin(11) - t812 * qJD(6) - t782 * t863 + (pkin(10) * t770 + pkin(9)) * qJD(1) * t853 + t866 * t676 + t867 * t663 + t862 * t662, -t811 * qJD(6) + t687 * pkin(4) + t705 * pkin(3) + t704 * pkin(11) - t739 * pkin(2) + t740 * t863 + t866 * (t713 * qJD(4) + t705 * t772 - t740 * t855) + t867 * (t811 * qJD(5) + t687 * t775 + t704 * t771) + t862 * (t687 * t771 - t704 * t775 + (t713 * t775 - t771 * t799) * qJD(5)), -(t721 * t851 - t724 * t775) * qJD(6) + t697 * pkin(11) + t867 * (t819 * t771 + t789 * t775) + t862 * (t789 * t771 - t819 * t775) + t721 * t794 - t801 * t698, -t681 * t866 + t783 * t808 + t793 * t880, -t881 * qJD(6) + t862 * t885 + t867 * t886, -t886; 0, -t807 * qJD(6) + t707 * pkin(4) + t730 * pkin(3) - t729 * pkin(11) + t866 * (t737 * qJD(4) + t730 * t772 - t776 * t825) + t867 * (t807 * qJD(5) + t707 * t775 - t729 * t771) + t862 * (t707 * t771 + t729 * t775 + (t737 * t775 + t748 * t771) * qJD(5)) + (-pkin(2) * t774 + t777 * t863) * t847, -(t741 * t851 + t742 * t775) * qJD(6) + t717 * pkin(11) + t867 * (t818 * t771 + t788 * t775) + t862 * (t788 * t771 - t818 * t775) + t741 * t794 + t801 * t716, t866 * t701 + t783 * t806 + t793 * (-t719 * qJD(4) - t717 * t772 + t776 * t826), t809 * qJD(6) + t862 * (t701 * t775 + t716 * t771 + (-t719 * t771 + t741 * t775) * qJD(5)) - t867 * t682, t682;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end