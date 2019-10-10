% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRP11
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
% Datum: 2019-10-10 13:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP11_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP11_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP11_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:13:34
	% EndTime: 2019-10-10 13:13:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:13:34
	% EndTime: 2019-10-10 13:13:34
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
	% StartTime: 2019-10-10 13:13:35
	% EndTime: 2019-10-10 13:13:35
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
	% StartTime: 2019-10-10 13:13:36
	% EndTime: 2019-10-10 13:13:36
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
	% StartTime: 2019-10-10 13:13:38
	% EndTime: 2019-10-10 13:13:39
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
	% StartTime: 2019-10-10 13:13:45
	% EndTime: 2019-10-10 13:13:48
	% DurationCPUTime: 3.16s
	% Computational Cost: add. (1927->216), mult. (5954->367), div. (0->0), fcn. (6464->14), ass. (0->129)
	t676 = sin(qJ(2));
	t679 = cos(qJ(2));
	t680 = cos(qJ(1));
	t750 = cos(pkin(6));
	t719 = t680 * t750;
	t752 = sin(qJ(1));
	t654 = t676 * t719 + t752 * t679;
	t710 = t750 * t752;
	t691 = t680 * t676 + t679 * t710;
	t641 = t691 * qJD(1) + t654 * qJD(2);
	t702 = t676 * t710;
	t725 = t752 * t676;
	t737 = t680 * t679;
	t642 = -qJD(1) * t702 - qJD(2) * t725 + (qJD(2) * t750 + qJD(1)) * t737;
	t653 = -t679 * t719 + t725;
	t670 = sin(pkin(7));
	t672 = cos(pkin(7));
	t675 = sin(qJ(3));
	t671 = sin(pkin(6));
	t741 = t671 * t680;
	t732 = t670 * t741;
	t753 = cos(qJ(3));
	t707 = t753 * t732;
	t721 = qJD(3) * t753;
	t712 = t672 * t721;
	t722 = qJD(1) * t752;
	t713 = t671 * t722;
	t600 = (-qJD(3) * t654 - t641 * t672 + t670 * t713) * t675 - qJD(3) * t707 + t642 * t753 - t653 * t712;
	t748 = t641 * t670;
	t630 = t672 * t713 + t748;
	t674 = sin(qJ(4));
	t678 = cos(qJ(4));
	t730 = t654 * t753;
	t626 = (t653 * t672 + t732) * t675 - t730;
	t645 = -t653 * t670 + t672 * t741;
	t706 = t626 * t674 - t645 * t678;
	t592 = -t706 * qJD(4) - t600 * t678 - t630 * t674;
	t745 = t670 * t671;
	t716 = t752 * t745;
	t699 = t753 * t716;
	t718 = t675 * t732;
	t729 = t672 * t753;
	t740 = t672 * t675;
	t601 = qJD(1) * t699 - t641 * t729 - t642 * t675 + (t653 * t740 + t718 - t730) * qJD(3);
	t673 = sin(qJ(5));
	t677 = cos(qJ(5));
	t772 = t592 * t673 - t601 * t677;
	t771 = t592 * t677 + t601 * t673;
	t613 = t626 * t678 + t645 * t674;
	t623 = t653 * t729 + t654 * t675 + t707;
	t770 = -t613 * t673 - t623 * t677;
	t769 = t613 * t677 - t623 * t673;
	t768 = t613 * qJD(4) - t600 * t674 + t630 * t678;
	t754 = r_i_i_C(3) + pkin(12);
	t692 = t702 - t737;
	t756 = t675 * t692 + t699;
	t628 = -t692 * t753 + (-t672 * t691 + t716) * t675;
	t728 = t672 * t752;
	t647 = t670 * t691 + t671 * t728;
	t755 = -t628 * t674 + t647 * t678;
	t700 = qJD(5) * (r_i_i_C(1) * t673 + r_i_i_C(2) * t677);
	t751 = pkin(10) * t670;
	t744 = t670 * t674;
	t743 = t670 * t678;
	t742 = t671 * t679;
	t739 = t675 * t676;
	t738 = t675 * t679;
	t736 = qJD(2) * t671;
	t735 = qJD(5) * t673;
	t734 = qJD(5) * t677;
	t733 = t676 * t745;
	t731 = t671 * t739;
	t727 = t753 * t676;
	t726 = t753 * t679;
	t723 = t670 * t736;
	t720 = t670 * t750;
	t717 = t672 * t726;
	t715 = t676 * t723;
	t714 = t679 * t723;
	t711 = t675 * t720;
	t708 = t671 * t717;
	t615 = t628 * t678 + t647 * t674;
	t694 = t672 * t738 + t727;
	t644 = t694 * t671 + t711;
	t652 = -t670 * t742 + t750 * t672;
	t621 = t644 * t678 + t652 * t674;
	t705 = -t644 * t674 + t652 * t678;
	t704 = r_i_i_C(1) * t677 - r_i_i_C(2) * t673 + pkin(4);
	t703 = t753 * t720;
	t636 = -t653 * t753 - t654 * t740;
	t616 = t636 * t678 + t654 * t744;
	t638 = -t691 * t753 + t692 * t740;
	t617 = t638 * t678 - t692 * t744;
	t695 = -t672 * t739 + t726;
	t651 = t695 * t671;
	t639 = t651 * t678 + t674 * t733;
	t697 = t653 * t675 - t654 * t729;
	t696 = t675 * t691 + t692 * t729;
	t693 = t672 * t727 + t738;
	t689 = t692 * qJD(2);
	t686 = -t754 * t674 - t704 * t678 - pkin(3);
	t685 = t678 * t700 + (t704 * t674 - t754 * t678) * qJD(4);
	t684 = t653 * qJD(1) + t689;
	t683 = t684 * t675;
	t682 = t684 * t753;
	t681 = t645 * qJD(1) - t670 * t689;
	t650 = t693 * t671;
	t643 = -t703 - t708 + t731;
	t640 = t654 * qJD(1) + t691 * qJD(2);
	t632 = (-t694 * qJD(2) - t693 * qJD(3)) * t671;
	t631 = -qJD(2) * t708 - t721 * t742 + (qJD(3) * t672 + qJD(2)) * t731;
	t627 = t691 * t729 - t756;
	t619 = qJD(3) * t703 + ((t717 - t739) * qJD(3) + t695 * qJD(2)) * t671;
	t618 = qJD(3) * t711 + (t693 * qJD(2) + t694 * qJD(3)) * t671;
	t610 = t674 * t714 + t632 * t678 + (-t651 * t674 + t678 * t733) * qJD(4);
	t608 = t697 * qJD(3) - t641 * t753 - t642 * t740;
	t607 = t636 * qJD(3) - t641 * t675 + t642 * t729;
	t606 = t696 * qJD(3) + t640 * t740 + t682;
	t605 = t638 * qJD(3) - t640 * t729 + t683;
	t604 = t705 * qJD(4) + t619 * t678 + t674 * t715;
	t598 = qJD(1) * t718 + t756 * qJD(3) - t640 * t753 + t672 * t683 - t691 * t712;
	t597 = -qJD(1) * t707 + t628 * qJD(3) - t640 * t675 - t672 * t682;
	t596 = t642 * t744 + t608 * t678 + (-t636 * t674 + t654 * t743) * qJD(4);
	t594 = -t640 * t744 + t606 * t678 + (-t638 * t674 - t692 * t743) * qJD(4);
	t588 = t755 * qJD(4) + t598 * t678 + t681 * t674;
	t587 = t615 * qJD(4) + t598 * t674 - t681 * t678;
	t586 = t588 * t677 + t597 * t673 + (-t615 * t673 + t627 * t677) * qJD(5);
	t585 = -t588 * t673 + t597 * t677 + (-t615 * t677 - t627 * t673) * qJD(5);
	t1 = [t771 * r_i_i_C(1) - t772 * r_i_i_C(2) + t592 * pkin(4) - t600 * pkin(3) + t601 * pkin(11) - t642 * pkin(2) - pkin(10) * t748 + t754 * t768 + (t770 * r_i_i_C(1) - t769 * r_i_i_C(2)) * qJD(5) + (-t680 * pkin(1) + (-t752 * pkin(9) - pkin(10) * t728) * t671) * qJD(1), (t594 * t677 + t605 * t673 + (-t617 * t673 - t677 * t696) * qJD(5)) * r_i_i_C(1) + (-t594 * t673 + t605 * t677 + (-t617 * t677 + t673 * t696) * qJD(5)) * r_i_i_C(2) + t594 * pkin(4) + t606 * pkin(3) + t605 * pkin(11) + t684 * pkin(2) - t640 * t751 + t754 * (t617 * qJD(4) + t606 * t674 + t640 * t743), (t598 * t673 + t628 * t734) * r_i_i_C(1) + (t598 * t677 - t628 * t735) * r_i_i_C(2) + t598 * pkin(11) + t686 * t597 + t685 * t627, -t704 * t587 + t754 * t588 - t755 * t700, r_i_i_C(1) * t585 - r_i_i_C(2) * t586, 0; -pkin(1) * t722 - t640 * pkin(2) + t598 * pkin(3) + t588 * pkin(4) + t597 * pkin(11) + t586 * r_i_i_C(1) + t585 * r_i_i_C(2) - t684 * t751 + (pkin(10) * t672 + pkin(9)) * qJD(1) * t741 + t754 * t587, (t596 * t677 + t607 * t673) * r_i_i_C(1) + (-t596 * t673 + t607 * t677) * r_i_i_C(2) + t596 * pkin(4) + t608 * pkin(3) + t607 * pkin(11) - t641 * pkin(2) + t642 * t751 + t754 * (t616 * qJD(4) + t608 * t674 - t642 * t743) + ((-t616 * t673 - t677 * t697) * r_i_i_C(1) + (-t616 * t677 + t673 * t697) * r_i_i_C(2)) * qJD(5), (t600 * t673 - t626 * t734) * r_i_i_C(1) + (t600 * t677 + t626 * t735) * r_i_i_C(2) + t600 * pkin(11) - t686 * t601 + t685 * t623, -t592 * t754 - t706 * t700 + t704 * t768, t772 * r_i_i_C(1) + t771 * r_i_i_C(2) + (t769 * r_i_i_C(1) + t770 * r_i_i_C(2)) * qJD(5), 0; 0, (t610 * t677 - t631 * t673) * r_i_i_C(1) + (-t610 * t673 - t631 * t677) * r_i_i_C(2) + t610 * pkin(4) + t632 * pkin(3) - t631 * pkin(11) + t754 * (t639 * qJD(4) + t632 * t674 - t678 * t714) + ((-t639 * t673 + t650 * t677) * r_i_i_C(1) + (-t639 * t677 - t650 * t673) * r_i_i_C(2)) * qJD(5) + (-pkin(2) * t676 + t679 * t751) * t736, (t619 * t673 + t644 * t734) * r_i_i_C(1) + (t619 * t677 - t644 * t735) * r_i_i_C(2) + t619 * pkin(11) + t686 * t618 + t685 * t643, t754 * t604 - t705 * t700 + t704 * (-t621 * qJD(4) - t619 * t674 + t678 * t715), (-t604 * t673 + t618 * t677) * r_i_i_C(1) + (-t604 * t677 - t618 * t673) * r_i_i_C(2) + ((-t621 * t677 - t643 * t673) * r_i_i_C(1) + (t621 * t673 - t643 * t677) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:13:44
	% EndTime: 2019-10-10 13:13:48
	% DurationCPUTime: 3.45s
	% Computational Cost: add. (2401->227), mult. (7333->371), div. (0->0), fcn. (8019->14), ass. (0->136)
	t682 = cos(qJ(2));
	t683 = cos(qJ(1));
	t773 = cos(pkin(6));
	t740 = t683 * t773;
	t679 = sin(qJ(2));
	t777 = sin(qJ(1));
	t745 = t777 * t679;
	t655 = -t682 * t740 + t745;
	t674 = cos(pkin(7));
	t678 = sin(qJ(3));
	t656 = t679 * t740 + t777 * t682;
	t778 = cos(qJ(3));
	t750 = t656 * t778;
	t672 = sin(pkin(7));
	t673 = sin(pkin(6));
	t763 = t673 * t683;
	t751 = t672 * t763;
	t629 = (t655 * t674 + t751) * t678 - t750;
	t647 = -t655 * t672 + t674 * t763;
	t677 = sin(qJ(4));
	t681 = cos(qJ(4));
	t615 = t629 * t681 + t647 * t677;
	t729 = t778 * t751;
	t749 = t674 * t778;
	t782 = -t655 * t749 - t729;
	t626 = t656 * t678 - t782;
	t676 = sin(qJ(5));
	t680 = cos(qJ(5));
	t792 = -t615 * t676 - t626 * t680;
	t726 = t615 * t680 - t626 * t676;
	t731 = t773 * t777;
	t696 = t683 * t679 + t682 * t731;
	t643 = t696 * qJD(1) + t656 * qJD(2);
	t717 = t679 * t731;
	t759 = t683 * t682;
	t644 = -qJD(1) * t717 - qJD(2) * t745 + (qJD(2) * t773 + qJD(1)) * t759;
	t742 = t777 * qJD(1);
	t733 = t673 * t742;
	t602 = (-qJD(3) * t656 - t643 * t674 + t672 * t733) * t678 + t644 * t778 + t782 * qJD(3);
	t769 = t643 * t672;
	t698 = t674 * t733 + t769;
	t791 = t615 * qJD(4) - t602 * t677 + t698 * t681;
	t722 = t629 * t677 - t647 * t681;
	t592 = t722 * qJD(4) + t602 * t681 + t698 * t677;
	t779 = pkin(5) + r_i_i_C(1);
	t714 = t680 * r_i_i_C(2) + t779 * t676;
	t706 = pkin(11) + t714;
	t699 = t714 * qJD(5);
	t774 = r_i_i_C(3) + qJ(6) + pkin(12);
	t697 = t717 - t759;
	t764 = t673 * t672;
	t736 = t777 * t764;
	t631 = -t697 * t778 + (-t696 * t674 + t736) * t678;
	t748 = t674 * t777;
	t649 = t696 * t672 + t673 * t748;
	t616 = -t631 * t677 + t649 * t681;
	t691 = t696 * t778;
	t710 = t778 * t736;
	t781 = -t674 * t691 + t678 * t697 + t710;
	t739 = t678 * t751;
	t762 = t674 * t678;
	t780 = qJD(1) * t710 - t643 * t749 - t644 * t678 + (t655 * t762 + t739 - t750) * qJD(3);
	t776 = pkin(10) * t672;
	t642 = t656 * qJD(1) + t696 * qJD(2);
	t694 = t697 * qJD(2);
	t688 = t655 * qJD(1) + t694;
	t686 = t688 * t778;
	t599 = -qJD(1) * t729 + t631 * qJD(3) - t642 * t678 - t674 * t686;
	t772 = t599 * t676;
	t640 = t697 * t762 - t691;
	t687 = t688 * t678;
	t607 = t640 * qJD(3) - t642 * t749 + t687;
	t771 = t607 * t676;
	t766 = t672 * t677;
	t765 = t672 * t681;
	t761 = t678 * t679;
	t760 = t678 * t682;
	t758 = qJD(2) * t673;
	t757 = qJD(5) * t676;
	t756 = qJD(5) * t680;
	t755 = pkin(5) * t757;
	t754 = pkin(5) * t756;
	t753 = t679 * t764;
	t752 = t673 * t761;
	t747 = t778 * t679;
	t746 = t778 * t682;
	t743 = t672 * t758;
	t741 = t672 * t773;
	t737 = t674 * t746;
	t735 = t679 * t743;
	t734 = t682 * t743;
	t732 = t678 * t741;
	t730 = t673 * t737;
	t728 = -t592 * t676 - t680 * t780;
	t702 = -t674 * t761 + t746;
	t718 = t778 * t741;
	t621 = qJD(3) * t718 + ((t737 - t761) * qJD(3) + t702 * qJD(2)) * t673;
	t701 = t674 * t760 + t747;
	t646 = t701 * t673 + t732;
	t654 = t773 * t674 - t682 * t764;
	t719 = -t646 * t677 + t654 * t681;
	t606 = t719 * qJD(4) + t621 * t681 + t677 * t735;
	t700 = t674 * t747 + t760;
	t620 = qJD(3) * t732 + (t700 * qJD(2) + t701 * qJD(3)) * t673;
	t727 = -t606 * t676 + t620 * t680;
	t623 = t646 * t681 + t654 * t677;
	t645 = -t718 - t730 + t752;
	t723 = -t623 * t680 - t645 * t676;
	t617 = t631 * t681 + t649 * t677;
	t671 = t680 * pkin(5) + pkin(4);
	t716 = t680 * r_i_i_C(1) - t676 * r_i_i_C(2) + t671;
	t638 = -t655 * t778 - t656 * t762;
	t713 = -t638 * t677 + t656 * t765;
	t618 = t638 * t681 + t656 * t766;
	t712 = -t640 * t677 - t697 * t765;
	t619 = t640 * t681 - t697 * t766;
	t653 = t702 * t673;
	t705 = -t653 * t677 + t681 * t753;
	t641 = t653 * t681 + t677 * t753;
	t703 = t655 * t678 - t656 * t749;
	t690 = -t774 * t677 - t716 * t681 - pkin(3);
	t600 = qJD(1) * t739 + t781 * qJD(3) - t642 * t778 + t674 * t687;
	t684 = t647 * qJD(1) - t672 * t694;
	t590 = t616 * qJD(4) + t600 * t681 + t684 * t677;
	t587 = -t590 * t676 + t599 * t680 + (-t617 * t680 + t676 * t781) * qJD(5);
	t639 = -t696 * t678 - t697 * t749;
	t685 = -t677 * qJD(6) + t681 * t699 + (t716 * t677 - t774 * t681) * qJD(4);
	t652 = t700 * t673;
	t634 = (-t701 * qJD(2) - t700 * qJD(3)) * t673;
	t610 = t703 * qJD(3) - t643 * t778 - t644 * t762;
	t608 = -t639 * qJD(3) + t642 * t762 + t686;
	t605 = t623 * qJD(4) + t621 * t677 - t681 * t735;
	t596 = t712 * qJD(4) + t608 * t681 - t642 * t766;
	t589 = t617 * qJD(4) + t600 * t677 - t684 * t681;
	t588 = t590 * t680 + t772 + (-t617 * t676 - t680 * t781) * qJD(5);
	t1 = [t722 * qJD(6) - t602 * pkin(3) - t644 * pkin(2) - pkin(10) * t769 - t716 * t592 + t706 * t780 + t774 * t791 + (-t683 * pkin(1) + (-t777 * pkin(9) - pkin(10) * t748) * t673) * qJD(1) + (-t726 * r_i_i_C(2) + t779 * t792) * qJD(5), (t596 * t680 + t771 + (-t619 * t676 + t639 * t680) * qJD(5)) * r_i_i_C(1) + (-t596 * t676 + t607 * t680 + (-t619 * t680 - t639 * t676) * qJD(5)) * r_i_i_C(2) + t596 * t671 - t619 * t755 - t712 * qJD(6) + pkin(5) * t771 + t639 * t754 + t608 * pkin(3) + t607 * pkin(11) + t688 * pkin(2) - t642 * t776 + t774 * (t619 * qJD(4) + t608 * t677 + t642 * t765), (t600 * t680 - t631 * t757) * r_i_i_C(2) + t600 * pkin(11) + t690 * t599 - t685 * t781 + t779 * (t600 * t676 + t631 * t756), qJD(6) * t617 - t589 * t716 + t590 * t774 - t616 * t699, -r_i_i_C(2) * t588 + t779 * t587, t589; -pkin(1) * t742 - t642 * pkin(2) + t600 * pkin(3) + pkin(5) * t772 + t599 * pkin(11) + t588 * r_i_i_C(1) + t587 * r_i_i_C(2) - t616 * qJD(6) + t590 * t671 - t617 * t755 - t781 * t754 - t688 * t776 + (pkin(10) * t674 + pkin(9)) * qJD(1) * t763 + t774 * t589, -t713 * qJD(6) + t610 * pkin(3) - t643 * pkin(2) + t644 * t776 + t716 * (qJD(4) * t713 + t610 * t681 + t644 * t766) + t706 * (qJD(3) * t638 - t643 * t678 + t644 * t749) + t774 * (qJD(4) * t618 + t610 * t677 - t644 * t765) + ((-t618 * t680 + t676 * t703) * r_i_i_C(2) + t779 * (-t618 * t676 - t680 * t703)) * qJD(5), (t602 * t680 + t629 * t757) * r_i_i_C(2) + t602 * pkin(11) - t690 * t780 + t685 * t626 + t779 * (t602 * t676 - t629 * t756), -qJD(6) * t615 + t592 * t774 - t699 * t722 + t716 * t791, t728 * r_i_i_C(1) + (-t592 * t680 + t676 * t780) * r_i_i_C(2) + (t726 * r_i_i_C(1) + t792 * r_i_i_C(2)) * qJD(5) + (qJD(5) * t726 + t728) * pkin(5), -t791; 0, -t705 * qJD(6) + t634 * pkin(3) + t716 * (qJD(4) * t705 + t634 * t681 + t677 * t734) - t706 * (-qJD(2) * t730 - t673 * qJD(3) * t746 + (qJD(3) * t674 + qJD(2)) * t752) + t774 * (qJD(4) * t641 + t634 * t677 - t681 * t734) + (-pkin(2) * t679 + t682 * t776) * t758 + ((-t641 * t680 - t652 * t676) * r_i_i_C(2) + t779 * (-t641 * t676 + t652 * t680)) * qJD(5), (t621 * t680 - t646 * t757) * r_i_i_C(2) + t621 * pkin(11) + t690 * t620 + t685 * t645 + t779 * (t621 * t676 + t646 * t756), qJD(6) * t623 - t605 * t716 + t606 * t774 - t699 * t719, t727 * r_i_i_C(1) + (-t606 * t680 - t620 * t676) * r_i_i_C(2) + (t723 * r_i_i_C(1) + (t623 * t676 - t645 * t680) * r_i_i_C(2)) * qJD(5) + (qJD(5) * t723 + t727) * pkin(5), t605;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end