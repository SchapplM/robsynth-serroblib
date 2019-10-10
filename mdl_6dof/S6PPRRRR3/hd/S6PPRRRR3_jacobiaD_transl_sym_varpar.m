% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:22
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobiaD_transl_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (23->19), mult. (86->45), div. (0->0), fcn. (88->10), ass. (0->22)
	t111 = sin(pkin(13));
	t117 = cos(pkin(6));
	t126 = t111 * t117;
	t112 = sin(pkin(7));
	t113 = sin(pkin(6));
	t125 = t112 * t113;
	t115 = cos(pkin(13));
	t124 = t115 * t117;
	t116 = cos(pkin(7));
	t118 = sin(qJ(3));
	t123 = t116 * t118;
	t119 = cos(qJ(3));
	t122 = t116 * t119;
	t121 = t118 * t125;
	t120 = t119 * t125;
	t114 = cos(pkin(14));
	t110 = sin(pkin(14));
	t109 = -t110 * t126 + t115 * t114;
	t108 = -t115 * t110 - t114 * t126;
	t107 = t110 * t124 + t111 * t114;
	t106 = -t111 * t110 + t114 * t124;
	t1 = [0, 0, ((-t108 * t123 - t109 * t119 - t111 * t121) * r_i_i_C(1) + (-t108 * t122 + t109 * t118 - t111 * t120) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, 0, ((-t106 * t123 - t107 * t119 + t115 * t121) * r_i_i_C(1) + (-t106 * t122 + t107 * t118 + t115 * t120) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, 0, ((-r_i_i_C(1) * t118 - r_i_i_C(2) * t119) * t117 * t112 + ((-t110 * t119 - t114 * t123) * r_i_i_C(1) + (t110 * t118 - t114 * t122) * r_i_i_C(2)) * t113) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:22
	% EndTime: 2019-10-09 21:22:23
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (257->72), mult. (859->144), div. (0->0), fcn. (972->14), ass. (0->56)
	t318 = sin(pkin(14));
	t322 = sin(pkin(6));
	t323 = cos(pkin(14));
	t329 = sin(qJ(3));
	t331 = cos(qJ(3));
	t326 = cos(pkin(7));
	t343 = t326 * t329;
	t321 = sin(pkin(7));
	t327 = cos(pkin(6));
	t349 = t321 * t327;
	t304 = t322 * (t318 * t331 + t323 * t343) + t329 * t349;
	t324 = cos(pkin(13));
	t319 = sin(pkin(13));
	t353 = t319 * t327;
	t313 = -t318 * t353 + t324 * t323;
	t312 = -t324 * t318 - t323 * t353;
	t350 = t321 * t322;
	t334 = t312 * t326 + t319 * t350;
	t300 = t313 * t331 + t329 * t334;
	t346 = t324 * t327;
	t311 = t318 * t346 + t319 * t323;
	t356 = t311 * t329;
	t355 = t311 * t331;
	t320 = sin(pkin(8));
	t328 = sin(qJ(4));
	t352 = t320 * t328;
	t330 = cos(qJ(4));
	t351 = t320 * t330;
	t348 = t322 * t323;
	t347 = t322 * t326;
	t325 = cos(pkin(8));
	t345 = t325 * t328;
	t344 = t325 * t330;
	t342 = qJD(3) * t329;
	t341 = qJD(3) * t331;
	t340 = t324 * t350;
	t338 = t321 * t341;
	t337 = t326 * t341;
	t336 = t330 * r_i_i_C(1) - t328 * r_i_i_C(2) + pkin(3);
	t310 = -t319 * t318 + t323 * t346;
	t335 = -t310 * t326 + t340;
	t332 = (t328 * r_i_i_C(1) + t330 * r_i_i_C(2)) * t325 + (-pkin(10) - r_i_i_C(3)) * t320;
	t309 = -t321 * t348 + t327 * t326;
	t306 = -t312 * t321 + t319 * t347;
	t305 = -t310 * t321 - t324 * t347;
	t303 = t331 * t349 + (t323 * t326 * t331 - t318 * t329) * t322;
	t302 = t304 * qJD(3);
	t301 = t318 * t322 * t342 - t327 * t338 - t337 * t348;
	t299 = -t313 * t329 + t331 * t334;
	t298 = t310 * t343 - t329 * t340 + t355;
	t297 = -t331 * t335 - t356;
	t296 = t300 * qJD(3);
	t295 = -t319 * t322 * t338 - t312 * t337 + t313 * t342;
	t294 = (t329 * t335 - t355) * qJD(3);
	t293 = -t310 * t337 + (t331 * t340 + t356) * qJD(3);
	t1 = [0, 0, -t336 * t296 + t332 * t295 + ((-t299 * t328 - t300 * t344) * r_i_i_C(1) + (-t299 * t330 + t300 * t345) * r_i_i_C(2)) * qJD(4), (t295 * t328 - t296 * t344) * r_i_i_C(1) + (t295 * t330 + t296 * t345) * r_i_i_C(2) + ((-t299 * t345 - t300 * t330 - t306 * t352) * r_i_i_C(1) + (-t299 * t344 + t300 * t328 - t306 * t351) * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, t336 * t294 + t332 * t293 + ((-t297 * t328 - t298 * t344) * r_i_i_C(1) + (-t297 * t330 + t298 * t345) * r_i_i_C(2)) * qJD(4), (t293 * t328 + t294 * t344) * r_i_i_C(1) + (t293 * t330 - t294 * t345) * r_i_i_C(2) + ((-t297 * t345 - t298 * t330 - t305 * t352) * r_i_i_C(1) + (-t297 * t344 + t298 * t328 - t305 * t351) * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, -t336 * t302 + t332 * t301 + ((-t303 * t328 - t304 * t344) * r_i_i_C(1) + (-t303 * t330 + t304 * t345) * r_i_i_C(2)) * qJD(4), (t301 * t328 - t302 * t344) * r_i_i_C(1) + (t301 * t330 + t302 * t345) * r_i_i_C(2) + ((-t303 * t345 - t304 * t330 - t309 * t352) * r_i_i_C(1) + (-t303 * t344 + t304 * t328 - t309 * t351) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:25
	% EndTime: 2019-10-09 21:22:26
	% DurationCPUTime: 0.94s
	% Computational Cost: add. (1005->131), mult. (3235->244), div. (0->0), fcn. (3808->16), ass. (0->82)
	t510 = sin(pkin(14));
	t514 = sin(pkin(6));
	t515 = cos(pkin(14));
	t522 = sin(qJ(3));
	t525 = cos(qJ(3));
	t518 = cos(pkin(7));
	t545 = t518 * t522;
	t513 = sin(pkin(7));
	t519 = cos(pkin(6));
	t551 = t513 * t519;
	t496 = (t510 * t525 + t515 * t545) * t514 + t522 * t551;
	t521 = sin(qJ(4));
	t524 = cos(qJ(4));
	t495 = t525 * t551 + (t515 * t518 * t525 - t510 * t522) * t514;
	t550 = t514 * t515;
	t501 = -t513 * t550 + t519 * t518;
	t512 = sin(pkin(8));
	t517 = cos(pkin(8));
	t534 = t495 * t517 + t501 * t512;
	t480 = t496 * t524 + t534 * t521;
	t516 = cos(pkin(13));
	t511 = sin(pkin(13));
	t555 = t511 * t519;
	t505 = -t510 * t555 + t516 * t515;
	t504 = -t516 * t510 - t515 * t555;
	t552 = t513 * t514;
	t531 = t504 * t518 + t511 * t552;
	t491 = t505 * t525 + t531 * t522;
	t490 = -t505 * t522 + t531 * t525;
	t549 = t514 * t518;
	t498 = -t504 * t513 + t511 * t549;
	t535 = t490 * t517 + t498 * t512;
	t474 = t491 * t524 + t535 * t521;
	t548 = t516 * t519;
	t502 = -t511 * t510 + t515 * t548;
	t540 = t516 * t552;
	t503 = t510 * t548 + t511 * t515;
	t557 = t503 * t525;
	t489 = t502 * t545 - t522 * t540 + t557;
	t532 = -t502 * t518 + t540;
	t558 = t503 * t522;
	t488 = -t532 * t525 - t558;
	t497 = -t502 * t513 - t516 * t549;
	t536 = t488 * t517 + t497 * t512;
	t472 = t489 * t524 + t536 * t521;
	t562 = r_i_i_C(3) + pkin(11);
	t520 = sin(qJ(5));
	t554 = t512 * t520;
	t523 = cos(qJ(5));
	t553 = t512 * t523;
	t547 = t517 * t521;
	t546 = t517 * t524;
	t544 = qJD(3) * t522;
	t543 = qJD(3) * t525;
	t542 = qJD(5) * t520;
	t541 = qJD(5) * t523;
	t538 = t513 * t543;
	t537 = t518 * t543;
	t533 = t523 * r_i_i_C(1) - t520 * r_i_i_C(2) + pkin(4);
	t477 = t488 * t524 - t489 * t547;
	t478 = t490 * t524 - t491 * t547;
	t483 = t495 * t524 - t496 * t547;
	t529 = qJD(5) * (-t520 * r_i_i_C(1) - t523 * r_i_i_C(2));
	t528 = -t489 * t521 + t536 * t524;
	t527 = -t491 * t521 + t535 * t524;
	t526 = -t496 * t521 + t534 * t524;
	t494 = t496 * qJD(3);
	t493 = t514 * t510 * t544 - t519 * t538 - t537 * t550;
	t492 = -t495 * t512 + t501 * t517;
	t487 = t491 * qJD(3);
	t486 = -t511 * t514 * t538 - t504 * t537 + t505 * t544;
	t485 = (t532 * t522 - t557) * qJD(3);
	t484 = -t502 * t537 + (t525 * t540 + t558) * qJD(3);
	t482 = -t490 * t512 + t498 * t517;
	t481 = -t488 * t512 + t497 * t517;
	t476 = t493 * t547 - t494 * t524 + (-t495 * t521 - t496 * t546) * qJD(4);
	t470 = t526 * qJD(4) - t493 * t524 - t494 * t547;
	t468 = t486 * t547 - t487 * t524 + (-t490 * t521 - t491 * t546) * qJD(4);
	t466 = t484 * t547 + t485 * t524 + (-t488 * t521 - t489 * t546) * qJD(4);
	t464 = t527 * qJD(4) - t486 * t524 - t487 * t547;
	t462 = t528 * qJD(4) - t484 * t524 + t485 * t547;
	t1 = [0, 0, (t468 * t523 - t478 * t542) * r_i_i_C(1) + (-t468 * t520 - t478 * t541) * r_i_i_C(2) + t468 * pkin(4) - t487 * pkin(3) + t562 * (t478 * qJD(4) - t486 * t546 - t487 * t521) + ((-t486 * t520 + t491 * t541) * r_i_i_C(1) + (-t486 * t523 - t491 * t542) * r_i_i_C(2) - t486 * pkin(10)) * t512, t562 * t464 + t527 * t529 + t533 * (-t474 * qJD(4) + t486 * t521 - t487 * t546), (-t464 * t520 + t487 * t553) * r_i_i_C(1) + (-t464 * t523 - t487 * t554) * r_i_i_C(2) + ((-t474 * t523 - t482 * t520) * r_i_i_C(1) + (t474 * t520 - t482 * t523) * r_i_i_C(2)) * qJD(5), 0; 0, 0, (t466 * t523 - t477 * t542) * r_i_i_C(1) + (-t466 * t520 - t477 * t541) * r_i_i_C(2) + t466 * pkin(4) + t485 * pkin(3) + t562 * (t477 * qJD(4) - t484 * t546 + t485 * t521) + ((-t484 * t520 + t489 * t541) * r_i_i_C(1) + (-t484 * t523 - t489 * t542) * r_i_i_C(2) - t484 * pkin(10)) * t512, t562 * t462 + t528 * t529 + t533 * (-t472 * qJD(4) + t484 * t521 + t485 * t546), (-t462 * t520 - t485 * t553) * r_i_i_C(1) + (-t462 * t523 + t485 * t554) * r_i_i_C(2) + ((-t472 * t523 - t481 * t520) * r_i_i_C(1) + (t472 * t520 - t481 * t523) * r_i_i_C(2)) * qJD(5), 0; 0, 0, (t476 * t523 - t483 * t542) * r_i_i_C(1) + (-t476 * t520 - t483 * t541) * r_i_i_C(2) + t476 * pkin(4) - t494 * pkin(3) + t562 * (t483 * qJD(4) - t493 * t546 - t494 * t521) + ((-t493 * t520 + t496 * t541) * r_i_i_C(1) + (-t493 * t523 - t496 * t542) * r_i_i_C(2) - t493 * pkin(10)) * t512, t562 * t470 + t526 * t529 + t533 * (-t480 * qJD(4) + t493 * t521 - t494 * t546), (-t470 * t520 + t494 * t553) * r_i_i_C(1) + (-t470 * t523 - t494 * t554) * r_i_i_C(2) + ((-t480 * t523 - t492 * t520) * r_i_i_C(1) + (t480 * t520 - t492 * t523) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:31
	% EndTime: 2019-10-09 21:22:33
	% DurationCPUTime: 1.77s
	% Computational Cost: add. (3022->198), mult. (9505->345), div. (0->0), fcn. (11498->18), ass. (0->111)
	t678 = sin(qJ(6));
	t682 = cos(qJ(6));
	t698 = qJD(6) * (t678 * r_i_i_C(1) + t682 * r_i_i_C(2));
	t670 = sin(pkin(14));
	t673 = sin(pkin(6));
	t674 = cos(pkin(14));
	t681 = sin(qJ(3));
	t684 = cos(qJ(3));
	t677 = cos(pkin(6));
	t727 = sin(pkin(7));
	t708 = t727 * t677;
	t676 = cos(pkin(7));
	t718 = t676 * t681;
	t658 = (t670 * t684 + t674 * t718) * t673 + t681 * t708;
	t728 = cos(pkin(13));
	t709 = t728 * t674;
	t671 = sin(pkin(13));
	t723 = t671 * t677;
	t665 = -t670 * t723 + t709;
	t710 = t728 * t670;
	t664 = -t674 * t723 - t710;
	t711 = t673 * t727;
	t692 = t664 * t676 + t671 * t711;
	t649 = t665 * t684 + t692 * t681;
	t731 = r_i_i_C(3) + pkin(12);
	t730 = cos(qJ(4));
	t672 = sin(pkin(8));
	t729 = pkin(10) * t672;
	t663 = t671 * t674 + t677 * t710;
	t726 = t663 * t681;
	t725 = t663 * t684;
	t679 = sin(qJ(5));
	t722 = t672 * t679;
	t683 = cos(qJ(5));
	t721 = t672 * t683;
	t720 = t673 * t676;
	t675 = cos(pkin(8));
	t680 = sin(qJ(4));
	t719 = t675 * t680;
	t717 = qJD(3) * t681;
	t716 = qJD(3) * t684;
	t715 = qJD(6) * t678;
	t714 = qJD(6) * t682;
	t713 = t675 * t730;
	t712 = t676 * t716;
	t706 = t727 * t716;
	t662 = -t671 * t670 + t677 * t709;
	t700 = t728 * t711;
	t690 = -t662 * t676 + t700;
	t646 = -t690 * t684 - t726;
	t647 = t662 * t718 - t681 * t700 + t725;
	t691 = -t662 * t727 - t728 * t720;
	t687 = t691 * t672;
	t620 = t647 * t730 + (t646 * t675 + t687) * t680;
	t634 = -t646 * t672 + t691 * t675;
	t608 = t620 * t683 + t634 * t679;
	t704 = -t620 * t679 + t634 * t683;
	t648 = -t665 * t681 + t692 * t684;
	t694 = -t664 * t727 + t671 * t720;
	t689 = t694 * t672;
	t622 = t649 * t730 + (t648 * t675 + t689) * t680;
	t635 = -t648 * t672 + t694 * t675;
	t610 = t622 * t683 + t635 * t679;
	t703 = -t622 * t679 + t635 * t683;
	t657 = t684 * t708 + (t674 * t676 * t684 - t670 * t681) * t673;
	t693 = -t674 * t711 + t677 * t676;
	t688 = t693 * t672;
	t633 = t658 * t730 + (t657 * t675 + t688) * t680;
	t650 = -t657 * t672 + t693 * t675;
	t624 = t633 * t683 + t650 * t679;
	t702 = -t633 * t679 + t650 * t683;
	t701 = t682 * r_i_i_C(1) - t678 * r_i_i_C(2) + pkin(5);
	t628 = t646 * t730 - t647 * t719;
	t615 = t628 * t683 + t647 * t722;
	t630 = t648 * t730 - t649 * t719;
	t616 = t630 * t683 + t649 * t722;
	t637 = t657 * t730 - t658 * t719;
	t631 = t637 * t683 + t658 * t722;
	t697 = -t646 * t680 - t647 * t713;
	t696 = -t648 * t680 - t649 * t713;
	t695 = -t657 * t680 - t658 * t713;
	t686 = -t731 * t679 - t701 * t683 - pkin(4);
	t621 = -t648 * t713 + t649 * t680 - t730 * t689;
	t632 = -t657 * t713 + t658 * t680 - t730 * t688;
	t619 = -t646 * t713 + t647 * t680 - t730 * t687;
	t685 = t683 * t698 + (t701 * t679 - t731 * t683) * qJD(5);
	t656 = t658 * qJD(3);
	t655 = -t677 * t706 + (t670 * t717 - t674 * t712) * t673;
	t645 = t649 * qJD(3);
	t644 = -t671 * t673 * t706 - t664 * t712 + t665 * t717;
	t643 = (t690 * t681 - t725) * qJD(3);
	t642 = -t662 * t712 + (t684 * t700 + t726) * qJD(3);
	t626 = t695 * qJD(4) + t655 * t719 - t656 * t730;
	t625 = t637 * qJD(4) - t655 * t713 - t656 * t680;
	t618 = -t632 * qJD(4) - t655 * t730 - t656 * t719;
	t617 = t633 * qJD(4) - t655 * t680 + t656 * t713;
	t614 = t696 * qJD(4) + t644 * t719 - t645 * t730;
	t613 = t630 * qJD(4) - t644 * t713 - t645 * t680;
	t612 = t697 * qJD(4) + t642 * t719 + t643 * t730;
	t611 = t628 * qJD(4) - t642 * t713 + t643 * t680;
	t606 = -t621 * qJD(4) - t644 * t730 - t645 * t719;
	t605 = t622 * qJD(4) - t644 * t680 + t645 * t713;
	t604 = -t619 * qJD(4) - t642 * t730 + t643 * t719;
	t603 = t620 * qJD(4) - t642 * t680 - t643 * t713;
	t602 = -t655 * t722 + t626 * t683 + (-t637 * t679 + t658 * t721) * qJD(5);
	t600 = t702 * qJD(5) + t618 * t683 + t656 * t722;
	t598 = -t644 * t722 + t614 * t683 + (-t630 * t679 + t649 * t721) * qJD(5);
	t596 = -t642 * t722 + t612 * t683 + (-t628 * t679 + t647 * t721) * qJD(5);
	t594 = t703 * qJD(5) + t606 * t683 + t645 * t722;
	t592 = t704 * qJD(5) + t604 * t683 - t643 * t722;
	t1 = [0, 0, (t598 * t682 + t613 * t678) * r_i_i_C(1) + (-t598 * t678 + t613 * t682) * r_i_i_C(2) + t598 * pkin(5) + t614 * pkin(4) + t613 * pkin(11) - t645 * pkin(3) - t644 * t729 + t731 * (t616 * qJD(5) + t614 * t679 + t644 * t721) + ((-t616 * t678 - t682 * t696) * r_i_i_C(1) + (-t616 * t682 + t678 * t696) * r_i_i_C(2)) * qJD(6), (t606 * t678 + t622 * t714) * r_i_i_C(1) + (t606 * t682 - t622 * t715) * r_i_i_C(2) + t606 * pkin(11) + t686 * t605 + t685 * t621, t731 * t594 - t703 * t698 + t701 * (-t610 * qJD(5) - t606 * t679 + t645 * t721), (-t594 * t678 + t605 * t682) * r_i_i_C(1) + (-t594 * t682 - t605 * t678) * r_i_i_C(2) + ((-t610 * t682 - t621 * t678) * r_i_i_C(1) + (t610 * t678 - t621 * t682) * r_i_i_C(2)) * qJD(6); 0, 0, (t596 * t682 + t611 * t678) * r_i_i_C(1) + (-t596 * t678 + t611 * t682) * r_i_i_C(2) + t596 * pkin(5) + t612 * pkin(4) + t611 * pkin(11) + t643 * pkin(3) - t642 * t729 + t731 * (t615 * qJD(5) + t612 * t679 + t642 * t721) + ((-t615 * t678 - t682 * t697) * r_i_i_C(1) + (-t615 * t682 + t678 * t697) * r_i_i_C(2)) * qJD(6), (t604 * t678 + t620 * t714) * r_i_i_C(1) + (t604 * t682 - t620 * t715) * r_i_i_C(2) + t604 * pkin(11) + t686 * t603 + t685 * t619, t731 * t592 - t704 * t698 + t701 * (-t608 * qJD(5) - t604 * t679 - t643 * t721), (-t592 * t678 + t603 * t682) * r_i_i_C(1) + (-t592 * t682 - t603 * t678) * r_i_i_C(2) + ((-t608 * t682 - t619 * t678) * r_i_i_C(1) + (t608 * t678 - t619 * t682) * r_i_i_C(2)) * qJD(6); 0, 0, (t602 * t682 + t625 * t678) * r_i_i_C(1) + (-t602 * t678 + t625 * t682) * r_i_i_C(2) + t602 * pkin(5) + t626 * pkin(4) + t625 * pkin(11) - t656 * pkin(3) - t655 * t729 + t731 * (t631 * qJD(5) + t626 * t679 + t655 * t721) + ((-t631 * t678 - t682 * t695) * r_i_i_C(1) + (-t631 * t682 + t678 * t695) * r_i_i_C(2)) * qJD(6), (t618 * t678 + t633 * t714) * r_i_i_C(1) + (t618 * t682 - t633 * t715) * r_i_i_C(2) + t618 * pkin(11) + t686 * t617 + t685 * t632, t731 * t600 - t702 * t698 + t701 * (-t624 * qJD(5) - t618 * t679 + t656 * t721), (-t600 * t678 + t617 * t682) * r_i_i_C(1) + (-t600 * t682 - t617 * t678) * r_i_i_C(2) + ((-t624 * t682 - t632 * t678) * r_i_i_C(1) + (t624 * t678 - t632 * t682) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end