% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRPRR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(12));
	t58 = sin(pkin(12));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:21
	% EndTime: 2019-10-09 22:35:21
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (78->41), mult. (285->88), div. (0->0), fcn. (301->10), ass. (0->35)
	t259 = sin(pkin(12));
	t262 = cos(pkin(12));
	t266 = sin(qJ(2));
	t264 = cos(pkin(6));
	t268 = cos(qJ(2));
	t283 = t264 * t268;
	t253 = -t259 * t266 + t262 * t283;
	t260 = sin(pkin(7));
	t261 = sin(pkin(6));
	t287 = t260 * t261;
	t263 = cos(pkin(7));
	t265 = sin(qJ(3));
	t286 = t263 * t265;
	t267 = cos(qJ(3));
	t285 = t263 * t267;
	t284 = t264 * t266;
	t282 = t265 * t266;
	t281 = t265 * t268;
	t280 = t266 * t267;
	t279 = t267 * t268;
	t277 = qJD(3) * t260 * t264;
	t276 = -t253 * t263 + t262 * t287;
	t274 = t259 * t283 + t262 * t266;
	t275 = -t259 * t287 + t263 * t274;
	t254 = t259 * t268 + t262 * t284;
	t273 = t259 * t284 - t262 * t268;
	t272 = -t263 * t279 + t282;
	t271 = -t263 * t280 - t281;
	t270 = -t263 * t281 - t280;
	t269 = t263 * t282 - t279;
	t252 = t273 * qJD(2);
	t251 = t274 * qJD(2);
	t250 = t254 * qJD(2);
	t249 = t253 * qJD(2);
	t1 = [0, t251 * t286 + t252 * t267 + (t265 * t274 + t273 * t285) * qJD(3), t252 * t285 + t251 * t265 + (t275 * t265 + t267 * t273) * qJD(3), 0, 0, 0; 0, -t249 * t286 - t250 * t267 + (-t253 * t265 - t254 * t285) * qJD(3), -t250 * t285 - t249 * t265 + (-t254 * t267 + t276 * t265) * qJD(3), 0, 0, 0; 0, (t270 * qJD(2) + t271 * qJD(3)) * t261, -t265 * t277 + (t271 * qJD(2) + t270 * qJD(3)) * t261, 0, 0, 0; 0, t251 * t285 - t252 * t265 + (t267 * t274 - t273 * t286) * qJD(3), -t252 * t286 + t251 * t267 + (-t265 * t273 + t275 * t267) * qJD(3), 0, 0, 0; 0, -t249 * t285 + t250 * t265 + (-t253 * t267 + t254 * t286) * qJD(3), t250 * t286 - t249 * t267 + (t254 * t265 + t276 * t267) * qJD(3), 0, 0, 0; 0, (t272 * qJD(2) + t269 * qJD(3)) * t261, -t267 * t277 + (t269 * qJD(2) + t272 * qJD(3)) * t261, 0, 0, 0; 0, -t251 * t260, 0, 0, 0, 0; 0, t249 * t260, 0, 0, 0, 0; 0, qJD(2) * t268 * t287, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:22
	% EndTime: 2019-10-09 22:35:22
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (124->49), mult. (453->105), div. (0->0), fcn. (477->12), ass. (0->46)
	t413 = sin(pkin(12));
	t417 = cos(pkin(12));
	t421 = sin(qJ(2));
	t419 = cos(pkin(6));
	t423 = cos(qJ(2));
	t439 = t419 * t423;
	t406 = -t413 * t421 + t417 * t439;
	t412 = sin(pkin(13));
	t414 = sin(pkin(7));
	t446 = t412 * t414;
	t415 = sin(pkin(6));
	t444 = t414 * t415;
	t416 = cos(pkin(13));
	t443 = t414 * t416;
	t418 = cos(pkin(7));
	t420 = sin(qJ(3));
	t442 = t418 * t420;
	t422 = cos(qJ(3));
	t441 = t418 * t422;
	t440 = t419 * t421;
	t438 = t420 * t421;
	t437 = t420 * t423;
	t436 = t421 * t422;
	t435 = t422 * t423;
	t433 = qJD(3) * t414 * t419;
	t432 = qJD(2) * t423 * t444;
	t431 = -t406 * t418 + t417 * t444;
	t429 = t413 * t439 + t417 * t421;
	t430 = t413 * t444 - t418 * t429;
	t407 = t413 * t423 + t417 * t440;
	t428 = t413 * t440 - t417 * t423;
	t427 = t418 * t435 - t438;
	t426 = -t418 * t436 - t437;
	t425 = -t418 * t437 - t436;
	t424 = -t418 * t438 + t435;
	t405 = t428 * qJD(2);
	t404 = t429 * qJD(2);
	t403 = t407 * qJD(2);
	t402 = t406 * qJD(2);
	t401 = (t425 * qJD(2) + t426 * qJD(3)) * t415;
	t400 = -t420 * t433 + (t426 * qJD(2) + t425 * qJD(3)) * t415;
	t399 = t404 * t442 + t405 * t422 + (t420 * t429 + t428 * t441) * qJD(3);
	t398 = -t402 * t442 - t403 * t422 + (-t406 * t420 - t407 * t441) * qJD(3);
	t397 = t405 * t441 + t404 * t420 + (-t430 * t420 + t422 * t428) * qJD(3);
	t396 = -t403 * t441 - t402 * t420 + (-t407 * t422 + t431 * t420) * qJD(3);
	t1 = [0, t399 * t416 - t404 * t446, t397 * t416, 0, 0, 0; 0, t398 * t416 + t402 * t446, t396 * t416, 0, 0, 0; 0, t401 * t416 + t412 * t432, t400 * t416, 0, 0, 0; 0, -t399 * t412 - t404 * t443, -t397 * t412, 0, 0, 0; 0, -t398 * t412 + t402 * t443, -t396 * t412, 0, 0, 0; 0, -t401 * t412 + t416 * t432, -t400 * t412, 0, 0, 0; 0, -t404 * t441 + t405 * t420 + (-t422 * t429 + t428 * t442) * qJD(3), t405 * t442 - t404 * t422 + (t420 * t428 + t430 * t422) * qJD(3), 0, 0, 0; 0, t402 * t441 - t403 * t420 + (t406 * t422 - t407 * t442) * qJD(3), -t403 * t442 + t402 * t422 + (-t407 * t420 - t431 * t422) * qJD(3), 0, 0, 0; 0, (t427 * qJD(2) + t424 * qJD(3)) * t415, t422 * t433 + (t424 * qJD(2) + t427 * qJD(3)) * t415, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:22
	% EndTime: 2019-10-09 22:35:23
	% DurationCPUTime: 0.63s
	% Computational Cost: add. (329->88), mult. (921->177), div. (0->0), fcn. (1021->12), ass. (0->69)
	t496 = sin(pkin(12));
	t499 = cos(pkin(12));
	t503 = sin(qJ(2));
	t501 = cos(pkin(6));
	t505 = cos(qJ(2));
	t526 = t501 * t505;
	t487 = -t496 * t503 + t499 * t526;
	t502 = sin(qJ(3));
	t504 = cos(qJ(3));
	t527 = t501 * t503;
	t510 = t496 * t527 - t499 * t505;
	t500 = cos(pkin(7));
	t511 = t496 * t526 + t499 * t503;
	t497 = sin(pkin(7));
	t498 = sin(pkin(6));
	t532 = t497 * t498;
	t512 = t496 * t532 - t500 * t511;
	t473 = t512 * t502 - t504 * t510;
	t488 = t496 * t505 + t499 * t527;
	t513 = -t487 * t500 + t499 * t532;
	t538 = -t488 * t504 + t513 * t502;
	t495 = pkin(13) + qJ(5);
	t493 = sin(t495);
	t535 = t493 * t497;
	t494 = cos(t495);
	t534 = t494 * t497;
	t531 = t497 * t501;
	t530 = t498 * t500;
	t529 = t500 * t502;
	t528 = t500 * t504;
	t525 = t502 * t503;
	t524 = t502 * t505;
	t523 = t503 * t504;
	t522 = t504 * t505;
	t521 = qJD(5) * t493;
	t520 = qJD(5) * t494;
	t519 = t503 * t532;
	t517 = qJD(2) * t532;
	t516 = qJD(3) * t531;
	t515 = t503 * t517;
	t514 = t505 * t517;
	t475 = t487 * t504 - t488 * t529;
	t476 = -t504 * t511 + t510 * t529;
	t509 = t500 * t522 - t525;
	t508 = -t500 * t523 - t524;
	t507 = t500 * t524 + t523;
	t506 = -t500 * t525 + t522;
	t470 = -t488 * t502 - t513 * t504;
	t472 = t502 * t510 + t512 * t504;
	t486 = t501 * t500 - t505 * t532;
	t485 = t510 * qJD(2);
	t484 = t511 * qJD(2);
	t483 = t488 * qJD(2);
	t482 = t487 * qJD(2);
	t481 = t506 * t498;
	t480 = t496 * t530 + t497 * t511;
	t479 = -t487 * t497 - t499 * t530;
	t478 = t507 * t498 + t502 * t531;
	t477 = t509 * t498 + t504 * t531;
	t474 = (-t507 * qJD(2) + t508 * qJD(3)) * t498;
	t469 = t504 * t516 + (t506 * qJD(2) + t509 * qJD(3)) * t498;
	t468 = -t502 * t516 + (t508 * qJD(2) - t507 * qJD(3)) * t498;
	t467 = t484 * t529 + t485 * t504 + (t502 * t511 + t510 * t528) * qJD(3);
	t466 = -t482 * t529 - t483 * t504 + (-t487 * t502 - t488 * t528) * qJD(3);
	t465 = t472 * qJD(3) - t484 * t504 + t485 * t529;
	t464 = -t473 * qJD(3) + t484 * t502 + t485 * t528;
	t463 = t470 * qJD(3) + t482 * t504 - t483 * t529;
	t462 = t538 * qJD(3) - t482 * t502 - t483 * t528;
	t1 = [0, -t484 * t535 + t467 * t494 + (-t476 * t493 - t510 * t534) * qJD(5), t464 * t494 - t472 * t521, 0, -t485 * t534 - t465 * t493 + (-t473 * t494 - t480 * t493) * qJD(5), 0; 0, t482 * t535 + t466 * t494 + (-t475 * t493 + t488 * t534) * qJD(5), t462 * t494 - t470 * t521, 0, t483 * t534 - t463 * t493 + (-t479 * t493 + t494 * t538) * qJD(5), 0; 0, t493 * t514 + t474 * t494 + (-t481 * t493 + t494 * t519) * qJD(5), t468 * t494 - t477 * t521, 0, t494 * t515 - t469 * t493 + (-t478 * t494 - t486 * t493) * qJD(5), 0; 0, -t484 * t534 - t467 * t493 + (-t476 * t494 + t510 * t535) * qJD(5), -t464 * t493 - t472 * t520, 0, t485 * t535 - t465 * t494 + (t473 * t493 - t480 * t494) * qJD(5), 0; 0, t482 * t534 - t466 * t493 + (-t475 * t494 - t488 * t535) * qJD(5), -t462 * t493 - t470 * t520, 0, -t483 * t535 - t463 * t494 + (-t479 * t494 - t493 * t538) * qJD(5), 0; 0, t494 * t514 - t474 * t493 + (-t481 * t494 - t493 * t519) * qJD(5), -t468 * t493 - t477 * t520, 0, -t493 * t515 - t469 * t494 + (t478 * t493 - t486 * t494) * qJD(5), 0; 0, t476 * qJD(3) - t484 * t528 + t485 * t502, t465, 0, 0, 0; 0, t475 * qJD(3) + t482 * t528 - t483 * t502, t463, 0, 0, 0; 0, (t509 * qJD(2) + t506 * qJD(3)) * t498, t469, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:26
	% EndTime: 2019-10-09 22:35:27
	% DurationCPUTime: 1.07s
	% Computational Cost: add. (940->148), mult. (2555->280), div. (0->0), fcn. (2925->14), ass. (0->111)
	t697 = cos(pkin(7));
	t703 = cos(qJ(3));
	t704 = cos(qJ(2));
	t733 = t703 * t704;
	t700 = sin(qJ(3));
	t701 = sin(qJ(2));
	t736 = t700 * t701;
	t711 = t697 * t733 - t736;
	t693 = sin(pkin(12));
	t696 = cos(pkin(12));
	t698 = cos(pkin(6));
	t737 = t698 * t704;
	t680 = -t693 * t701 + t696 * t737;
	t738 = t698 * t701;
	t681 = t693 * t704 + t696 * t738;
	t751 = t681 * t700;
	t750 = t681 * t703;
	t712 = t693 * t738 - t696 * t704;
	t749 = t712 * t700;
	t692 = pkin(13) + qJ(5);
	t690 = sin(t692);
	t694 = sin(pkin(7));
	t748 = t690 * t694;
	t691 = cos(t692);
	t747 = t691 * t694;
	t695 = sin(pkin(6));
	t745 = t694 * t695;
	t744 = t694 * t700;
	t743 = t694 * t703;
	t742 = t695 * t696;
	t741 = t695 * t697;
	t740 = t697 * t700;
	t739 = t697 * t703;
	t735 = t700 * t704;
	t734 = t701 * t703;
	t732 = qJD(5) * t690;
	t731 = qJD(5) * t691;
	t730 = qJD(6) * t691;
	t699 = sin(qJ(6));
	t729 = qJD(6) * t699;
	t702 = cos(qJ(6));
	t728 = qJD(6) * t702;
	t727 = t701 * t745;
	t726 = t695 * t743;
	t723 = t698 * t743;
	t722 = qJD(2) * t745;
	t721 = qJD(3) * t744;
	t720 = t701 * t722;
	t719 = t704 * t722;
	t675 = t680 * qJD(2);
	t676 = t681 * qJD(2);
	t715 = t680 * t697 - t694 * t742;
	t637 = -t676 * t740 + t675 * t703 + (t715 * t703 - t751) * qJD(3);
	t654 = -t680 * t739 + t696 * t726 + t751;
	t718 = t654 * t730 + t637;
	t713 = t693 * t737 + t696 * t701;
	t677 = t713 * qJD(2);
	t678 = t712 * qJD(2);
	t714 = t693 * t745 - t697 * t713;
	t639 = t678 * t740 - t677 * t703 + (t714 * t703 + t749) * qJD(3);
	t656 = -t693 * t726 + t713 * t739 - t749;
	t717 = t656 * t730 + t639;
	t708 = -t697 * t736 + t733;
	t653 = qJD(3) * t723 + (t708 * qJD(2) + t711 * qJD(3)) * t695;
	t665 = -t711 * t695 - t723;
	t716 = t665 * t730 + t653;
	t655 = t715 * t700 + t750;
	t667 = -t680 * t694 - t696 * t741;
	t641 = t655 * t691 + t667 * t690;
	t640 = -t655 * t690 + t667 * t691;
	t657 = t714 * t700 - t703 * t712;
	t668 = t693 * t741 + t694 * t713;
	t643 = t657 * t691 + t668 * t690;
	t642 = -t657 * t690 + t668 * t691;
	t709 = t697 * t735 + t734;
	t666 = t709 * t695 + t698 * t744;
	t679 = t698 * t697 - t704 * t745;
	t651 = t666 * t691 + t679 * t690;
	t650 = -t666 * t690 + t679 * t691;
	t661 = t680 * t703 - t681 * t740;
	t648 = t661 * t691 + t681 * t748;
	t663 = -t703 * t713 + t712 * t740;
	t649 = t663 * t691 - t712 * t748;
	t660 = t680 * t700 + t681 * t739;
	t662 = -t700 * t713 - t712 * t739;
	t710 = t697 * t734 + t735;
	t674 = t708 * t695;
	t664 = t674 * t691 + t690 * t727;
	t636 = t675 * t700 + t676 * t739 - t721 * t742 + (t680 * t740 + t750) * qJD(3);
	t707 = qJD(6) * t655 - t636 * t691 + t654 * t732;
	t638 = t657 * qJD(3) - t677 * t700 - t678 * t739;
	t706 = qJD(6) * t657 - t638 * t691 + t656 * t732;
	t652 = t698 * t721 + (t710 * qJD(2) + t709 * qJD(3)) * t695;
	t705 = qJD(6) * t666 - t652 * t691 + t665 * t732;
	t673 = t710 * t695;
	t659 = (-t709 * qJD(2) - t710 * qJD(3)) * t695;
	t658 = (t711 * qJD(2) + t708 * qJD(3)) * t695;
	t647 = -t662 * qJD(3) + t677 * t740 + t678 * t703;
	t646 = t663 * qJD(3) - t677 * t739 + t678 * t700;
	t645 = -t660 * qJD(3) - t675 * t740 - t676 * t703;
	t644 = t661 * qJD(3) + t675 * t739 - t676 * t700;
	t635 = t690 * t719 + t659 * t691 + (-t674 * t690 + t691 * t727) * qJD(5);
	t634 = t650 * qJD(5) + t653 * t691 + t690 * t720;
	t633 = -t651 * qJD(5) - t653 * t690 + t691 * t720;
	t632 = -t677 * t748 + t647 * t691 + (-t663 * t690 - t712 * t747) * qJD(5);
	t631 = t675 * t748 + t645 * t691 + (-t661 * t690 + t681 * t747) * qJD(5);
	t630 = t642 * qJD(5) + t639 * t691 - t678 * t748;
	t629 = -t643 * qJD(5) - t639 * t690 - t678 * t747;
	t628 = t640 * qJD(5) + t637 * t691 + t676 * t748;
	t627 = -t641 * qJD(5) - t637 * t690 + t676 * t747;
	t1 = [0, t632 * t702 + t646 * t699 + (-t649 * t699 + t662 * t702) * qJD(6), t717 * t699 + t706 * t702, 0, t629 * t702 - t642 * t729, -t630 * t699 + t638 * t702 + (-t643 * t702 - t656 * t699) * qJD(6); 0, t631 * t702 + t644 * t699 + (-t648 * t699 + t660 * t702) * qJD(6), t718 * t699 + t707 * t702, 0, t627 * t702 - t640 * t729, -t628 * t699 + t636 * t702 + (-t641 * t702 - t654 * t699) * qJD(6); 0, t635 * t702 + t658 * t699 + (-t664 * t699 + t673 * t702) * qJD(6), t716 * t699 + t705 * t702, 0, t633 * t702 - t650 * t729, -t634 * t699 + t652 * t702 + (-t651 * t702 - t665 * t699) * qJD(6); 0, -t632 * t699 + t646 * t702 + (-t649 * t702 - t662 * t699) * qJD(6), -t706 * t699 + t717 * t702, 0, -t629 * t699 - t642 * t728, -t630 * t702 - t638 * t699 + (t643 * t699 - t656 * t702) * qJD(6); 0, -t631 * t699 + t644 * t702 + (-t648 * t702 - t660 * t699) * qJD(6), -t707 * t699 + t718 * t702, 0, -t627 * t699 - t640 * t728, -t628 * t702 - t636 * t699 + (t641 * t699 - t654 * t702) * qJD(6); 0, -t635 * t699 + t658 * t702 + (-t664 * t702 - t673 * t699) * qJD(6), -t705 * t699 + t716 * t702, 0, -t633 * t699 - t650 * t728, -t634 * t702 - t652 * t699 + (t651 * t699 - t665 * t702) * qJD(6); 0, t649 * qJD(5) + t647 * t690 + t677 * t747, -t638 * t690 - t656 * t731, 0, t630, 0; 0, t648 * qJD(5) + t645 * t690 - t675 * t747, -t636 * t690 - t654 * t731, 0, t628, 0; 0, t664 * qJD(5) + t659 * t690 - t691 * t719, -t652 * t690 - t665 * t731, 0, t634, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end