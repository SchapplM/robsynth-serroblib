% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR8
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRPRR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
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
	% StartTime: 2019-10-09 22:39:11
	% EndTime: 2019-10-09 22:39:11
	% DurationCPUTime: 0.19s
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
	% StartTime: 2019-10-09 22:39:12
	% EndTime: 2019-10-09 22:39:12
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (78->41), mult. (285->88), div. (0->0), fcn. (301->10), ass. (0->35)
	t359 = sin(pkin(12));
	t362 = cos(pkin(12));
	t366 = sin(qJ(2));
	t364 = cos(pkin(6));
	t368 = cos(qJ(2));
	t383 = t364 * t368;
	t353 = -t359 * t366 + t362 * t383;
	t360 = sin(pkin(7));
	t361 = sin(pkin(6));
	t387 = t360 * t361;
	t363 = cos(pkin(7));
	t365 = sin(qJ(3));
	t386 = t363 * t365;
	t367 = cos(qJ(3));
	t385 = t363 * t367;
	t384 = t364 * t366;
	t382 = t365 * t366;
	t381 = t365 * t368;
	t380 = t366 * t367;
	t379 = t367 * t368;
	t377 = qJD(3) * t360 * t364;
	t376 = t353 * t363 - t362 * t387;
	t374 = t359 * t383 + t362 * t366;
	t375 = t359 * t387 - t363 * t374;
	t354 = t359 * t368 + t362 * t384;
	t373 = t359 * t384 - t362 * t368;
	t372 = t363 * t379 - t382;
	t371 = t363 * t380 + t381;
	t370 = t363 * t381 + t380;
	t369 = -t363 * t382 + t379;
	t352 = t373 * qJD(2);
	t351 = t374 * qJD(2);
	t350 = t354 * qJD(2);
	t349 = t353 * qJD(2);
	t1 = [0, -t351 * t360, 0, 0, 0, 0; 0, t349 * t360, 0, 0, 0, 0; 0, qJD(2) * t368 * t387, 0, 0, 0, 0; 0, -t351 * t386 - t352 * t367 + (-t365 * t374 - t373 * t385) * qJD(3), -t352 * t385 - t351 * t365 + (t375 * t365 - t367 * t373) * qJD(3), 0, 0, 0; 0, t349 * t386 + t350 * t367 + (t353 * t365 + t354 * t385) * qJD(3), t350 * t385 + t349 * t365 + (t354 * t367 + t376 * t365) * qJD(3), 0, 0, 0; 0, (t370 * qJD(2) + t371 * qJD(3)) * t361, t365 * t377 + (t371 * qJD(2) + t370 * qJD(3)) * t361, 0, 0, 0; 0, -t351 * t385 + t352 * t365 + (-t367 * t374 + t373 * t386) * qJD(3), t352 * t386 - t351 * t367 + (t365 * t373 + t375 * t367) * qJD(3), 0, 0, 0; 0, t349 * t385 - t350 * t365 + (t353 * t367 - t354 * t386) * qJD(3), -t350 * t386 + t349 * t367 + (-t354 * t365 + t376 * t367) * qJD(3), 0, 0, 0; 0, (t372 * qJD(2) + t369 * qJD(3)) * t361, t367 * t377 + (t369 * qJD(2) + t372 * qJD(3)) * t361, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:13
	% EndTime: 2019-10-09 22:39:13
	% DurationCPUTime: 0.60s
	% Computational Cost: add. (269->90), mult. (921->177), div. (0->0), fcn. (1021->12), ass. (0->68)
	t480 = sin(pkin(12));
	t483 = cos(pkin(12));
	t488 = sin(qJ(2));
	t485 = cos(pkin(6));
	t491 = cos(qJ(2));
	t512 = t485 * t491;
	t474 = -t480 * t488 + t483 * t512;
	t487 = sin(qJ(3));
	t490 = cos(qJ(3));
	t513 = t485 * t488;
	t496 = t480 * t513 - t483 * t491;
	t484 = cos(pkin(7));
	t497 = t480 * t512 + t483 * t488;
	t481 = sin(pkin(7));
	t482 = sin(pkin(6));
	t520 = t481 * t482;
	t498 = t480 * t520 - t484 * t497;
	t524 = t487 * t496 + t498 * t490;
	t475 = t480 * t491 + t483 * t513;
	t499 = -t474 * t484 + t483 * t520;
	t457 = t475 * t487 + t499 * t490;
	t519 = t481 * t485;
	t486 = sin(qJ(5));
	t518 = t481 * t486;
	t489 = cos(qJ(5));
	t517 = t481 * t489;
	t516 = t482 * t484;
	t515 = t484 * t487;
	t514 = t484 * t490;
	t511 = t487 * t488;
	t510 = t487 * t491;
	t509 = t488 * t490;
	t508 = t490 * t491;
	t507 = qJD(5) * t486;
	t506 = qJD(5) * t489;
	t505 = t488 * t520;
	t503 = qJD(2) * t520;
	t502 = qJD(3) * t519;
	t501 = t488 * t503;
	t500 = t491 * t503;
	t462 = t474 * t487 + t475 * t514;
	t463 = -t487 * t497 - t496 * t514;
	t495 = t484 * t508 - t511;
	t494 = t484 * t509 + t510;
	t493 = t484 * t510 + t509;
	t492 = -t484 * t511 + t508;
	t458 = t475 * t490 - t499 * t487;
	t460 = t498 * t487 - t490 * t496;
	t473 = t485 * t484 - t491 * t520;
	t472 = t496 * qJD(2);
	t471 = t497 * qJD(2);
	t470 = t475 * qJD(2);
	t469 = t474 * qJD(2);
	t468 = t494 * t482;
	t467 = t480 * t516 + t481 * t497;
	t466 = -t474 * t481 - t483 * t516;
	t465 = t493 * t482 + t487 * t519;
	t464 = -t495 * t482 - t490 * t519;
	t461 = (t495 * qJD(2) + t492 * qJD(3)) * t482;
	t456 = t490 * t502 + (t492 * qJD(2) + t495 * qJD(3)) * t482;
	t455 = t487 * t502 + (t494 * qJD(2) + t493 * qJD(3)) * t482;
	t454 = -t471 * t514 + t472 * t487 + (-t490 * t497 + t496 * t515) * qJD(3);
	t453 = t469 * t514 - t470 * t487 + (t474 * t490 - t475 * t515) * qJD(3);
	t452 = t524 * qJD(3) - t471 * t490 + t472 * t515;
	t451 = t460 * qJD(3) - t471 * t487 - t472 * t514;
	t450 = -t457 * qJD(3) + t469 * t490 - t470 * t515;
	t449 = t458 * qJD(3) + t469 * t487 + t470 * t514;
	t1 = [0, -t471 * t517 + t454 * t486 + (t463 * t489 + t496 * t518) * qJD(5), t452 * t486 + t460 * t506, 0, t472 * t518 + t451 * t489 + (-t467 * t489 + t486 * t524) * qJD(5), 0; 0, t469 * t517 + t453 * t486 + (t462 * t489 - t475 * t518) * qJD(5), t450 * t486 + t458 * t506, 0, -t470 * t518 + t449 * t489 + (-t457 * t486 - t466 * t489) * qJD(5), 0; 0, t489 * t500 + t461 * t486 + (t468 * t489 - t486 * t505) * qJD(5), t456 * t486 + t465 * t506, 0, -t486 * t501 + t455 * t489 + (-t464 * t486 - t473 * t489) * qJD(5), 0; 0, t471 * t518 + t454 * t489 + (-t463 * t486 + t496 * t517) * qJD(5), t452 * t489 - t460 * t507, 0, t472 * t517 - t451 * t486 + (t467 * t486 + t489 * t524) * qJD(5), 0; 0, -t469 * t518 + t453 * t489 + (-t462 * t486 - t475 * t517) * qJD(5), t450 * t489 - t458 * t507, 0, -t470 * t517 - t449 * t486 + (-t457 * t489 + t466 * t486) * qJD(5), 0; 0, -t486 * t500 + t461 * t489 + (-t468 * t486 - t489 * t505) * qJD(5), t456 * t489 - t465 * t507, 0, -t489 * t501 - t455 * t486 + (-t464 * t489 + t473 * t486) * qJD(5), 0; 0, -t463 * qJD(3) + t471 * t515 + t472 * t490, -t451, 0, 0, 0; 0, -t462 * qJD(3) - t469 * t515 - t470 * t490, -t449, 0, 0, 0; 0, (-t493 * qJD(2) - t494 * qJD(3)) * t482, -t455, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:16
	% EndTime: 2019-10-09 22:39:17
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (784->147), mult. (2555->280), div. (0->0), fcn. (2925->14), ass. (0->110)
	t702 = cos(pkin(7));
	t710 = cos(qJ(3));
	t711 = cos(qJ(2));
	t740 = t710 * t711;
	t706 = sin(qJ(3));
	t707 = sin(qJ(2));
	t743 = t706 * t707;
	t718 = t702 * t740 - t743;
	t698 = sin(pkin(12));
	t701 = cos(pkin(12));
	t703 = cos(pkin(6));
	t744 = t703 * t711;
	t688 = -t698 * t707 + t701 * t744;
	t745 = t703 * t707;
	t689 = t698 * t711 + t701 * t745;
	t758 = t689 * t706;
	t757 = t689 * t710;
	t719 = t698 * t745 - t701 * t711;
	t756 = t719 * t706;
	t699 = sin(pkin(7));
	t700 = sin(pkin(6));
	t754 = t699 * t700;
	t705 = sin(qJ(5));
	t753 = t699 * t705;
	t752 = t699 * t706;
	t709 = cos(qJ(5));
	t751 = t699 * t709;
	t750 = t699 * t710;
	t749 = t700 * t701;
	t748 = t700 * t702;
	t747 = t702 * t706;
	t746 = t702 * t710;
	t742 = t706 * t711;
	t741 = t707 * t710;
	t739 = qJD(5) * t705;
	t738 = qJD(5) * t709;
	t704 = sin(qJ(6));
	t737 = qJD(6) * t704;
	t736 = qJD(6) * t705;
	t708 = cos(qJ(6));
	t735 = qJD(6) * t708;
	t734 = t707 * t754;
	t733 = t700 * t750;
	t730 = t703 * t750;
	t729 = qJD(2) * t754;
	t728 = qJD(3) * t752;
	t727 = t707 * t729;
	t726 = t711 * t729;
	t683 = t688 * qJD(2);
	t684 = t689 * qJD(2);
	t643 = t683 * t706 + t684 * t746 - t728 * t749 + (t688 * t747 + t757) * qJD(3);
	t722 = t688 * t702 - t699 * t749;
	t661 = t722 * t706 + t757;
	t725 = -t661 * t736 - t643;
	t720 = t698 * t744 + t701 * t707;
	t721 = t698 * t754 - t702 * t720;
	t663 = t721 * t706 - t710 * t719;
	t685 = t720 * qJD(2);
	t686 = t719 * qJD(2);
	t645 = t663 * qJD(3) - t685 * t706 - t686 * t746;
	t724 = -t663 * t736 - t645;
	t716 = t702 * t742 + t741;
	t717 = t702 * t741 + t742;
	t658 = t703 * t728 + (t717 * qJD(2) + t716 * qJD(3)) * t700;
	t674 = t716 * t700 + t703 * t752;
	t723 = -t674 * t736 - t658;
	t660 = -t688 * t746 + t701 * t733 + t758;
	t675 = -t688 * t699 - t701 * t748;
	t652 = t660 * t709 - t675 * t705;
	t653 = t660 * t705 + t675 * t709;
	t662 = -t698 * t733 + t720 * t746 - t756;
	t676 = t698 * t748 + t699 * t720;
	t654 = t662 * t709 - t676 * t705;
	t655 = t662 * t705 + t676 * t709;
	t673 = -t718 * t700 - t730;
	t687 = t703 * t702 - t711 * t754;
	t664 = t673 * t709 - t687 * t705;
	t665 = t673 * t705 + t687 * t709;
	t668 = t688 * t706 + t689 * t746;
	t656 = t668 * t705 + t689 * t751;
	t670 = -t706 * t720 - t719 * t746;
	t657 = t670 * t705 - t719 * t751;
	t669 = t688 * t710 - t689 * t747;
	t671 = -t710 * t720 + t719 * t747;
	t715 = -t702 * t743 + t740;
	t681 = t717 * t700;
	t672 = t681 * t705 + t709 * t734;
	t644 = -t684 * t747 + t683 * t710 + (t722 * t710 - t758) * qJD(3);
	t714 = -qJD(6) * t660 + t644 * t705 + t661 * t738;
	t646 = t686 * t747 - t685 * t710 + (t721 * t710 + t756) * qJD(3);
	t713 = -qJD(6) * t662 + t646 * t705 + t663 * t738;
	t659 = qJD(3) * t730 + (t715 * qJD(2) + t718 * qJD(3)) * t700;
	t712 = -qJD(6) * t673 + t659 * t705 + t674 * t738;
	t682 = t715 * t700;
	t667 = (-t716 * qJD(2) - t717 * qJD(3)) * t700;
	t666 = (t718 * qJD(2) + t715 * qJD(3)) * t700;
	t651 = -t670 * qJD(3) + t685 * t747 + t686 * t710;
	t650 = t671 * qJD(3) - t685 * t746 + t686 * t706;
	t649 = -t668 * qJD(3) - t683 * t747 - t684 * t710;
	t648 = t669 * qJD(3) + t683 * t746 - t684 * t706;
	t647 = t709 * t726 + t666 * t705 + (t681 * t709 - t705 * t734) * qJD(5);
	t642 = t664 * qJD(5) + t658 * t705 + t709 * t727;
	t641 = -t665 * qJD(5) + t658 * t709 - t705 * t727;
	t640 = -t685 * t751 + t650 * t705 + (t670 * t709 + t719 * t753) * qJD(5);
	t639 = t683 * t751 + t648 * t705 + (t668 * t709 - t689 * t753) * qJD(5);
	t638 = t654 * qJD(5) + t645 * t705 - t686 * t751;
	t637 = -t655 * qJD(5) + t645 * t709 + t686 * t753;
	t636 = t652 * qJD(5) + t643 * t705 + t684 * t751;
	t635 = -t653 * qJD(5) + t643 * t709 - t684 * t753;
	t1 = [0, t640 * t708 + t651 * t704 + (-t657 * t704 + t671 * t708) * qJD(6), t724 * t704 + t713 * t708, 0, t637 * t708 - t654 * t737, -t638 * t704 + t646 * t708 + (-t655 * t708 - t663 * t704) * qJD(6); 0, t639 * t708 + t649 * t704 + (-t656 * t704 + t669 * t708) * qJD(6), t725 * t704 + t714 * t708, 0, t635 * t708 - t652 * t737, -t636 * t704 + t644 * t708 + (-t653 * t708 - t661 * t704) * qJD(6); 0, t647 * t708 + t667 * t704 + (-t672 * t704 + t682 * t708) * qJD(6), t723 * t704 + t712 * t708, 0, t641 * t708 - t664 * t737, -t642 * t704 + t659 * t708 + (-t665 * t708 - t674 * t704) * qJD(6); 0, -t640 * t704 + t651 * t708 + (-t657 * t708 - t671 * t704) * qJD(6), -t713 * t704 + t724 * t708, 0, -t637 * t704 - t654 * t735, -t638 * t708 - t646 * t704 + (t655 * t704 - t663 * t708) * qJD(6); 0, -t639 * t704 + t649 * t708 + (-t656 * t708 - t669 * t704) * qJD(6), -t714 * t704 + t725 * t708, 0, -t635 * t704 - t652 * t735, -t636 * t708 - t644 * t704 + (t653 * t704 - t661 * t708) * qJD(6); 0, -t647 * t704 + t667 * t708 + (-t672 * t708 - t682 * t704) * qJD(6), -t712 * t704 + t723 * t708, 0, -t641 * t704 - t664 * t735, -t642 * t708 - t659 * t704 + (t665 * t704 - t674 * t708) * qJD(6); 0, t657 * qJD(5) - t650 * t709 - t685 * t753, -t646 * t709 + t663 * t739, 0, t638, 0; 0, t656 * qJD(5) - t648 * t709 + t683 * t753, -t644 * t709 + t661 * t739, 0, t636, 0; 0, t672 * qJD(5) - t666 * t709 + t705 * t726, -t659 * t709 + t674 * t739, 0, t642, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end