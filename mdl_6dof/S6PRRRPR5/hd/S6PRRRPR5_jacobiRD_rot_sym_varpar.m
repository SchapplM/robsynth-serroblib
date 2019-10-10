% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRRPR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:08
	% EndTime: 2019-10-09 22:54:08
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
	% StartTime: 2019-10-09 22:54:09
	% EndTime: 2019-10-09 22:54:09
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
	% StartTime: 2019-10-09 22:54:11
	% EndTime: 2019-10-09 22:54:11
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (269->87), mult. (921->177), div. (0->0), fcn. (1021->12), ass. (0->68)
	t472 = sin(pkin(12));
	t475 = cos(pkin(12));
	t480 = sin(qJ(2));
	t477 = cos(pkin(6));
	t483 = cos(qJ(2));
	t504 = t477 * t483;
	t466 = -t472 * t480 + t475 * t504;
	t479 = sin(qJ(3));
	t482 = cos(qJ(3));
	t505 = t477 * t480;
	t488 = t472 * t505 - t475 * t483;
	t476 = cos(pkin(7));
	t489 = t472 * t504 + t475 * t480;
	t473 = sin(pkin(7));
	t474 = sin(pkin(6));
	t512 = t473 * t474;
	t490 = t472 * t512 - t476 * t489;
	t452 = t490 * t479 - t482 * t488;
	t467 = t472 * t483 + t475 * t505;
	t491 = -t466 * t476 + t475 * t512;
	t516 = -t467 * t482 + t491 * t479;
	t511 = t473 * t477;
	t478 = sin(qJ(4));
	t510 = t473 * t478;
	t481 = cos(qJ(4));
	t509 = t473 * t481;
	t508 = t474 * t476;
	t507 = t476 * t479;
	t506 = t476 * t482;
	t503 = t479 * t480;
	t502 = t479 * t483;
	t501 = t480 * t482;
	t500 = t482 * t483;
	t499 = qJD(4) * t478;
	t498 = qJD(4) * t481;
	t497 = t480 * t512;
	t495 = qJD(2) * t512;
	t494 = qJD(3) * t511;
	t493 = t480 * t495;
	t492 = t483 * t495;
	t454 = t466 * t482 - t467 * t507;
	t455 = -t482 * t489 + t488 * t507;
	t487 = t476 * t500 - t503;
	t486 = -t476 * t501 - t502;
	t485 = t476 * t502 + t501;
	t484 = -t476 * t503 + t500;
	t449 = -t467 * t479 - t491 * t482;
	t451 = t479 * t488 + t490 * t482;
	t465 = t477 * t476 - t483 * t512;
	t464 = t488 * qJD(2);
	t463 = t489 * qJD(2);
	t462 = t467 * qJD(2);
	t461 = t466 * qJD(2);
	t460 = t484 * t474;
	t459 = t472 * t508 + t473 * t489;
	t458 = -t466 * t473 - t475 * t508;
	t457 = t485 * t474 + t479 * t511;
	t456 = t487 * t474 + t482 * t511;
	t453 = (-t485 * qJD(2) + t486 * qJD(3)) * t474;
	t448 = t482 * t494 + (t484 * qJD(2) + t487 * qJD(3)) * t474;
	t447 = -t479 * t494 + (t486 * qJD(2) - t485 * qJD(3)) * t474;
	t446 = t463 * t507 + t464 * t482 + (t479 * t489 + t488 * t506) * qJD(3);
	t445 = -t461 * t507 - t462 * t482 + (-t466 * t479 - t467 * t506) * qJD(3);
	t444 = t451 * qJD(3) - t463 * t482 + t464 * t507;
	t443 = -t452 * qJD(3) + t463 * t479 + t464 * t506;
	t442 = t449 * qJD(3) + t461 * t482 - t462 * t507;
	t441 = t516 * qJD(3) - t461 * t479 - t462 * t506;
	t1 = [0, -t463 * t510 + t446 * t481 + (-t455 * t478 - t488 * t509) * qJD(4), t443 * t481 - t451 * t499, -t464 * t509 - t444 * t478 + (-t452 * t481 - t459 * t478) * qJD(4), 0, 0; 0, t461 * t510 + t445 * t481 + (-t454 * t478 + t467 * t509) * qJD(4), t441 * t481 - t449 * t499, t462 * t509 - t442 * t478 + (-t458 * t478 + t481 * t516) * qJD(4), 0, 0; 0, t478 * t492 + t453 * t481 + (-t460 * t478 + t481 * t497) * qJD(4), t447 * t481 - t456 * t499, t481 * t493 - t448 * t478 + (-t457 * t481 - t465 * t478) * qJD(4), 0, 0; 0, -t463 * t509 - t446 * t478 + (-t455 * t481 + t488 * t510) * qJD(4), -t443 * t478 - t451 * t498, t464 * t510 - t444 * t481 + (t452 * t478 - t459 * t481) * qJD(4), 0, 0; 0, t461 * t509 - t445 * t478 + (-t454 * t481 - t467 * t510) * qJD(4), -t441 * t478 - t449 * t498, -t462 * t510 - t442 * t481 + (-t458 * t481 - t478 * t516) * qJD(4), 0, 0; 0, t481 * t492 - t453 * t478 + (-t460 * t481 - t478 * t497) * qJD(4), -t447 * t478 - t456 * t498, -t478 * t493 - t448 * t481 + (t457 * t478 - t465 * t481) * qJD(4), 0, 0; 0, t455 * qJD(3) - t463 * t506 + t464 * t479, t444, 0, 0, 0; 0, t454 * qJD(3) + t461 * t506 - t462 * t479, t442, 0, 0, 0; 0, (t487 * qJD(2) + t484 * qJD(3)) * t474, t448, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:11
	% EndTime: 2019-10-09 22:54:11
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (329->88), mult. (921->177), div. (0->0), fcn. (1021->12), ass. (0->69)
	t502 = sin(pkin(12));
	t505 = cos(pkin(12));
	t509 = sin(qJ(2));
	t507 = cos(pkin(6));
	t511 = cos(qJ(2));
	t532 = t507 * t511;
	t493 = -t502 * t509 + t505 * t532;
	t508 = sin(qJ(3));
	t510 = cos(qJ(3));
	t533 = t507 * t509;
	t516 = t502 * t533 - t505 * t511;
	t506 = cos(pkin(7));
	t517 = t502 * t532 + t505 * t509;
	t503 = sin(pkin(7));
	t504 = sin(pkin(6));
	t538 = t503 * t504;
	t518 = t502 * t538 - t506 * t517;
	t479 = t518 * t508 - t510 * t516;
	t494 = t502 * t511 + t505 * t533;
	t519 = -t493 * t506 + t505 * t538;
	t544 = -t494 * t510 + t519 * t508;
	t501 = qJ(4) + pkin(13);
	t499 = sin(t501);
	t541 = t499 * t503;
	t500 = cos(t501);
	t540 = t500 * t503;
	t537 = t503 * t507;
	t536 = t504 * t506;
	t535 = t506 * t508;
	t534 = t506 * t510;
	t531 = t508 * t509;
	t530 = t508 * t511;
	t529 = t509 * t510;
	t528 = t510 * t511;
	t527 = qJD(4) * t499;
	t526 = qJD(4) * t500;
	t525 = t509 * t538;
	t523 = qJD(2) * t538;
	t522 = qJD(3) * t537;
	t521 = t509 * t523;
	t520 = t511 * t523;
	t481 = t493 * t510 - t494 * t535;
	t482 = -t510 * t517 + t516 * t535;
	t515 = t506 * t528 - t531;
	t514 = -t506 * t529 - t530;
	t513 = t506 * t530 + t529;
	t512 = -t506 * t531 + t528;
	t476 = -t494 * t508 - t519 * t510;
	t478 = t508 * t516 + t518 * t510;
	t492 = t507 * t506 - t511 * t538;
	t491 = t516 * qJD(2);
	t490 = t517 * qJD(2);
	t489 = t494 * qJD(2);
	t488 = t493 * qJD(2);
	t487 = t512 * t504;
	t486 = t502 * t536 + t503 * t517;
	t485 = -t493 * t503 - t505 * t536;
	t484 = t513 * t504 + t508 * t537;
	t483 = t515 * t504 + t510 * t537;
	t480 = (-t513 * qJD(2) + t514 * qJD(3)) * t504;
	t475 = t510 * t522 + (t512 * qJD(2) + t515 * qJD(3)) * t504;
	t474 = -t508 * t522 + (t514 * qJD(2) - t513 * qJD(3)) * t504;
	t473 = t490 * t535 + t491 * t510 + (t508 * t517 + t516 * t534) * qJD(3);
	t472 = -t488 * t535 - t489 * t510 + (-t493 * t508 - t494 * t534) * qJD(3);
	t471 = t478 * qJD(3) - t490 * t510 + t491 * t535;
	t470 = -t479 * qJD(3) + t490 * t508 + t491 * t534;
	t469 = t476 * qJD(3) + t488 * t510 - t489 * t535;
	t468 = t544 * qJD(3) - t488 * t508 - t489 * t534;
	t1 = [0, -t490 * t541 + t473 * t500 + (-t482 * t499 - t516 * t540) * qJD(4), t470 * t500 - t478 * t527, -t491 * t540 - t471 * t499 + (-t479 * t500 - t486 * t499) * qJD(4), 0, 0; 0, t488 * t541 + t472 * t500 + (-t481 * t499 + t494 * t540) * qJD(4), t468 * t500 - t476 * t527, t489 * t540 - t469 * t499 + (-t485 * t499 + t500 * t544) * qJD(4), 0, 0; 0, t499 * t520 + t480 * t500 + (-t487 * t499 + t500 * t525) * qJD(4), t474 * t500 - t483 * t527, t500 * t521 - t475 * t499 + (-t484 * t500 - t492 * t499) * qJD(4), 0, 0; 0, -t490 * t540 - t473 * t499 + (-t482 * t500 + t516 * t541) * qJD(4), -t470 * t499 - t478 * t526, t491 * t541 - t471 * t500 + (t479 * t499 - t486 * t500) * qJD(4), 0, 0; 0, t488 * t540 - t472 * t499 + (-t481 * t500 - t494 * t541) * qJD(4), -t468 * t499 - t476 * t526, -t489 * t541 - t469 * t500 + (-t485 * t500 - t499 * t544) * qJD(4), 0, 0; 0, t500 * t520 - t480 * t499 + (-t487 * t500 - t499 * t525) * qJD(4), -t474 * t499 - t483 * t526, -t499 * t521 - t475 * t500 + (t484 * t499 - t492 * t500) * qJD(4), 0, 0; 0, t482 * qJD(3) - t490 * t534 + t491 * t508, t471, 0, 0, 0; 0, t481 * qJD(3) + t488 * t534 - t489 * t508, t469, 0, 0, 0; 0, (t515 * qJD(2) + t512 * qJD(3)) * t504, t475, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:15
	% EndTime: 2019-10-09 22:54:16
	% DurationCPUTime: 1.07s
	% Computational Cost: add. (940->148), mult. (2555->280), div. (0->0), fcn. (2925->14), ass. (0->111)
	t702 = cos(pkin(7));
	t708 = cos(qJ(3));
	t709 = cos(qJ(2));
	t738 = t708 * t709;
	t705 = sin(qJ(3));
	t706 = sin(qJ(2));
	t741 = t705 * t706;
	t716 = t702 * t738 - t741;
	t698 = sin(pkin(12));
	t701 = cos(pkin(12));
	t703 = cos(pkin(6));
	t742 = t703 * t709;
	t685 = -t698 * t706 + t701 * t742;
	t743 = t703 * t706;
	t686 = t698 * t709 + t701 * t743;
	t756 = t686 * t705;
	t755 = t686 * t708;
	t717 = t698 * t743 - t701 * t709;
	t754 = t717 * t705;
	t697 = qJ(4) + pkin(13);
	t695 = sin(t697);
	t699 = sin(pkin(7));
	t753 = t695 * t699;
	t696 = cos(t697);
	t752 = t696 * t699;
	t700 = sin(pkin(6));
	t750 = t699 * t700;
	t749 = t699 * t705;
	t748 = t699 * t708;
	t747 = t700 * t701;
	t746 = t700 * t702;
	t745 = t702 * t705;
	t744 = t702 * t708;
	t740 = t705 * t709;
	t739 = t706 * t708;
	t737 = qJD(4) * t695;
	t736 = qJD(4) * t696;
	t735 = qJD(6) * t696;
	t704 = sin(qJ(6));
	t734 = qJD(6) * t704;
	t707 = cos(qJ(6));
	t733 = qJD(6) * t707;
	t732 = t706 * t750;
	t731 = t700 * t748;
	t728 = t703 * t748;
	t727 = qJD(2) * t750;
	t726 = qJD(3) * t749;
	t725 = t706 * t727;
	t724 = t709 * t727;
	t680 = t685 * qJD(2);
	t681 = t686 * qJD(2);
	t720 = t685 * t702 - t699 * t747;
	t642 = -t681 * t745 + t680 * t708 + (t720 * t708 - t756) * qJD(3);
	t659 = -t685 * t744 + t701 * t731 + t756;
	t723 = t659 * t735 + t642;
	t718 = t698 * t742 + t701 * t706;
	t682 = t718 * qJD(2);
	t683 = t717 * qJD(2);
	t719 = t698 * t750 - t702 * t718;
	t644 = t683 * t745 - t682 * t708 + (t719 * t708 + t754) * qJD(3);
	t661 = -t698 * t731 + t718 * t744 - t754;
	t722 = t661 * t735 + t644;
	t713 = -t702 * t741 + t738;
	t658 = qJD(3) * t728 + (t713 * qJD(2) + t716 * qJD(3)) * t700;
	t670 = -t716 * t700 - t728;
	t721 = t670 * t735 + t658;
	t660 = t720 * t705 + t755;
	t672 = -t685 * t699 - t701 * t746;
	t646 = t660 * t696 + t672 * t695;
	t645 = -t660 * t695 + t672 * t696;
	t662 = t719 * t705 - t708 * t717;
	t673 = t698 * t746 + t699 * t718;
	t648 = t662 * t696 + t673 * t695;
	t647 = -t662 * t695 + t673 * t696;
	t714 = t702 * t740 + t739;
	t671 = t714 * t700 + t703 * t749;
	t684 = t703 * t702 - t709 * t750;
	t656 = t671 * t696 + t684 * t695;
	t655 = -t671 * t695 + t684 * t696;
	t666 = t685 * t708 - t686 * t745;
	t653 = t666 * t696 + t686 * t753;
	t668 = -t708 * t718 + t717 * t745;
	t654 = t668 * t696 - t717 * t753;
	t665 = t685 * t705 + t686 * t744;
	t667 = -t705 * t718 - t717 * t744;
	t715 = t702 * t739 + t740;
	t679 = t713 * t700;
	t669 = t679 * t696 + t695 * t732;
	t641 = t680 * t705 + t681 * t744 - t726 * t747 + (t685 * t745 + t755) * qJD(3);
	t712 = qJD(6) * t660 - t641 * t696 + t659 * t737;
	t643 = t662 * qJD(3) - t682 * t705 - t683 * t744;
	t711 = qJD(6) * t662 - t643 * t696 + t661 * t737;
	t657 = t703 * t726 + (t715 * qJD(2) + t714 * qJD(3)) * t700;
	t710 = qJD(6) * t671 - t657 * t696 + t670 * t737;
	t678 = t715 * t700;
	t664 = (-t714 * qJD(2) - t715 * qJD(3)) * t700;
	t663 = (t716 * qJD(2) + t713 * qJD(3)) * t700;
	t652 = -t667 * qJD(3) + t682 * t745 + t683 * t708;
	t651 = t668 * qJD(3) - t682 * t744 + t683 * t705;
	t650 = -t665 * qJD(3) - t680 * t745 - t681 * t708;
	t649 = t666 * qJD(3) + t680 * t744 - t681 * t705;
	t640 = t695 * t724 + t664 * t696 + (-t679 * t695 + t696 * t732) * qJD(4);
	t639 = t655 * qJD(4) + t658 * t696 + t695 * t725;
	t638 = -t656 * qJD(4) - t658 * t695 + t696 * t725;
	t637 = -t682 * t753 + t652 * t696 + (-t668 * t695 - t717 * t752) * qJD(4);
	t636 = t680 * t753 + t650 * t696 + (-t666 * t695 + t686 * t752) * qJD(4);
	t635 = t647 * qJD(4) + t644 * t696 - t683 * t753;
	t634 = -t648 * qJD(4) - t644 * t695 - t683 * t752;
	t633 = t645 * qJD(4) + t642 * t696 + t681 * t753;
	t632 = -t646 * qJD(4) - t642 * t695 + t681 * t752;
	t1 = [0, t637 * t707 + t651 * t704 + (-t654 * t704 + t667 * t707) * qJD(6), t722 * t704 + t711 * t707, t634 * t707 - t647 * t734, 0, -t635 * t704 + t643 * t707 + (-t648 * t707 - t661 * t704) * qJD(6); 0, t636 * t707 + t649 * t704 + (-t653 * t704 + t665 * t707) * qJD(6), t723 * t704 + t712 * t707, t632 * t707 - t645 * t734, 0, -t633 * t704 + t641 * t707 + (-t646 * t707 - t659 * t704) * qJD(6); 0, t640 * t707 + t663 * t704 + (-t669 * t704 + t678 * t707) * qJD(6), t721 * t704 + t710 * t707, t638 * t707 - t655 * t734, 0, -t639 * t704 + t657 * t707 + (-t656 * t707 - t670 * t704) * qJD(6); 0, -t637 * t704 + t651 * t707 + (-t654 * t707 - t667 * t704) * qJD(6), -t711 * t704 + t722 * t707, -t634 * t704 - t647 * t733, 0, -t635 * t707 - t643 * t704 + (t648 * t704 - t661 * t707) * qJD(6); 0, -t636 * t704 + t649 * t707 + (-t653 * t707 - t665 * t704) * qJD(6), -t712 * t704 + t723 * t707, -t632 * t704 - t645 * t733, 0, -t633 * t707 - t641 * t704 + (t646 * t704 - t659 * t707) * qJD(6); 0, -t640 * t704 + t663 * t707 + (-t669 * t707 - t678 * t704) * qJD(6), -t710 * t704 + t721 * t707, -t638 * t704 - t655 * t733, 0, -t639 * t707 - t657 * t704 + (t656 * t704 - t670 * t707) * qJD(6); 0, t654 * qJD(4) + t652 * t695 + t682 * t752, -t643 * t695 - t661 * t736, t635, 0, 0; 0, t653 * qJD(4) + t650 * t695 - t680 * t752, -t641 * t695 - t659 * t736, t633, 0, 0; 0, t669 * qJD(4) + t664 * t695 - t696 * t724, -t657 * t695 - t670 * t736, t639, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end