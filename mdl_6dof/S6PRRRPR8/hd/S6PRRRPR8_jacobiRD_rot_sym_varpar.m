% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRRPR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:05
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
	% StartTime: 2019-10-09 23:00:06
	% EndTime: 2019-10-09 23:00:06
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
	% StartTime: 2019-10-09 23:00:08
	% EndTime: 2019-10-09 23:00:08
	% DurationCPUTime: 0.59s
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
	% StartTime: 2019-10-09 23:00:09
	% EndTime: 2019-10-09 23:00:10
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (269->87), mult. (921->177), div. (0->0), fcn. (1021->12), ass. (0->68)
	t527 = sin(pkin(12));
	t530 = cos(pkin(12));
	t535 = sin(qJ(2));
	t532 = cos(pkin(6));
	t538 = cos(qJ(2));
	t559 = t532 * t538;
	t521 = -t527 * t535 + t530 * t559;
	t534 = sin(qJ(3));
	t537 = cos(qJ(3));
	t560 = t532 * t535;
	t543 = t527 * t560 - t530 * t538;
	t531 = cos(pkin(7));
	t544 = t527 * t559 + t530 * t535;
	t528 = sin(pkin(7));
	t529 = sin(pkin(6));
	t567 = t528 * t529;
	t545 = t527 * t567 - t531 * t544;
	t507 = t545 * t534 - t537 * t543;
	t522 = t527 * t538 + t530 * t560;
	t546 = -t521 * t531 + t530 * t567;
	t571 = -t522 * t537 + t546 * t534;
	t566 = t528 * t532;
	t533 = sin(qJ(4));
	t565 = t528 * t533;
	t536 = cos(qJ(4));
	t564 = t528 * t536;
	t563 = t529 * t531;
	t562 = t531 * t534;
	t561 = t531 * t537;
	t558 = t534 * t535;
	t557 = t534 * t538;
	t556 = t535 * t537;
	t555 = t537 * t538;
	t554 = qJD(4) * t533;
	t553 = qJD(4) * t536;
	t552 = t535 * t567;
	t550 = qJD(2) * t567;
	t549 = qJD(3) * t566;
	t548 = t535 * t550;
	t547 = t538 * t550;
	t509 = t521 * t537 - t522 * t562;
	t510 = -t537 * t544 + t543 * t562;
	t542 = t531 * t555 - t558;
	t541 = -t531 * t556 - t557;
	t540 = t531 * t557 + t556;
	t539 = -t531 * t558 + t555;
	t504 = -t522 * t534 - t546 * t537;
	t506 = t534 * t543 + t545 * t537;
	t520 = t532 * t531 - t538 * t567;
	t519 = t543 * qJD(2);
	t518 = t544 * qJD(2);
	t517 = t522 * qJD(2);
	t516 = t521 * qJD(2);
	t515 = t539 * t529;
	t514 = t527 * t563 + t528 * t544;
	t513 = -t521 * t528 - t530 * t563;
	t512 = t540 * t529 + t534 * t566;
	t511 = t542 * t529 + t537 * t566;
	t508 = (-t540 * qJD(2) + t541 * qJD(3)) * t529;
	t503 = t537 * t549 + (t539 * qJD(2) + t542 * qJD(3)) * t529;
	t502 = -t534 * t549 + (t541 * qJD(2) - t540 * qJD(3)) * t529;
	t501 = t518 * t562 + t519 * t537 + (t534 * t544 + t543 * t561) * qJD(3);
	t500 = -t516 * t562 - t517 * t537 + (-t521 * t534 - t522 * t561) * qJD(3);
	t499 = t506 * qJD(3) - t518 * t537 + t519 * t562;
	t498 = -t507 * qJD(3) + t518 * t534 + t519 * t561;
	t497 = t504 * qJD(3) + t516 * t537 - t517 * t562;
	t496 = t571 * qJD(3) - t516 * t534 - t517 * t561;
	t1 = [0, t510 * qJD(3) - t518 * t561 + t519 * t534, t499, 0, 0, 0; 0, t509 * qJD(3) + t516 * t561 - t517 * t534, t497, 0, 0, 0; 0, (t542 * qJD(2) + t539 * qJD(3)) * t529, t503, 0, 0, 0; 0, t518 * t565 - t501 * t536 + (t510 * t533 + t543 * t564) * qJD(4), -t498 * t536 + t506 * t554, t519 * t564 + t499 * t533 + (t507 * t536 + t514 * t533) * qJD(4), 0, 0; 0, -t516 * t565 - t500 * t536 + (t509 * t533 - t522 * t564) * qJD(4), -t496 * t536 + t504 * t554, -t517 * t564 + t497 * t533 + (t513 * t533 - t536 * t571) * qJD(4), 0, 0; 0, -t533 * t547 - t508 * t536 + (t515 * t533 - t536 * t552) * qJD(4), -t502 * t536 + t511 * t554, -t536 * t548 + t503 * t533 + (t512 * t536 + t520 * t533) * qJD(4), 0, 0; 0, t518 * t564 + t501 * t533 + (t510 * t536 - t543 * t565) * qJD(4), t498 * t533 + t506 * t553, -t519 * t565 + t499 * t536 + (-t507 * t533 + t514 * t536) * qJD(4), 0, 0; 0, -t516 * t564 + t500 * t533 + (t509 * t536 + t522 * t565) * qJD(4), t496 * t533 + t504 * t553, t517 * t565 + t497 * t536 + (t513 * t536 + t533 * t571) * qJD(4), 0, 0; 0, -t536 * t547 + t508 * t533 + (t515 * t536 + t533 * t552) * qJD(4), t502 * t533 + t511 * t553, t533 * t548 + t503 * t536 + (-t512 * t533 + t520 * t536) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:12
	% EndTime: 2019-10-09 23:00:13
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (784->150), mult. (2555->280), div. (0->0), fcn. (2925->14), ass. (0->110)
	t723 = cos(pkin(7));
	t731 = cos(qJ(3));
	t732 = cos(qJ(2));
	t764 = t731 * t732;
	t727 = sin(qJ(3));
	t728 = sin(qJ(2));
	t767 = t727 * t728;
	t740 = t723 * t764 - t767;
	t719 = sin(pkin(12));
	t722 = cos(pkin(12));
	t724 = cos(pkin(6));
	t768 = t724 * t732;
	t709 = -t719 * t728 + t722 * t768;
	t769 = t724 * t728;
	t710 = t719 * t732 + t722 * t769;
	t782 = t710 * t727;
	t781 = t710 * t731;
	t741 = t719 * t769 - t722 * t732;
	t780 = t741 * t727;
	t720 = sin(pkin(7));
	t721 = sin(pkin(6));
	t778 = t720 * t721;
	t726 = sin(qJ(4));
	t777 = t720 * t726;
	t776 = t720 * t727;
	t730 = cos(qJ(4));
	t775 = t720 * t730;
	t774 = t720 * t731;
	t773 = t721 * t722;
	t772 = t721 * t723;
	t771 = t723 * t727;
	t770 = t723 * t731;
	t766 = t727 * t732;
	t765 = t728 * t731;
	t763 = qJD(4) * t726;
	t762 = qJD(4) * t730;
	t725 = sin(qJ(6));
	t761 = qJD(6) * t725;
	t760 = qJD(6) * t726;
	t729 = cos(qJ(6));
	t759 = qJD(6) * t729;
	t758 = t728 * t778;
	t757 = t721 * t774;
	t754 = t724 * t774;
	t753 = qJD(2) * t778;
	t752 = qJD(3) * t776;
	t751 = t728 * t753;
	t750 = t732 * t753;
	t704 = t709 * qJD(2);
	t705 = t710 * qJD(2);
	t744 = t709 * t723 - t720 * t773;
	t665 = -t705 * t771 + t704 * t731 + (t744 * t731 - t782) * qJD(3);
	t681 = -t709 * t770 + t722 * t757 + t782;
	t749 = t681 * t760 - t665;
	t742 = t719 * t768 + t722 * t728;
	t706 = t742 * qJD(2);
	t707 = t741 * qJD(2);
	t743 = t719 * t778 - t723 * t742;
	t667 = t707 * t771 - t706 * t731 + (t743 * t731 + t780) * qJD(3);
	t683 = -t719 * t757 + t742 * t770 - t780;
	t748 = t683 * t760 - t667;
	t737 = -t723 * t767 + t764;
	t680 = qJD(3) * t754 + (t737 * qJD(2) + t740 * qJD(3)) * t721;
	t694 = -t740 * t721 - t754;
	t747 = t694 * t760 - t680;
	t682 = t744 * t727 + t781;
	t696 = -t709 * t720 - t722 * t772;
	t674 = t682 * t730 + t696 * t726;
	t673 = t682 * t726 - t696 * t730;
	t684 = t743 * t727 - t731 * t741;
	t697 = t719 * t772 + t720 * t742;
	t676 = t684 * t730 + t697 * t726;
	t675 = t684 * t726 - t697 * t730;
	t738 = t723 * t766 + t765;
	t695 = t738 * t721 + t724 * t776;
	t708 = t724 * t723 - t732 * t778;
	t686 = t695 * t730 + t708 * t726;
	t685 = t695 * t726 - t708 * t730;
	t690 = t709 * t731 - t710 * t771;
	t746 = -t690 * t726 + t710 * t775;
	t692 = -t731 * t742 + t741 * t771;
	t745 = -t692 * t726 - t741 * t775;
	t689 = t709 * t727 + t710 * t770;
	t691 = -t727 * t742 - t741 * t770;
	t739 = t723 * t765 + t766;
	t703 = t737 * t721;
	t736 = -t703 * t726 + t730 * t758;
	t664 = t704 * t727 + t705 * t770 - t752 * t773 + (t709 * t771 + t781) * qJD(3);
	t735 = -qJD(6) * t682 - t664 * t726 - t681 * t762;
	t666 = t684 * qJD(3) - t706 * t727 - t707 * t770;
	t734 = -qJD(6) * t684 - t666 * t726 - t683 * t762;
	t679 = t724 * t752 + (t739 * qJD(2) + t738 * qJD(3)) * t721;
	t733 = -qJD(6) * t695 - t679 * t726 - t694 * t762;
	t702 = t739 * t721;
	t688 = (-t738 * qJD(2) - t739 * qJD(3)) * t721;
	t687 = (t740 * qJD(2) + t737 * qJD(3)) * t721;
	t672 = -t691 * qJD(3) + t706 * t771 + t707 * t731;
	t671 = t692 * qJD(3) - t706 * t770 + t707 * t727;
	t670 = -t689 * qJD(3) - t704 * t771 - t705 * t731;
	t669 = t690 * qJD(3) + t704 * t770 - t705 * t727;
	t668 = -t730 * t750 + t688 * t726 + (t703 * t730 + t726 * t758) * qJD(4);
	t663 = -t685 * qJD(4) + t680 * t730 + t726 * t751;
	t662 = t686 * qJD(4) + t680 * t726 - t730 * t751;
	t661 = t706 * t775 + t672 * t726 + (t692 * t730 - t741 * t777) * qJD(4);
	t660 = -t704 * t775 + t670 * t726 + (t690 * t730 + t710 * t777) * qJD(4);
	t659 = -t675 * qJD(4) + t667 * t730 - t707 * t777;
	t658 = t676 * qJD(4) + t667 * t726 + t707 * t775;
	t657 = -t673 * qJD(4) + t665 * t730 + t705 * t777;
	t656 = t674 * qJD(4) + t665 * t726 - t705 * t775;
	t1 = [0, t661 * t725 + t671 * t729 + (-t691 * t725 - t729 * t745) * qJD(6), t734 * t725 - t748 * t729, t659 * t725 + t676 * t759, 0, t658 * t729 - t666 * t725 + (-t675 * t725 - t683 * t729) * qJD(6); 0, t660 * t725 + t669 * t729 + (-t689 * t725 - t729 * t746) * qJD(6), t735 * t725 - t749 * t729, t657 * t725 + t674 * t759, 0, t656 * t729 - t664 * t725 + (-t673 * t725 - t681 * t729) * qJD(6); 0, t668 * t725 + t687 * t729 + (-t702 * t725 - t729 * t736) * qJD(6), t733 * t725 - t747 * t729, t663 * t725 + t686 * t759, 0, t662 * t729 - t679 * t725 + (-t685 * t725 - t694 * t729) * qJD(6); 0, t661 * t729 - t671 * t725 + (-t691 * t729 + t725 * t745) * qJD(6), t748 * t725 + t734 * t729, t659 * t729 - t676 * t761, 0, -t658 * t725 - t666 * t729 + (-t675 * t729 + t683 * t725) * qJD(6); 0, t660 * t729 - t669 * t725 + (-t689 * t729 + t725 * t746) * qJD(6), t749 * t725 + t735 * t729, t657 * t729 - t674 * t761, 0, -t656 * t725 - t664 * t729 + (-t673 * t729 + t681 * t725) * qJD(6); 0, t668 * t729 - t687 * t725 + (-t702 * t729 + t725 * t736) * qJD(6), t747 * t725 + t733 * t729, t663 * t729 - t686 * t761, 0, -t662 * t725 - t679 * t729 + (-t685 * t729 + t694 * t725) * qJD(6); 0, t745 * qJD(4) + t672 * t730 - t706 * t777, -t666 * t730 + t683 * t763, -t658, 0, 0; 0, t746 * qJD(4) + t670 * t730 + t704 * t777, -t664 * t730 + t681 * t763, -t656, 0, 0; 0, t736 * qJD(4) + t688 * t730 + t726 * t750, -t679 * t730 + t694 * t763, -t662, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end