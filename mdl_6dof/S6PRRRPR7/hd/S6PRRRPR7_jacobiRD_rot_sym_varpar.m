% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR7
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
% Datum: 2019-10-09 22:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRRPR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
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
	% StartTime: 2019-10-09 22:58:05
	% EndTime: 2019-10-09 22:58:05
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
	% StartTime: 2019-10-09 22:58:07
	% EndTime: 2019-10-09 22:58:08
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
	% StartTime: 2019-10-09 22:58:10
	% EndTime: 2019-10-09 22:58:10
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (437->102), mult. (1485->207), div. (0->0), fcn. (1639->14), ass. (0->82)
	t586 = sin(pkin(12));
	t590 = cos(pkin(12));
	t595 = sin(qJ(2));
	t592 = cos(pkin(6));
	t598 = cos(qJ(2));
	t622 = t592 * t598;
	t579 = -t586 * t595 + t590 * t622;
	t594 = sin(qJ(3));
	t597 = cos(qJ(3));
	t623 = t592 * t595;
	t606 = t586 * t623 - t590 * t598;
	t591 = cos(pkin(7));
	t607 = t586 * t622 + t590 * t595;
	t587 = sin(pkin(7));
	t588 = sin(pkin(6));
	t630 = t587 * t588;
	t608 = t586 * t630 - t591 * t607;
	t564 = t608 * t594 - t597 * t606;
	t580 = t586 * t598 + t590 * t623;
	t609 = -t579 * t591 + t590 * t630;
	t634 = -t580 * t597 + t609 * t594;
	t629 = t587 * t592;
	t593 = sin(qJ(4));
	t628 = t587 * t593;
	t596 = cos(qJ(4));
	t627 = t587 * t596;
	t626 = t588 * t591;
	t625 = t591 * t594;
	t624 = t591 * t597;
	t621 = t594 * t595;
	t620 = t594 * t598;
	t619 = t595 * t597;
	t618 = t597 * t598;
	t617 = qJD(4) * t593;
	t616 = qJD(4) * t596;
	t615 = t595 * t630;
	t613 = qJD(2) * t630;
	t612 = qJD(3) * t629;
	t611 = t595 * t613;
	t610 = t598 * t613;
	t567 = t579 * t597 - t580 * t625;
	t568 = -t597 * t607 + t606 * t625;
	t605 = t591 * t618 - t621;
	t604 = -t591 * t619 - t620;
	t603 = t591 * t620 + t619;
	t602 = -t591 * t621 + t618;
	t574 = t579 * qJD(2);
	t575 = t580 * qJD(2);
	t550 = t634 * qJD(3) - t574 * t594 - t575 * t624;
	t561 = -t580 * t594 - t609 * t597;
	t601 = -t550 * t596 + t561 * t617;
	t576 = t607 * qJD(2);
	t577 = t606 * qJD(2);
	t552 = -t564 * qJD(3) + t576 * t594 + t577 * t624;
	t563 = t594 * t606 + t608 * t597;
	t600 = -t552 * t596 + t563 * t617;
	t559 = -t594 * t612 + (t604 * qJD(2) - t603 * qJD(3)) * t588;
	t569 = t605 * t588 + t597 * t629;
	t599 = -t559 * t596 + t569 * t617;
	t589 = cos(pkin(13));
	t585 = sin(pkin(13));
	t578 = t592 * t591 - t598 * t630;
	t573 = t602 * t588;
	t572 = t586 * t626 + t587 * t607;
	t571 = -t579 * t587 - t590 * t626;
	t570 = t603 * t588 + t594 * t629;
	t566 = (-t603 * qJD(2) + t604 * qJD(3)) * t588;
	t565 = (t605 * qJD(2) + t602 * qJD(3)) * t588;
	t560 = t597 * t612 + (t602 * qJD(2) + t605 * qJD(3)) * t588;
	t558 = t576 * t625 + t577 * t597 + (t594 * t607 + t606 * t624) * qJD(3);
	t557 = t568 * qJD(3) - t576 * t624 + t577 * t594;
	t556 = -t574 * t625 - t575 * t597 + (-t579 * t594 - t580 * t624) * qJD(3);
	t555 = t567 * qJD(3) + t574 * t624 - t575 * t594;
	t554 = t593 * t610 + t566 * t596 + (-t573 * t593 + t596 * t615) * qJD(4);
	t553 = t563 * qJD(3) - t576 * t597 + t577 * t625;
	t551 = t561 * qJD(3) + t574 * t597 - t575 * t625;
	t549 = t596 * t611 - t560 * t593 + (-t570 * t596 - t578 * t593) * qJD(4);
	t548 = -t576 * t628 + t558 * t596 + (-t568 * t593 - t606 * t627) * qJD(4);
	t547 = t574 * t628 + t556 * t596 + (-t567 * t593 + t580 * t627) * qJD(4);
	t546 = -t577 * t627 - t553 * t593 + (-t564 * t596 - t572 * t593) * qJD(4);
	t545 = t575 * t627 - t551 * t593 + (-t571 * t593 + t596 * t634) * qJD(4);
	t1 = [0, t548 * t589 + t557 * t585, t553 * t585 - t600 * t589, t546 * t589, 0, 0; 0, t547 * t589 + t555 * t585, t551 * t585 - t601 * t589, t545 * t589, 0, 0; 0, t554 * t589 + t565 * t585, t560 * t585 - t599 * t589, t549 * t589, 0, 0; 0, -t548 * t585 + t557 * t589, t553 * t589 + t600 * t585, -t546 * t585, 0, 0; 0, -t547 * t585 + t555 * t589, t551 * t589 + t601 * t585, -t545 * t585, 0, 0; 0, -t554 * t585 + t565 * t589, t560 * t589 + t599 * t585, -t549 * t585, 0, 0; 0, t576 * t627 + t558 * t593 + (t568 * t596 - t606 * t628) * qJD(4), t552 * t593 + t563 * t616, -t577 * t628 + t553 * t596 + (-t564 * t593 + t572 * t596) * qJD(4), 0, 0; 0, -t574 * t627 + t556 * t593 + (t567 * t596 + t580 * t628) * qJD(4), t550 * t593 + t561 * t616, t575 * t628 + t551 * t596 + (t571 * t596 + t593 * t634) * qJD(4), 0, 0; 0, -t596 * t610 + t566 * t593 + (t573 * t596 + t593 * t615) * qJD(4), t559 * t593 + t569 * t616, t593 * t611 + t560 * t596 + (-t570 * t593 + t578 * t596) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:10
	% EndTime: 2019-10-09 22:58:12
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (874->148), mult. (2555->280), div. (0->0), fcn. (2925->14), ass. (0->111)
	t698 = cos(pkin(7));
	t704 = cos(qJ(3));
	t705 = cos(qJ(2));
	t734 = t704 * t705;
	t701 = sin(qJ(3));
	t702 = sin(qJ(2));
	t737 = t701 * t702;
	t712 = t698 * t734 - t737;
	t694 = sin(pkin(12));
	t697 = cos(pkin(12));
	t699 = cos(pkin(6));
	t738 = t699 * t705;
	t681 = -t694 * t702 + t697 * t738;
	t739 = t699 * t702;
	t682 = t694 * t705 + t697 * t739;
	t752 = t682 * t701;
	t751 = t682 * t704;
	t713 = t694 * t739 - t697 * t705;
	t750 = t713 * t701;
	t695 = sin(pkin(7));
	t696 = sin(pkin(6));
	t748 = t695 * t696;
	t700 = sin(qJ(4));
	t747 = t695 * t700;
	t746 = t695 * t701;
	t703 = cos(qJ(4));
	t745 = t695 * t703;
	t744 = t695 * t704;
	t743 = t696 * t697;
	t742 = t696 * t698;
	t741 = t698 * t701;
	t740 = t698 * t704;
	t736 = t701 * t705;
	t735 = t702 * t704;
	t733 = qJD(4) * t700;
	t732 = qJD(4) * t703;
	t693 = pkin(13) + qJ(6);
	t691 = sin(t693);
	t731 = qJD(6) * t691;
	t692 = cos(t693);
	t730 = qJD(6) * t692;
	t729 = qJD(6) * t703;
	t728 = t702 * t748;
	t727 = t696 * t744;
	t724 = t699 * t744;
	t723 = qJD(2) * t748;
	t722 = qJD(3) * t746;
	t721 = t702 * t723;
	t720 = t705 * t723;
	t676 = t681 * qJD(2);
	t677 = t682 * qJD(2);
	t716 = t681 * t698 - t695 * t743;
	t637 = -t677 * t741 + t676 * t704 + (t716 * t704 - t752) * qJD(3);
	t653 = -t681 * t740 + t697 * t727 + t752;
	t719 = t653 * t729 + t637;
	t714 = t694 * t738 + t697 * t702;
	t678 = t714 * qJD(2);
	t679 = t713 * qJD(2);
	t715 = t694 * t748 - t698 * t714;
	t639 = t679 * t741 - t678 * t704 + (t715 * t704 + t750) * qJD(3);
	t655 = -t694 * t727 + t714 * t740 - t750;
	t718 = t655 * t729 + t639;
	t709 = -t698 * t737 + t734;
	t652 = qJD(3) * t724 + (t709 * qJD(2) + t712 * qJD(3)) * t696;
	t666 = -t712 * t696 - t724;
	t717 = t666 * t729 + t652;
	t654 = t716 * t701 + t751;
	t668 = -t681 * t695 - t697 * t742;
	t646 = t654 * t703 + t668 * t700;
	t645 = -t654 * t700 + t668 * t703;
	t656 = t715 * t701 - t704 * t713;
	t669 = t694 * t742 + t695 * t714;
	t648 = t656 * t703 + t669 * t700;
	t647 = -t656 * t700 + t669 * t703;
	t710 = t698 * t736 + t735;
	t667 = t710 * t696 + t699 * t746;
	t680 = t699 * t698 - t705 * t748;
	t658 = t667 * t703 + t680 * t700;
	t657 = -t667 * t700 + t680 * t703;
	t662 = t681 * t704 - t682 * t741;
	t649 = t662 * t703 + t682 * t747;
	t664 = -t704 * t714 + t713 * t741;
	t650 = t664 * t703 - t713 * t747;
	t661 = t681 * t701 + t682 * t740;
	t663 = -t701 * t714 - t713 * t740;
	t711 = t698 * t735 + t736;
	t675 = t709 * t696;
	t665 = t675 * t703 + t700 * t728;
	t636 = t676 * t701 + t677 * t740 - t722 * t743 + (t681 * t741 + t751) * qJD(3);
	t708 = qJD(6) * t654 - t636 * t703 + t653 * t733;
	t638 = t656 * qJD(3) - t678 * t701 - t679 * t740;
	t707 = qJD(6) * t656 - t638 * t703 + t655 * t733;
	t651 = t699 * t722 + (t711 * qJD(2) + t710 * qJD(3)) * t696;
	t706 = qJD(6) * t667 - t651 * t703 + t666 * t733;
	t674 = t711 * t696;
	t660 = (-t710 * qJD(2) - t711 * qJD(3)) * t696;
	t659 = (t712 * qJD(2) + t709 * qJD(3)) * t696;
	t644 = -t663 * qJD(3) + t678 * t741 + t679 * t704;
	t643 = t664 * qJD(3) - t678 * t740 + t679 * t701;
	t642 = -t661 * qJD(3) - t676 * t741 - t677 * t704;
	t641 = t662 * qJD(3) + t676 * t740 - t677 * t701;
	t640 = t700 * t720 + t660 * t703 + (-t675 * t700 + t703 * t728) * qJD(4);
	t635 = t657 * qJD(4) + t652 * t703 + t700 * t721;
	t634 = -t658 * qJD(4) - t652 * t700 + t703 * t721;
	t633 = -t678 * t747 + t644 * t703 + (-t664 * t700 - t713 * t745) * qJD(4);
	t632 = t676 * t747 + t642 * t703 + (-t662 * t700 + t682 * t745) * qJD(4);
	t631 = t647 * qJD(4) + t639 * t703 - t679 * t747;
	t630 = -t648 * qJD(4) - t639 * t700 - t679 * t745;
	t629 = t645 * qJD(4) + t637 * t703 + t677 * t747;
	t628 = -t646 * qJD(4) - t637 * t700 + t677 * t745;
	t1 = [0, t633 * t692 + t643 * t691 + (-t650 * t691 + t663 * t692) * qJD(6), t718 * t691 + t707 * t692, t630 * t692 - t647 * t731, 0, -t631 * t691 + t638 * t692 + (-t648 * t692 - t655 * t691) * qJD(6); 0, t632 * t692 + t641 * t691 + (-t649 * t691 + t661 * t692) * qJD(6), t719 * t691 + t708 * t692, t628 * t692 - t645 * t731, 0, -t629 * t691 + t636 * t692 + (-t646 * t692 - t653 * t691) * qJD(6); 0, t640 * t692 + t659 * t691 + (-t665 * t691 + t674 * t692) * qJD(6), t717 * t691 + t706 * t692, t634 * t692 - t657 * t731, 0, -t635 * t691 + t651 * t692 + (-t658 * t692 - t666 * t691) * qJD(6); 0, -t633 * t691 + t643 * t692 + (-t650 * t692 - t663 * t691) * qJD(6), -t707 * t691 + t718 * t692, -t630 * t691 - t647 * t730, 0, -t631 * t692 - t638 * t691 + (t648 * t691 - t655 * t692) * qJD(6); 0, -t632 * t691 + t641 * t692 + (-t649 * t692 - t661 * t691) * qJD(6), -t708 * t691 + t719 * t692, -t628 * t691 - t645 * t730, 0, -t629 * t692 - t636 * t691 + (t646 * t691 - t653 * t692) * qJD(6); 0, -t640 * t691 + t659 * t692 + (-t665 * t692 - t674 * t691) * qJD(6), -t706 * t691 + t717 * t692, -t634 * t691 - t657 * t730, 0, -t635 * t692 - t651 * t691 + (t658 * t691 - t666 * t692) * qJD(6); 0, t650 * qJD(4) + t644 * t700 + t678 * t745, -t638 * t700 - t655 * t732, t631, 0, 0; 0, t649 * qJD(4) + t642 * t700 - t676 * t745, -t636 * t700 - t653 * t732, t629, 0, 0; 0, t665 * qJD(4) + t660 * t700 - t703 * t720, -t651 * t700 - t666 * t732, t635, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end