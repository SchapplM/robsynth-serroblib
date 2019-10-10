% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:19
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRRRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(13));
	t58 = sin(pkin(13));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:24
	% EndTime: 2019-10-09 23:19:25
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (78->41), mult. (285->88), div. (0->0), fcn. (301->10), ass. (0->35)
	t259 = sin(pkin(13));
	t262 = cos(pkin(13));
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
	% StartTime: 2019-10-09 23:19:26
	% EndTime: 2019-10-09 23:19:27
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (269->87), mult. (921->177), div. (0->0), fcn. (1021->12), ass. (0->68)
	t472 = sin(pkin(13));
	t475 = cos(pkin(13));
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
	% StartTime: 2019-10-09 23:19:27
	% EndTime: 2019-10-09 23:19:27
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (487->77), mult. (1225->156), div. (0->0), fcn. (1363->12), ass. (0->76)
	t536 = sin(pkin(13));
	t539 = cos(pkin(13));
	t543 = sin(qJ(2));
	t541 = cos(pkin(6));
	t545 = cos(qJ(2));
	t573 = t541 * t545;
	t526 = -t536 * t543 + t539 * t573;
	t542 = sin(qJ(3));
	t544 = cos(qJ(3));
	t574 = t541 * t543;
	t552 = t536 * t574 - t539 * t545;
	t540 = cos(pkin(7));
	t553 = t536 * t573 + t539 * t543;
	t537 = sin(pkin(7));
	t538 = sin(pkin(6));
	t579 = t537 * t538;
	t554 = t536 * t579 - t540 * t553;
	t587 = t554 * t542 - t544 * t552;
	t527 = t536 * t545 + t539 * t574;
	t555 = -t526 * t540 + t539 * t579;
	t586 = -t527 * t544 + t555 * t542;
	t535 = qJ(4) + qJ(5);
	t532 = sin(t535);
	t534 = qJD(4) + qJD(5);
	t583 = t532 * t534;
	t533 = cos(t535);
	t582 = t533 * t534;
	t581 = t534 * t537;
	t578 = t537 * t541;
	t577 = t538 * t540;
	t576 = t540 * t542;
	t575 = t540 * t544;
	t572 = t542 * t543;
	t571 = t542 * t545;
	t570 = t543 * t544;
	t569 = t544 * t545;
	t567 = qJD(2) * t579;
	t566 = qJD(3) * t578;
	t509 = -t527 * t542 - t555 * t544;
	t521 = t526 * qJD(2);
	t522 = t527 * qJD(2);
	t502 = t509 * qJD(3) + t521 * t544 - t522 * t576;
	t565 = -(-t526 * t537 - t539 * t577) * t534 - t502;
	t511 = t542 * t552 + t554 * t544;
	t523 = t553 * qJD(2);
	t524 = t552 * qJD(2);
	t504 = t511 * qJD(3) - t523 * t544 + t524 * t576;
	t564 = -(t536 * t577 + t537 * t553) * t534 - t504;
	t548 = -t540 * t572 + t569;
	t551 = t540 * t569 - t572;
	t508 = t544 * t566 + (t548 * qJD(2) + t551 * qJD(3)) * t538;
	t563 = -(t541 * t540 - t545 * t579) * t534 - t508;
	t562 = t527 * t581 - t521 * t576 - t522 * t544 + (-t526 * t542 - t527 * t575) * qJD(3);
	t561 = -t552 * t581 + t523 * t576 + t524 * t544 + (t542 * t553 + t552 * t575) * qJD(3);
	t560 = -t522 * t537 - t586 * t534;
	t559 = t524 * t537 + t587 * t534;
	t514 = t526 * t544 - t527 * t576;
	t558 = -t514 * t534 + t521 * t537;
	t515 = -t544 * t553 + t552 * t576;
	t557 = -t515 * t534 - t523 * t537;
	t549 = t540 * t571 + t570;
	t550 = -t540 * t570 - t571;
	t556 = t534 * t543 * t579 + (-t549 * qJD(2) + t550 * qJD(3)) * t538;
	t547 = -(t549 * t538 + t542 * t578) * t534 + t543 * t567;
	t546 = -t534 * t538 * t548 + t545 * t567;
	t516 = t551 * t538 + t544 * t578;
	t507 = -t542 * t566 + (t550 * qJD(2) - t549 * qJD(3)) * t538;
	t503 = -t587 * qJD(3) + t523 * t542 + t524 * t575;
	t501 = t586 * qJD(3) - t521 * t542 - t522 * t575;
	t500 = -t547 * t532 + t563 * t533;
	t499 = t563 * t532 + t547 * t533;
	t498 = t559 * t532 + t564 * t533;
	t497 = t564 * t532 - t559 * t533;
	t496 = t560 * t532 + t565 * t533;
	t495 = t565 * t532 - t560 * t533;
	t1 = [0, t557 * t532 + t561 * t533, t503 * t533 - t511 * t583, t497, t497, 0; 0, t558 * t532 + t562 * t533, t501 * t533 - t509 * t583, t495, t495, 0; 0, t546 * t532 + t556 * t533, t507 * t533 - t516 * t583, t499, t499, 0; 0, -t561 * t532 + t557 * t533, -t503 * t532 - t511 * t582, t498, t498, 0; 0, -t532 * t562 + t558 * t533, -t501 * t532 - t509 * t582, t496, t496, 0; 0, -t532 * t556 + t546 * t533, -t507 * t532 - t516 * t582, t500, t500, 0; 0, qJD(3) * t515 - t523 * t575 + t524 * t542, t504, 0, 0, 0; 0, qJD(3) * t514 + t521 * t575 - t522 * t542, t502, 0, 0, 0; 0, (qJD(2) * t551 + qJD(3) * t548) * t538, t508, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:32
	% EndTime: 2019-10-09 23:19:33
	% DurationCPUTime: 1.17s
	% Computational Cost: add. (1243->146), mult. (3133->278), div. (0->0), fcn. (3596->14), ass. (0->125)
	t773 = cos(pkin(7));
	t779 = cos(qJ(3));
	t780 = cos(qJ(2));
	t817 = t779 * t780;
	t776 = sin(qJ(3));
	t777 = sin(qJ(2));
	t820 = t776 * t777;
	t789 = t773 * t817 - t820;
	t769 = sin(pkin(13));
	t772 = cos(pkin(13));
	t774 = cos(pkin(6));
	t821 = t774 * t780;
	t755 = -t769 * t777 + t772 * t821;
	t822 = t774 * t777;
	t756 = t769 * t780 + t772 * t822;
	t837 = t756 * t776;
	t836 = t756 * t779;
	t790 = t769 * t822 - t772 * t780;
	t835 = t790 * t776;
	t768 = qJ(4) + qJ(5);
	t765 = sin(t768);
	t767 = qJD(4) + qJD(5);
	t834 = t765 * t767;
	t770 = sin(pkin(7));
	t833 = t765 * t770;
	t766 = cos(t768);
	t832 = t766 * t767;
	t831 = t767 * t770;
	t771 = sin(pkin(6));
	t829 = t770 * t771;
	t828 = t770 * t776;
	t827 = t770 * t779;
	t826 = t771 * t772;
	t825 = t771 * t773;
	t824 = t773 * t776;
	t823 = t773 * t779;
	t819 = t776 * t780;
	t818 = t777 * t779;
	t816 = qJD(6) * t766;
	t775 = sin(qJ(6));
	t815 = qJD(6) * t775;
	t778 = cos(qJ(6));
	t814 = qJD(6) * t778;
	t813 = t777 * t829;
	t812 = t771 * t827;
	t809 = t774 * t827;
	t808 = qJD(2) * t829;
	t807 = qJD(3) * t828;
	t750 = t755 * qJD(2);
	t751 = t756 * qJD(2);
	t793 = t755 * t773 - t770 * t826;
	t712 = -t751 * t824 + t750 * t779 + (t793 * t779 - t837) * qJD(3);
	t742 = -t755 * t770 - t772 * t825;
	t806 = t742 * t767 + t712;
	t791 = t769 * t821 + t772 * t777;
	t752 = t791 * qJD(2);
	t753 = t790 * qJD(2);
	t792 = t769 * t829 - t773 * t791;
	t714 = t753 * t824 - t752 * t779 + (t792 * t779 + t835) * qJD(3);
	t743 = t769 * t825 + t770 * t791;
	t805 = t743 * t767 + t714;
	t786 = -t773 * t820 + t817;
	t726 = qJD(3) * t809 + (t786 * qJD(2) + t789 * qJD(3)) * t771;
	t754 = t774 * t773 - t780 * t829;
	t804 = t754 * t767 + t726;
	t735 = t755 * t776 + t756 * t823;
	t803 = -t735 * qJD(3) - t750 * t824 - t751 * t779 + t756 * t831;
	t737 = -t776 * t791 - t790 * t823;
	t802 = -t737 * qJD(3) + t752 * t824 + t753 * t779 - t790 * t831;
	t729 = -t755 * t823 + t772 * t812 + t837;
	t801 = t729 * t816 + t712;
	t731 = -t769 * t812 + t791 * t823 - t835;
	t800 = t731 * t816 + t714;
	t740 = -t789 * t771 - t809;
	t799 = t740 * t816 + t726;
	t730 = t793 * t776 + t836;
	t798 = -t730 * t767 + t751 * t770;
	t732 = t792 * t776 - t779 * t790;
	t797 = -t732 * t767 - t753 * t770;
	t736 = t755 * t779 - t756 * t824;
	t796 = t736 * t767 - t750 * t770;
	t738 = -t779 * t791 + t790 * t824;
	t795 = t738 * t767 + t752 * t770;
	t787 = t773 * t819 + t818;
	t788 = t773 * t818 + t819;
	t794 = t767 * t813 + (-t787 * qJD(2) - t788 * qJD(3)) * t771;
	t741 = t787 * t771 + t774 * t828;
	t785 = -t741 * t767 + t777 * t808;
	t749 = t786 * t771;
	t784 = -t749 * t767 + t780 * t808;
	t711 = t750 * t776 + t751 * t823 - t807 * t826 + (t755 * t824 + t836) * qJD(3);
	t783 = qJD(6) * t730 - t711 * t766 + t729 * t834;
	t713 = t732 * qJD(3) - t752 * t776 - t753 * t823;
	t782 = qJD(6) * t732 - t713 * t766 + t731 * t834;
	t725 = t774 * t807 + (t788 * qJD(2) + t787 * qJD(3)) * t771;
	t781 = qJD(6) * t741 - t725 * t766 + t740 * t834;
	t748 = t788 * t771;
	t739 = t749 * t766 + t765 * t813;
	t733 = (t789 * qJD(2) + t786 * qJD(3)) * t771;
	t728 = t741 * t766 + t754 * t765;
	t727 = -t741 * t765 + t754 * t766;
	t724 = t738 * t766 - t790 * t833;
	t723 = t736 * t766 + t756 * t833;
	t721 = t738 * qJD(3) - t752 * t823 + t753 * t776;
	t719 = t736 * qJD(3) + t750 * t823 - t751 * t776;
	t718 = t732 * t766 + t743 * t765;
	t717 = -t732 * t765 + t743 * t766;
	t716 = t730 * t766 + t742 * t765;
	t715 = -t730 * t765 + t742 * t766;
	t710 = t784 * t765 + t794 * t766;
	t709 = t785 * t765 + t804 * t766;
	t708 = -t804 * t765 + t785 * t766;
	t707 = -t795 * t765 + t802 * t766;
	t706 = -t796 * t765 + t803 * t766;
	t705 = t708 * t778 - t727 * t815;
	t704 = -t708 * t775 - t727 * t814;
	t703 = t797 * t765 + t805 * t766;
	t702 = -t805 * t765 + t797 * t766;
	t701 = t798 * t765 + t806 * t766;
	t700 = -t806 * t765 + t798 * t766;
	t699 = t702 * t778 - t717 * t815;
	t698 = -t702 * t775 - t717 * t814;
	t697 = t700 * t778 - t715 * t815;
	t696 = -t700 * t775 - t715 * t814;
	t1 = [0, t707 * t778 + t721 * t775 + (-t724 * t775 + t737 * t778) * qJD(6), t800 * t775 + t782 * t778, t699, t699, -t703 * t775 + t713 * t778 + (-t718 * t778 - t731 * t775) * qJD(6); 0, t706 * t778 + t719 * t775 + (-t723 * t775 + t735 * t778) * qJD(6), t801 * t775 + t783 * t778, t697, t697, -t701 * t775 + t711 * t778 + (-t716 * t778 - t729 * t775) * qJD(6); 0, t710 * t778 + t733 * t775 + (-t739 * t775 + t748 * t778) * qJD(6), t799 * t775 + t781 * t778, t705, t705, -t709 * t775 + t725 * t778 + (-t728 * t778 - t740 * t775) * qJD(6); 0, -t707 * t775 + t721 * t778 + (-t724 * t778 - t737 * t775) * qJD(6), -t782 * t775 + t800 * t778, t698, t698, -t703 * t778 - t713 * t775 + (t718 * t775 - t731 * t778) * qJD(6); 0, -t706 * t775 + t719 * t778 + (-t723 * t778 - t735 * t775) * qJD(6), -t783 * t775 + t801 * t778, t696, t696, -t701 * t778 - t711 * t775 + (t716 * t775 - t729 * t778) * qJD(6); 0, -t710 * t775 + t733 * t778 + (-t739 * t778 - t748 * t775) * qJD(6), -t781 * t775 + t799 * t778, t704, t704, -t709 * t778 - t725 * t775 + (t728 * t775 - t740 * t778) * qJD(6); 0, t802 * t765 + t795 * t766, -t713 * t765 - t731 * t832, t703, t703, 0; 0, t803 * t765 + t796 * t766, -t711 * t765 - t729 * t832, t701, t701, 0; 0, t794 * t765 - t784 * t766, -t725 * t765 - t740 * t832, t709, t709, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end