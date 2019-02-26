% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:25
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR14_jacobiaD_rot_5_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_5_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_5_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_rot_5_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:25:35
% EndTime: 2019-01-03 10:25:40
% DurationCPUTime: 4.60s
% Computational Cost: add. (22947->221), mult. (69993->417), div. (705->12), fcn. (88776->19), ass. (0->197)
t498 = sin(pkin(14));
t616 = cos(pkin(6));
t617 = sin(qJ(2));
t574 = t616 * t617;
t618 = sin(qJ(1));
t620 = cos(qJ(2));
t621 = cos(qJ(1));
t534 = t621 * t574 + t618 * t620;
t614 = cos(pkin(14));
t528 = t534 * t614;
t575 = t616 * t620;
t546 = -t621 * t575 + t618 * t617;
t615 = cos(pkin(7));
t538 = t546 * t615;
t612 = sin(pkin(7));
t613 = sin(pkin(6));
t568 = t613 * t612;
t556 = t621 * t568;
t474 = -t528 + (t538 + t556) * t498;
t502 = sin(qJ(4));
t642 = t474 * t502;
t619 = cos(qJ(4));
t641 = t474 * t619;
t535 = t618 * t575 + t621 * t617;
t482 = t535 * qJD(1) + t534 * qJD(2);
t563 = -t618 * t574 + t621 * t620;
t483 = t563 * qJD(1) - t546 * qJD(2);
t551 = t618 * t568;
t458 = (qJD(1) * t551 - t615 * t482) * t498 + t483 * t614;
t640 = t458 * t502;
t499 = sin(pkin(8));
t570 = t615 * t613;
t557 = t621 * t570;
t624 = t546 * t612 - t557;
t639 = t624 * t499;
t500 = cos(pkin(8));
t529 = t535 * t614;
t567 = t612 * t614;
t549 = t613 * t567;
t542 = t618 * t549;
t517 = t563 * t498 + t615 * t529 - t542;
t552 = t618 * t570;
t523 = t535 * t612 + t552;
t638 = t523 * t499 - t517 * t500;
t571 = t615 * t614;
t459 = -qJD(1) * t542 + t482 * t571 + t483 * t498;
t533 = qJD(1) * t552 + t482 * t612;
t627 = t459 * t500 - t533 * t499;
t637 = t458 * t619 - t627 * t502;
t531 = t534 * t498;
t543 = t621 * t549;
t635 = t531 + t543;
t539 = t546 * t614;
t630 = (-t615 * t539 - t635) * t500 + t639;
t438 = -t630 * t619 - t642;
t436 = t438 ^ 2;
t569 = t614 * t613;
t550 = t615 * t569;
t582 = t498 * t613;
t564 = t617 * t582;
t526 = -t620 * t550 + t564;
t486 = t616 * t567 - t526;
t553 = t617 * t569;
t555 = t620 * t570;
t487 = t553 + (t612 * t616 + t555) * t498;
t554 = t620 * t568;
t493 = t616 * t615 - t554;
t586 = t500 * t619;
t587 = t499 * t619;
t536 = t486 * t586 - t487 * t502 + t493 * t587;
t453 = 0.1e1 / t536 ^ 2;
t428 = t436 * t453 + 0.1e1;
t422 = 0.1e1 / t428;
t440 = t630 * t502 - t641;
t411 = t440 * qJD(4) + t627 * t619 + t640;
t456 = t487 * t619 + (t486 * t500 + t493 * t499) * t502;
t491 = -t617 * t550 - t620 * t582;
t489 = t491 * qJD(2);
t492 = -t615 * t564 + t620 * t569;
t490 = t492 * qJD(2);
t547 = t617 * t499 * t568;
t540 = t619 * t547;
t432 = -qJD(2) * t540 + t456 * qJD(4) - t489 * t586 + t490 * t502;
t452 = 0.1e1 / t536;
t600 = t438 * t453;
t562 = t411 * t452 + t432 * t600;
t393 = t562 * t422;
t429 = atan2(-t438, -t536);
t416 = sin(t429);
t417 = cos(t429);
t566 = t416 * t536 - t417 * t438;
t388 = t566 * t393 - t411 * t416 + t417 * t432;
t405 = -t416 * t438 - t417 * t536;
t403 = 0.1e1 / t405 ^ 2;
t633 = t388 * t403;
t632 = t432 * t453;
t481 = t534 * qJD(1) + t535 * qJD(2);
t520 = t546 * qJD(1) - t563 * qJD(2);
t518 = t520 * t614;
t511 = qJD(1) * t543 + t481 * t498 + t615 * t518;
t513 = qJD(1) * t557 - t520 * t612;
t631 = t499 * t513 + t500 * t511;
t629 = t638 * t619;
t581 = t499 * t612;
t628 = (t546 * t498 - t615 * t528) * t500 + t534 * t581;
t475 = t563 * t614 + (-t535 * t615 + t551) * t498;
t443 = t475 * t502 - t629;
t437 = t443 ^ 2;
t399 = t403 * t437 + 0.1e1;
t397 = 0.1e1 / t399;
t402 = 0.1e1 / t405;
t444 = t475 * t619 + t502 * t638;
t457 = -t481 * t614 + (qJD(1) * t556 + t520 * t615) * t498;
t409 = t444 * qJD(4) + t457 * t502 - t631 * t619;
t604 = t409 * t403;
t610 = t402 * t633;
t611 = (-t437 * t610 + t443 * t604) / t399 ^ 2;
t626 = -t397 * t633 - 0.2e1 * t402 * t611;
t622 = 0.2e1 * t443;
t577 = t610 * t622;
t596 = 0.2e1 * t611;
t605 = t403 * t443;
t625 = t596 * t605 + (t577 - t604) * t397;
t463 = t517 * t499 + t523 * t500;
t501 = sin(qJ(5));
t503 = cos(qJ(5));
t427 = t444 * t503 + t463 * t501;
t419 = 0.1e1 / t427;
t420 = 0.1e1 / t427 ^ 2;
t623 = -0.2e1 * t438;
t598 = qJD(4) * t502;
t410 = t629 * qJD(4) + t457 * t619 - t475 * t598 + t631 * t502;
t434 = -t511 * t499 + t513 * t500;
t400 = t427 * qJD(5) + t410 * t501 - t434 * t503;
t426 = t444 * t501 - t463 * t503;
t418 = t426 ^ 2;
t408 = t418 * t420 + 0.1e1;
t603 = t420 * t426;
t597 = qJD(5) * t426;
t401 = t410 * t503 + t434 * t501 - t597;
t606 = t401 * t419 * t420;
t609 = (t400 * t603 - t418 * t606) / t408 ^ 2;
t602 = t452 * t632;
t608 = (t411 * t600 + t436 * t602) / t428 ^ 2;
t607 = t397 * t402;
t601 = t438 * t452;
t595 = -0.2e1 * t609;
t594 = 0.2e1 * t609;
t593 = -0.2e1 * t608;
t592 = t452 * t608;
t591 = t397 * t605;
t588 = t426 * t606;
t585 = qJD(2) * t619;
t583 = t498 * t615;
t579 = 0.2e1 * t588;
t578 = t602 * t623;
t473 = t614 * t538 + t635;
t442 = t641 + (t473 * t500 - t639) * t502;
t462 = -t473 * t499 - t500 * t624;
t425 = t442 * t503 + t462 * t501;
t424 = t442 * t501 - t462 * t503;
t545 = t563 * t615;
t479 = t535 * t498 - t614 * t545;
t480 = -t498 * t545 - t529;
t544 = t563 * t612;
t541 = t499 * t544;
t449 = t480 * t619 + (t479 * t500 + t541) * t502;
t467 = -t479 * t499 + t500 * t544;
t431 = t449 * t503 + t467 * t501;
t430 = t449 * t501 - t467 * t503;
t565 = t619 * t581;
t561 = -t501 * t419 + t503 * t603;
t560 = t440 * t452 + t456 * t600;
t478 = -t615 * t531 - t539;
t447 = t478 * t502 - t628 * t619;
t466 = -t491 * t586 + t492 * t502 - t540;
t559 = t447 * t452 + t466 * t600;
t548 = -t416 + (-t417 * t601 + t416) * t422;
t537 = t473 * t586 - t587 * t624 - t642;
t448 = -t479 * t586 + t480 * t502 - t619 * t541;
t465 = t481 * t583 + t518;
t464 = t481 * t571 - t520 * t498;
t446 = -t481 * t612 * t500 - t464 * t499;
t445 = (-t498 * t555 - t553) * qJD(2) * t502 + t492 * qJD(4) * t619 - t499 * t554 * t585 + t547 * t598 + (t491 * t598 - t526 * t585) * t500;
t435 = -t459 * t499 - t500 * t533;
t433 = t490 * t619 + (qJD(2) * t547 + t489 * t500) * t502 + t536 * qJD(4);
t415 = (-t482 * t614 - t483 * t583) * t502 - (t482 * t498 - t483 * t571) * t586 - t483 * t565 + (t478 * t619 + t628 * t502) * qJD(4);
t414 = t465 * t619 + (t464 * t500 - t481 * t581) * t502 - t448 * qJD(4);
t413 = t537 * qJD(4) - t637;
t412 = -t438 * qJD(4) + t637;
t406 = 0.1e1 / t408;
t396 = t559 * t422;
t394 = t560 * t422;
t389 = t566 * t394 - t416 * t440 + t417 * t456;
t387 = t559 * t593 + (-t466 * t578 + t415 * t452 + (t411 * t466 + t432 * t447 + t438 * t445) * t453) * t422;
t385 = t560 * t593 + (-t456 * t578 + t412 * t452 + (t411 * t456 + t432 * t440 + t433 * t438) * t453) * t422;
t1 = [-t592 * t622 + (t409 * t452 + t443 * t632) * t422, t387, 0, t385, 0, 0; (t442 * qJD(4) - t459 * t586 + t533 * t587 - t640) * t607 - (t548 * t409 + ((t393 * t422 * t601 + t593) * t416 + (-t592 * t623 - t393 + (t393 - t562) * t422) * t417) * t443) * t591 - t626 * t537 + t625 * t548 * t443 (t449 * qJD(4) - t464 * t586 + t465 * t502 + t481 * t565) * t607 - ((-t387 * t438 - t396 * t411 + t445 + (t396 * t536 - t447) * t393) * t417 + (t387 * t536 - t396 * t432 - t415 + (t396 * t438 - t466) * t393) * t416) * t591 + t626 * t448 + t625 * (t566 * t396 - t416 * t447 + t417 * t466) 0 (t389 * t605 - t402 * t444) * t596 + (t389 * t577 + t410 * t402 + (-t444 * t388 - t389 * t409 + (-(-t385 * t438 - t394 * t411 + t433 + (t394 * t536 - t440) * t393) * t417 - (t385 * t536 - t394 * t432 - t412 + (t394 * t438 - t456) * t393) * t416) * t443) * t403) * t397, 0, 0; (-t419 * t424 + t425 * t603) * t594 + ((t425 * qJD(5) + t413 * t501 - t435 * t503) * t419 + t425 * t579 + (-t424 * t401 - (-t424 * qJD(5) + t413 * t503 + t435 * t501) * t426 - t425 * t400) * t420) * t406 (-t419 * t430 + t431 * t603) * t594 + ((t431 * qJD(5) + t414 * t501 - t446 * t503) * t419 + t431 * t579 + (-t430 * t401 - (-t430 * qJD(5) + t414 * t503 + t446 * t501) * t426 - t431 * t400) * t420) * t406, 0, t561 * t443 * t595 + (t561 * t409 + ((-qJD(5) * t419 - 0.2e1 * t588) * t503 + (t400 * t503 + (t401 - t597) * t501) * t420) * t443) * t406, t595 + 0.2e1 * (t400 * t420 * t406 + (-t406 * t606 - t420 * t609) * t426) * t426, 0;];
JaD_rot  = t1;
