% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR10_jacobiaD_rot_4_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_rot_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_rot_4_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiaD_rot_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:19
% EndTime: 2018-11-23 11:27:23
% DurationCPUTime: 4.29s
% Computational Cost: add. (47732->253), mult. (52968->440), div. (706->12), fcn. (50696->29), ass. (0->187)
t550 = sin(pkin(8));
t629 = pkin(6) - qJ(2);
t608 = sin(t629);
t590 = t608 / 0.2e1;
t628 = pkin(6) + qJ(2);
t607 = sin(t628);
t510 = (t590 - t607 / 0.2e1) * qJD(2);
t555 = sin(qJ(1));
t559 = cos(qJ(1));
t595 = cos(t628) / 0.2e1;
t611 = cos(t629);
t524 = t611 / 0.2e1 + t595;
t655 = sin(qJ(2));
t577 = t524 * t555 + t559 * t655;
t558 = cos(qJ(2));
t634 = qJD(2) * t558;
t488 = qJD(1) * t577 - t510 * t559 + t555 * t634;
t551 = sin(pkin(6));
t654 = cos(pkin(7));
t612 = t551 * t654;
t597 = qJD(1) * t612;
t652 = sin(pkin(7));
t569 = t488 * t652 + t555 * t597;
t653 = cos(pkin(8));
t511 = t524 * qJD(2);
t589 = t607 / 0.2e1;
t518 = t589 - t608 / 0.2e1;
t616 = qJD(2) * t655;
t536 = t555 * t616;
t636 = qJD(1) * t555;
t489 = -t518 * t636 - t536 + (qJD(1) * t558 + t511) * t559;
t495 = t518 * t559 + t555 * t558;
t627 = pkin(7) - qJ(3);
t606 = sin(t627);
t588 = t606 / 0.2e1;
t626 = pkin(7) + qJ(3);
t605 = sin(t626);
t506 = (t588 - t605 / 0.2e1) * qJD(3);
t593 = cos(t626) / 0.2e1;
t610 = cos(t627);
t521 = t593 - t610 / 0.2e1;
t508 = t521 * qJD(3);
t587 = t605 / 0.2e1;
t515 = t587 + t588;
t522 = t610 / 0.2e1 + t593;
t554 = sin(qJ(3));
t557 = cos(qJ(3));
t632 = qJD(3) * t557;
t660 = -t524 * t559 + t555 * t655;
t656 = -t488 * t522 - t489 * t554 - t495 * t632 - t660 * t506 + (-t508 * t559 + t515 * t636) * t551;
t425 = t550 * t656 - t569 * t653;
t494 = t559 * t612 - t652 * t660;
t638 = t551 * t559;
t657 = -t495 * t554 - t515 * t638 - t522 * t660;
t666 = t494 * t653 + t550 * t657;
t456 = t666 ^ 2;
t517 = t589 + t590;
t523 = t595 - t611 / 0.2e1;
t552 = cos(pkin(6));
t470 = -(t515 * t552 + t517 * t522 + t523 * t554) * t550 + (-t517 * t652 + t552 * t654) * t653;
t468 = 0.1e1 / t470 ^ 2;
t447 = t456 * t468 + 0.1e1;
t640 = t666 * t468;
t467 = 0.1e1 / t470;
t509 = t517 * qJD(2);
t512 = t523 * qJD(2);
t584 = t653 * t652;
t455 = -(t506 * t517 + t508 * t552 - t509 * t554 + t512 * t522 + t523 * t632) * t550 - t512 * t584;
t662 = t455 * t468;
t642 = t467 * t662;
t649 = (t425 * t640 - t456 * t642) / t447 ^ 2;
t664 = -0.2e1 * t649;
t499 = t518 * t555 - t558 * t559;
t639 = t551 * t555;
t475 = t499 * t554 + t515 * t639 - t522 * t577;
t516 = t587 - t606 / 0.2e1;
t477 = t499 * t557 + t516 * t577 + t521 * t639;
t624 = pkin(8) + qJ(4);
t603 = sin(t624);
t585 = t603 / 0.2e1;
t625 = pkin(8) - qJ(4);
t604 = sin(t625);
t586 = t604 / 0.2e1;
t513 = t585 + t586;
t591 = cos(t624) / 0.2e1;
t609 = cos(t625);
t520 = t609 / 0.2e1 + t591;
t553 = sin(qJ(4));
t571 = -t555 * t612 - t577 * t652;
t440 = -t475 * t520 - t477 * t553 + t513 * t571;
t434 = t440 ^ 2;
t514 = t585 - t604 / 0.2e1;
t519 = t591 - t609 / 0.2e1;
t556 = cos(qJ(4));
t567 = t475 * t514 - t477 * t556 + t519 * t571;
t436 = 0.1e1 / t567 ^ 2;
t663 = t434 * t436;
t505 = t515 * qJD(3);
t507 = t522 * qJD(3);
t633 = qJD(3) * t554;
t659 = t660 * t507 + t488 * t516 + (t505 * t559 + t521 * t636) * t551 + t495 * t633 - t489 * t557;
t658 = -t495 * t557 + t516 * t660 - t521 * t638;
t485 = qJD(1) * t660 - t510 * t555 - t559 * t634;
t486 = qJD(1) * t495 + t511 * t555 + t559 * t616;
t635 = qJD(1) * t559;
t430 = -t485 * t516 + t486 * t557 + t577 * t507 - t499 * t633 + (-t505 * t555 + t521 * t635) * t551;
t448 = atan2(t666, t470);
t443 = sin(t448);
t444 = cos(t448);
t420 = t443 * t666 + t444 * t470;
t417 = 0.1e1 / t420;
t435 = 0.1e1 / t567;
t418 = 0.1e1 / t420 ^ 2;
t460 = t475 * t550 + t571 * t653;
t457 = t460 ^ 2;
t416 = t418 * t457 + 0.1e1;
t428 = t499 * t632 + t485 * t522 + t486 * t554 - t577 * t506 + (t508 * t555 + t515 * t635) * t551;
t570 = t485 * t652 - t559 * t597;
t424 = t428 * t550 + t570 * t653;
t646 = t424 * t418;
t445 = 0.1e1 / t447;
t582 = t425 * t467 - t455 * t640;
t408 = t582 * t445;
t583 = -t443 * t470 + t444 * t666;
t403 = t408 * t583 + t425 * t443 + t444 * t455;
t650 = t403 * t417 * t418;
t651 = (-t457 * t650 + t460 * t646) / t416 ^ 2;
t501 = t513 * qJD(4);
t503 = t520 * qJD(4);
t631 = qJD(4) * t553;
t413 = t428 * t514 - t430 * t556 + t475 * t503 + t477 * t631 - t501 * t571 + t519 * t570;
t437 = t435 * t436;
t648 = t413 * t437;
t647 = t418 * t460;
t645 = t436 * t440;
t644 = t443 * t460;
t643 = t444 * t460;
t641 = t666 * t467;
t630 = qJD(4) * t556;
t623 = -0.2e1 * t651;
t622 = 0.2e1 * t651;
t621 = -0.2e1 * t650;
t423 = 0.1e1 + t663;
t502 = (t586 - t603 / 0.2e1) * qJD(4);
t504 = t519 * qJD(4);
t412 = -t428 * t520 - t430 * t553 - t475 * t502 - t477 * t630 + t504 * t571 + t513 * t570;
t618 = t412 * t645;
t620 = 0.2e1 * (-t434 * t648 + t618) / t423 ^ 2;
t619 = 0.2e1 * t649;
t617 = t666 * t642;
t615 = t499 * t652;
t614 = t513 * t652;
t613 = t519 * t652;
t602 = t467 * t664;
t601 = 0.2e1 * t440 * t648;
t600 = t460 * t621;
t464 = (-t495 * t522 + t554 * t660) * t550 - t495 * t584;
t481 = -(-t517 * t554 + t522 * t523) * t550 - t523 * t584;
t581 = -t464 * t467 + t481 * t640;
t484 = -t516 * t517 + t521 * t552 + t523 * t557;
t580 = t467 * t658 + t484 * t640;
t573 = t443 + (t444 * t641 - t443) * t445;
t490 = qJD(1) * t499 - t511 * t559 + t536;
t483 = t499 * t516 - t557 * t577;
t482 = t499 * t522 + t554 * t577;
t465 = -t482 * t550 - t499 * t584;
t463 = -t505 * t552 - t507 * t517 - t509 * t557 - t512 * t516 - t523 * t633;
t462 = -(t506 * t523 - t509 * t522 - t512 * t554 - t517 * t632) * t550 + t509 * t584;
t454 = t475 * t556 + t477 * t514;
t453 = t475 * t553 - t477 * t520;
t452 = t482 * t514 + t483 * t556 + t499 * t613;
t451 = -t482 * t520 + t483 * t553 + t499 * t614;
t450 = t485 * t557 + t486 * t516 + t499 * t507 + t577 * t633;
t449 = -t485 * t554 + t486 * t522 + t499 * t506 + t577 * t632;
t439 = -t494 * t519 - t514 * t657 + t556 * t658;
t438 = -t494 * t513 + t520 * t657 + t553 * t658;
t426 = (t488 * t554 + t490 * t522 - t495 * t506 + t632 * t660) * t550 + t490 * t584;
t421 = 0.1e1 / t423;
t414 = 0.1e1 / t416;
t411 = t580 * t550 * t445;
t409 = t581 * t445;
t407 = t573 * t460;
t405 = (t443 * t658 - t444 * t484) * t550 + t583 * t411;
t404 = -t409 * t583 + t443 * t464 + t444 * t481;
t402 = (t580 * t664 + (-0.2e1 * t484 * t617 + t659 * t467 + (t425 * t484 - t455 * t658 + t463 * t666) * t468) * t445) * t550;
t401 = t581 * t619 + (0.2e1 * t481 * t617 + t426 * t467 + (-t425 * t481 - t455 * t464 - t462 * t666) * t468) * t445;
t1 = [t460 * t602 + (t424 * t467 - t460 * t662) * t445, t401, t402, 0, 0, 0; t666 * t417 * t623 + (t425 * t417 + (-t403 * t666 + t407 * t424) * t418) * t414 + ((t407 * t621 + t573 * t646) * t414 + (t407 * t623 + ((-t408 * t445 * t641 + t619) * t644 + (t666 * t602 + t408 + (-t408 + t582) * t445) * t643) * t414) * t418) * t460 (-t404 * t647 - t417 * t465) * t622 + ((-t449 * t550 - t486 * t584) * t417 + t404 * t600 + (-t465 * t403 + t404 * t424 + (t401 * t666 - t409 * t425 + t462 + (t409 * t470 + t464) * t408) * t643 + (-t401 * t470 + t409 * t455 + t426 + (t409 * t666 - t481) * t408) * t644) * t418) * t414 (t417 * t477 * t550 - t405 * t647) * t622 + ((t583 * t402 + (-t408 * t420 + t425 * t444 - t443 * t455) * t411) * t647 + (t600 + t646) * t405 + (-t430 * t417 + (t477 * t403 + (t659 * t443 - t444 * t463 + (t443 * t484 + t444 * t658) * t408) * t460) * t418) * t550) * t414, 0, 0, 0; (-t435 * t438 + t439 * t645) * t620 + ((-t494 * t504 + t502 * t657 + t513 * t569 + t520 * t656 + t553 * t659 + t630 * t658) * t435 + t439 * t601 + (-t438 * t413 - (t494 * t501 - t503 * t657 - t514 * t656 + t519 * t569 + t556 * t659 - t631 * t658) * t440 - t439 * t412) * t436) * t421 (-t435 * t451 + t452 * t645) * t620 + ((-t449 * t520 + t450 * t553 - t482 * t502 + t483 * t630 + t486 * t614 + t504 * t615) * t435 + t452 * t601 + (-t451 * t413 - (t449 * t514 + t450 * t556 + t482 * t503 - t483 * t631 + t486 * t613 - t501 * t615) * t440 - t452 * t412) * t436) * t421 (-t435 * t453 + t454 * t645) * t620 + ((t428 * t553 - t430 * t520 + t475 * t630 - t477 * t502) * t435 + t454 * t601 + (-t453 * t413 - (t428 * t556 + t430 * t514 - t475 * t631 + t477 * t503) * t440 - t454 * t412) * t436) * t421 (-t435 * t567 - t663) * t620 + (0.2e1 * t618 + (-0.2e1 * t434 * t437 - t436 * t567 + t435) * t413) * t421, 0, 0;];
JaD_rot  = t1;
