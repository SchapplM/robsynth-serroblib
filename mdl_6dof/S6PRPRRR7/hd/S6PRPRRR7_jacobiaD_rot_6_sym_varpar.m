% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiaD_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:32
% EndTime: 2019-02-26 19:57:38
% DurationCPUTime: 6.08s
% Computational Cost: add. (40668->248), mult. (122035->475), div. (816->12), fcn. (158495->21), ass. (0->208)
t657 = sin(pkin(13));
t663 = cos(pkin(6));
t618 = t663 * t657;
t661 = cos(pkin(13));
t664 = sin(qJ(2));
t665 = cos(qJ(2));
t548 = -t618 * t664 + t661 * t665;
t550 = sin(pkin(14));
t587 = t618 * t665 + t661 * t664;
t659 = sin(pkin(6));
t612 = t657 * t659;
t658 = sin(pkin(7));
t598 = t658 * t612;
t660 = cos(pkin(14));
t662 = cos(pkin(7));
t522 = t548 * t660 + (-t587 * t662 + t598) * t550;
t555 = sin(qJ(4));
t558 = cos(qJ(4));
t551 = sin(pkin(8));
t552 = cos(pkin(8));
t584 = t587 * t660;
t571 = t548 * t550 + t584 * t662 - t598 * t660;
t574 = t587 * t658 + t612 * t662;
t670 = -t574 * t551 + t571 * t552;
t503 = t522 * t555 + t558 * t670;
t619 = t663 * t661;
t585 = -t619 * t665 + t657 * t664;
t586 = t619 * t664 + t657 * t665;
t614 = t659 * t658;
t599 = t661 * t614;
t521 = t586 * t660 + (-t585 * t662 - t599) * t550;
t581 = t585 * t660;
t569 = t550 * t586 + t581 * t662 + t599 * t660;
t616 = t662 * t659;
t573 = t585 * t658 - t616 * t661;
t565 = t551 * t573 - t552 * t569;
t501 = -t521 * t555 + t558 * t565;
t543 = t585 * qJD(2);
t544 = t586 * qJD(2);
t628 = t550 * t662;
t524 = -t543 * t660 - t544 * t628;
t617 = t662 * t660;
t592 = t543 * t550 - t544 * t617;
t626 = t658 * t551;
t579 = t544 * t626 + t552 * t592;
t478 = qJD(4) * t501 + t524 * t558 + t555 * t579;
t502 = t521 * t558 + t555 * t565;
t554 = sin(qJ(5));
t557 = cos(qJ(5));
t566 = t551 * t569 + t552 * t573;
t487 = t502 * t557 + t554 * t566;
t625 = t658 * t552;
t580 = t544 * t625 - t551 * t592;
t452 = qJD(5) * t487 + t478 * t554 - t557 * t580;
t485 = t502 * t554 - t557 * t566;
t483 = t485 ^ 2;
t615 = t660 * t659;
t601 = t664 * t615;
t604 = t665 * t616;
t613 = t658 * t663;
t537 = t601 + (t604 + t613) * t550;
t600 = t662 * t615;
t627 = t550 * t659;
t609 = t664 * t627;
t577 = -t600 * t665 + t609;
t536 = t613 * t660 - t577;
t603 = t665 * t614;
t547 = t662 * t663 - t603;
t610 = t536 * t552 + t547 * t551;
t512 = t537 * t558 + t555 * t610;
t520 = -t536 * t551 + t547 * t552;
t499 = t512 * t554 - t520 * t557;
t497 = 0.1e1 / t499 ^ 2;
t468 = t483 * t497 + 0.1e1;
t466 = 0.1e1 / t468;
t511 = -t537 * t555 + t558 * t610;
t542 = -t609 * t662 + t615 * t665;
t540 = t542 * qJD(2);
t541 = -t600 * t664 - t627 * t665;
t539 = t541 * qJD(2);
t602 = t664 * t614;
t593 = qJD(2) * t602;
t588 = t539 * t552 + t551 * t593;
t494 = qJD(4) * t511 + t540 * t558 + t555 * t588;
t500 = t512 * t557 + t520 * t554;
t533 = -t539 * t551 + t552 * t593;
t470 = qJD(5) * t500 + t494 * t554 - t533 * t557;
t496 = 0.1e1 / t499;
t648 = t485 * t497;
t435 = (-t452 * t496 + t470 * t648) * t466;
t469 = atan2(-t485, t499);
t464 = sin(t469);
t465 = cos(t469);
t611 = -t464 * t499 - t465 * t485;
t430 = t435 * t611 - t452 * t464 + t465 * t470;
t448 = -t464 * t485 + t465 * t499;
t445 = 0.1e1 / t448;
t446 = 0.1e1 / t448 ^ 2;
t672 = t430 * t445 * t446;
t504 = t522 * t558 - t555 * t670;
t568 = t551 * t571 + t552 * t574;
t488 = t504 * t554 - t557 * t568;
t624 = 0.2e1 * t488 * t672;
t605 = -t496 * t501 + t511 * t648;
t671 = t554 * t605;
t546 = t548 * qJD(2);
t545 = t587 * qJD(2);
t591 = t545 * t550 - t546 * t617;
t669 = t546 * t626 + t591 * t552;
t668 = t541 * t552 + t551 * t602;
t649 = t470 * t496 * t497;
t667 = -0.2e1 * (t452 * t648 - t483 * t649) / t468 ^ 2;
t532 = -t548 * t628 - t584;
t531 = -t548 * t617 + t550 * t587;
t596 = t531 * t552 + t548 * t626;
t666 = -t532 * t555 + t558 * t596;
t489 = t504 * t557 + t554 * t568;
t553 = sin(qJ(6));
t556 = cos(qJ(6));
t463 = t489 * t556 + t503 * t553;
t459 = 0.1e1 / t463;
t460 = 0.1e1 / t463 ^ 2;
t526 = -t545 * t660 - t546 * t628;
t480 = -t503 * qJD(4) + t526 * t558 + t669 * t555;
t578 = t546 * t625 - t551 * t591;
t455 = -qJD(5) * t488 + t480 * t557 + t554 * t578;
t479 = t504 * qJD(4) + t526 * t555 - t558 * t669;
t443 = qJD(6) * t463 + t455 * t553 - t479 * t556;
t462 = t489 * t553 - t503 * t556;
t458 = t462 ^ 2;
t451 = t458 * t460 + 0.1e1;
t652 = t460 * t462;
t634 = qJD(6) * t462;
t444 = t455 * t556 + t479 * t553 - t634;
t655 = t444 * t459 * t460;
t656 = (t443 * t652 - t458 * t655) / t451 ^ 2;
t654 = t446 * t488;
t454 = qJD(5) * t489 + t480 * t554 - t557 * t578;
t653 = t454 * t446;
t651 = t464 * t488;
t650 = t465 * t488;
t647 = t503 * t554;
t646 = t503 * t557;
t640 = t552 * t555;
t639 = t553 * t459;
t638 = t556 * t462;
t636 = qJD(5) * t554;
t635 = qJD(5) * t557;
t484 = t488 ^ 2;
t442 = t446 * t484 + 0.1e1;
t633 = 0.2e1 * (-t484 * t672 + t488 * t653) / t442 ^ 2;
t632 = -0.2e1 * t656;
t631 = 0.2e1 * t656;
t629 = t462 * t655;
t623 = 0.2e1 * t629;
t622 = -0.2e1 * t485 * t649;
t620 = qJD(6) * t646 + t480;
t507 = t532 * t558 + t555 * t596;
t517 = -t531 * t551 + t548 * t625;
t492 = t507 * t557 + t517 * t554;
t473 = t492 * t556 - t553 * t666;
t472 = t492 * t553 + t556 * t666;
t491 = t507 * t554 - t517 * t557;
t608 = t460 * t638 - t639;
t607 = -t487 * t496 + t500 * t648;
t583 = t586 * t662;
t530 = -t550 * t583 - t581;
t529 = t550 * t585 - t583 * t660;
t582 = t586 * t658;
t575 = t529 * t552 + t551 * t582;
t505 = t530 * t558 + t555 * t575;
t576 = -t529 * t551 + t552 * t582;
t490 = t505 * t554 - t557 * t576;
t518 = t542 * t558 + t555 * t668;
t534 = -t541 * t551 + t552 * t602;
t508 = t518 * t554 - t534 * t557;
t606 = -t490 * t496 + t508 * t648;
t527 = t545 * t617 + t546 * t550;
t597 = -t527 * t552 + t545 * t626;
t594 = qJD(2) * t603;
t590 = qJD(6) * t504 - t479 * t557 + t503 * t636;
t538 = t577 * qJD(2);
t528 = t545 * t628 - t546 * t660;
t525 = t543 * t617 + t544 * t550;
t515 = -t527 * t551 - t545 * t625;
t493 = -qJD(4) * t512 - t540 * t555 + t588 * t558;
t482 = qJD(4) * t666 + t528 * t558 - t597 * t555;
t481 = qJD(4) * t507 + t528 * t555 + t558 * t597;
t477 = -qJD(4) * t502 - t524 * t555 + t579 * t558;
t476 = (t538 * t640 + (-qJD(4) * t542 + t551 * t594) * t555 + ((-t550 * t604 - t601) * qJD(2) + t668 * qJD(4)) * t558) * t554 + t518 * t635 - (-t538 * t551 + t552 * t594) * t557 + t534 * t636;
t475 = t504 * t553 - t556 * t646;
t474 = -t504 * t556 - t553 * t646;
t471 = -qJD(5) * t499 + t494 * t557 + t533 * t554;
t457 = -qJD(5) * t491 + t482 * t557 + t515 * t554;
t456 = (qJD(5) * t505 + t525 * t551 + t543 * t625) * t557 + ((t543 * t628 - t544 * t660) * t558 + t525 * t640 - t543 * t555 * t626 + t576 * qJD(5) + (-t530 * t555 + t558 * t575) * qJD(4)) * t554;
t453 = -qJD(5) * t485 + t478 * t557 + t554 * t580;
t449 = 0.1e1 / t451;
t440 = 0.1e1 / t442;
t439 = t466 * t671;
t438 = t606 * t466;
t437 = t607 * t466;
t433 = (-t464 * t501 + t465 * t511) * t554 + t611 * t439;
t432 = t438 * t611 - t464 * t490 + t465 * t508;
t431 = t437 * t611 - t464 * t487 + t465 * t500;
t429 = t606 * t667 + (t508 * t622 - t456 * t496 + (t452 * t508 + t470 * t490 + t476 * t485) * t497) * t466;
t427 = t607 * t667 + (t500 * t622 - t453 * t496 + (t452 * t500 + t470 * t487 + t471 * t485) * t497) * t466;
t426 = t667 * t671 + (t605 * t635 + (t511 * t622 - t477 * t496 + (t452 * t511 + t470 * t501 + t485 * t493) * t497) * t554) * t466;
t1 = [0, t429, 0, t426, t427, 0; 0 (t432 * t654 - t445 * t491) * t633 + ((qJD(5) * t492 + t482 * t554 - t515 * t557) * t445 + t432 * t624 + (-t491 * t430 - t432 * t454 - (-t429 * t485 - t438 * t452 + t476 + (-t438 * t499 - t490) * t435) * t650 - (-t429 * t499 - t438 * t470 - t456 + (t438 * t485 - t508) * t435) * t651) * t446) * t440, 0 (t433 * t654 + t445 * t647) * t633 + ((-t479 * t554 - t503 * t635) * t445 + (-t653 + t624) * t433 + (t647 * t430 - (t511 * t635 - t426 * t485 - t439 * t452 + t493 * t554 + (-t439 * t499 - t501 * t554) * t435) * t650 - (-t501 * t635 - t426 * t499 - t439 * t470 - t477 * t554 + (t439 * t485 - t511 * t554) * t435) * t651) * t446) * t440 (t431 * t654 - t445 * t489) * t633 + (t431 * t624 + t455 * t445 + (-t489 * t430 - t431 * t454 - (-t427 * t485 - t437 * t452 + t471 + (-t437 * t499 - t487) * t435) * t650 - (-t427 * t499 - t437 * t470 - t453 + (t437 * t485 - t500) * t435) * t651) * t446) * t440, 0; 0 (-t459 * t472 + t473 * t652) * t631 + ((qJD(6) * t473 + t457 * t553 - t481 * t556) * t459 + t473 * t623 + (-t472 * t444 - (-qJD(6) * t472 + t457 * t556 + t481 * t553) * t462 - t473 * t443) * t460) * t449, 0 (-t459 * t474 + t475 * t652) * t631 + (t475 * t623 - t620 * t459 * t556 + t590 * t639 + (-t462 * t553 * t620 - t475 * t443 - t474 * t444 - t590 * t638) * t460) * t449, t608 * t488 * t632 + (t608 * t454 + ((-qJD(6) * t459 - 0.2e1 * t629) * t556 + (t443 * t556 + (t444 - t634) * t553) * t460) * t488) * t449, t632 + 0.2e1 * (t443 * t460 * t449 + (-t449 * t655 - t460 * t656) * t462) * t462;];
JaD_rot  = t1;
