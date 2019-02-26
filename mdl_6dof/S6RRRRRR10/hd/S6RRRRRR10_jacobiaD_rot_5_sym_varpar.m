% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR10_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiaD_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:49
% EndTime: 2019-02-26 22:52:56
% DurationCPUTime: 7.44s
% Computational Cost: add. (33277->290), mult. (101676->518), div. (962->12), fcn. (127596->19), ass. (0->217)
t695 = sin(qJ(2));
t696 = sin(qJ(1));
t554 = t696 * t695;
t563 = cos(qJ(2));
t694 = cos(pkin(6));
t641 = t694 * t695;
t621 = t696 * t641;
t698 = cos(qJ(1));
t642 = t694 * t698;
t539 = -qJD(1) * t621 - qJD(2) * t554 + (t698 * qJD(1) + qJD(2) * t642) * t563;
t560 = sin(qJ(3));
t562 = cos(qJ(3));
t619 = -t563 * t642 + t554;
t556 = sin(pkin(6));
t692 = sin(pkin(7));
t652 = t556 * t692;
t633 = t698 * t652;
t693 = cos(pkin(7));
t593 = t619 * t693 + t633;
t658 = t696 * t563;
t609 = t694 * t658 + t698 * t695;
t610 = t698 * t641 + t658;
t538 = t609 * qJD(1) + t610 * qJD(2);
t630 = t696 * t652;
t604 = qJD(1) * t630 - t538 * t693;
t493 = (t610 * qJD(3) - t604) * t560 + (t593 * qJD(3) - t539) * t562;
t559 = sin(qJ(4));
t529 = t593 * t560 - t610 * t562;
t571 = qJD(3) * t529 - t539 * t560 + t604 * t562;
t557 = cos(pkin(8));
t697 = cos(qJ(4));
t660 = t557 * t697;
t722 = t493 * t559 + t571 * t660;
t528 = t610 * t560 + t593 * t562;
t555 = sin(pkin(8));
t653 = t556 * t693;
t634 = t698 * t653;
t702 = t619 * t692 - t634;
t592 = t702 * t555;
t498 = t529 * t697 + (t528 * t557 - t592) * t559;
t721 = t529 * t559;
t717 = t493 * t697;
t576 = t528 * t697;
t591 = t697 * t592;
t494 = t557 * t576 - t591 - t721;
t648 = t563 * t693;
t616 = t560 * t648 + t695 * t562;
t637 = t694 * t692;
t543 = t616 * t556 + t560 * t637;
t548 = -t563 * t652 + t694 * t693;
t615 = -t695 * t560 + t562 * t648;
t542 = t615 * t556 + t562 * t637;
t659 = t697 * t542;
t661 = t555 * t697;
t602 = -t543 * t559 + t548 * t661 + t557 * t659;
t474 = atan2(-t494, -t602);
t465 = sin(t474);
t466 = cos(t474);
t448 = -t465 * t494 - t466 * t602;
t446 = 0.1e1 / t448 ^ 2;
t590 = -t609 * t693 + t630;
t608 = -t698 * t563 + t621;
t530 = t590 * t560 - t562 * t608;
t677 = t608 * t560;
t585 = -t590 * t562 - t677;
t581 = t585 * t697;
t631 = t696 * t653;
t589 = t609 * t692 + t631;
t588 = t589 * t555;
t587 = t697 * t588;
t499 = t530 * t559 + t557 * t581 - t587;
t489 = t499 ^ 2;
t444 = t446 * t489 + 0.1e1;
t537 = t610 * qJD(1) + t609 * qJD(2);
t536 = qJD(1) * t619 + t608 * qJD(2);
t620 = qJD(1) * t633;
t605 = t536 * t693 + t620;
t490 = qJD(3) * t677 + t605 * t560 + (t590 * qJD(3) - t537) * t562;
t607 = qJD(1) * t634 - t536 * t692;
t599 = t555 * t607;
t655 = qJD(4) * t697;
t643 = t530 * t655;
t672 = qJD(4) * t559;
t707 = -t530 * qJD(3) + t537 * t560;
t574 = t605 * t562 + t707;
t583 = t585 * t559;
t709 = -qJD(4) * t583 - t574 * t697;
t449 = t490 * t559 + t557 * t709 + t588 * t672 - t697 * t599 + t643;
t684 = t446 * t499;
t488 = t494 ^ 2;
t510 = 0.1e1 / t602 ^ 2;
t473 = t488 * t510 + 0.1e1;
t471 = 0.1e1 / t473;
t577 = t528 * t559;
t606 = qJD(1) * t631 + t538 * t692;
t598 = t555 * t606;
t451 = -t557 * qJD(4) * t577 - t529 * t655 + t592 * t672 - t697 * t598 - t722;
t513 = t543 * t697 + (t542 * t557 + t548 * t555) * t559;
t640 = t693 * t695;
t613 = -t560 * t563 - t562 * t640;
t629 = qJD(3) * t637;
t523 = -t560 * t629 + (t613 * qJD(2) - t616 * qJD(3)) * t556;
t614 = t560 * t640 - t562 * t563;
t524 = t562 * t629 + (-t614 * qJD(2) + t615 * qJD(3)) * t556;
t622 = t555 * t695 * t652;
t617 = t697 * t622;
t475 = -qJD(2) * t617 + t513 * qJD(4) - t523 * t660 + t524 * t559;
t509 = 0.1e1 / t602;
t680 = t494 * t510;
t628 = t451 * t509 + t475 * t680;
t435 = t628 * t471;
t636 = t465 * t602 - t466 * t494;
t429 = t636 * t435 - t451 * t465 + t466 * t475;
t445 = 0.1e1 / t448;
t690 = t429 * t445 * t446;
t670 = 0.2e1 * (t449 * t684 - t489 * t690) / t444 ^ 2;
t714 = t475 * t510;
t710 = -qJD(4) * t576 + t571 * t559;
t708 = -qJD(4) * t581 + t574 * t559;
t597 = t610 * t693;
t586 = t619 * t560 - t562 * t597;
t654 = t555 * t692;
t706 = t557 * t586 + t610 * t654;
t442 = 0.1e1 / t444;
t685 = t442 * t446;
t705 = t429 * t685 + t445 * t670;
t699 = 0.2e1 * t499;
t646 = t690 * t699;
t704 = t442 * t646 - t449 * t685 + t670 * t684;
t500 = t530 * t697 + (-t585 * t557 + t588) * t559;
t516 = t585 * t555 + t589 * t557;
t558 = sin(qJ(5));
t561 = cos(qJ(5));
t470 = t500 * t561 + t516 * t558;
t462 = 0.1e1 / t470;
t463 = 0.1e1 / t470 ^ 2;
t700 = -0.2e1 * t494;
t656 = t530 * t672;
t450 = qJD(4) * t587 + t490 * t697 + t557 * t708 + t559 * t599 - t656;
t649 = t562 * t693;
t477 = t607 * t557 + (-t536 * t649 - t562 * t620 - t707) * t555;
t440 = t470 * qJD(5) + t450 * t558 - t477 * t561;
t469 = t500 * t558 - t516 * t561;
t461 = t469 ^ 2;
t456 = t461 * t463 + 0.1e1;
t683 = t463 * t469;
t671 = qJD(5) * t469;
t441 = t450 * t561 + t477 * t558 - t671;
t687 = t441 * t462 * t463;
t689 = (t440 * t683 - t461 * t687) / t456 ^ 2;
t682 = t509 * t714;
t688 = (t451 * t680 + t488 * t682) / t473 ^ 2;
t686 = t442 * t445;
t681 = t494 * t509;
t675 = t555 * t558;
t674 = t555 * t561;
t673 = t557 * t559;
t669 = -0.2e1 * t689;
t668 = 0.2e1 * t689;
t667 = -0.2e1 * t688;
t665 = t509 * t688;
t663 = t442 * t684;
t662 = t469 * t687;
t657 = t529 * t672;
t651 = t557 * t692;
t650 = t560 * t693;
t645 = 0.2e1 * t662;
t644 = t682 * t700;
t515 = -t528 * t555 - t557 * t702;
t468 = t498 * t561 + t515 * t558;
t467 = t498 * t558 - t515 * t561;
t534 = t609 * t560 + t608 * t649;
t535 = -t609 * t562 + t608 * t650;
t505 = t535 * t697 + (t534 * t557 - t608 * t654) * t559;
t519 = -t534 * t555 - t608 * t651;
t482 = t505 * t561 + t519 * t558;
t481 = t505 * t558 - t519 * t561;
t635 = t697 * t654;
t627 = -t558 * t462 + t561 * t683;
t626 = -t498 * t509 + t513 * t680;
t533 = -t560 * t597 - t619 * t562;
t503 = t533 * t559 - t697 * t706;
t545 = t613 * t556;
t546 = t614 * t556;
t518 = -t545 * t660 - t546 * t559 - t617;
t625 = t503 * t509 + t518 * t680;
t506 = -t529 * t660 - t577;
t517 = t542 * t559 + t543 * t660;
t624 = t506 * t509 + t517 * t680;
t508 = -t530 * t673 - t581;
t623 = -t508 * t558 + t530 * t674;
t484 = t508 * t561 + t530 * t675;
t618 = -t465 + (-t466 * t681 + t465) * t471;
t603 = t528 * t660 - t661 * t702 - t721;
t596 = t534 * t660 - t535 * t559 - t608 * t635;
t502 = t534 * qJD(3) + t536 * t562 + t537 * t650;
t501 = -t535 * qJD(3) - t536 * t560 + t537 * t649;
t485 = -t501 * t555 - t537 * t651;
t480 = -t546 * t655 + (t545 * t557 + t622) * t672 + ((-t616 * qJD(2) + t613 * qJD(3)) * t559 - (-t615 * qJD(2) + t614 * qJD(3)) * t660 - qJD(2) * t563 * t635) * t556;
t479 = t524 * t660 + t523 * t559 + (-t543 * t673 + t659) * qJD(4);
t478 = t555 * t571 - t557 * t606;
t476 = t524 * t697 + (qJD(2) * t622 + t523 * t557) * t559 + t602 * qJD(4);
t460 = (t586 * qJD(3) - t538 * t562 - t539 * t650) * t559 + t533 * t655 - (-t533 * qJD(3) + t538 * t560 - t539 * t649) * t660 - t539 * t635 + t706 * t672;
t459 = t502 * t697 + (t501 * t557 - t537 * t654) * t559 + t596 * qJD(4);
t458 = -t493 * t660 + t557 * t657 + t710;
t457 = -t490 * t673 - t557 * t643 - t709;
t454 = 0.1e1 / t456;
t453 = t717 + (-t557 * t571 - t598) * t559 + t603 * qJD(4);
t452 = qJD(4) * t591 + t557 * t710 + t559 * t598 + t657 - t717;
t439 = t625 * t471;
t438 = t624 * t471;
t437 = t626 * t471;
t430 = t636 * t437 + t465 * t498 + t466 * t513;
t428 = t625 * t667 + (-t518 * t644 + t460 * t509 + (t451 * t518 + t475 * t503 + t480 * t494) * t510) * t471;
t427 = t624 * t667 + (-t517 * t644 + t458 * t509 + (t451 * t517 + t475 * t506 + t479 * t494) * t510) * t471;
t425 = t626 * t667 + (-t513 * t644 + t452 * t509 + (t451 * t513 - t475 * t498 + t476 * t494) * t510) * t471;
t1 = [-t665 * t699 + (t449 * t509 + t499 * t714) * t471, t428, t427, t425, 0, 0; (t498 * qJD(4) + t606 * t661 + t722) * t686 - (t618 * t449 + ((t435 * t471 * t681 + t667) * t465 + (-t665 * t700 - t435 + (t435 - t628) * t471) * t466) * t499) * t663 + t705 * t603 + t704 * t618 * t499 (t505 * qJD(4) - t501 * t660 + t502 * t559 + t537 * t635) * t686 - ((-t428 * t494 - t439 * t451 + t480 + (t439 * t602 - t503) * t435) * t466 + (t428 * t602 - t439 * t475 - t460 + (t439 * t494 - t518) * t435) * t465) * t663 + t705 * t596 + t704 * (t636 * t439 - t465 * t503 + t466 * t518) (t490 * t660 - t557 * t656 + t708) * t686 - ((-t427 * t494 - t438 * t451 + t479 + (t438 * t602 - t506) * t435) * t466 + (t427 * t602 - t438 * t475 - t458 + (t438 * t494 - t517) * t435) * t465) * t663 - t705 * (t530 * t660 - t583) + t704 * (t636 * t438 - t465 * t506 + t466 * t517) (t430 * t684 - t445 * t500) * t670 + (t430 * t646 + t450 * t445 + (-t500 * t429 - t430 * t449 + (-(-t425 * t494 - t437 * t451 + t476 + (t437 * t602 + t498) * t435) * t466 - (t425 * t602 - t437 * t475 - t452 + (t437 * t494 - t513) * t435) * t465) * t499) * t446) * t442, 0, 0; (-t462 * t467 + t468 * t683) * t668 + ((t468 * qJD(5) + t453 * t558 - t478 * t561) * t462 + t468 * t645 + (-t467 * t441 - (-t467 * qJD(5) + t453 * t561 + t478 * t558) * t469 - t468 * t440) * t463) * t454 (-t462 * t481 + t482 * t683) * t668 + ((t482 * qJD(5) + t459 * t558 - t485 * t561) * t462 + t482 * t645 + (-t481 * t441 - (-t481 * qJD(5) + t459 * t561 + t485 * t558) * t469 - t482 * t440) * t463) * t454 (t462 * t623 + t484 * t683) * t668 + ((t484 * qJD(5) + t457 * t558 - t490 * t674) * t462 + t484 * t645 + (t623 * t441 - (t623 * qJD(5) + t457 * t561 + t490 * t675) * t469 - t484 * t440) * t463) * t454, t627 * t499 * t669 + (t627 * t449 + ((-qJD(5) * t462 - 0.2e1 * t662) * t561 + (t440 * t561 + (t441 - t671) * t558) * t463) * t499) * t454, t669 + 0.2e1 * (t440 * t454 * t463 + (-t454 * t687 - t463 * t689) * t469) * t469, 0;];
JaD_rot  = t1;
