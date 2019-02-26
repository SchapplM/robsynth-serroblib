% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR12_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobiaD_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:11
% EndTime: 2019-02-26 21:21:20
% DurationCPUTime: 8.14s
% Computational Cost: add. (54698->282), mult. (164092->508), div. (983->12), fcn. (210779->21), ass. (0->223)
t728 = sin(pkin(14));
t735 = cos(pkin(6));
t679 = t735 * t728;
t732 = cos(pkin(14));
t736 = sin(qJ(1));
t738 = cos(qJ(1));
t646 = -t736 * t679 + t738 * t732;
t586 = t646 * qJD(1);
t597 = sin(qJ(3));
t730 = sin(pkin(7));
t731 = sin(pkin(6));
t677 = t731 * t730;
t661 = t736 * t677;
t737 = cos(qJ(3));
t651 = t737 * t661;
t665 = t738 * t677;
t681 = t735 * t732;
t645 = t736 * t681 + t738 * t728;
t585 = t645 * qJD(1);
t734 = cos(pkin(7));
t694 = t585 * t734;
t746 = -t738 * t681 + t736 * t728;
t748 = t746 * t734;
t644 = t738 * t679 + t736 * t732;
t749 = t644 * t737;
t616 = qJD(1) * t651 + (-t749 + (t748 + t665) * t597) * qJD(3) - t586 * t597 - t737 * t694;
t678 = t734 * t731;
t662 = t736 * t678;
t647 = qJD(1) * t662 + t585 * t730;
t729 = sin(pkin(8));
t733 = cos(pkin(8));
t541 = t616 * t729 - t647 * t733;
t595 = sin(qJ(5));
t767 = t541 * t595;
t599 = cos(qJ(5));
t766 = t541 * t599;
t596 = sin(qJ(4));
t611 = t616 * t733 + t647 * t729;
t765 = t611 * t596;
t600 = cos(qJ(4));
t764 = t611 * t600;
t747 = t597 * t665 - t749;
t573 = t597 * t748 + t747;
t763 = t573 * t596;
t762 = t573 * t600;
t628 = -t738 * t678 + t730 * t746;
t761 = t628 * t729;
t760 = t628 * t733;
t632 = t645 * t734;
t620 = t646 * t597 + t737 * t632 - t651;
t627 = t645 * t730 + t662;
t759 = -t620 * t733 + t627 * t729;
t758 = t732 * t678 + t735 * t730;
t639 = t737 * t748;
t574 = t646 * t737 + (-t632 + t661) * t597;
t676 = t731 * t728;
t581 = t758 * t597 + t737 * t676;
t580 = -t597 * t676 + t758 * t737;
t587 = -t732 * t677 + t735 * t734;
t655 = t733 * t580 + t729 * t587;
t563 = -t581 * t596 + t655 * t600;
t577 = t580 * qJD(3);
t578 = t581 * qJD(3);
t692 = t596 * t733;
t538 = t563 * qJD(4) + t577 * t600 - t578 * t692;
t564 = t581 * t600 + t655 * t596;
t570 = -t580 * t729 + t587 * t733;
t546 = t564 * t599 + t570 * t595;
t691 = t599 * t729;
t516 = t546 * qJD(5) + t538 * t595 - t578 * t691;
t545 = t564 * t595 - t570 * t599;
t543 = 0.1e1 / t545 ^ 2;
t752 = t516 * t543;
t542 = 0.1e1 / t545;
t652 = t737 * t665;
t750 = t644 * t597;
t742 = -t750 - t639;
t622 = t652 - t742;
t615 = -t622 * t733 + t761;
t548 = t615 * t600 + t763;
t549 = t615 * t596 - t762;
t614 = t622 * t729 + t760;
t527 = t549 * t595 - t614 * t599;
t714 = t527 * t543;
t669 = -t542 * t548 + t563 * t714;
t751 = t595 * t669;
t631 = qJD(1) * t748;
t741 = -t750 - t652;
t609 = t741 * qJD(1) + qJD(3) * t574 - t737 * t631;
t626 = qJD(1) * t628;
t745 = t609 * t733 + t729 * t626;
t705 = qJD(4) * t574;
t744 = -t733 * t705 - t609;
t743 = t759 * t600;
t572 = t639 - t741;
t656 = t733 * t572 - t761;
t739 = t656 * t600 - t763;
t515 = atan2(-t527, t545);
t506 = sin(t515);
t507 = cos(t515);
t485 = -t506 * t527 + t507 * t545;
t482 = 0.1e1 / t485;
t553 = t574 * t600 + t759 * t596;
t613 = t620 * t729 + t627 * t733;
t533 = t553 * t599 + t613 * t595;
t552 = t574 * t596 - t743;
t594 = sin(qJ(6));
t598 = cos(qJ(6));
t511 = t533 * t598 + t552 * t594;
t503 = 0.1e1 / t511;
t483 = 0.1e1 / t485 ^ 2;
t504 = 0.1e1 / t511 ^ 2;
t612 = t613 * t599;
t532 = t553 * t595 - t612;
t526 = t532 ^ 2;
t481 = t483 * t526 + 0.1e1;
t555 = t747 * qJD(1) - t620 * qJD(3) + t597 * t631;
t497 = t555 * t600 + t743 * qJD(4) + (-t705 - t745) * t596;
t607 = t609 * t729 - t733 * t626;
t486 = t533 * qJD(5) + t497 * t595 - t607 * t599;
t720 = t486 * t483;
t525 = t527 ^ 2;
t514 = t525 * t543 + 0.1e1;
t512 = 0.1e1 / t514;
t642 = qJD(1) * t661 - t694;
t683 = qJD(3) * t652 - t586 * t737;
t556 = t742 * qJD(3) + t642 * t597 - t683;
t499 = t548 * qJD(4) + t556 * t600 + t765;
t529 = t549 * t599 + t614 * t595;
t488 = t529 * qJD(5) + t499 * t595 + t766;
t673 = -t488 * t542 + t516 * t714;
t472 = t673 * t512;
t675 = -t506 * t545 - t507 * t527;
t466 = t675 * t472 - t488 * t506 + t507 * t516;
t484 = t482 * t483;
t726 = t466 * t484;
t727 = 0.2e1 * (-t526 * t726 + t532 * t720) / t481 ^ 2;
t704 = qJD(5) * t595;
t487 = qJD(5) * t612 + t497 * t599 - t553 * t704 + t607 * t595;
t496 = t553 * qJD(4) + t555 * t596 + t600 * t745;
t475 = t511 * qJD(6) + t487 * t594 - t496 * t598;
t510 = t533 * t594 - t552 * t598;
t502 = t510 ^ 2;
t493 = t502 * t504 + 0.1e1;
t719 = t504 * t510;
t702 = qJD(6) * t510;
t476 = t487 * t598 + t496 * t594 - t702;
t722 = t476 * t503 * t504;
t724 = 0.2e1 * (t475 * t719 - t502 * t722) / t493 ^ 2;
t716 = t542 * t752;
t723 = 0.2e1 * (t488 * t714 - t525 * t716) / t514 ^ 2;
t721 = t483 * t532;
t718 = t506 * t532;
t717 = t507 * t532;
t715 = t527 * t542;
t713 = t552 * t595;
t712 = t552 * t599;
t708 = t594 * t503;
t707 = t598 * t510;
t703 = qJD(5) * t599;
t701 = 0.2e1 * t484 * t532;
t700 = t482 * t727;
t699 = t483 * t727;
t698 = t542 * t723;
t697 = t483 * t718;
t696 = t483 * t717;
t695 = t510 * t722;
t693 = t595 * t729;
t690 = t600 * t733;
t689 = t729 * t555;
t688 = t466 * t701;
t687 = 0.2e1 * t695;
t686 = -0.2e1 * t527 * t716;
t682 = qJD(6) * t712 + t497;
t551 = t656 * t596 + t762;
t567 = -t572 * t729 - t760;
t531 = t551 * t599 + t567 * t595;
t509 = t531 * t598 - t594 * t739;
t508 = t531 * t594 + t598 * t739;
t561 = -t574 * t692 - t620 * t600;
t536 = t561 * t599 + t574 * t693;
t560 = t574 * t690 - t620 * t596;
t524 = t536 * t598 + t560 * t594;
t523 = t536 * t594 - t560 * t598;
t530 = t551 * t595 - t567 * t599;
t672 = t504 * t707 - t708;
t671 = -t529 * t542 + t546 * t714;
t559 = t573 * t692 - t622 * t600;
t534 = t559 * t595 + t573 * t691;
t568 = t580 * t600 - t581 * t692;
t562 = t568 * t595 - t581 * t691;
t670 = -t534 * t542 + t562 * t714;
t659 = -t561 * t595 + t574 * t691;
t658 = -t506 + (t507 * t715 + t506) * t512;
t653 = qJD(6) * t553 - t496 * t599 + t552 * t704;
t619 = qJD(4) * t620;
t558 = qJD(3) * t639 + (t644 * qJD(3) - t642) * t597 + t683;
t537 = -t564 * qJD(4) - t577 * t596 - t578 * t690;
t522 = (-t577 * t692 - t578 * t600 + (-t580 * t596 - t581 * t690) * qJD(4)) * t595 + t568 * t703 - t577 * t691 + t581 * qJD(5) * t693;
t521 = t553 * t594 - t598 * t712;
t520 = -t553 * t598 - t594 * t712;
t519 = -t555 * t692 + t596 * t619 + t600 * t744;
t518 = t555 * t690 + t596 * t744 - t600 * t619;
t517 = -t545 * qJD(5) + t538 * t599 + t578 * t693;
t501 = qJD(4) * t739 + t558 * t600 - t765;
t500 = t551 * qJD(4) + t558 * t596 + t764;
t498 = -qJD(4) * t549 - t556 * t596 + t764;
t495 = (-t556 * t692 + (t573 * t690 + t622 * t596) * qJD(4) + t616 * t600) * t595 - t556 * t691 + (t559 * t599 - t573 * t693) * qJD(5);
t494 = t659 * qJD(5) + t519 * t599 + t595 * t689;
t491 = 0.1e1 / t493;
t490 = -t530 * qJD(5) + t501 * t599 + t767;
t489 = -t527 * qJD(5) + t499 * t599 - t767;
t479 = 0.1e1 / t481;
t478 = t512 * t751;
t477 = t670 * t512;
t474 = t671 * t512;
t471 = t658 * t532;
t469 = (-t506 * t548 + t507 * t563) * t595 + t675 * t478;
t467 = t675 * t474 - t506 * t529 + t507 * t546;
t465 = -t670 * t723 + (t562 * t686 - t495 * t542 + (t488 * t562 + t516 * t534 + t522 * t527) * t543) * t512;
t463 = -t671 * t723 + (t546 * t686 - t489 * t542 + (t488 * t546 + t516 * t529 + t517 * t527) * t543) * t512;
t462 = -t723 * t751 + (t669 * t703 + (t563 * t686 - t498 * t542 + (t488 * t563 + t516 * t548 + t527 * t537) * t543) * t595) * t512;
t1 = [t532 * t698 + (-t486 * t542 + t532 * t752) * t512, 0, t465, t462, t463, 0; -t530 * t700 + ((t531 * qJD(5) + t501 * t595 - t766) * t482 + (-t530 * t466 - t471 * t486) * t483) * t479 + (t471 * t699 + (0.2e1 * t471 * t726 - (-t472 * t512 * t715 - t723) * t697 - (-t527 * t698 - t472 + (t472 - t673) * t512) * t696 - t658 * t720) * t479) * t532, 0, t659 * t700 + ((t536 * qJD(5) + t519 * t595 - t599 * t689) * t482 + t659 * t483 * t466 - ((-t465 * t527 - t477 * t488 + t522 + (-t477 * t545 - t534) * t472) * t507 + (-t465 * t545 - t477 * t516 - t495 + (t477 * t527 - t562) * t472) * t506) * t721) * t479 + (t532 * t699 + (-t720 + t688) * t479) * (t675 * t477 - t506 * t534 + t507 * t562) (t469 * t721 + t482 * t713) * t727 + (-t469 * t720 + (-t496 * t595 - t552 * t703) * t482 + (t469 * t701 + t483 * t713) * t466 - (t563 * t703 - t462 * t527 - t478 * t488 + t537 * t595 + (-t478 * t545 - t548 * t595) * t472) * t696 - (-t548 * t703 - t462 * t545 - t478 * t516 - t498 * t595 + (t478 * t527 - t563 * t595) * t472) * t697) * t479 (t467 * t721 - t482 * t533) * t727 + (t467 * t688 + t487 * t482 + (-t533 * t466 - t467 * t486 - (-t463 * t527 - t474 * t488 + t517 + (-t474 * t545 - t529) * t472) * t717 - (-t463 * t545 - t474 * t516 - t489 + (t474 * t527 - t546) * t472) * t718) * t483) * t479, 0; (-t503 * t508 + t509 * t719) * t724 + ((t509 * qJD(6) + t490 * t594 - t500 * t598) * t503 + t509 * t687 + (-t508 * t476 - (-t508 * qJD(6) + t490 * t598 + t500 * t594) * t510 - t509 * t475) * t504) * t491, 0 (-t503 * t523 + t524 * t719) * t724 + ((t524 * qJD(6) + t494 * t594 - t518 * t598) * t503 + t524 * t687 + (-t523 * t476 - (-t523 * qJD(6) + t494 * t598 + t518 * t594) * t510 - t524 * t475) * t504) * t491 (-t503 * t520 + t521 * t719) * t724 + (t521 * t687 - t682 * t503 * t598 + t653 * t708 + (-t682 * t510 * t594 - t521 * t475 - t520 * t476 - t653 * t707) * t504) * t491, -t672 * t532 * t724 + (t672 * t486 + ((-qJD(6) * t503 - 0.2e1 * t695) * t598 + (t475 * t598 + (t476 - t702) * t594) * t504) * t532) * t491, -t724 + (0.2e1 * t475 * t504 * t491 + (-0.2e1 * t491 * t722 - t504 * t724) * t510) * t510;];
JaD_rot  = t1;
