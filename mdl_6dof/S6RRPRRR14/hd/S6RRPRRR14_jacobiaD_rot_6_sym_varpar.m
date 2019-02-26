% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR14_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:49
% EndTime: 2019-02-26 22:55:58
% DurationCPUTime: 9.07s
% Computational Cost: add. (54756->300), mult. (163156->566), div. (983->12), fcn. (208653->21), ass. (0->256)
t645 = sin(pkin(14));
t792 = cos(pkin(6));
t793 = sin(qJ(2));
t736 = t792 * t793;
t794 = sin(qJ(1));
t795 = cos(qJ(2));
t796 = cos(qJ(1));
t691 = t796 * t736 + t794 * t795;
t790 = cos(pkin(14));
t684 = t691 * t790;
t737 = t792 * t795;
t700 = -t796 * t737 + t794 * t793;
t791 = cos(pkin(7));
t693 = t700 * t791;
t788 = sin(pkin(7));
t789 = sin(pkin(6));
t731 = t789 * t788;
t714 = t796 * t731;
t615 = (t693 + t714) * t645 - t684;
t650 = sin(qJ(4));
t818 = t615 * t650;
t653 = cos(qJ(4));
t817 = t615 * t653;
t649 = sin(qJ(5));
t692 = t794 * t737 + t796 * t793;
t625 = qJD(1) * t692 + qJD(2) * t691;
t724 = -t794 * t736 + t796 * t795;
t626 = qJD(1) * t724 - qJD(2) * t700;
t730 = t788 * t790;
t706 = t789 * t730;
t695 = t794 * t706;
t734 = t791 * t790;
t599 = -qJD(1) * t695 + t625 * t734 + t626 * t645;
t646 = sin(pkin(8));
t647 = cos(pkin(8));
t733 = t791 * t789;
t709 = t794 * t733;
t689 = qJD(1) * t709 + t625 * t788;
t672 = t599 * t646 + t647 * t689;
t816 = t649 * t672;
t652 = cos(qJ(5));
t815 = t652 * t672;
t708 = t794 * t731;
t598 = t645 * (qJD(1) * t708 - t625 * t791) + t626 * t790;
t671 = -t599 * t647 + t646 * t689;
t814 = t598 * t653 + t650 * t671;
t813 = -t598 * t650 + t671 * t653;
t715 = t796 * t733;
t679 = t700 * t788 - t715;
t810 = t646 * t679;
t809 = t647 * t679;
t685 = t692 * t790;
t673 = t645 * t724 + t685 * t791 - t695;
t678 = t692 * t788 + t709;
t808 = t678 * t646 - t673 * t647;
t686 = t691 * t645;
t696 = t796 * t706;
t807 = t686 + t696;
t732 = t790 * t789;
t710 = t793 * t732;
t713 = t795 * t733;
t629 = t710 + (t788 * t792 + t713) * t645;
t707 = t791 * t732;
t746 = t645 * t789;
t725 = t793 * t746;
t682 = -t707 * t795 + t725;
t628 = t730 * t792 - t682;
t712 = t795 * t731;
t636 = t791 * t792 - t712;
t726 = t628 * t647 + t636 * t646;
t595 = -t629 * t650 + t653 * t726;
t635 = -t725 * t791 + t732 * t795;
t633 = t635 * qJD(2);
t634 = -t707 * t793 - t746 * t795;
t632 = t634 * qJD(2);
t711 = t793 * t731;
t701 = qJD(2) * t711;
t690 = t632 * t647 + t646 * t701;
t573 = qJD(4) * t595 + t633 * t653 + t650 * t690;
t596 = t629 * t653 + t650 * t726;
t612 = -t628 * t646 + t636 * t647;
t578 = t596 * t652 + t612 * t649;
t618 = -t632 * t646 + t647 * t701;
t551 = qJD(5) * t578 + t573 * t649 - t618 * t652;
t577 = t596 * t649 - t612 * t652;
t575 = 0.1e1 / t577 ^ 2;
t805 = t551 * t575;
t574 = 0.1e1 / t577;
t694 = t700 * t790;
t669 = t694 * t791 + t807;
t662 = -t647 * t669 + t810;
t583 = t653 * t662 + t818;
t584 = t650 * t662 - t817;
t663 = t646 * t669 + t809;
t562 = t584 * t649 - t652 * t663;
t774 = t562 * t575;
t719 = -t574 * t583 + t595 * t774;
t804 = t649 * t719;
t624 = qJD(1) * t691 + qJD(2) * t692;
t675 = qJD(1) * t700 - t724 * qJD(2);
t674 = t675 * t790;
t665 = qJD(1) * t696 + t624 * t645 + t674 * t791;
t667 = qJD(1) * t715 - t675 * t788;
t803 = t646 * t667 + t647 * t665;
t802 = t808 * t653;
t801 = t634 * t647 + t646 * t711;
t699 = t724 * t791;
t622 = -t645 * t699 - t685;
t621 = t645 * t692 - t699 * t790;
t698 = t724 * t788;
t687 = t621 * t647 + t646 * t698;
t800 = -t622 * t650 + t653 * t687;
t614 = t693 * t790 + t807;
t727 = t614 * t647 - t810;
t799 = t653 * t727 - t818;
t546 = atan2(-t562, t577);
t537 = sin(t546);
t538 = cos(t546);
t520 = -t537 * t562 + t538 * t577;
t517 = 0.1e1 / t520;
t616 = t724 * t790 + (-t692 * t791 + t708) * t645;
t588 = t616 * t653 + t808 * t650;
t661 = t646 * t673 + t647 * t678;
t568 = t588 * t652 + t649 * t661;
t587 = t616 * t650 - t802;
t648 = sin(qJ(6));
t651 = cos(qJ(6));
t550 = t568 * t651 + t587 * t648;
t543 = 0.1e1 / t550;
t518 = 0.1e1 / t520 ^ 2;
t544 = 0.1e1 / t550 ^ 2;
t798 = -0.2e1 * t562;
t567 = t588 * t649 - t652 * t661;
t797 = 0.2e1 * t567;
t561 = t567 ^ 2;
t516 = t518 * t561 + 0.1e1;
t597 = -t624 * t790 + (qJD(1) * t714 + t675 * t791) * t645;
t762 = qJD(4) * t650;
t532 = qJD(4) * t802 + t597 * t653 - t616 * t762 + t650 * t803;
t660 = -t646 * t665 + t647 * t667;
t521 = qJD(5) * t568 + t532 * t649 - t652 * t660;
t781 = t521 * t518;
t560 = t562 ^ 2;
t541 = t560 * t575 + 0.1e1;
t539 = 0.1e1 / t541;
t534 = qJD(4) * t583 + t814;
t564 = t584 * t652 + t649 * t663;
t523 = qJD(5) * t564 + t534 * t649 - t815;
t723 = -t523 * t574 + t551 * t774;
t507 = t723 * t539;
t729 = -t537 * t577 - t538 * t562;
t501 = t507 * t729 - t523 * t537 + t538 * t551;
t519 = t517 * t518;
t786 = t501 * t519;
t787 = (-t561 * t786 + t567 * t781) / t516 ^ 2;
t522 = -qJD(5) * t567 + t532 * t652 + t649 * t660;
t531 = qJD(4) * t588 + t597 * t650 - t653 * t803;
t511 = qJD(6) * t550 + t522 * t648 - t531 * t651;
t549 = t568 * t648 - t587 * t651;
t542 = t549 ^ 2;
t528 = t542 * t544 + 0.1e1;
t778 = t544 * t549;
t758 = qJD(6) * t549;
t512 = t522 * t651 + t531 * t648 - t758;
t783 = t512 * t543 * t544;
t785 = (t511 * t778 - t542 * t783) / t528 ^ 2;
t776 = t574 * t805;
t784 = (t523 * t774 - t560 * t776) / t541 ^ 2;
t782 = t518 * t567;
t780 = t537 * t567;
t779 = t538 * t567;
t777 = t549 * t651;
t775 = t562 * t574;
t773 = t587 * t649;
t772 = t587 * t652;
t764 = t647 * t650;
t763 = t648 * t543;
t760 = qJD(5) * t649;
t759 = qJD(5) * t652;
t757 = 0.2e1 * t787;
t756 = -0.2e1 * t785;
t755 = 0.2e1 * t785;
t754 = -0.2e1 * t784;
t753 = t519 * t797;
t752 = t574 * t784;
t751 = t518 * t780;
t750 = t518 * t779;
t749 = t549 * t783;
t747 = t645 * t791;
t745 = t646 * t788;
t744 = t647 * t788;
t743 = t501 * t753;
t742 = 0.2e1 * t749;
t741 = t776 * t798;
t735 = qJD(6) * t772 + t532;
t586 = t650 * t727 + t817;
t603 = -t614 * t646 - t809;
t566 = t586 * t652 + t603 * t649;
t548 = t566 * t651 - t648 * t799;
t547 = t566 * t648 + t651 * t799;
t592 = t622 * t653 + t650 * t687;
t609 = -t621 * t646 + t647 * t698;
t571 = t592 * t652 + t609 * t649;
t556 = t571 * t651 - t648 * t800;
t555 = t571 * t648 + t651 * t800;
t565 = t586 * t649 - t603 * t652;
t570 = t592 * t649 - t609 * t652;
t722 = t544 * t777 - t763;
t721 = -t564 * t574 + t578 * t774;
t620 = -t686 * t791 - t694;
t619 = t645 * t700 - t684 * t791;
t683 = t691 * t788;
t680 = t619 * t647 + t646 * t683;
t590 = t620 * t653 + t650 * t680;
t681 = -t619 * t646 + t647 * t683;
t569 = t590 * t649 - t652 * t681;
t608 = t635 * t653 + t650 * t801;
t623 = -t634 * t646 + t647 * t711;
t593 = t608 * t649 - t623 * t652;
t720 = -t569 * t574 + t593 * t774;
t604 = t624 * t734 - t645 * t675;
t705 = -t604 * t647 + t624 * t745;
t704 = -t537 + (t538 * t775 + t537) * t539;
t702 = qJD(2) * t712;
t697 = qJD(6) * t588 - t531 * t652 + t587 * t760;
t631 = t682 * qJD(2);
t606 = t625 * t645 - t626 * t734;
t605 = t624 * t747 + t674;
t589 = -t604 * t646 - t624 * t744;
t572 = -qJD(4) * t596 - t633 * t650 + t690 * t653;
t559 = (t646 * t650 * t702 + t631 * t764 - t635 * t762 + ((-t645 * t713 - t710) * qJD(2) + t801 * qJD(4)) * t653) * t649 + t608 * t759 - (-t631 * t646 + t647 * t702) * t652 + t623 * t760;
t558 = t588 * t648 - t651 * t772;
t557 = -t588 * t651 - t648 * t772;
t554 = qJD(4) * t800 + t605 * t653 - t705 * t650;
t553 = qJD(4) * t592 + t605 * t650 + t653 * t705;
t552 = -qJD(5) * t577 + t573 * t652 + t618 * t649;
t536 = qJD(4) * t799 - t814;
t535 = qJD(4) * t586 + t813;
t533 = -qJD(4) * t584 + t813;
t530 = (qJD(5) * t590 + t606 * t646 - t626 * t744) * t652 + ((-t625 * t790 - t626 * t747) * t653 + t606 * t764 + t626 * t650 * t745 + t681 * qJD(5) + (-t620 * t650 + t653 * t680) * qJD(4)) * t649;
t529 = -qJD(5) * t570 + t554 * t652 + t589 * t649;
t526 = 0.1e1 / t528;
t525 = -qJD(5) * t565 + t536 * t652 - t816;
t524 = -qJD(5) * t562 + t534 * t652 + t816;
t514 = 0.1e1 / t516;
t513 = t539 * t804;
t510 = t720 * t539;
t509 = t721 * t539;
t506 = t704 * t567;
t504 = (-t537 * t583 + t538 * t595) * t649 + t729 * t513;
t503 = t510 * t729 - t537 * t569 + t538 * t593;
t502 = t509 * t729 - t537 * t564 + t538 * t578;
t500 = t720 * t754 + (t593 * t741 - t530 * t574 + (t523 * t593 + t551 * t569 + t559 * t562) * t575) * t539;
t498 = t721 * t754 + (t578 * t741 - t524 * t574 + (t523 * t578 + t551 * t564 + t552 * t562) * t575) * t539;
t497 = t754 * t804 + (t719 * t759 + (t595 * t741 - t533 * t574 + (t523 * t595 + t551 * t583 + t562 * t572) * t575) * t649) * t539;
t1 = [t752 * t797 + (-t521 * t574 + t567 * t805) * t539, t500, 0, t497, t498, 0; -0.2e1 * t565 * t517 * t787 + ((qJD(5) * t566 + t536 * t649 + t815) * t517 + (-t565 * t501 - t506 * t521) * t518) * t514 + (t506 * t518 * t757 + (0.2e1 * t506 * t786 - (-t507 * t539 * t775 + t754) * t751 - (t752 * t798 - t507 + (t507 - t723) * t539) * t750 - t704 * t781) * t514) * t567 (t503 * t782 - t517 * t570) * t757 + ((qJD(5) * t571 + t554 * t649 - t589 * t652) * t517 + t503 * t743 + (-t570 * t501 - t503 * t521 - (-t500 * t562 - t510 * t523 + t559 + (-t510 * t577 - t569) * t507) * t779 - (-t500 * t577 - t510 * t551 - t530 + (t510 * t562 - t593) * t507) * t780) * t518) * t514, 0 (t504 * t782 + t517 * t773) * t757 + (-t504 * t781 + (-t531 * t649 - t587 * t759) * t517 + (t504 * t753 + t518 * t773) * t501 - (t595 * t759 - t497 * t562 - t513 * t523 + t572 * t649 + (-t513 * t577 - t583 * t649) * t507) * t750 - (-t583 * t759 - t497 * t577 - t513 * t551 - t533 * t649 + (t513 * t562 - t595 * t649) * t507) * t751) * t514 (t502 * t782 - t517 * t568) * t757 + (t502 * t743 + t522 * t517 + (-t568 * t501 - t502 * t521 - (-t498 * t562 - t509 * t523 + t552 + (-t509 * t577 - t564) * t507) * t779 - (-t498 * t577 - t509 * t551 - t524 + (t509 * t562 - t578) * t507) * t780) * t518) * t514, 0; (-t543 * t547 + t548 * t778) * t755 + ((qJD(6) * t548 + t525 * t648 - t535 * t651) * t543 + t548 * t742 + (-t547 * t512 - (-qJD(6) * t547 + t525 * t651 + t535 * t648) * t549 - t548 * t511) * t544) * t526 (-t543 * t555 + t556 * t778) * t755 + ((qJD(6) * t556 + t529 * t648 - t553 * t651) * t543 + t556 * t742 + (-t555 * t512 - (-qJD(6) * t555 + t529 * t651 + t553 * t648) * t549 - t556 * t511) * t544) * t526, 0 (-t543 * t557 + t558 * t778) * t755 + (t558 * t742 - t735 * t543 * t651 + t697 * t763 + (-t549 * t648 * t735 - t558 * t511 - t557 * t512 - t697 * t777) * t544) * t526, t722 * t567 * t756 + (t722 * t521 + ((-qJD(6) * t543 - 0.2e1 * t749) * t651 + (t511 * t651 + (t512 - t758) * t648) * t544) * t567) * t526, t756 + 0.2e1 * (t511 * t544 * t526 + (-t526 * t783 - t544 * t785) * t549) * t549;];
JaD_rot  = t1;
