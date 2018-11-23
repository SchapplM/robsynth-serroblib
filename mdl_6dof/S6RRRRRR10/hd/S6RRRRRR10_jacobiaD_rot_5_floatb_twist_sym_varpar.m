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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR10_jacobiaD_rot_5_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_rot_5_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_jacobiaD_rot_5_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiaD_rot_5_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:20
% EndTime: 2018-11-23 11:27:30
% DurationCPUTime: 10.12s
% Computational Cost: add. (135037->348), mult. (138954->564), div. (962->12), fcn. (127596->31), ass. (0->245)
t822 = pkin(7) + qJ(3);
t803 = sin(t822);
t783 = t803 / 0.2e1;
t823 = pkin(7) - qJ(3);
t804 = sin(t823);
t685 = t783 - t804 / 0.2e1;
t850 = sin(qJ(2));
t851 = sin(qJ(1));
t710 = t851 * t850;
t824 = pkin(6) + qJ(2);
t791 = cos(t824) / 0.2e1;
t825 = pkin(6) - qJ(2);
t809 = cos(t825);
t766 = t809 / 0.2e1 + t791;
t853 = cos(qJ(1));
t751 = -t853 * t766 + t710;
t789 = cos(t822) / 0.2e1;
t808 = cos(t823);
t688 = t789 - t808 / 0.2e1;
t847 = sin(pkin(6));
t794 = t853 * t847;
t777 = t688 * t794;
t718 = cos(qJ(3));
t805 = sin(t824);
t785 = t805 / 0.2e1;
t806 = sin(t825);
t786 = -t806 / 0.2e1;
t763 = t785 + t786;
t852 = cos(qJ(2));
t795 = t851 * t852;
t743 = t763 * t853 + t795;
t861 = t743 * t718;
t642 = -t685 * t751 + t777 + t861;
t715 = sin(qJ(4));
t820 = pkin(8) + qJ(4);
t801 = sin(t820);
t781 = t801 / 0.2e1;
t821 = pkin(8) - qJ(4);
t802 = sin(t821);
t782 = t802 / 0.2e1;
t761 = t781 + t782;
t787 = cos(t820) / 0.2e1;
t807 = cos(t821);
t764 = t807 / 0.2e1 + t787;
t848 = cos(pkin(7));
t779 = t848 * t847;
t770 = t853 * t779;
t846 = sin(pkin(7));
t855 = -t751 * t846 + t770;
t784 = t804 / 0.2e1;
t762 = t783 + t784;
t755 = t762 * t847;
t750 = t853 * t755;
t765 = t808 / 0.2e1 + t789;
t849 = sin(qJ(3));
t860 = t743 * t849;
t865 = t751 * t765 + t750 + t860;
t604 = t642 * t715 + t761 * t855 + t764 * t865;
t706 = t806 / 0.2e1;
t686 = t785 + t706;
t689 = t791 - t809 / 0.2e1;
t713 = cos(pkin(6));
t652 = t686 * t765 + t689 * t849 + t713 * t762;
t653 = t686 * t685 - t713 * t688 - t689 * t718;
t669 = -t686 * t846 + t713 * t848;
t623 = -t652 * t764 + t653 * t715 - t669 * t761;
t575 = atan2(-t604, t623);
t570 = sin(t575);
t571 = cos(t575);
t555 = -t570 * t604 + t571 * t623;
t553 = 0.1e1 / t555 ^ 2;
t797 = t853 * t852;
t745 = -t851 * t763 + t797;
t741 = t745 * t718;
t796 = t853 * t850;
t744 = t766 * t851 + t796;
t793 = t851 * t847;
t776 = t688 * t793;
t645 = -t685 * t744 + t741 - t776;
t737 = t745 * t849;
t749 = t851 * t755;
t728 = t744 * t765 + t737 - t749;
t769 = t851 * t779;
t731 = -t744 * t846 - t769;
t609 = t645 * t715 + t728 * t764 + t731 * t761;
t603 = t609 ^ 2;
t551 = t553 * t603 + 0.1e1;
t705 = -t805 / 0.2e1;
t830 = t706 + t705;
t681 = t830 * qJD(2);
t657 = qJD(1) * t751 - qJD(2) * t797 - t851 * t681;
t682 = t766 * qJD(2);
t760 = -qJD(2) * t796 - t682 * t851;
t659 = -qJD(1) * t743 + t760;
t677 = (t784 - t803 / 0.2e1) * qJD(3);
t679 = t688 * qJD(3);
t596 = qJD(1) * t750 - qJD(3) * t741 + t657 * t765 - t659 * t849 - t677 * t744 + t679 * t793;
t676 = t762 * qJD(3);
t678 = t765 * qJD(3);
t597 = -qJD(1) * t777 - qJD(3) * t737 + t657 * t685 + t659 * t718 + t676 * t793 - t678 * t744;
t674 = (t782 - t801 / 0.2e1) * qJD(4);
t753 = -qJD(1) * t770 + t657 * t846;
t687 = t787 - t807 / 0.2e1;
t768 = t687 * qJD(4);
t717 = cos(qJ(4));
t827 = qJD(4) * t717;
t556 = -t596 * t764 + t597 * t715 + t645 * t827 + t674 * t728 + t731 * t768 + t753 * t761;
t839 = t556 * t553;
t602 = t604 ^ 2;
t621 = 0.1e1 / t623 ^ 2;
t574 = t602 * t621 + 0.1e1;
t572 = 0.1e1 / t574;
t660 = qJD(1) * t744 + qJD(2) * t795 - t681 * t853;
t780 = qJD(2) * t710 - t682 * t853;
t661 = qJD(1) * t745 - t780;
t600 = -qJD(1) * t749 + qJD(3) * t861 + t660 * t765 + t661 * t849 + t677 * t751 + t679 * t794;
t601 = qJD(1) * t776 + qJD(3) * t860 + t660 * t685 - t661 * t718 + t676 * t794 + t678 * t751;
t752 = qJD(1) * t769 + t660 * t846;
t558 = t600 * t764 - t601 * t715 + t642 * t827 + t674 * t865 - t752 * t761 + t768 * t855;
t680 = (t705 + t786) * qJD(2);
t683 = t689 * qJD(2);
t829 = qJD(3) * t718;
t630 = t686 * t677 + t713 * t679 + t680 * t849 + t683 * t765 + t689 * t829;
t812 = qJD(3) * t849;
t747 = -t713 * t676 - t686 * t678 + t680 * t718 - t683 * t685 - t689 * t812;
t754 = t761 * t846;
t568 = -t630 * t764 - t652 * t674 + t653 * t827 - t669 * t768 + t683 * t754 - t715 * t747;
t620 = 0.1e1 / t623;
t833 = t604 * t621;
t775 = -t558 * t620 + t568 * t833;
t542 = t775 * t572;
t778 = -t570 * t623 - t571 * t604;
t536 = t542 * t778 - t558 * t570 + t568 * t571;
t552 = 0.1e1 / t555;
t844 = t536 * t552 * t553;
t819 = 0.2e1 * (-t603 * t844 + t609 * t839) / t551 ^ 2;
t684 = t781 - t802 / 0.2e1;
t606 = t642 * t717 - t684 * t865 + t687 * t855;
t673 = t761 * qJD(4);
t675 = t764 * qJD(4);
t828 = qJD(4) * t715;
t560 = t600 * t684 + t601 * t717 + t642 * t828 + t673 * t855 + t675 * t865 + t687 * t752;
t862 = t568 * t621;
t711 = sin(pkin(8));
t712 = cos(pkin(8));
t629 = t711 * t728 - t712 * t731;
t714 = sin(qJ(5));
t716 = cos(qJ(5));
t725 = t645 * t717 - t684 * t728 + t687 * t731;
t585 = t629 * t714 + t716 * t725;
t579 = 0.1e1 / t585;
t580 = 0.1e1 / t585 ^ 2;
t854 = -0.2e1 * t604;
t557 = t596 * t684 + t597 * t717 - t645 * t828 - t673 * t731 - t675 * t728 + t687 * t753;
t586 = -t596 * t711 - t712 * t753;
t547 = qJD(5) * t585 + t557 * t714 - t586 * t716;
t584 = -t629 * t716 + t714 * t725;
t578 = t584 ^ 2;
t563 = t578 * t580 + 0.1e1;
t835 = t580 * t584;
t826 = qJD(5) * t584;
t548 = t557 * t716 + t586 * t714 - t826;
t841 = t548 * t579 * t580;
t843 = (t547 * t835 - t578 * t841) / t563 ^ 2;
t838 = t620 * t862;
t842 = (t558 * t833 - t602 * t838) / t574 ^ 2;
t840 = t553 * t609;
t837 = t570 * t609;
t836 = t571 * t609;
t834 = t604 * t620;
t832 = t711 * t714;
t831 = t711 * t716;
t818 = 0.2e1 * t844;
t817 = -0.2e1 * t843;
t816 = 0.2e1 * t843;
t815 = -0.2e1 * t842;
t814 = t620 * t842;
t813 = t584 * t841;
t811 = t687 * t846;
t810 = t712 * t846;
t800 = t609 * t818;
t799 = t838 * t854;
t798 = 0.2e1 * t813;
t628 = -t711 * t865 + t712 * t855;
t583 = -t606 * t716 + t628 * t714;
t582 = -t606 * t714 - t628 * t716;
t671 = -t830 * t851 - t797;
t736 = t744 * t849;
t650 = t671 * t765 + t736;
t740 = t744 * t718;
t651 = t671 * t685 - t740;
t616 = t650 * t684 + t651 * t717 + t671 * t811;
t634 = -t650 * t711 - t671 * t810;
t589 = t616 * t716 + t634 * t714;
t588 = t616 * t714 - t634 * t716;
t774 = -t714 * t579 + t716 * t835;
t624 = t652 * t684 + t653 * t717 - t669 * t687;
t773 = -t606 * t620 + t624 * t833;
t748 = t751 * t718;
t759 = -t830 * t853 + t795;
t649 = -t685 * t759 - t748;
t746 = t751 * t849;
t732 = t759 * t765 - t746;
t614 = t649 * t715 + t732 * t764 - t754 * t759;
t666 = -t686 * t849 + t689 * t765;
t667 = t689 * t685 + t686 * t718;
t627 = -t666 * t764 + t667 * t715 + t689 * t754;
t772 = -t614 * t620 + t627 * t833;
t617 = t642 * t764 - t715 * t865;
t625 = t652 * t715 + t653 * t764;
t771 = -t617 * t620 + t625 * t833;
t619 = -t645 * t684 - t717 * t728;
t590 = t619 * t714 - t645 * t831;
t591 = t619 * t716 + t645 * t832;
t767 = -t570 + (t571 * t834 + t570) * t572;
t758 = t768 * t846;
t727 = qJD(4) * t728;
t662 = qJD(1) * t671 + t780;
t658 = qJD(1) * t759 - t760;
t618 = t645 * t764 - t715 * t728;
t615 = -t650 * t764 + t651 * t715 + t671 * t754;
t613 = qJD(3) * t736 + t657 * t718 + t658 * t685 + t671 * t678;
t612 = qJD(3) * t740 - t657 * t849 + t658 * t765 + t671 * t677;
t592 = -t612 * t711 - t658 * t810;
t587 = -t600 * t711 - t712 * t752;
t577 = (t689 * t678 + t680 * t685 + t683 * t718 - t686 * t812) * t715 + t667 * t827 - (t689 * t677 + t680 * t765 - t683 * t849 - t686 * t829) * t764 - t666 * t674 + t680 * t754 + t689 * t758;
t576 = t630 * t715 + t652 * t827 + t653 * t674 - t747 * t764;
t569 = t630 * t684 + t652 * t675 - t653 * t828 + t669 * t673 + t683 * t811 - t717 * t747;
t567 = (qJD(3) * t746 - t660 * t718 + t662 * t685 - t678 * t759) * t715 + t649 * t827 - (qJD(3) * t748 + t660 * t849 + t662 * t765 - t677 * t759) * t764 + t732 * t674 + t662 * t754 - t759 * t758;
t566 = -t671 * t673 * t846 + t612 * t684 + t613 * t717 + t650 * t675 - t651 * t828 + t658 * t811;
t565 = -t600 * t715 - t601 * t764 + t642 * t674 - t827 * t865;
t564 = t596 * t717 - t597 * t684 - t645 * t675 + t715 * t727;
t561 = 0.1e1 / t563;
t549 = 0.1e1 / t551;
t546 = t772 * t572;
t545 = t771 * t572;
t544 = t773 * t572;
t541 = t767 * t609;
t539 = t546 * t778 - t570 * t614 + t571 * t627;
t538 = t545 * t778 - t570 * t617 + t571 * t625;
t537 = t544 * t778 - t570 * t606 + t571 * t624;
t534 = t772 * t815 + (t627 * t799 - t567 * t620 + (t558 * t627 + t568 * t614 + t577 * t604) * t621) * t572;
t533 = t771 * t815 + (t625 * t799 - t565 * t620 + (t558 * t625 + t568 * t617 + t576 * t604) * t621) * t572;
t532 = t773 * t815 + (t624 * t799 + t560 * t620 + (t558 * t624 + t568 * t606 + t569 * t604) * t621) * t572;
t1 = [0.2e1 * t609 * t814 + (-t556 * t620 + t609 * t862) * t572, t534, t533, t532, 0, 0; t604 * t552 * t819 + (-t558 * t552 + (t536 * t604 - t541 * t556) * t553) * t549 + ((t541 * t818 - t767 * t839) * t549 + (t541 * t819 + (-(-t542 * t572 * t834 + t815) * t837 - (t814 * t854 - t542 + (t542 - t775) * t572) * t836) * t549) * t553) * t609 (t539 * t840 - t552 * t615) * t819 + ((-t612 * t764 + t613 * t715 - t650 * t674 + t651 * t827 + t658 * t754 + t671 * t758) * t552 + t539 * t800 + (-t615 * t536 - t539 * t556 - (-t534 * t604 - t546 * t558 + t577 + (-t546 * t623 - t614) * t542) * t836 - (-t534 * t623 - t546 * t568 - t567 + (t546 * t604 - t627) * t542) * t837) * t553) * t549 (t538 * t840 - t552 * t618) * t819 + ((t596 * t715 + t597 * t764 + t645 * t674 - t717 * t727) * t552 + t538 * t800 + (-t618 * t536 - t538 * t556 - (-t533 * t604 - t545 * t558 + t576 + (-t545 * t623 - t617) * t542) * t836 - (-t533 * t623 - t545 * t568 - t565 + (t545 * t604 - t625) * t542) * t837) * t553) * t549 (t537 * t840 - t552 * t725) * t819 + (t557 * t552 + t537 * t800 + (-t725 * t536 - t537 * t556 - (-t532 * t604 - t544 * t558 + t569 + (-t544 * t623 - t606) * t542) * t836 - (-t532 * t623 - t544 * t568 + t560 + (t544 * t604 - t624) * t542) * t837) * t553) * t549, 0, 0; (-t579 * t582 + t583 * t835) * t816 + ((qJD(5) * t583 + t560 * t714 - t587 * t716) * t579 + t583 * t798 + (-t582 * t548 - (-qJD(5) * t582 + t560 * t716 + t587 * t714) * t584 - t583 * t547) * t580) * t561 (-t579 * t588 + t589 * t835) * t816 + ((qJD(5) * t589 + t566 * t714 - t592 * t716) * t579 + t589 * t798 + (-t588 * t548 - (-qJD(5) * t588 + t566 * t716 + t592 * t714) * t584 - t589 * t547) * t580) * t561 (-t579 * t590 + t591 * t835) * t816 + ((qJD(5) * t591 + t564 * t714 - t597 * t831) * t579 + t591 * t798 + (-t590 * t548 - (-qJD(5) * t590 + t564 * t716 + t597 * t832) * t584 - t591 * t547) * t580) * t561, t774 * t609 * t817 + (t774 * t556 + ((-qJD(5) * t579 - 0.2e1 * t813) * t716 + (t547 * t716 + (t548 - t826) * t714) * t580) * t609) * t561, t817 + 0.2e1 * (t547 * t580 * t561 + (-t561 * t841 - t580 * t843) * t584) * t584, 0;];
JaD_rot  = t1;
