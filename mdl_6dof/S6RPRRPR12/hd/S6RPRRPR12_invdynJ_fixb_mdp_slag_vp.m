% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR12_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR12_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR12_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:53:30
% EndTime: 2019-03-09 05:53:53
% DurationCPUTime: 18.48s
% Computational Cost: add. (14477->724), mult. (46067->947), div. (0->0), fcn. (40261->14), ass. (0->293)
t707 = sin(pkin(12));
t709 = cos(pkin(12));
t715 = cos(qJ(1));
t886 = cos(pkin(6));
t826 = t715 * t886;
t896 = sin(qJ(1));
t660 = t707 * t826 + t709 * t896;
t712 = sin(qJ(3));
t884 = sin(pkin(7));
t828 = t712 * t884;
t708 = sin(pkin(6));
t864 = t708 * t715;
t897 = cos(qJ(3));
t769 = t896 * t707 - t709 * t826;
t885 = cos(pkin(7));
t908 = t769 * t885;
t586 = -t660 * t897 + t712 * t908 + t828 * t864;
t711 = sin(qJ(4));
t714 = cos(qJ(4));
t831 = t708 * t885;
t904 = -t715 * t831 + t769 * t884;
t555 = t586 * t711 + t714 * t904;
t805 = t884 * t897;
t583 = t660 * t712 + t805 * t864 + t897 * t908;
t710 = sin(qJ(6));
t713 = cos(qJ(6));
t921 = t555 * t710 - t583 * t713;
t920 = t555 * t713 + t583 * t710;
t556 = t586 * t714 - t711 * t904;
t829 = t712 * t885;
t743 = t708 * (-t707 * t829 + t709 * t897);
t773 = qJD(3) * t805;
t917 = qJD(1) * t743 - t773;
t752 = t707 * t897 + t709 * t829;
t790 = t886 * t884;
t771 = t712 * t790;
t633 = t708 * t752 + t771;
t623 = t633 * qJD(3);
t758 = t897 * t790;
t806 = t885 * t897;
t779 = t709 * t806;
t766 = t708 * t779;
t865 = t707 * t712;
t632 = t708 * t865 - t758 - t766;
t569 = qJD(1) * t623 + t632 * qJDD(1);
t791 = t886 * t885;
t830 = t708 * t884;
t659 = t709 * t830 - t791;
t746 = qJDD(1) * t659 - qJDD(3);
t916 = -t569 * MDP(11) - t746 * MDP(12);
t624 = t633 * qJD(1);
t804 = qJD(1) * t830;
t910 = qJD(1) * t791 - t709 * t804 + qJD(3);
t577 = t714 * t624 + t711 * t910;
t574 = qJD(6) + t577;
t693 = t709 * t708 * qJ(2);
t824 = qJD(1) * t886;
t813 = pkin(1) * t824;
t657 = qJD(1) * t693 + t707 * t813;
t733 = (t709 * t831 + t790) * pkin(9);
t605 = qJD(1) * t733 + t657;
t691 = t709 * t813;
t866 = t707 * t708;
t735 = t886 * pkin(2) + (-pkin(9) * t885 - qJ(2)) * t866;
t614 = qJD(1) * t735 + t691;
t913 = t712 * t605 - t614 * t806;
t912 = t707 * t806 + t709 * t712;
t679 = qJD(1) * t758;
t852 = qJD(1) * t708;
t835 = t707 * t852;
t620 = -qJD(1) * t766 + t712 * t835 - t679;
t612 = qJD(4) + t620;
t650 = (-pkin(9) * t707 * t884 - pkin(2) * t709 - pkin(1)) * t708;
t640 = qJD(1) * t650 + qJD(2);
t571 = -t614 * t884 + t885 * t640;
t523 = t620 * pkin(3) - t624 * pkin(10) + t571;
t544 = t897 * t605 + t614 * t829 + t640 * t828;
t526 = pkin(10) * t910 + t544;
t487 = -t714 * t523 + t526 * t711;
t843 = qJD(5) + t487;
t568 = qJDD(4) + t569;
t705 = t708 ^ 2;
t909 = t705 * (t707 ^ 2 + t709 ^ 2);
t665 = t711 * t885 + t714 * t828;
t777 = t707 * t804;
t858 = -qJD(4) * t665 + t711 * t917 - t714 * t777;
t837 = pkin(1) * t886;
t663 = t707 * t837 + t693;
t626 = t733 + t663;
t698 = t709 * t837;
t634 = t698 + t735;
t731 = -t712 * t626 + t634 * t806 + t650 * t805;
t664 = t711 * t828 - t714 * t885;
t857 = -qJD(4) * t664 - t711 * t777 - t714 * t917;
t751 = t779 - t865;
t719 = (qJD(1) * qJD(3) * t751 + qJDD(1) * t752) * t708 + qJD(3) * t679 + qJDD(1) * t771;
t850 = qJD(4) * t577;
t514 = t711 * t719 + t714 * t746 + t850;
t807 = t886 * t896;
t755 = t715 * t707 + t709 * t807;
t907 = t755 * t885 - t896 * t830;
t802 = qJD(3) * t828;
t767 = -t852 * t912 + t802;
t849 = qJD(4) * t711;
t867 = t620 * t711;
t906 = -qJD(5) * t711 - t544 + (t849 + t867) * pkin(4);
t844 = pkin(5) * t577 + t843;
t661 = -t707 * t807 + t715 * t709;
t587 = t661 * t712 + t897 * t907;
t761 = g(1) * t587 + g(2) * t583 + g(3) * t632;
t488 = t711 * t523 + t714 * t526;
t485 = -qJ(5) * t612 - t488;
t742 = t714 * t910;
t575 = t624 * t711 - t742;
t893 = pkin(5) * t575;
t475 = -t485 - t893;
t513 = -qJD(4) * t742 + t624 * t849 + t711 * t746 - t714 * t719;
t512 = -qJDD(6) + t513;
t899 = pkin(4) + pkin(11);
t905 = t899 * t512 + (t475 - t488 + t893) * t574;
t901 = t577 ^ 2;
t900 = t612 ^ 2;
t898 = pkin(5) + pkin(10);
t895 = pkin(1) * t705;
t894 = pkin(4) * t568;
t892 = pkin(10) * t568;
t581 = t633 * t711 + t659 * t714;
t888 = t581 * pkin(4);
t887 = pkin(10) * qJD(4);
t883 = pkin(1) * qJDD(1);
t882 = qJ(5) * t575;
t881 = qJ(5) * t714;
t846 = qJD(6) * t713;
t839 = t710 * t514 + t713 * t568 + t575 * t846;
t847 = qJD(6) * t710;
t480 = -t612 * t847 + t839;
t880 = t480 * t713;
t879 = t485 * t612;
t878 = t488 * t612;
t877 = t512 * t710;
t868 = t612 * t710;
t545 = -t713 * t575 + t868;
t876 = t545 * t574;
t875 = t545 * t612;
t547 = t575 * t710 + t612 * t713;
t874 = t547 * t574;
t873 = t547 * t612;
t565 = t568 * qJ(5);
t872 = t575 * t612;
t871 = t577 * t575;
t870 = t577 * t612;
t582 = t633 * t714 - t659 * t711;
t869 = t582 * qJ(5);
t862 = t710 * t711;
t861 = t711 * t713;
t509 = t713 * t512;
t579 = -t634 * t884 + t885 * t650;
t535 = t632 * pkin(3) - t633 * pkin(10) + t579;
t838 = t897 * t626 + t634 * t829 + t650 * t828;
t542 = -pkin(10) * t659 + t838;
t860 = t711 * t535 + t714 * t542;
t543 = t640 * t805 - t913;
t570 = pkin(3) * t624 + pkin(10) * t620;
t859 = t714 * t543 + t711 * t570;
t848 = qJD(4) * t714;
t856 = -qJ(5) * t848 - t620 * t881 + t906;
t498 = -qJ(5) * t624 - t859;
t855 = -pkin(5) * t867 - t898 * t849 + t498;
t836 = t708 * t896;
t854 = t715 * pkin(1) + qJ(2) * t836;
t842 = qJD(1) * qJD(2);
t841 = g(1) * t896;
t687 = t898 * t714;
t812 = qJDD(1) * t837;
t833 = t708 * t842;
t642 = qJDD(1) * t693 + t707 * t812 + t709 * t833;
t834 = qJD(3) * t897;
t832 = -qJ(5) * t711 - pkin(3);
t717 = qJD(1) ^ 2;
t825 = t717 * t886;
t593 = qJDD(1) * t733 + t642;
t689 = t709 * t812;
t594 = qJDD(1) * t735 - t707 * t833 + t689;
t636 = qJDD(1) * t650 + qJDD(2);
t749 = qJD(3) * t913 - t897 * t593 - t594 * t829 - t636 * t828 - t640 * t773;
t494 = -pkin(10) * t746 - t749;
t561 = -t594 * t884 + t885 * t636;
t504 = t569 * pkin(3) - pkin(10) * t719 + t561;
t815 = -t711 * t494 + t714 * t504 - t523 * t849 - t526 * t848;
t782 = qJDD(5) - t815;
t461 = -pkin(5) * t513 - t568 * t899 + t782;
t803 = qJD(3) * t829;
t768 = t712 * t593 - t594 * t806 + t605 * t834 + t614 * t803 - t636 * t805 + t640 * t802;
t495 = pkin(3) * t746 + t768;
t722 = t513 * qJ(5) - t577 * qJD(5) + t495;
t465 = t514 * t899 + t722;
t823 = t713 * t461 - t710 * t465;
t822 = -t713 * t514 + t568 * t710;
t821 = t535 * t714 - t711 * t542;
t820 = t612 * t714;
t819 = t574 * t710;
t818 = t574 * t713;
t816 = -t714 * t494 - t711 * t504 - t523 * t848 + t526 * t849;
t808 = -pkin(1) * t896 + qJ(2) * t864;
t588 = t661 * t897 - t712 * t907;
t724 = -t755 * t884 - t831 * t896;
t557 = t588 * t711 + t714 * t724;
t801 = g(1) * t555 + g(2) * t557;
t558 = t588 * t714 - t711 * t724;
t800 = -g(1) * t556 - g(2) * t558;
t799 = -g(1) * t583 + g(2) * t587;
t492 = -qJ(5) * t632 - t860;
t798 = qJD(2) * t824;
t796 = g(2) * t864 - g(3) * t886;
t562 = t620 * t861 + t624 * t710;
t795 = t713 * t849 + t562;
t563 = -t620 * t862 + t624 * t713;
t794 = t710 * t849 - t563;
t539 = t711 * t543;
t666 = -t714 * t899 + t832;
t793 = qJD(6) * t666 + t539 + (-pkin(5) * t620 - t570) * t714 - t899 * t624 - qJD(4) * t687;
t686 = t898 * t711;
t792 = -qJD(6) * t686 - t906 - t612 * (pkin(11) * t711 - t881);
t530 = qJD(2) * t708 * t912 + t626 * t834 + t634 * t803 + t650 * t802;
t788 = t710 * t461 + t713 * t465;
t474 = -t612 * t899 + t844;
t525 = -pkin(3) * t910 - t543;
t720 = -t577 * qJ(5) + t525;
t483 = t575 * t899 + t720;
t466 = t474 * t713 - t483 * t710;
t467 = t474 * t710 + t483 * t713;
t477 = pkin(5) * t582 - t632 * t899 - t821;
t541 = t659 * pkin(3) - t731;
t725 = t541 - t869;
t489 = t581 * t899 + t725;
t787 = t477 * t713 - t489 * t710;
t786 = t477 * t710 + t489 * t713;
t784 = t581 * t713 - t632 * t710;
t552 = t581 * t710 + t632 * t713;
t783 = (-qJ(2) * t835 + t691) * t707 - t657 * t709;
t776 = t707 * qJD(2) * t830;
t775 = pkin(4) * t714 - t832;
t529 = qJD(2) * t743 + qJD(3) * t731;
t622 = (t708 * t751 + t758) * qJD(3);
t560 = t623 * pkin(3) - t622 * pkin(10) + t776;
t774 = -t711 * t529 - t535 * t849 - t542 * t848 + t560 * t714;
t772 = g(1) * t715 + g(2) * t896;
t606 = t612 * qJD(5);
t463 = -t565 - t606 + t816;
t765 = t525 * t612 - t892;
t497 = t575 * pkin(4) + t720;
t764 = -t497 * t612 + t892;
t763 = t714 * t529 + t535 * t848 - t542 * t849 + t711 * t560;
t762 = -g(1) * t558 + g(2) * t556 - g(3) * t582;
t760 = -g(1) * t588 + g(2) * t586 - g(3) * t633;
t754 = t713 * t664 + t710 * t805;
t753 = -t710 * t664 + t713 * t805;
t745 = -t513 + t872;
t741 = t612 * t887 - t761;
t550 = -t633 * t849 + (-qJD(4) * t659 + t622) * t714;
t740 = -qJ(5) * t550 - qJD(5) * t582 + t530;
t738 = g(1) * t557 - g(2) * t555 + g(3) * t581 + t815;
t737 = t762 - t816;
t468 = t514 * pkin(4) + t722;
t736 = t468 + t741;
t471 = -qJ(5) * t623 - qJD(5) * t632 - t763;
t462 = -pkin(5) * t514 - t463;
t730 = t462 + (t574 * t899 + t882) * t574 + t762;
t727 = t497 * t577 + qJDD(5) - t738;
t692 = -t708 * t883 + qJDD(2);
t662 = -qJ(2) * t866 + t698;
t641 = t689 + (-qJ(2) * qJDD(1) - t842) * t866;
t549 = qJD(4) * t582 + t622 * t711;
t533 = pkin(4) * t577 + t882;
t516 = t557 * t710 + t587 * t713;
t515 = t557 * t713 - t587 * t710;
t506 = qJD(6) * t784 + t549 * t710 + t623 * t713;
t505 = qJD(6) * t552 - t549 * t713 + t623 * t710;
t503 = t725 + t888;
t500 = -pkin(4) * t624 - t570 * t714 + t539;
t493 = -pkin(4) * t632 - t821;
t484 = -pkin(4) * t612 + t843;
t482 = -pkin(5) * t581 - t492;
t481 = qJD(6) * t547 + t822;
t476 = pkin(4) * t549 + t740;
t473 = t549 * t899 + t740;
t472 = -pkin(4) * t623 - t774;
t470 = -pkin(5) * t549 - t471;
t469 = pkin(5) * t550 - t623 * t899 - t774;
t464 = t782 - t894;
t459 = -t467 * qJD(6) + t823;
t458 = t466 * qJD(6) + t788;
t1 = [(-t642 * t886 - g(1) * t769 + g(2) * t755 + (t692 * t707 - t709 * t798) * t708 + (-t663 * t886 - t707 * t895) * qJDD(1)) * MDP(5) + (t641 * t886 + g(1) * t660 - g(2) * t661 + (-t692 * t709 - t707 * t798) * t708 + (t662 * t886 + t709 * t895) * qJDD(1)) * MDP(4) + (t624 * t622 + t633 * t719) * MDP(8) + (-t633 * t569 - t622 * t620 - t624 * t623 - t632 * t719) * MDP(9) + (t463 * t581 + t464 * t582 + t471 * t575 + t472 * t577 + t484 * t550 + t485 * t549 + t492 * t514 - t493 * t513 - t799) * MDP(22) + (t464 * t632 - t468 * t581 + t472 * t612 - t476 * t575 + t484 * t623 + t493 * t568 - t497 * t549 - t503 * t514 - t800) * MDP(23) + (-t463 * t632 - t468 * t582 - t471 * t612 - t476 * t577 - t485 * t623 - t492 * t568 - t497 * t550 + t503 * t513 - t801) * MDP(24) + (-t488 * t623 + t495 * t582 - t541 * t513 + t525 * t550 + t530 * t577 - t568 * t860 - t612 * t763 + t632 * t816 + t801) * MDP(21) + (-t481 * t582 - t505 * t574 - t512 * t784 - t545 * t550) * MDP(29) + (t480 * t784 - t481 * t552 - t505 * t547 - t506 * t545) * MDP(27) + (t842 * t909 + (-t641 * t707 + t642 * t709 + (-t662 * t707 + t663 * t709) * qJDD(1) - t772) * t708) * MDP(6) + (t642 * t663 + t641 * t662 - g(1) * t808 - g(2) * t854 + (-t692 * pkin(1) - qJD(2) * t783) * t708) * MDP(7) + (t468 * t503 + t497 * t476 + t463 * t492 + t485 * t471 + t464 * t493 + t484 * t472 - g(1) * (-t660 * pkin(2) + t586 * pkin(3) + t556 * pkin(4) - pkin(10) * t583 + t555 * qJ(5) + t808) - g(2) * (t661 * pkin(2) + t588 * pkin(3) + t558 * pkin(4) + t587 * pkin(10) + t557 * qJ(5) + t854) + (g(1) * t904 + g(2) * t724) * pkin(9)) * MDP(25) + (-g(1) * t586 - g(2) * t588 - t530 * t910 + t561 * t632 + t579 * t569 + t571 * t623 + t620 * t776 - t731 * t746) * MDP(13) + (t622 * t910 - t633 * t746) * MDP(10) + (-t529 * t910 + t561 * t633 + t571 * t622 + t579 * t719 + t624 * t776 + t746 * t838 + t799) * MDP(14) + (-t623 * t910 + t632 * t746) * MDP(11) + (-g(2) * t715 + t841) * MDP(2) + (-(qJD(6) * t787 + t469 * t710 + t473 * t713) * t574 + t786 * t512 - t458 * t582 - t467 * t550 + t470 * t547 + t482 * t480 + t462 * t552 + t475 * t506 - g(1) * t920 - g(2) * t515) * MDP(32) + ((-qJD(6) * t786 + t469 * t713 - t473 * t710) * t574 - t787 * t512 + t459 * t582 + t466 * t550 + t470 * t545 + t482 * t481 - t462 * t784 + t475 * t505 - g(1) * t921 - g(2) * t516) * MDP(31) + (-t487 * t623 + t495 * t581 + t541 * t514 + t525 * t549 + t530 * t575 + t568 * t821 + t612 * t774 + t632 * t815 + t800) * MDP(20) + qJDD(1) * MDP(1) + t772 * MDP(3) + (-t719 * MDP(10) + t768 * MDP(13) - t749 * MDP(14) - t916) * t659 + (-t514 * t632 - t549 * t612 - t568 * t581 - t575 * t623) * MDP(18) + (-t513 * t632 + t550 * t612 + t568 * t582 + t577 * t623) * MDP(17) + (t568 * t632 + t612 * t623) * MDP(19) + (t513 * t581 - t514 * t582 - t549 * t577 - t550 * t575) * MDP(16) + (-t512 * t582 + t550 * t574) * MDP(30) + (t480 * t582 + t506 * t574 - t512 * t552 + t547 * t550) * MDP(28) + (-t513 * t582 + t550 * t577) * MDP(15) + (t480 * t552 + t506 * t547) * MDP(26); -t717 * MDP(6) * t909 + (qJDD(2) + t796) * MDP(7) + (t569 * t885 - t620 * t777 - t746 * t805 - t767 * t910) * MDP(13) + (-t624 * t777 + t719 * t885 + t746 * t828 + t910 * t917) * MDP(14) + (-t513 * t664 - t514 * t665 - t575 * t857 - t577 * t858) * MDP(22) + (-g(1) * t836 - t463 * t665 + t464 * t664 - t468 * t805 - t484 * t858 - t485 * t857 + t497 * t767 + t796) * MDP(25) + (t665 * t481 - t754 * t512 + (t753 * qJD(6) - t710 * t767 - t713 * t858) * t574 + t857 * t545) * MDP(31) + (t665 * t480 - t753 * t512 + (-t754 * qJD(6) + t710 * t858 - t713 * t767) * t574 + t857 * t547) * MDP(32) + (MDP(20) - MDP(23)) * (-t514 * t805 - t664 * t568 + t575 * t767 + t612 * t858) + (-MDP(24) + MDP(21)) * (t513 * t805 - t665 * t568 + t577 * t767 - t612 * t857) + ((-qJDD(1) * t709 + t707 * t825) * MDP(4) + (qJDD(1) * t707 + t709 * t825) * MDP(5) + (qJD(1) * t783 - t841 - t883) * MDP(7)) * t708; (t545 * t563 + t547 * t562 + (-t545 * t710 + t547 * t713) * t849 + (-t880 + t481 * t710 + (t545 * t713 + t547 * t710) * qJD(6)) * t714) * MDP(27) + t916 + ((-t513 - t872) * t714 + (-t514 - t870) * t711) * MDP(16) + (-t481 * t711 + t795 * t574 + (t574 * t847 + t509 - t875) * t714) * MDP(29) + (t568 * t711 + t612 * t820) * MDP(17) + (t480 * t711 + t794 * t574 + (-t574 * t846 + t873 + t877) * t714) * MDP(28) + (-t498 * t575 - t500 * t577 + (-t463 + t612 * t484 + (-t514 + t850) * pkin(10)) * t714 + (t464 + t879 + (qJD(4) * t575 - t513) * pkin(10)) * t711 + t760) * MDP(22) + (-t500 * t612 + t514 * t775 - t575 * t856 + t711 * t764 + t714 * t736) * MDP(23) + (t498 * t612 - t513 * t775 - t577 * t856 - t711 * t736 + t714 * t764) * MDP(24) + (-pkin(3) * t514 + t539 * t612 - t544 * t575 + t765 * t711 + (-t495 + (-t570 - t887) * t612 + t761) * t714) * MDP(20) + (pkin(3) * t513 + t859 * t612 - t544 * t577 + t765 * t714 + (t495 + t741) * t711) * MDP(21) + ((t666 * t713 + t686 * t710) * t512 - t458 * t711 + t687 * t480 - g(1) * (-t587 * t861 - t588 * t710) - g(2) * (-t583 * t861 + t586 * t710) - g(3) * (-t632 * t861 - t633 * t710) + (t710 * t793 + t713 * t792) * t574 + t855 * t547 + t794 * t475 + (-t462 * t710 - t467 * t612 - t475 * t846) * t714) * MDP(32) + (-(-t666 * t710 + t686 * t713) * t512 + t459 * t711 + t687 * t481 - g(1) * (-t587 * t862 + t588 * t713) - g(2) * (-t583 * t862 - t586 * t713) - g(3) * (-t632 * t862 + t633 * t713) + (t710 * t792 - t713 * t793) * t574 + t855 * t545 - t795 * t475 + (t462 * t713 + t466 * t612 - t475 * t847) * t714) * MDP(31) + (t620 * t910 + t719) * MDP(10) + (t543 * t910 + t571 * t620 + t749 - t760) * MDP(14) + (t544 * t910 + t761 - t768) * MDP(13) - t620 ^ 2 * MDP(9) + (-t484 * t500 - t485 * t498 + t856 * t497 + (-t463 * t714 + t464 * t711 + (t484 * t714 + t485 * t711) * qJD(4) + t760) * pkin(10) + (-t468 + t761) * t775) * MDP(25) + (-t513 * t711 + t577 * t820) * MDP(15) + (-t512 * t711 + t574 * t820) * MDP(30) + (-t480 * t710 * t714 + (-t714 * t846 + t794) * t547) * MDP(26) + (MDP(11) * t910 - t571 * MDP(13) - t577 * MDP(17) + t575 * MDP(18) - t612 * MDP(19) + t487 * MDP(20) + t488 * MDP(21) - t484 * MDP(23) + t485 * MDP(24) + t620 * MDP(8) + MDP(9) * t624) * t624 + (t568 * t714 - t711 * t900) * MDP(18); MDP(15) * t871 + (-t575 ^ 2 + t901) * MDP(16) + t745 * MDP(17) + (t870 - t514) * MDP(18) + t568 * MDP(19) + (-t525 * t577 + t738 + t878) * MDP(20) + (-t487 * t612 + t525 * t575 - t737) * MDP(21) + (pkin(4) * t513 - qJ(5) * t514 + (-t485 - t488) * t577 + (t484 - t843) * t575) * MDP(22) + (t533 * t575 + t727 - t878 - 0.2e1 * t894) * MDP(23) + (-t497 * t575 + t533 * t577 + t612 * t843 + 0.2e1 * t565 + t606 + t737) * MDP(24) + (-t463 * qJ(5) - t464 * pkin(4) - t497 * t533 - t484 * t488 - g(1) * (-pkin(4) * t557 + qJ(5) * t558) - g(2) * (pkin(4) * t555 - qJ(5) * t556) - g(3) * (t869 - t888) - t843 * t485) * MDP(25) + (-t547 * t819 + t880) * MDP(26) + ((-t481 - t874) * t713 + (-t480 + t876) * t710) * MDP(27) + (t547 * t575 - t574 * t819 - t509) * MDP(28) + (-t545 * t575 - t574 * t818 + t877) * MDP(29) + t574 * t575 * MDP(30) + (qJ(5) * t481 + t466 * t575 + t844 * t545 + t730 * t710 + t713 * t905) * MDP(31) + (qJ(5) * t480 - t467 * t575 + t844 * t547 - t710 * t905 + t730 * t713) * MDP(32); t745 * MDP(22) + (-t900 - t901) * MDP(24) + (t727 + t879 - t894) * MDP(25) + (-t509 - t875) * MDP(31) + (-t873 + t877) * MDP(32) + (-MDP(31) * t819 - MDP(32) * t818) * t574 + (-t871 + t568) * MDP(23); t547 * t545 * MDP(26) + (-t545 ^ 2 + t547 ^ 2) * MDP(27) + (t839 + t876) * MDP(28) + (-t822 + t874) * MDP(29) - t512 * MDP(30) + (-g(1) * t515 + g(2) * t920 - g(3) * t784 + t467 * t574 - t475 * t547 + t823) * MDP(31) + (g(1) * t516 - g(2) * t921 + g(3) * t552 + t466 * t574 + t475 * t545 - t788) * MDP(32) + (-MDP(28) * t868 - MDP(29) * t547 - MDP(31) * t467 - MDP(32) * t466) * qJD(6);];
tau  = t1;
