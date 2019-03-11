% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:41:57
% EndTime: 2019-03-09 16:42:12
% DurationCPUTime: 14.12s
% Computational Cost: add. (13017->657), mult. (28770->796), div. (0->0), fcn. (21237->14), ass. (0->283)
t731 = sin(qJ(2));
t885 = cos(qJ(3));
t808 = qJD(1) * t885;
t730 = sin(qJ(3));
t733 = cos(qJ(2));
t843 = t730 * t733;
t650 = -qJD(1) * t843 - t731 * t808;
t722 = qJD(2) + qJD(3);
t726 = sin(pkin(10));
t727 = cos(pkin(10));
t623 = t650 * t727 - t722 * t726;
t729 = sin(qJ(5));
t884 = cos(qJ(5));
t637 = t726 * t650;
t898 = t722 * t727 + t637;
t776 = t884 * t898;
t568 = t729 * t623 + t776;
t916 = t568 ^ 2;
t827 = qJD(1) * t731;
t648 = t730 * t827 - t733 * t808;
t642 = qJD(5) + t648;
t915 = t568 * t642;
t667 = t731 * t885 + t843;
t617 = t722 * t667;
t811 = t885 * t733;
t822 = qJDD(1) * t731;
t577 = qJD(1) * t617 - qJDD(1) * t811 + t730 * t822;
t574 = qJDD(5) + t577;
t784 = t729 * t898;
t912 = -t623 * t884 + t784;
t886 = t912 ^ 2;
t914 = t642 * t912;
t735 = -pkin(8) - pkin(7);
t823 = qJD(1) * qJD(2);
t804 = t733 * t823;
t619 = qJDD(2) * pkin(2) - t735 * (-t804 - t822);
t805 = t731 * t823;
t821 = qJDD(1) * t733;
t626 = t735 * (-t805 + t821);
t684 = t735 * t731;
t669 = qJD(1) * t684;
t875 = qJD(2) * pkin(2);
t657 = t669 + t875;
t685 = t735 * t733;
t671 = qJD(1) * t685;
t807 = qJD(3) * t885;
t826 = qJD(3) * t730;
t799 = -t885 * t619 - t730 * t626 + t657 * t826 - t671 * t807;
t720 = qJDD(2) + qJDD(3);
t878 = t720 * pkin(3);
t530 = qJDD(4) + t799 - t878;
t725 = qJ(2) + qJ(3);
t718 = cos(t725);
t710 = g(3) * t718;
t717 = sin(t725);
t734 = cos(qJ(1));
t850 = t717 * t734;
t732 = sin(qJ(1));
t851 = t717 * t732;
t907 = g(1) * t850 + g(2) * t851;
t888 = t710 - t907;
t913 = t530 + t888;
t665 = t726 * t884 + t729 * t727;
t647 = t665 * qJD(5);
t909 = t665 * t648 + t647;
t767 = -t729 * t726 + t884 * t727;
t806 = qJD(5) * t884;
t825 = qJD(5) * t729;
t891 = -t726 * t825 + t727 * t806;
t908 = t767 * t648 + t891;
t703 = pkin(2) * t805;
t772 = -t730 * t731 + t811;
t616 = t722 * t772;
t742 = t616 * qJD(1);
t876 = t733 * pkin(2);
t714 = pkin(1) + t876;
t895 = -qJ(4) * t667 - t714;
t513 = t577 * pkin(3) - qJ(4) * t742 + t650 * qJD(4) + qJDD(1) * t895 + t703;
t753 = t730 * t619 - t626 * t885 + t657 * t807 + t671 * t826;
t526 = t720 * qJ(4) + t722 * qJD(4) + t753;
t483 = t727 * t513 - t726 * t526;
t484 = t726 * t513 + t727 * t526;
t783 = -t483 * t726 + t484 * t727;
t906 = -t530 - t710;
t608 = -pkin(3) * t772 + t895;
t628 = t730 * t684 - t685 * t885;
t557 = t727 * t608 - t628 * t726;
t857 = t667 * t727;
t536 = -pkin(4) * t772 - pkin(9) * t857 + t557;
t558 = t726 * t608 + t727 * t628;
t858 = t667 * t726;
t546 = -pkin(9) * t858 + t558;
t905 = t729 * t536 + t884 * t546;
t605 = -pkin(3) * t650 + qJ(4) * t648;
t591 = pkin(2) * t827 + t605;
t653 = t730 * t671;
t614 = t669 * t885 + t653;
t549 = t726 * t591 + t727 * t614;
t697 = pkin(2) * t807 + qJD(4);
t904 = -t697 * t727 + t549;
t548 = t727 * t591 - t614 * t726;
t859 = t648 * t727;
t797 = -t650 * pkin(4) + pkin(9) * t859;
t521 = t548 + t797;
t860 = t648 * t726;
t820 = pkin(9) * t860;
t535 = t820 + t549;
t706 = pkin(2) * t730 + qJ(4);
t655 = (-pkin(9) - t706) * t726;
t719 = t727 * pkin(9);
t855 = t706 * t727;
t656 = t719 + t855;
t769 = t655 * t884 - t729 * t656;
t903 = -qJD(5) * t769 + t729 * t521 + t884 * t535 - t697 * t767;
t602 = t729 * t655 + t656 * t884;
t902 = -qJD(5) * t602 - t521 * t884 + t729 * t535 - t665 * t697;
t635 = pkin(4) * t860;
t901 = pkin(5) * t909 - qJ(6) * t908 - qJD(6) * t665 + t635;
t609 = t657 * t885 + t653;
t550 = t727 * t605 - t609 * t726;
t525 = t550 + t797;
t551 = t726 * t605 + t727 * t609;
t537 = t820 + t551;
t728 = -pkin(9) - qJ(4);
t681 = t728 * t726;
t873 = qJ(4) * t727;
t682 = t719 + t873;
t768 = t681 * t884 - t729 * t682;
t900 = -qJD(4) * t767 - qJD(5) * t768 + t729 * t525 + t884 * t537;
t622 = t729 * t681 + t682 * t884;
t899 = -qJD(4) * t665 - qJD(5) * t622 - t525 * t884 + t729 * t537;
t896 = t718 * pkin(3) + t717 * qJ(4);
t741 = t667 * qJDD(1) + t742;
t740 = t726 * t720 + t727 * t741;
t572 = t726 * t741;
t802 = t720 * t727 - t572;
t894 = t726 * t740 + t727 * t802;
t892 = -t623 * t726 + t727 * t898;
t890 = t885 * t684 + t730 * t685;
t654 = t885 * t671;
t613 = t730 * t669 - t654;
t889 = -pkin(2) * t826 + t613;
t791 = g(1) * t734 + g(2) * t732;
t818 = t731 * t875;
t547 = pkin(3) * t617 - qJ(4) * t616 - qJD(4) * t667 + t818;
t813 = qJD(2) * t735;
t670 = t731 * t813;
t672 = t733 * t813;
t569 = qJD(3) * t890 + t885 * t670 + t730 * t672;
t505 = t727 * t547 - t569 * t726;
t863 = t616 * t727;
t492 = pkin(4) * t617 - pkin(9) * t863 + t505;
t506 = t726 * t547 + t727 * t569;
t864 = t616 * t726;
t501 = -pkin(9) * t864 + t506;
t887 = -qJD(5) * t905 + t492 * t884 - t729 * t501;
t883 = pkin(2) * t731;
t882 = pkin(5) * t574;
t709 = g(3) * t717;
t879 = t650 * pkin(5);
t877 = t727 * pkin(4);
t872 = qJ(6) * t574;
t683 = t714 * qJD(1);
t587 = pkin(3) * t648 + qJ(4) * t650 - t683;
t610 = t730 * t657 - t654;
t593 = qJ(4) * t722 + t610;
t539 = t727 * t587 - t593 * t726;
t514 = pkin(4) * t648 + pkin(9) * t623 + t539;
t540 = t726 * t587 + t727 * t593;
t519 = pkin(9) * t898 + t540;
t481 = t729 * t514 + t519 * t884;
t871 = t481 * t642;
t869 = t912 * t568;
t868 = t577 * t726;
t867 = t769 * t574;
t866 = t602 * t574;
t865 = t610 * t722;
t862 = t768 * t574;
t861 = t622 * t574;
t721 = pkin(10) + qJ(5);
t715 = sin(t721);
t854 = t715 * t718;
t716 = cos(t721);
t853 = t716 * t718;
t852 = t717 * t728;
t707 = pkin(3) + t877;
t676 = t718 * t707;
t849 = t718 * t728;
t848 = t718 * t732;
t847 = t718 * t734;
t842 = t732 * t716;
t841 = t734 * t715;
t639 = t650 * qJ(6);
t837 = t639 - t903;
t836 = -t879 - t902;
t835 = -t889 + t901;
t834 = -t610 + t901;
t833 = t639 - t900;
t832 = -t879 - t899;
t723 = t731 ^ 2;
t828 = -t733 ^ 2 + t723;
t480 = t514 * t884 - t729 * t519;
t824 = qJD(6) - t480;
t819 = t885 * pkin(2);
t816 = -g(3) * t854 + t715 * t907;
t815 = pkin(5) * t853 + qJ(6) * t854 + t676;
t814 = g(1) * t847 + g(2) * t848 + t709;
t803 = pkin(4) * t726 - t735;
t470 = t577 * pkin(4) - pkin(9) * t740 + t483;
t474 = pkin(9) * t802 + t484;
t800 = -t884 * t470 + t729 * t474 + t514 * t825 + t519 * t806;
t713 = -t819 - pkin(3);
t796 = -g(1) * t851 + g(2) * t850;
t795 = t635 - t889;
t794 = -pkin(3) * t717 - t883;
t629 = t715 * t848 + t716 * t734;
t631 = t718 * t841 - t842;
t793 = -g(1) * t629 + g(2) * t631;
t630 = t718 * t842 - t841;
t632 = t715 * t732 + t716 * t847;
t792 = g(1) * t630 - g(2) * t632;
t790 = g(1) * t732 - g(2) * t734;
t789 = t539 * t650 + t727 * t907;
t782 = t676 - t852;
t781 = -g(3) * t853 + t716 * t907;
t589 = pkin(4) * t858 - t890;
t779 = pkin(5) * t716 + qJ(6) * t715 + t707;
t778 = t790 * t718;
t777 = -0.2e1 * pkin(1) * t823 - pkin(7) * qJDD(2);
t773 = t536 * t884 - t729 * t546;
t570 = t730 * t670 - t672 * t885 + t684 * t826 - t685 * t807;
t764 = t729 * t470 + t884 * t474 + t514 * t806 - t519 * t825;
t763 = t729 * t492 + t884 * t501 + t536 * t806 - t546 * t825;
t497 = -qJD(5) * t776 - t623 * t825 - t729 * t802 - t884 * t740;
t538 = pkin(4) * t864 + t570;
t604 = -pkin(5) * t767 - t665 * qJ(6) - t707;
t590 = -t722 * pkin(3) + qJD(4) - t609;
t761 = -t539 * t859 - t540 * t860 + t783 - t814;
t736 = qJD(2) ^ 2;
t759 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t736 + t790;
t737 = qJD(1) ^ 2;
t758 = pkin(1) * t737 - pkin(7) * qJDD(1) + t791;
t757 = g(1) * t631 + g(2) * t629 + t715 * t709 - t800;
t756 = -t683 * t650 - t799 - t888;
t498 = qJD(5) * t784 - t623 * t806 + t729 * t740 - t884 * t802;
t504 = -pkin(4) * t802 + t530;
t466 = t498 * pkin(5) + t497 * qJ(6) - qJD(6) * t912 + t504;
t478 = t642 * qJ(6) + t481;
t554 = -pkin(4) * t898 + t590;
t499 = -pkin(5) * t568 - qJ(6) * t912 + t554;
t755 = -t466 * t665 + t478 * t650 - t499 * t908 + t816;
t754 = -t481 * t650 + t504 * t665 + t554 * t908 - t816;
t458 = qJD(6) * t642 + t764 + t872;
t460 = qJDD(6) + t800 - t882;
t477 = -t642 * pkin(5) + t824;
t752 = t458 * t767 + t460 * t665 + t477 * t908 - t478 * t909 - t814;
t751 = -t466 * t767 - t477 * t650 + t499 * t909 + t781;
t750 = t480 * t650 - t504 * t767 + t554 * t909 + t781;
t748 = t499 * t912 + qJDD(6) - t757;
t747 = -g(1) * t632 - g(2) * t630 - t709 * t716 + t764;
t746 = (-t497 * t767 - t498 * t665 + t568 * t908 - t909 * t912) * MDP(23) + (-t497 * t665 + t908 * t912) * MDP(22) + (t574 * t665 + t642 * t908 + t650 * t912) * MDP(24) + (t568 * t650 + t574 * t767 - t642 * t909) * MDP(25) + (t648 * t722 + t741) * MDP(13) + (-t650 * t722 - t577) * MDP(14) + (-t648 ^ 2 + t650 ^ 2) * MDP(12) + t720 * MDP(15) + (-MDP(11) * t648 + MDP(26) * t642) * t650;
t745 = -t540 * t650 + t590 * t859 + t726 * t913;
t744 = -t683 * t648 - t753 + t814;
t689 = t734 * t714;
t687 = qJ(4) * t847;
t686 = qJ(4) * t848;
t680 = t713 - t877;
t643 = -qJDD(1) * t714 + t703;
t597 = t767 * t667;
t596 = t665 * t667;
t592 = -t819 + t604;
t573 = -t635 + t610;
t524 = t616 * t665 + t667 * t891;
t523 = -t616 * t767 + t647 * t667;
t517 = t596 * pkin(5) - t597 * qJ(6) + t589;
t515 = pkin(5) * t912 - qJ(6) * t568;
t496 = pkin(5) * t772 - t773;
t495 = -qJ(6) * t772 + t905;
t479 = -t497 - t915;
t471 = t524 * pkin(5) + t523 * qJ(6) - t597 * qJD(6) + t538;
t465 = -t617 * pkin(5) - t887;
t462 = qJ(6) * t617 - qJD(6) * t772 + t763;
t1 = [(t458 * t495 + t478 * t462 + t466 * t517 + t499 * t471 + t460 * t496 + t477 * t465 - g(1) * (-pkin(5) * t630 - qJ(6) * t629) - g(2) * (pkin(5) * t632 + qJ(6) * t631 + t689) + (-g(1) * t803 - g(2) * t782) * t734 + (-g(1) * (-t714 - t782) - g(2) * t803) * t732) * MDP(32) + (qJDD(1) * t723 + 0.2e1 * t731 * t804) * MDP(4) + (-t667 * t577 - t616 * t648 + t650 * t617 + t741 * t772) * MDP(12) + (t480 * t617 + t589 * t498 + t504 * t596 + t554 * t524 - t538 * t568 + t773 * t574 + t642 * t887 + t772 * t800 + t792) * MDP(27) + (t460 * t772 - t465 * t642 + t466 * t596 - t471 * t568 - t477 * t617 - t496 * t574 + t498 * t517 + t499 * t524 + t792) * MDP(29) + (-t458 * t596 + t460 * t597 + t462 * t568 + t465 * t912 - t477 * t523 - t478 * t524 - t495 * t498 - t496 * t497 - t796) * MDP(30) + (t497 * t596 - t498 * t597 - t523 * t568 - t524 * t912) * MDP(23) + (t498 * t772 - t524 * t642 + t568 * t617 - t574 * t596) * MDP(25) + t790 * MDP(2) + t791 * MDP(3) + (t731 * t777 + t733 * t759) * MDP(9) + (-t731 * t759 + t733 * t777) * MDP(10) + 0.2e1 * (t731 * t821 - t823 * t828) * MDP(5) + (-t458 * t772 + t462 * t642 - t466 * t597 - t471 * t912 + t478 * t617 + t495 * t574 + t497 * t517 + t499 * t523 - t793) * MDP(31) + (t497 * t772 - t523 * t642 + t574 * t597 + t617 * t912) * MDP(24) + (-t481 * t617 - t589 * t497 + t504 * t597 - t554 * t523 + t538 * t912 - t574 * t905 - t642 * t763 + t764 * t772 + t793) * MDP(28) + (-t497 * t597 - t523 * t912) * MDP(22) + (-t574 * t772 + t617 * t642) * MDP(26) + (-t617 * t722 + t720 * t772) * MDP(14) + (-t483 * t857 - t484 * t858 + t505 * t623 + t506 * t898 - t539 * t863 - t540 * t864 - t557 * t740 + t558 * t802 - t796) * MDP(20) + (-t650 * t616 + t667 * t741) * MDP(11) + qJDD(1) * MDP(1) + (qJDD(2) * t731 + t733 * t736) * MDP(6) + (qJDD(2) * t733 - t731 * t736) * MDP(7) + (-t569 * t722 - t683 * t616 - t628 * t720 + t643 * t667 - t650 * t818 - t714 * t741 + t796) * MDP(17) + (-t506 * t648 - t558 * t577 + t484 * t772 - t540 * t617 - t570 * t623 - t890 * t740 + t530 * t857 + t590 * t863 - g(1) * (t726 * t848 + t727 * t734) - g(2) * (-t726 * t847 + t727 * t732)) * MDP(19) + (-t570 * t722 - t577 * t714 - t617 * t683 - t643 * t772 + t648 * t818 + t720 * t890 + t778) * MDP(16) + (-g(2) * t689 + t483 * t557 + t484 * t558 + t539 * t505 + t540 * t506 - t530 * t890 + t590 * t570 + (g(1) * t735 - g(2) * t896) * t734 + (-g(1) * (-t714 - t896) + g(2) * t735) * t732) * MDP(21) + (t505 * t648 + t557 * t577 - t483 * t772 + t539 * t617 - t570 * t898 + t890 * t802 + t727 * t778 + (t530 * t667 + t590 * t616 - t791) * t726) * MDP(18) + (t616 * t722 + t667 * t720) * MDP(13); (-t548 * t623 - t549 * t898 + t697 * t892 + t706 * t894 + t761) * MDP(20) + (t498 * t592 - t568 * t835 - t642 * t836 + t751 + t867) * MDP(29) + (-g(3) * t733 + t731 * t758) * MDP(9) + (t530 * t713 - g(1) * (t734 * t794 + t687) - g(2) * (t732 * t794 + t686) - g(3) * (t896 + t876) + t783 * t706 - t889 * t590 - t904 * t540 + (-t697 * t726 - t548) * t539) * MDP(21) + t746 + (t458 * t602 + t466 * t592 - t460 * t769 - g(3) * (t815 - t852 + t876) + t835 * t499 + t837 * t478 + t836 * t477 + t791 * (t717 * t779 + t849 + t883)) * MDP(32) + (t613 * t722 + (-t648 * t827 + t720 * t885 - t722 * t826) * pkin(2) + t756) * MDP(16) + qJDD(2) * MDP(8) + (t614 * t722 + (t650 * t827 - t720 * t730 - t722 * t807) * pkin(2) + t744) * MDP(17) + (-t577 * t855 + t623 * t889 + t648 * t904 + t713 * t740 + t745) * MDP(19) + (g(3) * t731 + t733 * t758) * MDP(10) + MDP(7) * t821 + (t497 * t769 - t498 * t602 + t568 * t837 + t836 * t912 + t752) * MDP(30) + (-t680 * t497 + t642 * t903 + t795 * t912 + t754 - t866) * MDP(28) + (-t706 * t868 - t713 * t802 + (-t548 + (t590 - t697) * t726) * t648 + t789 + t889 * t898 + t906 * t727) * MDP(18) + (t497 * t592 + t642 * t837 - t835 * t912 + t755 + t866) * MDP(31) + MDP(6) * t822 + (t680 * t498 - t568 * t795 + t642 * t902 + t750 + t867) * MDP(27) + (-t731 * t733 * MDP(4) + MDP(5) * t828) * t737; (t498 * t604 - t568 * t834 - t642 * t832 + t751 + t862) * MDP(29) + (-t707 * t498 + t568 * t573 + t642 * t899 + t750 + t862) * MDP(27) + (t458 * t622 + t466 * t604 - t460 * t768 - g(3) * t815 + t791 * t849 + t834 * t499 + t833 * t478 + t832 * t477 + (g(3) * t728 + t779 * t791) * t717) * MDP(32) + (qJ(4) * t894 + qJD(4) * t892 - t550 * t623 - t551 * t898 + t761) * MDP(20) + t746 + (-t530 * pkin(3) - t540 * t551 - t539 * t550 - t590 * t610 - g(1) * (-pkin(3) * t850 + t687) - g(2) * (-pkin(3) * t851 + t686) - g(3) * t896 + (-t539 * t726 + t540 * t727) * qJD(4) + t783 * qJ(4)) * MDP(21) + (t609 * t722 + t744) * MDP(17) + (-pkin(3) * t740 - qJD(4) * t859 + t551 * t648 - t577 * t873 + t610 * t623 + t745) * MDP(19) + (t497 * t768 - t498 * t622 + t568 * t833 + t832 * t912 + t752) * MDP(30) + (t707 * t497 - t573 * t912 + t642 * t900 + t754 - t861) * MDP(28) + (-qJ(4) * t868 - pkin(3) * t572 + t610 * t637 + (t865 + t878 + t906) * t727 + (-t550 + (-qJD(4) + t590) * t726) * t648 + t789) * MDP(18) + (t497 * t604 + t642 * t833 - t834 * t912 + t755 + t861) * MDP(31) + (t756 + t865) * MDP(16); (-t623 * t648 - t802) * MDP(18) + (t648 * t898 + t740) * MDP(19) + (-t623 ^ 2 - t898 ^ 2) * MDP(20) + (-t539 * t623 - t540 * t898 + t913) * MDP(21) + (-t886 - t916) * MDP(30) + (-t477 * t912 - t478 * t568 + t466 + t888) * MDP(32) + (-MDP(28) + MDP(31)) * (t497 - t915) + (MDP(27) + MDP(29)) * (t498 + t914); -MDP(22) * t869 + (t886 - t916) * MDP(23) + t479 * MDP(24) + (-t498 + t914) * MDP(25) + t574 * MDP(26) + (-t554 * t912 + t757 + t871) * MDP(27) + (t480 * t642 - t554 * t568 - t747) * MDP(28) + (t515 * t568 - t748 + t871 + 0.2e1 * t882) * MDP(29) + (pkin(5) * t497 - qJ(6) * t498 + (t478 - t481) * t912 - (t477 - t824) * t568) * MDP(30) + (0.2e1 * t872 + t499 * t568 + t515 * t912 + (0.2e1 * qJD(6) - t480) * t642 + t747) * MDP(31) + (t458 * qJ(6) - t460 * pkin(5) - t499 * t515 - t477 * t481 - g(1) * (-pkin(5) * t631 + qJ(6) * t632) - g(2) * (-pkin(5) * t629 + qJ(6) * t630) - (-pkin(5) * t715 + qJ(6) * t716) * t709 + t824 * t478) * MDP(32); t479 * MDP(30) + (-t642 ^ 2 - t886) * MDP(31) + (-t478 * t642 + t748 - t882) * MDP(32) + (-t869 - t574) * MDP(29);];
tau  = t1;
