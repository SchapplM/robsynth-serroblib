% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:50:42
% EndTime: 2019-03-09 13:51:14
% DurationCPUTime: 20.49s
% Computational Cost: add. (34155->684), mult. (63300->866), div. (0->0), fcn. (70626->8), ass. (0->404)
t487 = sin(qJ(6));
t490 = cos(qJ(6));
t749 = Ifges(7,4) * t490;
t454 = Ifges(7,1) * t487 + t749;
t654 = t490 * t454;
t750 = Ifges(7,4) * t487;
t453 = Ifges(7,2) * t490 + t750;
t667 = t487 * t453;
t489 = sin(qJ(2));
t761 = sin(qJ(4));
t763 = cos(qJ(4));
t764 = cos(qJ(2));
t438 = t489 * t763 - t761 * t764;
t488 = sin(qJ(5));
t650 = t489 * t761 + t764 * t763;
t762 = cos(qJ(5));
t807 = t488 * t438 + t762 * t650;
t871 = Ifges(6,5) * t807;
t744 = Ifges(7,2) * t487;
t578 = -t744 + t749;
t358 = -t762 * t438 + t488 * t650;
t739 = Ifges(7,6) * t358;
t847 = -t578 * t807 - t739;
t894 = t490 * t847;
t752 = Ifges(7,1) * t490;
t579 = -t750 + t752;
t746 = Ifges(7,5) * t358;
t846 = -t579 * t807 - t746;
t896 = t487 * t846;
t711 = t490 * Ifges(7,6);
t717 = t487 * Ifges(7,5);
t451 = t711 + t717;
t897 = t358 * t451;
t905 = Ifges(6,6) * t358;
t914 = -t897 / 0.4e1 + t905 / 0.2e1 + t894 / 0.4e1 + t896 / 0.4e1;
t644 = t764 * pkin(7);
t547 = -pkin(8) * t764 + t644;
t756 = pkin(7) * t489;
t599 = -t489 * pkin(8) + t756;
t384 = -t761 * t547 + t763 * t599;
t315 = -t438 * pkin(9) + t384;
t808 = t763 * t547 + t599 * t761;
t851 = -t650 * pkin(9) + t808;
t922 = -t762 * t315 + t488 * t851;
t928 = t922 * mrSges(6,2);
t714 = t490 * mrSges(7,1);
t719 = t487 * mrSges(7,2);
t580 = t714 - t719;
t809 = t315 * t488 + t762 * t851;
t934 = t809 * t580;
t936 = t809 * mrSges(6,1);
t917 = t934 / 0.2e1 + t936 / 0.2e1 - t928 / 0.2e1;
t803 = -t917 + t914 - t871 / 0.2e1;
t945 = t803 - (t654 / 0.4e1 - t667 / 0.4e1) * t807;
t474 = t489 * qJ(3);
t823 = -t764 * pkin(2) - t474;
t446 = -pkin(1) + t823;
t411 = t764 * pkin(3) - t446;
t375 = pkin(4) * t650 + t411;
t921 = m(6) * t375 + mrSges(6,1) * t807 - mrSges(6,2) * t358;
t713 = t490 * mrSges(7,2);
t720 = t487 * mrSges(7,1);
t449 = t713 + t720;
t222 = t449 * t358;
t790 = t222 / 0.2e1;
t834 = mrSges(7,2) * t358;
t864 = t487 * t807;
t230 = -mrSges(7,3) * t864 - t834;
t926 = t490 * t230;
t835 = mrSges(7,1) * t358;
t863 = t490 * t807;
t234 = -mrSges(7,3) * t863 + t835;
t927 = t487 * t234;
t819 = t926 / 0.2e1 - t927 / 0.2e1;
t843 = t358 / 0.2e1;
t911 = -t358 / 0.2e1;
t944 = t819 + t790 + (t911 + t843) * mrSges(6,3);
t865 = t449 * t807;
t943 = -t809 * t222 - t922 * t865;
t743 = Ifges(5,6) * t438;
t833 = Ifges(5,5) * t650;
t870 = t808 * mrSges(5,1);
t903 = t384 * mrSges(5,2);
t912 = -t934 - t936 + t905 - t897 / 0.2e1 + t894 / 0.2e1 + t896 / 0.2e1 + t928;
t942 = t912 - t743 - t833 - t870 - t871 - t903;
t797 = m(6) * pkin(4);
t900 = t922 * t488;
t940 = (t762 * t809 + t900) * t797;
t437 = -t488 * t763 - t761 * t762;
t697 = t922 * t437;
t436 = t488 * t761 - t762 * t763;
t700 = t436 * t809;
t939 = -t700 - t697;
t158 = pkin(5) * t807 + pkin(10) * t358 + t375;
t78 = t158 * t490 - t487 * t809;
t788 = t807 / 0.2e1;
t79 = t158 * t487 + t490 * t809;
t731 = t807 * Ifges(7,5);
t132 = -t358 * t579 + t731;
t662 = t490 * t132;
t730 = t807 * Ifges(7,6);
t129 = -t358 * t578 + t730;
t676 = t487 * t129;
t815 = t375 * mrSges(6,2) + Ifges(6,1) * t911 - Ifges(6,4) * t807 - t676 / 0.2e1 + t662 / 0.2e1;
t477 = Ifges(7,6) * t487;
t478 = Ifges(7,5) * t490;
t824 = t478 - t477;
t841 = -t650 / 0.2e1;
t938 = -t79 * t230 - t78 * t234 - (t788 * t824 + t815) * t807 + (t411 * mrSges(5,1) - Ifges(5,4) * t438 + 0.2e1 * (Ifges(5,1) - Ifges(5,2)) * t841) * t438 + (-t411 * mrSges(5,2) - 0.2e1 * Ifges(5,4) * t841) * t650 + t943;
t799 = m(7) / 0.2e1;
t935 = t799 * t809;
t933 = t809 * t922;
t931 = t927 - t926;
t888 = -m(7) * t809 + t865;
t491 = -pkin(2) - pkin(3);
t444 = -qJ(3) * t761 + t763 * t491;
t443 = -pkin(4) + t444;
t445 = qJ(3) * t763 + t491 * t761;
t621 = t762 * t445;
t369 = t488 * t443 + t621;
t902 = t922 * t369;
t901 = t922 * t487;
t899 = t922 * t490;
t728 = t807 * mrSges(6,3);
t484 = t487 ^ 2;
t485 = t490 ^ 2;
t649 = t484 + t485;
t595 = t649 * t436;
t765 = t490 / 0.2e1;
t767 = -t487 / 0.2e1;
t839 = Ifges(6,4) * t843;
t920 = t765 * t846 + t767 * t847 + t839;
t918 = m(5) * t411 + mrSges(5,1) * t650 + t438 * mrSges(5,2);
t916 = t833 / 0.2e1 + t870 / 0.2e1 + t903 / 0.2e1;
t729 = t358 * mrSges(6,3);
t576 = t487 * t78 - t490 * t79;
t892 = (t809 + t576) * t799;
t643 = t762 * pkin(4);
t472 = -t643 - pkin(5);
t768 = t472 / 0.2e1;
t891 = t865 * t768;
t665 = t488 * t445;
t368 = t443 * t762 - t665;
t365 = pkin(5) - t368;
t782 = t365 / 0.2e1;
t890 = t865 * t782;
t371 = t444 * t762 - t665;
t367 = t371 * mrSges(6,2);
t370 = t444 * t488 + t621;
t684 = t370 * t580;
t726 = t370 * mrSges(6,1);
t592 = t649 * t371;
t810 = mrSges(7,3) * t592;
t889 = -t445 * mrSges(5,1) - t444 * mrSges(5,2) - t367 - t684 - t726 + t810;
t712 = t490 * mrSges(7,3);
t233 = mrSges(7,1) * t807 + t358 * t712;
t718 = t487 * mrSges(7,3);
t642 = t358 * t718;
t229 = -mrSges(7,2) * t807 + t642;
t658 = t490 * t229;
t766 = t487 / 0.2e1;
t555 = t233 * t766 - t658 / 0.2e1;
t887 = t555 - t865 / 0.2e1;
t448 = -t764 * mrSges(4,1) - t489 * mrSges(4,3);
t882 = m(4) * t446 + t448;
t241 = -pkin(5) * t358 + pkin(10) * t807;
t881 = m(5) * (-t384 * t761 + t763 * t808);
t876 = -t807 / 0.2e1;
t878 = (t876 + t788) * mrSges(6,3) + t887;
t877 = t449 / 0.2e1;
t872 = -Ifges(4,5) + Ifges(3,4);
t596 = t650 * mrSges(5,3);
t866 = t365 * t809;
t678 = t437 * t580;
t633 = -t478 / 0.2e1;
t559 = t633 + t477 / 0.2e1;
t862 = t559 * t807;
t769 = t454 / 0.4e1;
t858 = t807 * t769;
t759 = pkin(4) * t438;
t192 = t241 + t759;
t92 = t192 * t490 + t901;
t93 = t192 * t487 - t899;
t574 = -t487 * t92 + t490 * t93;
t628 = t712 / 0.2e1;
t631 = -t718 / 0.2e1;
t798 = -pkin(5) / 0.2e1;
t854 = -pkin(5) * t935 - t865 * t798 + t92 * t631 + t93 * t628 + (t574 * t799 - t819) * pkin(10) + t945;
t554 = t667 / 0.2e1 - t654 / 0.2e1;
t475 = t764 * qJ(3);
t433 = t489 * t491 + t475;
t380 = t433 - t759;
t159 = -t241 + t380;
t80 = t159 * t490 - t901;
t81 = t159 * t487 + t899;
t575 = t80 * t487 - t81 * t490;
t827 = t490 * t578;
t828 = t487 * t579;
t825 = -t827 / 0.4e1 - t828 / 0.4e1;
t507 = -t828 / 0.2e1 - t827 / 0.2e1 + t554;
t840 = -t678 / 0.2e1;
t838 = -t763 / 0.2e1;
t305 = t436 * t369;
t651 = t368 * t437 - t305;
t562 = t437 * mrSges(6,1) + t436 * mrSges(6,2) - mrSges(7,3) * t595;
t515 = -mrSges(5,1) * t761 - mrSges(5,2) * t763 + t562;
t822 = t515 + t678;
t219 = t580 * t358;
t758 = pkin(4) * t488;
t471 = pkin(10) + t758;
t634 = t749 / 0.2e1;
t770 = t453 / 0.4e1;
t513 = -t487 * (t634 + t769 - t744 / 0.4e1) - t490 * (t770 - t752 / 0.4e1);
t224 = t453 * t358;
t225 = t454 * t358;
t804 = t490 * t224 / 0.4e1 + t662 / 0.4e1 + t487 * t225 / 0.4e1 - t676 / 0.4e1 + t922 * t877;
t630 = t718 / 0.2e1;
t820 = (t630 + t631) * t79;
t533 = t804 + t820;
t656 = t490 * t233;
t672 = t487 * t229;
t556 = -t672 / 0.2e1 - t656 / 0.2e1;
t583 = mrSges(7,3) * (-t485 / 0.2e1 - t484 / 0.2e1);
t692 = t807 * t824;
t620 = t692 / 0.4e1;
t496 = t556 * t471 - (t471 * t583 + t513) * t358 + t620 - t219 * t768 + t533;
t794 = mrSges(7,2) / 0.2e1;
t795 = -mrSges(7,1) / 0.2e1;
t535 = Ifges(7,3) * t843 + t794 * t93 + t795 * t92;
t16 = t496 + t535 - t862;
t405 = t472 * t449;
t235 = -t405 + t507;
t528 = (Ifges(7,1) / 0.2e1 - Ifges(7,2) / 0.2e1) * t490 - t453 / 0.2e1 - t750 / 0.2e1;
t581 = t454 / 0.2e1 + t634;
t776 = -t371 / 0.2e1;
t783 = -t365 / 0.2e1;
t59 = (t783 + t768) * t449 + (mrSges(7,2) * t776 + t581) * t490 + (mrSges(7,1) * t776 + t528) * t487;
t558 = t713 / 0.2e1 + t720 / 0.2e1;
t257 = (-t449 / 0.2e1 + t558) * t436;
t646 = t257 * qJD(3);
t821 = t16 * qJD(1) - t59 * qJD(2) - t235 * qJD(4) - t646;
t757 = pkin(5) * t449;
t501 = t405 / 0.2e1 - t757 / 0.2e1 - t507;
t625 = t487 * t762;
t584 = -t625 / 0.2e1;
t624 = t490 * t762;
t518 = (mrSges(7,1) * t584 - mrSges(7,2) * t624 / 0.2e1) * pkin(4);
t149 = t518 - t501;
t506 = pkin(10) * t583 + t513;
t510 = pkin(10) * t556 - t219 * t798 + t533;
t98 = t241 * t490 + t901;
t99 = t241 * t487 - t899;
t560 = t794 * t99 + t795 * t98;
t18 = (t824 / 0.4e1 - t559) * t807 - (-Ifges(7,3) / 0.2e1 + t506) * t358 + t510 + t560;
t247 = t507 + t757;
t779 = -t368 / 0.2e1;
t61 = (t783 + t798) * t449 + (mrSges(7,2) * t779 + t581) * t490 + (mrSges(7,1) * t779 + t528) * t487;
t816 = -t18 * qJD(1) + t61 * qJD(2) + t149 * qJD(4) + t247 * qJD(5) + t646;
t573 = -t98 * t487 + t99 * t490;
t106 = m(7) * (-0.1e1 + t649) * t437 * t436;
t647 = t106 * qJD(3);
t227 = t718 * t807 + t834;
t660 = t490 * t227;
t231 = t712 * t807 - t835;
t671 = t487 * t231;
t814 = t660 / 0.2e1 - t671 / 0.2e1;
t813 = mrSges(7,3) * t649;
t366 = -pkin(10) + t369;
t629 = -t712 / 0.2e1;
t698 = t922 * t370;
t777 = t370 / 0.2e1;
t801 = m(6) / 0.2e1;
t805 = t92 * t630 + t93 * t629 + (-t902 + t698 + (-t368 + t371) * t809) * t801 + (t366 * t574 - t371 * t576 + t698 + t866) * t799 - t890 - t222 * t777 - t803 + (t369 * t843 - t358 * t777 + (t776 - t779) * t807) * mrSges(6,3) + t916;
t800 = -m(7) / 0.2e1;
t796 = m(7) * pkin(4);
t793 = Ifges(7,3) / 0.2e1;
t792 = t98 / 0.2e1;
t791 = -t99 / 0.2e1;
t781 = -t366 / 0.2e1;
t780 = t366 / 0.2e1;
t778 = t368 / 0.2e1;
t772 = t437 / 0.2e1;
t760 = m(7) * (pkin(5) * t437 - pkin(10) * t595);
t20 = (t878 + t892) * t436 + ((-t922 - t574) * t799 + t944) * t437;
t22 = (t790 + (-t922 - t573) * t799 - t814) * t437 + (t887 + t892) * t436;
t754 = t20 * qJD(4) + t22 * qJD(5);
t561 = t375 * mrSges(6,1) + t824 * t911 + t839 + (Ifges(6,2) + Ifges(7,3)) * t788;
t639 = Ifges(6,2) / 0.2e1 + t793;
t524 = t639 * t807 + t561;
t3 = (-pkin(1) * mrSges(3,1) + t446 * mrSges(4,1) - t872 * t489) * t489 + t81 * t229 + t80 * t233 + m(7) * (t78 * t80 + t79 * t81 - t933) + t882 * (pkin(2) * t489 - t475) + ((Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t489 - mrSges(4,3) * t446 - pkin(1) * mrSges(3,2) + t872 * t764) * t764 + t918 * t433 + (-Ifges(6,1) * t788 - Ifges(6,4) * t911 - (Ifges(7,1) * t863 + t746) * t765 - (-Ifges(7,2) * t864 + t739) * t767 + t863 * t750 + t524) * t358 - t938 + t921 * t380;
t732 = t3 * qJD(1);
t727 = t369 * mrSges(6,1);
t33 = (t656 + m(7) * (t487 * t79 + t490 * t78) + t672 - t882 + t918 + t921) * t489;
t705 = qJD(1) * t33;
t632 = t730 / 0.2e1;
t10 = t79 * t233 + t922 * t219 - ((t132 / 0.2e1 + t224 / 0.2e1 + t731 / 0.2e1) * t487 + (-t225 / 0.2e1 + t129 / 0.2e1 + t632 + t79 * mrSges(7,3)) * t490) * t358 + (t642 - t229) * t78;
t704 = t10 * qJD(1);
t687 = t365 * t449;
t686 = t369 * t580;
t685 = t370 * t436;
t683 = t436 * t219;
t682 = t436 * t449;
t679 = t436 * t488;
t544 = t558 * t436;
t259 = t682 / 0.2e1 + t544;
t653 = t259 * qJD(6) + t647;
t652 = -t257 * qJD(6) - t647;
t645 = qJD(4) + qJD(5);
t641 = -t758 / 0.2e1;
t640 = t758 / 0.2e1;
t637 = t764 * mrSges(4,2);
t627 = Ifges(6,6) / 0.2e1 - t451 / 0.4e1;
t626 = -t365 * t437 - t366 * t595;
t616 = t371 * t229 / 0.2e1;
t615 = t233 * t776;
t613 = t437 * t766;
t610 = t864 / 0.2e1;
t607 = t667 / 0.4e1;
t604 = -t863 / 0.2e1;
t602 = -t654 / 0.4e1;
t594 = t649 * t437;
t593 = t649 * t368;
t591 = t649 * t471;
t590 = mrSges(6,3) * t643;
t589 = mrSges(7,3) * t643;
t586 = t643 / 0.2e1;
t495 = -t472 * t935 + t891 + t940 / 0.2e1 + t641 * t729 + t590 * t876 + (-t575 * t799 + t819) * t471 + t916;
t2 = -t366 * t819 + t487 * t615 + t490 * t616 + t81 * t629 + t80 * t630 - t495 + t803 + t805;
t37 = -m(7) * (t365 * t370 + t366 * t592) - m(6) * (-t368 * t370 + t369 * t371) + t889;
t572 = t2 * qJD(1) - t37 * qJD(2);
t514 = (pkin(5) * t809 - pkin(10) * t575) * t800 + pkin(5) * t865 / 0.2e1;
t519 = -t627 * t358 - t890 - t369 * t222 / 0.2e1;
t538 = -t846 / 0.4e1 + t231 * t781 + t233 * t779;
t539 = -t847 / 0.4e1 + t227 * t780 + t229 * t778;
t565 = t866 + t902;
t4 = (-pkin(10) * t230 / 0.2e1 + t539) * t490 + (pkin(10) * t234 / 0.2e1 + t538) * t487 + (t366 * t573 - t368 * t576 + t565) * t799 + (-t711 / 0.4e1 - t717 / 0.4e1 + t627) * t358 + t825 * t807 + ((t791 - t81 / 0.2e1) * t490 + (t80 / 0.2e1 + t792) * t487) * mrSges(7,3) + t514 + t519;
t364 = t368 * mrSges(6,2);
t520 = mrSges(7,3) * t593 - t364 - t686 - t727;
t41 = -m(7) * (t365 * t369 + t366 * t593) + t520;
t571 = t4 * qJD(1) - t41 * qJD(2);
t6 = t93 * t229 + t92 * t233 + m(7) * (t78 * t92 + t79 * t93 + t933) - (Ifges(6,1) * t876 + t524 + t920) * t358 + t921 * t759 + t938;
t570 = t6 * qJD(1) + t20 * qJD(3);
t9 = t98 * t233 + t78 * t231 + m(7) * (t78 * t98 + t79 * t99 + t933) + t99 * t229 + t79 * t227 - (t561 + t920) * t358 + (t862 - (-Ifges(6,1) / 0.2e1 + t639) * t358 - t815) * t807 + t943;
t569 = t9 * qJD(1) + t22 * qJD(3);
t499 = -t881 / 0.2e1 - m(6) * t939 / 0.2e1 + (t437 * t575 - t700) * t800;
t500 = t881 / 0.2e1 + t939 * t801 + (t576 * t436 - t697) * t799;
t14 = t878 * t436 + t437 * t944 + t499 + t500;
t53 = -mrSges(4,3) - m(4) * qJ(3) - m(7) * t626 - m(6) * t651 - m(5) * (-t444 * t761 + t445 * t763) + t822;
t568 = -t14 * qJD(1) + t53 * qJD(2);
t503 = t840 + (-t371 * t437 + t651 + t685) * t801 + (-t371 * t594 + t626 + t685) * t799;
t542 = -t472 * t437 - t471 * t595;
t532 = m(7) * t542;
t536 = (t437 * t762 - t679) * t797;
t511 = t532 / 0.2e1 + t536 / 0.2e1;
t34 = t840 + t503 - t511 - t515;
t567 = t20 * qJD(1) + t34 * qJD(2);
t504 = t840 + (-t368 * t594 + t305 + t626) * t799;
t36 = t840 - t760 / 0.2e1 + t504 - t562;
t566 = t22 * qJD(1) + t36 * qJD(2);
t497 = t556 * t366 - (t366 * t583 - t513) * t358 - t692 / 0.4e1 - t219 * t782 - t804 + t820;
t534 = Ifges(7,3) * t911 + t794 * t81 + t795 * t80;
t12 = t497 + t534 + t862;
t160 = t687 + t507;
t564 = t12 * qJD(1) - t160 * qJD(2);
t258 = (t877 + t558) * t436;
t545 = (t714 / 0.2e1 - t719 / 0.2e1) * t489;
t30 = t683 / 0.2e1 + t545 + (-t358 * t583 + t556) * t437;
t563 = t30 * qJD(1) + t258 * qJD(2);
t548 = t678 + t562;
t546 = t649 * t762;
t531 = t687 / 0.2e1 + 0.2e1 * t607 + 0.2e1 * t602 + 0.2e1 * t825;
t461 = t484 * t589;
t462 = t485 * t589;
t525 = -mrSges(6,2) * t643 - t758 * (t580 + mrSges(6,1)) + t461 + t462;
t202 = -(t471 * t546 + t472 * t488) * t796 - t525;
t494 = -t364 / 0.2e1 - t461 / 0.2e1 - t462 / 0.2e1 + (t472 * t369 + t368 * t591 + (t365 * t488 + t366 * t546) * pkin(4)) * t799 - t727 / 0.2e1 - t686 / 0.2e1 + mrSges(6,1) * t640 - t580 * t641 + mrSges(6,2) * t586 + t778 * t813;
t502 = t367 / 0.2e1 + (-pkin(5) * t370 + pkin(10) * t592) * t800 + t726 / 0.2e1 + t684 / 0.2e1 - t810 / 0.2e1;
t29 = t494 + t502;
t261 = t760 / 0.2e1;
t505 = m(7) * ((-t437 * t546 + t679) * pkin(4) + t542);
t65 = t261 - t505 / 0.2e1;
t516 = Ifges(6,5) * t788 + t917;
t492 = (t472 * t809 + (t624 * t79 - t625 * t78 + t900) * pkin(4)) * t799 - t891 - t222 * t640 + t98 * t631 + t99 * t628 + pkin(4) * t233 * t584 + t586 * t658 - t516 + (t573 * t799 + t814) * t471 + (t607 + t602) * t807 + t914;
t8 = -t492 + t854;
t526 = -t8 * qJD(1) + t29 * qJD(2) - t65 * qJD(3) - t202 * qJD(4);
t517 = t81 * t628 + t80 * t631 - t945;
t389 = t678 / 0.2e1;
t260 = -t682 / 0.2e1 + t544;
t150 = t518 + t501;
t62 = t757 / 0.2e1 - t558 * t368 + t531;
t60 = -t405 / 0.2e1 - t558 * t371 + t531;
t54 = t505 / 0.2e1 + t261 + t548;
t38 = t261 + t389 + t504;
t35 = t389 + t503 + t511;
t31 = -t683 / 0.2e1 + t229 * t613 + t545 + (-t358 * t813 + t656) * t772;
t28 = t494 - t502;
t17 = t807 * t633 + t487 * t632 + t510 - t560 + t620 - (t506 + t793) * t358;
t15 = Ifges(7,5) * t604 + Ifges(7,6) * t610 + t496 - t535;
t13 = m(4) * t644 + t637 + t234 * t613 + t729 * t772 + t596 * t838 + (t761 * t438 + t650 * t838) * mrSges(5,3) + (mrSges(6,3) * t788 + t555) * t436 - t499 + t500 + (t865 + t728) * t436 / 0.2e1 + (-t926 / 0.2e1 + mrSges(6,3) * t843 + t790) * t437;
t11 = t497 + Ifges(7,5) * t863 / 0.2e1 - Ifges(7,6) * t864 / 0.2e1 - t534;
t7 = t492 + t854;
t5 = -t514 + t516 + t517 + (mrSges(7,3) * t791 + (t366 * t99 + t368 * t79) * t799 + t858 + t539) * t490 + ((-t366 * t98 - t368 * t78) * t799 - t807 * t453 / 0.4e1 + mrSges(7,3) * t792 + t538) * t487 + t565 * t799 + t519 + t819 * pkin(10);
t1 = t495 + (-t230 * t780 + t616 + t858) * t490 + (-t234 * t781 - t770 * t807 + t615) * t487 + t743 + t517 + t805;
t19 = [qJD(2) * t3 + qJD(3) * t33 + qJD(4) * t6 + qJD(5) * t9 - qJD(6) * t10, t13 * qJD(3) + t1 * qJD(4) + t5 * qJD(5) + t11 * qJD(6) + t732 + (t575 * mrSges(7,3) + t554 * t807 + (Ifges(3,5) + Ifges(4,4)) * t764 + m(5) * (-t384 * t445 + t444 * t808) + t445 * t438 * mrSges(5,3) + mrSges(3,2) * t756 - t368 * t728 - t369 * t729 - mrSges(4,2) * t474 - mrSges(3,1) * t644 - pkin(2) * t637 + (m(4) * t823 + t448) * pkin(7) - t444 * t596 + (-m(7) * t575 - t931) * t366 + t888 * t365 + m(6) * (t368 * t809 + t902) + (-Ifges(3,6) + Ifges(4,6)) * t489 + t942) * qJD(2), qJD(2) * t13 + qJD(6) * t31 + t705 + t754, t1 * qJD(2) + (t453 * t610 + t454 * t604 - t940 + t758 * t729 + t807 * t590 - t888 * t472 + (m(7) * t574 + t931) * t471 + t574 * mrSges(7,3) + t942) * qJD(4) + t7 * qJD(5) + t15 * qJD(6) + t570, t5 * qJD(2) + t7 * qJD(4) + t17 * qJD(6) + t569 + ((-Ifges(6,5) + t554) * t807 + t888 * pkin(5) + (m(7) * t573 + t660 - t671) * pkin(10) + t573 * mrSges(7,3) + t912) * qJD(5), -t704 + t11 * qJD(2) + t31 * qJD(3) + t15 * qJD(4) + t17 * qJD(5) + (-mrSges(7,1) * t79 - mrSges(7,2) * t78 + t897) * qJD(6); qJD(3) * t14 + qJD(4) * t2 + qJD(5) * t4 + qJD(6) * t12 - t732, -qJD(3) * t53 - qJD(4) * t37 - qJD(5) * t41 - qJD(6) * t160, t35 * qJD(4) + t38 * qJD(5) + t260 * qJD(6) - t568 + t647, t35 * qJD(3) + (m(7) * (t370 * t472 + t371 * t591) + (-t370 * t762 + t371 * t488) * t797 + t889) * qJD(4) + t28 * qJD(5) + t60 * qJD(6) + t572, t38 * qJD(3) + t28 * qJD(4) + (m(7) * (-pkin(5) * t369 + pkin(10) * t593) + t520) * qJD(5) + t62 * qJD(6) + t571, t260 * qJD(3) + t60 * qJD(4) + t62 * qJD(5) + (-t366 * t580 - t824) * qJD(6) + t564; -qJD(2) * t14 - qJD(6) * t30 - t705 + t754, qJD(4) * t34 + qJD(5) * t36 - qJD(6) * t258 + t568, t645 * t106 (t532 + t536 + t822) * qJD(4) + t54 * qJD(5) + t567 + t653, t54 * qJD(4) + (t548 + t760) * qJD(5) + t566 + t653, qJD(6) * t678 + t259 * t645 - t563; -qJD(2) * t2 - qJD(5) * t8 + qJD(6) * t16 - t570, -qJD(3) * t34 + qJD(5) * t29 - qJD(6) * t59 - t572, -qJD(5) * t65 - t567 + t652, -qJD(5) * t202 - qJD(6) * t235 ((-pkin(5) * t488 + pkin(10) * t546) * t796 + t525) * qJD(5) + t150 * qJD(6) + t526, t150 * qJD(5) + (-t471 * t580 + t824) * qJD(6) + t821; -qJD(2) * t4 + qJD(4) * t8 + qJD(6) * t18 - t569, -qJD(3) * t36 - qJD(4) * t29 - qJD(6) * t61 - t571, qJD(4) * t65 - t566 + t652, -qJD(6) * t149 - t526, -t247 * qJD(6) (-pkin(10) * t580 + t824) * qJD(6) - t816; -qJD(2) * t12 + qJD(3) * t30 - qJD(4) * t16 - qJD(5) * t18 + t704, qJD(3) * t258 + qJD(4) * t59 + qJD(5) * t61 - t564, t257 * t645 + t563, qJD(5) * t149 - t821, t816, 0;];
Cq  = t19;
