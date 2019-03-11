% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:00:36
% EndTime: 2019-03-09 07:01:15
% DurationCPUTime: 25.62s
% Computational Cost: add. (42258->862), mult. (86209->1140), div. (0->0), fcn. (90827->10), ass. (0->452)
t888 = qJD(4) + qJD(5);
t547 = sin(qJ(5));
t548 = sin(qJ(4));
t551 = cos(qJ(5));
t552 = cos(qJ(4));
t666 = t551 * t552;
t494 = -t547 * t548 + t666;
t495 = -t547 * t552 - t551 * t548;
t546 = sin(qJ(6));
t550 = cos(qJ(6));
t394 = t494 * t546 - t495 * t550;
t624 = t550 * t494 + t495 * t546;
t906 = Ifges(7,5) * t624 - Ifges(7,6) * t394;
t837 = -pkin(9) - pkin(8);
t517 = t837 * t548;
t518 = t837 * t552;
t419 = t517 * t547 - t551 * t518;
t361 = pkin(10) * t494 + t419;
t866 = t551 * t517 + t547 * t518;
t889 = t495 * pkin(10) + t866;
t935 = -t361 * t546 + t550 * t889;
t943 = t935 * mrSges(7,2);
t216 = t361 * t550 + t546 * t889;
t946 = t216 * mrSges(7,1);
t950 = -t946 / 0.2e1 - t943 / 0.2e1;
t955 = 0.2e1 * t950 + t906;
t957 = t955 * qJD(6);
t762 = Ifges(7,4) * t394;
t245 = Ifges(7,2) * t624 + t762;
t402 = -mrSges(6,1) * t495 + mrSges(6,2) * t494;
t489 = Ifges(6,4) * t494;
t405 = Ifges(6,2) * t495 + t489;
t408 = -Ifges(6,1) * t495 + t489;
t536 = -pkin(4) * t552 - pkin(3);
t817 = -t394 / 0.2e1;
t878 = t624 / 0.2e1;
t899 = t394 / 0.2e1;
t386 = Ifges(7,4) * t624;
t248 = Ifges(7,1) * t394 + t386;
t880 = -Ifges(7,2) * t394 + t386;
t905 = t880 + t248;
t908 = Ifges(7,1) * t624 - t762;
t436 = -pkin(5) * t494 + t536;
t882 = mrSges(7,1) * t394 + mrSges(7,2) * t624;
t913 = t436 * t882;
t956 = t913 + t536 * t402 + (t408 / 0.2e1 + t405 / 0.2e1) * t494 + t245 * t817 + t908 * t899 + t905 * t878;
t553 = cos(qJ(3));
t549 = sin(qJ(3));
t673 = t548 * t549;
t462 = t547 * t673 - t549 * t666;
t463 = t495 * t549;
t625 = t462 * t546 + t550 * t463;
t735 = t625 * mrSges(7,3);
t285 = mrSges(7,2) * t553 + t735;
t832 = t285 / 0.2e1;
t343 = t462 * t550 - t463 * t546;
t739 = t343 * mrSges(7,3);
t287 = -mrSges(7,1) * t553 + t739;
t834 = -t216 / 0.2e1;
t930 = t287 * t834;
t954 = t935 * t832 + t930;
t241 = -mrSges(7,1) * t624 + mrSges(7,2) * t394;
t953 = m(7) * t436 + t241;
t404 = Ifges(6,5) * t494 + Ifges(6,6) * t495;
t620 = t404 + t906;
t891 = t866 * mrSges(6,2);
t918 = t419 * mrSges(6,1);
t936 = -t943 - t946;
t952 = t620 - t891 - t918 + t936;
t835 = -t935 / 0.2e1;
t944 = (t835 + t935 / 0.2e1) * mrSges(7,2) + (t834 + t216 / 0.2e1) * mrSges(7,1);
t951 = qJD(6) * t944;
t892 = t343 * mrSges(7,1);
t641 = t892 / 0.2e1;
t826 = -t625 / 0.2e1;
t830 = -t892 / 0.2e1;
t877 = t625 / 0.2e1;
t579 = t641 + t830 + (t826 + t877) * mrSges(7,2);
t949 = qJD(2) * t579 + qJD(3) * t944;
t870 = t624 * mrSges(7,3);
t945 = -t216 * t550 + t546 * t935;
t843 = m(7) * pkin(5);
t920 = -t843 / 0.2e1;
t938 = t882 / 0.2e1;
t763 = Ifges(7,4) * t343;
t188 = Ifges(7,2) * t625 - t553 * Ifges(7,6) - t763;
t532 = sin(pkin(11)) * pkin(1) + pkin(7);
t513 = t549 * t532;
t472 = pkin(4) * t673 + t513;
t372 = -pkin(5) * t463 + t472;
t869 = t625 * mrSges(7,2);
t904 = -t892 + t869;
t907 = Ifges(7,1) * t625 + t763;
t934 = -(-t907 / 0.4e1 + t188 / 0.4e1) * t394 + t372 * t938 + t436 * t904 / 0.2e1;
t633 = -cos(pkin(11)) * pkin(1) - pkin(2);
t484 = -pkin(3) * t553 - t549 * pkin(8) + t633;
t466 = t552 * t484;
t671 = t549 * t552;
t648 = pkin(9) * t671;
t349 = -t648 + t466 + (-t532 * t548 - pkin(4)) * t553;
t514 = t553 * t532;
t397 = t548 * t484 + t514 * t552;
t371 = -pkin(9) * t673 + t397;
t350 = t547 * t371;
t217 = t551 * t349 - t350;
t454 = t462 * pkin(10);
t155 = t217 + t454;
t782 = pkin(5) * t553;
t151 = -t782 + t155;
t352 = t551 * t371;
t218 = t349 * t547 + t352;
t780 = pkin(10) * t463;
t156 = t218 + t780;
t706 = t156 * t550;
t88 = t151 * t546 + t706;
t91 = -t155 * t546 - t706;
t932 = t88 + t91;
t929 = t372 * t904;
t396 = -t514 * t548 + t466;
t370 = t396 - t648;
t222 = t551 * t370 - t350;
t928 = t217 - t222;
t464 = t495 * t553;
t465 = t494 * t553;
t345 = t464 * t550 - t465 * t546;
t348 = t464 * t546 + t465 * t550;
t722 = t465 * mrSges(6,2);
t723 = t464 * mrSges(6,1);
t598 = -t723 / 0.2e1 + t722 / 0.2e1;
t733 = t348 * mrSges(7,2);
t737 = t345 * mrSges(7,1);
t659 = t737 / 0.2e1 - t733 / 0.2e1;
t622 = t659 - t598;
t926 = (t345 * t550 + t348 * t546) * t920 - t622;
t925 = t908 / 0.4e1 - t245 / 0.4e1;
t599 = -t891 / 0.2e1 - t918 / 0.2e1;
t724 = t463 * mrSges(6,3);
t615 = mrSges(6,2) * t553 + t724;
t725 = t462 * mrSges(6,3);
t616 = -mrSges(6,1) * t553 + t725;
t718 = t495 * mrSges(6,3);
t635 = t718 / 0.2e1;
t638 = -t870 / 0.2e1;
t640 = mrSges(7,3) * t817;
t642 = -t735 / 0.2e1;
t643 = t739 / 0.2e1;
t788 = -t553 / 0.4e1;
t808 = t419 / 0.2e1;
t809 = t866 / 0.2e1;
t810 = -t419 / 0.2e1;
t707 = t156 * t546;
t87 = t151 * t550 - t707;
t442 = Ifges(6,4) * t463;
t337 = -Ifges(6,1) * t462 - t553 * Ifges(6,5) + t442;
t355 = Ifges(6,2) * t462 + t442;
t890 = t355 + t337;
t895 = -t866 / 0.2e1;
t324 = Ifges(7,4) * t625;
t190 = -Ifges(7,1) * t343 - t553 * Ifges(7,5) + t324;
t765 = Ifges(6,4) * t462;
t335 = Ifges(6,2) * t463 - t553 * Ifges(6,6) - t765;
t356 = Ifges(6,1) * t463 + t765;
t764 = Ifges(6,4) * t495;
t406 = Ifges(6,2) * t494 - t764;
t407 = Ifges(6,1) * t494 + t764;
t439 = t462 * mrSges(6,1);
t627 = t463 * mrSges(6,2) - t439;
t881 = Ifges(7,2) * t343 + t324;
t902 = (t407 / 0.4e1 - t406 / 0.4e1) * t462 - (t335 / 0.4e1 - t356 / 0.4e1) * t495 - t536 * t627 / 0.2e1 - t472 * t402 / 0.2e1 - (t408 + t405) * t463 / 0.4e1 - (t190 + t881) * t624 / 0.4e1 - t905 * t625 / 0.4e1 + t925 * t343 - t934;
t923 = t616 * t810 + t615 * t809 + t87 * t638 + t88 * t640 + t935 * t642 + t216 * t643 + t724 * t895 + t725 * t808 + t218 * t635 + t890 * t494 / 0.4e1 + t620 * t788 - t902;
t849 = -m(7) / 0.2e1;
t921 = -Ifges(7,3) / 0.2e1;
t900 = -t343 / 0.2e1;
t916 = t287 * t900;
t221 = -t370 * t547 - t352;
t910 = t218 + t221;
t601 = -Ifges(7,5) * t348 / 0.2e1 - Ifges(7,6) * t345 / 0.2e1;
t777 = t549 * pkin(3);
t781 = pkin(8) * t553;
t519 = t777 - t781;
t420 = t548 * t513 + t552 * t519;
t664 = t552 * t553;
t376 = t549 * pkin(4) - pkin(9) * t664 + t420;
t421 = t548 * t519 - t532 * t671;
t672 = t548 * t553;
t395 = -pkin(9) * t672 + t421;
t232 = t551 * t376 - t395 * t547;
t170 = pkin(5) * t549 - pkin(10) * t465 + t232;
t233 = t547 * t376 + t551 * t395;
t174 = pkin(10) * t464 + t233;
t100 = t170 * t546 + t174 * t550;
t99 = t170 * t550 - t174 * t546;
t864 = t100 * mrSges(7,2) / 0.2e1 - t99 * mrSges(7,1) / 0.2e1;
t617 = -t549 * t921 - t601 - t864;
t871 = Ifges(7,5) * t625;
t894 = Ifges(7,6) * t343;
t658 = t871 + t894;
t630 = t894 / 0.2e1 + t871 / 0.2e1;
t603 = -t420 * t548 + t421 * t552;
t785 = pkin(4) * t551;
t535 = pkin(5) + t785;
t677 = t546 * t547;
t478 = -pkin(4) * t677 + t535 * t550;
t675 = t547 * t550;
t479 = pkin(4) * t675 + t535 * t546;
t485 = (-t546 * t551 - t675) * pkin(4);
t486 = (t550 * t551 - t677) * pkin(4);
t887 = ((t479 + t485) * t935 + (-t478 + t486) * t216) * t849 - t950;
t165 = t221 - t780;
t166 = t454 + t222;
t96 = t165 * t546 + t166 * t550;
t771 = t96 * mrSges(7,2);
t95 = t165 * t550 - t166 * t546;
t772 = t95 * mrSges(7,1);
t768 = t772 / 0.2e1 - t771 / 0.2e1;
t886 = (t546 * t96 + t550 * t95) * t920 - t768;
t543 = t548 ^ 2;
t544 = t552 ^ 2;
t875 = mrSges(5,3) * (t544 / 0.2e1 + t543 / 0.2e1);
t540 = Ifges(5,5) * t552;
t758 = Ifges(5,6) * t548;
t868 = Ifges(4,4) - t540 / 0.2e1 + t758 / 0.2e1;
t541 = Ifges(5,4) * t552;
t865 = -Ifges(5,2) * t548 + t541;
t511 = Ifges(5,1) * t548 + t541;
t621 = Ifges(6,5) * t463 + Ifges(6,6) * t462 + t658;
t715 = t548 * mrSges(5,2);
t863 = -t552 * mrSges(5,1) + t715;
t740 = t233 * mrSges(6,2);
t741 = t232 * mrSges(6,1);
t859 = -t740 / 0.2e1 + t741 / 0.2e1;
t858 = -Ifges(6,6) * t464 / 0.2e1 - Ifges(6,5) * t465 / 0.2e1;
t857 = qJD(6) * t579;
t582 = 0.2e1 * t641 - t869;
t856 = t582 * qJD(6);
t775 = t88 * mrSges(7,1);
t776 = t87 * mrSges(7,2);
t855 = -t775 / 0.2e1 - t776 / 0.2e1 + t630;
t853 = 0.2e1 * qJD(3);
t852 = m(5) / 0.2e1;
t851 = -m(6) / 0.2e1;
t850 = m(6) / 0.2e1;
t848 = m(7) / 0.2e1;
t845 = -pkin(8) / 0.2e1;
t844 = m(6) * pkin(4);
t842 = -mrSges(5,1) / 0.2e1;
t841 = mrSges(5,2) / 0.2e1;
t840 = -t87 / 0.2e1;
t839 = -t88 / 0.2e1;
t838 = -t95 / 0.2e1;
t193 = -mrSges(7,1) * t625 - mrSges(7,2) * t343;
t836 = t193 / 0.2e1;
t833 = -t285 / 0.2e1;
t824 = t345 / 0.2e1;
t822 = t343 / 0.2e1;
t820 = t348 / 0.2e1;
t353 = -mrSges(6,1) * t463 - mrSges(6,2) * t462;
t819 = t353 / 0.2e1;
t816 = -t624 / 0.2e1;
t812 = t406 / 0.2e1;
t784 = pkin(5) * t462;
t415 = pkin(4) * t671 - t784;
t811 = t415 / 0.2e1;
t807 = -t462 / 0.2e1;
t806 = t463 / 0.2e1;
t805 = t464 / 0.2e1;
t804 = t465 / 0.2e1;
t803 = -t478 / 0.2e1;
t802 = -t479 / 0.2e1;
t801 = t485 / 0.2e1;
t800 = t486 / 0.2e1;
t798 = -t495 / 0.2e1;
t712 = t552 * mrSges(5,2);
t716 = t548 * mrSges(5,1);
t508 = t712 + t716;
t797 = t508 / 0.2e1;
t796 = -t546 / 0.2e1;
t795 = -t548 / 0.2e1;
t794 = t548 / 0.2e1;
t793 = t549 / 0.2e1;
t792 = -t550 / 0.2e1;
t791 = t550 / 0.2e1;
t790 = t552 / 0.2e1;
t789 = -t553 / 0.2e1;
t787 = t553 / 0.2e1;
t786 = pkin(4) * t548;
t783 = pkin(5) * t495;
t778 = t546 * pkin(5);
t774 = t91 * mrSges(7,1);
t92 = t155 * t550 - t707;
t773 = t92 * mrSges(7,2);
t769 = t774 / 0.2e1 - t773 / 0.2e1;
t766 = Ifges(5,4) * t548;
t754 = pkin(4) * qJD(4);
t753 = pkin(5) * qJD(5);
t745 = t217 * mrSges(6,2);
t744 = t218 * mrSges(6,1);
t743 = t221 * mrSges(6,1);
t742 = t222 * mrSges(6,2);
t730 = t394 * mrSges(7,3);
t721 = t485 * mrSges(7,1);
t720 = t486 * mrSges(7,2);
t719 = t494 * mrSges(6,3);
t717 = t547 * mrSges(6,1);
t714 = t551 * mrSges(6,2);
t711 = t553 * Ifges(5,5);
t710 = t553 * Ifges(5,6);
t709 = -mrSges(4,1) + t863;
t501 = mrSges(5,2) * t553 - mrSges(5,3) * t673;
t503 = -mrSges(5,1) * t553 - mrSges(5,3) * t671;
t568 = (t343 * t822 + t625 * t877) * mrSges(7,3) + t916 + t625 * t833;
t595 = t830 - t439 / 0.2e1 + t869 / 0.2e1;
t613 = t87 * t343 + t625 * t88;
t649 = t844 / 0.2e1;
t15 = (-t343 * t96 + t625 * t95 + t613) * t849 + (t462 * t928 + t910 * t463) * t851 + (t501 * t794 + t503 * t790 + t549 * t875) * t549 + (m(7) * t811 + t439 / 0.2e1 + (-t715 / 0.2e1 + (t649 + mrSges(5,1) / 0.2e1) * t552) * t549 + t595) * t553 + t568;
t708 = t15 * qJD(1);
t18 = (-t343 * t92 + t625 * t91 + t613) * t849 + ((t920 + mrSges(6,1) / 0.2e1) * t462 + t595) * t553 + t568;
t705 = t18 * qJD(1);
t26 = t904 * t787 + t916 + t285 * t826 + (t625 ^ 2 / 0.2e1 + t343 ^ 2 / 0.2e1) * mrSges(7,3);
t703 = t26 * qJD(1);
t689 = t462 * t241;
t685 = t478 * t625;
t684 = t478 * t624;
t683 = t479 * t343;
t682 = t479 * t394;
t681 = t495 * t193;
t679 = t546 * t343;
t678 = t546 * t394;
t422 = -mrSges(6,2) * t549 + t464 * mrSges(6,3);
t676 = t547 * t422;
t504 = t549 * mrSges(5,1) - mrSges(5,3) * t664;
t674 = t548 * t504;
t670 = t549 * t553;
t669 = t550 * t625;
t668 = t550 * t624;
t423 = mrSges(6,1) * t549 - t465 * mrSges(6,3);
t667 = t551 * t423;
t502 = -t549 * mrSges(5,2) - mrSges(5,3) * t672;
t665 = t552 * t502;
t660 = t478 * t343 + t479 * t625;
t529 = pkin(4) * t672;
t473 = t514 + t529;
t652 = t543 + t544;
t651 = qJD(3) * t553;
t150 = -t479 * mrSges(7,1) - t478 * mrSges(7,2);
t650 = qJD(6) * t150;
t639 = t870 / 0.2e1;
t637 = -t730 / 0.2e1;
t636 = -t718 / 0.2e1;
t634 = t711 / 0.2e1;
t632 = t882 * t789;
t631 = -t217 / 0.2e1 + t222 / 0.2e1;
t623 = t87 * mrSges(7,3) - t190 / 0.2e1;
t619 = -t778 / 0.2e1 + t802;
t618 = pkin(5) * t792 + t803;
t512 = Ifges(5,1) * t552 - t766;
t509 = Ifges(5,2) * t552 + t766;
t194 = t733 - t737;
t286 = -mrSges(7,2) * t549 + mrSges(7,3) * t345;
t288 = mrSges(7,1) * t549 - mrSges(7,3) * t348;
t354 = t722 - t723;
t373 = -pkin(5) * t464 + t473;
t481 = t553 * t508;
t597 = -t712 / 0.2e1 - t716 / 0.2e1;
t14 = t422 * t807 + t287 * t824 + t288 * t877 + t285 * t820 + t286 * t900 + t423 * t806 + (t462 * t805 + t463 * t804) * mrSges(6,3) + (t501 * t790 + t503 * t795 - t194 / 0.2e1 - t354 / 0.2e1 - t481 / 0.2e1 + t598) * t553 + (t665 / 0.2e1 - t674 / 0.2e1 + t836 + t819 - t597 * t549) * t549 + (-t100 * t343 + t87 * t345 + t88 * t348 + t372 * t549 - t373 * t553 + t625 * t99) * t848 + (t217 * t464 + t218 * t465 + t232 * t463 - t233 * t462 + t472 * t549 - t473 * t553) * t850 + ((-t396 * t548 + t397 * t552 - t514) * t553 + (t603 + t513) * t549) * t852;
t189 = Ifges(7,4) * t348 + Ifges(7,2) * t345 + t549 * Ifges(7,6);
t191 = Ifges(7,1) * t348 + Ifges(7,4) * t345 + t549 * Ifges(7,5);
t336 = Ifges(6,4) * t465 + Ifges(6,2) * t464 + t549 * Ifges(6,6);
t338 = Ifges(6,1) * t465 + Ifges(6,4) * t464 + t549 * Ifges(6,5);
t457 = t549 * t865 - t710;
t458 = Ifges(5,6) * t549 + t553 * t865;
t459 = t549 * t512 - t711;
t460 = Ifges(5,5) * t549 + t512 * t553;
t571 = t601 + t858;
t6 = m(5) * (t396 * t420 + t397 * t421) + t421 * t501 + t397 * t502 + t420 * t503 + t396 * t504 + t472 * t354 + t473 * t353 + t217 * t423 + t218 * t422 + t372 * t194 + t373 * t193 + t99 * t287 + t87 * t288 + t100 * t285 + t88 * t286 + (t232 * t462 + t233 * t463) * mrSges(6,3) + t189 * t877 + t190 * t820 + t188 * t824 + t337 * t804 + t335 * t805 + t336 * t806 + t338 * t807 + (t633 * mrSges(4,1) + Ifges(6,5) * t807 + Ifges(7,5) * t900 + Ifges(6,6) * t806 + Ifges(7,6) * t877 + t458 * t795 + t460 * t790 + t532 * t481 - t868 * t549) * t549 + t191 * t900 + m(6) * (t217 * t232 + t218 * t233 + t472 * t473) + m(7) * (t100 * t88 + t372 * t373 + t87 * t99) + ((-Ifges(5,3) + Ifges(4,1) - Ifges(4,2) - Ifges(6,3) - Ifges(7,3) + (m(5) * t532 + t508) * t532) * t549 + t633 * mrSges(4,2) + t457 * t795 + t459 * t790 + t571 + t740 - t741 + t868 * t553) * t553;
t610 = t6 * qJD(1) + t14 * qJD(2);
t482 = t549 * t509;
t483 = t511 * t549;
t557 = -t929 - t218 * t725 - t472 * t627 - t88 * t739 + t188 * t900 + t881 * t826 + t907 * t822 + t335 * t807 + t462 * t356 / 0.2e1 + t621 * t787;
t7 = t222 * t615 + t221 * t616 - t217 * t724 + t396 * t501 - t397 * t503 - t87 * t735 + t415 * t193 + t190 * t877 + t95 * t287 + t96 * t285 + ((t634 - t459 / 0.2e1 + t482 / 0.2e1 + t396 * mrSges(5,3) - mrSges(5,2) * t513) * t548 + (t710 / 0.2e1 - t483 / 0.2e1 - t457 / 0.2e1 - t397 * mrSges(5,3) + mrSges(5,1) * t513 + (m(6) * t472 + t353) * pkin(4)) * t552) * t549 - t557 + m(6) * (t217 * t221 + t218 * t222) + m(7) * (t372 * t415 + t87 * t95 + t88 * t96) + t890 * t806;
t609 = t7 * qJD(1) - t15 * qJD(2);
t10 = t218 * t616 - t217 * t615 + t193 * t784 - t91 * t287 - t92 * t285 - m(7) * (-t372 * t784 + t87 * t91 + t88 * t92) + t557 + t623 * t625 + (-t337 / 0.2e1 - t355 / 0.2e1 + t217 * mrSges(6,3)) * t463;
t608 = -t10 * qJD(1) - t18 * qJD(2);
t13 = t929 + t87 * t285 - t88 * t287 + t658 * t789 - (-t88 * mrSges(7,3) - t188 / 0.2e1 + t907 / 0.2e1) * t343 + (t881 / 0.2e1 - t623) * t625;
t607 = t13 * qJD(1) - t26 * qJD(2);
t66 = m(7) * (-t343 * t348 + t345 * t625 - t670) + m(6) * (-t462 * t465 + t463 * t464 - t670) + m(5) * (-0.1e1 + t652) * t670;
t606 = t14 * qJD(1) + t66 * qJD(2);
t604 = -t372 * t495 - t436 * t462;
t596 = -t482 / 0.4e1 + t459 / 0.4e1 + t503 * t845;
t594 = t833 * t935 - t930;
t593 = t683 / 0.2e1 - t685 / 0.2e1;
t592 = -t678 / 0.2e1 - t668 / 0.2e1;
t589 = -t904 - t627;
t588 = (t343 * t550 + t546 * t625) * t843;
t586 = t945 * t848;
t403 = -mrSges(6,1) * t494 - mrSges(6,2) * t495;
t461 = -t783 + t786;
t19 = -pkin(3) * t508 + (-t509 / 0.2e1 + t512 / 0.2e1 + pkin(4) * t403) * t548 + (t812 - t407 / 0.2e1) * t495 + m(6) * t536 * t786 + (t865 / 0.2e1 + t511 / 0.2e1) * t552 + t953 * t461 + t956;
t561 = -t461 * t553 * t849 - t529 * t851;
t564 = (t345 * t478 + t348 * t479) * t848 + (t464 * t551 + t465 * t547) * t649 + t622;
t23 = (t878 + t816) * t735 + (t938 + t402 / 0.2e1 + t797 + t597) * t553 + t561 + t564;
t554 = (-t419 * t928 + t472 * t786 + t910 * t866) * t850 + t221 * t635 + t461 * t836 + t241 * t811 + (t372 * t461 + t415 * t436 + (t88 + t95) * t935 + (-t87 + t96) * t216) * t848 + t540 * t788 + (t394 * t838 + t878 * t96) * mrSges(7,3) + t631 * t719 + t923 + t954;
t556 = t532 * t797 - pkin(8) * t875 + (pkin(3) * t841 - t511 / 0.4e1 - t865 / 0.4e1) * t548 + (pkin(3) * t842 + t512 / 0.4e1 - t509 / 0.4e1 + (t536 * t850 + t403 / 0.2e1) * pkin(4)) * t552;
t560 = (t100 * t479 + t478 * t99) * t849 + t420 * t842 + t421 * t841 + t288 * t803 + t286 * t802 - (t232 * t551 + t233 * t547) * t844 / 0.2e1;
t580 = -t457 / 0.4e1 - t483 / 0.4e1 + pkin(4) * t819 + t501 * t845;
t3 = (-Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1 + t921 + t556) * t549 + t554 + t560 + (-t676 / 0.2e1 - t667 / 0.2e1) * pkin(4) + (0.3e1 / 0.4e1 * t710 + t580) * t548 + (-t711 / 0.2e1 + t596) * t552 + t571 - t859 + t864;
t585 = t3 * qJD(1) - t23 * qJD(2) + t19 * qJD(3);
t20 = t407 * t798 + t495 * t812 - t783 * t953 + t956;
t581 = (t882 + t402) * t789 + (t638 + t639) * t625 + (t635 + t636) * t462 + (t637 - t640) * t343;
t562 = t495 * t782 * t848 + t581;
t28 = t562 + t926;
t569 = t546 * t286 / 0.2e1 + t288 * t791 + (t100 * t546 + t550 * t99) * t848;
t577 = t932 * t935 + (-t87 + t92) * t216;
t578 = Ifges(6,3) * t793 + t617 - t858 + t859;
t5 = t594 + t578 + (t681 / 0.2e1 + t689 / 0.2e1 + t604 * t849 + t569) * pkin(5) + (t216 * t900 + t816 * t92 + t87 * t878 + t877 * t935 + t899 * t932) * mrSges(7,3) + t577 * t849 + (-t355 / 0.4e1 - t337 / 0.4e1) * t494 + (t404 / 0.4e1 + t906 / 0.4e1 + t599) * t553 + t902;
t584 = -t5 * qJD(1) + t28 * qJD(2) + t20 * qJD(3);
t34 = t913 + (t908 / 0.2e1 - t245 / 0.2e1) * t394 + (t248 / 0.2e1 + t880 / 0.2e1) * t624;
t44 = t632 - t659;
t555 = (t881 / 0.4e1 + t190 / 0.4e1) * t624 + (mrSges(7,3) * t835 + t248 / 0.4e1 + t880 / 0.4e1) * t625 - (mrSges(7,3) * t834 + t925) * t343 + t906 * t788 + t934 + t954;
t9 = t555 - t617;
t583 = t9 * qJD(1) + t44 * qJD(2) + t34 * qJD(3);
t558 = (t714 / 0.2e1 + t717 / 0.2e1) * pkin(4) * t553 + (t478 * t91 + t479 * t92 + t485 * t87 + t486 * t88) * t848 + t287 * t801 + t285 * t800 + t769;
t12 = t631 * mrSges(6,2) + (-t218 / 0.2e1 - t221 / 0.2e1) * mrSges(6,1) + ((-t679 / 0.2e1 + t669 / 0.2e1) * pkin(5) + t593) * mrSges(7,3) + t558 + t886;
t572 = t721 + (-t714 - t717) * pkin(4) - t720;
t148 = m(7) * (t478 * t485 + t479 * t486) + t572;
t27 = (t809 + t895) * mrSges(6,2) + (t808 + t810) * mrSges(6,1) + pkin(5) * t586 + (t682 / 0.2e1 + t684 / 0.2e1 + t486 * t816 + t394 * t801 + t592 * pkin(5)) * mrSges(7,3) + t950 + t887;
t567 = (-t343 * t486 + t485 * t625 + t660) * t848;
t64 = -t588 / 0.2e1 + t567;
t576 = t12 * qJD(1) + t64 * qJD(2) - t27 * qJD(3) + t148 * qJD(4);
t566 = (-t343 * t802 + t625 * t803) * mrSges(7,3) + t478 * t832 + t287 * t802 + t630;
t17 = (t840 + t96 / 0.2e1) * mrSges(7,2) + (t839 + t838) * mrSges(7,1) + t566 - t630;
t575 = t17 * qJD(1) + t150 * qJD(4) + t949;
t123 = (t800 + t618) * mrSges(7,2) + (-t485 / 0.2e1 + t619) * mrSges(7,1);
t565 = (t285 * t791 + t287 * t796 + (-t343 * t796 + t625 * t792) * mrSges(7,3)) * pkin(5) + t630;
t22 = (t840 + t92 / 0.2e1) * mrSges(7,2) + (t839 - t91 / 0.2e1) * mrSges(7,1) + t565 - t630;
t506 = (mrSges(7,1) * t546 + mrSges(7,2) * t550) * pkin(5);
t570 = -qJD(1) * t22 - qJD(4) * t123 + qJD(5) * t506 - t949;
t493 = t506 * qJD(6);
t124 = -t720 / 0.2e1 + t721 / 0.2e1 + t618 * mrSges(7,2) + t619 * mrSges(7,1);
t50 = t588 / 0.2e1 + t567 + t589;
t43 = t632 + t659;
t29 = t562 - t926;
t25 = t479 * t640 + t478 * t638 + t486 * t639 + t485 * t637 + (mrSges(7,3) * t592 + t586) * pkin(5) + 0.2e1 * t599 + t620 + t950 - t887;
t24 = t508 * t789 + t553 * t597 - t561 + t564 + t581;
t21 = t565 + t769 + t855;
t16 = t566 + t768 + t855;
t11 = t743 / 0.2e1 - t742 / 0.2e1 - t745 / 0.2e1 - t744 / 0.2e1 + t558 + t643 * t778 + t593 * mrSges(7,3) + t550 * pkin(5) * t642 + t621 - t886;
t8 = t555 + t617;
t4 = t578 + t218 * t636 + t92 * t639 + t91 * t637 + t569 * pkin(5) + (pkin(5) * t604 + t577) * t848 - t594 - (t681 + t689) * pkin(5) / 0.2e1 + t923;
t2 = t556 * t549 + t554 - t560 + t578 + Ifges(5,3) * t793 + (t710 / 0.4e1 + t580) * t548 - Ifges(5,6) * t672 / 0.2e1 + (t676 + t667) * pkin(4) / 0.2e1 + (t596 + t634) * t552;
t1 = qJD(3) * t14 - qJD(4) * t15 - qJD(5) * t18 - qJD(6) * t26;
t30 = [qJD(3) * t6 + qJD(4) * t7 - qJD(5) * t10 + qJD(6) * t13, t1, t2 * qJD(4) + t4 * qJD(5) + t8 * qJD(6) + (t509 * t795 + t511 * t790 + t532 * t709 + Ifges(4,5)) * t651 + ((t232 * t866 + t233 * t419 + t473 * t536) * t850 + (t100 * t216 + t373 * t436 + t935 * t99) * t848 + (-pkin(3) * t514 + pkin(8) * t603) * t852) * t853 + t610 + (mrSges(4,2) * t513 - Ifges(4,6) * t549 + t603 * mrSges(5,3) + t536 * t354 + t494 * t336 / 0.2e1 - pkin(3) * t481 + t473 * t403 + t436 * t194 + t419 * t422 + t373 * t241 + t216 * t286 - t99 * t730 + t866 * t423 + t100 * t870 + t189 * t878 + t248 * t820 + t245 * t824 + t408 * t804 + t406 * t805 + t458 * t790 + t460 * t794 + t338 * t798 + (t665 - t674) * pkin(8) + t191 * t899 + (Ifges(5,5) * t548 - Ifges(6,5) * t495 + Ifges(7,5) * t394 + Ifges(5,6) * t552 + Ifges(6,6) * t494 + Ifges(7,6) * t624) * t793 + t232 * t718 + t233 * t719 + t935 * t288) * qJD(3), t2 * qJD(3) + (m(7) * (t478 * t95 + t479 * t96) - t771 + t772 + t743 - t742 - t396 * mrSges(5,2) - t397 * mrSges(5,1) - Ifges(5,5) * t673 - Ifges(5,6) * t671 + t621 + (t683 - t685) * mrSges(7,3)) * qJD(4) + t11 * qJD(5) + t16 * qJD(6) + (m(6) * (t221 * t551 + t222 * t547) + (t462 * t547 - t463 * t551) * mrSges(6,3)) * t754 + t609, t4 * qJD(3) + t11 * qJD(4) + (t621 - t744 - t745 - t773 + t774) * qJD(5) + t21 * qJD(6) + (m(7) * (t546 * t92 + t550 * t91) + (-t669 + t679) * mrSges(7,3)) * t753 + t608, t8 * qJD(3) + t16 * qJD(4) + t21 * qJD(5) + (t658 - t775 - t776) * qJD(6) + t607; t1, t66 * qJD(3), t24 * qJD(4) + t29 * qJD(5) + t43 * qJD(6) + (mrSges(5,3) * t652 - mrSges(4,2)) * t651 + ((t216 * t348 + t345 * t935 + t436 * t549) * t848 + (t419 * t465 + t464 * t866 + t536 * t549) * t850 + (t652 * t781 - t777) * t852) * t853 + t606 + (-t345 * t730 + t348 * t870 + t464 * t718 + t465 * t719 + (t241 + t403 + t709) * t549) * qJD(3), -t708 + t24 * qJD(3) + (-mrSges(5,1) * t671 + mrSges(5,2) * t673 + t589) * qJD(4) + t50 * qJD(5) + 0.2e1 * (t660 * t848 + (pkin(4) * t463 * t547 + t462 * t785) * t850) * qJD(4) + t856, -t705 + t29 * qJD(3) + t50 * qJD(4) + (t588 + t589) * qJD(5) + t856, t43 * qJD(3) - qJD(6) * t904 + t582 * t888 - t703; qJD(4) * t3 - qJD(5) * t5 + qJD(6) * t9 - t610, -qJD(4) * t23 + qJD(5) * t28 + qJD(6) * t44 - t606, qJD(4) * t19 + qJD(5) * t20 + qJD(6) * t34 (m(7) * (-t216 * t478 + t479 * t935) + t540 - t758 + t863 * pkin(8) + (-t682 - t684) * mrSges(7,3) + t952) * qJD(4) + t25 * qJD(5) + t957 + (m(6) * (-t419 * t551 + t547 * t866) + (-t494 * t551 + t495 * t547) * mrSges(6,3)) * t754 + t585, t25 * qJD(4) + t952 * qJD(5) + t957 + (m(7) * t945 + (-t668 - t678) * mrSges(7,3)) * t753 + t584 (t906 + t936) * qJD(6) + t583 + t888 * t955; -qJD(3) * t3 + qJD(5) * t12 + qJD(6) * t17 - t609, qJD(3) * t23 + qJD(5) * t64 + t708 + t857, -qJD(5) * t27 - t585 + t951, qJD(5) * t148 + t650 ((t485 * t550 + t486 * t546) * t843 + t572) * qJD(5) + t124 * qJD(6) + t576, t124 * qJD(5) + t575 + t650; qJD(3) * t5 - qJD(4) * t12 + qJD(6) * t22 - t608, -qJD(3) * t28 - qJD(4) * t64 + t705 + t857, qJD(4) * t27 - t584 + t951, qJD(6) * t123 - t576, -t493, -t493 - t570; -qJD(3) * t9 - qJD(4) * t17 - qJD(5) * t22 - t607, -t44 * qJD(3) - t579 * t888 + t703, -t888 * t944 - t583, -qJD(5) * t123 - t575, t570, 0;];
Cq  = t30;
