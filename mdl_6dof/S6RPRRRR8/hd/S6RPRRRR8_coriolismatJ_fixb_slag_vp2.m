% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:19:50
% EndTime: 2019-03-09 07:20:16
% DurationCPUTime: 15.00s
% Computational Cost: add. (37238->752), mult. (69284->953), div. (0->0), fcn. (76065->8), ass. (0->440)
t867 = qJD(3) + qJD(4);
t488 = sin(qJ(4));
t489 = sin(qJ(3));
t770 = cos(qJ(4));
t771 = cos(qJ(3));
t453 = t488 * t489 - t770 * t771;
t487 = sin(qJ(5));
t484 = t487 ^ 2;
t490 = cos(qJ(5));
t485 = t490 ^ 2;
t852 = t484 + t485;
t572 = t852 * t453;
t486 = sin(qJ(6));
t769 = cos(qJ(6));
t560 = t770 * t769;
t611 = t490 * t770;
t405 = (-t486 * t611 - t487 * t560) * pkin(3);
t612 = t487 * t770;
t406 = (-t486 * t612 + t490 * t560) * pkin(3);
t642 = t770 * pkin(3);
t566 = -t642 / 0.2e1;
t556 = t487 * t566;
t818 = m(7) * pkin(5);
t652 = t818 / 0.2e1;
t754 = mrSges(6,2) * t490;
t905 = -mrSges(7,2) / 0.2e1;
t906 = mrSges(7,1) / 0.2e1;
t856 = t405 * t906 + t406 * t905;
t920 = (t405 * t769 + t406 * t486) * t652 + mrSges(6,1) * t556 + t566 * t754 + t856;
t546 = t486 * t490 + t487 * t769;
t715 = t490 * mrSges(6,1);
t717 = t487 * mrSges(6,2);
t844 = -t486 * t487 + t769 * t490;
t719 = t546 * mrSges(7,2);
t721 = t844 * mrSges(7,1);
t854 = t721 / 0.2e1 - t719 / 0.2e1;
t919 = -t717 / 0.2e1 + t715 / 0.2e1 + (t486 * t546 + t769 * t844) * t652 + t854;
t454 = t488 * t771 + t489 * t770;
t372 = -pkin(4) * t453 + t454 * pkin(9);
t814 = pkin(1) + pkin(7);
t629 = t489 * t814;
t458 = -pkin(8) * t489 - t629;
t561 = t771 * t814;
t536 = -pkin(8) * t771 - t561;
t373 = t458 * t488 - t770 * t536;
t875 = t487 * t373;
t221 = t490 * t372 + t875;
t680 = t454 * t490;
t559 = -t453 * pkin(5) + pkin(10) * t680;
t134 = t221 + t559;
t874 = t490 * t373;
t222 = t487 * t372 - t874;
t681 = t454 * t487;
t649 = pkin(10) * t681;
t177 = t649 + t222;
t100 = t486 * t134 + t177 * t769;
t816 = -mrSges(6,2) / 0.2e1;
t817 = mrSges(6,1) / 0.2e1;
t99 = t134 * t769 - t486 * t177;
t862 = t100 * t905 + t906 * t99;
t918 = -t221 * t817 - t222 * t816 - (t100 * t486 + t769 * t99) * t652 - t862;
t643 = t771 * pkin(3);
t350 = t643 + t372;
t205 = t490 * t350 + t875;
t206 = t487 * t350 - t874;
t133 = t205 + t559;
t165 = t649 + t206;
t95 = t133 * t769 - t486 * t165;
t96 = t486 * t133 + t165 * t769;
t863 = t96 * t905 + t906 * t95;
t917 = -t205 * t817 - t206 * t816 - (t486 * t96 + t769 * t95) * t652 - t863;
t898 = t844 * t454;
t233 = -mrSges(7,1) * t453 + mrSges(7,3) * t898;
t644 = pkin(5) * t769;
t763 = pkin(5) * t486;
t784 = -t453 / 0.2e1;
t803 = -t898 / 0.2e1;
t859 = t546 * t454;
t805 = t859 / 0.2e1;
t549 = Ifges(7,5) * t803 + Ifges(7,6) * t805;
t783 = t453 / 0.2e1;
t841 = Ifges(7,3) * t783 - t549;
t231 = mrSges(7,2) * t453 + mrSges(7,3) * t859;
t888 = t231 / 0.2e1;
t916 = t841 - t233 * t644 / 0.2e1 - t763 * t888 - Ifges(6,3) * t784;
t767 = pkin(3) * t488;
t475 = pkin(9) + t767;
t755 = pkin(10) + t475;
t446 = t755 * t487;
t447 = t755 * t490;
t343 = -t486 * t446 + t447 * t769;
t571 = -t769 * t446 - t447 * t486;
t915 = -t343 * mrSges(7,1) - t571 * mrSges(7,2);
t813 = -pkin(10) - pkin(9);
t469 = t813 * t487;
t470 = t813 * t490;
t393 = t486 * t469 - t470 * t769;
t570 = t769 * t469 + t470 * t486;
t914 = -t393 * mrSges(7,1) - t570 * mrSges(7,2);
t913 = t862 - t841;
t912 = t863 - t841;
t147 = -Ifges(7,4) * t898 + Ifges(7,2) * t859 - t453 * Ifges(7,6);
t149 = -Ifges(7,1) * t898 + Ifges(7,4) * t859 - t453 * Ifges(7,5);
t482 = Ifges(6,4) * t490;
t557 = Ifges(6,2) * t487 - t482;
t261 = -Ifges(6,6) * t453 + t454 * t557;
t751 = Ifges(6,4) * t487;
t466 = Ifges(6,1) * t490 - t751;
t263 = -Ifges(6,5) * t453 - t454 * t466;
t463 = Ifges(6,2) * t490 + t751;
t465 = Ifges(6,1) * t487 + t482;
t594 = -t680 / 0.2e1;
t595 = t681 / 0.2e1;
t772 = t490 / 0.2e1;
t774 = t487 / 0.2e1;
t785 = t546 / 0.2e1;
t787 = t844 / 0.2e1;
t438 = Ifges(7,4) * t844;
t371 = Ifges(7,1) * t546 + t438;
t792 = t371 / 0.2e1;
t749 = Ifges(7,4) * t546;
t369 = Ifges(7,2) * t844 + t749;
t795 = t369 / 0.2e1;
t524 = -Ifges(5,5) * t454 + Ifges(5,6) * t453 + t147 * t787 + t149 * t785 + t261 * t772 + t263 * t774 + t859 * t795 - t898 * t792 + t463 * t595 + t465 * t594 + (Ifges(6,5) * t487 + Ifges(7,5) * t546 + Ifges(6,6) * t490 + Ifges(7,6) * t844) * t784;
t461 = -t715 + t717;
t851 = t770 * t458 + t488 * t536;
t873 = t851 * t461;
t880 = t851 * mrSges(5,1);
t882 = t373 * mrSges(5,2);
t367 = t719 - t721;
t872 = -pkin(5) * t681 + t851;
t897 = t872 * t367;
t911 = t524 + t873 + t882 - t880 + t897;
t328 = t546 * t453;
t609 = t769 * t328;
t736 = t328 * mrSges(7,1);
t325 = t844 * t453;
t737 = t325 * mrSges(7,2);
t660 = t736 / 0.2e1 + t737 / 0.2e1;
t682 = t453 * t490;
t683 = t453 * t487;
t698 = t325 * t486;
t512 = (t609 - t698) * t652 + t683 * t817 + mrSges(6,2) * t682 / 0.2e1 + t660;
t365 = mrSges(7,1) * t546 + mrSges(7,2) * t844;
t685 = t453 * t365;
t718 = t546 * mrSges(7,3);
t618 = -t718 / 0.2e1;
t619 = t718 / 0.2e1;
t850 = (-t618 - t619) * t898;
t651 = pkin(5) * t683;
t885 = m(7) * t651;
t462 = t487 * mrSges(6,1) + t754;
t684 = t453 * t462;
t887 = -t684 / 0.2e1;
t894 = -t885 / 0.2e1 - t685 / 0.2e1 + t887 + t512 + t850;
t910 = qJD(2) * t894;
t909 = -t873 / 0.2e1 + t880 / 0.2e1 - t882 / 0.2e1 - t897 / 0.2e1;
t614 = t685 / 0.2e1 + t850;
t555 = t614 + t660;
t389 = t453 * t454;
t697 = t859 * t328;
t699 = t898 * t325;
t62 = m(7) * (t389 - t697 - t699) + m(6) * (-t454 * t572 + t389);
t663 = t62 * qJD(2);
t893 = t885 / 0.2e1 + t684 / 0.2e1 + t512 + t614;
t908 = t893 * qJD(5) + t555 * qJD(6) + t663;
t554 = t614 - t660;
t907 = -qJD(5) * t894 + qJD(6) * t554 - t663;
t471 = t489 * pkin(3) + qJ(2);
t765 = pkin(4) * t454;
t349 = pkin(9) * t453 + t471 + t765;
t203 = t490 * t349 - t487 * t851;
t162 = pkin(10) * t682 + t203;
t204 = t487 * t349 + t490 * t851;
t163 = pkin(10) * t683 + t204;
t671 = t486 * t163;
t105 = t162 * t769 - t671;
t132 = t454 * pkin(5) + t162;
t93 = t132 * t769 - t671;
t904 = t105 - t93;
t266 = t373 - t651;
t903 = t266 * t872;
t476 = -t642 - pkin(4);
t761 = t490 * pkin(5);
t459 = t476 - t761;
t902 = t459 * t872;
t477 = -pkin(4) - t761;
t901 = t477 * t872;
t654 = Ifges(7,5) * t844 - Ifges(7,6) * t546;
t66 = t654 + t915;
t900 = t66 * qJD(6);
t73 = t654 + t914;
t899 = t73 * qJD(6);
t178 = -mrSges(7,1) * t325 + mrSges(7,2) * t328;
t232 = -mrSges(7,2) * t454 + t328 * mrSges(7,3);
t806 = -t859 / 0.2e1;
t234 = mrSges(7,1) * t454 + t325 * mrSges(7,3);
t808 = t234 / 0.2e1;
t503 = t178 * t783 + t232 * t806 - t898 * t808 + (t699 / 0.2e1 + t697 / 0.2e1) * mrSges(7,3);
t575 = -t454 * mrSges(5,1) + t453 * mrSges(5,2);
t896 = m(5) * t471 - t575;
t574 = -mrSges(7,1) * t898 + mrSges(7,2) * t859;
t895 = qJD(6) * t574;
t750 = Ifges(7,4) * t325;
t148 = Ifges(7,2) * t328 + t454 * Ifges(7,6) - t750;
t315 = Ifges(7,4) * t328;
t150 = -Ifges(7,1) * t325 + Ifges(7,5) * t454 + t315;
t179 = -mrSges(7,1) * t859 - mrSges(7,2) * t898;
t180 = -t736 - t737;
t345 = t462 * t454;
t353 = mrSges(6,2) * t453 + mrSges(6,3) * t681;
t355 = -t453 * mrSges(6,1) + mrSges(6,3) * t680;
t748 = Ifges(6,5) * t454;
t264 = -t453 * t466 + t748;
t666 = t490 * t264;
t745 = Ifges(6,6) * t454;
t262 = t453 * t557 + t745;
t669 = t487 * t262;
t804 = t328 / 0.2e1;
t807 = -t325 / 0.2e1;
t481 = Ifges(6,5) * t490;
t744 = Ifges(6,6) * t487;
t879 = Ifges(5,4) - t481 / 0.2e1 + t744 / 0.2e1;
t610 = t769 * t163;
t94 = t486 * t132 + t610;
t892 = t872 * t180 + (t669 / 0.2e1 - t666 / 0.2e1 + t549 + t879 * t454) * t454 - t851 * t684 + t147 * t804 + t148 * t805 + t149 * t807 + t150 * t803 + t266 * t179 + t203 * t355 + t204 * t353 + t94 * t231 + t93 * t233 - t373 * t345 + t471 * (-mrSges(5,1) * t453 - mrSges(5,2) * t454);
t77 = t93 * t328;
t891 = t266 * t454 - t94 * t325 + t453 * t872 + t77;
t884 = pkin(4) * t851;
t104 = -t486 * t162 - t610;
t878 = t104 + t94;
t877 = t373 * t488;
t693 = t373 * t851;
t876 = t476 * t851;
t344 = t461 * t453;
t822 = m(7) / 0.2e1;
t827 = t453 ^ 2;
t871 = (-t761 * t827 - t878 * t859 + t898 * t904) * t822 + t344 * t783 + t503;
t370 = Ifges(7,1) * t844 - t749;
t580 = t369 / 0.4e1 - t370 / 0.4e1;
t368 = -Ifges(7,2) * t546 + t438;
t581 = t368 / 0.4e1 + t371 / 0.4e1;
t866 = t373 * t462 / 0.2e1 + t581 * t328 + t580 * t325;
t865 = -t204 * t682 + t373 * t454 + (t203 * t487 + t851) * t453;
t427 = t453 * t767;
t855 = -t454 * t642 - t427;
t853 = t465 - t557;
t356 = t454 * mrSges(6,1) + mrSges(6,3) * t682;
t664 = t490 * t356;
t354 = -mrSges(6,2) * t454 + mrSges(6,3) * t683;
t668 = t487 * t354;
t848 = -t664 / 0.2e1 - t668 / 0.2e1;
t665 = t490 * t353;
t667 = t487 * t355;
t847 = t665 / 0.2e1 - t667 / 0.2e1;
t686 = t393 * t325;
t688 = t570 * t328;
t846 = -t688 / 0.2e1 + t686 / 0.2e1;
t196 = t571 * t328;
t197 = t343 * t325;
t845 = -t196 / 0.2e1 + t197 / 0.2e1;
t840 = -t489 * mrSges(4,1) - t771 * mrSges(4,2);
t348 = t453 * t465;
t811 = t180 / 0.2e1;
t838 = pkin(5) * t811 + t348 / 0.4e1;
t577 = t463 / 0.4e1 - t466 / 0.4e1;
t837 = qJD(2) * t554;
t753 = mrSges(6,3) * t453;
t833 = t848 + t852 * t753 / 0.2e1;
t775 = t477 / 0.2e1;
t591 = t178 * t775;
t789 = t570 / 0.2e1;
t821 = -pkin(4) / 0.2e1;
t832 = t232 * t789 + t344 * t821 - t393 * t808 + t591;
t781 = t459 / 0.2e1;
t593 = t178 * t781;
t776 = t476 / 0.2e1;
t801 = t571 / 0.2e1;
t831 = t232 * t801 - t343 * t808 + t344 * t776 + t593;
t773 = -t490 / 0.2e1;
t830 = (-Ifges(6,3) - Ifges(5,2) - Ifges(7,3) + Ifges(5,1)) * t454 + t263 * t773 + t261 * t774;
t182 = Ifges(7,1) * t328 + t750;
t782 = t454 / 0.4e1;
t828 = t266 * t365 / 0.2e1 + t654 * t782 - (-t182 / 0.4e1 + t148 / 0.4e1) * t546;
t826 = 2 * qJD(4);
t824 = m(6) / 0.2e1;
t823 = -m(7) / 0.2e1;
t820 = -pkin(9) / 0.2e1;
t819 = m(6) * pkin(3);
t815 = mrSges(7,3) / 0.2e1;
t809 = t232 / 0.2e1;
t802 = -t571 / 0.2e1;
t800 = -t343 / 0.2e1;
t347 = t453 * t463;
t799 = t347 / 0.4e1;
t797 = -t367 / 0.2e1;
t790 = -t570 / 0.2e1;
t788 = -t393 / 0.2e1;
t780 = -t463 / 0.2e1;
t777 = -t475 / 0.2e1;
t766 = pkin(4) * t345;
t764 = pkin(4) * t462;
t762 = pkin(5) * t487;
t760 = t93 * mrSges(7,2);
t759 = t94 * mrSges(7,1);
t747 = Ifges(7,5) * t325;
t743 = Ifges(7,6) * t328;
t739 = t104 * mrSges(7,1);
t738 = t105 * mrSges(7,2);
t720 = t844 * mrSges(7,3);
t6 = m(7) * (t93 * t95 + t94 * t96 + t903) + m(6) * (t203 * t205 + t204 * t206 + t693) + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t489) * t489 - Ifges(5,4) * t827 + t206 * t354 + t205 * t356 + t95 * t234 + t96 * t232 + (-Ifges(6,5) * t682 + Ifges(6,6) * t683 + t743 - t747) * t784 + ((-Ifges(4,1) + Ifges(4,2)) * t489 + qJ(2) * mrSges(4,1) - Ifges(4,4) * t771) * t771 + t830 * t453 + t896 * t643 + t892;
t714 = t6 * qJD(1);
t8 = m(7) * (t100 * t94 + t93 * t99 + t903) + m(6) * (t203 * t221 + t204 * t222 + t693) + (t747 / 0.2e1 - t743 / 0.2e1 - t879 * t453 + t830) * t453 + t222 * t354 + t221 * t356 + t99 * t234 + t100 * t232 + t892;
t713 = t8 * qJD(1);
t181 = Ifges(7,2) * t325 + t315;
t659 = Ifges(7,5) * t328 + Ifges(7,6) * t325;
t510 = t266 * t178 + (t150 / 0.2e1 + t181 / 0.2e1) * t328 + (-t182 / 0.2e1 + t94 * mrSges(7,3) + t148 / 0.2e1) * t325 - mrSges(7,3) * t77 + t454 * t659 / 0.2e1;
t626 = t748 / 0.2e1;
t9 = m(7) * (t93 * t104 + t94 * t105) + t105 * t232 + t104 * t234 + t373 * t344 - t204 * t356 + t203 * t354 + ((-t203 * mrSges(6,3) + t626 + t264 / 0.2e1 + t347 / 0.2e1) * t487 + (t204 * mrSges(6,3) + t745 / 0.2e1 - t348 / 0.2e1 + t262 / 0.2e1 + (-m(7) * t266 - t180) * pkin(5)) * t490) * t453 + t510;
t712 = t9 * qJD(1);
t500 = (t179 / 0.2e1 + t356 * t774 + t354 * t773 - t345 / 0.2e1) * t453 + (t811 + t887 + t847) * t454 + t898 * t888 + t232 * t807 + t233 * t806 + t234 * t804;
t705 = t206 * t490;
t706 = t205 * t487;
t553 = t705 - t706;
t494 = (t553 * t454 + t865) * t824 + (-t859 * t95 + t898 * t96 + t891) * t822 + t500;
t544 = m(7) * (t343 * t546 + t571 * t844);
t15 = -t544 / 0.2e1 + t494;
t711 = t15 * qJD(1);
t702 = t222 * t490;
t703 = t221 * t487;
t552 = t702 - t703;
t493 = (t552 * t454 + t865) * t824 + (t100 * t898 - t859 * t99 + t891) * t822 + t500;
t543 = m(7) * (t393 * t546 + t570 * t844);
t17 = -t543 / 0.2e1 + t493;
t710 = t17 * qJD(1);
t18 = t93 * t232 - t94 * t234 + t510;
t709 = t18 * qJD(1);
t20 = t454 * t833 + t871 - t919;
t708 = t20 * qJD(1);
t22 = t503 - t854;
t704 = t22 * qJD(1);
t696 = t571 * t233;
t695 = t343 * t231;
t37 = t546 * t232 + t844 * t234 + t668 + t664 + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(7) * (t546 * t94 + t844 * t93) + m(6) * (t203 * t490 + t204 * t487) - t840 + t896;
t694 = t37 * qJD(1);
t689 = t570 * t233;
t687 = t393 * t231;
t402 = t454 * t461;
t679 = t459 * t179;
t678 = t459 * t365;
t677 = t476 * t345;
t676 = t476 * t462;
t675 = t477 * t179;
t674 = t477 * t365;
t653 = mrSges(7,3) * t763;
t650 = pkin(5) * t682;
t647 = pkin(9) * t667;
t646 = pkin(9) * t665;
t641 = mrSges(6,3) * t706;
t640 = mrSges(6,3) * t705;
t639 = mrSges(6,3) * t703;
t638 = mrSges(6,3) * t702;
t622 = -t720 / 0.2e1;
t621 = t720 / 0.2e1;
t617 = -t93 / 0.2e1 + t105 / 0.2e1;
t616 = t94 / 0.2e1 + t104 / 0.2e1;
t615 = t459 * t454 + t196 - t197;
t613 = t476 * t454 - t475 * t572;
t584 = t264 / 0.4e1 + t799;
t583 = t802 + t801;
t582 = t800 + t343 / 0.2e1;
t579 = t790 + t789;
t578 = t788 + t393 / 0.2e1;
t573 = t481 - t744;
t568 = mrSges(7,3) * t644;
t558 = (t781 + t775) * t365;
t551 = (t484 / 0.2e1 + t485 / 0.2e1) * t753;
t542 = t852 * t770;
t491 = t828 - t669 / 0.4e1 + t105 * t621 + t93 * t622 + t853 * t683 / 0.4e1 + (t181 + t150) * t844 / 0.4e1 + t666 / 0.4e1 - Ifges(6,6) * t681 / 0.2e1 + t577 * t682 + t878 * t618 + (t626 + t799) * t490 + t838 * t487 + t573 * t782 + t650 * t797 + t866 + t916;
t251 = t266 * t762;
t528 = t343 * t904 + t878 * t571 + t251;
t2 = t491 + (-t459 * t650 + t528) * t822 + t845 * mrSges(7,3) + t833 * t475 + t831 + t917;
t444 = t459 * t762;
t529 = -(-t370 / 0.2e1 + t795) * t546 + (t792 + t368 / 0.2e1) * t844;
t507 = (t465 / 0.2e1 - t557 / 0.2e1) * t490 + (pkin(5) * t367 + t466 / 0.2e1 + t780) * t487 + t529;
t42 = m(7) * t444 + t507 + t676 + t678;
t541 = t2 * qJD(1) + t42 * qJD(3) - t910;
t499 = (pkin(3) * t454 * t542 + t427 + t613) * t824 + (-t405 * t859 + t406 * t898 + t427 + t615) * t822;
t511 = (t454 * t477 - t686 + t688) * t822 + (-pkin(9) * t572 - t765) * t824;
t43 = t499 - t511;
t501 = -t405 * t718 + t406 * t720 + (-mrSges(5,1) + t367 + t461) * t767 + (mrSges(6,3) * t852 - mrSges(5,2)) * t642;
t64 = m(7) * (t343 * t406 + t405 * t571 + t459 * t767) + (t475 * t542 + t476 * t488) * t819 + t501;
t492 = (t876 + (-t203 * t612 + t204 * t611 + t877) * pkin(3)) * t824 + (t100 * t343 + t266 * t767 + t405 * t93 + t406 * t94 + t571 * t99 + t902) * t822 + t696 / 0.2e1 + t695 / 0.2e1 + t405 * t808 + t406 * t809 + t679 / 0.2e1 - t677 / 0.2e1 + t100 * t621 - t639 / 0.2e1 + t638 / 0.2e1 + t99 * t618 + t356 * t556 + pkin(3) * t354 * t611 / 0.2e1 + (t180 - t684) * t767 / 0.2e1 + (t552 * t824 + t847) * t475 - t909;
t495 = -m(6) * (pkin(9) * t553 - t884) / 0.2e1 + (t393 * t96 + t570 * t95 + t901) * t823 - t766 / 0.2e1 - t689 / 0.2e1 - t687 / 0.2e1 - t675 / 0.2e1 + t641 / 0.2e1 - t640 / 0.2e1 + t647 / 0.2e1 - t646 / 0.2e1 + t95 * t619 + t96 * t622 + t909;
t7 = t495 + t492;
t540 = t7 * qJD(1) + t43 * qJD(2) + t64 * qJD(3);
t513 = (t150 / 0.4e1 + t181 / 0.4e1) * t844 + t828;
t496 = (mrSges(7,3) * t802 + t581) * t328 + (t343 * t815 + t580) * t325 + t571 * t809 + t234 * t800 + t593 + t513;
t11 = t496 - t912;
t52 = t529 + t678;
t538 = t11 * qJD(1) + t52 * qJD(3) + t837;
t460 = t477 * t762;
t519 = (t444 + t460) * t822;
t25 = t519 + t466 * t774 + t487 * t780 + t676 / 0.2e1 + t674 / 0.2e1 + t678 / 0.2e1 - t764 / 0.2e1 - t546 * t795 + t370 * t785 + t367 * t762 + (t371 + t368) * t787 + t853 * t772 + (t622 + t621) * (t570 + t571) - t920;
t527 = t393 * t904 + t878 * t570 + t251;
t4 = t491 + (-t477 * t650 + t527) * t822 + t846 * mrSges(7,3) + t833 * pkin(9) + t832 + t918;
t46 = m(7) * t460 + t507 + t674 - t764;
t532 = t4 * qJD(1) + t25 * qJD(3) + t46 * qJD(4) - t910;
t531 = -t546 * t616 + t617 * t844;
t497 = (t393 * t815 + t580) * t325 + (mrSges(7,3) * t790 + t581) * t328 + t570 * t809 + t234 * t788 + t591 + t513;
t13 = t497 - t913;
t515 = t558 + t529;
t38 = t515 - t856;
t57 = t529 + t674;
t530 = t13 * qJD(1) + t38 * qJD(3) + t57 * qJD(4) + t837;
t526 = -t546 * t653 - t568 * t844 + t573 + t654;
t525 = -mrSges(6,3) * t572 - t325 * t720 - t328 * t718 + t454 * t367 + t402 + t575;
t521 = (t465 / 0.4e1 - t557 / 0.4e1) * t453 - t262 / 0.4e1 - t745 / 0.4e1 + t838;
t514 = (t769 * t809 - t486 * t234 / 0.2e1 + (t698 / 0.2e1 - t609 / 0.2e1) * mrSges(7,3)) * pkin(5);
t24 = -mrSges(7,1) * t616 + mrSges(7,2) * t617 + t514;
t457 = (mrSges(7,1) * t486 + mrSges(7,2) * t769) * pkin(5);
t70 = mrSges(7,1) * t582 + mrSges(7,2) * t583;
t79 = mrSges(7,1) * t578 + mrSges(7,2) * t579;
t520 = -t24 * qJD(1) - t70 * qJD(3) - t79 * qJD(4) + t457 * qJD(5);
t498 = Ifges(6,5) * t594 + Ifges(6,6) * t595 + t481 * t782 + t513 + t866 - t916;
t450 = t457 * qJD(6);
t39 = t515 + t856;
t35 = t499 + t511 + t525;
t26 = t519 + (t776 + t821) * t462 + t558 + (-(-t578 - t582) * t546 + (t579 + t583) * t844) * mrSges(7,3) + t507 + t920;
t23 = t503 + t854;
t21 = -t760 / 0.2e1 - t759 / 0.2e1 - t738 / 0.2e1 + t739 / 0.2e1 + t514 + t659;
t19 = (t551 + t848) * t454 + t871 + t919;
t16 = t543 / 0.2e1 + t493;
t14 = t544 / 0.2e1 + t494;
t12 = t497 + t913;
t10 = t496 + t912;
t5 = t498 + (t356 * t820 + ((t477 * t823 + t797) * pkin(5) + t577) * t453 + t584) * t490 + t527 * t822 + (t354 * t820 + t521) * t487 + (t531 + t846) * mrSges(7,3) + pkin(9) * t551 + t832 - t918;
t3 = t475 * t551 + t498 + (t356 * t777 + ((t459 * t823 + t797) * pkin(5) + t577) * t453 + t584) * t490 + t528 * t822 + (t354 * t777 + t521) * t487 + (t531 + t845) * mrSges(7,3) + t831 - t917;
t1 = -t495 + t492 + t524;
t27 = [qJD(2) * t37 + qJD(3) * t6 + qJD(4) * t8 + qJD(5) * t9 + qJD(6) * t18, t694 + m(7) * (t546 * t898 - t844 * t859) * qJD(2) + t14 * qJD(3) + t16 * qJD(4) + t19 * qJD(5) + t23 * qJD(6), t714 + t14 * qJD(2) + (m(7) * (t343 * t96 + t571 * t95 + t902) - t677 - t95 * t718 + t96 * t720 + t695 + t696 + mrSges(4,1) * t629 - t641 + t640 + m(5) * (-t770 * t851 - t877) * pkin(3) + t679 - Ifges(4,6) * t771 + mrSges(4,2) * t561 + (m(6) * t553 + t665 - t667) * t475 - Ifges(4,5) * t489 - t855 * mrSges(5,3) + m(6) * t876 + t911) * qJD(3) + t1 * qJD(4) + t3 * qJD(5) + t10 * qJD(6), t713 + t16 * qJD(2) + t1 * qJD(3) + (t100 * t720 - t718 * t99 + t638 - t639 + t646 - t647 + t675 + t687 + t689 + t766 + t911) * qJD(4) + t5 * qJD(5) + t12 * qJD(6) + ((t100 * t393 + t570 * t99 + t901) * t822 + (pkin(9) * t552 - t884) * t824) * t826, t712 + t19 * qJD(2) + t3 * qJD(3) + t5 * qJD(4) + (Ifges(6,5) * t683 + Ifges(6,6) * t682 + (t104 * t769 + t105 * t486) * t818 + t325 * t653 - t328 * t568 - t738 + t739 - t203 * mrSges(6,2) - t204 * mrSges(6,1) + t659) * qJD(5) + t21 * qJD(6), t709 + t23 * qJD(2) + t10 * qJD(3) + t12 * qJD(4) + t21 * qJD(5) + (t659 - t759 - t760) * qJD(6); qJD(3) * t15 + qJD(4) * t17 + qJD(5) * t20 + qJD(6) * t22 - t694, t867 * t62, t35 * qJD(4) + t711 + (m(5) * t855 + 0.2e1 * t613 * t824 + 0.2e1 * t615 * t822 + t525 + t840) * qJD(3) + t908, t35 * qJD(3) + qJD(4) * t525 + t511 * t826 + t710 + t908, t708 + (m(7) * (-t644 * t898 - t763 * t859) + t402 + t574) * qJD(5) + t895 + t867 * t893, qJD(5) * t574 + t555 * t867 + t704 + t895; -qJD(2) * t15 + qJD(4) * t7 + qJD(5) * t2 + qJD(6) * t11 - t714, qJD(4) * t43 - t711 + t907, qJD(4) * t64 + qJD(5) * t42 + qJD(6) * t52 (m(7) * (t393 * t406 + t405 * t570 + t477 * t767) + (-pkin(4) * t488 + pkin(9) * t542) * t819 + t501) * qJD(4) + t26 * qJD(5) + t39 * qJD(6) + t540, t26 * qJD(4) + ((-t343 * t769 + t486 * t571) * t818 + t526 + t461 * t475 + t915) * qJD(5) + t900 + t541, t39 * qJD(4) + t66 * qJD(5) + t538 + t900; -qJD(2) * t17 - qJD(3) * t7 + qJD(5) * t4 + qJD(6) * t13 - t713, -qJD(3) * t43 - t710 + t907, qJD(5) * t25 + qJD(6) * t38 - t540, qJD(5) * t46 + qJD(6) * t57 ((-t393 * t769 + t486 * t570) * t818 + t526 + t461 * pkin(9) + t914) * qJD(5) + t899 + t532, t73 * qJD(5) + t530 + t899; -qJD(2) * t20 - qJD(3) * t2 - qJD(4) * t4 + qJD(6) * t24 - t712, t867 * t894 - t708, -qJD(4) * t25 + qJD(6) * t70 - t541, qJD(6) * t79 - t532, -t450, -t450 - t520; -qJD(2) * t22 - qJD(3) * t11 - qJD(4) * t13 - qJD(5) * t24 - t709, -t554 * t867 - t704, -qJD(4) * t38 - qJD(5) * t70 - t538, -qJD(5) * t79 - t530, t520, 0;];
Cq  = t27;
