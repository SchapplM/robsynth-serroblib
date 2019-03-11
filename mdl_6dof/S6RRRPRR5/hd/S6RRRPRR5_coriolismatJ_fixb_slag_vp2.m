% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:20:43
% EndTime: 2019-03-09 18:21:14
% DurationCPUTime: 15.51s
% Computational Cost: add. (38182->733), mult. (71059->935), div. (0->0), fcn. (78680->8), ass. (0->439)
t923 = qJ(4) / 0.2e1;
t496 = sin(qJ(5));
t493 = t496 ^ 2;
t499 = cos(qJ(5));
t494 = t499 ^ 2;
t659 = t493 + t494;
t497 = sin(qJ(3));
t498 = sin(qJ(2));
t796 = cos(qJ(3));
t797 = cos(qJ(2));
t454 = t497 * t498 - t796 * t797;
t771 = mrSges(6,3) * t454;
t922 = -t771 * t659 / 0.2e1;
t495 = sin(qJ(6));
t795 = cos(qJ(6));
t552 = -t495 * t496 + t499 * t795;
t789 = pkin(2) * t497;
t398 = t552 * t789;
t553 = t495 * t499 + t795 * t496;
t399 = t553 * t789;
t643 = -t789 / 0.2e1;
t581 = t496 * t643;
t772 = mrSges(6,1) * t499;
t634 = t772 / 0.2e1;
t827 = m(7) * pkin(5);
t657 = t827 / 0.2e1;
t916 = -mrSges(7,2) / 0.2e1;
t917 = mrSges(7,1) / 0.2e1;
t856 = t398 * t917 + t399 * t916;
t921 = (t398 * t795 + t399 * t495) * t657 + mrSges(6,2) * t581 + t634 * t789 + t856;
t322 = t552 * t454;
t324 = t553 * t454;
t179 = mrSges(7,1) * t324 + mrSges(7,2) * t322;
t490 = t496 * pkin(10);
t823 = pkin(3) + pkin(9);
t674 = t496 * t823;
t459 = -t490 - t674;
t460 = (-pkin(10) - t823) * t499;
t366 = t459 * t795 + t495 * t460;
t758 = t324 * mrSges(7,3);
t760 = t322 * mrSges(7,3);
t585 = -t459 * t495 + t795 * t460;
t809 = -t585 / 0.2e1;
t455 = t497 * t797 + t498 * t796;
t228 = mrSges(7,1) * t455 - t758;
t820 = -t228 / 0.2e1;
t491 = t496 * pkin(5);
t479 = qJ(4) + t491;
t915 = t479 / 0.2e1;
t920 = t179 * t915 + t760 * t809 + (-t758 / 0.2e1 + t820) * t366;
t325 = t553 * t455;
t816 = t325 / 0.2e1;
t858 = t552 * t455;
t818 = t858 / 0.2e1;
t559 = Ifges(7,5) * t816 + Ifges(7,6) * t818;
t850 = Ifges(7,3) * t454 / 0.2e1 - t559;
t355 = pkin(3) * t455 + qJ(4) * t454;
t492 = t498 * pkin(2);
t342 = t355 + t492;
t447 = t455 * pkin(9);
t254 = t342 + t447;
t785 = pkin(10) * t455;
t472 = (-pkin(8) - pkin(7)) * t498;
t653 = t797 * pkin(7);
t473 = pkin(8) * t797 + t653;
t855 = t497 * t472 + t796 * t473;
t880 = -t454 * pkin(4) + t855;
t911 = t499 * t880;
t918 = -t454 * pkin(5) + t911;
t114 = (-t254 - t785) * t496 + t918;
t912 = t496 * t880;
t141 = t499 * t254 + t912;
t690 = t455 * t499;
t424 = pkin(10) * t690;
t120 = t424 + t141;
t78 = t114 * t795 - t495 * t120;
t79 = t495 * t114 + t120 * t795;
t870 = t78 * t917 + t79 * t916;
t840 = t850 - t870;
t280 = t355 + t447;
t115 = (-t280 - t785) * t496 + t918;
t162 = t499 * t280 + t912;
t129 = t424 + t162;
t82 = t115 * t795 - t495 * t129;
t83 = t495 * t115 + t129 * t795;
t869 = t82 * t917 + t83 * t916;
t919 = t869 - t850;
t385 = -t796 * t472 + t473 * t497;
t914 = pkin(2) * (t385 * t497 + t796 * t855);
t307 = pkin(4) * t455 + t385;
t288 = t499 * t307;
t485 = -pkin(2) * t797 - pkin(1);
t694 = t455 * qJ(4);
t551 = t485 - t694;
t253 = t454 * t823 + t551;
t588 = pkin(10) * t454 + t253;
t118 = -t496 * t588 + t288;
t113 = pkin(5) * t455 + t118;
t714 = t307 * t496;
t119 = t499 * t588 + t714;
t681 = t495 * t119;
t76 = t113 * t795 - t681;
t88 = t118 * t795 - t681;
t913 = -t76 + t88;
t357 = mrSges(7,1) * t553 + mrSges(7,2) * t552;
t488 = t499 * mrSges(6,2);
t464 = t496 * mrSges(6,1) + t488;
t573 = -t464 - t357;
t570 = -mrSges(5,3) + t573;
t784 = t454 * pkin(3);
t340 = t551 + t784;
t857 = m(5) * t340 - mrSges(5,2) * t454 - mrSges(5,3) * t455;
t439 = Ifges(7,4) * t553;
t360 = -Ifges(7,2) * t552 - t439;
t363 = Ifges(7,1) * t552 - t439;
t910 = t360 + t363;
t226 = -mrSges(7,2) * t455 + t760;
t337 = t464 * t454;
t696 = t454 * t496;
t343 = mrSges(6,1) * t455 - mrSges(6,3) * t696;
t695 = t454 * t499;
t345 = -t455 * mrSges(6,2) + mrSges(6,3) * t695;
t597 = -t674 / 0.2e1;
t668 = t499 * t823;
t786 = pkin(5) * t499;
t619 = -pkin(4) - t786;
t224 = t454 * t619 + t855;
t717 = t224 * t499;
t829 = -m(7) / 0.2e1;
t615 = t795 * t119;
t77 = t495 * t113 + t615;
t87 = -t495 * t118 - t615;
t895 = t87 + t77;
t909 = -((t479 * t696 + t717) * pkin(5) + t895 * t585 + t913 * t366) * t829 + t337 * t923 - t226 * t809 - t343 * t597 - t345 * t668 / 0.2e1 - t823 * t922 + t920;
t652 = t796 * pkin(2);
t484 = -t652 - pkin(3);
t476 = -pkin(9) + t484;
t687 = t476 * t496;
t428 = -t490 + t687;
t429 = (-pkin(10) + t476) * t499;
t329 = t428 * t795 + t495 * t429;
t586 = -t428 * t495 + t795 * t429;
t908 = -t329 * mrSges(7,1) - t586 * mrSges(7,2);
t907 = -t366 * mrSges(7,1) - t585 * mrSges(7,2);
t906 = qJD(6) * t357;
t538 = t455 * t619 - t385;
t859 = t538 * t357;
t886 = t307 * t464;
t905 = t859 / 0.2e1 - t886 / 0.2e1;
t759 = t324 * mrSges(7,2);
t761 = t322 * mrSges(7,1);
t180 = t759 - t761;
t422 = Ifges(6,4) * t695;
t338 = -Ifges(6,2) * t696 + t422;
t770 = Ifges(6,4) * t496;
t469 = Ifges(6,1) * t499 - t770;
t339 = t454 * t469;
t463 = -mrSges(6,2) * t496 + t772;
t765 = Ifges(6,6) * t499;
t766 = Ifges(6,5) * t496;
t465 = -t765 - t766;
t568 = Ifges(6,2) * t499 + t770;
t735 = t455 * Ifges(6,6);
t249 = t454 * t568 + t735;
t672 = t499 * t249;
t251 = Ifges(6,1) * t696 + Ifges(6,5) * t455 + t422;
t678 = t496 * t251;
t788 = pkin(5) * t357;
t854 = t469 - t568;
t899 = t455 / 0.4e1;
t904 = t696 * t788 / 0.2e1 + t180 * t786 / 0.2e1 + t499 * t339 / 0.4e1 - t672 / 0.4e1 - t496 * t338 / 0.4e1 - t678 / 0.4e1 + t465 * t899 + t880 * t463 / 0.2e1 + t854 * t695 / 0.4e1;
t148 = Ifges(7,4) * t325 + Ifges(7,2) * t858 - t454 * Ifges(7,6);
t150 = Ifges(7,1) * t325 + Ifges(7,4) * t858 - t454 * Ifges(7,5);
t250 = -Ifges(6,6) * t454 + t455 * t568;
t769 = Ifges(6,4) * t499;
t569 = Ifges(6,1) * t496 + t769;
t252 = -Ifges(6,5) * t454 + t455 * t569;
t467 = -Ifges(6,2) * t496 + t769;
t600 = t690 / 0.2e1;
t691 = t455 * t496;
t601 = t691 / 0.2e1;
t798 = t499 / 0.2e1;
t799 = -t496 / 0.2e1;
t805 = -t553 / 0.2e1;
t808 = t552 / 0.2e1;
t811 = t363 / 0.2e1;
t767 = Ifges(7,4) * t552;
t361 = -Ifges(7,2) * t553 + t767;
t812 = t361 / 0.2e1;
t526 = t148 * t805 + t150 * t808 + t250 * t799 + t252 * t798 + t858 * t812 + t325 * t811 + t467 * t600 + t469 * t601 - (Ifges(6,5) * t499 + Ifges(7,5) * t552 - Ifges(6,6) * t496 - Ifges(7,6) * t553) * t454 / 0.2e1 + (-Ifges(4,6) + Ifges(5,5)) * t455 + (-Ifges(4,5) + Ifges(5,4)) * t454;
t889 = t855 * mrSges(5,2);
t890 = t855 * mrSges(4,1);
t891 = t385 * mrSges(5,3);
t892 = t385 * mrSges(4,2);
t903 = t526 + t859 + t889 + t892 - t886 - t890 - t891;
t902 = t889 / 0.2e1 - t890 / 0.2e1 - t891 / 0.2e1 + t892 / 0.2e1;
t901 = -m(6) / 0.2e1;
t825 = -mrSges(6,2) / 0.2e1;
t739 = t552 * mrSges(7,3);
t630 = t739 / 0.2e1;
t897 = -t858 / 0.2e1;
t66 = t552 * t77;
t888 = qJ(4) * t307;
t887 = t307 * t880;
t478 = qJ(4) + t789;
t884 = t478 * t307;
t461 = t478 + t491;
t139 = t253 * t499 + t714;
t629 = -t739 / 0.2e1;
t768 = Ifges(7,4) * t324;
t147 = Ifges(7,2) * t322 + t455 * Ifges(7,6) + t768;
t313 = Ifges(7,4) * t322;
t149 = Ifges(7,1) * t324 + t455 * Ifges(7,5) + t313;
t182 = -Ifges(7,2) * t324 + t313;
t183 = Ifges(7,1) * t322 - t768;
t353 = mrSges(7,1) * t552 - mrSges(7,2) * t553;
t359 = -Ifges(7,5) * t553 - Ifges(7,6) * t552;
t362 = -Ifges(7,1) * t553 - t767;
t835 = t359 * t899 + t224 * t353 / 0.2e1 + t910 * t322 / 0.4e1 - (t182 + t149) * t553 / 0.4e1 + (t183 / 0.4e1 - t147 / 0.4e1) * t552 + (t362 / 0.4e1 - t361 / 0.4e1) * t324;
t531 = t77 * t629 + t835;
t605 = t696 / 0.4e1;
t606 = -t696 / 0.4e1;
t623 = mrSges(6,3) * t798;
t738 = t553 * mrSges(7,3);
t626 = -t738 / 0.2e1;
t627 = t738 / 0.2e1;
t726 = t139 * t499;
t512 = t87 * t629 + t88 * t626 + t76 * t627 + t467 * t606 + t139 * t623 - t569 * t605 + t531 - mrSges(6,3) * t726 / 0.2e1 + t904;
t815 = -t586 / 0.2e1;
t817 = -t324 / 0.2e1;
t557 = t322 * t815 + t329 * t817;
t800 = t478 / 0.2e1;
t819 = t228 / 0.2e1;
t828 = m(7) / 0.2e1;
t802 = t461 / 0.2e1;
t813 = t586 / 0.2e1;
t877 = t179 * t802 + t226 * t813;
t883 = t512 + t337 * t800 + ((t461 * t696 + t717) * pkin(5) + t895 * t586 + t913 * t329) * t828 - t329 * t819 + t877 + mrSges(7,3) * t557;
t669 = t499 * t345;
t676 = t496 * t343;
t555 = -t676 / 0.2e1 + t669 / 0.2e1;
t881 = t555 + t922;
t879 = -mrSges(6,3) * t659 - mrSges(4,1) + mrSges(5,2);
t632 = -t766 / 0.2e1;
t878 = Ifges(5,6) + Ifges(4,4) + t632 - t765 / 0.2e1;
t876 = -t478 * t385 + t484 * t855;
t138 = -t253 * t496 + t288;
t181 = -mrSges(7,1) * t858 + mrSges(7,2) * t325;
t227 = mrSges(7,2) * t454 + mrSges(7,3) * t858;
t229 = -mrSges(7,1) * t454 - mrSges(7,3) * t325;
t335 = t463 * t454;
t336 = t463 * t455;
t344 = -mrSges(6,1) * t454 - mrSges(6,3) * t691;
t346 = mrSges(6,2) * t454 + mrSges(6,3) * t690;
t874 = t307 * t335 + (Ifges(7,5) * t817 - Ifges(7,6) * t322 / 0.2e1 + t250 * t798 + t496 * t252 / 0.2e1 + t878 * t454) * t454 + t538 * t180 + t138 * t344 + t139 * t346 + t224 * t181 - t880 * t336 + t340 * (-mrSges(5,2) * t455 + mrSges(5,3) * t454) + t485 * (mrSges(4,1) * t455 - mrSges(4,2) * t454) + t76 * t229 + t77 * t227 + t322 * t148 / 0.2e1 + t147 * t818 + t324 * t150 / 0.2e1 + t149 * t816;
t833 = m(5) / 0.2e1;
t873 = (-pkin(3) * t855 - qJ(4) * t385) * t833;
t701 = t399 * t553;
t702 = t398 * t552;
t562 = -t701 - t702;
t866 = t562 * mrSges(7,3);
t863 = t224 * t538;
t861 = t461 * t538;
t860 = t479 * t538;
t721 = t162 * t496;
t161 = -t280 * t496 + t911;
t722 = t161 * t499;
t564 = t721 + t722;
t831 = m(6) / 0.2e1;
t849 = t552 * t82 + t553 * t83;
t853 = t564 * t831 + t849 * t828;
t723 = t141 * t496;
t140 = -t254 * t496 + t911;
t724 = t140 * t499;
t565 = t723 + t724;
t733 = t79 * t553;
t734 = t78 * t552;
t848 = t733 + t734;
t847 = m(6) / 0.4e1 + m(5) / 0.4e1;
t845 = t565 * t831 + t848 * t828 + t833 * t855;
t844 = mrSges(7,2) * t897 - mrSges(7,1) * t325 / 0.2e1;
t843 = -(m(6) + m(5)) * qJ(4) + t570;
t841 = (-mrSges(4,2) - t570) * t796;
t806 = t553 / 0.2e1;
t807 = -t552 / 0.2e1;
t839 = t226 * t808 + t228 * t805 + (t322 * t807 - t324 * t806) * mrSges(7,3);
t683 = t479 * t181;
t703 = t366 * t227;
t704 = t585 * t229;
t728 = qJ(4) * t336;
t838 = t873 + (-t565 * t823 - t888) * t831 + (t366 * t79 + t585 * t78 + t860) * t828 - t728 / 0.2e1 + t704 / 0.2e1 + t703 / 0.2e1 + t683 / 0.2e1 + (-t724 / 0.2e1 - t723 / 0.2e1) * mrSges(6,3) + (-t733 / 0.2e1 - t734 / 0.2e1) * mrSges(7,3) + t905;
t837 = t738 * t897 + t325 * t630 + t696 * t825 - t759 / 0.2e1 + t761 / 0.2e1 + t880 * t901;
t834 = 2 * qJD(3);
t826 = mrSges(6,1) / 0.2e1;
t822 = mrSges(7,3) * pkin(5);
t821 = -t226 / 0.2e1;
t814 = -t329 / 0.2e1;
t810 = t585 / 0.2e1;
t787 = pkin(5) * t495;
t656 = t552 * t787;
t418 = mrSges(7,3) * t656;
t401 = t418 / 0.2e1;
t803 = t455 / 0.2e1;
t801 = -t467 / 0.2e1;
t793 = m(5) * t855;
t782 = t76 * mrSges(7,2);
t781 = t77 * mrSges(7,1);
t776 = t87 * mrSges(7,1);
t775 = t88 * mrSges(7,2);
t764 = Ifges(6,3) * t454;
t511 = t672 / 0.2e1 + t678 / 0.2e1 - t878 * t455 + (-Ifges(6,3) - Ifges(4,1) + Ifges(5,3) - Ifges(7,3) + Ifges(4,2) - Ifges(5,2)) * t454 + t559;
t2 = (mrSges(4,2) * t492 + t511) * t455 + m(7) * (t76 * t78 + t77 * t79 + t863) + m(6) * (t138 * t140 + t139 * t141 - t887) + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t498 + (m(4) * t485 + mrSges(4,1) * t454) * pkin(2)) * t498 + t140 * t343 + t141 * t345 + t78 * t228 + t79 * t226 + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) - Ifges(3,2)) * t498 + Ifges(3,4) * t797) * t797 + t874 + t857 * t342;
t762 = t2 * qJD(1);
t737 = t454 * mrSges(5,1);
t736 = t455 * mrSges(5,1);
t8 = t511 * t455 + m(7) * (t76 * t82 + t77 * t83 + t863) + m(6) * (t138 * t161 + t139 * t162 - t887) + t161 * t343 + t162 * t345 + t82 * t228 + t83 * t226 + t857 * t355 + t874;
t732 = t8 * qJD(1);
t421 = Ifges(6,5) * t695;
t667 = Ifges(7,5) * t322 - Ifges(7,6) * t324;
t518 = t224 * t179 + (t182 / 0.2e1 + t149 / 0.2e1) * t322 - (t77 * mrSges(7,3) + t147 / 0.2e1 - t183 / 0.2e1) * t324 - t76 * t760 + t667 * t803;
t624 = -t735 / 0.2e1;
t9 = t87 * t228 + t88 * t226 + m(7) * (t76 * t87 + t77 * t88) + t138 * t345 - t139 * t343 + t421 * t803 + t880 * t337 + ((-t138 * mrSges(6,3) + t251 / 0.2e1 + t338 / 0.2e1) * t499 + (-t139 * mrSges(6,3) + t624 + t339 / 0.2e1 - t249 / 0.2e1 + (m(7) * t224 + t180) * pkin(5)) * t496) * t454 + t518;
t729 = t9 * qJD(1);
t727 = qJ(4) * t463;
t14 = t76 * t226 - t77 * t228 + t518;
t725 = t14 * qJD(1);
t507 = t839 + t844;
t716 = t507 * qJD(1);
t27 = -t226 * t858 + t228 * t325 + ((t553 * t76 - t66) * m(7) + m(6) * (t138 * t496 - t726) - t669 + t676 - t857) * t455;
t715 = t27 * qJD(1);
t706 = t586 * t229;
t705 = t329 * t227;
t689 = t461 * t181;
t688 = t461 * t353;
t686 = t478 * t336;
t685 = t478 * t463;
t682 = t479 * t353;
t680 = t495 * t324;
t675 = t496 * t346;
t670 = t499 * t344;
t660 = -t418 / 0.2e1 + t401;
t654 = pkin(5) * t795;
t651 = mrSges(6,3) * t721;
t648 = t478 * t736;
t642 = -t787 / 0.2e1;
t641 = t787 / 0.2e1;
t636 = t476 * t675;
t635 = t476 * t670;
t628 = t66 / 0.2e1;
t625 = t737 / 0.2e1;
t616 = t796 * t478;
t614 = t795 * t322;
t602 = -t694 / 0.2e1;
t596 = t670 / 0.2e1;
t594 = t815 + t813;
t593 = t814 + t329 / 0.2e1;
t591 = t809 + t810;
t590 = -t493 / 0.2e1 - t494 / 0.2e1;
t587 = t659 * t497;
t583 = mrSges(7,3) * t641;
t582 = t454 * t652;
t423 = t553 * t654;
t580 = mrSges(7,3) * t628;
t579 = t654 / 0.2e1;
t574 = t359 + t660;
t572 = (t802 + t915) * t353;
t571 = t614 * t822;
t501 = mrSges(5,1) * t582 / 0.2e1 + t82 * t630 + t484 * t625 + t83 * t627 + t161 * t623 - t636 / 0.2e1 - t705 / 0.2e1 - t689 / 0.2e1 - m(5) * (t876 + t914) / 0.2e1 - (-t335 + t180) * t652 / 0.2e1 + (t343 * t499 + t736) * t643 - t905 - t635 / 0.2e1 + (-t884 + t564 * t476 + (t796 * t880 + (t138 * t499 + t139 * t496) * t497) * pkin(2)) * t901 + t686 / 0.2e1 + t651 / 0.2e1 - t706 / 0.2e1 + (t224 * t652 + t329 * t83 + t398 * t76 + t399 * t77 + t586 * t82 + t861) * t829 + t345 * t581 + t648 / 0.2e1 - t902 + t398 * t820 + t399 * t821;
t3 = mrSges(5,1) * t602 + pkin(3) * t625 + t346 * t597 - t823 * t596 + t501 + t838 + t902;
t48 = m(7) * (t329 * t399 + t398 * t586 + t461 * t652) + t866 + (m(6) * (t476 * t587 + t616) + m(5) * (t484 * t497 + t616) + t841) * pkin(2) + t879 * t789;
t567 = -t3 * qJD(1) + t48 * qJD(2);
t510 = t531 + (t628 + t557) * mrSges(7,3) + t228 * t814 + t877;
t11 = t510 + t840;
t532 = -(-t362 / 0.2e1 + t812) * t552 - (t811 + t360 / 0.2e1) * t553;
t46 = t532 + t688;
t566 = t11 * qJD(1) + t46 * qJD(2);
t188 = m(7) * t461 + 0.4e1 * t478 * t847 - t570;
t554 = t596 + t675 / 0.2e1;
t528 = t227 * t806 + t229 * t808 + t554;
t503 = t454 * t634 + t528 + t837;
t530 = -pkin(5) * t695 + t880;
t506 = -t793 / 0.2e1 + (t325 * t586 - t329 * t858 + t530) * t829;
t20 = t503 + t506 + t845;
t561 = qJD(1) * t20 - qJD(2) * t188;
t430 = t461 * t786;
t515 = (-t469 / 0.2e1 + t568 / 0.2e1) * t496 + (t788 - t569 / 0.2e1 + t801) * t499 + t532;
t37 = m(7) * t430 + t515 + t685 + t688;
t523 = t140 * t826 + t141 * t825 + (t495 * t79 + t78 * t795) * t657;
t5 = t764 / 0.2e1 + t499 * t624 + t455 * t632 + t227 * t642 - t229 * t654 / 0.2e1 - t523 + t881 * t476 + t840 + t883;
t546 = t5 * qJD(1) + t37 * qJD(2);
t514 = t226 * t810 + t531 + t580 + t920;
t13 = t514 - t919;
t522 = t572 + t532;
t35 = t522 - t856;
t47 = t532 + t682;
t545 = t13 * qJD(1) + t35 * qJD(2) + t47 * qJD(3);
t520 = mrSges(6,2) * t600 + (mrSges(6,1) + (t495 ^ 2 + t795 ^ 2) * t827) * t601;
t536 = (t87 * t552 + t553 * t913 + t66) * t828;
t16 = t536 + t507 - t520 + t881;
t544 = t16 * qJD(1);
t542 = Ifges(6,5) * t601 + Ifges(6,6) * t600 + t227 * t641 + t229 * t579 - t764 / 0.2e1 - t850;
t104 = (t702 / 0.2e1 + t701 / 0.2e1 - t479) * m(7) + (t829 + (-0.1e1 / 0.2e1 - t590) * m(6)) * t789 + t843;
t225 = m(7) * t479 - t843;
t509 = m(7) * (t325 * t585 - t366 * t858 + t530);
t24 = t503 - t509 / 0.2e1 + t853;
t541 = -qJD(1) * t24 - qJD(2) * t104 + qJD(3) * t225;
t462 = t479 * t786;
t524 = (t430 + t462) * t828;
t29 = -t569 * t798 + t499 * t801 + t682 / 0.2e1 + t685 / 0.2e1 + t688 / 0.2e1 + t727 / 0.2e1 + t361 * t807 + t362 * t808 + t357 * t786 + t524 + t910 * t805 + t854 * t799 + (t627 + t626) * (t585 + t586) - t921;
t39 = m(7) * t462 + t515 + t682 + t727;
t513 = t161 * t826 + t162 * t825 + (t495 * t83 + t795 * t82) * t657 + t542 + t869;
t6 = t467 * t605 - t569 * t606 + t76 * t626 + t88 * t627 + t87 * t630 + t513 + t580 - t835 - t904 - t909;
t533 = -t6 * qJD(1) + t29 * qJD(2) + t39 * qJD(3);
t529 = t590 * t771 + t555;
t527 = mrSges(7,3) * t423 + t359 - t418 + t465;
t18 = (-t88 / 0.2e1 + t76 / 0.2e1) * mrSges(7,2) + (t87 / 0.2e1 + t77 / 0.2e1) * mrSges(7,1) + (t795 * t821 + t495 * t819 + (t680 / 0.2e1 + t614 / 0.2e1) * mrSges(7,3)) * pkin(5);
t458 = (mrSges(7,1) * t495 + mrSges(7,2) * t795) * pkin(5);
t52 = mrSges(7,1) * t593 + mrSges(7,2) * t594 + t660;
t86 = mrSges(7,2) * t591 - t552 * t583 + t401;
t525 = t18 * qJD(1) - t52 * qJD(2) + t86 * qJD(3) + t458 * qJD(5);
t521 = t839 - t844;
t517 = -mrSges(6,1) * t695 / 0.2e1 + t528 - t737 - t837;
t486 = qJ(4) * t652;
t450 = t458 * qJD(6);
t105 = m(7) * t491 + t562 * t829 + (t659 * t831 + t833) * t789 - t570 + (m(7) / 0.4e1 + t847) * (0.4e1 * qJ(4) + 0.2e1 * t789);
t62 = t574 + t907;
t49 = t574 + t908;
t36 = t522 + t856;
t30 = t524 + (t800 + t923) * t463 + t572 + (-(-t591 + t594) * t553 + t593 * t552) * mrSges(7,3) + t515 + t921;
t21 = t509 / 0.2e1 + t793 + t517 + t853;
t19 = -t506 + t517 + t845;
t17 = -t781 / 0.2e1 - t782 / 0.2e1 + t226 * t579 - t324 * t583 + t228 * t642 - t571 / 0.2e1 + t776 / 0.2e1 - t775 / 0.2e1 + t667;
t15 = t536 + t520 + t521 + t529;
t12 = t514 + t919;
t10 = t510 - t840;
t7 = t512 + t513 + t909;
t4 = t476 * t529 + t523 + t542 + t870 + t883;
t1 = -(-mrSges(4,2) / 0.2e1 + mrSges(5,3) / 0.2e1) * t385 + (mrSges(5,2) / 0.2e1 - mrSges(4,1) / 0.2e1) * t855 + (t602 + t784 / 0.2e1) * mrSges(5,1) + t526 - t501 - t554 * t823 + t838;
t22 = [qJD(2) * t2 + qJD(3) * t8 + qJD(4) * t27 + qJD(5) * t9 + qJD(6) * t14, t762 + (-m(4) * t914 + t636 + t705 + t689 + (mrSges(3,2) * pkin(7) - Ifges(3,6)) * t498 + m(6) * (t476 * t565 - t884) + m(5) * t876 + (-t455 * t789 + t582) * mrSges(4,3) + t635 - t686 + t706 + m(7) * (t329 * t79 + t586 * t78 + t861) - t848 * mrSges(7,3) - t565 * mrSges(6,3) + Ifges(3,5) * t797 - t484 * t737 - t648 + t903 - mrSges(3,1) * t653) * qJD(2) + t1 * qJD(3) + t19 * qJD(4) + t4 * qJD(5) + t10 * qJD(6), t732 + t1 * qJD(2) + (-mrSges(5,1) * t694 - mrSges(6,3) * t722 - mrSges(7,3) * t849 + pkin(3) * t737 - t344 * t668 - t346 * t674 - t651 + t683 + t703 + t704 - t728 + t903) * qJD(3) + t21 * qJD(4) + t7 * qJD(5) + t12 * qJD(6) + ((-t564 * t823 - t888) * t831 + (t366 * t83 + t585 * t82 + t860) * t828 + t873) * t834, t19 * qJD(2) + t21 * qJD(3) + t15 * qJD(5) + qJD(6) * t521 + t715, t729 + t4 * qJD(2) + t7 * qJD(3) + t15 * qJD(4) + (t421 - Ifges(6,6) * t696 + t776 - t775 + (t495 * t88 + t795 * t87) * t827 - t680 * t822 - t571 - t138 * mrSges(6,2) - t139 * mrSges(6,1) + t667) * qJD(5) + t17 * qJD(6), t725 + t10 * qJD(2) + t12 * qJD(3) + t521 * qJD(4) + t17 * qJD(5) + (t667 - t781 - t782) * qJD(6); -qJD(3) * t3 - qJD(4) * t20 + qJD(5) * t5 + qJD(6) * t11 - t762, qJD(3) * t48 + qJD(4) * t188 + qJD(5) * t37 + qJD(6) * t46, t105 * qJD(4) + t30 * qJD(5) + t36 * qJD(6) + ((t366 * t399 + t398 * t585 + t479 * t652) * t828 + (-pkin(2) * t587 * t823 + t486) * t831 + (-pkin(3) * t789 + t486) * t833) * t834 + t567 + (t866 + (t497 * t879 + t841) * pkin(2)) * qJD(3), qJD(3) * t105 - t561, t30 * qJD(3) + ((-t329 * t795 + t495 * t586) * t827 - t476 * t488 - mrSges(6,1) * t687 + t527 + t908) * qJD(5) + t49 * qJD(6) + t546, t36 * qJD(3) + t49 * qJD(5) + (t359 + t908) * qJD(6) + t566; qJD(2) * t3 - qJD(4) * t24 - qJD(5) * t6 + qJD(6) * t13 - t732, -qJD(4) * t104 + qJD(5) * t29 + qJD(6) * t35 - t567, qJD(4) * t225 + qJD(5) * t39 + qJD(6) * t47, t541 ((-t366 * t795 + t495 * t585) * t827 + mrSges(6,2) * t668 + mrSges(6,1) * t674 + t527 + t907) * qJD(5) + t62 * qJD(6) + t533, t62 * qJD(5) + (t359 + t907) * qJD(6) + t545; qJD(2) * t20 + qJD(3) * t24 + qJD(5) * t16 + qJD(6) * t507 - t715, qJD(3) * t104 + t561, -t541, 0 (m(7) * (-t423 + t656) + t573) * qJD(5) - t906 + t544, -qJD(5) * t357 + t716 - t906; -qJD(2) * t5 + qJD(3) * t6 - qJD(4) * t16 - qJD(6) * t18 - t729, -qJD(3) * t29 + qJD(6) * t52 - t546, -qJD(6) * t86 - t533, -t544, -t450, -t450 - t525; -qJD(2) * t11 - qJD(3) * t13 - qJD(4) * t507 + qJD(5) * t18 - t725, -qJD(3) * t35 - qJD(5) * t52 - t566, qJD(5) * t86 - t545, -t716, t525, 0;];
Cq  = t22;
