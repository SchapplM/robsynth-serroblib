% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:19
% EndTime: 2019-03-09 00:13:51
% DurationCPUTime: 17.34s
% Computational Cost: add. (21832->887), mult. (51634->1168), div. (0->0), fcn. (53692->10), ass. (0->449)
t918 = Ifges(7,4) + Ifges(6,5);
t924 = Ifges(6,3) + Ifges(7,2);
t532 = cos(qJ(4));
t508 = -pkin(4) * t532 - pkin(3);
t786 = cos(qJ(5));
t648 = t786 * t532;
t528 = sin(qJ(5));
t529 = sin(qJ(4));
t697 = t528 * t529;
t577 = t648 - t697;
t649 = t786 * t529;
t578 = t528 * t532 + t649;
t594 = -pkin(5) * t577 - qJ(6) * t578;
t306 = t508 + t594;
t339 = -mrSges(7,1) * t577 - mrSges(7,3) * t578;
t905 = m(7) * t306 + t339;
t530 = sin(qJ(3));
t533 = cos(qJ(3));
t492 = pkin(3) * t530 - pkin(9) * t533;
t695 = t529 * t530;
t412 = pkin(8) * t695 + t532 * t492;
t691 = t532 * t533;
t323 = pkin(4) * t530 - pkin(10) * t691 + t412;
t693 = t530 * t532;
t413 = -pkin(8) * t693 + t529 * t492;
t694 = t529 * t533;
t366 = -pkin(10) * t694 + t413;
t187 = t528 * t323 + t786 * t366;
t159 = qJ(6) * t530 + t187;
t184 = t323 * t786 - t528 * t366;
t160 = -t530 * pkin(5) - t184;
t434 = t578 * t533;
t376 = -mrSges(7,2) * t434 + mrSges(7,3) * t530;
t379 = -mrSges(6,2) * t530 - mrSges(6,3) * t434;
t435 = t577 * t533;
t382 = mrSges(6,1) * t530 - mrSges(6,3) * t435;
t782 = pkin(4) * t528;
t503 = qJ(6) + t782;
t672 = t786 * pkin(4);
t507 = -t672 - pkin(5);
t607 = t672 / 0.2e1;
t723 = t533 * Ifges(5,5);
t651 = t723 / 0.2e1;
t668 = t782 / 0.2e1;
t846 = m(6) * pkin(4);
t676 = t846 / 0.2e1;
t790 = t530 / 0.2e1;
t793 = t503 / 0.2e1;
t725 = t530 * mrSges(7,1);
t731 = t435 * mrSges(7,2);
t383 = -t725 + t731;
t814 = t383 / 0.2e1;
t843 = -mrSges(5,2) / 0.2e1;
t844 = mrSges(5,1) / 0.2e1;
t847 = m(7) / 0.2e1;
t809 = t435 / 0.2e1;
t909 = -mrSges(7,3) / 0.2e1;
t910 = mrSges(6,2) / 0.2e1;
t922 = mrSges(7,1) / 0.2e1;
t855 = t187 * t910 - t184 * mrSges(6,1) / 0.2e1 + t160 * t922 + t159 * t909;
t921 = t434 / 0.2e1;
t853 = Ifges(7,6) * t921 - Ifges(6,6) * t434 / 0.2e1 + t918 * t809 + t924 * t790 - t855;
t923 = (t159 * t503 + t160 * t507) * t847 + Ifges(5,3) * t790 + t412 * t844 + t413 * t843 + t376 * t793 + t507 * t814 + (t184 * t786 + t187 * t528) * t676 + t532 * t651 - Ifges(5,6) * t694 / 0.2e1 + t379 * t668 + t382 * t607 + t853;
t523 = t529 ^ 2;
t525 = t532 ^ 2;
t869 = t525 + t523;
t919 = pkin(9) * t869;
t888 = mrSges(7,1) + mrSges(6,1);
t917 = Ifges(6,6) - Ifges(7,6);
t527 = sin(pkin(6));
t531 = sin(qJ(2));
t701 = t527 * t531;
t720 = cos(pkin(6));
t448 = t530 * t701 - t533 * t720;
t266 = t577 * t448;
t662 = t909 + t910;
t916 = t662 * t266;
t730 = t577 * mrSges(7,2);
t653 = -t730 / 0.2e1;
t654 = t730 / 0.2e1;
t325 = t654 + t653;
t915 = qJD(3) * t325;
t914 = qJD(6) * t325;
t519 = t530 * pkin(8);
t477 = pkin(4) * t695 + t519;
t433 = t578 * t530;
t432 = t528 * t695 - t530 * t648;
t719 = qJ(6) * t432;
t595 = pkin(5) * t433 + t719;
t234 = t595 + t477;
t290 = -mrSges(6,1) * t432 - mrSges(6,2) * t433;
t521 = t532 * pkin(9);
t679 = t532 * pkin(10) + t521;
t840 = -pkin(10) - pkin(9);
t372 = t528 * t679 - t840 * t649;
t603 = t918 * t577 - t917 * t578;
t515 = t533 * mrSges(7,3);
t733 = t433 * mrSges(7,2);
t377 = -t515 - t733;
t732 = t433 * mrSges(6,3);
t378 = mrSges(6,2) * t533 - t732;
t622 = t378 / 0.2e1 + t377 / 0.2e1;
t787 = -t533 / 0.4e1;
t792 = t508 / 0.2e1;
t338 = mrSges(6,1) * t578 + mrSges(6,2) * t577;
t824 = t338 / 0.2e1;
t337 = mrSges(7,1) * t578 - mrSges(7,3) * t577;
t825 = t337 / 0.2e1;
t289 = -mrSges(7,1) * t432 + mrSges(7,3) * t433;
t828 = t289 / 0.2e1;
t913 = t234 * t825 + t290 * t792 + t306 * t828 - t622 * t372 + t477 * t824 + t603 * t787;
t769 = Ifges(7,5) * t577;
t341 = Ifges(7,3) * t578 + t769;
t463 = Ifges(6,4) * t577;
t343 = -Ifges(6,2) * t578 + t463;
t460 = Ifges(7,5) * t578;
t345 = Ifges(7,1) * t577 + t460;
t346 = Ifges(7,1) * t578 - t769;
t773 = Ifges(6,4) * t578;
t347 = Ifges(6,1) * t577 - t773;
t348 = Ifges(6,1) * t578 + t463;
t342 = -Ifges(7,3) * t577 + t460;
t344 = Ifges(6,2) * t577 + t773;
t859 = t342 / 0.2e1 - t344 / 0.2e1;
t912 = -(-t348 / 0.2e1 - t343 / 0.2e1 - t346 / 0.2e1 + t341 / 0.2e1) * t577 - (-t347 / 0.2e1 - t345 / 0.2e1 - t859) * t578 + t306 * t337 + t508 * t338;
t481 = -pkin(3) * t533 - pkin(9) * t530 - pkin(2);
t673 = pkin(8) * t691;
t360 = t673 + (-pkin(10) * t530 + t481) * t529;
t650 = t786 * t360;
t465 = t532 * t481;
t601 = -pkin(10) * t693 + t465;
t321 = (-pkin(8) * t529 - pkin(4)) * t533 + t601;
t699 = t528 * t321;
t174 = t650 + t699;
t690 = t533 * qJ(6);
t144 = t174 - t690;
t449 = t530 * t720 + t533 * t701;
t534 = cos(qJ(2));
t700 = t527 * t534;
t349 = t449 * t532 - t529 * t700;
t580 = -t449 * t529 - t532 * t700;
t186 = t349 * t528 - t786 * t580;
t288 = -pkin(5) * t432 + qJ(6) * t433;
t675 = pkin(4) * t693;
t256 = t288 + t675;
t472 = mrSges(5,2) * t533 - mrSges(5,3) * t695;
t555 = t349 * t786 + t528 * t580;
t566 = t580 * t529;
t674 = pkin(8) * t694;
t359 = t601 - t674;
t195 = t359 * t528 + t650;
t698 = t528 * t360;
t196 = t359 * t786 - t698;
t616 = t186 * t195 + t196 * t555;
t906 = t432 * mrSges(6,3);
t380 = -mrSges(6,1) * t533 + t906;
t734 = t432 * mrSges(7,2);
t381 = mrSges(7,1) * t533 - t734;
t621 = -t380 / 0.2e1 + t381 / 0.2e1;
t713 = t349 * t532;
t777 = mrSges(5,3) * t530;
t474 = -mrSges(5,1) * t533 - mrSges(5,3) * t693;
t797 = -t474 / 0.2e1;
t811 = -t433 / 0.2e1;
t812 = t432 / 0.2e1;
t849 = m(6) / 0.2e1;
t173 = t321 * t786 - t698;
t882 = t173 * t555;
t151 = t533 * pkin(5) - t173;
t884 = t151 * t555;
t886 = mrSges(6,3) + mrSges(7,2);
t891 = t580 / 0.2e1;
t911 = t886 * (t186 * t811 + t555 * t812) + (-t174 * t186 + t448 * t675 + t616 - t882) * t849 + (-t144 * t186 + t256 * t448 + t616 + t884) * t847 + t349 * t797 + t472 * t891 + (-t713 / 0.2e1 + t566 / 0.2e1) * t777 + t621 * t555 - t622 * t186;
t265 = t578 * t448;
t907 = t265 / 0.2e1;
t887 = mrSges(6,2) - mrSges(7,3);
t605 = t888 * t907 + t916;
t706 = t448 * t529;
t634 = t706 / 0.2e1;
t705 = t448 * t532;
t903 = (-t265 * t507 - t266 * t503) * t847 + (t265 * t786 - t266 * t528) * t676 + mrSges(5,1) * t634 + mrSges(5,2) * t705 / 0.2e1 + t605;
t688 = t533 * t534;
t401 = (-t529 * t688 + t531 * t532) * t527;
t402 = (t529 * t531 + t532 * t688) * t527;
t236 = -t401 * t786 + t402 * t528;
t237 = t528 * t401 + t402 * t786;
t606 = -t888 * t236 / 0.2e1 + (mrSges(7,3) / 0.2e1 - mrSges(6,2) / 0.2e1) * t237;
t902 = (t236 * t507 + t237 * t503) * t847 + t401 * t844 + t402 * t843 + (-t236 * t786 + t237 * t528) * t676 + t606;
t514 = m(7) * qJ(6) + mrSges(7,3);
t901 = qJD(5) * t514;
t900 = t514 * qJD(6);
t336 = pkin(5) * t578 - qJ(6) * t577;
t781 = pkin(4) * t529;
t310 = t336 + t781;
t862 = -mrSges(5,1) * t532 + mrSges(5,2) * t529;
t453 = t862 * t530;
t858 = t679 * t786 + t840 * t697;
t614 = t195 * t372 + t196 * t858;
t823 = t339 / 0.2e1;
t291 = mrSges(7,1) * t433 + mrSges(7,3) * t432;
t827 = t291 / 0.2e1;
t899 = (-t144 * t372 + t151 * t858 + t234 * t310 + t256 * t306 + t614) * t847 + pkin(3) * t453 / 0.2e1 + t256 * t823 + t310 * t827 + t621 * t858 + t913;
t897 = -pkin(5) * t555 - qJ(6) * t186;
t404 = -t733 / 0.2e1;
t894 = 0.2e1 * t404;
t893 = m(7) + m(6);
t892 = t555 / 0.2e1;
t516 = Ifges(5,5) * t532;
t766 = Ifges(5,6) * t529;
t885 = Ifges(4,4) - t516 / 0.2e1 + t766 / 0.2e1;
t883 = t160 * t847;
t485 = t530 * mrSges(4,1) + t533 * mrSges(4,2);
t881 = -t485 * t700 / 0.2e1;
t880 = t507 * t555;
t340 = -mrSges(6,1) * t577 + mrSges(6,2) * t578;
t875 = t339 + t340;
t874 = t348 + t346;
t873 = t377 + t378;
t872 = t379 + t376;
t871 = -t380 + t381;
t517 = Ifges(5,4) * t532;
t870 = -Ifges(5,2) * t529 + t517;
t488 = Ifges(5,1) * t529 + t517;
t868 = t372 * t811 + t812 * t858;
t775 = Ifges(5,4) * t529;
t486 = Ifges(5,2) * t532 + t775;
t456 = t530 * t486;
t865 = pkin(9) * t797 - t456 / 0.4e1;
t604 = t917 * t432 - t918 * t433;
t864 = -t412 * t529 + t413 * t532;
t292 = mrSges(6,1) * t433 - mrSges(6,2) * t432;
t484 = mrSges(5,1) * t529 + mrSges(5,2) * t532;
t454 = t484 * t530;
t861 = t291 + t292 + t454;
t457 = t530 * t488;
t798 = -t472 / 0.2e1;
t860 = pkin(9) * t798 + pkin(4) * t292 / 0.2e1 - t457 / 0.4e1;
t852 = 0.2e1 * m(7);
t851 = 0.2e1 * qJD(3);
t850 = m(5) / 0.2e1;
t848 = -m(7) / 0.2e1;
t845 = m(7) * pkin(4);
t842 = -mrSges(7,2) / 0.2e1;
t841 = -mrSges(6,3) / 0.2e1;
t839 = -qJ(6) / 0.2e1;
t838 = t144 / 0.2e1;
t837 = -t151 / 0.2e1;
t836 = t173 / 0.2e1;
t835 = t174 / 0.2e1;
t832 = -t195 / 0.2e1;
t831 = -t196 / 0.2e1;
t830 = Ifges(7,5) * t809 + Ifges(7,6) * t790 + Ifges(7,3) * t921;
t774 = Ifges(6,4) * t432;
t269 = -Ifges(6,2) * t433 - t533 * Ifges(6,6) - t774;
t829 = t269 / 0.4e1;
t770 = Ifges(7,5) * t433;
t295 = -Ifges(7,3) * t432 - t770;
t826 = t295 / 0.4e1;
t822 = t340 / 0.2e1;
t821 = t341 / 0.4e1;
t820 = t344 / 0.4e1;
t819 = t858 / 0.2e1;
t808 = t448 / 0.2e1;
t806 = -t453 / 0.2e1;
t803 = t577 / 0.2e1;
t801 = t578 / 0.2e1;
t799 = -t578 / 0.2e1;
t796 = t484 / 0.2e1;
t795 = -t486 / 0.4e1;
t489 = Ifges(5,1) * t532 - t775;
t794 = t489 / 0.4e1;
t791 = -t529 / 0.2e1;
t789 = t532 / 0.2e1;
t788 = -t533 / 0.2e1;
t785 = m(7) * t555;
t784 = m(7) * t195;
t783 = m(7) * t858;
t522 = t533 * pkin(8);
t779 = mrSges(4,2) * t530;
t760 = t173 * mrSges(6,2);
t759 = t173 * mrSges(7,3);
t758 = t174 * mrSges(6,1);
t757 = t174 * mrSges(7,1);
t755 = t555 * mrSges(6,1);
t754 = t555 * mrSges(7,1);
t752 = t186 * mrSges(6,2);
t751 = t186 * mrSges(7,3);
t750 = t195 * mrSges(6,1);
t749 = t195 * mrSges(7,1);
t748 = t196 * mrSges(6,2);
t747 = t196 * mrSges(7,3);
t738 = t858 * mrSges(6,1);
t737 = t858 * mrSges(7,1);
t736 = t372 * mrSges(6,2);
t735 = t372 * mrSges(7,3);
t729 = t577 * mrSges(6,3);
t728 = t578 * mrSges(7,2);
t727 = t578 * mrSges(6,3);
t722 = t533 * Ifges(5,6);
t721 = -mrSges(4,1) + t862;
t318 = t448 * t449;
t24 = m(5) * (t318 + (t566 - t713) * t448) + t893 * (-t186 * t265 - t266 * t555 + t318);
t718 = t24 * qJD(1);
t659 = t530 * t700;
t375 = t448 * t659;
t25 = m(5) * (t349 * t402 + t401 * t580 + t375) + m(4) * (t448 * t530 + t449 * t533 - t701) * t700 + t893 * (t186 * t236 + t237 * t555 + t375);
t717 = t25 * qJD(1);
t711 = t401 * t529;
t710 = t402 * t532;
t707 = t432 * t448;
t704 = t477 * t529;
t703 = t503 * t578;
t702 = t507 * t577;
t427 = t530 * t870 - t722;
t696 = t529 * t427;
t429 = t489 * t530 - t723;
t692 = t532 * t429;
t689 = t533 * t555;
t687 = -t144 + t174;
t682 = 0.2e1 * t654;
t478 = pkin(4) * t694 + t522;
t677 = mrSges(6,3) * t782;
t671 = t503 * t734;
t669 = -t782 / 0.2e1;
t667 = pkin(8) * t796;
t664 = mrSges(6,1) / 0.2e1 + t922;
t663 = t841 + t842;
t661 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t660 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t652 = -t728 / 0.2e1;
t633 = t448 * t799;
t629 = -t555 / 0.2e1 + t892;
t417 = Ifges(7,5) * t432;
t267 = -t533 * Ifges(7,6) + Ifges(7,3) * t433 - t417;
t628 = t267 / 0.2e1 - t269 / 0.2e1;
t271 = -Ifges(7,1) * t432 - Ifges(7,4) * t533 + t770;
t420 = Ifges(6,4) * t433;
t273 = -Ifges(6,1) * t432 - Ifges(6,5) * t533 - t420;
t627 = t271 / 0.2e1 + t273 / 0.2e1;
t626 = t828 + t290 / 0.2e1;
t625 = t825 + t824;
t620 = -t503 + t782;
t619 = t516 - t766;
t615 = t432 * t677;
t613 = -t265 * t372 - t266 * t858;
t610 = mrSges(6,3) * t672;
t598 = t433 * t610;
t235 = pkin(5) * t434 - qJ(6) * t435 + t478;
t293 = mrSges(7,1) * t434 - mrSges(7,3) * t435;
t294 = mrSges(6,1) * t434 + mrSges(6,2) * t435;
t410 = t465 - t674;
t411 = t481 * t529 + t673;
t455 = t484 * t533;
t473 = -mrSges(5,2) * t530 - mrSges(5,3) * t694;
t475 = mrSges(5,1) * t530 - mrSges(5,3) * t691;
t535 = (t413 * t349 + t412 * t580 + t449 * t519 + (t410 * t529 - t411 * t532 + t522) * t448) * t850 + (t187 * t555 + t448 * t478 + t449 * t477) * t849 + (t159 * t555 + t234 * t449 + t235 * t448) * t847 + t349 * t473 / 0.2e1 + t475 * t891 + t474 * t634 + t705 * t798 + t881 + t872 * t892 - (t144 * t847 + t174 * t849 + t622) * t266 - (t151 * t847 - t173 * t849 + t621) * t265 + (-t184 * t849 + t883 - t382 / 0.2e1 + t814) * t186 + (t293 + t294 + t455) * t808 + t861 * t449 / 0.2e1;
t589 = t372 * t236 + t237 * t858;
t537 = (-pkin(3) * t659 + (t710 - t711) * pkin(9)) * t850 + (t508 * t659 + t589) * t849 + (t306 * t659 + t589) * t847 + t881 + (-t711 / 0.2e1 + t710 / 0.2e1) * mrSges(5,3) + (t862 + t875) * t659 / 0.2e1 + t886 * (t236 * t801 + t237 * t803);
t6 = -t537 + t535;
t270 = Ifges(6,4) * t435 - Ifges(6,2) * t434 + t530 * Ifges(6,6);
t272 = Ifges(7,1) * t435 + t530 * Ifges(7,4) + Ifges(7,5) * t434;
t274 = Ifges(6,1) * t435 - Ifges(6,4) * t434 + t530 * Ifges(6,5);
t428 = Ifges(5,6) * t530 + t533 * t870;
t430 = t530 * Ifges(5,5) + t489 * t533;
t561 = t434 * t660 + t435 * t661;
t7 = m(5) * (t410 * t412 + t411 * t413) + (-t272 / 0.2e1 - t274 / 0.2e1) * t432 + (-t270 / 0.2e1 + t830) * t433 + t627 * t435 + t628 * t434 + m(6) * (t173 * t184 + t174 * t187 + t477 * t478) + m(7) * (t144 * t159 + t151 * t160 + t234 * t235) + (t430 * t789 + t428 * t791 + pkin(8) * t455 - t885 * t530 - t660 * t433 + t661 * t432 + (m(5) * pkin(8) ^ 2 + Ifges(4,1) - Ifges(4,2) - Ifges(5,3) - t924) * t533) * t530 + t477 * t294 + t478 * t292 - pkin(2) * t485 + t413 * t472 + t411 * t473 + t412 * t474 + t410 * t475 + t159 * t377 + t187 * t378 + t174 * t379 + t184 * t380 + t160 * t381 + t173 * t382 + t151 * t383 + t144 * t376 + t235 * t291 + t234 * t293 + (t692 / 0.2e1 - t696 / 0.2e1 + pkin(8) * t454 + t561 + t885 * t533) * t533;
t593 = t6 * qJD(1) + t7 * qJD(2);
t296 = Ifges(6,2) * t432 - t420;
t297 = -Ifges(7,1) * t433 - t417;
t298 = -Ifges(6,1) * t433 + t774;
t542 = t144 * t734 + t173 * t732 + t234 * t289 + t477 * t290 + (-t298 / 0.2e1 + t174 * mrSges(6,3) - t297 / 0.2e1 - t628) * t432 + (t295 / 0.2e1 - t151 * mrSges(7,2) - t296 / 0.2e1 - t627) * t433 + t604 * t788;
t10 = m(6) * (-t173 * t195 + t174 * t196) + m(7) * (t144 * t196 + t151 * t195 + t234 * t256) + t410 * t472 - t411 * t474 + t256 * t291 + t542 + (-pkin(8) * t453 + (t410 * mrSges(5,3) + t651 + t456 / 0.2e1 - t429 / 0.2e1) * t529 + (-t411 * mrSges(5,3) + t722 / 0.2e1 - t457 / 0.2e1 - t427 / 0.2e1 + (m(6) * t477 + t292) * pkin(4)) * t532) * t530 + t871 * t195 + t873 * t196;
t9 = t448 * t806 + (t290 + t289) * t808 - t902 + t911;
t592 = t9 * qJD(1) + t10 * qJD(2);
t549 = t626 * t448 + (-t432 * t663 + t621) * t555;
t564 = t433 * t663 - t622;
t583 = m(7) * (-pkin(5) * t236 + qJ(6) * t237);
t585 = t288 * t448 + t882 + t884;
t12 = t662 * t237 + t664 * t236 + t564 * t186 + (t186 * t687 + t585) * t847 - t583 / 0.2e1 + t549;
t13 = t542 + (m(7) * t234 + t291) * t288 + (m(7) * t151 + t871) * t174 + (m(7) * t144 + t873) * t173;
t591 = t12 * qJD(1) + t13 * qJD(2);
t48 = m(7) * (-t144 * t533 + t234 * t432) + t432 * t291 - t533 * t377;
t57 = (t236 / 0.4e1 + t689 / 0.4e1 - t707 / 0.4e1) * t852;
t590 = qJD(1) * t57 - qJD(2) * t48;
t587 = t883 - t725 / 0.2e1;
t584 = m(7) * t897;
t582 = m(7) * (pkin(5) * t265 - qJ(6) * t266);
t581 = m(7) * (-pkin(5) * t858 - qJ(6) * t372);
t576 = t829 - t298 / 0.4e1 - t297 / 0.4e1 - t267 / 0.4e1;
t575 = -t296 / 0.4e1 - t273 / 0.4e1 - t271 / 0.4e1 + t826;
t574 = t821 - t343 / 0.4e1 - t346 / 0.4e1 - t348 / 0.4e1;
t573 = -t342 / 0.4e1 + t820 - t345 / 0.4e1 - t347 / 0.4e1;
t546 = pkin(4) * t706 * t849 + t310 * t448 * t847;
t15 = t448 * t796 + t546 + (t338 + t337) * t808 - t903 + t886 * (t799 + t801) * t555;
t563 = -t173 * t858 - t174 * t372 + t614;
t2 = -t923 + t899 - (t347 + t345 + t342) * t432 / 0.4e1 + t530 * t667 - t696 / 0.4e1 + t692 / 0.4e1 + (t273 + t271 + t296) * t577 / 0.4e1 - t577 * t826 + (t298 + t297 + t267) * t578 / 0.4e1 - t578 * t829 - (t488 + t870) * t695 / 0.4e1 + ((t508 * t693 + t704) * pkin(4) + t563) * t849 - t727 * t835 - t729 * t836 + t432 * t820 + t433 * t821 + t675 * t822 + t619 * t787 + t693 * t794 + t693 * t795 - t777 * t919 / 0.2e1 + t886 * (t195 * t801 + t196 * t803 + t868) + t144 * t652 + t151 * t654 + t860 * t529 + t865 * t532 - (t343 + t874) * t433 / 0.4e1;
t20 = -pkin(3) * t484 + (t870 / 0.2e1 + t488 / 0.2e1) * t532 + (-t486 / 0.2e1 + t489 / 0.2e1 + pkin(4) * t340) * t529 + m(6) * t508 * t781 + t912 + t905 * t310;
t568 = t15 * qJD(1) + t2 * qJD(2) + t20 * qJD(3);
t544 = (t336 * t847 + t625) * t448;
t17 = -t916 - t664 * t265 - t582 / 0.2e1 + t544;
t23 = t336 * t905 + t912;
t536 = -t576 * t578 - t575 * t577 + (-(-t174 / 0.2e1 + t838) * t578 - (t837 - t173 / 0.2e1) * t577 + t868) * mrSges(7,2) + (t372 * t841 + t574) * t433 + t573 * t432 + (t234 * t336 + t288 * t306 + t372 * t687) * t847 + t288 * t823 + t336 * t827 + (t621 + t906 / 0.2e1 + (t151 + t173) * t847) * t858 + t913;
t554 = (-pkin(5) * t160 + qJ(6) * t159) * t848 + pkin(5) * t814 + t376 * t839;
t4 = (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t530 + t536 + t554 + t561 + t855;
t567 = t17 * qJD(1) + t4 * qJD(2) + t23 * qJD(3);
t130 = (t633 + t907) * m(7);
t552 = (-t234 * t578 + t306 * t432 - t533 * t858) * t848 - t432 * t339 / 0.2e1 + t291 * t801;
t44 = (-t577 * t788 + t809) * mrSges(7,2) + t552 + t587;
t84 = t905 * t578;
t565 = -qJD(1) * t130 + qJD(2) * t44 + qJD(3) * t84;
t538 = (t503 * t173 + t507 * t174 + (t144 * t786 + t151 * t528) * pkin(4)) * t847 - t760 / 0.2e1 + t759 / 0.2e1 - t758 / 0.2e1 - t757 / 0.2e1 + t671 / 0.2e1 + t507 * t404 + t380 * t669 + t381 * t668 + t615 / 0.2e1 + t598 / 0.2e1 + t873 * t607;
t543 = (-pkin(5) * t195 + qJ(6) * t196) * t848 + t750 / 0.2e1 + t749 / 0.2e1 + t748 / 0.2e1 - t747 / 0.2e1 + pkin(5) * t404 + t719 * t842;
t19 = t538 + t543;
t551 = (t186 * t620 + t555 * t672 + t880) * t847;
t22 = -t584 / 0.2e1 + t551 + t888 * t629;
t556 = -t672 * t887 - t888 * t782;
t255 = -(t503 * t786 + t507 * t528) * t845 - t556;
t550 = (t372 * t620 + (t507 + t672) * t858) * t847;
t29 = -t581 / 0.2e1 + (-(t793 + t839 + t669) * t578 - (-t507 / 0.2e1 - pkin(5) / 0.2e1 - t672 / 0.2e1) * t577) * mrSges(7,2) + t550 + t888 * (t819 - t858 / 0.2e1);
t560 = t22 * qJD(1) + t19 * qJD(2) + t29 * qJD(3) - t255 * qJD(4);
t476 = m(7) * t503 + mrSges(7,3);
t553 = -t515 + ((-qJ(6) - t503) * t533 + t174) * t847;
t53 = -t784 / 0.2e1 + t553;
t558 = -qJD(2) * t53 - qJD(4) * t476 + t915;
t55 = -t515 + (t650 / 0.4e1 - t690 / 0.2e1 + t699 / 0.4e1 - t174 / 0.4e1) * t852;
t557 = qJD(2) * t55 + qJD(4) * t514 + t901 + t915;
t526 = t533 ^ 2;
t524 = t530 ^ 2;
t482 = t524 * pkin(8) * t700;
t447 = mrSges(7,3) + (qJ(6) + 0.2e1 * t668) * m(7);
t197 = t682 + t783;
t129 = m(7) * t633 - t265 * t847;
t103 = t783 / 0.2e1 + m(7) * t819 + t682;
t66 = t785 / 0.2e1 + m(7) * t892;
t58 = (-t689 + t707 + t236) * t847;
t54 = (t174 - 0.2e1 * t690) * t847 - t515 + t894 + m(7) * t835;
t52 = t894 + t784 / 0.2e1 + t553;
t45 = t533 * t653 + t731 / 0.2e1 - t552 + t587;
t26 = t736 / 0.2e1 - t735 / 0.2e1 - t738 / 0.2e1 - t737 / 0.2e1 + t581 / 0.2e1 + pkin(5) * t653 + qJ(6) * t652 - t664 * t858 + t662 * t372 + (t702 / 0.2e1 - t703 / 0.2e1 + (t528 * t801 + t786 * t803) * pkin(4)) * mrSges(7,2) + t550 + t603;
t21 = t584 / 0.2e1 - t755 / 0.2e1 - t754 / 0.2e1 - t751 / 0.2e1 + t752 / 0.2e1 - t664 * t555 + t662 * t186 + t551;
t18 = t538 - t543 + t604;
t16 = t582 / 0.2e1 + t544 + t605;
t14 = (t796 + t625) * t448 + t546 + t886 * t578 * t629 + t903;
t11 = t583 / 0.2e1 + t585 * t847 + (t687 * t847 + t564) * t186 + t549 + t606;
t8 = (t806 + t626) * t448 + t902 + t911;
t5 = t537 + t535;
t3 = t536 - t554 + t853;
t1 = (t372 * t663 + t574) * t433 + (-t663 * t858 + t573) * t432 + (t722 / 0.4e1 - t427 / 0.4e1 + t860) * t529 + (t429 / 0.4e1 + t865) * t532 + (t667 + (-t488 / 0.4e1 - t870 / 0.4e1) * t529 + (-t525 / 0.2e1 - t523 / 0.2e1) * pkin(9) * mrSges(5,3) + (t794 + t795 + (m(6) * t792 + t822) * pkin(4)) * t532) * t530 - ((t831 + t836) * mrSges(6,3) + (t837 + t831) * mrSges(7,2) + t575) * t577 - ((t835 + t832) * mrSges(6,3) + (t838 + t832) * mrSges(7,2) + t576) * t578 + (pkin(4) * t704 + t563) * t849 + t516 * t787 + t899 + t923;
t27 = [qJD(2) * t25 + qJD(3) * t24, t5 * qJD(3) + t8 * qJD(4) + t11 * qJD(5) + t58 * qJD(6) + t717 + (-t236 * t380 + t236 * t381 + t237 * t377 + t237 * t378 + t401 * t474 + t402 * t472 + ((-mrSges(4,1) * t533 - mrSges(3,1) + t779) * t531 + (-mrSges(3,2) + (t524 + t526) * mrSges(4,3) + t861 * t530) * t534) * t527 + 0.2e1 * (t144 * t237 + t151 * t236 + t234 * t659) * t847 + 0.2e1 * (-t173 * t236 + t174 * t237 + t477 * t659) * t849 + 0.2e1 * (t401 * t410 + t402 * t411 + t482) * t850 + m(4) * (t482 + (pkin(8) * t526 * t534 - pkin(2) * t531) * t527)) * qJD(2), t718 + t5 * qJD(2) + t14 * qJD(4) + t16 * qJD(5) + t129 * qJD(6) + ((t306 * t449 + t613) * t847 + (t449 * t508 + t613) * t849 + (-pkin(3) * t449 - t448 * t919) * t850) * t851 + ((-mrSges(5,3) * t869 + mrSges(4,2)) * t448 + (t721 + t875) * t449 + t886 * (-t265 * t578 - t266 * t577)) * qJD(3), t8 * qJD(2) + t14 * qJD(3) + (m(7) * (-t186 * t503 + t880) - t755 - t754 - t751 + t752 + (-t186 * t528 - t555 * t786) * t846 - t580 * mrSges(5,2) - t349 * mrSges(5,1)) * qJD(4) + t21 * qJD(5) + t66 * qJD(6), t11 * qJD(2) + t16 * qJD(3) + t21 * qJD(4) + (t186 * t887 - t555 * t888) * qJD(5) + (t897 * qJD(5) / 0.2e1 + qJD(6) * t892) * t852, t58 * qJD(2) + t129 * qJD(3) + t66 * qJD(4) + qJD(5) * t785; qJD(3) * t6 + qJD(4) * t9 + qJD(5) * t12 - qJD(6) * t57 - t717, qJD(3) * t7 + qJD(4) * t10 + qJD(5) * t13 + qJD(6) * t48, t1 * qJD(4) + t3 * qJD(5) + t45 * qJD(6) + ((-pkin(3) * t522 + pkin(9) * t864) * t850 + (-t184 * t372 + t187 * t858 + t478 * t508) * t849 + (t159 * t858 + t160 * t372 + t235 * t306) * t847) * t851 + t593 + (t160 * t728 - t577 * t830 + t270 * t803 + t428 * t789 + pkin(8) * t779 + t473 * t521 - Ifges(4,6) * t530 + t508 * t294 + t478 * t340 - pkin(3) * t455 + t235 * t339 + t306 * t293 + (pkin(8) * t721 + t486 * t791 + t488 * t789 + Ifges(4,5)) * t533 - t184 * t727 + t187 * t729 + t159 * t730 + t874 * t809 + (t274 + t272) * t801 + (t430 / 0.2e1 - pkin(9) * t475) * t529 + t859 * t434 + t872 * t858 + (-t382 + t383) * t372 + t864 * mrSges(5,3) + (Ifges(5,5) * t529 + Ifges(5,6) * t532 + t577 * t917 + t578 * t918) * t790) * qJD(3), t1 * qJD(3) + (-Ifges(5,5) * t695 - Ifges(5,6) * t693 + m(7) * (t195 * t507 + t196 * t503) + t747 + t671 - t748 - t750 - t749 + (-t195 * t786 + t196 * t528) * t846 - t507 * t733 - t410 * mrSges(5,2) - t411 * mrSges(5,1) + t598 + t615 + t604) * qJD(4) + t18 * qJD(5) + t52 * qJD(6) + t592, t3 * qJD(3) + t18 * qJD(4) + (m(7) * (-pkin(5) * t174 + qJ(6) * t173) + t759 - t760 - t758 - t757 + t595 * mrSges(7,2) + t604) * qJD(5) + t54 * qJD(6) + t591, qJD(3) * t45 + qJD(4) * t52 + qJD(5) * t54 - t590; -qJD(2) * t6 + qJD(4) * t15 + qJD(5) * t17 + qJD(6) * t130 - t718, qJD(4) * t2 + qJD(5) * t4 - qJD(6) * t44 - t593, qJD(4) * t20 + qJD(5) * t23 - qJD(6) * t84 (m(7) * (-t372 * t503 + t507 * t858) + t736 - t735 + (-t372 * t528 - t786 * t858) * t846 - t738 - t737 - t577 * t610 - t578 * t677 + t603 + t619 + t862 * pkin(9) + (t702 - t703) * mrSges(7,2)) * qJD(4) + t26 * qJD(5) + t103 * qJD(6) + t568, t26 * qJD(4) + (t594 * mrSges(7,2) + t603 + (-m(7) * pkin(5) - t888) * t858 + (mrSges(6,2) - t514) * t372) * qJD(5) + t197 * qJD(6) + t567, qJD(4) * t103 + qJD(5) * t197 - t565; -qJD(2) * t9 - qJD(3) * t15 + qJD(5) * t22, -qJD(3) * t2 + qJD(5) * t19 + qJD(6) * t53 - t592, qJD(5) * t29 - t568 - t914, -qJD(5) * t255 + qJD(6) * t476 ((-pkin(5) * t528 + qJ(6) * t786) * t845 + t556) * qJD(5) + t447 * qJD(6) + t560, qJD(5) * t447 - t558; -qJD(2) * t12 - qJD(3) * t17 - qJD(4) * t22, -qJD(3) * t4 - qJD(4) * t19 + qJD(6) * t55 - t591, -qJD(4) * t29 - t567 + t914, -t560 + t900, t900, t557; qJD(2) * t57 - qJD(3) * t130, qJD(3) * t44 - qJD(4) * t53 - qJD(5) * t55 + t590, t565 + (qJD(4) - qJD(5)) * t325, t558 - t901, -t557, 0;];
Cq  = t27;
