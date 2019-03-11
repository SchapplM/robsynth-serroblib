% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRRPP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:11:15
% EndTime: 2019-03-09 21:11:54
% DurationCPUTime: 20.58s
% Computational Cost: add. (22937->945), mult. (48698->1164), div. (0->0), fcn. (47361->6), ass. (0->455)
t895 = Ifges(5,5) + Ifges(7,5);
t832 = Ifges(6,4) - t895;
t896 = Ifges(7,4) + Ifges(6,5);
t831 = Ifges(5,6) - t896;
t477 = sin(qJ(4));
t478 = sin(qJ(3));
t480 = cos(qJ(4));
t481 = cos(qJ(3));
t418 = t477 * t481 + t478 * t480;
t479 = sin(qJ(2));
t386 = t418 * t479;
t898 = -t386 / 0.2e1;
t897 = Ifges(5,1) + Ifges(7,3);
t871 = Ifges(7,2) + Ifges(6,3);
t894 = -Ifges(6,6) + Ifges(7,6);
t749 = pkin(3) * t481;
t458 = -pkin(2) - t749;
t672 = qJ(5) * t418;
t534 = t458 - t672;
t648 = t480 * t481;
t417 = t477 * t478 - t648;
t476 = pkin(4) + qJ(6);
t655 = t476 * t417;
t208 = t534 + t655;
t293 = -mrSges(7,2) * t418 + mrSges(7,3) * t417;
t586 = m(7) * t208 + t293;
t802 = -pkin(9) - pkin(8);
t436 = t802 * t478;
t437 = t802 * t481;
t318 = -t480 * t436 - t437 * t477;
t257 = t418 * pkin(5) + t318;
t527 = t417 * t832 - t831 * t418;
t842 = t436 * t477 - t480 * t437;
t856 = -pkin(5) * t417 + t842;
t893 = -t257 * mrSges(7,2) - t856 * mrSges(7,3) + t527;
t687 = t417 * mrSges(7,1);
t611 = t687 / 0.2e1;
t612 = -t687 / 0.2e1;
t279 = t612 + t611;
t892 = qJD(2) * t279;
t891 = qJD(5) * t279;
t652 = t478 * t479;
t388 = -t477 * t652 + t479 * t648;
t482 = cos(qJ(2));
t328 = t388 * mrSges(6,1) - mrSges(6,2) * t482;
t691 = t388 * mrSges(5,3);
t331 = -mrSges(5,1) * t482 - t691;
t589 = t331 / 0.2e1 - t328 / 0.2e1;
t847 = mrSges(5,3) + mrSges(6,1);
t890 = -t847 * (t318 * t898 - t842 * t388 / 0.2e1) + t589 * t842;
t469 = Ifges(4,4) * t481;
t433 = Ifges(4,1) * t478 + t469;
t760 = t481 / 0.2e1;
t888 = t433 * t760;
t747 = pkin(4) * t417;
t264 = t534 + t747;
t296 = -mrSges(6,2) * t417 - mrSges(6,3) * t418;
t843 = m(6) * t264 + t296;
t737 = mrSges(7,2) + mrSges(6,3);
t274 = qJ(5) * (-m(7) - m(6)) - t737;
t885 = qJD(4) * t274;
t441 = m(7) * t476 + mrSges(7,3);
t884 = qJD(4) * t441;
t540 = Ifges(5,2) - Ifges(6,2) + t871 - t897;
t633 = Ifges(5,4) - t894;
t882 = (t458 * mrSges(5,1) - t264 * mrSges(6,2) + t208 * mrSges(7,3) - t418 * t633) * t418 + (-t458 * mrSges(5,2) + t208 * mrSges(7,2) + t264 * mrSges(6,3) + t417 * t633 + t418 * t540) * t417;
t750 = pkin(3) * t480;
t457 = -pkin(4) - t750;
t447 = qJ(6) - t457;
t751 = pkin(3) * t477;
t455 = qJ(5) + t751;
t881 = -t257 * t455 - t856 * t447;
t679 = t482 * mrSges(7,2);
t695 = t386 * mrSges(7,1);
t327 = -t679 - t695;
t786 = -t327 / 0.2e1;
t467 = t482 * mrSges(7,3);
t692 = t388 * mrSges(7,1);
t325 = t467 + t692;
t788 = -t325 / 0.2e1;
t880 = -t257 * t786 + t856 * t788;
t879 = -qJD(5) * t274 + qJD(6) * t441;
t817 = m(7) / 0.2e1;
t820 = m(6) / 0.2e1;
t498 = (-pkin(4) * t842 - qJ(5) * t318) * t820 + (-qJ(5) * t257 - t476 * t856) * t817;
t740 = t479 * pkin(8);
t430 = -pkin(2) * t482 - pkin(1) - t740;
t415 = t481 * t430;
t650 = t479 * t481;
t570 = -pkin(9) * t650 + t415;
t506 = (-pkin(7) * t478 - pkin(3)) * t482 + t570;
t646 = t481 * t482;
t631 = pkin(7) * t646;
t309 = t631 + (-pkin(9) * t479 + t430) * t478;
t649 = t480 * t309;
t136 = t477 * t506 + t649;
t671 = qJ(5) * t482;
t120 = -t136 + t671;
t743 = t386 * pkin(5);
t79 = -t120 - t743;
t869 = t257 * t79;
t694 = t386 * mrSges(5,3);
t329 = mrSges(5,2) * t482 - t694;
t468 = t482 * mrSges(6,3);
t696 = t386 * mrSges(6,1);
t326 = t468 + t696;
t787 = t326 / 0.2e1;
t590 = -t329 / 0.2e1 + t787;
t864 = t590 * t318;
t387 = t482 * t418;
t651 = t478 * t482;
t389 = -t477 * t651 + t480 * t646;
t863 = t387 * t871 + t389 * t894 + t479 * t896;
t862 = t895 * t479 + t897 * t389 + (-Ifges(5,4) + Ifges(7,6)) * t387;
t717 = Ifges(7,6) * t418;
t721 = Ifges(6,6) * t418;
t861 = t871 * t417 + t717 - t721;
t719 = Ifges(7,6) * t388;
t723 = Ifges(6,6) * t388;
t860 = t871 * t386 + t719 - t723;
t718 = Ifges(7,6) * t417;
t729 = Ifges(5,4) * t417;
t859 = t418 * t897 + t718 - t729;
t720 = Ifges(7,6) * t386;
t731 = Ifges(5,4) * t386;
t858 = t388 * t897 + t720 - t731;
t654 = t477 * t309;
t135 = -t480 * t506 + t654;
t474 = t482 * pkin(4);
t123 = t135 + t474;
t472 = t479 * pkin(7);
t428 = pkin(3) * t652 + t472;
t673 = qJ(5) * t388;
t567 = t428 - t673;
t656 = t476 * t386;
t124 = t567 + t656;
t748 = pkin(4) * t386;
t191 = t567 + t748;
t741 = t479 * pkin(2);
t438 = -pkin(8) * t482 + t741;
t361 = pkin(7) * t652 + t481 * t438;
t276 = t479 * pkin(3) - pkin(9) * t646 + t361;
t362 = -pkin(7) * t650 + t478 * t438;
t313 = -pkin(9) * t651 + t362;
t147 = t477 * t276 + t480 * t313;
t130 = -qJ(5) * t479 - t147;
t146 = t276 * t480 - t477 * t313;
t131 = -pkin(4) * t479 - t146;
t571 = Ifges(5,6) / 0.2e1 - Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1;
t851 = -Ifges(6,4) / 0.2e1;
t572 = Ifges(7,5) / 0.2e1 + t851 + Ifges(5,5) / 0.2e1;
t73 = pkin(5) * t389 - t476 * t479 - t146;
t807 = -t73 / 0.2e1;
t810 = mrSges(6,2) / 0.2e1;
t811 = -mrSges(5,2) / 0.2e1;
t85 = -pkin(5) * t387 - t130;
t491 = -t571 * t387 + t572 * t389 - t130 * mrSges(6,3) / 0.2e1 + t131 * t810 + t146 * mrSges(5,1) / 0.2e1 + t147 * t811 + mrSges(7,3) * t807 + t85 * mrSges(7,2) / 0.2e1;
t722 = Ifges(6,6) * t417;
t549 = -Ifges(6,2) * t418 + t722;
t724 = Ifges(6,6) * t386;
t550 = -Ifges(6,2) * t388 + t724;
t728 = Ifges(5,4) * t418;
t553 = -Ifges(5,2) * t417 + t728;
t730 = Ifges(5,4) * t388;
t554 = -Ifges(5,2) * t386 + t730;
t688 = t417 * mrSges(6,1);
t613 = -t688 / 0.2e1;
t615 = -t695 / 0.2e1;
t667 = t856 * t388;
t669 = t136 * t418;
t685 = t418 * mrSges(7,1);
t686 = t418 * mrSges(6,1);
t742 = t388 * pkin(5);
t93 = -t135 - t742;
t517 = -t474 + t93;
t71 = qJ(6) * t482 - t517;
t809 = -mrSges(5,3) / 0.2e1;
t827 = -t458 * (mrSges(5,1) * t388 - mrSges(5,2) * t386) / 0.2e1 - t264 * (-mrSges(6,2) * t388 + mrSges(6,3) * t386) / 0.2e1 - t208 * (mrSges(7,2) * t386 + mrSges(7,3) * t388) / 0.2e1 + t527 * t482 / 0.4e1;
t855 = -t827 - t79 * t685 / 0.2e1 + t120 * t686 / 0.2e1 - (Ifges(6,2) * t386 - Ifges(5,6) * t482 + t554 + t723) * t418 / 0.4e1 - (Ifges(6,2) * t417 + t553 + t721) * t388 / 0.4e1 - mrSges(7,1) * t667 / 0.2e1 - t890 + t669 * t809 + t257 * t615 + t71 * t612 + t123 * t613 + t428 * (mrSges(5,1) * t418 - mrSges(5,2) * t417) / 0.2e1 + t191 * (-mrSges(6,2) * t418 + mrSges(6,3) * t417) / 0.2e1 + t124 * (mrSges(7,2) * t417 + mrSges(7,3) * t418) / 0.2e1 + (-t417 * t897 + t717 - t728 + t861) * t388 / 0.4e1 + (-t386 * t897 - t482 * t896 + t719 - t730 + t860) * t418 / 0.4e1 + (t871 * t418 + t549 - t718 + t722) * t386 / 0.4e1 + (-Ifges(6,4) * t482 + t871 * t388 + t550 - t720 + t724) * t417 / 0.4e1 - (-Ifges(5,2) * t388 - t482 * t895 - t731 + t858) * t417 / 0.4e1 - (-Ifges(5,2) * t418 - t729 + t859) * t386 / 0.4e1 + t491;
t854 = -t318 * t455 + t842 * t457;
t853 = t318 * t477 + t480 * t842;
t779 = -t387 / 0.2e1;
t773 = t389 / 0.2e1;
t850 = mrSges(7,1) + mrSges(6,1);
t849 = -mrSges(5,2) + mrSges(6,3);
t848 = -mrSges(6,2) + mrSges(5,1);
t725 = Ifges(4,6) * t478;
t727 = Ifges(4,5) * t481;
t529 = t727 / 0.2e1 - t725 / 0.2e1;
t846 = Ifges(3,4) - t529;
t845 = t633 * t386;
t844 = t633 * t388;
t841 = -Ifges(4,2) * t478 + t469;
t840 = (-mrSges(7,3) - t848) * t477 + (-mrSges(5,2) + t737) * t480;
t659 = t455 * t388;
t839 = -t659 / 0.2e1 - t673 / 0.2e1;
t837 = t615 - t696;
t836 = -t361 * t478 + t362 * t481;
t835 = -mrSges(4,1) * t481 + mrSges(4,2) * t478;
t816 = m(7) / 0.4e1;
t833 = t816 + m(6) / 0.4e1;
t528 = t386 * t832 - t388 * t831;
t625 = -t750 / 0.2e1;
t626 = t751 / 0.2e1;
t658 = t455 * t418;
t830 = t417 * t625 + t418 * t626 - t658 / 0.2e1 - t672 / 0.2e1;
t680 = t479 * mrSges(7,3);
t689 = t389 * mrSges(7,1);
t321 = -t680 + t689;
t682 = t479 * mrSges(6,2);
t690 = t389 * mrSges(6,1);
t324 = t682 + t690;
t573 = Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1 + Ifges(5,3) / 0.2e1;
t322 = mrSges(6,1) * t387 - mrSges(6,3) * t479;
t681 = t479 * mrSges(7,2);
t693 = t387 * mrSges(7,1);
t323 = t681 - t693;
t591 = -t322 / 0.2e1 + t323 / 0.2e1;
t829 = (qJ(5) * t85 - t476 * t73) * t817 + (-pkin(4) * t131 - qJ(5) * t130) * t820 + t573 * t479 + t591 * qJ(5) - t476 * t321 / 0.2e1 - pkin(4) * t324 / 0.2e1;
t764 = t457 / 0.2e1;
t765 = -t447 / 0.2e1;
t812 = -mrSges(4,2) / 0.2e1;
t814 = mrSges(4,1) / 0.2e1;
t828 = (-t130 * t455 + t131 * t457) * t820 + (-t447 * t73 + t455 * t85) * t817 + t591 * t455 + t324 * t764 + t321 * t765 + t362 * t812 + t361 * t814;
t825 = 2 * qJD(3);
t824 = 0.2e1 * qJD(4);
t823 = -m(5) / 0.2e1;
t822 = m(5) / 0.2e1;
t821 = -m(6) / 0.2e1;
t818 = -m(7) / 0.2e1;
t815 = pkin(8) / 0.2e1;
t813 = mrSges(6,1) / 0.2e1;
t808 = t71 / 0.2e1;
t806 = t79 / 0.2e1;
t805 = t93 / 0.2e1;
t94 = t136 - t743;
t804 = t94 / 0.2e1;
t632 = pkin(7) * t651;
t308 = t570 - t632;
t155 = t308 * t477 + t649;
t106 = t155 - t743;
t801 = -t106 / 0.2e1;
t156 = t308 * t480 - t654;
t107 = t156 - t742;
t800 = t107 / 0.2e1;
t799 = -t120 / 0.2e1;
t798 = t123 / 0.2e1;
t797 = -t136 / 0.2e1;
t796 = -t155 / 0.2e1;
t795 = t156 / 0.2e1;
t794 = Ifges(6,2) * t773 + Ifges(6,6) * t779 + t479 * t851;
t249 = -mrSges(7,2) * t388 + mrSges(7,3) * t386;
t793 = -t249 / 0.2e1;
t253 = -mrSges(6,2) * t386 - mrSges(6,3) * t388;
t792 = -t253 / 0.2e1;
t791 = -t293 / 0.2e1;
t790 = -t296 / 0.2e1;
t789 = -t318 / 0.2e1;
t781 = t386 / 0.2e1;
t778 = t387 / 0.2e1;
t774 = -t389 / 0.2e1;
t772 = -t417 / 0.2e1;
t770 = t417 / 0.2e1;
t767 = t418 / 0.2e1;
t763 = -t478 / 0.2e1;
t762 = t478 / 0.2e1;
t758 = t482 / 0.2e1;
t756 = m(6) * t842;
t755 = m(7) * t107;
t754 = m(7) * t257;
t753 = m(7) * t856;
t745 = pkin(8) * t478;
t744 = pkin(8) * t481;
t471 = t478 * pkin(3);
t473 = t482 * pkin(7);
t739 = t93 * mrSges(7,2);
t738 = t94 * mrSges(7,3);
t475 = m(7) * qJD(5);
t736 = m(7) * qJD(6);
t732 = Ifges(4,4) * t478;
t716 = pkin(3) * qJD(3);
t714 = t106 * mrSges(7,3);
t713 = t107 * mrSges(7,2);
t712 = t124 * mrSges(7,2);
t711 = t124 * mrSges(7,3);
t710 = t135 * mrSges(5,2);
t709 = t135 * mrSges(6,3);
t708 = t136 * mrSges(5,1);
t707 = t136 * mrSges(6,2);
t706 = t191 * mrSges(6,2);
t705 = t191 * mrSges(6,3);
t700 = t318 * mrSges(5,2);
t699 = t318 * mrSges(6,3);
t698 = t842 * mrSges(5,1);
t697 = t842 * mrSges(6,2);
t684 = t428 * mrSges(5,1);
t683 = t428 * mrSges(5,2);
t678 = t482 * Ifges(4,5);
t677 = t482 * Ifges(4,6);
t429 = pkin(3) * t651 + t473;
t566 = -qJ(5) * t389 + t429;
t125 = t387 * t476 + t566;
t192 = pkin(4) * t387 + t566;
t229 = Ifges(5,4) * t389 - Ifges(5,2) * t387 + Ifges(5,6) * t479;
t251 = -mrSges(7,2) * t389 + mrSges(7,3) * t387;
t252 = mrSges(5,1) * t386 + mrSges(5,2) * t388;
t254 = mrSges(5,1) * t387 + mrSges(5,2) * t389;
t255 = -mrSges(6,2) * t387 - mrSges(6,3) * t389;
t330 = -mrSges(5,2) * t479 - mrSges(5,3) * t387;
t332 = mrSges(5,1) * t479 - mrSges(5,3) * t389;
t359 = t415 - t632;
t360 = t478 * t430 + t631;
t379 = Ifges(4,6) * t479 + t482 * t841;
t559 = Ifges(4,1) * t481 - t732;
t381 = Ifges(4,5) * t479 + t482 * t559;
t560 = mrSges(4,1) * t478 + mrSges(4,2) * t481;
t405 = t560 * t482;
t424 = mrSges(4,2) * t482 - mrSges(4,3) * t652;
t425 = -t479 * mrSges(4,2) - mrSges(4,3) * t651;
t426 = -mrSges(4,1) * t482 - mrSges(4,3) * t650;
t427 = t479 * mrSges(4,1) - mrSges(4,3) * t646;
t519 = t559 * t479;
t380 = t519 - t678;
t647 = t481 * t380;
t378 = t479 * t841 - t677;
t653 = t478 * t378;
t5 = (t647 / 0.2e1 - t653 / 0.2e1 - pkin(1) * mrSges(3,2) + t832 * t389 + t831 * t387 + (-Ifges(3,2) + Ifges(3,1) - Ifges(6,1) - Ifges(7,1) - Ifges(4,3) - Ifges(5,3) + (m(4) * pkin(7) + t560) * pkin(7)) * t479 + t846 * t482) * t482 + (-pkin(1) * mrSges(3,1) + pkin(7) * t405 + t379 * t763 + t381 * t760 - t571 * t386 - t846 * t479) * t479 + t229 * t898 + (t572 * t479 + t794 + t862 / 0.2e1) * t388 + t554 * t779 + t550 * t774 + m(4) * (t359 * t361 + t360 * t362) + t362 * t424 + t360 * t425 + t361 * t426 + t359 * t427 + t428 * t254 + t429 * t252 + t131 * t328 + t147 * t329 + t136 * t330 + t146 * t331 - t135 * t332 + t79 * t323 + t123 * t324 + t73 * t325 + t130 * t326 + t85 * t327 + t71 * t321 + t120 * t322 + t192 * t253 + t191 * t255 + t125 * t249 + t124 * t251 + t858 * t773 + t860 * t778 + t863 * t781 + m(5) * (-t135 * t146 + t136 * t147 + t428 * t429) + m(6) * (t120 * t130 + t123 * t131 + t191 * t192) + m(7) * (t124 * t125 + t71 * t73 + t79 * t85);
t676 = t5 * qJD(1);
t250 = pkin(4) * t388 + qJ(5) * t386;
t451 = pkin(3) * t650;
t211 = t250 + t451;
t363 = t388 * qJ(6);
t137 = t211 + t363;
t432 = Ifges(4,2) * t481 + t732;
t485 = (-t71 * mrSges(7,1) + t482 * t572 - t683 + t705 + t712 + t845) * t386 + (t120 * mrSges(6,1) - t136 * mrSges(5,3) + t386 * t540 + t482 * t571 + t684 - t706 + t711 - t844) * t388 - t79 * t692 - t123 * t696 - t135 * t694 - t528 * t482 / 0.2e1;
t551 = Ifges(4,5) * t478 + Ifges(4,6) * t481;
t640 = -t331 + t328;
t641 = t329 - t326;
t6 = (t551 * t758 + t380 * t763 - t481 * t378 / 0.2e1 + (-pkin(7) * t835 + t432 * t762 - t888) * t479 + (m(5) * t428 + t252) * t749 + (t359 * t478 - t360 * t481) * mrSges(4,3)) * t479 + t485 + t641 * t156 + t640 * t155 + m(5) * (t135 * t155 + t136 * t156) + m(6) * (-t120 * t156 + t123 * t155 + t191 * t211) + m(7) * (t106 * t71 + t107 * t79 + t124 * t137) + t359 * t424 - t360 * t426 + t106 * t325 + t107 * t327 + t211 * t253 + t137 * t249;
t675 = t6 * qJD(1);
t166 = t250 + t363;
t7 = m(6) * (t120 * t135 + t123 * t136 + t191 * t250) + m(7) * (t124 * t166 + t71 * t94 + t79 * t93) + t485 + t640 * t136 - t641 * t135 + t94 * t325 + t93 * t327 + t250 * t253 + t166 * t249;
t674 = t7 * qJD(1);
t670 = t135 * t417;
t20 = (t326 - t327) * t482 + (-t249 - t253) * t388 + m(7) * (-t124 * t388 - t482 * t79) + m(6) * (t120 * t482 - t191 * t388);
t668 = t20 * qJD(1);
t31 = m(7) * (t124 * t386 + t482 * t71) + t482 * t325 + t386 * t249;
t666 = t31 * qJD(1);
t661 = t428 * t478;
t657 = t455 * t480;
t645 = t482 * t417;
t630 = m(7) * t805;
t629 = mrSges(4,3) * t815;
t628 = t457 * t696;
t627 = -t751 / 0.2e1;
t624 = t750 / 0.2e1;
t622 = t813 + mrSges(5,3) / 0.2e1;
t621 = t810 - mrSges(5,1) / 0.2e1;
t620 = mrSges(6,3) / 0.2e1 + t811;
t619 = Ifges(4,1) / 0.4e1 - Ifges(4,2) / 0.4e1;
t618 = t764 + pkin(4) / 0.2e1;
t614 = t695 / 0.2e1;
t610 = t685 / 0.2e1;
t596 = t418 * t758;
t595 = t155 / 0.2e1 + t797;
t594 = t795 + t135 / 0.2e1;
t588 = t765 + t476 / 0.2e1;
t587 = t455 / 0.2e1 - qJ(5) / 0.2e1;
t219 = t257 * t818;
t584 = t219 - t685;
t583 = t468 + t679;
t578 = -mrSges(4,3) * t740 / 0.2e1;
t569 = 0.2e1 * t612 - t688;
t568 = t136 - 0.2e1 * t671;
t563 = t587 * t388;
t552 = -t725 + t727;
t294 = pkin(4) * t418 + qJ(5) * t417;
t295 = mrSges(5,1) * t417 + mrSges(5,2) * t418;
t483 = (t667 / 0.2e1 + t257 * t781) * mrSges(7,1) + t491 + t827;
t270 = t294 + t471;
t406 = t418 * qJ(6);
t209 = t270 + t406;
t537 = t155 * t318 + t156 * t842;
t486 = (t120 * t318 + t123 * t842 + t191 * t270 + t211 * t264 + t537) * t821 + (t106 * t257 + t124 * t209 + t137 * t208 - t869 + (t107 + t71) * t856) * t818 + t137 * t791 + t209 * t793 + t211 * t790 + t270 * t792 + t880;
t507 = Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(7,3) / 0.2e1 - Ifges(5,2) / 0.2e1 - Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1;
t488 = t844 + (-Ifges(5,6) / 0.4e1 + Ifges(7,4) / 0.4e1 + Ifges(6,5) / 0.4e1) * t482 + t507 * t386 - t711 / 0.2e1 + t706 / 0.2e1 - t684 / 0.2e1;
t489 = -t845 + (-Ifges(7,5) / 0.4e1 + Ifges(6,4) / 0.4e1 - Ifges(5,5) / 0.4e1) * t482 + t507 * t388 - t712 / 0.2e1 - t705 / 0.2e1 + t683 / 0.2e1;
t505 = t135 * t842 - t136 * t318 + t537;
t512 = (t146 * t480 + t147 * t477) * t822;
t524 = t477 * t330 / 0.2e1 + t480 * t332 / 0.2e1;
t531 = Ifges(4,3) / 0.2e1 + t573;
t1 = t828 + ((pkin(2) * t812 + t433 / 0.4e1 + t469 / 0.4e1 - pkin(7) * mrSges(4,1) / 0.2e1 + (t629 + t619) * t478) * t478 + (pkin(2) * t814 + 0.3e1 / 0.4e1 * t732 + pkin(7) * t812 + t432 / 0.4e1 + (t629 - t619) * t481 + (-t295 / 0.2e1 + t458 * t823) * pkin(3)) * t481 + t531) * t479 + (-t380 / 0.4e1 + 0.3e1 / 0.4e1 * t678 + t426 * t815) * t481 + (-0.3e1 / 0.4e1 * t677 + t378 / 0.4e1 + t424 * t815 - pkin(3) * t252 / 0.2e1) * t478 + (-t595 * mrSges(5,3) + (t806 + t801) * mrSges(7,1) + (t799 + t796) * mrSges(6,1) + t488) * t418 + (t594 * mrSges(5,3) + (t808 + t800) * mrSges(7,1) + (t798 + t795) * mrSges(6,1) + t489) * t417 + t486 + t483 - t864 + pkin(3) * t512 + t524 * pkin(3) + (pkin(3) * t661 + t505) * t823 + t890;
t523 = t432 * t763 + t888;
t10 = -pkin(2) * t560 + t559 * t762 + t841 * t760 + t523 + (m(5) * t458 + t295) * t471 + t843 * t270 + t882 + t586 * t209;
t542 = -t1 * qJD(1) + t10 * qJD(2);
t243 = t294 + t406;
t11 = t243 * t586 + t843 * t294 + t882;
t487 = (t191 * t294 + t250 * t264 + (t123 - t135) * t842 + (t120 + t136) * t318) * t821 + (t124 * t243 + t166 * t208 + t257 * t94 - t869 + (t71 + t93) * t856) * t818 + t166 * t791 + t243 * t793 + t250 * t790 + t294 * t792 + t880;
t4 = t487 + t483 + ((t808 + t805) * mrSges(7,1) + (-t135 / 0.2e1 + t798) * mrSges(6,1) + t489) * t417 + ((-t94 / 0.2e1 + t806) * mrSges(7,1) + (t799 + t797) * mrSges(6,1) + t488) * t418 + (t386 * t622 - t590) * t318 + (t388 * t622 + t589) * t842 + t829;
t541 = -t4 * qJD(1) + t11 * qJD(2);
t490 = (t793 + t792) * t418 + (t791 + t790) * t388 + (-t418 * t191 - t264 * t388 - t482 * t842) * t820 + (-t418 * t124 - t208 * t388 - t482 * t856) * t817;
t535 = m(7) * t807 + t131 * t821;
t15 = (mrSges(7,3) / 0.2e1 - mrSges(6,2) / 0.2e1) * t479 + t490 + t535 + t850 * (t645 / 0.2e1 + t774);
t39 = (-t586 - t843) * t418;
t539 = -qJD(1) * t15 - qJD(2) * t39;
t497 = (t417 * t124 + t208 * t386 + t257 * t482) * t817 + t293 * t781 + t249 * t770;
t532 = t85 * t818 - t681 / 0.2e1;
t22 = (t596 + t778) * mrSges(7,1) + t497 + t532;
t56 = t586 * t417;
t538 = qJD(1) * t22 + qJD(2) * t56;
t533 = m(6) * t796 + m(7) * t801;
t526 = (-qJ(5) - t455) * t482 + t136;
t518 = m(7) * t804 + t136 * t820 + t615;
t513 = (t627 + t587) * t418;
t492 = (pkin(3) * t853 + t854) * t821 + ((t257 * t477 + t480 * t856) * pkin(3) + t881) * t818;
t13 = (t513 + (t624 + t588) * t417) * mrSges(7,1) + (t513 + (t624 + t618) * t417) * mrSges(6,1) + t492 + t498 + t849 * (t789 + t318 / 0.2e1);
t493 = (-t135 * t455 + t136 * t457 + (-t120 * t480 + t123 * t477) * pkin(3)) * t821 + (-t447 * t94 + t455 * t93 + (t477 * t71 + t480 * t79) * pkin(3)) * t818;
t499 = (-pkin(4) * t155 + qJ(5) * t156) * t820 + (qJ(5) * t107 - t106 * t476) * t817;
t9 = (t801 + t804) * mrSges(7,3) + (t800 - t93 / 0.2e1) * mrSges(7,2) + (t386 * t588 + t563) * mrSges(7,1) + (t386 * t618 + t563) * mrSges(6,1) + ((-t694 / 0.2e1 + t786 + t590) * t480 + (t691 / 0.2e1 + t788 + t589) * t477) * pkin(3) + t493 + t499 - t848 * t595 + t849 * t594;
t99 = (m(7) * (-t447 * t477 + t657) + m(6) * (t457 * t477 + t657) + t840) * pkin(3);
t510 = -t9 * qJD(1) - t13 * qJD(2) + t99 * qJD(3);
t494 = t615 + t526 * t820 + (t526 - t743) * t817 - t583;
t26 = t614 + t494 + t533;
t338 = 0.4e1 * t455 * t833 + t737;
t509 = qJD(1) * t26 + qJD(3) * t338 - t892;
t501 = -t467 + ((-qJ(6) - t447) * t482 + t517) * t817;
t35 = -t755 / 0.2e1 + t501;
t416 = m(7) * t447 + mrSges(7,3);
t87 = t219 + t754 / 0.2e1;
t508 = qJD(1) * t35 + qJD(2) * t87 + qJD(3) * t416;
t500 = t467 + ((-qJ(6) - t476) * t482 + t517) * t818;
t36 = t630 + t500;
t504 = -qJD(1) * t36 + qJD(3) * t441 + t884;
t495 = t568 * t821 + (t568 - t743) * t818 + t583;
t28 = t614 + t495 + t518;
t503 = -qJD(1) * t28 - qJD(3) * t274 - t885 + t892;
t466 = qJ(5) * t750;
t423 = (qJD(1) * t482 - qJD(3) - qJD(4)) * m(7);
t396 = (0.4e1 * pkin(4) + 0.4e1 * qJ(6) + 0.2e1 * t750) * t816 + mrSges(7,3) + m(7) * t624;
t275 = (t817 + t820) * t751 + t737 + t833 * (0.4e1 * qJ(5) + 0.2e1 * t751);
t61 = -t754 / 0.2e1 + t584;
t60 = -t257 * t817 + t584;
t52 = t569 + t753 + t756;
t41 = t753 / 0.2e1 + t756 / 0.2e1 + t842 * t820 + t856 * t817 + t569;
t33 = -t500 + t630 - t692;
t32 = -t692 + t755 / 0.2e1 + t501;
t24 = -t495 + t518 + t837;
t23 = t494 - t533 + t837;
t21 = mrSges(7,1) * t596 - t693 / 0.2e1 + t497 - t532;
t14 = -t680 / 0.2e1 + t689 / 0.2e1 + t682 / 0.2e1 + t690 / 0.2e1 + (t813 + mrSges(7,1) / 0.2e1) * t645 + t490 - t535;
t12 = -t620 * t318 + t621 * t842 + t457 * t613 + t447 * t611 + t498 - t698 / 0.2e1 + t697 / 0.2e1 - t699 / 0.2e1 + t700 / 0.2e1 - t492 + (t655 / 0.2e1 + t830) * mrSges(7,1) + (t747 / 0.2e1 + t830) * mrSges(6,1) + t893;
t8 = -t628 / 0.2e1 + t499 + t707 / 0.2e1 + t528 + t713 / 0.2e1 + t710 / 0.2e1 - t708 / 0.2e1 - t709 / 0.2e1 + t620 * t156 + t621 * t155 + t447 * t614 + t326 * t625 - t714 / 0.2e1 - t738 / 0.2e1 + t739 / 0.2e1 - t493 + (t691 + t331) * t627 + (t328 + t325) * t626 + (t656 / 0.2e1 + t839) * mrSges(7,1) + (t748 / 0.2e1 + t839) * mrSges(6,1) + (t694 + t329 + t327) * t624;
t3 = t610 * t94 + t612 * t93 + t670 * t813 - t487 + t847 * t669 / 0.2e1 + t864 + t829 + t855;
t2 = t828 + t855 - t426 * t744 / 0.2e1 - t424 * t745 / 0.2e1 + t318 * t787 + t560 * t472 / 0.2e1 + t252 * t471 / 0.2e1 + (-t552 / 0.4e1 + t529) * t482 + (t519 / 0.4e1 + t578 * t481) * t481 - (t841 + 0.2e1 * t433) * t652 / 0.4e1 + t835 * t741 / 0.2e1 - t432 * t650 / 0.2e1 + t478 ^ 2 * t578 - t486 + t295 * t451 / 0.2e1 - t653 / 0.4e1 + t670 * t809 + t329 * t789 + t106 * t610 + t107 * t612 + t647 / 0.4e1 + t847 * (t155 * t767 + t156 * t772) + t531 * t479 + (t512 + t524) * pkin(3) + ((t458 * t650 + t661) * pkin(3) + t505) * t822;
t16 = [qJD(2) * t5 + qJD(3) * t6 + qJD(4) * t7 + qJD(5) * t20 + qJD(6) * t31, t2 * qJD(3) + t3 * qJD(4) + t14 * qJD(5) + t21 * qJD(6) + t676 + (-t427 * t745 - t85 * t687 + t842 * t330 - t842 * t322 + 0.2e1 * (-t130 * t842 + t131 * t318 + t192 * t264) * t820 + 0.2e1 * (-t146 * t318 + t147 * t842 + t429 * t458) * t822 + (Ifges(3,5) + (-mrSges(3,1) + t835) * pkin(7) + t523) * t482 + t856 * t323 + 0.2e1 * (t125 * t208 + t257 * t73 + t85 * t856) * t817 + m(4) * (-pkin(2) * t473 + pkin(8) * t836) + t836 * mrSges(4,3) + t553 * t779 + t418 * t794 + t229 * t772 + t549 * t774 + t379 * t760 + t381 * t762 + mrSges(3,2) * t472 + t425 * t744 + t73 * t685 + t131 * t686 + t130 * t688 - Ifges(3,6) * t479 + t458 * t254 + t429 * t295 - pkin(2) * t405 - t318 * t332 + t318 * t324 + t257 * t321 + t192 * t296 + t125 * t293 + t264 * t255 + t208 * t251 + (-t146 * t418 - t147 * t417) * mrSges(5,3) + t863 * t770 + t859 * t773 + t861 * t778 + t862 * t767 + (-t831 * t417 - t832 * t418 + t551) * t479 / 0.2e1) * qJD(2), t675 + t2 * qJD(2) + (-t360 * mrSges(4,1) - t359 * mrSges(4,2) - Ifges(4,5) * t652 - Ifges(4,6) * t650 - t155 * t848 + t156 * t849 + t447 * t695 - t659 * t850 + t528 - t628 + t713 - t714) * qJD(3) + t8 * qJD(4) + t23 * qJD(5) + t32 * qJD(6) + ((-t106 * t447 + t107 * t455) * t817 + (t155 * t457 + t156 * t455) * t820) * t825 + (m(5) * (-t155 * t480 + t156 * t477) + (t386 * t480 - t388 * t477) * mrSges(5,3)) * t716, t674 + t3 * qJD(2) + t8 * qJD(3) + (-mrSges(6,1) * t673 + pkin(4) * t696 + t528 + t707 - t708 - t709 + t710 - t738 + t739 + (t656 - t673) * mrSges(7,1)) * qJD(4) + t24 * qJD(5) + t33 * qJD(6) + ((qJ(5) * t93 - t476 * t94) * t817 + (-pkin(4) * t136 - qJ(5) * t135) * t820) * t824, qJD(2) * t14 + qJD(3) * t23 + qJD(4) * t24 + t668, qJD(2) * t21 + qJD(3) * t32 + qJD(4) * t33 + t666; -qJD(3) * t1 - qJD(4) * t4 + qJD(5) * t15 + qJD(6) * t22 - t676, qJD(3) * t10 + qJD(4) * t11 + qJD(5) * t39 + qJD(6) * t56 (pkin(8) * t835 - t318 * t849 + t447 * t687 - t457 * t688 - t658 * t850 - t842 * t848 + t552 + t893) * qJD(3) + t12 * qJD(4) + t41 * qJD(5) + t61 * qJD(6) + (t817 * t881 + t854 * t820) * t825 + (-m(5) * t853 + (t417 * t480 - t418 * t477) * mrSges(5,3)) * t716 + t542, t12 * qJD(3) + (-mrSges(6,1) * t672 + pkin(4) * t688 + t697 - t698 - t699 + t700 + (t655 - t672) * mrSges(7,1) + t893) * qJD(4) + t52 * qJD(5) + t60 * qJD(6) + t498 * t824 + t541, qJD(3) * t41 + qJD(4) * t52 - t539, qJD(3) * t61 + qJD(4) * t60 + t538; qJD(2) * t1 - qJD(4) * t9 + qJD(5) * t26 + qJD(6) * t35 - t675, -qJD(4) * t13 + qJD(6) * t87 - t542 - t891, qJD(4) * t99 + qJD(5) * t338 + qJD(6) * t416, t275 * qJD(5) + t396 * qJD(6) + ((-t476 * t751 + t466) * t817 + (-pkin(4) * t751 + t466) * t820) * t824 + t510 + t840 * pkin(3) * qJD(4), qJD(4) * t275 + t509, qJD(4) * t396 + t508; qJD(2) * t4 + qJD(3) * t9 - qJD(5) * t28 - qJD(6) * t36 - t674, qJD(3) * t13 - t541 + t891, -t510 + t879, t879, t503, t504; -t15 * qJD(2) - t26 * qJD(3) + t28 * qJD(4) + t482 * t736 - t668, t539 + (qJD(3) - qJD(4)) * t279, -t509 - t736 + t885, -t503 - t736, 0, t423; -t22 * qJD(2) - t35 * qJD(3) + t36 * qJD(4) - t475 * t482 - t666, -qJD(3) * t87 - t538, t475 - t508 - t884, t475 - t504, -t423, 0;];
Cq  = t16;
