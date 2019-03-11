% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRRPP3
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
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRRPP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:53:55
% EndTime: 2019-03-09 20:54:22
% DurationCPUTime: 14.06s
% Computational Cost: add. (21470->844), mult. (42947->1036), div. (0->0), fcn. (43440->6), ass. (0->443)
t751 = sin(qJ(3));
t752 = sin(qJ(2));
t753 = cos(qJ(3));
t754 = cos(qJ(2));
t445 = t751 * t752 - t753 * t754;
t446 = -t751 * t754 - t752 * t753;
t512 = cos(qJ(4));
t504 = Ifges(5,4) * t512;
t511 = sin(qJ(4));
t589 = Ifges(5,2) * t511 - t504;
t203 = -Ifges(5,6) * t446 + t445 * t589;
t730 = Ifges(6,6) * t511;
t585 = Ifges(6,2) * t512 - t730;
t214 = -Ifges(6,4) * t446 + t445 * t585;
t729 = Ifges(6,6) * t512;
t466 = -Ifges(6,2) * t511 - t729;
t736 = Ifges(5,4) * t511;
t468 = Ifges(5,2) * t512 + t736;
t684 = t445 * t512;
t626 = t684 / 0.2e1;
t627 = -t684 / 0.2e1;
t685 = t445 * t511;
t628 = t685 / 0.2e1;
t629 = -t685 / 0.2e1;
t755 = t512 / 0.2e1;
t756 = -t512 / 0.2e1;
t757 = t511 / 0.2e1;
t758 = -t511 / 0.2e1;
t764 = -t446 / 0.2e1;
t462 = -Ifges(6,3) * t512 - t730;
t499 = Ifges(7,6) * t511;
t465 = -Ifges(7,2) * t512 + t499;
t819 = t462 + t465;
t728 = Ifges(7,6) * t512;
t461 = Ifges(7,3) * t511 - t728;
t470 = Ifges(5,1) * t511 + t504;
t820 = t461 + t470;
t584 = -Ifges(6,3) * t511 + t729;
t210 = -Ifges(6,5) * t446 + t445 * t584;
t464 = Ifges(7,2) * t511 + t728;
t212 = -Ifges(7,4) * t446 - t445 * t464;
t822 = t212 + t210;
t471 = Ifges(5,1) * t512 - t736;
t205 = -Ifges(5,5) * t446 - t445 * t471;
t583 = -Ifges(7,3) * t512 - t499;
t208 = -Ifges(7,5) * t446 + t445 * t583;
t823 = t208 + t205;
t835 = Ifges(5,5) + Ifges(7,5);
t846 = Ifges(6,5) + Ifges(7,4);
t532 = -Ifges(4,5) * t445 + Ifges(4,6) * t446 + t203 * t755 + t214 * t758 + t466 * t626 + t468 * t628 + t823 * t757 + t822 * t756 + t819 * t629 + t820 * t627 + ((Ifges(5,6) - t846) * t512 + (-Ifges(6,4) + t835) * t511) * t764;
t454 = t512 * mrSges(6,2) - t511 * mrSges(6,3);
t653 = t752 * pkin(7);
t479 = -pkin(8) * t752 - t653;
t656 = t754 * pkin(7);
t481 = pkin(8) * t754 + t656;
t805 = t751 * t479 + t753 * t481;
t807 = -pkin(4) * t685 + qJ(5) * t684 + t805;
t824 = t807 * t454;
t455 = -t512 * mrSges(5,1) + t511 * mrSges(5,2);
t825 = t805 * t455;
t833 = t805 * mrSges(4,1);
t359 = -t753 * t479 + t481 * t751;
t834 = t359 * mrSges(4,2);
t457 = -t511 * mrSges(7,2) - mrSges(7,3) * t512;
t673 = t511 * qJ(6);
t821 = -t445 * t673 + t807;
t841 = t821 * t457;
t849 = t532 + t824 + t825 + t834 - t833 + t841;
t848 = t824 / 0.2e1 + t834 / 0.2e1 + t841 / 0.2e1;
t401 = mrSges(6,1) * t684;
t718 = t446 * mrSges(6,2);
t321 = -t401 - t718;
t847 = t321 / 0.2e1;
t683 = t446 * t511;
t564 = -pkin(4) * t683 + t359;
t710 = qJ(5) * t512;
t576 = -t673 + t710;
t126 = t446 * t576 + t564;
t845 = t126 * t821;
t510 = pkin(4) + qJ(6);
t674 = t511 * qJ(5);
t421 = -t510 * t512 - pkin(3) - t674;
t655 = t753 * pkin(2);
t386 = -t655 + t421;
t844 = t386 * t821;
t843 = t421 * t821;
t721 = t805 * mrSges(4,3);
t842 = t446 * t721;
t839 = t825 / 0.2e1 - t833 / 0.2e1;
t496 = -pkin(2) * t754 - pkin(1);
t304 = t445 * pkin(3) + t446 * pkin(9) + t496;
t691 = t805 * t512;
t149 = t511 * t304 + t691;
t425 = t445 * qJ(5);
t118 = -t149 - t425;
t148 = -t512 * t304 + t511 * t805;
t119 = -t445 * pkin(4) + t148;
t682 = t446 * t512;
t180 = qJ(5) * t682 + t564;
t590 = mrSges(7,2) * t512 - t511 * mrSges(7,3);
t287 = t590 * t445;
t289 = t590 * t446;
t591 = t511 * mrSges(6,2) + mrSges(6,3) * t512;
t292 = t591 * t445;
t459 = t511 * mrSges(5,1) + mrSges(5,2) * t512;
t293 = t459 * t445;
t295 = t591 * t446;
t310 = mrSges(5,2) * t446 + mrSges(5,3) * t685;
t312 = -t446 * mrSges(5,1) + mrSges(5,3) * t684;
t651 = mrSges(7,1) * t684;
t716 = t446 * mrSges(7,3);
t318 = -t651 + t716;
t319 = -mrSges(6,1) * t685 + mrSges(6,3) * t446;
t650 = mrSges(7,1) * t685;
t717 = t446 * mrSges(7,2);
t320 = t650 - t717;
t560 = t459 * t446;
t112 = pkin(5) * t682 - t148;
t810 = t445 * t510;
t70 = -t112 - t810;
t113 = t691 + (pkin(5) * t446 + t304) * t511;
t90 = t425 + t113;
t838 = -t821 * t289 - t807 * t295 + t805 * t560 - t118 * t319 - t119 * t321 - t126 * t287 + t148 * t312 - t149 * t310 - t180 * t292 + t359 * t293 - t70 * t318 - t90 * t320 - t842 - t496 * (-mrSges(4,1) * t446 - mrSges(4,2) * t445);
t836 = pkin(3) * t805;
t832 = t180 * t807;
t831 = t359 * t751;
t693 = t359 * t805;
t582 = -t512 * pkin(4) - t674;
t452 = -pkin(3) + t582;
t424 = -t655 + t452;
t830 = t424 * t807;
t829 = t452 * t807;
t495 = -t655 - pkin(3);
t828 = t495 * t805;
t827 = t511 * t359;
t826 = t512 * t359;
t508 = t511 ^ 2;
t509 = t512 ^ 2;
t818 = (t508 + t509) * t655;
t500 = Ifges(7,5) * t512;
t501 = Ifges(6,5) * t511;
t502 = Ifges(5,5) * t512;
t503 = Ifges(7,4) * t511;
t817 = Ifges(6,4) * t512 - t500 - t501 - t502 - t503;
t816 = t510 * t318 / 0.2e1 + pkin(4) * t847;
t789 = m(7) / 0.2e1;
t814 = 0.2e1 * t789;
t791 = m(6) / 0.2e1;
t813 = 0.2e1 * t791;
t812 = m(6) + m(7);
t811 = mrSges(6,1) + mrSges(5,3);
t741 = mrSges(7,2) + mrSges(6,3);
t809 = t289 + t295;
t665 = t386 + t421;
t622 = t457 * t758;
t808 = -t717 / 0.2e1 + t446 * t622;
t663 = t424 + t452;
t661 = t454 + t457;
t806 = (t622 + mrSges(7,2) / 0.2e1) * t446 - t650;
t744 = t445 * pkin(9);
t323 = -t446 * pkin(3) + t744;
t654 = t752 * pkin(2);
t305 = t654 + t323;
t162 = t305 * t512 + t827;
t163 = t511 * t305 - t826;
t573 = -t162 * t511 + t163 * t512;
t426 = t446 * qJ(5);
t121 = t426 - t163;
t743 = t446 * pkin(4);
t122 = -t162 + t743;
t575 = -t121 * t512 + t122 * t511;
t761 = t457 / 0.2e1;
t762 = t454 / 0.2e1;
t619 = t761 + t762;
t803 = t401 + (-mrSges(7,3) / 0.2e1 + mrSges(6,2) / 0.2e1 + t619 * t512) * t446 + t651;
t802 = Ifges(5,6) * t511 + t817;
t801 = (t319 / 0.2e1 - t320 / 0.2e1) * qJ(5) + t816;
t772 = t163 / 0.2e1;
t783 = mrSges(7,1) / 0.2e1;
t785 = -mrSges(6,1) / 0.2e1;
t598 = pkin(5) * t685 - t426;
t92 = t598 + t163;
t800 = mrSges(5,3) * t772 + t121 * t785 + t92 * t783;
t773 = -t162 / 0.2e1;
t592 = -pkin(5) * t684 + t446 * t510;
t74 = -t162 + t592;
t781 = t74 / 0.2e1;
t784 = mrSges(6,1) / 0.2e1;
t799 = mrSges(7,1) * t781 + mrSges(5,3) * t773 + t122 * t784;
t506 = t511 * pkin(4);
t442 = t506 - t576;
t749 = m(7) * t442;
t456 = t506 - t710;
t750 = m(6) * t456;
t798 = -t665 * t749 / 0.2e1 - t663 * t750 / 0.2e1;
t797 = -Ifges(4,1) + Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t742 = mrSges(6,1) + mrSges(7,1);
t796 = qJ(5) * t742 + Ifges(5,6);
t765 = t445 / 0.2e1;
t795 = mrSges(7,1) * t627 + t716 / 0.2e1 - t718 / 0.2e1 + (mrSges(7,1) * t765 + t446 * t619) * t512;
t727 = Ifges(5,3) * t446;
t738 = Ifges(7,1) * t446;
t739 = Ifges(6,1) * t446;
t768 = -t319 / 0.2e1;
t794 = Ifges(5,6) * t628 + Ifges(6,4) * t626 - t727 / 0.2e1 - t738 / 0.2e1 - t739 / 0.2e1 - t816 + t846 * t629 + t835 * t627 + (t768 + t320 / 0.2e1) * qJ(5);
t793 = m(5) / 0.2e1;
t792 = -m(6) / 0.2e1;
t790 = -m(7) / 0.2e1;
t786 = m(7) * pkin(2);
t782 = mrSges(6,3) / 0.2e1;
t174 = t323 * t512 + t827;
t80 = -t174 + t592;
t780 = -t80 / 0.2e1;
t175 = t511 * t323 - t826;
t98 = t598 + t175;
t779 = -t98 / 0.2e1;
t778 = m(7) * t92;
t777 = m(7) * t98;
t776 = t112 / 0.2e1;
t775 = -t122 / 0.2e1;
t129 = -t174 + t743;
t774 = -t129 / 0.2e1;
t290 = t454 * t446;
t771 = t290 / 0.2e1;
t294 = t446 * t457;
t770 = t294 / 0.2e1;
t316 = mrSges(7,1) * t683 + t445 * mrSges(7,2);
t769 = -t316 / 0.2e1;
t507 = t512 * pkin(5);
t652 = t751 * pkin(2);
t494 = t652 + pkin(9);
t678 = t494 * t512;
t444 = t507 + t678;
t767 = t444 / 0.2e1;
t766 = -t445 / 0.2e1;
t763 = -t590 / 0.2e1;
t760 = -t591 / 0.2e1;
t745 = pkin(9) * t512;
t480 = t507 + t745;
t759 = t480 / 0.2e1;
t748 = pkin(3) * t293;
t747 = pkin(3) * t459;
t737 = Ifges(4,4) * t446;
t735 = Ifges(6,4) * t445;
t731 = Ifges(5,6) * t445;
t204 = t446 * t589 + t731;
t206 = Ifges(5,5) * t445 - t446 * t471;
t207 = Ifges(7,5) * t445 + t446 * t583;
t209 = Ifges(6,5) * t445 + t446 * t584;
t402 = Ifges(7,6) * t682;
t211 = Ifges(7,4) * t445 - Ifges(7,2) * t683 - t402;
t403 = Ifges(6,6) * t683;
t213 = Ifges(6,2) * t682 - t403 + t735;
t311 = -mrSges(5,2) * t445 + mrSges(5,3) * t683;
t313 = t445 * mrSges(5,1) + mrSges(5,3) * t682;
t649 = mrSges(7,1) * t682;
t719 = t445 * mrSges(7,3);
t314 = -t649 - t719;
t315 = -mrSges(6,1) * t683 - t445 * mrSges(6,3);
t317 = -mrSges(6,1) * t682 + t445 * mrSges(6,2);
t429 = Ifges(4,4) * t445;
t4 = t446 * (-Ifges(4,2) * t445 - t737) / 0.2e1 + m(5) * (-t148 * t162 + t149 * t163 + t693) + (-0.2e1 * t429 + (-Ifges(4,1) + Ifges(4,2)) * t446) * t766 + (Ifges(3,1) - Ifges(3,2)) * t754 * t752 + t74 * t314 + t121 * t315 + t92 * t316 + t122 * t317 + t163 * t311 + t162 * t313 + (t445 * t797 + t446 * t802 + t737) * t764 + (t445 * t802 - t727 - t738 - t739) * t765 - t842 + (m(4) * t496 + mrSges(4,1) * t445 - mrSges(4,2) * t446) * t654 + t204 * t628 + t213 * t626 + (t211 + t209) * t629 + (t207 + t206) * t627 + (-t752 ^ 2 + t754 ^ 2) * Ifges(3,4) + m(6) * (t118 * t121 + t119 * t122 + t832) - pkin(1) * (mrSges(3,1) * t752 + mrSges(3,2) * t754) + m(7) * (t70 * t74 + t90 * t92 + t845) - t838 + (t214 / 0.2e1 - t823 / 0.2e1) * t682 + (t203 / 0.2e1 - t822 / 0.2e1) * t683;
t720 = t4 * qJD(1);
t715 = t511 * mrSges(7,1);
t498 = t512 * mrSges(6,1);
t497 = t512 * mrSges(7,1);
t128 = t426 - t175;
t593 = t206 / 0.2e1 + t207 / 0.2e1 - t213 / 0.2e1;
t594 = -t204 / 0.2e1 + t209 / 0.2e1 + t211 / 0.2e1;
t644 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t599 = -Ifges(5,6) / 0.2e1 + t644;
t643 = Ifges(5,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t600 = Ifges(6,4) / 0.2e1 - t643;
t8 = -t128 * t315 - t98 * t316 - t129 * t317 - t175 * t311 - t174 * t313 - t80 * t314 - m(5) * (-t148 * t174 + t149 * t175 + t693) - m(6) * (t118 * t128 + t119 * t129 + t832) - m(7) * (t70 * t80 + t90 * t98 + t845) + (t721 + t737 + (-t214 / 0.2e1 + t208 / 0.2e1 + t205 / 0.2e1 + t600 * t446) * t512 + (t212 / 0.2e1 + t210 / 0.2e1 - t203 / 0.2e1 - t599 * t446) * t511 + (Ifges(4,2) + t797) * t445) * t446 + (-t429 + (-t445 * t600 + t593) * t512 + (t445 * t599 + t594) * t511) * t445 + t838;
t714 = t8 * qJD(1);
t288 = -pkin(4) * t682 - t446 * t674;
t196 = -qJ(6) * t682 + t288;
t291 = t455 * t446;
t296 = Ifges(7,3) * t683 - t402;
t297 = -Ifges(6,3) * t682 - t403;
t298 = t446 * t465;
t299 = t446 * t466;
t300 = t446 * t468;
t301 = t446 * t470;
t405 = Ifges(7,4) * t682;
t664 = -Ifges(6,4) * t683 - Ifges(6,5) * t682;
t9 = t405 * t766 + t664 * t765 + t359 * t291 + t113 * t314 + t112 * t316 + t180 * t290 + t126 * t294 + t288 * t295 + t196 * t289 + (t317 - t313) * t149 + (t315 - t311) * t148 + m(6) * (t118 * t148 + t119 * t149 + t180 * t288) + m(7) * (t112 * t90 + t113 * t70 + t126 * t196) + ((t731 / 0.2e1 + t149 * mrSges(5,3) - t118 * mrSges(6,1) + t90 * mrSges(7,1) - t296 / 0.2e1 + t299 / 0.2e1 - t301 / 0.2e1 - t594) * t512 + (t148 * mrSges(5,3) + t119 * mrSges(6,1) + t70 * mrSges(7,1) + t300 / 0.2e1 - t297 / 0.2e1 - t298 / 0.2e1 + t643 * t445 + t593) * t511) * t446;
t713 = t9 * qJD(1);
t712 = t112 + t70;
t711 = t113 - t90;
t19 = (-m(6) * t180 - m(7) * t126 - t809) * t682 + (m(6) * t118 - m(7) * t90 + t315 - t316) * t445;
t709 = qJD(1) * t19;
t672 = t511 * t126;
t34 = -t445 * t314 - t289 * t683 + m(7) * (-t445 * t70 - t446 * t672);
t708 = qJD(1) * t34;
t704 = t126 * t512;
t702 = t128 * t512;
t701 = t129 * t511;
t698 = t174 * t511;
t697 = t175 * t512;
t690 = t386 * t287;
t689 = t421 * t287;
t688 = t424 * t292;
t443 = (pkin(5) + t494) * t511;
t687 = t443 * t318;
t686 = t444 * t320;
t681 = t452 * t292;
t478 = (pkin(5) + pkin(9)) * t511;
t680 = t478 * t318;
t679 = t480 * t320;
t677 = t495 * t293;
t676 = t495 * t459;
t671 = t511 * t180;
t670 = t511 * t312;
t669 = t511 * t321;
t668 = t512 * t310;
t667 = t512 * t319;
t666 = -t442 * t457 - t456 * t454;
t662 = t818 * t494;
t660 = t818 * pkin(9);
t659 = t497 + t498;
t658 = -t753 / 0.2e1;
t657 = t753 / 0.2e1;
t648 = m(7) * t776;
t642 = t494 * t670;
t641 = t494 * t669;
t640 = t494 * t668;
t639 = t494 * t667;
t638 = -t735 / 0.4e1;
t637 = -t731 / 0.4e1;
t632 = t511 * t753;
t631 = t512 * t753;
t621 = t289 * t756;
t620 = -t386 / 0.2e1 - t421 / 0.2e1;
t617 = m(6) * t663;
t616 = m(7) * t665;
t615 = -t591 + t750;
t614 = m(7) * t386 + t457;
t613 = m(7) * t421 + t457;
t612 = -t590 + t749;
t609 = pkin(2) * t632;
t607 = t655 / 0.2e1;
t601 = Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1 + Ifges(5,3) / 0.2e1;
t597 = -t616 / 0.2e1;
t596 = (-t289 / 0.2e1 - t295 / 0.2e1) * t511;
t536 = -t666 + (-t585 + t468) * t758 + (-t584 + t464 + t466) * t756 + (-t589 + t820) * t755 + (-t583 + t471 + t819) * t757;
t26 = t386 * t612 + t424 * t615 + t536 + t676;
t140 = t456 * t180;
t520 = (t501 / 0.4e1 + t502 / 0.4e1 + t503 / 0.4e1 + t500 / 0.4e1) * t445 + t126 * t763 + t180 * t760 + t196 * t761 + t288 * t762 + t359 * t459 / 0.2e1 + t442 * t289 / 0.2e1 + t456 * t295 / 0.2e1;
t99 = t442 * t126;
t514 = (t424 * t288 + t140) * t791 + (t196 * t386 + t443 * t711 + t444 * t712 + t99) * t789 + t386 * t770 + t424 * t771 + t443 * t769 + t314 * t767 + t495 * t291 / 0.2e1 + t520;
t519 = (-pkin(4) * t122 - qJ(5) * t121) * t792 + (qJ(5) * t92 - t510 * t74) * t790 + t121 * t782 + mrSges(6,2) * t775 + mrSges(5,1) * t773 + mrSges(5,2) * t772 + mrSges(7,3) * t781 - t92 * mrSges(7,2) / 0.2e1;
t528 = (t113 / 0.2e1 - t90 / 0.2e1) * mrSges(7,1) + (t149 / 0.2e1 + t118 / 0.2e1) * mrSges(6,1) - t204 / 0.4e1 + t209 / 0.4e1 + t211 / 0.4e1 + t296 / 0.4e1 - t299 / 0.4e1 + t301 / 0.4e1;
t558 = (t118 + t149) * t791 + t315 / 0.2e1 - t311 / 0.2e1;
t521 = t494 * t558 + t528;
t529 = (t776 + t70 / 0.2e1) * mrSges(7,1) + (-t148 / 0.2e1 + t119 / 0.2e1) * mrSges(6,1) + t206 / 0.4e1 + t207 / 0.4e1 - t213 / 0.4e1 - t297 / 0.4e1 - t298 / 0.4e1 + t300 / 0.4e1;
t557 = (t119 - t148) * t791 + t317 / 0.2e1 - t313 / 0.2e1;
t522 = t494 * t557 + t529;
t551 = t583 / 0.4e1 - t462 / 0.4e1 - t465 / 0.4e1 - t585 / 0.4e1 + t468 / 0.4e1 - t471 / 0.4e1;
t541 = mrSges(7,1) * t767 + t551;
t550 = t461 / 0.4e1 + t584 / 0.4e1 - t464 / 0.4e1 - t466 / 0.4e1 - t589 / 0.4e1 + t470 / 0.4e1;
t542 = t443 * t783 + t550;
t546 = t811 * (t508 / 0.2e1 + t509 / 0.2e1);
t543 = t546 * t494;
t562 = (-0.3e1 / 0.4e1 * Ifges(5,6) + t644) * t445;
t563 = (-0.3e1 / 0.4e1 * Ifges(6,4) + t643) * t445;
t3 = (t543 + t601) * t446 + t514 + (t446 * t541 + t522 + t563) * t512 + (t446 * t542 + t521 + t562) * t511 + t519 + t801;
t579 = t3 * qJD(1) + t26 * qJD(2);
t517 = -mrSges(4,2) * t655 + (-mrSges(4,1) + t455 + t661) * t652 + t818 * (mrSges(7,1) + t811);
t38 = (t386 * t751 + t443 * t632 + t444 * t631) * t786 + t517 + m(6) * (t424 * t652 + t662) + m(5) * (t495 * t652 + t662);
t572 = t697 - t698;
t574 = t701 - t702;
t513 = -t640 / 0.2e1 + t639 / 0.2e1 + t677 / 0.2e1 + t497 * t779 - (-t560 + t809) * t652 / 0.2e1 - t848 + t715 * t780 + t702 * t784 + t701 * t785 - ((t317 + t314) * t511 + (t316 + t311) * t512) * t655 / 0.2e1 - t839 - t690 / 0.2e1 - t688 / 0.2e1 - t686 / 0.2e1 - t687 / 0.2e1 - m(5) * (t828 + t572 * t494 + (t148 * t632 + t149 * t631 + t831) * pkin(2)) / 0.2e1 + (t830 + t574 * t494 + (-t118 * t631 + t119 * t632 + t180 * t751) * pkin(2)) * t792 + (t844 + t443 * t80 + t444 * t98 + (t126 * t751 + t631 * t90 + t632 * t70) * pkin(2)) * t790 + (t313 * t511 + t315 * t512) * t607 + (-t697 / 0.2e1 + t698 / 0.2e1) * mrSges(5,3) - t641 / 0.2e1 + t642 / 0.2e1;
t516 = (pkin(9) * t573 - t836) * t793 + (pkin(9) * t575 + t829) * t791 + (t478 * t74 + t480 * t92 + t843) * t789 + t748 / 0.2e1 + t689 / 0.2e1 + t681 / 0.2e1 + t680 / 0.2e1 + t679 / 0.2e1 + t848;
t7 = t516 + t513 - (t667 + t670) * pkin(9) / 0.2e1 + (t668 + t669) * pkin(9) / 0.2e1 + t800 * t512 + t799 * t511 + t839;
t578 = -t7 * qJD(1) + t38 * qJD(2);
t302 = m(7) * t445;
t577 = m(7) * qJD(4) + qJD(1) * t302;
t525 = t596 + (-t671 + (t424 * t446 + t445 * t494) * t512) * t791 + (t386 * t682 + t444 * t445 - t672) * t789;
t567 = m(6) * t775 + t74 * t790;
t12 = t525 + t567 + t803;
t143 = (m(6) * t424 + t454 + t614) * t511;
t571 = -qJD(1) * t12 + qJD(2) * t143;
t281 = t614 * t512;
t535 = (-t386 * t683 - t443 * t445 - t704) * t789 + t621;
t30 = -t778 / 0.2e1 + t535 + t806;
t570 = -qJD(1) * t30 + qJD(2) * t281;
t527 = t741 * t445 + (-t118 + t425) * t791 + (0.2e1 * t425 + t113) * t789;
t565 = t113 * t790 + t149 * t792;
t27 = t527 + t565;
t486 = qJ(5) * t812 + t741;
t569 = qJD(1) * t27 + qJD(4) * t486;
t533 = (t112 + 0.2e1 * t810) * t790 - t719;
t39 = t648 + t533;
t485 = m(7) * t510 + mrSges(7,3);
t568 = -qJD(1) * t39 + qJD(4) * t485;
t566 = m(6) * t774 + m(7) * t780;
t559 = t796 * qJD(4) * t511;
t487 = qJ(5) * pkin(2) * t631;
t531 = (-pkin(4) * t609 + t487) * t791 + (-t510 * t609 + t487) * t789;
t537 = (mrSges(5,2) * t658 + t657 * t741) * pkin(2);
t538 = (mrSges(6,2) * t657 + (mrSges(5,1) + mrSges(7,3)) * t658) * pkin(2);
t10 = (-t495 / 0.2e1 + pkin(3) / 0.2e1) * t459 - (-t452 / 0.2e1 - t424 / 0.2e1) * t591 - t620 * t590 + (-t471 / 0.2e1 + t468 / 0.2e1 - t585 / 0.2e1 - t465 / 0.2e1 - t462 / 0.2e1 + t583 / 0.2e1 + t538) * t511 + (-t470 / 0.2e1 + t589 / 0.2e1 + t466 / 0.2e1 + t464 / 0.2e1 - t584 / 0.2e1 - t461 / 0.2e1 + t537) * t512 + t531 + t666 + t798;
t28 = t421 * t612 + t452 * t615 + t536 - t747;
t515 = (t452 * t288 + t140) * t791 + (t196 * t421 + t478 * t711 + t480 * t712 + t99) * t789 - pkin(3) * t291 / 0.2e1 + t421 * t770 + t452 * t771 + t478 * t769 + t314 * t759 + t520;
t518 = (-pkin(4) * t129 - qJ(5) * t128) * t792 + (qJ(5) * t98 - t510 * t80) * t790 + t128 * t782 + mrSges(6,2) * t774 - t174 * mrSges(5,1) / 0.2e1 + t175 * mrSges(5,2) / 0.2e1 + t80 * mrSges(7,3) / 0.2e1 + mrSges(7,2) * t779;
t523 = pkin(9) * t558 + t528;
t524 = pkin(9) * t557 + t529;
t539 = mrSges(7,1) * t759 + t551;
t540 = t478 * t783 + t550;
t544 = t546 * pkin(9);
t5 = (t544 + t601) * t446 + t515 + (t446 * t540 + t523 + t562) * t511 + (t446 * t539 + t524 + t563) * t512 + t518 + t801;
t556 = t5 * qJD(1) - t10 * qJD(2) + t28 * qJD(3);
t526 = t596 + (-t671 + (t446 * t452 + t744) * t512) * t791 + (t421 * t682 + t480 * t445 - t672) * t789;
t14 = t526 + t566 + t803;
t178 = (m(6) * t452 + t454 + t613) * t511;
t552 = t812 * t657 * pkin(2);
t60 = (t552 + t617 / 0.2e1 + t616 / 0.2e1 + t661) * t511;
t555 = -qJD(1) * t14 + qJD(2) * t60 + qJD(3) * t178;
t146 = (t457 + (t607 - t620) * m(7)) * t512;
t306 = t613 * t512;
t534 = (-t421 * t683 - t478 * t445 - t704) * t789 + t621;
t31 = -t777 / 0.2e1 + t534 + t806;
t554 = -qJD(1) * t31 + qJD(2) * t146 + qJD(3) * t306;
t545 = -pkin(4) * t498 - t510 * t497 - t817;
t530 = qJD(4) * (m(6) * t582 + t454 + t455);
t387 = -m(7) * t478 - t715;
t373 = -m(7) * t443 - t715;
t332 = m(6) * t745 + m(7) * t480 + t659;
t264 = m(6) * t678 + m(7) * t444 + t659;
t147 = (m(7) * t607 - t457 + t597) * t512;
t61 = -t617 * t757 + (t597 + t552 - t661) * t511;
t37 = -t533 + t648 + t649;
t36 = t777 / 0.2e1 + t534 + t808;
t32 = t778 / 0.2e1 + t535 + t808;
t24 = t683 * t742 + t527 - t565;
t16 = t526 - t566 + t795;
t13 = t525 - t567 + t795;
t11 = t511 * t538 + t512 * t537 + t536 + t676 / 0.2e1 - t747 / 0.2e1 + t531 + t665 * t763 + t663 * t760 - t798;
t6 = (t637 + t523) * t511 + (t638 + t524) * t512 + t794 + t515 - t518 + (t511 * t540 + t512 * t539 + t544) * t446;
t2 = (t637 + t521) * t511 + (t638 + t522) * t512 + t794 + t514 - t519 + (t511 * t542 + t512 * t541 + t543) * t446;
t1 = ((t847 - t312 / 0.2e1) * pkin(9) + t799) * t511 + ((t768 + t310 / 0.2e1) * pkin(9) + t800) * t512 + t532 + t516 - t513 + (t455 / 0.2e1 - mrSges(4,1) / 0.2e1) * t805;
t15 = [qJD(2) * t4 - qJD(3) * t8 + qJD(4) * t9 - qJD(5) * t19 + qJD(6) * t34, t720 + (-mrSges(3,1) * t656 + t640 - t639 - t677 + m(4) * (-t753 * t805 - t831) * pkin(2) + t575 * mrSges(6,1) + t573 * mrSges(5,3) + t92 * t497 + t74 * t715 + mrSges(3,2) * t653 + t690 + t688 + t686 + t687 + m(6) * (t494 * t575 + t830) + m(5) * (t494 * t573 + t828) + t849 - Ifges(3,6) * t752 + Ifges(3,5) * t754 + m(7) * (t443 * t74 + t444 * t92 + t844) + t641 - t642 + (t445 * t655 + t446 * t652) * mrSges(4,3)) * qJD(2) + t1 * qJD(3) + t2 * qJD(4) + t13 * qJD(5) + t32 * qJD(6), t1 * qJD(2) + t6 * qJD(4) + t16 * qJD(5) + t36 * qJD(6) - t714 + (t679 + t680 + t681 + t689 + t748 + (-t128 * mrSges(6,1) + t98 * mrSges(7,1) + t175 * mrSges(5,3) + (t310 - t319) * pkin(9)) * t512 + (t129 * mrSges(6,1) + t80 * mrSges(7,1) - t174 * mrSges(5,3) + (-t312 + t321) * pkin(9)) * t511 + (pkin(9) * t574 + t829) * t813 + 0.2e1 * (pkin(9) * t572 - t836) * t793 + (t478 * t80 + t480 * t98 + t843) * t814 + t849) * qJD(3), t2 * qJD(2) + t6 * qJD(3) + t24 * qJD(5) + t37 * qJD(6) + t713 + (t112 * mrSges(7,2) - t113 * mrSges(7,3) - t405 + t664 + (qJ(5) * t112 - t113 * t510) * t814 + (t796 * t512 + (-mrSges(6,1) * pkin(4) - mrSges(7,1) * t510 + t835) * t511) * t446 + (-pkin(4) * t813 - mrSges(5,1) + mrSges(6,2)) * t149 + (-qJ(5) * t813 + mrSges(5,2) - mrSges(6,3)) * t148) * qJD(4), qJD(2) * t13 + qJD(3) * t16 + qJD(4) * t24 - t709, qJD(2) * t32 + qJD(3) * t36 + qJD(4) * t37 + t708; -qJD(3) * t7 + qJD(4) * t3 + qJD(5) * t12 + qJD(6) * t30 - t720, qJD(3) * t38 + qJD(4) * t26 - qJD(5) * t143 - qJD(6) * t281 ((t421 * t751 + t478 * t632 + t480 * t631) * t786 + t517 + m(6) * (t452 * t652 + t660) + m(5) * (-pkin(3) * t652 + t660)) * qJD(3) + t11 * qJD(4) + t61 * qJD(5) + t147 * qJD(6) + t578, t11 * qJD(3) + (m(7) * (-qJ(5) * t443 - t444 * t510) - t443 * mrSges(7,2) - t444 * mrSges(7,3) + t545) * qJD(4) + t264 * qJD(5) + t373 * qJD(6) - t559 + t494 * t530 + t579, qJD(3) * t61 + qJD(4) * t264 - t571, qJD(3) * t147 + qJD(4) * t373 - t570; qJD(2) * t7 + qJD(4) * t5 + qJD(5) * t14 + qJD(6) * t31 + t714, -qJD(4) * t10 - qJD(5) * t60 - qJD(6) * t146 - t578, qJD(4) * t28 - qJD(5) * t178 - qJD(6) * t306 (m(7) * (-qJ(5) * t478 - t480 * t510) - t480 * mrSges(7,3) - t478 * mrSges(7,2) + t545) * qJD(4) + t332 * qJD(5) + t387 * qJD(6) - t559 + pkin(9) * t530 + t556, qJD(4) * t332 - t555, qJD(4) * t387 - t554; -qJD(2) * t3 - qJD(3) * t5 + qJD(5) * t27 - qJD(6) * t39 - t713, qJD(3) * t10 - t579, -t556, qJD(5) * t486 + qJD(6) * t485, t569, t568; -qJD(2) * t12 - qJD(3) * t14 - qJD(4) * t27 - qJD(6) * t302 + t709, qJD(3) * t60 + t571, t555, -m(7) * qJD(6) - t569, 0, -t577; -qJD(2) * t30 - qJD(3) * t31 + qJD(4) * t39 + qJD(5) * t302 - t708, qJD(3) * t146 + t570, t554, m(7) * qJD(5) - t568, t577, 0;];
Cq  = t15;
