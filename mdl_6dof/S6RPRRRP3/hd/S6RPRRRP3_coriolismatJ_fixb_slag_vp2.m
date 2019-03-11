% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:03:07
% EndTime: 2019-03-09 06:03:33
% DurationCPUTime: 13.88s
% Computational Cost: add. (20032->728), mult. (41387->940), div. (0->0), fcn. (39922->8), ass. (0->378)
t458 = cos(qJ(3));
t806 = t458 / 0.4e1;
t800 = Ifges(7,4) + Ifges(6,5);
t798 = Ifges(6,6) - Ifges(7,6);
t456 = sin(qJ(3));
t454 = sin(qJ(5));
t457 = cos(qJ(4));
t455 = sin(qJ(4));
t712 = cos(qJ(5));
t575 = t712 * t455;
t501 = t454 * t457 + t575;
t375 = t501 * t456;
t729 = -t375 / 0.2e1;
t790 = pkin(9) + pkin(8);
t602 = t790 * t457;
t621 = t454 * t455;
t754 = t602 * t712 - t621 * t790;
t736 = t754 / 0.2e1;
t705 = pkin(4) * t457;
t439 = -pkin(3) - t705;
t574 = t712 * t457;
t500 = t574 - t621;
t520 = -pkin(5) * t500 - qJ(6) * t501;
t245 = t439 + t520;
t272 = -mrSges(7,1) * t500 - mrSges(7,3) * t501;
t780 = m(7) * t245 + t272;
t433 = sin(pkin(10)) * pkin(1) + pkin(7);
t416 = t456 * t433;
t619 = t455 * t456;
t384 = pkin(4) * t619 + t416;
t374 = t454 * t619 - t456 * t574;
t653 = qJ(6) * t374;
t523 = pkin(5) * t375 + t653;
t143 = t523 + t384;
t697 = Ifges(6,4) * t374;
t205 = -Ifges(6,2) * t375 - t458 * Ifges(6,6) - t697;
t229 = -t374 * mrSges(7,1) + t375 * mrSges(7,3);
t230 = -t374 * mrSges(6,1) - t375 * mrSges(6,2);
t270 = mrSges(7,1) * t501 - mrSges(7,3) * t500;
t271 = mrSges(6,1) * t501 + mrSges(6,2) * t500;
t663 = t500 * Ifges(7,5);
t274 = Ifges(7,3) * t501 + t663;
t399 = Ifges(7,5) * t501;
t275 = -Ifges(7,3) * t500 + t399;
t402 = Ifges(6,4) * t500;
t278 = -Ifges(6,2) * t501 + t402;
t660 = t501 * Ifges(6,4);
t279 = Ifges(6,2) * t500 + t660;
t280 = Ifges(7,1) * t500 + t399;
t282 = Ifges(6,1) * t500 - t660;
t302 = t454 * t602 + t575 * t790;
t694 = Ifges(7,5) * t375;
t524 = -Ifges(7,3) * t374 - t694;
t537 = t800 * t500 - t798 * t501;
t770 = mrSges(7,2) + mrSges(6,3);
t207 = -Ifges(7,1) * t374 - Ifges(7,4) * t458 + t694;
t353 = Ifges(6,4) * t375;
t209 = -Ifges(6,1) * t374 - Ifges(6,5) * t458 - t353;
t793 = Ifges(6,2) * t374 + t207 + t209 - t353;
t350 = Ifges(7,5) * t374;
t203 = -t458 * Ifges(7,6) + Ifges(7,3) * t375 - t350;
t794 = t203 - t350 + t697 + (-Ifges(6,1) - Ifges(7,1)) * t375;
t281 = Ifges(7,1) * t501 - t663;
t283 = Ifges(6,1) * t501 + t402;
t797 = t283 + t281;
t805 = -t770 * (t302 * t729 + t374 * t736) - t439 * t230 / 0.2e1 - t384 * t271 / 0.2e1 - t375 * t274 / 0.4e1 - t374 * t279 / 0.4e1 - t245 * t229 / 0.2e1 - t143 * t270 / 0.2e1 + t537 * t806 + (t282 + t280 + t275) * t374 / 0.4e1 + (t278 + t797) * t375 / 0.4e1 - (-t205 / 0.4e1 + t794 / 0.4e1) * t501 - (-t524 / 0.4e1 + t793 / 0.4e1) * t500;
t804 = -mrSges(6,1) / 0.2e1;
t803 = mrSges(6,2) / 0.2e1;
t802 = -mrSges(7,3) / 0.2e1;
t799 = Ifges(7,2) + Ifges(6,3);
t664 = t500 * mrSges(7,2);
t580 = -t664 / 0.2e1;
t581 = t664 / 0.2e1;
t258 = t581 + t580;
t796 = qJD(3) * t258;
t795 = qJD(6) * t258;
t791 = -(-t281 / 0.2e1 + t274 / 0.2e1 - t283 / 0.2e1 - t278 / 0.2e1) * t500 - (-t280 / 0.2e1 - t275 / 0.2e1 - t282 / 0.2e1 + t279 / 0.2e1) * t501 + t245 * t270 + t439 * t271;
t713 = t458 / 0.2e1;
t784 = -t754 / 0.2e1;
t452 = t457 ^ 2;
t601 = t455 ^ 2 + t452;
t782 = mrSges(5,3) * t601;
t577 = -cos(pkin(10)) * pkin(1) - pkin(2);
t702 = t456 * pkin(8);
t396 = -pkin(3) * t458 + t577 - t702;
t380 = t457 * t396;
t617 = t456 * t457;
t592 = pkin(9) * t617;
t225 = -t592 + t380 + (-t433 * t455 - pkin(4)) * t458;
t417 = t458 * t433;
t620 = t455 * t396;
t503 = t417 * t457 + t620;
t238 = -pkin(9) * t619 + t503;
t622 = t454 * t238;
t102 = t225 * t712 - t622;
t701 = t458 * pkin(5);
t92 = -t102 + t701;
t781 = t102 + t92;
t445 = m(7) * qJ(6) + mrSges(7,3);
t778 = qJD(5) * t445;
t777 = t445 * qJD(6);
t703 = t456 * pkin(3);
t704 = pkin(8) * t458;
t420 = t703 - t704;
t305 = t455 * t416 + t457 * t420;
t615 = t457 * t458;
t243 = t456 * pkin(4) - pkin(9) * t615 + t305;
t306 = t455 * t420 - t433 * t617;
t618 = t455 * t458;
t257 = -pkin(9) * t618 + t306;
t122 = t454 * t243 + t712 * t257;
t114 = qJ(6) * t456 + t122;
t121 = t243 * t712 - t454 * t257;
t115 = -t456 * pkin(5) - t121;
t775 = t122 * t803 + t121 * t804 + t115 * mrSges(7,1) / 0.2e1 + t114 * t802;
t671 = t375 * mrSges(7,2);
t334 = -t671 / 0.2e1;
t774 = 0.2e1 * t334;
t773 = -t272 / 0.2e1;
t707 = pkin(4) * t454;
t587 = t707 / 0.2e1;
t771 = mrSges(7,1) + mrSges(6,1);
t491 = mrSges(5,3) * t503;
t576 = t712 * t238;
t623 = t454 * t225;
t103 = t576 + t623;
t614 = t458 * qJ(6);
t91 = t103 - t614;
t769 = t229 + t230;
t377 = t501 * t458;
t307 = -mrSges(7,2) * t377 + mrSges(7,3) * t456;
t310 = -mrSges(6,2) * t456 - mrSges(6,3) * t377;
t768 = t307 + t310;
t446 = t458 * mrSges(7,3);
t308 = -t446 - t671;
t670 = t375 * mrSges(6,3);
t309 = t458 * mrSges(6,2) - t670;
t767 = t308 + t309;
t447 = Ifges(5,4) * t457;
t693 = Ifges(5,2) * t455;
t766 = t447 - t693;
t589 = mrSges(5,3) * t619;
t508 = mrSges(5,2) * t458 - t589;
t509 = -mrSges(5,1) * t458 - mrSges(5,3) * t617;
t715 = t457 / 0.2e1;
t718 = t455 / 0.2e1;
t765 = t508 * t718 + t509 * t715;
t591 = t712 * pkin(4);
t764 = t374 * t591 - t375 * t707;
t259 = -t417 * t455 + t380;
t237 = t259 - t592;
t111 = t454 * t237 + t576;
t112 = t237 * t712 - t622;
t721 = -t501 / 0.2e1;
t723 = -t500 / 0.2e1;
t762 = t111 * t721 + t112 * t723;
t538 = t374 * t798 - t375 * t800;
t511 = -t305 * t455 + t306 * t457;
t529 = -t455 * Ifges(5,1) - t447;
t759 = -mrSges(5,1) * t457 + mrSges(5,2) * t455;
t758 = -t308 / 0.2e1 - t309 / 0.2e1;
t757 = (t270 + t271) * t713;
t231 = mrSges(7,1) * t375 + mrSges(7,3) * t374;
t756 = -m(7) * t143 - t231;
t755 = t782 / 0.2e1;
t753 = -m(7) * t91 - t767;
t672 = t374 * mrSges(6,3);
t311 = -t458 * mrSges(6,1) + t672;
t673 = t374 * mrSges(7,2);
t312 = t458 * mrSges(7,1) - t673;
t752 = -m(7) * t92 + t311 - t312;
t717 = t456 / 0.2e1;
t379 = t500 * t458;
t725 = t379 / 0.2e1;
t727 = t377 / 0.2e1;
t728 = -t377 / 0.2e1;
t751 = Ifges(6,6) * t728 + Ifges(7,6) * t727 + t717 * t799 + t725 * t800 - t775;
t748 = t456 ^ 2;
t747 = 0.2e1 * qJD(3);
t746 = -m(6) / 0.2e1;
t745 = m(6) / 0.2e1;
t744 = -m(7) / 0.2e1;
t743 = m(7) / 0.2e1;
t742 = m(6) * pkin(4);
t741 = m(7) * pkin(4);
t738 = -qJ(6) / 0.2e1;
t737 = -t231 / 0.2e1;
t658 = t456 * mrSges(7,1);
t666 = t379 * mrSges(7,2);
t314 = -t658 + t666;
t733 = t314 / 0.2e1;
t732 = -t374 / 0.2e1;
t731 = t374 / 0.2e1;
t730 = t375 / 0.2e1;
t726 = -t379 / 0.2e1;
t724 = t500 / 0.2e1;
t722 = t501 / 0.2e1;
t434 = qJ(6) + t707;
t720 = t434 / 0.2e1;
t719 = -t455 / 0.2e1;
t714 = -t458 / 0.2e1;
t711 = m(5) * t433;
t710 = m(7) * t111;
t709 = m(7) * t754;
t708 = m(7) * t374;
t706 = pkin(4) * t455;
t698 = Ifges(5,4) * t455;
t696 = Ifges(5,5) * t456;
t695 = Ifges(5,5) * t458;
t692 = Ifges(5,6) * t456;
t691 = t102 * mrSges(6,2);
t690 = t102 * mrSges(7,3);
t689 = t103 * mrSges(6,1);
t688 = t103 * mrSges(7,1);
t687 = t111 * mrSges(6,1);
t686 = t111 * mrSges(7,1);
t685 = t112 * mrSges(6,2);
t684 = t112 * mrSges(7,3);
t679 = t754 * mrSges(6,1);
t678 = t754 * mrSges(7,1);
t677 = t302 * mrSges(6,2);
t676 = t302 * mrSges(7,3);
t675 = t305 * mrSges(5,1);
t674 = t306 * mrSges(5,2);
t669 = t377 * mrSges(6,1);
t668 = t377 * mrSges(7,1);
t667 = t379 * mrSges(6,2);
t665 = t379 * mrSges(7,3);
t662 = t501 * mrSges(7,2);
t661 = t501 * mrSges(6,3);
t657 = t457 * mrSges(5,2);
t656 = t103 - t91;
t655 = qJD(6) * t708;
t654 = t759 - mrSges(4,1);
t465 = t311 * t731 + t312 * t732 + t769 * t714 + t767 * t729 + t770 * (-t375 ^ 2 / 0.2e1 - t374 ^ 2 / 0.2e1);
t490 = t759 * t717;
t593 = pkin(4) * t617;
t545 = t458 * t593;
t12 = -t465 + (-t545 + (t111 - t91 - t614) * t375 + (-t112 - t92 + t701) * t374) * t744 + (-t545 + (-t103 + t111) * t375 + (t102 - t112) * t374) * t746 - t458 * t490 + t748 * t755 + t765 * t456;
t652 = t12 * qJD(1);
t522 = pkin(5) * t374 - qJ(6) * t375;
t13 = (-t781 * t374 + t656 * t375 + t522 * t458) * t743 + t465;
t651 = t13 * qJD(1);
t40 = m(7) * (t143 * t374 - t458 * t91) - t458 * t308 + t374 * t231;
t631 = t40 * qJD(1);
t626 = t434 * t501;
t438 = -t591 - pkin(5);
t625 = t438 * t500;
t616 = t456 * t458;
t606 = -t438 * t374 - t434 * t375;
t603 = 0.2e1 * t581;
t431 = pkin(4) * t618;
t385 = t417 + t431;
t600 = qJD(3) * t458;
t597 = -t748 / 0.2e1;
t596 = mrSges(6,3) * t707;
t595 = t742 / 0.2e1;
t590 = t434 * t673;
t588 = -t707 / 0.2e1;
t586 = -mrSges(7,1) / 0.2e1 + t804;
t585 = t803 + t802;
t584 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t583 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t579 = -t662 / 0.2e1;
t578 = -t661 / 0.2e1;
t573 = t102 * t724;
t572 = t103 * t722;
t565 = t501 * t714;
t414 = t457 * Ifges(5,2) + t698;
t564 = t414 * t719;
t563 = t619 / 0.2e1;
t561 = -t618 / 0.2e1;
t560 = -t617 / 0.2e1;
t559 = t617 / 0.2e1;
t558 = -t615 / 0.2e1;
t557 = t500 * t713;
t551 = t374 * t596;
t548 = t302 * t377 + t379 * t754;
t546 = mrSges(6,3) * t591;
t541 = t591 / 0.2e1;
t540 = mrSges(6,3) * t572;
t539 = -t668 / 0.2e1 - t669 / 0.2e1 + t665 / 0.2e1 - t667 / 0.2e1;
t536 = t375 * t546;
t534 = mrSges(5,1) * t455 + t657;
t533 = mrSges(6,1) * t375 - mrSges(6,2) * t374;
t532 = t667 + t669;
t531 = -t665 + t668;
t530 = Ifges(5,1) * t457 - t698;
t526 = Ifges(5,5) * t457 - Ifges(5,6) * t455;
t525 = Ifges(5,5) * t455 + Ifges(5,6) * t457;
t521 = -pkin(5) * t377 + qJ(6) * t379;
t269 = pkin(5) * t501 - qJ(6) * t500;
t144 = -t521 + t385;
t204 = Ifges(7,5) * t379 + t456 * Ifges(7,6) + Ifges(7,3) * t377;
t206 = Ifges(6,4) * t379 - Ifges(6,2) * t377 + t456 * Ifges(6,6);
t208 = Ifges(7,1) * t379 + Ifges(7,4) * t456 + Ifges(7,5) * t377;
t210 = Ifges(6,1) * t379 - Ifges(6,4) * t377 + Ifges(6,5) * t456;
t313 = mrSges(6,1) * t456 - mrSges(6,3) * t379;
t370 = t458 * t766 + t692;
t371 = t458 * t530 + t696;
t479 = t377 * t583 + t379 * t584;
t6 = t102 * t313 + t103 * t310 + t114 * t308 + t115 * t312 + t121 * t311 + t122 * t309 + t144 * t231 + t91 * t307 + t92 * t314 + (-t143 * mrSges(7,3) + t207 / 0.2e1 + t209 / 0.2e1 + t384 * mrSges(6,2)) * t379 + (t143 * mrSges(7,1) + t384 * mrSges(6,1) + t203 / 0.2e1 - t205 / 0.2e1) * t377 + (t385 * mrSges(6,1) - t206 / 0.2e1 + t204 / 0.2e1) * t375 + (-t385 * mrSges(6,2) - t208 / 0.2e1 - t210 / 0.2e1) * t374 + m(5) * (t259 * t305 + t306 * t620) + m(6) * (t102 * t121 + t103 * t122 + t384 * t385) + m(7) * (t114 * t91 + t115 * t92 + t143 * t144) + (t577 * mrSges(4,1) + t259 * mrSges(5,1) - Ifges(4,4) * t456 + (t696 / 0.2e1 + t371 / 0.2e1 - t305 * mrSges(5,3)) * t457 - t583 * t375 + t584 * t374 + (-t396 * mrSges(5,2) - t306 * mrSges(5,3) - t370 / 0.2e1 - t692 / 0.2e1) * t455) * t456 + (Ifges(4,4) * t458 - t675 + t674 + t577 * mrSges(4,2) + (-t259 * mrSges(5,3) + t306 * t711 - t695) * t457 + (-Ifges(5,3) + Ifges(5,1) * t452 / 0.2e1 - Ifges(4,2) + Ifges(4,1) + (t657 + t711) * t433 - t799) * t456 + ((-mrSges(5,3) * t433 * t457 + Ifges(5,6)) * t458 - mrSges(5,3) * t620 + (0.2e1 * t433 * mrSges(5,1) - t447 + t693 / 0.2e1) * t456) * t455 + t479) * t458;
t495 = t456 * t533;
t9 = t509 * t618 / 0.2e1 + (t456 * mrSges(5,1) - mrSges(5,3) * t615) * t563 + t508 * t558 + (-t456 * mrSges(5,2) - mrSges(5,3) * t618) * t560 + (-t114 * t374 + t115 * t375 + t143 * t456 - t144 * t458 + t92 * t377 + t91 * t379) * t744 + t456 * t737 - t495 / 0.2e1 + (-t102 * t377 + t103 * t379 - t121 * t375 - t122 * t374 + t384 * t456 - t385 * t458) * t746 + t311 * t727 + t313 * t730 + t312 * t728 + t314 * t729 - m(5) * (((-0.1e1 + t452) * t417 + (-t259 + t380) * t455) * t458 + (t511 + t416) * t456) / 0.2e1 + t768 * t731 + t767 * t726 + (t531 + t532) * t713 + (t597 + t458 ^ 2 / 0.2e1) * t534;
t519 = t6 * qJD(1) - t9 * qJD(2);
t462 = -t102 * t670 - t103 * t672 - t143 * t229 + t205 * t732 - t384 * t230 + t524 * t729 + t538 * t713 + t92 * t671 - t91 * t673 + t730 * t793 + t731 * t794;
t480 = -Ifges(5,6) * t458 + t456 * t766;
t494 = t530 * t456;
t481 = t494 - t695;
t484 = -t522 + t593;
t7 = t462 + t480 * t559 + t481 * t563 - t525 * t616 / 0.2e1 - t495 * t705 + t457 * t529 * t597 + t503 * t509 - m(6) * t384 * t593 + t491 * t617 + (t433 * t759 + t564) * t748 + t756 * t484 + (-t508 - t589) * t259 + (-m(6) * t103 + t753) * t112 + (m(6) * t102 + t752) * t111;
t516 = -t7 * qJD(1) - t12 * qJD(2);
t8 = t102 * t753 + t103 * t752 - t522 * t756 + t462;
t515 = -t8 * qJD(1) + t13 * qJD(2);
t47 = m(5) * (-0.1e1 + t601) * t616 + (m(7) + m(6)) * (-t374 * t379 + t375 * t377 - t616);
t514 = -t9 * qJD(1) + t47 * qJD(2);
t513 = t302 * t111 + t112 * t754;
t510 = t115 * t744 + t658 / 0.2e1;
t507 = -t258 * t375 + t374 * t578 - t661 * t732 - t757;
t506 = m(7) * (-pkin(5) * t754 - qJ(6) * t302);
t505 = m(7) * t522;
t504 = m(7) * t521;
t250 = t269 + t706;
t273 = -mrSges(6,1) * t500 + mrSges(6,2) * t501;
t459 = pkin(4) * t273 * t560 + t484 * t773 - t533 * t706 / 0.2e1 - t534 * t416 / 0.2e1 + t455 * t480 / 0.4e1 + t414 * t559 - pkin(3) * t490 + t526 * t806 + t312 * t784 + t311 * t736 + t250 * t737 + (t250 * t143 + t245 * t484 + t754 * t92 + t513) * t744 + (-t102 * t754 + (t384 * t455 + t439 * t617) * pkin(4) + t513) * t746 - (t494 + t481) * t457 / 0.4e1 + t702 * t755 - (t103 * t746 + t744 * t91 + t758) * t302 + t765 * pkin(8) + t762 * mrSges(7,2) + (t573 + t762) * mrSges(6,3) + (t719 + t718) * t491 + (-t529 / 0.2e1 + t766 / 0.4e1) * t619;
t463 = t751 + (t121 * t712 + t122 * t454) * t595 + Ifges(5,5) * t615 / 0.2e1 + t310 * t587 + Ifges(5,6) * t561 + t313 * t541 - t674 / 0.2e1 + t675 / 0.2e1 + Ifges(5,3) * t717 + t307 * t720 + t438 * t733 + (t114 * t434 + t115 * t438) * t743;
t1 = t463 + t459 + t91 * t662 / 0.2e1 + t92 * t580 + t540 + t805;
t475 = -pkin(3) * t534 - t529 * t715 + t564;
t14 = t530 * t718 + t766 * t715 + t475 + t791 + (m(6) * t439 + t273) * t706 + t780 * t250;
t464 = -t250 * t458 * t744 - t431 * t746 + t534 * t713;
t466 = (t377 * t438 + t379 * t434) * t743 + (-t377 * t712 + t379 * t454) * t595 + mrSges(5,1) * t561 + mrSges(5,2) * t558 + t539;
t16 = t464 + t466 + t757;
t489 = -t1 * qJD(1) - t16 * qJD(2) + t14 * qJD(3);
t15 = t780 * t269 + t791;
t468 = -t269 * t458 * t743 + t507;
t19 = t585 * t379 - t586 * t377 - t504 / 0.2e1 + t468;
t467 = t103 * t578 + t91 * t579 + t92 * t581 - t805;
t460 = t467 + t758 * t302 + (t572 + t573) * mrSges(7,2) + (t143 * t269 - t245 * t522 + t302 * t656) * t743 + t269 * t231 / 0.2e1 + t522 * t773 + t540 + (-t311 / 0.2e1 + t312 / 0.2e1 + t781 * t743) * t754;
t473 = (-pkin(5) * t115 + qJ(6) * t114) * t744 + pkin(5) * t733 + t307 * t738;
t4 = t460 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t456 + t473 + t479 + t775;
t488 = t4 * qJD(1) + t19 * qJD(2) + t15 * qJD(3);
t227 = (-t565 + t728) * m(7);
t471 = (-t143 * t501 + t245 * t374 - t458 * t754) * t743 + t272 * t731 + t231 * t721;
t35 = (-t557 + t726) * mrSges(7,2) + t471 + t510;
t89 = t780 * t501;
t486 = qJD(1) * t35 + qJD(2) * t227 - qJD(3) * t89;
t411 = m(7) * t434 + mrSges(7,3);
t472 = -t446 + ((-qJ(6) - t434) * t458 + t103) * t743;
t42 = -t710 / 0.2e1 + t472;
t485 = -qJD(1) * t42 - qJD(4) * t411 + t796;
t461 = (t434 * t102 + t438 * t103 + (t454 * t92 + t712 * t91) * pkin(4)) * t743 - t691 / 0.2e1 + t690 / 0.2e1 - t689 / 0.2e1 - t688 / 0.2e1 + t590 / 0.2e1 + t438 * t334 + t311 * t588 + t312 * t587 + t551 / 0.2e1 + t536 / 0.2e1 + t767 * t541;
t469 = (-pkin(5) * t111 + qJ(6) * t112) * t744 + t687 / 0.2e1 + t686 / 0.2e1 + t685 / 0.2e1 - t684 / 0.2e1 + pkin(5) * t334 - mrSges(7,2) * t653 / 0.2e1;
t11 = t461 + t469;
t476 = -t771 * t707 + (mrSges(7,3) - mrSges(6,2)) * t591;
t174 = -(t434 * t712 + t438 * t454) * t741 - t476;
t470 = ((-t434 + t707) * t302 + (t438 + t591) * t754) * t743;
t23 = -t506 / 0.2e1 + (-(t738 + t720 + t588) * t501 - (-pkin(5) / 0.2e1 - t438 / 0.2e1 - t591 / 0.2e1) * t500) * mrSges(7,2) + t470 + t771 * (t736 + t784);
t474 = (t606 - t764) * t743;
t48 = -t505 / 0.2e1 + t474;
t478 = t11 * qJD(1) + t48 * qJD(2) + t23 * qJD(3) - t174 * qJD(4);
t44 = -t446 + 0.2e1 * (t576 / 0.4e1 - t614 / 0.2e1 + t623 / 0.4e1 - t103 / 0.4e1) * m(7);
t477 = qJD(1) * t44 + qJD(4) * t445 + t778 + t796;
t393 = mrSges(7,3) + (qJ(6) + 0.2e1 * t587) * m(7);
t226 = (-t565 + t727) * m(7);
t132 = t603 + t709;
t94 = t709 / 0.2e1 + m(7) * t736 + t603;
t43 = 0.2e1 * t743 * t91 - t446 + t774;
t41 = t774 + t710 / 0.2e1 + t472;
t39 = t505 / 0.2e1 + t474 - t769;
t34 = -mrSges(7,2) * t557 + t666 / 0.2e1 + t471 - t510;
t20 = -t678 / 0.2e1 - t676 / 0.2e1 + t677 / 0.2e1 - t679 / 0.2e1 + t506 / 0.2e1 + pkin(5) * t580 + qJ(6) * t579 + t586 * t754 + t585 * t302 + (t625 / 0.2e1 - t626 / 0.2e1 + (t454 * t722 + t712 * t724) * pkin(4)) * mrSges(7,2) + t470 + t537;
t18 = t504 / 0.2e1 + t468 + t539;
t17 = -t464 + t466 + t507;
t10 = t461 - t469 + t538;
t5 = t460 - t473 + t751;
t3 = -t9 * qJD(3) - t12 * qJD(4) + t13 * qJD(5);
t2 = t467 + t463 - t459;
t21 = [qJD(3) * t6 - qJD(4) * t7 - qJD(5) * t8 + qJD(6) * t40, t3, t2 * qJD(4) + t5 * qJD(5) + t34 * qJD(6) + ((-t121 * t302 + t122 * t754 + t385 * t439) * t745 + (t114 * t754 + t115 * t302 + t144 * t245) * t743) * t747 + (Ifges(4,5) + (-m(5) * pkin(3) + t654) * t433 + t475) * t600 + t519 + (t114 * t664 + t115 * t662 - Ifges(4,6) * t456 + t385 * t273 + t245 * t531 + t439 * t532 + t144 * t272 + (m(5) * t511 - t456 * t534) * pkin(8) + t122 * t500 * mrSges(6,3) + mrSges(4,2) * t416 - t121 * t661 + t370 * t715 + t371 * t718 + t204 * t723 + t206 * t724 + t275 * t727 + t279 * t728 + t797 * t725 + (t210 + t208) * t722 + t768 * t754 + (-t313 + t314) * t302 + t511 * mrSges(5,3) + (t500 * t798 + t501 * t800 + t525) * t717) * qJD(3), t2 * qJD(3) + (t551 + t536 - Ifges(5,5) * t619 - Ifges(5,6) * t617 - t438 * t671 + m(7) * (t111 * t438 + t112 * t434) + t590 + t684 + (-t111 * t712 + t112 * t454) * t742 - t686 - t685 - t687 - t503 * mrSges(5,1) - t259 * mrSges(5,2) + t538) * qJD(4) + t10 * qJD(5) + t41 * qJD(6) + t516, t5 * qJD(3) + t10 * qJD(4) + (m(7) * (-pkin(5) * t103 + qJ(6) * t102) + t690 - t688 - t691 - t689 + t523 * mrSges(7,2) + t538) * qJD(5) + t43 * qJD(6) + t515, t34 * qJD(3) + t41 * qJD(4) + t43 * qJD(5) + t631; t3, qJD(3) * t47, t17 * qJD(4) + t18 * qJD(5) + t226 * qJD(6) + (-mrSges(4,2) + t782) * t600 + ((t245 * t456 + t548) * t743 + (t439 * t456 + t548) * t745 + m(5) * (t601 * t704 - t703) / 0.2e1) * t747 + t514 + ((t272 + t273 + t654) * t456 + t770 * (t377 * t501 + t379 * t500)) * qJD(3), -t652 + t17 * qJD(3) + (-mrSges(5,1) * t617 + mrSges(5,2) * t619 - t769) * qJD(4) + t39 * qJD(5) + 0.2e1 * (t606 * t743 + t745 * t764) * qJD(4) - t655, t651 + t18 * qJD(3) + t39 * qJD(4) + (t505 - t769) * qJD(5) - t655, t226 * qJD(3) + 0.2e1 * (-qJD(4) / 0.2e1 - qJD(5) / 0.2e1) * t708; -qJD(4) * t1 + qJD(5) * t4 + qJD(6) * t35 - t519, -qJD(4) * t16 + qJD(5) * t19 + qJD(6) * t227 - t514, qJD(4) * t14 + qJD(5) * t15 - qJD(6) * t89 (m(7) * (-t302 * t434 + t438 * t754) - t678 - t676 + t677 - t679 + (-t302 * t454 - t712 * t754) * t742 - t501 * t596 - t500 * t546 + t526 + t537 + t759 * pkin(8) + (t625 - t626) * mrSges(7,2)) * qJD(4) + t20 * qJD(5) + t94 * qJD(6) + t489, t20 * qJD(4) + (t520 * mrSges(7,2) + t537 + (-m(7) * pkin(5) - t771) * t754 + (mrSges(6,2) - t445) * t302) * qJD(5) + t132 * qJD(6) + t488, qJD(4) * t94 + qJD(5) * t132 + t486; qJD(3) * t1 + qJD(5) * t11 + qJD(6) * t42 - t516, qJD(3) * t16 + qJD(5) * t48 + t652, qJD(5) * t23 - t489 - t795, -qJD(5) * t174 + qJD(6) * t411 ((-pkin(5) * t454 + qJ(6) * t712) * t741 + t476) * qJD(5) + t393 * qJD(6) + t478, qJD(5) * t393 - t485; -qJD(3) * t4 - qJD(4) * t11 + qJD(6) * t44 - t515, -qJD(3) * t19 - qJD(4) * t48 - t651, -qJD(4) * t23 - t488 + t795, -t478 + t777, t777, t477; -t35 * qJD(3) - t42 * qJD(4) - t44 * qJD(5) - t631, -t227 * qJD(3), -t486 + (qJD(4) - qJD(5)) * t258, t485 - t778, -t477, 0;];
Cq  = t21;
