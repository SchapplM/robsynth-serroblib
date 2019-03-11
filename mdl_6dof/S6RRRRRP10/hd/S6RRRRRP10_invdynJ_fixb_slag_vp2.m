% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP10_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP10_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:16:26
% EndTime: 2019-03-10 02:18:18
% DurationCPUTime: 70.22s
% Computational Cost: add. (29756->1101), mult. (70699->1453), div. (0->0), fcn. (56592->14), ass. (0->473)
t456 = sin(qJ(2));
t623 = cos(pkin(6));
t559 = pkin(1) * t623;
t438 = t456 * t559;
t455 = sin(qJ(3));
t459 = cos(qJ(3));
t517 = pkin(3) * t455 - pkin(10) * t459;
t452 = sin(pkin(6));
t460 = cos(qJ(2));
t604 = t452 * t460;
t832 = -(t438 + (pkin(8) + t517) * t604) * qJD(1) + t517 * qJD(3);
t454 = sin(qJ(4));
t458 = cos(qJ(4));
t581 = qJD(1) * t452;
t595 = t459 * t460;
t321 = (-t454 * t595 + t456 * t458) * t581;
t578 = qJD(3) * t459;
t551 = t454 * t578;
t575 = qJD(4) * t458;
t742 = t455 * t575 + t321 + t551;
t538 = t623 * qJD(1);
t525 = pkin(1) * t538;
t555 = t456 * t581;
t363 = -pkin(8) * t555 + t460 * t525;
t488 = (pkin(2) * t456 - pkin(9) * t460) * t452;
t364 = qJD(1) * t488;
t255 = t459 * t363 + t455 * t364;
t230 = pkin(10) * t555 + t255;
t404 = -pkin(3) * t459 - pkin(10) * t455 - pkin(2);
t577 = qJD(4) * t454;
t579 = qJD(3) * t455;
t744 = -t458 * t230 + t404 * t575 + (-t458 * t579 - t459 * t577) * pkin(9) + t832 * t454;
t567 = pkin(9) * t579;
t831 = t832 * t458 + (t230 + t567) * t454;
t451 = qJ(4) + qJ(5);
t447 = sin(t451);
t448 = cos(t451);
t823 = mrSges(6,2) - mrSges(7,3);
t815 = m(7) * qJ(6) - t823;
t824 = mrSges(6,1) + mrSges(7,1);
t817 = m(7) * pkin(5) + t824;
t830 = t447 * t815 + t448 * t817;
t583 = pkin(8) * t604 + t438;
t366 = t583 * qJD(1);
t430 = t538 + qJD(2);
t309 = t430 * pkin(9) + t366;
t323 = (-pkin(2) * t460 - pkin(9) * t456 - pkin(1)) * t581;
t200 = -t455 * t309 + t459 * t323;
t331 = t430 * t459 - t455 * t555;
t526 = t459 * t555;
t332 = t430 * t455 + t526;
t244 = pkin(3) * t332 - pkin(10) * t331;
t133 = -t200 * t454 + t458 * t244;
t461 = -pkin(11) - pkin(10);
t558 = qJD(4) * t461;
t617 = t331 * t458;
t829 = -pkin(4) * t332 + pkin(11) * t617 + t458 * t558 - t133;
t134 = t458 * t200 + t454 * t244;
t618 = t331 * t454;
t828 = -pkin(11) * t618 - t454 * t558 + t134;
t601 = t454 * t456;
t322 = (t458 * t595 + t601) * t581;
t596 = t458 * t459;
t440 = pkin(9) * t596;
t554 = t460 * t581;
t527 = t455 * t554;
t827 = -pkin(4) * t527 + pkin(11) * t322 + (pkin(4) * t455 - pkin(11) * t596) * qJD(3) + (-t440 + (pkin(11) * t455 - t404) * t454) * qJD(4) + t831;
t826 = pkin(11) * t742 - t744;
t410 = qJD(3) - t554;
t172 = -t410 * pkin(3) - t200;
t493 = t332 * t454 - t410 * t458;
t135 = pkin(4) * t493 + t172;
t325 = qJD(4) - t331;
t317 = qJD(5) + t325;
t457 = cos(qJ(5));
t453 = sin(qJ(5));
t308 = -t430 * pkin(2) - t363;
t167 = -t331 * pkin(3) - t332 * pkin(10) + t308;
t201 = t459 * t309 + t455 * t323;
t173 = pkin(10) * t410 + t201;
t111 = t454 * t167 + t458 * t173;
t88 = -pkin(11) * t493 + t111;
t626 = t453 * t88;
t110 = t458 * t167 - t173 * t454;
t260 = t332 * t458 + t410 * t454;
t87 = -pkin(11) * t260 + t110;
t74 = pkin(4) * t325 + t87;
t31 = t457 * t74 - t626;
t747 = qJD(6) - t31;
t27 = -pkin(5) * t317 + t747;
t150 = t453 * t260 + t457 * t493;
t147 = Ifges(6,4) * t150;
t466 = t457 * t260 - t453 * t493;
t633 = Ifges(7,5) * t150;
t776 = Ifges(7,4) + Ifges(6,5);
t778 = Ifges(6,1) + Ifges(7,1);
t767 = t317 * t776 + t466 * t778 - t147 + t633;
t799 = mrSges(7,2) * t27 - mrSges(6,3) * t31 + t767 / 0.2e1;
t825 = t135 * mrSges(6,2) + t799;
t777 = -Ifges(6,4) + Ifges(7,5);
t775 = -Ifges(6,6) + Ifges(7,6);
t774 = Ifges(6,3) + Ifges(7,2);
t394 = t453 * t458 + t454 * t457;
t369 = t394 * t455;
t393 = t453 * t454 - t457 * t458;
t731 = qJD(4) + qJD(5);
t185 = -t369 * t731 - t393 * t578;
t218 = t321 * t453 + t322 * t457;
t822 = t185 - t218;
t576 = qJD(4) * t455;
t474 = -t454 * t576 + t458 * t578;
t574 = qJD(5) * t453;
t600 = t455 * t458;
t602 = t454 * t455;
t186 = -t574 * t602 + (t600 * t731 + t551) * t457 + t474 * t453;
t217 = -t457 * t321 + t322 * t453;
t821 = t186 - t217;
t254 = -t455 * t363 + t364 * t459;
t229 = -pkin(3) * t555 - t254;
t446 = pkin(9) * t578;
t820 = t446 - t229;
t739 = t527 - t579;
t657 = cos(qJ(1));
t519 = t623 * t657;
t656 = sin(qJ(1));
t380 = t456 * t519 + t460 * t656;
t557 = t452 * t657;
t296 = t380 * t459 - t455 * t557;
t379 = t456 * t656 - t460 * t519;
t804 = -t296 * t454 + t379 * t458;
t627 = t332 * Ifges(4,4);
t756 = t410 * Ifges(4,6);
t193 = t331 * Ifges(4,2) + t627 + t756;
t630 = t201 * mrSges(4,3);
t669 = t317 / 0.2e1;
t684 = t466 / 0.2e1;
t687 = t150 / 0.2e1;
t688 = -t150 / 0.2e1;
t755 = t493 * Ifges(5,6);
t726 = Ifges(5,5) * t260 + Ifges(5,3) * t325 + t150 * t775 + t317 * t774 + t466 * t776 - t755;
t819 = t774 * t669 + t776 * t684 - t630 - t193 / 0.2e1 + Ifges(7,6) * t687 + Ifges(6,6) * t688 + t726 / 0.2e1;
t62 = t150 * pkin(5) - qJ(6) * t466 + t135;
t781 = t62 * mrSges(7,3);
t818 = -t781 + t825;
t146 = Ifges(7,5) * t466;
t78 = Ifges(7,6) * t317 + Ifges(7,3) * t150 + t146;
t634 = Ifges(6,4) * t466;
t81 = -Ifges(6,2) * t150 + Ifges(6,6) * t317 + t634;
t816 = -t78 / 0.2e1 - t62 * mrSges(7,1) + t81 / 0.2e1;
t790 = -m(7) - m(6);
t812 = -m(4) + t790;
t392 = t458 * t404;
t280 = -pkin(11) * t600 + t392 + (-pkin(9) * t454 - pkin(4)) * t459;
t342 = t454 * t404 + t440;
t303 = -pkin(11) * t602 + t342;
t573 = qJD(5) * t457;
t769 = t280 * t573 - t303 * t574 + t453 * t827 - t457 * t826;
t743 = t453 * t280 + t457 * t303;
t768 = -qJD(5) * t743 + t453 * t826 + t457 * t827;
t513 = -mrSges(5,1) * t458 + mrSges(5,2) * t454;
t480 = m(5) * pkin(3) - t513;
t724 = t480 + t830;
t811 = mrSges(4,1) + t724;
t566 = m(5) * pkin(10) + mrSges(5,3);
t779 = mrSges(6,3) + mrSges(7,2);
t728 = -t566 - t779;
t810 = mrSges(4,2) + t728;
t38 = t457 * t87 - t626;
t809 = pkin(4) * t573 - t38;
t413 = t461 * t454;
t414 = t461 * t458;
t312 = t413 * t453 - t414 * t457;
t749 = -qJD(5) * t312 + t828 * t453 + t457 * t829;
t492 = t457 * t413 + t414 * t453;
t748 = qJD(5) * t492 + t453 * t829 - t828 * t457;
t642 = mrSges(6,3) * t466;
t131 = mrSges(6,1) * t317 - t642;
t132 = -mrSges(7,1) * t317 + mrSges(7,2) * t466;
t808 = t131 - t132;
t740 = pkin(4) * t742 + t820;
t515 = t459 * mrSges(4,1) - t455 * mrSges(4,2);
t807 = -t455 * t566 - t459 * t480 - t515;
t571 = qJDD(1) * t452;
t806 = pkin(8) * t571 + qJD(2) * t525;
t533 = t623 * qJDD(1);
t572 = qJD(1) * qJD(2);
t805 = -pkin(8) * t452 * t572 + pkin(1) * t533;
t738 = -t201 + (t577 - t618) * pkin(4);
t518 = t623 * t656;
t382 = -t456 * t518 + t460 * t657;
t556 = t452 * t656;
t300 = t382 * t459 + t455 * t556;
t381 = t456 * t657 + t460 * t518;
t231 = -t300 * t454 + t381 * t458;
t624 = t457 * t88;
t32 = t453 * t74 + t624;
t28 = qJ(6) * t317 + t32;
t803 = mrSges(7,2) * t28 + mrSges(6,3) * t32;
t372 = (qJDD(1) * t456 + t460 * t572) * t452;
t429 = t533 + qJDD(2);
t213 = qJD(3) * t331 + t372 * t459 + t429 * t455;
t371 = (-qJDD(1) * t460 + t456 * t572) * t452;
t359 = qJDD(3) + t371;
t481 = t493 * qJD(4);
t124 = t213 * t458 + t359 * t454 - t481;
t125 = -qJD(4) * t260 - t213 * t454 + t359 * t458;
t46 = -qJD(5) * t150 + t457 * t124 + t453 * t125;
t47 = qJD(5) * t466 + t453 * t124 - t457 * t125;
t261 = t456 * t805 + t460 * t806;
t239 = pkin(9) * t429 + t261;
t565 = pkin(1) * t571;
t247 = pkin(2) * t371 - pkin(9) * t372 - t565;
t108 = -t455 * t239 + t247 * t459 - t309 * t578 - t323 * t579;
t93 = -pkin(3) * t359 - t108;
t59 = -pkin(4) * t125 + t93;
t11 = pkin(5) * t47 - qJ(6) * t46 - qJD(6) * t466 + t59;
t214 = -qJD(3) * t526 - t455 * t372 + t429 * t459 - t430 * t579;
t199 = qJDD(4) - t214;
t188 = qJDD(5) + t199;
t262 = -t456 * t806 + t460 * t805;
t240 = -t429 * pkin(2) - t262;
t105 = -t214 * pkin(3) - t213 * pkin(10) + t240;
t107 = t459 * t239 + t455 * t247 - t309 * t579 + t323 * t578;
t92 = pkin(10) * t359 + t107;
t26 = -qJD(4) * t111 + t458 * t105 - t454 * t92;
t22 = pkin(4) * t199 - pkin(11) * t124 + t26;
t25 = t454 * t105 + t167 * t575 - t173 * t577 + t458 * t92;
t24 = pkin(11) * t125 + t25;
t6 = -qJD(5) * t32 + t22 * t457 - t24 * t453;
t3 = -pkin(5) * t188 + qJDD(6) - t6;
t681 = t188 / 0.2e1;
t700 = t47 / 0.2e1;
t701 = -t47 / 0.2e1;
t702 = t46 / 0.2e1;
t772 = t188 * t776 + t46 * t778 + t47 * t777;
t801 = mrSges(6,2) * t59 + mrSges(7,2) * t3 - mrSges(6,3) * t6 - mrSges(7,3) * t11 + Ifges(6,4) * t701 + Ifges(7,5) * t700 + t776 * t681 + t778 * t702 + t772 / 0.2e1;
t759 = t135 * mrSges(6,1);
t800 = t759 - t803;
t95 = pkin(5) * t466 + qJ(6) * t150;
t219 = t296 * t447 - t379 * t448;
t220 = t296 * t448 + t379 * t447;
t614 = t379 * t454;
t798 = t296 * t458 + t614;
t712 = -Ifges(6,2) * t688 + Ifges(7,3) * t687 + t669 * t775 + t684 * t777 - t816;
t797 = t712 + t800;
t796 = -t31 * mrSges(6,1) + t27 * mrSges(7,1) + t32 * mrSges(6,2) - t28 * mrSges(7,3);
t670 = -t317 / 0.2e1;
t685 = -t466 / 0.2e1;
t795 = -Ifges(6,4) * t687 - Ifges(7,5) * t688 - t776 * t670 - t778 * t685 + t818;
t794 = Ifges(6,4) * t688 + Ifges(7,5) * t687 + t776 * t669 + t778 * t684;
t512 = t454 * mrSges(5,1) + t458 * mrSges(5,2);
t780 = mrSges(3,2) - mrSges(4,3);
t793 = pkin(4) * t454 * t790 - m(5) * pkin(9) - t447 * t817 + t448 * t815 - t512 + t780;
t445 = pkin(4) * t458 + pkin(3);
t791 = mrSges(3,1) - t807 + (t461 * t790 + t779) * t455 + (-t445 * t790 + t830) * t459;
t691 = t124 / 0.2e1;
t690 = t125 / 0.2e1;
t679 = t199 / 0.2e1;
t678 = t213 / 0.2e1;
t677 = t214 / 0.2e1;
t662 = t359 / 0.2e1;
t789 = -t493 / 0.2e1;
t788 = t493 / 0.2e1;
t773 = t188 * t774 + t46 * t776 + t47 * t775;
t771 = -qJ(6) * t739 - qJD(6) * t459 + t769;
t770 = pkin(5) * t739 - t768;
t370 = t393 * t455;
t766 = pkin(5) * t821 - qJ(6) * t822 + qJD(6) * t370 + t740;
t757 = t410 * Ifges(4,5);
t215 = t394 * t331;
t216 = t393 * t331;
t287 = t731 * t393;
t288 = t731 * t394;
t754 = -qJD(6) * t394 + t738 + (-t216 + t287) * qJ(6) + (-t215 + t288) * pkin(5);
t753 = qJD(6) + t809;
t752 = -qJ(6) * t332 + t748;
t751 = pkin(5) * t332 - t749;
t605 = t452 * t456;
t387 = -pkin(8) * t605 + t460 * t559;
t354 = -pkin(2) * t623 - t387;
t377 = t455 * t605 - t459 * t623;
t378 = t455 * t623 + t459 * t605;
t226 = t377 * pkin(3) - t378 * pkin(10) + t354;
t355 = pkin(9) * t623 + t583;
t584 = pkin(2) * t604 + pkin(9) * t605;
t651 = pkin(1) * t452;
t356 = -t584 - t651;
t246 = t459 * t355 + t455 * t356;
t228 = -pkin(10) * t604 + t246;
t127 = t458 * t226 - t228 * t454;
t487 = -t378 * t458 + t454 * t604;
t104 = pkin(4) * t377 + pkin(11) * t487 + t127;
t128 = t454 * t226 + t458 * t228;
t293 = -t378 * t454 - t458 * t604;
t112 = pkin(11) * t293 + t128;
t746 = t453 * t104 + t457 * t112;
t745 = -qJD(4) * t342 + t831;
t741 = t322 - t474;
t277 = t378 * t447 + t448 * t604;
t562 = t447 * t604;
t278 = t378 * t448 - t562;
t737 = t277 * t824 + t823 * t278;
t223 = t300 * t447 - t381 * t448;
t224 = t300 * t448 + t381 * t447;
t736 = t223 * t824 + t823 * t224;
t735 = t219 * t824 + t823 * t220;
t258 = Ifges(5,4) * t493;
t138 = t260 * Ifges(5,1) + t325 * Ifges(5,5) - t258;
t324 = Ifges(4,4) * t331;
t194 = t332 * Ifges(4,1) + t324 + t757;
t734 = t458 * t138 + t194;
t733 = t107 * t459 - t108 * t455;
t732 = t25 * t458 - t26 * t454;
t641 = Ifges(3,4) * t456;
t706 = t452 ^ 2;
t730 = (t456 * (Ifges(3,1) * t460 - t641) / 0.2e1 - pkin(1) * (mrSges(3,1) * t456 + mrSges(3,2) * t460)) * t706;
t54 = Ifges(5,5) * t124 + Ifges(5,6) * t125 + Ifges(5,3) * t199;
t729 = t54 + t773;
t580 = qJD(2) * t460;
t552 = t452 * t580;
t292 = -qJD(3) * t377 + t459 * t552;
t553 = qJD(2) * t605;
t166 = qJD(4) * t293 + t292 * t458 + t454 * t553;
t291 = qJD(3) * t378 + t455 * t552;
t365 = qJD(2) * t488;
t367 = t387 * qJD(2);
t144 = -t355 * t579 + t356 * t578 + t455 * t365 + t459 * t367;
t140 = pkin(10) * t553 + t144;
t368 = t583 * qJD(2);
t160 = t291 * pkin(3) - t292 * pkin(10) + t368;
t61 = -qJD(4) * t128 - t140 * t454 + t458 * t160;
t30 = pkin(4) * t291 - pkin(11) * t166 + t61;
t165 = qJD(4) * t487 - t292 * t454 + t458 * t553;
t60 = t458 * t140 + t454 * t160 + t226 * t575 - t228 * t577;
t40 = pkin(11) * t165 + t60;
t10 = -qJD(5) * t746 + t30 * t457 - t40 * t453;
t532 = mrSges(3,3) * t555;
t719 = -m(4) * t308 + mrSges(3,1) * t430 + mrSges(4,1) * t331 - mrSges(4,2) * t332 - t532;
t718 = t26 * mrSges(5,1) - t25 * mrSges(5,2);
t5 = t453 * t22 + t457 * t24 + t74 * t573 - t574 * t88;
t2 = qJ(6) * t188 + qJD(6) * t317 + t5;
t715 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3);
t714 = -t781 + t794;
t713 = -Ifges(6,2) * t687 + Ifges(7,3) * t688 + t670 * t775 + t685 * t777 + t816;
t709 = t713 + t803;
t708 = t111 * mrSges(5,2) - t110 * mrSges(5,1) - t308 * mrSges(4,1) + t756 / 0.2e1 + t755 / 0.2e1 + t796;
t707 = mrSges(6,1) * t59 + mrSges(7,1) * t11 - mrSges(7,2) * t2 - mrSges(6,3) * t5 + 0.2e1 * Ifges(7,3) * t700 - t46 * Ifges(6,4) / 0.2e1 - t188 * Ifges(6,6) / 0.2e1 + (t777 + Ifges(7,5)) * t702 + (t775 + Ifges(7,6)) * t681 + (-t701 + t700) * Ifges(6,2);
t55 = t124 * Ifges(5,4) + t125 * Ifges(5,2) + t199 * Ifges(5,6);
t699 = t55 / 0.2e1;
t698 = Ifges(5,1) * t691 + Ifges(5,4) * t690 + Ifges(5,5) * t679;
t692 = Ifges(4,1) * t678 + Ifges(4,4) * t677 + Ifges(4,5) * t662;
t689 = t138 / 0.2e1;
t674 = -t260 / 0.2e1;
t673 = t260 / 0.2e1;
t668 = -t325 / 0.2e1;
t667 = t325 / 0.2e1;
t665 = t331 / 0.2e1;
t663 = t332 / 0.2e1;
t650 = pkin(4) * t260;
t477 = pkin(4) * t293;
t649 = pkin(4) * t453;
t647 = pkin(4) * t457;
t646 = pkin(9) * t459;
t449 = t455 * pkin(9);
t643 = mrSges(6,3) * t150;
t640 = Ifges(3,4) * t460;
t639 = Ifges(4,4) * t455;
t638 = Ifges(4,4) * t459;
t637 = Ifges(5,4) * t260;
t636 = Ifges(5,4) * t454;
t635 = Ifges(5,4) * t458;
t632 = t111 * mrSges(5,3);
t631 = t200 * mrSges(4,3);
t625 = t455 * t93;
t611 = t381 * t454;
t137 = -t493 * Ifges(5,2) + Ifges(5,6) * t325 + t637;
t603 = t454 * t137;
t599 = t455 * t460;
t399 = pkin(4) * t602 + t449;
t582 = t657 * pkin(1) + pkin(8) * t556;
t563 = t452 * t599;
t561 = Ifges(4,5) * t213 + Ifges(4,6) * t214 + Ifges(4,3) * t359;
t560 = Ifges(3,5) * t372 - Ifges(3,6) * t371 + Ifges(3,3) * t429;
t547 = t605 / 0.2e1;
t546 = -t603 / 0.2e1;
t541 = -t576 / 0.2e1;
t35 = -t188 * mrSges(7,1) + t46 * mrSges(7,2);
t537 = -t219 * pkin(5) + qJ(6) * t220;
t536 = -t223 * pkin(5) + qJ(6) * t224;
t535 = -t277 * pkin(5) + qJ(6) * t278;
t245 = -t455 * t355 + t356 * t459;
t534 = -t380 * t455 - t459 * t557;
t531 = mrSges(3,3) * t554;
t523 = t804 * pkin(4);
t522 = t231 * pkin(4);
t520 = -pkin(1) * t656 + pkin(8) * t557;
t227 = pkin(3) * t604 - t245;
t516 = mrSges(4,1) * t377 + mrSges(4,2) * t378;
t514 = mrSges(5,1) * t293 + mrSges(5,2) * t487;
t509 = Ifges(4,1) * t459 - t639;
t508 = Ifges(5,1) * t458 - t636;
t507 = Ifges(5,1) * t454 + t635;
t506 = -Ifges(4,2) * t455 + t638;
t505 = -Ifges(5,2) * t454 + t635;
t504 = Ifges(5,2) * t458 + t636;
t503 = Ifges(4,5) * t459 - Ifges(4,6) * t455;
t502 = Ifges(5,5) * t458 - Ifges(5,6) * t454;
t501 = Ifges(5,5) * t454 + Ifges(5,6) * t458;
t52 = t104 * t457 - t112 * t453;
t174 = t280 * t457 - t303 * t453;
t494 = t457 * t293 + t453 * t487;
t180 = t293 * t453 - t457 * t487;
t490 = t382 * pkin(2) + pkin(9) * t381 + t582;
t145 = -t355 * t578 - t356 * t579 + t365 * t459 - t455 * t367;
t9 = t104 * t573 - t112 * t574 + t453 * t30 + t457 * t40;
t486 = t172 * t512;
t484 = (t460 * Ifges(3,2) + t641) * t452;
t483 = t493 * mrSges(5,3);
t153 = t227 - t477;
t476 = t430 * t452 * (Ifges(3,5) * t460 - Ifges(3,6) * t456);
t472 = -t380 * pkin(2) - t379 * pkin(9) + t520;
t141 = -pkin(3) * t553 - t145;
t464 = t715 + t773;
t94 = -pkin(4) * t165 + t141;
t444 = -pkin(5) - t647;
t442 = qJ(6) + t649;
t423 = Ifges(3,4) * t554;
t383 = (-mrSges(3,1) * t460 + mrSges(3,2) * t456) * t452;
t375 = t381 * pkin(2);
t373 = t379 * pkin(2);
t362 = -mrSges(3,2) * t430 + t531;
t341 = -t454 * t646 + t392;
t306 = Ifges(3,1) * t555 + t430 * Ifges(3,5) + t423;
t305 = t430 * Ifges(3,6) + qJD(1) * t484;
t299 = t382 * t455 - t459 * t556;
t271 = pkin(5) * t393 - qJ(6) * t394 - t445;
t265 = mrSges(4,1) * t410 - mrSges(4,3) * t332;
t264 = -mrSges(4,2) * t410 + mrSges(4,3) * t331;
t233 = pkin(5) * t369 + qJ(6) * t370 + t399;
t232 = t300 * t458 + t611;
t192 = t332 * Ifges(4,5) + t331 * Ifges(4,6) + t410 * Ifges(4,3);
t178 = mrSges(5,1) * t325 - mrSges(5,3) * t260;
t177 = -t325 * mrSges(5,2) - t483;
t170 = pkin(5) * t459 - t174;
t169 = -qJ(6) * t459 + t743;
t159 = -mrSges(4,2) * t359 + mrSges(4,3) * t214;
t158 = mrSges(4,1) * t359 - mrSges(4,3) * t213;
t157 = mrSges(5,1) * t493 + t260 * mrSges(5,2);
t130 = -mrSges(6,2) * t317 - t643;
t129 = -mrSges(7,2) * t150 + mrSges(7,3) * t317;
t126 = -mrSges(4,1) * t214 + mrSges(4,2) * t213;
t116 = t213 * Ifges(4,4) + t214 * Ifges(4,2) + t359 * Ifges(4,6);
t97 = mrSges(6,1) * t150 + mrSges(6,2) * t466;
t96 = mrSges(7,1) * t150 - mrSges(7,3) * t466;
t90 = -mrSges(5,2) * t199 + mrSges(5,3) * t125;
t89 = mrSges(5,1) * t199 - mrSges(5,3) * t124;
t75 = t650 + t95;
t72 = -pkin(5) * t494 - qJ(6) * t180 + t153;
t69 = qJD(5) * t180 - t457 * t165 + t166 * t453;
t68 = qJD(5) * t494 + t165 * t453 + t166 * t457;
t67 = -mrSges(5,1) * t125 + mrSges(5,2) * t124;
t49 = -pkin(5) * t377 - t52;
t48 = qJ(6) * t377 + t746;
t37 = t453 * t87 + t624;
t36 = -mrSges(6,2) * t188 - mrSges(6,3) * t47;
t34 = mrSges(6,1) * t188 - mrSges(6,3) * t46;
t33 = -mrSges(7,2) * t47 + mrSges(7,3) * t188;
t21 = mrSges(6,1) * t47 + mrSges(6,2) * t46;
t20 = mrSges(7,1) * t47 - mrSges(7,3) * t46;
t18 = pkin(5) * t69 - qJ(6) * t68 - qJD(6) * t180 + t94;
t8 = -pkin(5) * t291 - t10;
t7 = qJ(6) * t291 + qJD(6) * t377 + t9;
t1 = [(Ifges(5,1) * t166 + Ifges(5,4) * t165) * t673 + (t476 / 0.2e1 + t192 * t547) * qJD(2) + (t261 * t604 - t262 * t605 - t363 * t552 - t366 * t553 - t371 * t583 - t372 * t387) * mrSges(3,3) + (t262 * t623 - t371 * t651 + t387 * t429) * mrSges(3,1) + (Ifges(3,4) * t372 - Ifges(3,2) * t371 + Ifges(3,6) * t429) * t604 / 0.2e1 + (Ifges(3,1) * t372 - Ifges(3,4) * t371 + Ifges(3,5) * t429) * t547 + (-t261 * t623 - t372 * t651 - t429 * t583) * mrSges(3,2) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t706 + t261 * t583 + t262 * t387 + t366 * t367) + (t706 * qJD(1) * (-Ifges(3,2) * t456 + t640) + t452 * t306) * t580 / 0.2e1 + (Ifges(4,1) * t378 - Ifges(4,5) * t604) * t678 + (Ifges(4,1) * t292 + Ifges(4,5) * t553) * t663 + (t794 + t818) * t68 + (-Ifges(4,4) * t663 + Ifges(5,5) * t673 - Ifges(4,2) * t665 + Ifges(5,3) * t667 - t708 + t819) * t291 + t797 * t69 - t93 * t514 + t240 * t516 + (Ifges(5,5) * t166 + Ifges(5,6) * t165) * t667 + (Ifges(4,5) * t378 - Ifges(4,3) * t604) * t662 + t410 * (Ifges(4,5) * t292 + Ifges(4,3) * t553) / 0.2e1 - t383 * t565 + t801 * t180 + t623 * t560 / 0.2e1 + t429 * (Ifges(3,3) * t623 + (Ifges(3,5) * t456 + Ifges(3,6) * t460) * t452) / 0.2e1 - t561 * t604 / 0.2e1 + (-Ifges(4,2) * t677 - Ifges(4,4) * t678 + Ifges(5,3) * t679 - t107 * mrSges(4,3) - Ifges(4,6) * t662 - t116 / 0.2e1 + Ifges(5,6) * t690 + Ifges(5,5) * t691 + Ifges(7,6) * t700 + Ifges(6,6) * t701 + t776 * t702 + t774 * t681 + t715 + t718 + t729 / 0.2e1) * t377 + t7 * t129 + t9 * t130 + t10 * t131 + t8 * t132 + t372 * (Ifges(3,5) * t623 + (t456 * Ifges(3,1) + t640) * t452) / 0.2e1 + (Ifges(4,4) * t378 - Ifges(4,6) * t604) * t677 + (Ifges(4,4) * t292 + Ifges(4,6) * t553) * t665 - t371 * (Ifges(3,6) * t623 + t484) / 0.2e1 + m(5) * (t110 * t61 + t111 * t60 + t127 * t26 + t128 * t25 + t141 * t172 + t227 * t93) + m(7) * (t11 * t72 + t18 * t62 + t2 * t48 + t27 * t8 + t28 * t7 + t3 * t49) + t108 * (-mrSges(4,1) * t604 - mrSges(4,3) * t378) + m(6) * (t10 * t31 + t135 * t94 + t153 * t59 + t32 * t9 + t5 * t746 + t52 * t6) + t746 * t36 - t707 * t494 + (-Ifges(5,5) * t487 + Ifges(5,6) * t293) * t679 + (-t110 * t166 + t111 * t165 + t25 * t293 + t26 * t487) * mrSges(5,3) + (-Ifges(5,4) * t487 + Ifges(5,2) * t293) * t690 + (-Ifges(5,1) * t487 + Ifges(5,4) * t293) * t691 - t487 * t698 + (-m(3) * t363 - t719) * t368 + (-m(5) * (-pkin(3) * t296 + t472) + t798 * mrSges(5,1) + t804 * mrSges(5,2) + mrSges(2,1) * t656 + mrSges(2,2) * t657 - m(3) * t520 + t380 * mrSges(3,1) - mrSges(3,3) * t557 - m(4) * t472 + t296 * mrSges(4,1) + t790 * (-pkin(4) * t614 - t296 * t445 - t461 * t534 + t472) - t780 * t379 + t817 * t220 + t815 * t219 + t810 * t534) * g(1) + (-mrSges(2,1) * t657 + mrSges(2,2) * t656 - m(3) * t582 - t382 * mrSges(3,1) - mrSges(3,3) * t556 - m(5) * (pkin(3) * t300 + t490) - t232 * mrSges(5,1) - t231 * mrSges(5,2) - m(4) * t490 - t300 * mrSges(4,1) + t790 * (pkin(4) * t611 - t299 * t461 + t300 * t445 + t490) + t780 * t381 - t817 * t224 - t815 * t223 + t810 * t299) * g(2) + t18 * t96 + t94 * t97 + t72 * t20 + m(4) * (t107 * t246 + t108 * t245 + t144 * t201 + t145 * t200 + t240 * t354) + Ifges(2,3) * qJDD(1) + t52 * t34 + t48 * t33 + t49 * t35 + t367 * t362 + t354 * t126 + t292 * t194 / 0.2e1 + t166 * t689 + t378 * t692 + t153 * t21 + t141 * t157 + t165 * t137 / 0.2e1 + t172 * (-mrSges(5,1) * t165 + mrSges(5,2) * t166) + t127 * t89 + t128 * t90 + t60 * t177 + t61 * t178 + t730 * t572 + t293 * t699 + t227 * t67 + t245 * t158 + t246 * t159 + (t107 * t604 - t201 * t553 + t292 * t308) * mrSges(4,2) + t144 * t264 + t145 * t265 + (Ifges(5,4) * t166 + Ifges(5,2) * t165) * t789 - t305 * t553 / 0.2e1 + t200 * (mrSges(4,1) * t553 - mrSges(4,3) * t292); t6 * (-mrSges(6,1) * t459 + mrSges(6,3) * t370) + t3 * (mrSges(7,1) * t459 - mrSges(7,2) * t370) + (t32 * t527 - t370 * t59 + t459 * t5) * mrSges(6,2) + (-Ifges(6,4) * t370 - Ifges(6,6) * t459) * t701 + (t11 * t370 - t2 * t459 + t218 * t62 - t28 * t527) * mrSges(7,3) + (-Ifges(7,5) * t370 - Ifges(7,6) * t459) * t700 + (t331 * t506 + t332 * t509 + t410 * t503) * qJD(3) / 0.2e1 + (t712 - t803) * t186 + (-t796 + t819) * t579 + t707 * t369 + (t531 - t362) * t363 + (mrSges(6,1) * t821 + mrSges(6,2) * t822) * t135 + t410 * t308 * (mrSges(4,1) * t455 + mrSges(4,2) * t459) + (Ifges(6,4) * t218 + Ifges(6,6) * t527) * t687 + (-t158 + t67) * t449 + t26 * (-mrSges(5,1) * t459 - mrSges(5,3) * t600) + (Ifges(4,2) * t459 + t639) * t677 + (Ifges(4,1) * t455 + t638) * t678 + (-Ifges(5,3) * t459 + t455 * t502) * t679 - t240 * t515 + t820 * t157 + t193 * t527 / 0.2e1 + t709 * t217 + t560 + (-t200 * (mrSges(4,1) * t456 - mrSges(4,3) * t595) - t201 * (-mrSges(4,2) * t456 - mrSges(4,3) * t599)) * t581 - t55 * t602 / 0.2e1 + t25 * (mrSges(5,2) * t459 - mrSges(5,3) * t602) - t27 * (-mrSges(7,1) * t527 + mrSges(7,2) * t218) - t31 * (mrSges(6,1) * t527 - mrSges(6,3) * t218) + (-t476 / 0.2e1 + t305 * t547 - t730 * qJD(1)) * qJD(1) + (-t567 - t255) * t264 + (t714 + t799) * t185 - ((-Ifges(3,2) * t555 + t459 * t194 + t455 * t726 + t306 + t423) * t460 + t410 * (Ifges(4,3) * t456 + t460 * t503) + t332 * (Ifges(4,5) * t456 + t460 * t509) + t331 * (Ifges(4,6) * t456 + t460 * t506) + t456 * t192) * t581 / 0.2e1 + (t135 * t740 + t174 * t6 + t31 * t768 + t32 * t769 + t399 * t59 + t5 * t743) * m(6) + t743 * t36 + (Ifges(7,5) * t218 + Ifges(7,6) * t527) * t688 + (-t501 * t576 + (Ifges(5,3) * t455 + t459 * t502) * qJD(3)) * t667 + (Ifges(5,5) * t322 + Ifges(5,6) * t321 + Ifges(5,3) * t527) * t668 + (-t507 * t576 + (Ifges(5,5) * t455 + t459 * t508) * qJD(3)) * t673 + (Ifges(5,1) * t322 + Ifges(5,4) * t321 + Ifges(5,5) * t527) * t674 + (-t631 + t546) * t578 + t159 * t646 + (Ifges(4,5) * t455 + Ifges(4,6) * t459) * t662 + t512 * t625 + (t532 + t719) * t366 + (t383 + (-m(4) - m(5)) * t584 - t779 * t563 + t790 * (-t461 * t563 + t584) - t815 * (-t448 * t605 + t459 * t562) + (t790 * (pkin(4) * t601 + t445 * t595) + t807 * t460 + (-mrSges(4,3) - t512) * t456 - t817 * (t447 * t456 + t448 * t595)) * t452) * g(3) + (-t322 / 0.2e1 + t454 * t541) * t138 + (-t321 / 0.2e1 + t458 * t541) * t137 + t459 * t116 / 0.2e1 + t399 * t21 + t341 * t89 + t342 * t90 + (-t370 * t778 - t459 * t776) * t702 + (t218 * t776 + t527 * t774) * t670 + (-t370 * t776 - t459 * t774) * t681 + (t218 * t778 + t527 * t776) * t685 + t766 * t96 - t767 * t218 / 0.2e1 + t768 * t131 + t769 * t130 + t770 * t132 + t771 * t129 + (t11 * t233 + t169 * t2 + t170 * t3 + t27 * t770 + t28 * t771 + t62 * t766) * m(7) - t772 * t370 / 0.2e1 + (-Ifges(5,6) * t459 + t455 * t505) * t690 + (-Ifges(5,5) * t459 + t455 * t508) * t691 + t455 * t692 + t169 * t33 + t170 * t35 + t174 * t34 - pkin(2) * t126 + (-t446 - t254) * t265 - t729 * t459 / 0.2e1 + (-t200 * t254 - t201 * t255 - pkin(2) * t240 + ((-t200 * t459 - t201 * t455) * qJD(3) + t733) * pkin(9)) * m(4) + t733 * mrSges(4,3) + t734 * t578 / 0.2e1 + t740 * t97 + (-mrSges(5,1) * t739 + mrSges(5,3) * t741) * t110 + (mrSges(5,1) * t742 - mrSges(5,2) * t741) * t172 + (mrSges(5,2) * t739 - mrSges(5,3) * t742) * t111 + t744 * t177 + t745 * t178 + (t25 * t342 + t26 * t341 + (t172 * t578 + t625) * pkin(9) - t172 * t229 + t744 * t111 + t745 * t110) * m(5) + t600 * t698 + (m(5) * t375 + t812 * (t382 * pkin(9) - t375) + t793 * t382 + t791 * t381) * g(1) + (m(5) * t373 + t812 * (t380 * pkin(9) - t373) + t793 * t380 + t791 * t379) * g(2) + t233 * t20 - t261 * mrSges(3,2) + t262 * mrSges(3,1) + (Ifges(5,4) * t322 + Ifges(5,2) * t321 + Ifges(5,6) * t527) * t788 + (-t504 * t576 + (Ifges(5,6) * t455 + t459 * t505) * qJD(3)) * t789; (-pkin(3) * t93 - t110 * t133 - t111 * t134 - t172 * t201) * m(5) + (t260 * t508 + t325 * t502) * qJD(4) / 0.2e1 + t707 * t393 + (t790 * (-t299 * t445 - t300 * t461) + t810 * t300 + t811 * t299) * g(1) + (t790 * (-t296 * t461 + t445 * t534) + t810 * t296 - t811 * t534) * g(2) + t797 * t288 + (t546 + t486) * qJD(4) + t501 * t679 + t93 * t513 + t795 * t216 - t505 * t481 / 0.2e1 + t801 * t394 + t561 - (t714 + t825) * t287 - (t35 - t34) * t492 + (t135 * t738 + t31 * t749 + t312 * t5 + t32 * t748 - t445 * t59 + t492 * t6) * m(6) + (t11 * t271 + t2 * t312 + t27 * t751 + t28 * t752 - t3 * t492 + t62 * t754) * m(7) + (t36 + t33) * t312 + t603 * t665 - t577 * t632 + t193 * t663 - t107 * mrSges(4,2) + t108 * mrSges(4,1) - pkin(3) * t67 - t445 * t21 + (t516 + t790 * (-t377 * t445 - t378 * t461) + t728 * t378 + t724 * t377) * g(3) + (-t486 + t505 * t788 + t502 * t668 + t508 * t674 + t631 - t757 / 0.2e1 - t308 * mrSges(4,2)) * t331 + (Ifges(5,5) * t674 + Ifges(6,6) * t687 + Ifges(7,6) * t688 + Ifges(5,3) * t668 + t774 * t670 + t776 * t685 + t630 + t708) * t332 + t575 * t689 + t504 * t690 + t507 * t691 - t134 * t177 - t133 * t178 - (Ifges(4,1) * t331 - t627 + t726) * t332 / 0.2e1 + (-t178 * t575 - t177 * t577 + m(5) * ((-t110 * t458 - t111 * t454) * qJD(4) + t732) - t454 * t89 + t458 * t90) * pkin(10) + (t111 * t618 + (-t575 + t617) * t110 + t732) * mrSges(5,3) - (-Ifges(4,2) * t332 + t324 + t734) * t331 / 0.2e1 + t738 * t97 + t454 * t698 + t458 * t699 + t748 * t130 + t749 * t131 + t751 * t132 + t752 * t129 + t754 * t96 - t200 * t264 + t271 * t20 + (t709 - t759) * t215 + (-t157 + t265) * t201; (-t177 - t483) * t110 + (-t804 * mrSges(5,1) + t798 * mrSges(5,2) - m(7) * (t523 + t537) - m(6) * t523 + t735) * g(2) + (-t135 * t650 + t31 * t37 - t32 * t38 + (t453 * t5 + t457 * t6 + (-t31 * t453 + t32 * t457) * qJD(5)) * pkin(4)) * m(6) - t97 * t650 + t795 * t150 + (t713 - t800) * t466 + (t2 * t442 + t28 * t753 + t3 * t444 - t62 * t75) * m(7) + t54 + (-Ifges(5,5) * t493 - Ifges(5,6) * t260) * t668 + t137 * t673 + (-Ifges(5,1) * t493 - t637) * t674 + t34 * t647 + t36 * t649 + t260 * t632 - t75 * t96 + t442 * t33 + t444 * t35 + t493 * t689 - t172 * (t260 * mrSges(5,1) - mrSges(5,2) * t493) + (-m(7) * t27 + t808) * (-pkin(4) * t574 + t37) + t718 + t111 * t178 + t809 * t130 + (-m(7) * (t522 + t536) - m(6) * t522 - mrSges(5,1) * t231 + mrSges(5,2) * t232 + t736) * g(1) + (-t514 - m(7) * (t477 + t535) - m(6) * t477 + t737) * g(3) + t753 * t129 + t464 + (-Ifges(5,2) * t260 - t258) * t788; t81 * t684 + (t150 * t27 + t28 * t466) * mrSges(7,2) + qJD(6) * t129 + (t642 + t808) * t32 + (-t129 - t130 - t643) * t31 + t737 * g(3) + t736 * g(1) + t735 * g(2) - t95 * t96 - pkin(5) * t35 + qJ(6) * t33 + (Ifges(7,3) * t466 - t633) * t688 - t62 * (mrSges(7,1) * t466 + mrSges(7,3) * t150) - t135 * (mrSges(6,1) * t466 - mrSges(6,2) * t150) + t464 + (-t150 * t776 + t466 * t775) * t670 + (-pkin(5) * t3 - t536 * g(1) - t537 * g(2) - t535 * g(3) + qJ(6) * t2 - t27 * t32 + t28 * t747 - t62 * t95) * m(7) + (-Ifges(6,2) * t466 - t147 + t767) * t687 + (-t150 * t778 + t146 - t634 + t78) * t685; -t317 * t129 + t466 * t96 + (-g(1) * t223 - g(2) * t219 - g(3) * t277 - t28 * t317 + t466 * t62 + t3) * m(7) + t35;];
tau  = t1;
