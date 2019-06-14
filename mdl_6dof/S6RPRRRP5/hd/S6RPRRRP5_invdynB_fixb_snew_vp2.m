% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRP5
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 01:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:29:02
% EndTime: 2019-05-06 01:29:15
% DurationCPUTime: 12.22s
% Computational Cost: add. (170502->341), mult. (409052->414), div. (0->0), fcn. (318370->10), ass. (0->140)
t755 = Ifges(6,1) + Ifges(7,1);
t748 = Ifges(6,4) - Ifges(7,5);
t754 = -Ifges(6,5) - Ifges(7,4);
t753 = Ifges(6,2) + Ifges(7,3);
t746 = Ifges(6,6) - Ifges(7,6);
t752 = -Ifges(6,3) - Ifges(7,2);
t715 = qJD(1) ^ 2;
t751 = cos(qJ(5));
t707 = cos(pkin(10));
t750 = pkin(2) * t707;
t749 = -mrSges(6,3) - mrSges(7,2);
t706 = sin(pkin(10));
t745 = mrSges(3,2) * t706;
t704 = t707 ^ 2;
t744 = t704 * t715;
t711 = sin(qJ(1));
t714 = cos(qJ(1));
t692 = -t714 * g(1) - t711 * g(2);
t688 = -t715 * pkin(1) + qJDD(1) * qJ(2) + t692;
t736 = qJD(1) * qJD(2);
t733 = -t707 * g(3) - 0.2e1 * t706 * t736;
t662 = (-pkin(7) * qJDD(1) + t715 * t750 - t688) * t706 + t733;
t679 = -t706 * g(3) + (t688 + 0.2e1 * t736) * t707;
t735 = qJDD(1) * t707;
t663 = -pkin(2) * t744 + pkin(7) * t735 + t679;
t710 = sin(qJ(3));
t713 = cos(qJ(3));
t642 = t713 * t662 - t710 * t663;
t722 = t706 * t713 + t707 * t710;
t721 = -t706 * t710 + t707 * t713;
t686 = t721 * qJD(1);
t737 = t686 * qJD(3);
t677 = t722 * qJDD(1) + t737;
t687 = t722 * qJD(1);
t611 = (-t677 + t737) * pkin(8) + (t686 * t687 + qJDD(3)) * pkin(3) + t642;
t643 = t710 * t662 + t713 * t663;
t676 = -t687 * qJD(3) + t721 * qJDD(1);
t682 = qJD(3) * pkin(3) - t687 * pkin(8);
t685 = t686 ^ 2;
t618 = -t685 * pkin(3) + t676 * pkin(8) - qJD(3) * t682 + t643;
t709 = sin(qJ(4));
t712 = cos(qJ(4));
t609 = t709 * t611 + t712 * t618;
t670 = t709 * t686 + t712 * t687;
t637 = -t670 * qJD(4) + t712 * t676 - t709 * t677;
t669 = t712 * t686 - t709 * t687;
t653 = -t669 * mrSges(5,1) + t670 * mrSges(5,2);
t705 = qJD(3) + qJD(4);
t660 = t705 * mrSges(5,1) - t670 * mrSges(5,3);
t702 = qJDD(3) + qJDD(4);
t654 = -t669 * pkin(4) - t670 * pkin(9);
t701 = t705 ^ 2;
t605 = -t701 * pkin(4) + t702 * pkin(9) + t669 * t654 + t609;
t703 = t706 ^ 2;
t691 = t711 * g(1) - t714 * g(2);
t726 = qJDD(2) - t691;
t675 = (-pkin(1) - t750) * qJDD(1) + (-qJ(2) + (-t703 - t704) * pkin(7)) * t715 + t726;
t630 = -t676 * pkin(3) - t685 * pkin(8) + t687 * t682 + t675;
t638 = t669 * qJD(4) + t709 * t676 + t712 * t677;
t607 = (-t669 * t705 - t638) * pkin(9) + (t670 * t705 - t637) * pkin(4) + t630;
t708 = sin(qJ(5));
t602 = t751 * t605 + t708 * t607;
t656 = t751 * t670 + t708 * t705;
t614 = t656 * qJD(5) + t708 * t638 - t751 * t702;
t636 = qJDD(5) - t637;
t665 = qJD(5) - t669;
t646 = t665 * mrSges(6,1) - t656 * mrSges(6,3);
t655 = t708 * t670 - t751 * t705;
t639 = t655 * pkin(5) - t656 * qJ(6);
t664 = t665 ^ 2;
t598 = -t664 * pkin(5) + t636 * qJ(6) + 0.2e1 * qJD(6) * t665 - t655 * t639 + t602;
t647 = -t665 * mrSges(7,1) + t656 * mrSges(7,2);
t734 = m(7) * t598 + t636 * mrSges(7,3) + t665 * t647;
t640 = t655 * mrSges(7,1) - t656 * mrSges(7,3);
t739 = -t655 * mrSges(6,1) - t656 * mrSges(6,2) - t640;
t593 = m(6) * t602 - t636 * mrSges(6,2) + t749 * t614 - t665 * t646 + t739 * t655 + t734;
t601 = -t708 * t605 + t751 * t607;
t615 = -t655 * qJD(5) + t751 * t638 + t708 * t702;
t645 = -t665 * mrSges(6,2) - t655 * mrSges(6,3);
t599 = -t636 * pkin(5) - t664 * qJ(6) + t656 * t639 + qJDD(6) - t601;
t644 = -t655 * mrSges(7,2) + t665 * mrSges(7,3);
t727 = -m(7) * t599 + t636 * mrSges(7,1) + t665 * t644;
t595 = m(6) * t601 + t636 * mrSges(6,1) + t749 * t615 + t665 * t645 + t739 * t656 + t727;
t728 = t751 * t593 - t708 * t595;
t585 = m(5) * t609 - t702 * mrSges(5,2) + t637 * mrSges(5,3) + t669 * t653 - t705 * t660 + t728;
t608 = t712 * t611 - t709 * t618;
t659 = -t705 * mrSges(5,2) + t669 * mrSges(5,3);
t604 = -t702 * pkin(4) - t701 * pkin(9) + t670 * t654 - t608;
t600 = -0.2e1 * qJD(6) * t656 + (t655 * t665 - t615) * qJ(6) + (t656 * t665 + t614) * pkin(5) + t604;
t596 = m(7) * t600 + t614 * mrSges(7,1) - t615 * mrSges(7,3) + t655 * t644 - t656 * t647;
t717 = -m(6) * t604 - t614 * mrSges(6,1) - t615 * mrSges(6,2) - t655 * t645 - t656 * t646 - t596;
t590 = m(5) * t608 + t702 * mrSges(5,1) - t638 * mrSges(5,3) - t670 * t653 + t705 * t659 + t717;
t580 = t709 * t585 + t712 * t590;
t673 = -t686 * mrSges(4,1) + t687 * mrSges(4,2);
t680 = -qJD(3) * mrSges(4,2) + t686 * mrSges(4,3);
t578 = m(4) * t642 + qJDD(3) * mrSges(4,1) - t677 * mrSges(4,3) + qJD(3) * t680 - t687 * t673 + t580;
t681 = qJD(3) * mrSges(4,1) - t687 * mrSges(4,3);
t729 = t712 * t585 - t709 * t590;
t579 = m(4) * t643 - qJDD(3) * mrSges(4,2) + t676 * mrSges(4,3) - qJD(3) * t681 + t686 * t673 + t729;
t572 = t713 * t578 + t710 * t579;
t678 = -t706 * t688 + t733;
t720 = mrSges(3,3) * qJDD(1) + t715 * (-mrSges(3,1) * t707 + t745);
t570 = m(3) * t678 - t720 * t706 + t572;
t730 = -t710 * t578 + t713 * t579;
t571 = m(3) * t679 + t720 * t707 + t730;
t731 = -t706 * t570 + t707 * t571;
t565 = m(2) * t692 - t715 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t731;
t684 = -qJDD(1) * pkin(1) - t715 * qJ(2) + t726;
t588 = t708 * t593 + t751 * t595;
t719 = m(5) * t630 - t637 * mrSges(5,1) + t638 * mrSges(5,2) - t669 * t659 + t670 * t660 + t588;
t718 = m(4) * t675 - t676 * mrSges(4,1) + t677 * mrSges(4,2) - t686 * t680 + t687 * t681 + t719;
t716 = -m(3) * t684 + mrSges(3,1) * t735 - t718 + (t703 * t715 + t744) * mrSges(3,3);
t582 = m(2) * t691 + t716 + (mrSges(2,1) - t745) * qJDD(1) - t715 * mrSges(2,2);
t743 = t711 * t565 + t714 * t582;
t566 = t707 * t570 + t706 * t571;
t742 = t753 * t655 - t748 * t656 - t746 * t665;
t741 = t746 * t655 + t754 * t656 + t752 * t665;
t740 = -t748 * t655 + t755 * t656 - t754 * t665;
t723 = Ifges(3,5) * t706 + Ifges(3,6) * t707;
t738 = t715 * t723;
t732 = t714 * t565 - t711 * t582;
t725 = Ifges(3,1) * t706 + Ifges(3,4) * t707;
t724 = Ifges(3,4) * t706 + Ifges(3,2) * t707;
t668 = Ifges(4,1) * t687 + Ifges(4,4) * t686 + Ifges(4,5) * qJD(3);
t667 = Ifges(4,4) * t687 + Ifges(4,2) * t686 + Ifges(4,6) * qJD(3);
t666 = Ifges(4,5) * t687 + Ifges(4,6) * t686 + Ifges(4,3) * qJD(3);
t650 = Ifges(5,1) * t670 + Ifges(5,4) * t669 + Ifges(5,5) * t705;
t649 = Ifges(5,4) * t670 + Ifges(5,2) * t669 + Ifges(5,6) * t705;
t648 = Ifges(5,5) * t670 + Ifges(5,6) * t669 + Ifges(5,3) * t705;
t587 = mrSges(6,2) * t604 + mrSges(7,2) * t599 - mrSges(6,3) * t601 - mrSges(7,3) * t600 - qJ(6) * t596 - t748 * t614 + t755 * t615 - t636 * t754 + t741 * t655 + t742 * t665;
t586 = -mrSges(6,1) * t604 - mrSges(7,1) * t600 + mrSges(7,2) * t598 + mrSges(6,3) * t602 - pkin(5) * t596 - t753 * t614 + t748 * t615 + t746 * t636 + t741 * t656 + t740 * t665;
t574 = Ifges(5,4) * t638 + Ifges(5,2) * t637 + Ifges(5,6) * t702 - t670 * t648 + t705 * t650 - mrSges(5,1) * t630 + mrSges(5,3) * t609 - mrSges(6,1) * t601 + mrSges(6,2) * t602 + mrSges(7,1) * t599 - mrSges(7,3) * t598 - pkin(5) * t727 - qJ(6) * t734 - pkin(4) * t588 + (pkin(5) * t640 + t742) * t656 + (qJ(6) * t640 - t740) * t655 + t752 * t636 + (pkin(5) * mrSges(7,2) + t754) * t615 + (qJ(6) * mrSges(7,2) + t746) * t614;
t573 = mrSges(5,2) * t630 - mrSges(5,3) * t608 + Ifges(5,1) * t638 + Ifges(5,4) * t637 + Ifges(5,5) * t702 - pkin(9) * t588 - t708 * t586 + t751 * t587 + t669 * t648 - t705 * t649;
t562 = mrSges(4,2) * t675 - mrSges(4,3) * t642 + Ifges(4,1) * t677 + Ifges(4,4) * t676 + Ifges(4,5) * qJDD(3) - pkin(8) * t580 - qJD(3) * t667 + t712 * t573 - t709 * t574 + t686 * t666;
t561 = -mrSges(4,1) * t675 + mrSges(4,3) * t643 + Ifges(4,4) * t677 + Ifges(4,2) * t676 + Ifges(4,6) * qJDD(3) - pkin(3) * t719 + pkin(8) * t729 + qJD(3) * t668 + t709 * t573 + t712 * t574 - t687 * t666;
t560 = -pkin(9) * t728 + (Ifges(2,6) - t723) * qJDD(1) - t751 * t586 + t686 * t668 - t687 * t667 + mrSges(2,3) * t692 - Ifges(4,5) * t677 - mrSges(3,1) * t678 + mrSges(3,2) * t679 + t669 * t650 - t670 * t649 - Ifges(4,6) * t676 - Ifges(5,5) * t638 - mrSges(4,1) * t642 + mrSges(4,2) * t643 - Ifges(5,6) * t637 + mrSges(5,2) * t609 - mrSges(5,1) * t608 - pkin(1) * t566 - Ifges(4,3) * qJDD(3) + mrSges(2,1) * g(3) - t708 * t587 - pkin(4) * t717 - Ifges(5,3) * t702 - pkin(3) * t580 - pkin(2) * t572 + (-t706 * t724 + t707 * t725 + Ifges(2,5)) * t715;
t559 = mrSges(3,2) * t684 - mrSges(3,3) * t678 - pkin(7) * t572 + t725 * qJDD(1) - t710 * t561 + t713 * t562 + t707 * t738;
t558 = -mrSges(3,1) * t684 + mrSges(3,3) * t679 - pkin(2) * t718 + pkin(7) * t730 + t724 * qJDD(1) + t713 * t561 + t710 * t562 - t706 * t738;
t557 = -mrSges(2,2) * g(3) - mrSges(2,3) * t691 + Ifges(2,5) * qJDD(1) - t715 * Ifges(2,6) - qJ(2) * t566 - t706 * t558 + t707 * t559;
t1 = [-m(1) * g(1) + t732; -m(1) * g(2) + t743; (-m(1) - m(2)) * g(3) + t566; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t743 + t714 * t557 - t711 * t560; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t732 + t711 * t557 + t714 * t560; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t691 - mrSges(2,2) * t692 + t706 * t559 + t707 * t558 + pkin(1) * (-qJDD(1) * t745 + t716) + qJ(2) * t731;];
tauB  = t1;
