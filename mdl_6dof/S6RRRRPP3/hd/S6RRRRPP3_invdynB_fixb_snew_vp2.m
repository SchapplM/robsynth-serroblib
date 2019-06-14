% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPP3
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 18:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:08:19
% EndTime: 2019-05-07 18:08:31
% DurationCPUTime: 6.70s
% Computational Cost: add. (75309->344), mult. (150896->405), div. (0->0), fcn. (103808->8), ass. (0->131)
t746 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t745 = -Ifges(6,1) - Ifges(5,3) - Ifges(7,1);
t732 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t731 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t744 = Ifges(5,2) + Ifges(7,2) + Ifges(6,3);
t730 = Ifges(5,6) - Ifges(6,5) - Ifges(7,4);
t743 = -2 * qJD(5);
t742 = cos(qJ(4));
t710 = qJD(1) ^ 2;
t741 = pkin(2) * t710;
t740 = -mrSges(7,1) - mrSges(5,3);
t704 = sin(qJ(3));
t705 = sin(qJ(2));
t707 = cos(qJ(3));
t708 = cos(qJ(2));
t681 = (t704 * t708 + t705 * t707) * qJD(1);
t701 = qJD(2) + qJD(3);
t703 = sin(qJ(4));
t668 = t681 * t703 - t742 * t701;
t734 = qJD(1) * t708;
t735 = qJD(1) * t705;
t680 = -t704 * t735 + t707 * t734;
t676 = qJD(4) - t680;
t739 = t668 * t676;
t706 = sin(qJ(1));
t709 = cos(qJ(1));
t694 = -g(1) * t709 - g(2) * t706;
t683 = -pkin(1) * t710 + qJDD(1) * pkin(7) + t694;
t738 = t683 * t705;
t733 = qJD(1) * qJD(2);
t688 = qJDD(1) * t705 + t708 * t733;
t641 = qJDD(2) * pkin(2) - pkin(8) * t688 - t738 + (pkin(8) * t733 + t705 * t741 - g(3)) * t708;
t671 = -g(3) * t705 + t708 * t683;
t689 = qJDD(1) * t708 - t705 * t733;
t692 = qJD(2) * pkin(2) - pkin(8) * t735;
t702 = t708 ^ 2;
t642 = pkin(8) * t689 - qJD(2) * t692 - t702 * t741 + t671;
t608 = t704 * t641 + t707 * t642;
t649 = -qJD(3) * t681 - t688 * t704 + t689 * t707;
t664 = -mrSges(4,1) * t680 + mrSges(4,2) * t681;
t673 = mrSges(4,1) * t701 - mrSges(4,3) * t681;
t700 = qJDD(2) + qJDD(3);
t650 = qJD(3) * t680 + t688 * t707 + t689 * t704;
t693 = g(1) * t706 - t709 * g(2);
t719 = -qJDD(1) * pkin(1) - t693;
t651 = -pkin(2) * t689 + t692 * t735 + (-pkin(8) * t702 - pkin(7)) * t710 + t719;
t602 = (-t680 * t701 - t650) * pkin(9) + (t681 * t701 - t649) * pkin(3) + t651;
t665 = -pkin(3) * t680 - pkin(9) * t681;
t699 = t701 ^ 2;
t606 = -pkin(3) * t699 + pkin(9) * t700 + t665 * t680 + t608;
t599 = t742 * t602 - t703 * t606;
t615 = -t668 * qJD(4) + t742 * t650 + t703 * t700;
t669 = t742 * t681 + t703 * t701;
t637 = -mrSges(7,2) * t669 + mrSges(7,3) * t668;
t639 = mrSges(5,1) * t668 + mrSges(5,2) * t669;
t648 = qJDD(4) - t649;
t654 = mrSges(6,1) * t668 - mrSges(6,3) * t676;
t657 = -mrSges(5,2) * t676 - mrSges(5,3) * t668;
t638 = pkin(4) * t668 - qJ(5) * t669;
t675 = t676 ^ 2;
t597 = -t648 * pkin(4) - t675 * qJ(5) + t669 * t638 + qJDD(5) - t599;
t640 = -mrSges(6,2) * t668 - mrSges(6,3) * t669;
t591 = -0.2e1 * qJD(6) * t676 + (t668 * t669 - t648) * qJ(6) + (t615 + t739) * pkin(5) + t597;
t655 = -mrSges(7,1) * t668 + mrSges(7,2) * t676;
t721 = m(7) * t591 - t648 * mrSges(7,3) - t676 * t655;
t717 = -m(6) * t597 - t615 * mrSges(6,1) - t669 * t640 - t721;
t585 = m(5) * t599 + (-t654 + t657) * t676 + (-t637 - t639) * t669 + (mrSges(5,1) - mrSges(6,2)) * t648 + t740 * t615 + t717;
t600 = t703 * t602 + t742 * t606;
t614 = qJD(4) * t669 + t650 * t703 - t742 * t700;
t658 = mrSges(5,1) * t676 - mrSges(5,3) * t669;
t715 = -pkin(4) * t675 + qJ(5) * t648 - t638 * t668 + t600;
t596 = t676 * t743 - t715;
t656 = mrSges(6,1) * t669 + mrSges(6,2) * t676;
t652 = pkin(5) * t669 - qJ(6) * t676;
t667 = t668 ^ 2;
t593 = -pkin(5) * t614 - qJ(6) * t667 + qJDD(6) + ((2 * qJD(5)) + t652) * t676 + t715;
t653 = mrSges(7,1) * t669 - mrSges(7,3) * t676;
t729 = m(7) * t593 + t648 * mrSges(7,2) + t676 * t653;
t718 = -m(6) * t596 + t648 * mrSges(6,3) + t676 * t656 + t729;
t736 = -t637 - t640;
t588 = m(5) * t600 - mrSges(5,2) * t648 - t658 * t676 + (-t639 + t736) * t668 + (-mrSges(6,1) + t740) * t614 + t718;
t722 = -t585 * t703 + t742 * t588;
t580 = m(4) * t608 - mrSges(4,2) * t700 + mrSges(4,3) * t649 + t664 * t680 - t673 * t701 + t722;
t607 = t641 * t707 - t704 * t642;
t672 = -mrSges(4,2) * t701 + mrSges(4,3) * t680;
t605 = -pkin(3) * t700 - pkin(9) * t699 + t681 * t665 - t607;
t713 = (-t615 + t739) * qJ(5) + t605 + (pkin(4) * t676 + t743) * t669;
t598 = pkin(4) * t614 + t713;
t595 = -pkin(5) * t667 + 0.2e1 * qJD(6) * t668 - t652 * t669 + (pkin(4) + qJ(6)) * t614 + t713;
t720 = m(7) * t595 - t615 * mrSges(7,2) + t614 * mrSges(7,3) - t669 * t653 + t668 * t655;
t716 = -m(6) * t598 + t614 * mrSges(6,2) + t668 * t654 - t720;
t712 = -m(5) * t605 - t614 * mrSges(5,1) - t668 * t657 + (t656 - t658) * t669 + (-mrSges(5,2) + mrSges(6,3)) * t615 + t716;
t583 = m(4) * t607 + mrSges(4,1) * t700 - mrSges(4,3) * t650 - t664 * t681 + t672 * t701 + t712;
t574 = t704 * t580 + t707 * t583;
t670 = -g(3) * t708 - t738;
t687 = (-mrSges(3,1) * t708 + mrSges(3,2) * t705) * qJD(1);
t691 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t734;
t572 = m(3) * t670 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t688 + qJD(2) * t691 - t687 * t735 + t574;
t690 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t735;
t723 = t707 * t580 - t583 * t704;
t573 = m(3) * t671 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t689 - qJD(2) * t690 + t687 * t734 + t723;
t724 = -t572 * t705 + t708 * t573;
t565 = m(2) * t694 - mrSges(2,1) * t710 - qJDD(1) * mrSges(2,2) + t724;
t682 = -pkin(7) * t710 + t719;
t581 = t742 * t585 + t703 * t588;
t714 = m(4) * t651 - t649 * mrSges(4,1) + mrSges(4,2) * t650 - t680 * t672 + t673 * t681 + t581;
t711 = -m(3) * t682 + t689 * mrSges(3,1) - mrSges(3,2) * t688 - t690 * t735 + t691 * t734 - t714;
t577 = m(2) * t693 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t710 + t711;
t737 = t706 * t565 + t709 * t577;
t566 = t708 * t572 + t705 * t573;
t728 = t668 * t730 - t669 * t731 + t676 * t745;
t727 = t668 * t744 - t669 * t732 - t676 * t730;
t726 = -t732 * t668 + t669 * t746 + t731 * t676;
t725 = t709 * t565 - t577 * t706;
t679 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t705 + Ifges(3,4) * t708) * qJD(1);
t678 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t705 + Ifges(3,2) * t708) * qJD(1);
t677 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t705 + Ifges(3,6) * t708) * qJD(1);
t662 = Ifges(4,1) * t681 + Ifges(4,4) * t680 + Ifges(4,5) * t701;
t661 = Ifges(4,4) * t681 + Ifges(4,2) * t680 + Ifges(4,6) * t701;
t660 = Ifges(4,5) * t681 + Ifges(4,6) * t680 + Ifges(4,3) * t701;
t590 = mrSges(7,1) * t615 + t637 * t669 + t721;
t589 = -mrSges(6,3) * t615 - t656 * t669 - t716;
t575 = mrSges(6,1) * t597 + mrSges(7,1) * t591 + mrSges(5,2) * t605 - mrSges(7,2) * t595 - mrSges(5,3) * t599 - mrSges(6,3) * t598 + pkin(5) * t590 - qJ(5) * t589 - t732 * t614 + t615 * t746 + t731 * t648 + t728 * t668 + t727 * t676;
t568 = -mrSges(5,1) * t605 + mrSges(5,3) * t600 - mrSges(6,1) * t596 + mrSges(6,2) * t598 + mrSges(7,1) * t593 - mrSges(7,3) * t595 - pkin(5) * (t637 * t668 - t729) - qJ(6) * t720 - pkin(4) * t589 + t726 * t676 + t728 * t669 + t730 * t648 + t732 * t615 + (-mrSges(7,1) * pkin(5) - t744) * t614;
t567 = Ifges(4,6) * t700 + t701 * t662 - t681 * t660 - pkin(4) * (-t654 * t676 + t717) + Ifges(4,2) * t649 + Ifges(4,4) * t650 - mrSges(4,1) * t651 + mrSges(4,3) * t608 + mrSges(5,2) * t600 - mrSges(5,1) * t599 + mrSges(6,3) * t596 - mrSges(6,2) * t597 - mrSges(7,2) * t593 + qJ(6) * t590 + mrSges(7,3) * t591 - pkin(3) * t581 - qJ(5) * t718 + (pkin(4) * t637 + t727) * t669 + (pkin(4) * mrSges(6,2) + t745) * t648 + (pkin(4) * mrSges(7,1) - t731) * t615 + (-qJ(5) * t736 - t726) * t668 + (-qJ(5) * (-mrSges(6,1) - mrSges(7,1)) + t730) * t614;
t562 = mrSges(4,2) * t651 - mrSges(4,3) * t607 + Ifges(4,1) * t650 + Ifges(4,4) * t649 + Ifges(4,5) * t700 - pkin(9) * t581 - t703 * t568 + t742 * t575 + t680 * t660 - t701 * t661;
t561 = mrSges(3,2) * t682 - mrSges(3,3) * t670 + Ifges(3,1) * t688 + Ifges(3,4) * t689 + Ifges(3,5) * qJDD(2) - pkin(8) * t574 - qJD(2) * t678 + t562 * t707 - t567 * t704 + t677 * t734;
t560 = -pkin(1) * t566 + mrSges(2,3) * t694 - pkin(2) * t574 - Ifges(3,5) * t688 - Ifges(3,6) * t689 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t670 + mrSges(3,2) * t671 - t742 * t568 - pkin(3) * t712 - pkin(9) * t722 - Ifges(4,5) * t650 - Ifges(4,6) * t649 - Ifges(4,3) * t700 - mrSges(4,1) * t607 + mrSges(4,2) * t608 - t703 * t575 + mrSges(2,1) * g(3) + t710 * Ifges(2,5) - t681 * t661 + t680 * t662 + Ifges(2,6) * qJDD(1) + (-t705 * t678 + t708 * t679) * qJD(1);
t559 = -mrSges(3,1) * t682 + mrSges(3,3) * t671 + Ifges(3,4) * t688 + Ifges(3,2) * t689 + Ifges(3,6) * qJDD(2) - pkin(2) * t714 + pkin(8) * t723 + qJD(2) * t679 + t704 * t562 + t707 * t567 - t677 * t735;
t558 = -mrSges(2,2) * g(3) - mrSges(2,3) * t693 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t710 - pkin(7) * t566 - t559 * t705 + t561 * t708;
t1 = [-m(1) * g(1) + t725; -m(1) * g(2) + t737; (-m(1) - m(2)) * g(3) + t566; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t737 + t709 * t558 - t706 * t560; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t725 + t706 * t558 + t709 * t560; -mrSges(1,1) * g(2) + mrSges(2,1) * t693 + mrSges(1,2) * g(1) - mrSges(2,2) * t694 + Ifges(2,3) * qJDD(1) + pkin(1) * t711 + pkin(7) * t724 + t708 * t559 + t705 * t561;];
tauB  = t1;
