% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-05-05 17:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPPR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:22:25
% EndTime: 2019-05-05 17:22:29
% DurationCPUTime: 2.50s
% Computational Cost: add. (12187->314), mult. (24258->348), div. (0->0), fcn. (10744->6), ass. (0->131)
t697 = sin(qJ(3));
t700 = cos(qJ(3));
t732 = qJD(1) * qJD(3);
t660 = qJDD(1) * t697 + t700 * t732;
t731 = (qJD(2) * qJD(1));
t680 = 2 * t731;
t678 = t697 * t732;
t661 = qJDD(1) * t700 - t678;
t703 = qJD(1) ^ 2;
t698 = sin(qJ(1));
t701 = cos(qJ(1));
t673 = -t701 * g(1) - t698 * g(2);
t721 = t703 * pkin(1) - qJDD(1) * qJ(2) - t673;
t719 = -t703 * pkin(7) - t721;
t716 = t660 * pkin(3) - t661 * qJ(4) + t719;
t720 = pkin(3) * t700 + qJ(4) * t697;
t734 = qJD(4) * t700;
t608 = t680 + (qJD(3) * t720 - 0.2e1 * t734) * qJD(1) + t716;
t736 = qJD(1) * t700;
t668 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t736;
t670 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t736;
t737 = qJD(1) * t697;
t671 = -mrSges(5,2) * t737 + qJD(3) * mrSges(5,3);
t667 = -qJD(3) * pkin(4) - qJ(5) * t736;
t681 = -2 * t731;
t743 = t697 ^ 2 * t703;
t707 = -qJ(5) * t743 + 0.2e1 * qJD(1) * t734 + t667 * t736 + qJDD(5) + t681 - t716;
t747 = -pkin(5) - qJ(4);
t756 = -pkin(4) - pkin(8);
t757 = -pkin(3) - pkin(8);
t600 = t707 + t756 * t660 + (t697 * t747 + t700 * t757) * t732 + t661 * pkin(5);
t672 = t698 * g(1) - t701 * g(2);
t715 = -t703 * qJ(2) + qJDD(2) - t672;
t627 = (-pkin(1) - pkin(7)) * qJDD(1) + t715;
t659 = (pkin(5) * t700 - pkin(8) * t697) * qJD(1);
t702 = qJD(3) ^ 2;
t754 = t697 * g(3);
t758 = -2 * qJD(5);
t656 = (pkin(3) * t697 - qJ(4) * t700) * qJD(1);
t763 = t656 * t736 + qJDD(4);
t709 = t700 * t703 * t697 * pkin(4) + t736 * t758 + (-t661 - t678) * qJ(5) - t754 + t763;
t603 = t747 * t702 + (-qJD(1) * t659 - t627) * t700 + (-pkin(3) + t756) * qJDD(3) + t709;
t696 = sin(qJ(6));
t699 = cos(qJ(6));
t598 = t600 * t699 - t603 * t696;
t653 = -qJD(3) * t699 - t696 * t737;
t618 = qJD(6) * t653 - qJDD(3) * t696 + t660 * t699;
t654 = -qJD(3) * t696 + t699 * t737;
t619 = -mrSges(7,1) * t653 + mrSges(7,2) * t654;
t675 = qJD(6) + t736;
t622 = -mrSges(7,2) * t675 + mrSges(7,3) * t653;
t651 = qJDD(6) + t661;
t596 = m(7) * t598 + mrSges(7,1) * t651 - mrSges(7,3) * t618 - t619 * t654 + t622 * t675;
t599 = t600 * t696 + t603 * t699;
t617 = -qJD(6) * t654 - qJDD(3) * t699 - t660 * t696;
t623 = mrSges(7,1) * t675 - mrSges(7,3) * t654;
t597 = m(7) * t599 - mrSges(7,2) * t651 + mrSges(7,3) * t617 + t619 * t653 - t623 * t675;
t587 = t699 * t596 + t696 * t597;
t604 = -t660 * pkin(4) - t720 * t732 + t707;
t717 = -m(6) * t604 - t661 * mrSges(6,1) - t587;
t665 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t737;
t744 = t665 * t697;
t704 = -((t668 + t670) * t700 + t744) * qJD(1) + m(5) * t608 + t660 * mrSges(5,1) - t661 * mrSges(5,3) + t671 * t737 + t717;
t585 = -t660 * mrSges(6,2) + t704;
t626 = t680 + t719;
t666 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t737;
t669 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t736;
t760 = -m(4) * t626 - t661 * mrSges(4,2) - t666 * t737 - t669 * t736;
t768 = t660 * mrSges(4,1) + t585 - t760;
t767 = Ifges(4,1) + Ifges(5,1) + Ifges(6,2);
t730 = Ifges(4,5) + Ifges(5,4) + Ifges(6,6);
t729 = Ifges(5,5) - Ifges(6,4) - Ifges(4,4);
t728 = Ifges(4,6) - Ifges(5,6) + Ifges(6,5);
t766 = Ifges(5,3) + Ifges(6,1) + Ifges(4,2);
t765 = -Ifges(6,3) - Ifges(4,3) - Ifges(5,2);
t621 = -t700 * g(3) + t697 * t627;
t764 = -qJDD(3) * qJ(4) - t621;
t628 = t681 + t721;
t762 = -m(3) * t628 + t703 * mrSges(3,2) + qJDD(1) * mrSges(3,3);
t753 = t702 * pkin(3);
t752 = mrSges(2,1) - mrSges(3,2);
t751 = mrSges(4,1) - mrSges(6,2);
t750 = -mrSges(4,3) - mrSges(5,2);
t749 = Ifges(2,5) - Ifges(3,4);
t748 = -Ifges(2,6) + Ifges(3,5);
t742 = t700 * t627;
t741 = t702 * qJ(4);
t620 = t742 + t754;
t611 = -qJDD(3) * pkin(3) - t620 - t741 + t763;
t606 = -t741 - t742 + (-pkin(3) - pkin(4)) * qJDD(3) + t709;
t655 = (mrSges(6,1) * t700 + mrSges(6,2) * t697) * qJD(1);
t739 = -t696 * t596 + t699 * t597;
t718 = -m(6) * t606 + t655 * t736 - t739;
t711 = -m(5) * t611 + qJDD(3) * mrSges(5,1) + qJD(3) * t671 + t718;
t657 = (mrSges(5,1) * t697 - mrSges(5,3) * t700) * qJD(1);
t722 = qJD(1) * (-t657 - (mrSges(4,1) * t697 + mrSges(4,2) * t700) * qJD(1));
t584 = m(4) * t620 + t751 * qJDD(3) + (-t665 + t666) * qJD(3) + t700 * t722 + (mrSges(6,3) + t750) * t661 + t711;
t679 = 0.2e1 * qJD(4) * qJD(3);
t610 = -t656 * t737 + t679 - t753 - t764;
t710 = pkin(4) * t743 - t660 * qJ(5) + t764;
t733 = t758 + t656;
t605 = t753 + (-0.2e1 * qJD(4) - t667) * qJD(3) + t733 * t737 + t710;
t602 = qJDD(3) * pkin(5) + qJD(3) * t667 + t679 + t757 * t702 + (t659 - t733) * t737 - t710;
t714 = -m(7) * t602 + t617 * mrSges(7,1) - t618 * mrSges(7,2) + t653 * t622 - t654 * t623;
t708 = -m(6) * t605 + qJDD(3) * mrSges(6,1) + t660 * mrSges(6,3) + qJD(3) * t668 + t655 * t737 - t714;
t706 = m(5) * t610 + qJDD(3) * mrSges(5,3) + qJD(3) * t670 + t708;
t589 = m(4) * t621 - qJDD(3) * mrSges(4,2) - qJD(3) * t669 + t660 * t750 + t697 * t722 + t706;
t580 = t700 * t584 + t697 * t589;
t629 = -qJDD(1) * pkin(1) + t715;
t712 = -m(3) * t629 + t703 * mrSges(3,3) - t580;
t578 = m(2) * t672 - t703 * mrSges(2,2) + qJDD(1) * t752 + t712;
t583 = m(2) * t673 - t703 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t660 * t751 + t704 - t760 + t762;
t740 = t701 * t578 + t698 * t583;
t735 = qJD(3) * t665;
t727 = t765 * qJD(3) + (t728 * t697 - t730 * t700) * qJD(1);
t726 = -t728 * qJD(3) + (t766 * t697 + t729 * t700) * qJD(1);
t725 = t730 * qJD(3) + (t729 * t697 + t767 * t700) * qJD(1);
t724 = -t578 * t698 + t701 * t583;
t723 = -t697 * t584 + t700 * t589;
t614 = Ifges(7,1) * t654 + Ifges(7,4) * t653 + Ifges(7,5) * t675;
t613 = Ifges(7,4) * t654 + Ifges(7,2) * t653 + Ifges(7,6) * t675;
t612 = Ifges(7,5) * t654 + Ifges(7,6) * t653 + Ifges(7,3) * t675;
t591 = mrSges(7,2) * t602 - mrSges(7,3) * t598 + Ifges(7,1) * t618 + Ifges(7,4) * t617 + Ifges(7,5) * t651 + t612 * t653 - t613 * t675;
t590 = -mrSges(7,1) * t602 + mrSges(7,3) * t599 + Ifges(7,4) * t618 + Ifges(7,2) * t617 + Ifges(7,6) * t651 - t612 * t654 + t614 * t675;
t586 = qJDD(3) * mrSges(6,2) - t661 * mrSges(6,3) - t718 + t735;
t579 = -m(3) * g(3) + t723;
t576 = -qJ(4) * t585 - qJ(5) * t586 + pkin(5) * t587 + mrSges(7,1) * t598 - mrSges(7,2) * t599 + mrSges(6,1) * t604 - mrSges(6,3) * t606 - mrSges(5,3) * t608 + mrSges(5,2) * t611 + Ifges(7,6) * t617 + Ifges(7,5) * t618 - mrSges(4,3) * t620 + mrSges(4,2) * t626 + Ifges(7,3) * t651 - t653 * t614 + t654 * t613 + t767 * t661 + t729 * t660 + t730 * qJDD(3) + t726 * qJD(3) + t727 * t737;
t575 = -pkin(3) * t585 + pkin(8) * t587 - mrSges(6,2) * t604 + mrSges(6,3) * t605 - mrSges(5,1) * t608 + mrSges(5,2) * t610 + mrSges(4,3) * t621 - mrSges(4,1) * t626 - qJ(5) * t708 + t696 * t590 - t699 * t591 - pkin(4) * t717 - t729 * t661 + (mrSges(6,2) * pkin(4) - t766) * t660 + t728 * qJDD(3) + t725 * qJD(3) + (pkin(4) * t744 + (pkin(4) * t668 + t727) * t700) * qJD(1);
t574 = qJ(4) * t706 + pkin(3) * (t711 - t735) - qJ(2) * t579 + pkin(2) * t580 - pkin(4) * t586 - pkin(8) * t739 - mrSges(6,1) * t605 + mrSges(6,2) * t606 + mrSges(5,3) * t610 - mrSges(5,1) * t611 + mrSges(4,1) * t620 - mrSges(4,2) * t621 + mrSges(3,1) * t629 - pkin(5) * t714 - mrSges(2,3) * t672 - t696 * t591 - t699 * t590 + t748 * t703 + t749 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (-mrSges(5,2) * qJ(4) - t728) * t660 + (-mrSges(6,2) * pkin(3) - t765) * qJDD(3) + (pkin(3) * (-mrSges(5,2) + mrSges(6,3)) + t730) * t661 + ((-pkin(3) * t657 - t726) * t700 + (-qJ(4) * t657 + t725) * t697) * qJD(1);
t573 = -mrSges(3,1) * t628 + mrSges(2,3) * t673 - pkin(1) * t579 + pkin(2) * t768 - pkin(7) * t723 + t752 * g(3) - t748 * qJDD(1) - t700 * t575 - t697 * t576 + t749 * t703;
t1 = [-m(1) * g(1) + t724; -m(1) * g(2) + t740; (-m(1) - m(2) - m(3)) * g(3) + t723; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t740 - t698 * t573 + t701 * t574; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t724 + t701 * t573 + t698 * t574; pkin(1) * t712 - pkin(7) * t580 + mrSges(2,1) * t672 - mrSges(2,2) * t673 + t700 * t576 - t697 * t575 + mrSges(3,2) * t629 - mrSges(3,3) * t628 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (t762 + t768) * qJ(2);];
tauB  = t1;
