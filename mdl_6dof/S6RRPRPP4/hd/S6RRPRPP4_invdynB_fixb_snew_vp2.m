% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-05-06 12:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:35:57
% EndTime: 2019-05-06 12:36:05
% DurationCPUTime: 5.63s
% Computational Cost: add. (60855->345), mult. (127613->405), div. (0->0), fcn. (76675->8), ass. (0->134)
t783 = -2 * qJD(3);
t782 = Ifges(3,1) + Ifges(4,2);
t781 = Ifges(6,1) + Ifges(7,1);
t773 = Ifges(3,4) + Ifges(4,6);
t772 = Ifges(6,4) - Ifges(7,5);
t771 = Ifges(7,4) + Ifges(6,5);
t770 = Ifges(3,5) - Ifges(4,4);
t780 = Ifges(3,2) + Ifges(4,3);
t779 = Ifges(6,2) + Ifges(7,3);
t778 = -Ifges(7,2) - Ifges(6,3);
t769 = Ifges(3,6) - Ifges(4,5);
t768 = Ifges(6,6) - Ifges(7,6);
t777 = (Ifges(3,3) + Ifges(4,1));
t731 = sin(qJ(1));
t734 = cos(qJ(1));
t716 = -g(1) * t734 - g(2) * t731;
t736 = qJD(1) ^ 2;
t691 = -pkin(1) * t736 + qJDD(1) * pkin(7) + t716;
t730 = sin(qJ(2));
t733 = cos(qJ(2));
t677 = -t730 * g(3) + t733 * t691;
t703 = (-pkin(2) * t733 - qJ(3) * t730) * qJD(1);
t735 = qJD(2) ^ 2;
t756 = qJD(1) * t733;
t654 = t735 * pkin(2) - qJDD(2) * qJ(3) + (qJD(2) * t783) - t703 * t756 - t677;
t776 = -2 * qJD(5);
t775 = t736 * pkin(7);
t774 = -mrSges(6,3) - mrSges(7,2);
t767 = cos(pkin(9));
t676 = -t733 * g(3) - t730 * t691;
t704 = (mrSges(4,2) * t733 - mrSges(4,3) * t730) * qJD(1);
t705 = (-mrSges(3,1) * t733 + mrSges(3,2) * t730) * qJD(1);
t755 = qJD(1) * qJD(2);
t752 = t733 * t755;
t706 = qJDD(1) * t730 + t752;
t711 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t756;
t712 = -mrSges(4,1) * t756 - qJD(2) * mrSges(4,3);
t751 = t730 * t755;
t707 = qJDD(1) * t733 - t751;
t757 = qJD(1) * t730;
t714 = pkin(3) * t757 - (qJD(2) * pkin(8));
t727 = t733 ^ 2;
t715 = t731 * g(1) - t734 * g(2);
t746 = -qJDD(1) * pkin(1) - t715;
t741 = pkin(2) * t751 + t757 * t783 + (-t706 - t752) * qJ(3) + t746;
t630 = -t714 * t757 + (-pkin(3) * t727 - pkin(7)) * t736 + (-pkin(2) - pkin(8)) * t707 + t741;
t655 = -qJDD(2) * pkin(2) - t735 * qJ(3) + t703 * t757 + qJDD(3) - t676;
t637 = (-t730 * t733 * t736 - qJDD(2)) * pkin(8) + (t706 - t752) * pkin(3) + t655;
t729 = sin(qJ(4));
t732 = cos(qJ(4));
t625 = -t729 * t630 + t732 * t637;
t701 = -qJD(2) * t729 - t732 * t756;
t669 = qJD(4) * t701 + qJDD(2) * t732 - t707 * t729;
t700 = qJDD(4) + t706;
t702 = qJD(2) * t732 - t729 * t756;
t719 = qJD(4) + t757;
t622 = (t701 * t719 - t669) * qJ(5) + (t701 * t702 + t700) * pkin(4) + t625;
t626 = t732 * t630 + t729 * t637;
t668 = -qJD(4) * t702 - qJDD(2) * t729 - t707 * t732;
t674 = pkin(4) * t719 - qJ(5) * t702;
t699 = t701 ^ 2;
t624 = -pkin(4) * t699 + qJ(5) * t668 - t674 * t719 + t626;
t728 = sin(pkin(9));
t670 = -t767 * t701 + t728 * t702;
t618 = t728 * t622 + t767 * t624 + t670 * t776;
t644 = -t767 * t668 + t728 * t669;
t671 = t728 * t701 + t767 * t702;
t658 = mrSges(6,1) * t719 - mrSges(6,3) * t671;
t648 = pkin(5) * t670 - qJ(6) * t671;
t717 = t719 ^ 2;
t615 = -pkin(5) * t717 + qJ(6) * t700 + 0.2e1 * qJD(6) * t719 - t648 * t670 + t618;
t659 = -mrSges(7,1) * t719 + mrSges(7,2) * t671;
t753 = m(7) * t615 + t700 * mrSges(7,3) + t719 * t659;
t649 = mrSges(7,1) * t670 - mrSges(7,3) * t671;
t761 = -mrSges(6,1) * t670 - mrSges(6,2) * t671 - t649;
t607 = m(6) * t618 - t700 * mrSges(6,2) + t774 * t644 - t719 * t658 + t761 * t670 + t753;
t744 = t767 * t622 - t728 * t624;
t617 = t671 * t776 + t744;
t645 = t728 * t668 + t767 * t669;
t656 = -mrSges(6,2) * t719 - mrSges(6,3) * t670;
t616 = -t700 * pkin(5) - t717 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t648) * t671 - t744;
t657 = -mrSges(7,2) * t670 + mrSges(7,3) * t719;
t747 = -m(7) * t616 + t700 * mrSges(7,1) + t719 * t657;
t609 = m(6) * t617 + t700 * mrSges(6,1) + t774 * t645 + t719 * t656 + t761 * t671 + t747;
t602 = t728 * t607 + t767 * t609;
t672 = -mrSges(5,1) * t701 + mrSges(5,2) * t702;
t673 = -mrSges(5,2) * t719 + mrSges(5,3) * t701;
t600 = m(5) * t625 + mrSges(5,1) * t700 - mrSges(5,3) * t669 - t672 * t702 + t673 * t719 + t602;
t675 = mrSges(5,1) * t719 - mrSges(5,3) * t702;
t748 = t767 * t607 - t609 * t728;
t601 = m(5) * t626 - mrSges(5,2) * t700 + mrSges(5,3) * t668 + t672 * t701 - t675 * t719 + t748;
t597 = t732 * t600 + t729 * t601;
t742 = -m(4) * t655 - t706 * mrSges(4,1) - t597;
t595 = m(3) * t676 - t706 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t711 - t712) * qJD(2) + (-t704 - t705) * t757 + t742;
t710 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t757;
t713 = mrSges(4,1) * t757 + qJD(2) * mrSges(4,2);
t636 = -t727 * t736 * pkin(8) + t707 * pkin(3) + qJD(2) * t714 - t654;
t628 = -t668 * pkin(4) - t699 * qJ(5) + t702 * t674 + qJDD(5) + t636;
t620 = t628 + (t671 * t719 + t644) * pkin(5) + (t670 * t719 - t645) * qJ(6) - 0.2e1 * qJD(6) * t671;
t613 = m(7) * t620 + t644 * mrSges(7,1) - t645 * mrSges(7,3) + t670 * t657 - t671 * t659;
t740 = m(6) * t628 + t644 * mrSges(6,1) + t645 * mrSges(6,2) + t670 * t656 + t671 * t658 + t613;
t738 = -m(5) * t636 + t668 * mrSges(5,1) - t669 * mrSges(5,2) + t701 * t673 - t702 * t675 - t740;
t737 = -m(4) * t654 + qJDD(2) * mrSges(4,3) + qJD(2) * t713 + t704 * t756 - t738;
t612 = t737 - qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t707 + t705 * t756 + m(3) * t677 - qJD(2) * t710;
t749 = -t595 * t730 + t733 * t612;
t590 = m(2) * t716 - mrSges(2,1) * t736 - qJDD(1) * mrSges(2,2) + t749;
t690 = t746 - t775;
t651 = -t707 * pkin(2) + t741 - t775;
t765 = -t729 * t600 + t732 * t601;
t745 = -m(4) * t651 - t707 * mrSges(4,2) + t713 * t757 - t765;
t739 = -m(3) * t690 + t711 * t756 + t707 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t706 + (-t710 * t730 - t712 * t733) * qJD(1) + t745;
t593 = m(2) * t715 + qJDD(1) * mrSges(2,1) - t736 * mrSges(2,2) + t739;
t766 = t731 * t590 + t734 * t593;
t591 = t733 * t595 + t730 * t612;
t764 = t779 * t670 - t772 * t671 - t768 * t719;
t763 = t768 * t670 - t771 * t671 + t778 * t719;
t762 = -t772 * t670 + t781 * t671 + t771 * t719;
t760 = (t777 * qJD(2)) + (t770 * t730 + t769 * t733) * qJD(1);
t759 = -t769 * qJD(2) + (-t773 * t730 - t780 * t733) * qJD(1);
t758 = t770 * qJD(2) + (t782 * t730 + t773 * t733) * qJD(1);
t750 = t734 * t590 - t593 * t731;
t663 = Ifges(5,1) * t702 + Ifges(5,4) * t701 + Ifges(5,5) * t719;
t662 = Ifges(5,4) * t702 + Ifges(5,2) * t701 + Ifges(5,6) * t719;
t661 = Ifges(5,5) * t702 + Ifges(5,6) * t701 + Ifges(5,3) * t719;
t604 = mrSges(6,2) * t628 + mrSges(7,2) * t616 - mrSges(6,3) * t617 - mrSges(7,3) * t620 - qJ(6) * t613 - t772 * t644 + t781 * t645 + t763 * t670 + t771 * t700 + t764 * t719;
t603 = -mrSges(6,1) * t628 - mrSges(7,1) * t620 + mrSges(7,2) * t615 + mrSges(6,3) * t618 - pkin(5) * t613 - t779 * t644 + t772 * t645 + t763 * t671 + t768 * t700 + t762 * t719;
t596 = -t706 * mrSges(4,3) + t712 * t756 - t745;
t587 = mrSges(5,2) * t636 - mrSges(5,3) * t625 + Ifges(5,1) * t669 + Ifges(5,4) * t668 + Ifges(5,5) * t700 - qJ(5) * t602 - t728 * t603 + t767 * t604 + t701 * t661 - t719 * t662;
t586 = -mrSges(5,1) * t636 + mrSges(5,3) * t626 + Ifges(5,4) * t669 + Ifges(5,2) * t668 + Ifges(5,6) * t700 - pkin(4) * t740 + qJ(5) * t748 + t767 * t603 + t728 * t604 - t702 * t661 + t719 * t663;
t585 = pkin(4) * t602 + (Ifges(5,3) - t778) * t700 + t782 * t706 + (-pkin(5) * t649 - t764) * t671 + t759 * qJD(2) + t760 * t756 + (-qJ(6) * t649 + t762) * t670 + pkin(3) * t597 - qJ(3) * t596 + qJ(6) * t753 + pkin(5) * t747 + (-mrSges(7,2) * qJ(6) - t768) * t644 + t770 * qJDD(2) + (-mrSges(7,2) * pkin(5) + t771) * t645 + t773 * t707 + mrSges(7,3) * t615 - mrSges(7,1) * t616 + mrSges(6,1) * t617 - mrSges(6,2) * t618 + mrSges(5,1) * t625 - mrSges(5,2) * t626 - mrSges(4,3) * t651 + mrSges(4,1) * t655 + Ifges(5,6) * t668 + Ifges(5,5) * t669 - mrSges(3,3) * t676 + mrSges(3,2) * t690 - t701 * t663 + t702 * t662;
t584 = -mrSges(3,1) * t690 - mrSges(4,1) * t654 + mrSges(4,2) * t651 + mrSges(3,3) * t677 - pkin(2) * t596 - pkin(3) * t738 - pkin(8) * t765 + t758 * qJD(2) + t769 * qJDD(2) - t732 * t586 - t729 * t587 + t773 * t706 + t780 * t707 - t760 * t757;
t583 = -pkin(1) * t591 + mrSges(2,3) * t716 - pkin(2) * (-qJD(2) * t712 + t742) - qJ(3) * t737 + t729 * t586 + pkin(8) * t597 - mrSges(3,1) * t676 + mrSges(3,2) * t677 - t732 * t587 - mrSges(4,2) * t655 + mrSges(4,3) * t654 + mrSges(2,1) * g(3) + t736 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-mrSges(4,1) * qJ(3) - t769) * t707 - t770 * t706 + (mrSges(4,2) * pkin(2) - t777) * qJDD(2) + (t758 * t733 + (pkin(2) * t704 + t759) * t730) * qJD(1);
t582 = -mrSges(2,2) * g(3) - mrSges(2,3) * t715 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t736 - pkin(7) * t591 - t584 * t730 + t585 * t733;
t1 = [-m(1) * g(1) + t750; -m(1) * g(2) + t766; (-m(1) - m(2)) * g(3) + t591; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t766 + t734 * t582 - t731 * t583; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t750 + t731 * t582 + t734 * t583; -mrSges(1,1) * g(2) + mrSges(2,1) * t715 + mrSges(1,2) * g(1) - mrSges(2,2) * t716 + Ifges(2,3) * qJDD(1) + pkin(1) * t739 + pkin(7) * t749 + t733 * t584 + t730 * t585;];
tauB  = t1;
