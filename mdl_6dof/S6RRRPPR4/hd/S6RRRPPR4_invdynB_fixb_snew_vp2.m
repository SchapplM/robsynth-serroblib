% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 04:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:41:53
% EndTime: 2019-05-07 04:42:15
% DurationCPUTime: 10.47s
% Computational Cost: add. (151254->364), mult. (314362->441), div. (0->0), fcn. (220278->10), ass. (0->142)
t786 = -2 * qJD(4);
t785 = Ifges(5,1) + Ifges(6,1);
t778 = Ifges(5,4) - Ifges(6,5);
t777 = Ifges(5,5) + Ifges(6,4);
t784 = Ifges(5,2) + Ifges(6,3);
t783 = -Ifges(6,2) - Ifges(5,3);
t782 = Ifges(6,6) - Ifges(5,6);
t741 = sin(qJ(1));
t745 = cos(qJ(1));
t729 = t741 * g(1) - t745 * g(2);
t747 = qJD(1) ^ 2;
t710 = -qJDD(1) * pkin(1) - t747 * pkin(7) - t729;
t740 = sin(qJ(2));
t744 = cos(qJ(2));
t765 = qJD(1) * qJD(2);
t763 = t744 * t765;
t724 = t740 * qJDD(1) + t763;
t764 = t740 * t765;
t725 = t744 * qJDD(1) - t764;
t672 = (-t724 - t763) * pkin(8) + (-t725 + t764) * pkin(2) + t710;
t730 = -t745 * g(1) - t741 * g(2);
t711 = -t747 * pkin(1) + qJDD(1) * pkin(7) + t730;
t702 = -t740 * g(3) + t744 * t711;
t723 = (-pkin(2) * t744 - pkin(8) * t740) * qJD(1);
t746 = qJD(2) ^ 2;
t766 = t744 * qJD(1);
t677 = -t746 * pkin(2) + qJDD(2) * pkin(8) + t723 * t766 + t702;
t739 = sin(qJ(3));
t743 = cos(qJ(3));
t645 = t743 * t672 - t739 * t677;
t768 = qJD(1) * t740;
t720 = t743 * qJD(2) - t739 * t768;
t693 = t720 * qJD(3) + t739 * qJDD(2) + t743 * t724;
t719 = qJDD(3) - t725;
t721 = t739 * qJD(2) + t743 * t768;
t733 = qJD(3) - t766;
t634 = (t720 * t733 - t693) * qJ(4) + (t720 * t721 + t719) * pkin(3) + t645;
t646 = t739 * t672 + t743 * t677;
t692 = -t721 * qJD(3) + t743 * qJDD(2) - t739 * t724;
t699 = t733 * pkin(3) - t721 * qJ(4);
t718 = t720 ^ 2;
t637 = -t718 * pkin(3) + t692 * qJ(4) - t733 * t699 + t646;
t737 = sin(pkin(10));
t775 = cos(pkin(10));
t696 = t737 * t720 + t775 * t721;
t625 = t775 * t634 - t737 * t637 + t696 * t786;
t661 = t737 * t692 + t775 * t693;
t701 = -t744 * g(3) - t740 * t711;
t753 = qJDD(2) * pkin(2) + t746 * pkin(8) - t723 * t768 + t701;
t751 = t692 * pkin(3) + t718 * qJ(4) - t721 * t699 - qJDD(4) + t753;
t695 = -t775 * t720 + t737 * t721;
t774 = t695 * t733;
t781 = (-t661 + t774) * qJ(5) - t751;
t780 = 2 * qJD(5);
t779 = -mrSges(5,3) - mrSges(6,2);
t722 = (-mrSges(3,1) * t744 + mrSges(3,2) * t740) * qJD(1);
t726 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t768;
t626 = t737 * t634 + t775 * t637 + t695 * t786;
t660 = -t775 * t692 + t737 * t693;
t679 = t733 * mrSges(5,1) - t696 * mrSges(5,3);
t667 = t695 * pkin(4) - t696 * qJ(5);
t732 = t733 ^ 2;
t623 = -t732 * pkin(4) + t719 * qJ(5) - t695 * t667 + t733 * t780 + t626;
t680 = -t733 * mrSges(6,1) + t696 * mrSges(6,2);
t624 = -t719 * pkin(4) - t732 * qJ(5) + t696 * t667 + qJDD(5) - t625;
t618 = (-t661 - t774) * pkin(9) + (t695 * t696 - t719) * pkin(5) + t624;
t682 = -t733 * pkin(5) - t696 * pkin(9);
t694 = t695 ^ 2;
t619 = -t694 * pkin(5) + t660 * pkin(9) + t733 * t682 + t623;
t738 = sin(qJ(6));
t742 = cos(qJ(6));
t616 = t742 * t618 - t738 * t619;
t665 = t742 * t695 - t738 * t696;
t632 = t665 * qJD(6) + t738 * t660 + t742 * t661;
t666 = t738 * t695 + t742 * t696;
t643 = -t665 * mrSges(7,1) + t666 * mrSges(7,2);
t731 = qJD(6) - t733;
t649 = -t731 * mrSges(7,2) + t665 * mrSges(7,3);
t715 = qJDD(6) - t719;
t614 = m(7) * t616 + t715 * mrSges(7,1) - t632 * mrSges(7,3) - t666 * t643 + t731 * t649;
t617 = t738 * t618 + t742 * t619;
t631 = -t666 * qJD(6) + t742 * t660 - t738 * t661;
t650 = t731 * mrSges(7,1) - t666 * mrSges(7,3);
t615 = m(7) * t617 - t715 * mrSges(7,2) + t631 * mrSges(7,3) + t665 * t643 - t731 * t650;
t758 = -t738 * t614 + t742 * t615;
t754 = m(6) * t623 + t719 * mrSges(6,3) + t733 * t680 + t758;
t668 = t695 * mrSges(6,1) - t696 * mrSges(6,3);
t769 = -t695 * mrSges(5,1) - t696 * mrSges(5,2) - t668;
t604 = m(5) * t626 - t719 * mrSges(5,2) + t779 * t660 - t733 * t679 + t769 * t695 + t754;
t678 = -t733 * mrSges(5,2) - t695 * mrSges(5,3);
t607 = t742 * t614 + t738 * t615;
t681 = -t695 * mrSges(6,2) + t733 * mrSges(6,3);
t752 = -m(6) * t624 + t719 * mrSges(6,1) + t733 * t681 - t607;
t606 = m(5) * t625 + t719 * mrSges(5,1) + t779 * t661 + t733 * t678 + t769 * t696 + t752;
t601 = t737 * t604 + t775 * t606;
t697 = -t720 * mrSges(4,1) + t721 * mrSges(4,2);
t698 = -t733 * mrSges(4,2) + t720 * mrSges(4,3);
t599 = m(4) * t645 + t719 * mrSges(4,1) - t693 * mrSges(4,3) - t721 * t697 + t733 * t698 + t601;
t700 = t733 * mrSges(4,1) - t721 * mrSges(4,3);
t759 = t775 * t604 - t737 * t606;
t600 = m(4) * t646 - t719 * mrSges(4,2) + t692 * mrSges(4,3) + t720 * t697 - t733 * t700 + t759;
t760 = -t739 * t599 + t743 * t600;
t594 = m(3) * t702 - qJDD(2) * mrSges(3,2) + t725 * mrSges(3,3) - qJD(2) * t726 + t722 * t766 + t760;
t727 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t766;
t628 = -0.2e1 * qJD(5) * t696 + (t696 * t733 + t660) * pkin(4) + t781;
t621 = -t694 * pkin(9) + (-pkin(4) - pkin(5)) * t660 + (-pkin(4) * t733 + t682 + t780) * t696 - t781;
t755 = -m(7) * t621 + t631 * mrSges(7,1) - t632 * mrSges(7,2) + t665 * t649 - t666 * t650;
t612 = m(6) * t628 + t660 * mrSges(6,1) - t661 * mrSges(6,3) - t696 * t680 + t695 * t681 + t755;
t749 = -m(5) * t751 + t660 * mrSges(5,1) + t661 * mrSges(5,2) + t695 * t678 + t696 * t679 + t612;
t748 = m(4) * t753 + t692 * mrSges(4,1) - t693 * mrSges(4,2) + t720 * t698 - t721 * t700 - t749;
t611 = m(3) * t701 + qJDD(2) * mrSges(3,1) - t724 * mrSges(3,3) + qJD(2) * t727 - t722 * t768 + t748;
t761 = t744 * t594 - t740 * t611;
t588 = m(2) * t730 - t747 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t761;
t595 = t743 * t599 + t739 * t600;
t750 = -m(3) * t710 + t725 * mrSges(3,1) - t724 * mrSges(3,2) - t726 * t768 + t727 * t766 - t595;
t591 = m(2) * t729 + qJDD(1) * mrSges(2,1) - t747 * mrSges(2,2) + t750;
t773 = t741 * t588 + t745 * t591;
t589 = t740 * t594 + t744 * t611;
t772 = t784 * t695 - t778 * t696 + t782 * t733;
t771 = -t782 * t695 - t777 * t696 + t783 * t733;
t770 = -t778 * t695 + t785 * t696 + t777 * t733;
t762 = t745 * t588 - t741 * t591;
t709 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t740 + Ifges(3,4) * t744) * qJD(1);
t708 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t740 + Ifges(3,2) * t744) * qJD(1);
t707 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t740 + Ifges(3,6) * t744) * qJD(1);
t686 = Ifges(4,1) * t721 + Ifges(4,4) * t720 + Ifges(4,5) * t733;
t685 = Ifges(4,4) * t721 + Ifges(4,2) * t720 + Ifges(4,6) * t733;
t684 = Ifges(4,5) * t721 + Ifges(4,6) * t720 + Ifges(4,3) * t733;
t640 = Ifges(7,1) * t666 + Ifges(7,4) * t665 + Ifges(7,5) * t731;
t639 = Ifges(7,4) * t666 + Ifges(7,2) * t665 + Ifges(7,6) * t731;
t638 = Ifges(7,5) * t666 + Ifges(7,6) * t665 + Ifges(7,3) * t731;
t609 = mrSges(7,2) * t621 - mrSges(7,3) * t616 + Ifges(7,1) * t632 + Ifges(7,4) * t631 + Ifges(7,5) * t715 + t665 * t638 - t731 * t639;
t608 = -mrSges(7,1) * t621 + mrSges(7,3) * t617 + Ifges(7,4) * t632 + Ifges(7,2) * t631 + Ifges(7,6) * t715 - t666 * t638 + t731 * t640;
t597 = -mrSges(5,2) * t751 + mrSges(6,2) * t624 - mrSges(5,3) * t625 - mrSges(6,3) * t628 - pkin(9) * t607 - qJ(5) * t612 - t738 * t608 + t742 * t609 - t778 * t660 + t785 * t661 + t771 * t695 + t777 * t719 + t772 * t733;
t596 = mrSges(5,1) * t751 - mrSges(6,1) * t628 + mrSges(6,2) * t623 + mrSges(5,3) * t626 - pkin(4) * t612 - pkin(5) * t755 - pkin(9) * t758 - t742 * t608 - t738 * t609 - t784 * t660 + t778 * t661 + t771 * t696 - t719 * t782 + t770 * t733;
t585 = -mrSges(4,2) * t753 - mrSges(4,3) * t645 + Ifges(4,1) * t693 + Ifges(4,4) * t692 + Ifges(4,5) * t719 - qJ(4) * t601 - t737 * t596 + t775 * t597 + t720 * t684 - t733 * t685;
t584 = mrSges(4,1) * t753 + mrSges(4,3) * t646 + Ifges(4,4) * t693 + Ifges(4,2) * t692 + Ifges(4,6) * t719 - pkin(3) * t749 + qJ(4) * t759 + t775 * t596 + t737 * t597 - t721 * t684 + t733 * t686;
t583 = (pkin(4) * mrSges(6,2) - t777) * t661 + (-Ifges(4,3) + t783) * t719 + (qJ(5) * t668 - t770) * t695 + (pkin(4) * t668 + t772) * t696 + (qJ(5) * mrSges(6,2) - t782) * t660 + Ifges(3,6) * qJDD(2) - pkin(3) * t601 - t665 * t640 + t666 * t639 + mrSges(3,3) * t702 + qJD(2) * t709 - mrSges(3,1) * t710 - Ifges(4,6) * t692 - Ifges(4,5) * t693 - pkin(2) * t595 + Ifges(7,3) * t715 + t720 * t686 - mrSges(7,2) * t617 - mrSges(4,1) * t645 + mrSges(4,2) * t646 + Ifges(7,6) * t631 + Ifges(7,5) * t632 + mrSges(6,1) * t624 - mrSges(5,1) * t625 + mrSges(5,2) * t626 + pkin(5) * t607 - t721 * t685 + Ifges(3,4) * t724 + Ifges(3,2) * t725 - mrSges(6,3) * t623 - pkin(4) * t752 - qJ(5) * t754 + mrSges(7,1) * t616 - t707 * t768;
t582 = mrSges(3,2) * t710 - mrSges(3,3) * t701 + Ifges(3,1) * t724 + Ifges(3,4) * t725 + Ifges(3,5) * qJDD(2) - pkin(8) * t595 - qJD(2) * t708 - t739 * t584 + t743 * t585 + t707 * t766;
t581 = Ifges(2,6) * qJDD(1) + t747 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t730 - Ifges(3,5) * t724 - Ifges(3,6) * t725 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t701 + mrSges(3,2) * t702 - t739 * t585 - t743 * t584 - pkin(2) * t748 - pkin(8) * t760 - pkin(1) * t589 + (-t740 * t708 + t744 * t709) * qJD(1);
t580 = -mrSges(2,2) * g(3) - mrSges(2,3) * t729 + Ifges(2,5) * qJDD(1) - t747 * Ifges(2,6) - pkin(7) * t589 + t744 * t582 - t740 * t583;
t1 = [-m(1) * g(1) + t762; -m(1) * g(2) + t773; (-m(1) - m(2)) * g(3) + t589; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t773 + t745 * t580 - t741 * t581; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t762 + t741 * t580 + t745 * t581; -mrSges(1,1) * g(2) + mrSges(2,1) * t729 + mrSges(1,2) * g(1) - mrSges(2,2) * t730 + Ifges(2,3) * qJDD(1) + pkin(1) * t750 + pkin(7) * t761 + t740 * t582 + t744 * t583;];
tauB  = t1;
