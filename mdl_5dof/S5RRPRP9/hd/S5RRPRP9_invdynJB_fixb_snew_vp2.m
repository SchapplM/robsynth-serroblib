% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRP9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP9_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP9_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP9_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:06:08
% EndTime: 2019-12-31 20:06:15
% DurationCPUTime: 4.89s
% Computational Cost: add. (50718->288), mult. (108789->352), div. (0->0), fcn. (71284->8), ass. (0->114)
t755 = Ifges(5,1) + Ifges(6,1);
t748 = Ifges(5,4) - Ifges(6,5);
t747 = -Ifges(5,5) - Ifges(6,4);
t754 = Ifges(5,2) + Ifges(6,3);
t746 = Ifges(5,6) - Ifges(6,6);
t753 = -Ifges(5,3) - Ifges(6,2);
t719 = sin(qJ(1));
t721 = cos(qJ(1));
t708 = -t721 * g(1) - t719 * g(2);
t723 = qJD(1) ^ 2;
t691 = -t723 * pkin(1) + qJDD(1) * pkin(6) + t708;
t718 = sin(qJ(2));
t720 = cos(qJ(2));
t674 = -t720 * g(3) - t718 * t691;
t700 = (-pkin(2) * t720 - qJ(3) * t718) * qJD(1);
t722 = qJD(2) ^ 2;
t739 = qJD(1) * t718;
t659 = -qJDD(2) * pkin(2) - t722 * qJ(3) + t700 * t739 + qJDD(3) - t674;
t737 = qJD(1) * qJD(2);
t735 = t720 * t737;
t702 = t718 * qJDD(1) + t735;
t715 = sin(pkin(8));
t716 = cos(pkin(8));
t678 = t716 * qJDD(2) - t715 * t702;
t697 = t715 * qJD(2) + t716 * t739;
t738 = t720 * qJD(1);
t680 = -pkin(3) * t738 - t697 * pkin(7);
t696 = t716 * qJD(2) - t715 * t739;
t695 = t696 ^ 2;
t637 = -t678 * pkin(3) - t695 * pkin(7) + t697 * t680 + t659;
t717 = sin(qJ(4));
t750 = cos(qJ(4));
t671 = t717 * t696 + t750 * t697;
t679 = t715 * qJDD(2) + t716 * t702;
t640 = t671 * qJD(4) - t750 * t678 + t717 * t679;
t670 = -t750 * t696 + t717 * t697;
t641 = -t670 * qJD(4) + t717 * t678 + t750 * t679;
t710 = qJD(4) - t738;
t628 = -0.2e1 * qJD(5) * t671 + (t670 * t710 - t641) * qJ(5) + (t671 * t710 + t640) * pkin(4) + t637;
t663 = -t710 * mrSges(6,1) + t671 * mrSges(6,2);
t664 = -t670 * mrSges(6,2) + t710 * mrSges(6,3);
t624 = m(6) * t628 + t640 * mrSges(6,1) - t641 * mrSges(6,3) - t671 * t663 + t670 * t664;
t707 = t719 * g(1) - t721 * g(2);
t690 = -qJDD(1) * pkin(1) - t723 * pkin(6) - t707;
t711 = t718 * t737;
t703 = t720 * qJDD(1) - t711;
t655 = (-t702 - t735) * qJ(3) + (-t703 + t711) * pkin(2) + t690;
t675 = -t718 * g(3) + t720 * t691;
t660 = -t722 * pkin(2) + qJDD(2) * qJ(3) + t700 * t738 + t675;
t635 = -0.2e1 * qJD(3) * t697 + t716 * t655 - t715 * t660;
t632 = (-t696 * t738 - t679) * pkin(7) + (t696 * t697 - t703) * pkin(3) + t635;
t636 = 0.2e1 * qJD(3) * t696 + t715 * t655 + t716 * t660;
t634 = -t695 * pkin(3) + t678 * pkin(7) + t680 * t738 + t636;
t630 = t717 * t632 + t750 * t634;
t650 = t670 * pkin(4) - t671 * qJ(5);
t699 = qJDD(4) - t703;
t709 = t710 ^ 2;
t626 = -t709 * pkin(4) + t699 * qJ(5) + 0.2e1 * qJD(5) * t710 - t670 * t650 + t630;
t741 = -t748 * t670 + t755 * t671 - t747 * t710;
t742 = t670 * t746 + t671 * t747 + t710 * t753;
t611 = -mrSges(5,1) * t637 - mrSges(6,1) * t628 + mrSges(6,2) * t626 + mrSges(5,3) * t630 - pkin(4) * t624 - t640 * t754 + t748 * t641 + t742 * t671 + t746 * t699 + t741 * t710;
t629 = t750 * t632 - t717 * t634;
t627 = -t699 * pkin(4) - t709 * qJ(5) + t671 * t650 + qJDD(5) - t629;
t743 = t754 * t670 - t671 * t748 - t746 * t710;
t612 = mrSges(5,2) * t637 + mrSges(6,2) * t627 - mrSges(5,3) * t629 - mrSges(6,3) * t628 - qJ(5) * t624 - t748 * t640 + t755 * t641 + t742 * t670 - t747 * t699 + t743 * t710;
t665 = Ifges(4,5) * t697 + Ifges(4,6) * t696 - Ifges(4,3) * t738;
t667 = Ifges(4,1) * t697 + Ifges(4,4) * t696 - Ifges(4,5) * t738;
t661 = -t710 * mrSges(5,2) - t670 * mrSges(5,3);
t662 = t710 * mrSges(5,1) - t671 * mrSges(5,3);
t725 = m(5) * t637 + t640 * mrSges(5,1) + t641 * mrSges(5,2) + t670 * t661 + t671 * t662 + t624;
t736 = m(6) * t626 + t699 * mrSges(6,3) + t710 * t663;
t651 = t670 * mrSges(6,1) - t671 * mrSges(6,3);
t740 = -t670 * mrSges(5,1) - t671 * mrSges(5,2) - t651;
t749 = -mrSges(5,3) - mrSges(6,2);
t616 = m(5) * t630 - t699 * mrSges(5,2) + t749 * t640 - t710 * t662 + t740 * t670 + t736;
t731 = -m(6) * t627 + t699 * mrSges(6,1) + t710 * t664;
t618 = m(5) * t629 + t699 * mrSges(5,1) + t749 * t641 + t710 * t661 + t740 * t671 + t731;
t732 = t750 * t616 - t717 * t618;
t593 = -mrSges(4,1) * t659 + mrSges(4,3) * t636 + Ifges(4,4) * t679 + Ifges(4,2) * t678 - Ifges(4,6) * t703 - pkin(3) * t725 + pkin(7) * t732 + t750 * t611 + t717 * t612 - t697 * t665 - t667 * t738;
t613 = t717 * t616 + t750 * t618;
t666 = Ifges(4,4) * t697 + Ifges(4,2) * t696 - Ifges(4,6) * t738;
t594 = mrSges(4,2) * t659 - mrSges(4,3) * t635 + Ifges(4,1) * t679 + Ifges(4,4) * t678 - Ifges(4,5) * t703 - pkin(7) * t613 - t717 * t611 + t750 * t612 + t696 * t665 + t666 * t738;
t672 = -t696 * mrSges(4,1) + t697 * mrSges(4,2);
t729 = mrSges(4,2) * t738 + t696 * mrSges(4,3);
t609 = m(4) * t635 - t703 * mrSges(4,1) - t679 * mrSges(4,3) - t697 * t672 - t729 * t738 + t613;
t677 = -mrSges(4,1) * t738 - t697 * mrSges(4,3);
t610 = m(4) * t636 + t703 * mrSges(4,2) + t678 * mrSges(4,3) + t696 * t672 + t677 * t738 + t732;
t607 = -t715 * t609 + t716 * t610;
t621 = m(4) * t659 - t678 * mrSges(4,1) + t679 * mrSges(4,2) + t697 * t677 - t696 * t729 + t725;
t688 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t718 + Ifges(3,2) * t720) * qJD(1);
t689 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t718 + Ifges(3,4) * t720) * qJD(1);
t752 = mrSges(3,1) * t674 - mrSges(3,2) * t675 + Ifges(3,5) * t702 + Ifges(3,6) * t703 + Ifges(3,3) * qJDD(2) - pkin(2) * t621 + qJ(3) * t607 + t716 * t593 + t715 * t594 + (t718 * t688 - t720 * t689) * qJD(1);
t623 = t641 * mrSges(6,2) + t671 * t651 - t731;
t751 = t746 * t640 + t747 * t641 - t741 * t670 + t743 * t671 + t753 * t699 - mrSges(5,1) * t629 + mrSges(6,1) * t627 + mrSges(5,2) * t630 - mrSges(6,3) * t626 + pkin(4) * t623 - qJ(5) * (-t640 * mrSges(6,2) - t670 * t651 + t736);
t701 = (-mrSges(3,1) * t720 + mrSges(3,2) * t718) * qJD(1);
t705 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t739;
t605 = m(3) * t675 - qJDD(2) * mrSges(3,2) + t703 * mrSges(3,3) - qJD(2) * t705 + t701 * t738 + t607;
t706 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t738;
t620 = m(3) * t674 + qJDD(2) * mrSges(3,1) - t702 * mrSges(3,3) + qJD(2) * t706 - t701 * t739 - t621;
t733 = t720 * t605 - t718 * t620;
t597 = m(2) * t708 - t723 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t733;
t606 = t716 * t609 + t715 * t610;
t726 = -m(3) * t690 + t703 * mrSges(3,1) - t702 * mrSges(3,2) - t705 * t739 + t706 * t738 - t606;
t601 = m(2) * t707 + qJDD(1) * mrSges(2,1) - t723 * mrSges(2,2) + t726;
t744 = t719 * t597 + t721 * t601;
t599 = t718 * t605 + t720 * t620;
t734 = t721 * t597 - t719 * t601;
t687 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t718 + Ifges(3,6) * t720) * qJD(1);
t590 = mrSges(3,2) * t690 - mrSges(3,3) * t674 + Ifges(3,1) * t702 + Ifges(3,4) * t703 + Ifges(3,5) * qJDD(2) - qJ(3) * t606 - qJD(2) * t688 - t715 * t593 + t716 * t594 + t687 * t738;
t592 = qJD(2) * t689 - mrSges(3,1) * t690 + t696 * t667 - t697 * t666 + mrSges(3,3) * t675 - Ifges(4,6) * t678 - Ifges(4,5) * t679 + Ifges(3,4) * t702 + t751 - mrSges(4,1) * t635 + mrSges(4,2) * t636 + Ifges(3,6) * qJDD(2) - pkin(3) * t613 - pkin(2) * t606 - t687 * t739 + (Ifges(4,3) + Ifges(3,2)) * t703;
t728 = mrSges(2,1) * t707 - mrSges(2,2) * t708 + Ifges(2,3) * qJDD(1) + pkin(1) * t726 + pkin(6) * t733 + t718 * t590 + t720 * t592;
t588 = mrSges(2,1) * g(3) + mrSges(2,3) * t708 + t723 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t599 - t752;
t587 = -mrSges(2,2) * g(3) - mrSges(2,3) * t707 + Ifges(2,5) * qJDD(1) - t723 * Ifges(2,6) - pkin(6) * t599 + t720 * t590 - t718 * t592;
t1 = [-m(1) * g(1) + t734; -m(1) * g(2) + t744; (-m(1) - m(2)) * g(3) + t599; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t744 + t721 * t587 - t719 * t588; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t734 + t719 * t587 + t721 * t588; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t728; t728; t752; t621; -t751; t623;];
tauJB = t1;
