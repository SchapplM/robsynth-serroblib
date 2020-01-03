% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRP9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRP9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP9_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP9_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP9_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:04:19
% EndTime: 2019-12-31 22:04:28
% DurationCPUTime: 5.56s
% Computational Cost: add. (58716->289), mult. (116426->351), div. (0->0), fcn. (77410->8), ass. (0->117)
t776 = Ifges(5,1) + Ifges(6,1);
t766 = Ifges(5,4) - Ifges(6,5);
t774 = Ifges(6,4) + Ifges(5,5);
t775 = Ifges(5,2) + Ifges(6,3);
t773 = Ifges(5,6) - Ifges(6,6);
t772 = -Ifges(5,3) - Ifges(6,2);
t739 = sin(qJ(3));
t742 = cos(qJ(3));
t740 = sin(qJ(2));
t762 = qJD(1) * t740;
t720 = t742 * qJD(2) - t739 * t762;
t721 = t739 * qJD(2) + t742 * t762;
t738 = sin(qJ(4));
t768 = cos(qJ(4));
t694 = -t768 * t720 + t738 * t721;
t695 = t738 * t720 + t768 * t721;
t743 = cos(qJ(2));
t761 = t743 * qJD(1);
t733 = qJD(3) - t761;
t732 = qJD(4) + t733;
t771 = t775 * t694 - t766 * t695 - t773 * t732;
t770 = -t766 * t694 + t776 * t695 + t774 * t732;
t741 = sin(qJ(1));
t744 = cos(qJ(1));
t730 = -t744 * g(1) - t741 * g(2);
t746 = qJD(1) ^ 2;
t714 = -t746 * pkin(1) + qJDD(1) * pkin(6) + t730;
t699 = -t743 * g(3) - t740 * t714;
t723 = (-pkin(2) * t743 - pkin(7) * t740) * qJD(1);
t745 = qJD(2) ^ 2;
t679 = -qJDD(2) * pkin(2) - t745 * pkin(7) + t723 * t762 - t699;
t760 = qJD(1) * qJD(2);
t758 = t743 * t760;
t724 = t740 * qJDD(1) + t758;
t692 = -t721 * qJD(3) + t742 * qJDD(2) - t739 * t724;
t701 = t733 * pkin(3) - t721 * pkin(8);
t718 = t720 ^ 2;
t647 = -t692 * pkin(3) - t718 * pkin(8) + t721 * t701 + t679;
t693 = t720 * qJD(3) + t739 * qJDD(2) + t742 * t724;
t657 = t695 * qJD(4) - t768 * t692 + t738 * t693;
t658 = -t694 * qJD(4) + t738 * t692 + t768 * t693;
t639 = -0.2e1 * qJD(5) * t695 + (t694 * t732 - t658) * qJ(5) + (t695 * t732 + t657) * pkin(4) + t647;
t681 = -t694 * mrSges(6,2) + t732 * mrSges(6,3);
t684 = -t732 * mrSges(6,1) + t695 * mrSges(6,2);
t629 = m(6) * t639 + t657 * mrSges(6,1) - t658 * mrSges(6,3) + t694 * t681 - t695 * t684;
t729 = t741 * g(1) - t744 * g(2);
t713 = -qJDD(1) * pkin(1) - t746 * pkin(6) - t729;
t734 = t740 * t760;
t725 = t743 * qJDD(1) - t734;
t675 = (-t724 - t758) * pkin(7) + (-t725 + t734) * pkin(2) + t713;
t700 = -t740 * g(3) + t743 * t714;
t680 = -t745 * pkin(2) + qJDD(2) * pkin(7) + t723 * t761 + t700;
t659 = t742 * t675 - t739 * t680;
t719 = qJDD(3) - t725;
t644 = (t720 * t733 - t693) * pkin(8) + (t720 * t721 + t719) * pkin(3) + t659;
t660 = t739 * t675 + t742 * t680;
t646 = -t718 * pkin(3) + t692 * pkin(8) - t733 * t701 + t660;
t642 = t738 * t644 + t768 * t646;
t670 = t694 * pkin(4) - t695 * qJ(5);
t715 = qJDD(4) + t719;
t731 = t732 ^ 2;
t636 = -t731 * pkin(4) + t715 * qJ(5) + 0.2e1 * qJD(5) * t732 - t694 * t670 + t642;
t764 = t773 * t694 - t774 * t695 + t772 * t732;
t615 = -mrSges(5,1) * t647 - mrSges(6,1) * t639 + mrSges(6,2) * t636 + mrSges(5,3) * t642 - pkin(4) * t629 - t775 * t657 + t766 * t658 + t764 * t695 + t773 * t715 + t770 * t732;
t641 = t768 * t644 - t738 * t646;
t637 = -t715 * pkin(4) - t731 * qJ(5) + t695 * t670 + qJDD(5) - t641;
t616 = mrSges(5,2) * t647 + mrSges(6,2) * t637 - mrSges(5,3) * t641 - mrSges(6,3) * t639 - qJ(5) * t629 - t766 * t657 + t776 * t658 + t764 * t694 + t774 * t715 + t771 * t732;
t686 = Ifges(4,5) * t721 + Ifges(4,6) * t720 + Ifges(4,3) * t733;
t688 = Ifges(4,1) * t721 + Ifges(4,4) * t720 + Ifges(4,5) * t733;
t682 = -t732 * mrSges(5,2) - t694 * mrSges(5,3);
t683 = t732 * mrSges(5,1) - t695 * mrSges(5,3);
t751 = m(5) * t647 + t657 * mrSges(5,1) + t658 * mrSges(5,2) + t694 * t682 + t695 * t683 + t629;
t759 = m(6) * t636 + t715 * mrSges(6,3) + t732 * t684;
t671 = t694 * mrSges(6,1) - t695 * mrSges(6,3);
t763 = -t694 * mrSges(5,1) - t695 * mrSges(5,2) - t671;
t767 = -mrSges(5,3) - mrSges(6,2);
t626 = m(5) * t642 - t715 * mrSges(5,2) + t767 * t657 - t732 * t683 + t763 * t694 + t759;
t754 = -m(6) * t637 + t715 * mrSges(6,1) + t732 * t681;
t628 = m(5) * t641 + t715 * mrSges(5,1) + t767 * t658 + t732 * t682 + t763 * t695 + t754;
t755 = t768 * t626 - t738 * t628;
t600 = -mrSges(4,1) * t679 + mrSges(4,3) * t660 + Ifges(4,4) * t693 + Ifges(4,2) * t692 + Ifges(4,6) * t719 - pkin(3) * t751 + pkin(8) * t755 + t768 * t615 + t738 * t616 - t721 * t686 + t733 * t688;
t620 = t738 * t626 + t768 * t628;
t687 = Ifges(4,4) * t721 + Ifges(4,2) * t720 + Ifges(4,6) * t733;
t601 = mrSges(4,2) * t679 - mrSges(4,3) * t659 + Ifges(4,1) * t693 + Ifges(4,4) * t692 + Ifges(4,5) * t719 - pkin(8) * t620 - t738 * t615 + t768 * t616 + t720 * t686 - t733 * t687;
t696 = -t720 * mrSges(4,1) + t721 * mrSges(4,2);
t697 = -t733 * mrSges(4,2) + t720 * mrSges(4,3);
t618 = m(4) * t659 + t719 * mrSges(4,1) - t693 * mrSges(4,3) - t721 * t696 + t733 * t697 + t620;
t698 = t733 * mrSges(4,1) - t721 * mrSges(4,3);
t619 = m(4) * t660 - t719 * mrSges(4,2) + t692 * mrSges(4,3) + t720 * t696 - t733 * t698 + t755;
t614 = -t739 * t618 + t742 * t619;
t623 = -m(4) * t679 + t692 * mrSges(4,1) - t693 * mrSges(4,2) + t720 * t697 - t721 * t698 - t751;
t711 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t740 + Ifges(3,2) * t743) * qJD(1);
t712 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t740 + Ifges(3,4) * t743) * qJD(1);
t769 = mrSges(3,1) * t699 - mrSges(3,2) * t700 + Ifges(3,5) * t724 + Ifges(3,6) * t725 + Ifges(3,3) * qJDD(2) + pkin(2) * t623 + pkin(7) * t614 + t742 * t600 + t739 * t601 + (t740 * t711 - t743 * t712) * qJD(1);
t722 = (-mrSges(3,1) * t743 + mrSges(3,2) * t740) * qJD(1);
t727 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t762;
t612 = m(3) * t700 - qJDD(2) * mrSges(3,2) + t725 * mrSges(3,3) - qJD(2) * t727 + t722 * t761 + t614;
t728 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t761;
t622 = m(3) * t699 + qJDD(2) * mrSges(3,1) - t724 * mrSges(3,3) + qJD(2) * t728 - t722 * t762 + t623;
t756 = t743 * t612 - t740 * t622;
t604 = m(2) * t730 - t746 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t756;
t613 = t742 * t618 + t739 * t619;
t750 = -m(3) * t713 + t725 * mrSges(3,1) - t724 * mrSges(3,2) - t727 * t762 + t728 * t761 - t613;
t608 = m(2) * t729 + qJDD(1) * mrSges(2,1) - t746 * mrSges(2,2) + t750;
t765 = t741 * t604 + t744 * t608;
t606 = t740 * t612 + t743 * t622;
t757 = t744 * t604 - t741 * t608;
t710 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t740 + Ifges(3,6) * t743) * qJD(1);
t597 = mrSges(3,2) * t713 - mrSges(3,3) * t699 + Ifges(3,1) * t724 + Ifges(3,4) * t725 + Ifges(3,5) * qJDD(2) - pkin(7) * t613 - qJD(2) * t711 - t739 * t600 + t742 * t601 + t710 * t761;
t633 = t658 * mrSges(6,2) + t695 * t671 - t754;
t749 = -mrSges(5,1) * t641 + mrSges(6,1) * t637 + mrSges(5,2) * t642 - mrSges(6,3) * t636 + pkin(4) * t633 - qJ(5) * t759 + t772 * t715 + t771 * t695 + (qJ(5) * t671 - t770) * t694 - t774 * t658 + (qJ(5) * mrSges(6,2) + t773) * t657;
t747 = mrSges(4,1) * t659 - mrSges(4,2) * t660 + Ifges(4,5) * t693 + Ifges(4,6) * t692 + Ifges(4,3) * t719 + pkin(3) * t620 + t721 * t687 - t720 * t688 - t749;
t599 = -mrSges(3,1) * t713 + mrSges(3,3) * t700 + Ifges(3,4) * t724 + Ifges(3,2) * t725 + Ifges(3,6) * qJDD(2) - pkin(2) * t613 + qJD(2) * t712 - t710 * t762 - t747;
t752 = mrSges(2,1) * t729 - mrSges(2,2) * t730 + Ifges(2,3) * qJDD(1) + pkin(1) * t750 + pkin(6) * t756 + t740 * t597 + t743 * t599;
t595 = mrSges(2,1) * g(3) + mrSges(2,3) * t730 + t746 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t606 - t769;
t594 = -mrSges(2,2) * g(3) - mrSges(2,3) * t729 + Ifges(2,5) * qJDD(1) - t746 * Ifges(2,6) - pkin(6) * t606 + t743 * t597 - t740 * t599;
t1 = [-m(1) * g(1) + t757; -m(1) * g(2) + t765; (-m(1) - m(2)) * g(3) + t606; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t765 + t744 * t594 - t741 * t595; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t757 + t741 * t594 + t744 * t595; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t752; t752; t769; t747; -t749; t633;];
tauJB = t1;
