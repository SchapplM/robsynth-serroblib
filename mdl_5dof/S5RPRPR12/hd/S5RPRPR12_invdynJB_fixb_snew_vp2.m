% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR12
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR12_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR12_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR12_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:56
% EndTime: 2019-12-31 18:30:03
% DurationCPUTime: 6.76s
% Computational Cost: add. (88309->291), mult. (216936->370), div. (0->0), fcn. (157130->10), ass. (0->127)
t742 = qJD(1) ^ 2;
t772 = cos(qJ(3));
t735 = cos(pkin(8));
t771 = pkin(2) * t735;
t733 = sin(pkin(8));
t770 = mrSges(3,2) * t733;
t730 = t735 ^ 2;
t769 = t730 * t742;
t738 = sin(qJ(1));
t740 = cos(qJ(1));
t720 = -t740 * g(1) - t738 * g(2);
t715 = -t742 * pkin(1) + qJDD(1) * qJ(2) + t720;
t764 = qJD(1) * qJD(2);
t760 = -t735 * g(3) - 0.2e1 * t733 * t764;
t683 = (-pkin(6) * qJDD(1) + t742 * t771 - t715) * t733 + t760;
t702 = -t733 * g(3) + (t715 + 0.2e1 * t764) * t735;
t762 = qJDD(1) * t735;
t689 = -pkin(2) * t769 + pkin(6) * t762 + t702;
t737 = sin(qJ(3));
t669 = t737 * t683 + t772 * t689;
t761 = t735 * t772;
t765 = t733 * qJD(1);
t713 = -qJD(1) * t761 + t737 * t765;
t749 = t772 * t733 + t735 * t737;
t714 = t749 * qJD(1);
t694 = t713 * pkin(3) - t714 * qJ(4);
t741 = qJD(3) ^ 2;
t659 = -t741 * pkin(3) + qJDD(3) * qJ(4) - t713 * t694 + t669;
t729 = t733 ^ 2;
t719 = t738 * g(1) - t740 * g(2);
t755 = qJDD(2) - t719;
t698 = (-pkin(1) - t771) * qJDD(1) + (-qJ(2) + (-t729 - t730) * pkin(6)) * t742 + t755;
t763 = qJDD(1) * t733;
t766 = t714 * qJD(3);
t699 = -qJDD(1) * t761 + t737 * t763 + t766;
t767 = t713 * qJD(3);
t700 = t749 * qJDD(1) - t767;
t662 = (-t700 + t767) * qJ(4) + (t699 + t766) * pkin(3) + t698;
t732 = sin(pkin(9));
t734 = cos(pkin(9));
t707 = t732 * qJD(3) + t734 * t714;
t651 = -0.2e1 * qJD(4) * t707 - t732 * t659 + t734 * t662;
t688 = t732 * qJDD(3) + t734 * t700;
t706 = t734 * qJD(3) - t732 * t714;
t649 = (t713 * t706 - t688) * pkin(7) + (t706 * t707 + t699) * pkin(4) + t651;
t652 = 0.2e1 * qJD(4) * t706 + t734 * t659 + t732 * t662;
t686 = t713 * pkin(4) - t707 * pkin(7);
t687 = t734 * qJDD(3) - t732 * t700;
t705 = t706 ^ 2;
t650 = -t705 * pkin(4) + t687 * pkin(7) - t713 * t686 + t652;
t736 = sin(qJ(5));
t739 = cos(qJ(5));
t647 = t739 * t649 - t736 * t650;
t676 = t739 * t706 - t736 * t707;
t658 = t676 * qJD(5) + t736 * t687 + t739 * t688;
t677 = t736 * t706 + t739 * t707;
t667 = -t676 * mrSges(6,1) + t677 * mrSges(6,2);
t711 = qJD(5) + t713;
t670 = -t711 * mrSges(6,2) + t676 * mrSges(6,3);
t697 = qJDD(5) + t699;
t644 = m(6) * t647 + t697 * mrSges(6,1) - t658 * mrSges(6,3) - t677 * t667 + t711 * t670;
t648 = t736 * t649 + t739 * t650;
t657 = -t677 * qJD(5) + t739 * t687 - t736 * t688;
t671 = t711 * mrSges(6,1) - t677 * mrSges(6,3);
t645 = m(6) * t648 - t697 * mrSges(6,2) + t657 * mrSges(6,3) + t676 * t667 - t711 * t671;
t636 = t739 * t644 + t736 * t645;
t678 = -t706 * mrSges(5,1) + t707 * mrSges(5,2);
t754 = -t713 * mrSges(5,2) + t706 * mrSges(5,3);
t634 = m(5) * t651 + t699 * mrSges(5,1) - t688 * mrSges(5,3) - t707 * t678 + t713 * t754 + t636;
t685 = t713 * mrSges(5,1) - t707 * mrSges(5,3);
t756 = -t736 * t644 + t739 * t645;
t635 = m(5) * t652 - t699 * mrSges(5,2) + t687 * mrSges(5,3) + t706 * t678 - t713 * t685 + t756;
t630 = -t732 * t634 + t734 * t635;
t695 = t713 * mrSges(4,1) + t714 * mrSges(4,2);
t709 = qJD(3) * mrSges(4,1) - t714 * mrSges(4,3);
t628 = m(4) * t669 - qJDD(3) * mrSges(4,2) - t699 * mrSges(4,3) - qJD(3) * t709 - t713 * t695 + t630;
t668 = t772 * t683 - t737 * t689;
t656 = -qJDD(3) * pkin(3) - t741 * qJ(4) + t714 * t694 + qJDD(4) - t668;
t653 = -t687 * pkin(4) - t705 * pkin(7) + t707 * t686 + t656;
t747 = m(6) * t653 - t657 * mrSges(6,1) + t658 * mrSges(6,2) - t676 * t670 + t677 * t671;
t646 = m(5) * t656 - t687 * mrSges(5,1) + t688 * mrSges(5,2) + t707 * t685 - t706 * t754 + t747;
t708 = -qJD(3) * mrSges(4,2) - t713 * mrSges(4,3);
t640 = m(4) * t668 + qJDD(3) * mrSges(4,1) - t700 * mrSges(4,3) + qJD(3) * t708 - t714 * t695 - t646;
t619 = t737 * t628 + t772 * t640;
t701 = -t733 * t715 + t760;
t750 = mrSges(3,3) * qJDD(1) + t742 * (-mrSges(3,1) * t735 + t770);
t617 = m(3) * t701 - t750 * t733 + t619;
t757 = t772 * t628 - t737 * t640;
t618 = m(3) * t702 + t750 * t735 + t757;
t758 = -t733 * t617 + t735 * t618;
t610 = m(2) * t720 - t742 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t758;
t712 = -qJDD(1) * pkin(1) - t742 * qJ(2) + t755;
t629 = t734 * t634 + t732 * t635;
t746 = m(4) * t698 + t699 * mrSges(4,1) + t700 * mrSges(4,2) + t713 * t708 + t714 * t709 + t629;
t745 = -m(3) * t712 + mrSges(3,1) * t762 - t746 + (t729 * t742 + t769) * mrSges(3,3);
t623 = t745 + (mrSges(2,1) - t770) * qJDD(1) - t742 * mrSges(2,2) + m(2) * t719;
t768 = t738 * t610 + t740 * t623;
t612 = t735 * t617 + t733 * t618;
t759 = t740 * t610 - t738 * t623;
t753 = Ifges(3,1) * t733 + Ifges(3,4) * t735;
t752 = Ifges(3,4) * t733 + Ifges(3,2) * t735;
t751 = Ifges(3,5) * t733 + Ifges(3,6) * t735;
t663 = Ifges(6,5) * t677 + Ifges(6,6) * t676 + Ifges(6,3) * t711;
t665 = Ifges(6,1) * t677 + Ifges(6,4) * t676 + Ifges(6,5) * t711;
t637 = -mrSges(6,1) * t653 + mrSges(6,3) * t648 + Ifges(6,4) * t658 + Ifges(6,2) * t657 + Ifges(6,6) * t697 - t677 * t663 + t711 * t665;
t664 = Ifges(6,4) * t677 + Ifges(6,2) * t676 + Ifges(6,6) * t711;
t638 = mrSges(6,2) * t653 - mrSges(6,3) * t647 + Ifges(6,1) * t658 + Ifges(6,4) * t657 + Ifges(6,5) * t697 + t676 * t663 - t711 * t664;
t672 = Ifges(5,5) * t707 + Ifges(5,6) * t706 + Ifges(5,3) * t713;
t674 = Ifges(5,1) * t707 + Ifges(5,4) * t706 + Ifges(5,5) * t713;
t620 = -mrSges(5,1) * t656 + mrSges(5,3) * t652 + Ifges(5,4) * t688 + Ifges(5,2) * t687 + Ifges(5,6) * t699 - pkin(4) * t747 + pkin(7) * t756 + t739 * t637 + t736 * t638 - t707 * t672 + t713 * t674;
t673 = Ifges(5,4) * t707 + Ifges(5,2) * t706 + Ifges(5,6) * t713;
t621 = mrSges(5,2) * t656 - mrSges(5,3) * t651 + Ifges(5,1) * t688 + Ifges(5,4) * t687 + Ifges(5,5) * t699 - pkin(7) * t636 - t736 * t637 + t739 * t638 + t706 * t672 - t713 * t673;
t690 = Ifges(4,5) * t714 - Ifges(4,6) * t713 + Ifges(4,3) * qJD(3);
t691 = Ifges(4,4) * t714 - Ifges(4,2) * t713 + Ifges(4,6) * qJD(3);
t607 = mrSges(4,2) * t698 - mrSges(4,3) * t668 + Ifges(4,1) * t700 - Ifges(4,4) * t699 + Ifges(4,5) * qJDD(3) - qJ(4) * t629 - qJD(3) * t691 - t732 * t620 + t734 * t621 - t713 * t690;
t692 = Ifges(4,1) * t714 - Ifges(4,4) * t713 + Ifges(4,5) * qJD(3);
t744 = mrSges(6,1) * t647 - mrSges(6,2) * t648 + Ifges(6,5) * t658 + Ifges(6,6) * t657 + Ifges(6,3) * t697 + t677 * t664 - t676 * t665;
t613 = -t744 + Ifges(4,6) * qJDD(3) + (-Ifges(4,2) - Ifges(5,3)) * t699 - t714 * t690 + Ifges(4,4) * t700 + t706 * t674 - t707 * t673 + qJD(3) * t692 - mrSges(4,1) * t698 - Ifges(5,6) * t687 - Ifges(5,5) * t688 + mrSges(4,3) * t669 + mrSges(5,2) * t652 - mrSges(5,1) * t651 - pkin(4) * t636 - pkin(3) * t629;
t717 = t751 * qJD(1);
t603 = -mrSges(3,1) * t712 + mrSges(3,3) * t702 - pkin(2) * t746 + pkin(6) * t757 + t752 * qJDD(1) + t737 * t607 + t772 * t613 - t717 * t765;
t606 = t735 * qJD(1) * t717 + mrSges(3,2) * t712 - mrSges(3,3) * t701 - pkin(6) * t619 + t753 * qJDD(1) + t772 * t607 - t737 * t613;
t625 = mrSges(3,2) * t763 - t745;
t748 = mrSges(2,1) * t719 - mrSges(2,2) * t720 + Ifges(2,3) * qJDD(1) - pkin(1) * t625 + qJ(2) * t758 + t735 * t603 + t733 * t606;
t743 = mrSges(4,1) * t668 - mrSges(4,2) * t669 + Ifges(4,5) * t700 - Ifges(4,6) * t699 + Ifges(4,3) * qJDD(3) - pkin(3) * t646 + qJ(4) * t630 + t734 * t620 + t732 * t621 + t714 * t691 + t713 * t692;
t604 = -t743 + mrSges(2,1) * g(3) + (Ifges(2,6) - t751) * qJDD(1) + mrSges(2,3) * t720 - mrSges(3,1) * t701 + mrSges(3,2) * t702 - pkin(2) * t619 - pkin(1) * t612 + (-t733 * t752 + t735 * t753 + Ifges(2,5)) * t742;
t601 = -mrSges(2,2) * g(3) - mrSges(2,3) * t719 + Ifges(2,5) * qJDD(1) - t742 * Ifges(2,6) - qJ(2) * t612 - t733 * t603 + t735 * t606;
t1 = [-m(1) * g(1) + t759; -m(1) * g(2) + t768; (-m(1) - m(2)) * g(3) + t612; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t768 + t740 * t601 - t738 * t604; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t759 + t738 * t601 + t740 * t604; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t748; t748; t625; t743; t646; t744;];
tauJB = t1;
