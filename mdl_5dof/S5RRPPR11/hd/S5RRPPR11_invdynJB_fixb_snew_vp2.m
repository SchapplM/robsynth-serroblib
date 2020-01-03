% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPPR11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR11_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR11_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR11_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:47
% EndTime: 2019-12-31 19:46:53
% DurationCPUTime: 3.83s
% Computational Cost: add. (36905->294), mult. (80581->357), div. (0->0), fcn. (45778->8), ass. (0->118)
t800 = -2 * qJD(3);
t799 = Ifges(3,1) + Ifges(4,2);
t792 = Ifges(3,4) + Ifges(4,6);
t791 = Ifges(3,5) - Ifges(4,4);
t798 = Ifges(3,2) + Ifges(4,3);
t790 = Ifges(3,6) - Ifges(4,5);
t797 = Ifges(3,3) + Ifges(4,1);
t760 = cos(qJ(2));
t757 = sin(qJ(2));
t781 = qJD(1) * qJD(2);
t779 = t757 * t781;
t732 = t760 * qJDD(1) - t779;
t782 = t757 * qJD(1);
t738 = pkin(3) * t782 - qJD(2) * qJ(4);
t753 = t760 ^ 2;
t763 = qJD(1) ^ 2;
t778 = t760 * t781;
t731 = t757 * qJDD(1) + t778;
t758 = sin(qJ(1));
t761 = cos(qJ(1));
t741 = t758 * g(1) - t761 * g(2);
t774 = -qJDD(1) * pkin(1) - t741;
t767 = pkin(2) * t779 + t782 * t800 + (-t731 - t778) * qJ(3) + t774;
t667 = -t738 * t782 + (-pkin(3) * t753 - pkin(6)) * t763 + (-pkin(2) - qJ(4)) * t732 + t767;
t742 = -t761 * g(1) - t758 * g(2);
t717 = -t763 * pkin(1) + qJDD(1) * pkin(6) + t742;
t695 = -t760 * g(3) - t757 * t717;
t728 = (-pkin(2) * t760 - qJ(3) * t757) * qJD(1);
t762 = qJD(2) ^ 2;
t685 = -qJDD(2) * pkin(2) - t762 * qJ(3) + t728 * t782 + qJDD(3) - t695;
t680 = (-t757 * t760 * t763 - qJDD(2)) * qJ(4) + (t731 - t778) * pkin(3) + t685;
t754 = sin(pkin(8));
t755 = cos(pkin(8));
t783 = qJD(1) * t760;
t723 = t755 * qJD(2) - t754 * t783;
t661 = -0.2e1 * qJD(4) * t723 - t754 * t667 + t755 * t680;
t701 = t755 * qJDD(2) - t754 * t732;
t722 = -t754 * qJD(2) - t755 * t783;
t659 = (t722 * t782 - t701) * pkin(7) + (t722 * t723 + t731) * pkin(4) + t661;
t662 = 0.2e1 * qJD(4) * t722 + t755 * t667 + t754 * t680;
t700 = -t754 * qJDD(2) - t755 * t732;
t702 = pkin(4) * t782 - t723 * pkin(7);
t721 = t722 ^ 2;
t660 = -t721 * pkin(4) + t700 * pkin(7) - t702 * t782 + t662;
t756 = sin(qJ(5));
t759 = cos(qJ(5));
t658 = t756 * t659 + t759 * t660;
t696 = -t757 * g(3) + t760 * t717;
t684 = t762 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t800 - t728 * t783 - t696;
t676 = -t753 * t763 * qJ(4) + t732 * pkin(3) + qJD(2) * t738 + qJDD(4) - t684;
t664 = -t700 * pkin(4) - t721 * pkin(7) + t723 * t702 + t676;
t693 = t756 * t722 + t759 * t723;
t671 = -t693 * qJD(5) + t759 * t700 - t756 * t701;
t692 = t759 * t722 - t756 * t723;
t672 = t692 * qJD(5) + t756 * t700 + t759 * t701;
t744 = qJD(5) + t782;
t677 = Ifges(6,5) * t693 + Ifges(6,6) * t692 + Ifges(6,3) * t744;
t679 = Ifges(6,1) * t693 + Ifges(6,4) * t692 + Ifges(6,5) * t744;
t727 = qJDD(5) + t731;
t644 = -mrSges(6,1) * t664 + mrSges(6,3) * t658 + Ifges(6,4) * t672 + Ifges(6,2) * t671 + Ifges(6,6) * t727 - t693 * t677 + t744 * t679;
t657 = t759 * t659 - t756 * t660;
t678 = Ifges(6,4) * t693 + Ifges(6,2) * t692 + Ifges(6,6) * t744;
t645 = mrSges(6,2) * t664 - mrSges(6,3) * t657 + Ifges(6,1) * t672 + Ifges(6,4) * t671 + Ifges(6,5) * t727 + t692 * t677 - t744 * t678;
t688 = Ifges(5,5) * t723 + Ifges(5,6) * t722 + Ifges(5,3) * t782;
t690 = Ifges(5,1) * t723 + Ifges(5,4) * t722 + Ifges(5,5) * t782;
t686 = -t744 * mrSges(6,2) + t692 * mrSges(6,3);
t687 = t744 * mrSges(6,1) - t693 * mrSges(6,3);
t771 = m(6) * t664 - t671 * mrSges(6,1) + t672 * mrSges(6,2) - t692 * t686 + t693 * t687;
t682 = -t692 * mrSges(6,1) + t693 * mrSges(6,2);
t652 = m(6) * t657 + t727 * mrSges(6,1) - t672 * mrSges(6,3) - t693 * t682 + t744 * t686;
t653 = m(6) * t658 - t727 * mrSges(6,2) + t671 * mrSges(6,3) + t692 * t682 - t744 * t687;
t775 = -t756 * t652 + t759 * t653;
t629 = -mrSges(5,1) * t676 + mrSges(5,3) * t662 + Ifges(5,4) * t701 + Ifges(5,2) * t700 + Ifges(5,6) * t731 - pkin(4) * t771 + pkin(7) * t775 + t759 * t644 + t756 * t645 - t723 * t688 + t690 * t782;
t643 = t759 * t652 + t756 * t653;
t689 = Ifges(5,4) * t723 + Ifges(5,2) * t722 + Ifges(5,6) * t782;
t630 = mrSges(5,2) * t676 - mrSges(5,3) * t661 + Ifges(5,1) * t701 + Ifges(5,4) * t700 + Ifges(5,5) * t731 - pkin(7) * t643 - t756 * t644 + t759 * t645 + t722 * t688 - t689 * t782;
t729 = (mrSges(4,2) * t760 - mrSges(4,3) * t757) * qJD(1);
t739 = -mrSges(4,1) * t783 - qJD(2) * mrSges(4,3);
t694 = -t722 * mrSges(5,1) + t723 * mrSges(5,2);
t698 = -mrSges(5,2) * t782 + t722 * mrSges(5,3);
t641 = m(5) * t661 + t731 * mrSges(5,1) - t701 * mrSges(5,3) - t723 * t694 + t698 * t782 + t643;
t699 = mrSges(5,1) * t782 - t723 * mrSges(5,3);
t642 = m(5) * t662 - t731 * mrSges(5,2) + t700 * mrSges(5,3) + t722 * t694 - t699 * t782 + t775;
t638 = t755 * t641 + t754 * t642;
t770 = -m(4) * t685 - t731 * mrSges(4,1) - t638;
t637 = qJDD(2) * mrSges(4,2) + qJD(2) * t739 + t729 * t782 - t770;
t655 = m(5) * t676 - t700 * mrSges(5,1) + t701 * mrSges(5,2) - t722 * t698 + t723 * t699 + t771;
t740 = mrSges(4,1) * t782 + qJD(2) * mrSges(4,2);
t766 = -m(4) * t684 + qJDD(2) * mrSges(4,3) + qJD(2) * t740 + t729 * t783 + t655;
t784 = t791 * qJD(2) + (t799 * t757 + t792 * t760) * qJD(1);
t785 = t790 * qJD(2) + (t792 * t757 + t798 * t760) * qJD(1);
t796 = (t785 * t757 - t784 * t760) * qJD(1) + t797 * qJDD(2) + t791 * t731 + t790 * t732 + mrSges(3,1) * t695 - mrSges(3,2) * t696 + mrSges(4,2) * t685 - mrSges(4,3) * t684 - pkin(2) * t637 + qJ(3) * (t732 * mrSges(4,1) + t766) - qJ(4) * t638 - t754 * t629 + t755 * t630;
t793 = t763 * pkin(6);
t730 = (-mrSges(3,1) * t760 + mrSges(3,2) * t757) * qJD(1);
t737 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t783;
t635 = m(3) * t695 - t731 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t737 - t739) * qJD(2) + (-t729 - t730) * t782 + t770;
t736 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t782;
t648 = t730 * t783 - qJD(2) * t736 + m(3) * t696 + t766 - qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t732;
t776 = -t757 * t635 + t760 * t648;
t626 = m(2) * t742 - t763 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t776;
t716 = t774 - t793;
t683 = -t732 * pkin(2) + t767 - t793;
t787 = -t754 * t641 + t755 * t642;
t773 = -m(4) * t683 - t732 * mrSges(4,2) + t740 * t782 - t787;
t765 = -m(3) * t716 + t737 * t783 + t732 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t731 + (-t736 * t757 - t739 * t760) * qJD(1) + t773;
t632 = m(2) * t741 + qJDD(1) * mrSges(2,1) - t763 * mrSges(2,2) + t765;
t788 = t758 * t626 + t761 * t632;
t628 = t760 * t635 + t757 * t648;
t786 = t797 * qJD(2) + (t791 * t757 + t790 * t760) * qJD(1);
t777 = t761 * t626 - t758 * t632;
t636 = -t731 * mrSges(4,3) + t739 * t783 - t773;
t621 = -mrSges(3,1) * t716 - mrSges(4,1) * t684 + mrSges(4,2) * t683 + mrSges(3,3) * t696 - pkin(2) * t636 + pkin(3) * t655 - qJ(4) * t787 + t784 * qJD(2) + t790 * qJDD(2) - t755 * t629 - t754 * t630 + t792 * t731 + t798 * t732 - t786 * t782;
t768 = mrSges(6,1) * t657 - mrSges(6,2) * t658 + Ifges(6,5) * t672 + Ifges(6,6) * t671 + Ifges(6,3) * t727 + t693 * t678 - t692 * t679;
t623 = t786 * t783 + t792 * t732 + (Ifges(5,3) + t799) * t731 - t785 * qJD(2) - t722 * t690 + t723 * t689 + Ifges(5,6) * t700 + Ifges(5,5) * t701 + mrSges(3,2) * t716 - mrSges(3,3) * t695 - mrSges(4,3) * t683 + mrSges(4,1) * t685 + mrSges(5,1) * t661 - mrSges(5,2) * t662 + pkin(4) * t643 + pkin(3) * t638 - qJ(3) * t636 + t768 + t791 * qJDD(2);
t769 = mrSges(2,1) * t741 - mrSges(2,2) * t742 + Ifges(2,3) * qJDD(1) + pkin(1) * t765 + pkin(6) * t776 + t760 * t621 + t757 * t623;
t619 = mrSges(2,1) * g(3) + mrSges(2,3) * t742 + t763 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t628 - t796;
t618 = -mrSges(2,2) * g(3) - mrSges(2,3) * t741 + Ifges(2,5) * qJDD(1) - t763 * Ifges(2,6) - pkin(6) * t628 - t757 * t621 + t760 * t623;
t1 = [-m(1) * g(1) + t777; -m(1) * g(2) + t788; (-m(1) - m(2)) * g(3) + t628; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t788 + t761 * t618 - t758 * t619; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t777 + t758 * t618 + t761 * t619; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t769; t769; t796; t637; t655; t768;];
tauJB = t1;
