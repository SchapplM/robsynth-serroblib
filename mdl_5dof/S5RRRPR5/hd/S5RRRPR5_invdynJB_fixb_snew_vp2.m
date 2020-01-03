% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:13:37
% EndTime: 2019-12-31 21:13:45
% DurationCPUTime: 8.20s
% Computational Cost: add. (121344->315), mult. (271991->400), div. (0->0), fcn. (192957->10), ass. (0->127)
t759 = sin(qJ(2));
t763 = cos(qJ(2));
t781 = qJD(1) * qJD(2);
t739 = t759 * qJDD(1) + t763 * t781;
t760 = sin(qJ(1));
t764 = cos(qJ(1));
t746 = -t764 * g(1) - t760 * g(2);
t765 = qJD(1) ^ 2;
t734 = -t765 * pkin(1) + qJDD(1) * pkin(6) + t746;
t785 = t759 * t734;
t786 = pkin(2) * t765;
t699 = qJDD(2) * pkin(2) - t739 * pkin(7) - t785 + (pkin(7) * t781 + t759 * t786 - g(3)) * t763;
t722 = -t759 * g(3) + t763 * t734;
t740 = t763 * qJDD(1) - t759 * t781;
t783 = qJD(1) * t759;
t744 = qJD(2) * pkin(2) - pkin(7) * t783;
t754 = t763 ^ 2;
t700 = t740 * pkin(7) - qJD(2) * t744 - t754 * t786 + t722;
t758 = sin(qJ(3));
t762 = cos(qJ(3));
t676 = t762 * t699 - t758 * t700;
t731 = (-t758 * t759 + t762 * t763) * qJD(1);
t707 = t731 * qJD(3) + t762 * t739 + t758 * t740;
t732 = (t758 * t763 + t759 * t762) * qJD(1);
t751 = qJDD(2) + qJDD(3);
t752 = qJD(2) + qJD(3);
t662 = (t731 * t752 - t707) * qJ(4) + (t731 * t732 + t751) * pkin(3) + t676;
t677 = t758 * t699 + t762 * t700;
t706 = -t732 * qJD(3) - t758 * t739 + t762 * t740;
t724 = t752 * pkin(3) - t732 * qJ(4);
t727 = t731 ^ 2;
t664 = -t727 * pkin(3) + t706 * qJ(4) - t752 * t724 + t677;
t755 = sin(pkin(9));
t756 = cos(pkin(9));
t718 = t756 * t731 - t755 * t732;
t787 = 2 * qJD(4);
t659 = t755 * t662 + t756 * t664 + t718 * t787;
t683 = t756 * t706 - t755 * t707;
t719 = t755 * t731 + t756 * t732;
t693 = -t718 * mrSges(5,1) + t719 * mrSges(5,2);
t710 = t752 * mrSges(5,1) - t719 * mrSges(5,3);
t694 = -t718 * pkin(4) - t719 * pkin(8);
t750 = t752 ^ 2;
t656 = -t750 * pkin(4) + t751 * pkin(8) + t718 * t694 + t659;
t745 = t760 * g(1) - t764 * g(2);
t773 = -qJDD(1) * pkin(1) - t745;
t708 = -t740 * pkin(2) + t744 * t783 + (-pkin(7) * t754 - pkin(6)) * t765 + t773;
t669 = -t706 * pkin(3) - t727 * qJ(4) + t732 * t724 + qJDD(4) + t708;
t684 = t755 * t706 + t756 * t707;
t660 = (-t718 * t752 - t684) * pkin(8) + (t719 * t752 - t683) * pkin(4) + t669;
t757 = sin(qJ(5));
t761 = cos(qJ(5));
t653 = -t757 * t656 + t761 * t660;
t704 = -t757 * t719 + t761 * t752;
t667 = t704 * qJD(5) + t761 * t684 + t757 * t751;
t682 = qJDD(5) - t683;
t705 = t761 * t719 + t757 * t752;
t685 = -t704 * mrSges(6,1) + t705 * mrSges(6,2);
t712 = qJD(5) - t718;
t686 = -t712 * mrSges(6,2) + t704 * mrSges(6,3);
t649 = m(6) * t653 + t682 * mrSges(6,1) - t667 * mrSges(6,3) - t705 * t685 + t712 * t686;
t654 = t761 * t656 + t757 * t660;
t666 = -t705 * qJD(5) - t757 * t684 + t761 * t751;
t687 = t712 * mrSges(6,1) - t705 * mrSges(6,3);
t650 = m(6) * t654 - t682 * mrSges(6,2) + t666 * mrSges(6,3) + t704 * t685 - t712 * t687;
t776 = -t757 * t649 + t761 * t650;
t635 = m(5) * t659 - t751 * mrSges(5,2) + t683 * mrSges(5,3) + t718 * t693 - t752 * t710 + t776;
t775 = -t756 * t662 + t755 * t664;
t658 = -0.2e1 * qJD(4) * t719 - t775;
t709 = -t752 * mrSges(5,2) + t718 * mrSges(5,3);
t655 = -t751 * pkin(4) - t750 * pkin(8) + (t787 + t694) * t719 + t775;
t771 = -m(6) * t655 + t666 * mrSges(6,1) - t667 * mrSges(6,2) + t704 * t686 - t705 * t687;
t645 = m(5) * t658 + t751 * mrSges(5,1) - t684 * mrSges(5,3) - t719 * t693 + t752 * t709 + t771;
t629 = t755 * t635 + t756 * t645;
t720 = -t731 * mrSges(4,1) + t732 * mrSges(4,2);
t723 = -t752 * mrSges(4,2) + t731 * mrSges(4,3);
t626 = m(4) * t676 + t751 * mrSges(4,1) - t707 * mrSges(4,3) - t732 * t720 + t752 * t723 + t629;
t725 = t752 * mrSges(4,1) - t732 * mrSges(4,3);
t777 = t756 * t635 - t755 * t645;
t627 = m(4) * t677 - t751 * mrSges(4,2) + t706 * mrSges(4,3) + t731 * t720 - t752 * t725 + t777;
t620 = t762 * t626 + t758 * t627;
t721 = -t763 * g(3) - t785;
t729 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t759 + Ifges(3,2) * t763) * qJD(1);
t730 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t759 + Ifges(3,4) * t763) * qJD(1);
t670 = Ifges(6,5) * t705 + Ifges(6,6) * t704 + Ifges(6,3) * t712;
t672 = Ifges(6,1) * t705 + Ifges(6,4) * t704 + Ifges(6,5) * t712;
t642 = -mrSges(6,1) * t655 + mrSges(6,3) * t654 + Ifges(6,4) * t667 + Ifges(6,2) * t666 + Ifges(6,6) * t682 - t705 * t670 + t712 * t672;
t671 = Ifges(6,4) * t705 + Ifges(6,2) * t704 + Ifges(6,6) * t712;
t643 = mrSges(6,2) * t655 - mrSges(6,3) * t653 + Ifges(6,1) * t667 + Ifges(6,4) * t666 + Ifges(6,5) * t682 + t704 * t670 - t712 * t671;
t689 = Ifges(5,4) * t719 + Ifges(5,2) * t718 + Ifges(5,6) * t752;
t690 = Ifges(5,1) * t719 + Ifges(5,4) * t718 + Ifges(5,5) * t752;
t714 = Ifges(4,4) * t732 + Ifges(4,2) * t731 + Ifges(4,6) * t752;
t715 = Ifges(4,1) * t732 + Ifges(4,4) * t731 + Ifges(4,5) * t752;
t768 = -mrSges(4,1) * t676 - mrSges(5,1) * t658 + mrSges(4,2) * t677 + mrSges(5,2) * t659 - pkin(3) * t629 - pkin(4) * t771 - pkin(8) * t776 - t761 * t642 - t757 * t643 - t719 * t689 + t718 * t690 + t731 * t715 - Ifges(5,6) * t683 - Ifges(5,5) * t684 - t732 * t714 - Ifges(4,6) * t706 - Ifges(4,5) * t707 + (-Ifges(4,3) - Ifges(5,3)) * t751;
t788 = mrSges(3,1) * t721 - mrSges(3,2) * t722 + Ifges(3,5) * t739 + Ifges(3,6) * t740 + Ifges(3,3) * qJDD(2) + pkin(2) * t620 + (t759 * t729 - t763 * t730) * qJD(1) - t768;
t738 = (-mrSges(3,1) * t763 + mrSges(3,2) * t759) * qJD(1);
t782 = qJD(1) * t763;
t743 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t782;
t618 = m(3) * t721 + qJDD(2) * mrSges(3,1) - t739 * mrSges(3,3) + qJD(2) * t743 - t738 * t783 + t620;
t742 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t783;
t778 = -t758 * t626 + t762 * t627;
t619 = m(3) * t722 - qJDD(2) * mrSges(3,2) + t740 * mrSges(3,3) - qJD(2) * t742 + t738 * t782 + t778;
t779 = -t759 * t618 + t763 * t619;
t611 = m(2) * t746 - t765 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t779;
t733 = -t765 * pkin(6) + t773;
t638 = t761 * t649 + t757 * t650;
t636 = m(5) * t669 - t683 * mrSges(5,1) + t684 * mrSges(5,2) - t718 * t709 + t719 * t710 + t638;
t770 = m(4) * t708 - t706 * mrSges(4,1) + t707 * mrSges(4,2) - t731 * t723 + t732 * t725 + t636;
t767 = -m(3) * t733 + t740 * mrSges(3,1) - t739 * mrSges(3,2) - t742 * t783 + t743 * t782 - t770;
t631 = m(2) * t745 + qJDD(1) * mrSges(2,1) - t765 * mrSges(2,2) + t767;
t784 = t760 * t611 + t764 * t631;
t613 = t763 * t618 + t759 * t619;
t780 = t764 * t611 - t760 * t631;
t688 = Ifges(5,5) * t719 + Ifges(5,6) * t718 + Ifges(5,3) * t752;
t621 = mrSges(5,2) * t669 - mrSges(5,3) * t658 + Ifges(5,1) * t684 + Ifges(5,4) * t683 + Ifges(5,5) * t751 - pkin(8) * t638 - t757 * t642 + t761 * t643 + t718 * t688 - t752 * t689;
t769 = mrSges(6,1) * t653 - mrSges(6,2) * t654 + Ifges(6,5) * t667 + Ifges(6,6) * t666 + Ifges(6,3) * t682 + t705 * t671 - t704 * t672;
t622 = -mrSges(5,1) * t669 + mrSges(5,3) * t659 + Ifges(5,4) * t684 + Ifges(5,2) * t683 + Ifges(5,6) * t751 - pkin(4) * t638 - t719 * t688 + t752 * t690 - t769;
t713 = Ifges(4,5) * t732 + Ifges(4,6) * t731 + Ifges(4,3) * t752;
t608 = -mrSges(4,1) * t708 + mrSges(4,3) * t677 + Ifges(4,4) * t707 + Ifges(4,2) * t706 + Ifges(4,6) * t751 - pkin(3) * t636 + qJ(4) * t777 + t755 * t621 + t756 * t622 - t732 * t713 + t752 * t715;
t614 = mrSges(4,2) * t708 - mrSges(4,3) * t676 + Ifges(4,1) * t707 + Ifges(4,4) * t706 + Ifges(4,5) * t751 - qJ(4) * t629 + t756 * t621 - t755 * t622 + t731 * t713 - t752 * t714;
t728 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t759 + Ifges(3,6) * t763) * qJD(1);
t604 = -mrSges(3,1) * t733 + mrSges(3,3) * t722 + Ifges(3,4) * t739 + Ifges(3,2) * t740 + Ifges(3,6) * qJDD(2) - pkin(2) * t770 + pkin(7) * t778 + qJD(2) * t730 + t762 * t608 + t758 * t614 - t728 * t783;
t606 = mrSges(3,2) * t733 - mrSges(3,3) * t721 + Ifges(3,1) * t739 + Ifges(3,4) * t740 + Ifges(3,5) * qJDD(2) - pkin(7) * t620 - qJD(2) * t729 - t758 * t608 + t762 * t614 + t728 * t782;
t772 = mrSges(2,1) * t745 - mrSges(2,2) * t746 + Ifges(2,3) * qJDD(1) + pkin(1) * t767 + pkin(6) * t779 + t763 * t604 + t759 * t606;
t607 = mrSges(2,1) * g(3) + mrSges(2,3) * t746 + t765 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t613 - t788;
t602 = -mrSges(2,2) * g(3) - mrSges(2,3) * t745 + Ifges(2,5) * qJDD(1) - t765 * Ifges(2,6) - pkin(6) * t613 - t759 * t604 + t763 * t606;
t1 = [-m(1) * g(1) + t780; -m(1) * g(2) + t784; (-m(1) - m(2)) * g(3) + t613; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t784 + t764 * t602 - t760 * t607; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t780 + t760 * t602 + t764 * t607; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t772; t772; t788; -t768; t636; t769;];
tauJB = t1;
