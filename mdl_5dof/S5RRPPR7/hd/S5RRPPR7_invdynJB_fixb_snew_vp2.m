% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:35:20
% EndTime: 2019-12-31 19:35:26
% DurationCPUTime: 3.91s
% Computational Cost: add. (35262->296), mult. (81123->357), div. (0->0), fcn. (51335->8), ass. (0->121)
t800 = -2 * qJD(3);
t799 = Ifges(4,1) + Ifges(5,2);
t798 = -Ifges(5,1) - Ifges(4,3);
t793 = Ifges(4,4) + Ifges(5,6);
t792 = Ifges(4,5) - Ifges(5,4);
t797 = -Ifges(4,2) - Ifges(5,3);
t791 = Ifges(4,6) - Ifges(5,5);
t755 = sin(qJ(2));
t758 = cos(qJ(2));
t778 = qJD(1) * qJD(2);
t740 = qJDD(1) * t755 + t758 * t778;
t756 = sin(qJ(1));
t759 = cos(qJ(1));
t747 = -g(1) * t759 - g(2) * t756;
t761 = qJD(1) ^ 2;
t735 = -pkin(1) * t761 + qJDD(1) * pkin(6) + t747;
t789 = t755 * t735;
t794 = pkin(2) * t761;
t683 = qJDD(2) * pkin(2) - t740 * qJ(3) - t789 + (qJ(3) * t778 + t755 * t794 - g(3)) * t758;
t715 = -g(3) * t755 + t758 * t735;
t741 = qJDD(1) * t758 - t755 * t778;
t782 = qJD(1) * t755;
t743 = qJD(2) * pkin(2) - qJ(3) * t782;
t752 = t758 ^ 2;
t684 = qJ(3) * t741 - qJD(2) * t743 - t752 * t794 + t715;
t753 = sin(pkin(8));
t790 = cos(pkin(8));
t729 = (t753 * t758 + t755 * t790) * qJD(1);
t668 = t683 * t790 - t684 * t753 + t729 * t800;
t781 = qJD(1) * t758;
t728 = t753 * t782 - t781 * t790;
t701 = mrSges(4,1) * t728 + mrSges(4,2) * t729;
t710 = t740 * t790 + t741 * t753;
t716 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t728;
t718 = mrSges(5,1) * t728 - qJD(2) * mrSges(5,3);
t700 = pkin(3) * t728 - qJ(4) * t729;
t760 = qJD(2) ^ 2;
t665 = -qJDD(2) * pkin(3) - t760 * qJ(4) + t729 * t700 + qJDD(4) - t668;
t780 = qJD(2) * t728;
t660 = (t728 * t729 - qJDD(2)) * pkin(7) + (t710 + t780) * pkin(4) + t665;
t709 = t740 * t753 - t741 * t790;
t720 = pkin(4) * t729 - qJD(2) * pkin(7);
t727 = t728 ^ 2;
t746 = t756 * g(1) - t759 * g(2);
t771 = -qJDD(1) * pkin(1) - t746;
t688 = -t741 * pkin(2) + qJDD(3) + t743 * t782 + (-qJ(3) * t752 - pkin(6)) * t761 + t771;
t795 = -2 * qJD(4);
t764 = (-t710 + t780) * qJ(4) + t688 + (pkin(3) * qJD(2) + t795) * t729;
t663 = -t727 * pkin(4) - t729 * t720 + (pkin(3) + pkin(7)) * t709 + t764;
t754 = sin(qJ(5));
t757 = cos(qJ(5));
t658 = t660 * t757 - t663 * t754;
t711 = -qJD(2) * t754 + t728 * t757;
t679 = t711 * qJD(5) + qJDD(2) * t757 + t709 * t754;
t712 = qJD(2) * t757 + t728 * t754;
t685 = -mrSges(6,1) * t711 + mrSges(6,2) * t712;
t726 = qJD(5) + t729;
t689 = -mrSges(6,2) * t726 + mrSges(6,3) * t711;
t708 = qJDD(5) + t710;
t655 = m(6) * t658 + mrSges(6,1) * t708 - mrSges(6,3) * t679 - t685 * t712 + t689 * t726;
t659 = t660 * t754 + t663 * t757;
t678 = -t712 * qJD(5) - qJDD(2) * t754 + t709 * t757;
t690 = mrSges(6,1) * t726 - mrSges(6,3) * t712;
t656 = m(6) * t659 - mrSges(6,2) * t708 + mrSges(6,3) * t678 + t685 * t711 - t690 * t726;
t646 = t757 * t655 + t754 * t656;
t702 = -mrSges(5,2) * t728 - mrSges(5,3) * t729;
t767 = -m(5) * t665 - t710 * mrSges(5,1) - t729 * t702 - t646;
t642 = m(4) * t668 - t710 * mrSges(4,3) - t729 * t701 + (mrSges(4,1) - mrSges(5,2)) * qJDD(2) + (t716 - t718) * qJD(2) + t767;
t723 = t728 * t800;
t786 = t683 * t753 + t684 * t790;
t669 = t723 + t786;
t717 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t729;
t770 = t760 * pkin(3) - qJDD(2) * qJ(4) - t786;
t664 = qJD(2) * t795 + ((2 * qJD(3)) + t700) * t728 + t770;
t719 = mrSges(5,1) * t729 + qJD(2) * mrSges(5,2);
t662 = -t709 * pkin(4) - t727 * pkin(7) - t728 * t700 + t723 + ((2 * qJD(4)) + t720) * qJD(2) - t770;
t769 = -m(6) * t662 + t678 * mrSges(6,1) - t679 * mrSges(6,2) + t711 * t689 - t712 * t690;
t765 = -m(5) * t664 + qJDD(2) * mrSges(5,3) + qJD(2) * t719 - t769;
t651 = m(4) * t669 - qJDD(2) * mrSges(4,2) - qJD(2) * t717 + (-t701 - t702) * t728 + (-mrSges(4,3) - mrSges(5,1)) * t709 + t765;
t637 = t642 * t790 + t651 * t753;
t645 = qJDD(2) * mrSges(5,2) + qJD(2) * t718 - t767;
t671 = Ifges(6,5) * t712 + Ifges(6,6) * t711 + Ifges(6,3) * t726;
t673 = Ifges(6,1) * t712 + Ifges(6,4) * t711 + Ifges(6,5) * t726;
t647 = -mrSges(6,1) * t662 + mrSges(6,3) * t659 + Ifges(6,4) * t679 + Ifges(6,2) * t678 + Ifges(6,6) * t708 - t671 * t712 + t673 * t726;
t672 = Ifges(6,4) * t712 + Ifges(6,2) * t711 + Ifges(6,6) * t726;
t648 = mrSges(6,2) * t662 - mrSges(6,3) * t658 + Ifges(6,1) * t679 + Ifges(6,4) * t678 + Ifges(6,5) * t708 + t671 * t711 - t672 * t726;
t714 = -t758 * g(3) - t789;
t731 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t755 + Ifges(3,2) * t758) * qJD(1);
t732 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t755 + Ifges(3,4) * t758) * qJD(1);
t783 = qJD(2) * t792 - t728 * t793 + t729 * t799;
t784 = qJD(2) * t791 + t728 * t797 + t729 * t793;
t796 = (t731 * t755 - t732 * t758) * qJD(1) + (Ifges(3,3) - t798) * qJDD(2) - t791 * t709 + t792 * t710 + t783 * t728 + t784 * t729 + mrSges(3,1) * t714 + mrSges(4,1) * t668 - mrSges(3,2) * t715 - mrSges(4,2) * t669 + mrSges(5,2) * t665 - mrSges(5,3) * t664 + Ifges(3,5) * t740 + Ifges(3,6) * t741 + pkin(2) * t637 - pkin(3) * t645 - pkin(7) * t646 + qJ(4) * (-t709 * mrSges(5,1) - t728 * t702 + t765) - t754 * t647 + t757 * t648;
t739 = (-mrSges(3,1) * t758 + mrSges(3,2) * t755) * qJD(1);
t745 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t781;
t635 = m(3) * t714 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t740 + qJD(2) * t745 - t739 * t782 + t637;
t744 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t782;
t774 = -t642 * t753 + t651 * t790;
t636 = m(3) * t715 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t741 - qJD(2) * t744 + t739 * t781 + t774;
t775 = -t635 * t755 + t636 * t758;
t628 = m(2) * t747 - mrSges(2,1) * t761 - qJDD(1) * mrSges(2,2) + t775;
t667 = t709 * pkin(3) + t764;
t787 = -t655 * t754 + t656 * t757;
t644 = m(5) * t667 - t709 * mrSges(5,2) - t710 * mrSges(5,3) - t728 * t718 - t729 * t719 + t787;
t643 = m(4) * t688 + t709 * mrSges(4,1) + t710 * mrSges(4,2) + t728 * t716 + t729 * t717 + t644;
t734 = -t761 * pkin(6) + t771;
t763 = -m(3) * t734 + t741 * mrSges(3,1) - mrSges(3,2) * t740 - t744 * t782 + t745 * t781 - t643;
t639 = m(2) * t746 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t761 + t763;
t788 = t628 * t756 + t639 * t759;
t630 = t635 * t758 + t636 * t755;
t785 = qJD(2) * t798 + t728 * t791 - t729 * t792;
t776 = t628 * t759 - t639 * t756;
t625 = -mrSges(4,1) * t688 - mrSges(5,1) * t664 + mrSges(5,2) * t667 + mrSges(4,3) * t669 - pkin(3) * t644 - pkin(4) * t769 - pkin(7) * t787 + t783 * qJD(2) + t791 * qJDD(2) - t757 * t647 - t754 * t648 + t709 * t797 + t793 * t710 + t785 * t729;
t766 = mrSges(6,1) * t658 - mrSges(6,2) * t659 + Ifges(6,5) * t679 + Ifges(6,6) * t678 + Ifges(6,3) * t708 + t712 * t672 - t711 * t673;
t631 = mrSges(5,1) * t665 + mrSges(4,2) * t688 - mrSges(4,3) * t668 - mrSges(5,3) * t667 + pkin(4) * t646 - qJ(4) * t644 - t784 * qJD(2) + t792 * qJDD(2) - t793 * t709 + t710 * t799 + t785 * t728 + t766;
t730 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t755 + Ifges(3,6) * t758) * qJD(1);
t621 = -mrSges(3,1) * t734 + mrSges(3,3) * t715 + Ifges(3,4) * t740 + Ifges(3,2) * t741 + Ifges(3,6) * qJDD(2) - pkin(2) * t643 + qJ(3) * t774 + qJD(2) * t732 + t625 * t790 + t753 * t631 - t730 * t782;
t624 = mrSges(3,2) * t734 - mrSges(3,3) * t714 + Ifges(3,1) * t740 + Ifges(3,4) * t741 + Ifges(3,5) * qJDD(2) - qJ(3) * t637 - qJD(2) * t731 - t625 * t753 + t631 * t790 + t730 * t781;
t768 = mrSges(2,1) * t746 - mrSges(2,2) * t747 + Ifges(2,3) * qJDD(1) + pkin(1) * t763 + pkin(6) * t775 + t621 * t758 + t624 * t755;
t622 = mrSges(2,1) * g(3) + mrSges(2,3) * t747 + t761 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t630 - t796;
t619 = -mrSges(2,2) * g(3) - mrSges(2,3) * t746 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t761 - pkin(6) * t630 - t621 * t755 + t624 * t758;
t1 = [-m(1) * g(1) + t776; -m(1) * g(2) + t788; (-m(1) - m(2)) * g(3) + t630; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t788 + t619 * t759 - t622 * t756; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t776 + t756 * t619 + t759 * t622; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t768; t768; t796; t643; t645; t766;];
tauJB = t1;
