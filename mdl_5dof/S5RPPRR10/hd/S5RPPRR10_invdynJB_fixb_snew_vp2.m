% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRR10
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR10_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR10_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR10_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:13
% EndTime: 2019-12-31 18:04:16
% DurationCPUTime: 3.17s
% Computational Cost: add. (29148->251), mult. (70171->314), div. (0->0), fcn. (46705->8), ass. (0->109)
t753 = sin(qJ(1));
t756 = cos(qJ(1));
t722 = -g(1) * t756 - g(2) * t753;
t757 = qJD(1) ^ 2;
t795 = -pkin(1) * t757 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t722;
t749 = sin(pkin(8));
t750 = cos(pkin(8));
t794 = (Ifges(3,6) - Ifges(4,6)) * t750 + (Ifges(4,4) + Ifges(3,5)) * t749;
t742 = t749 ^ 2;
t743 = t750 ^ 2;
t784 = t743 * t757;
t791 = t742 * t757 + t784;
t788 = Ifges(3,4) - Ifges(4,5);
t790 = t788 * t749;
t721 = t753 * g(1) - t756 * g(2);
t710 = -qJDD(1) * pkin(1) - t757 * qJ(2) + qJDD(2) - t721;
t775 = qJDD(1) * t750;
t776 = qJDD(1) * t749;
t779 = t749 * qJD(1);
t697 = -pkin(2) * t775 - qJ(3) * t776 - 0.2e1 * qJD(3) * t779 + t710;
t700 = -t750 * g(3) - t749 * t795;
t789 = Ifges(3,1) + Ifges(4,1);
t787 = Ifges(3,2) + Ifges(4,3);
t786 = mrSges(3,2) * t749;
t716 = (-mrSges(4,1) * t750 - mrSges(4,3) * t749) * qJD(1);
t717 = (-mrSges(3,1) * t750 + t786) * qJD(1);
t715 = (-pkin(2) * t750 - qJ(3) * t749) * qJD(1);
t680 = t715 * t779 + qJDD(3) - t700;
t674 = (-pkin(3) * t750 * t757 - pkin(6) * qJDD(1)) * t749 + t680;
t701 = -t749 * g(3) + t750 * t795;
t778 = t750 * qJD(1);
t682 = t715 * t778 + t701;
t676 = -pkin(3) * t784 - pkin(6) * t775 + t682;
t752 = sin(qJ(4));
t755 = cos(qJ(4));
t657 = t674 * t755 - t752 * t676;
t767 = t749 * t755 - t750 * t752;
t766 = -t749 * t752 - t750 * t755;
t712 = t766 * qJD(1);
t780 = t712 * qJD(4);
t699 = qJDD(1) * t767 + t780;
t713 = t767 * qJD(1);
t652 = (-t699 + t780) * pkin(7) + (t712 * t713 + qJDD(4)) * pkin(4) + t657;
t658 = t674 * t752 + t676 * t755;
t698 = -t713 * qJD(4) + qJDD(1) * t766;
t704 = qJD(4) * pkin(4) - pkin(7) * t713;
t711 = t712 ^ 2;
t653 = -pkin(4) * t711 + pkin(7) * t698 - qJD(4) * t704 + t658;
t751 = sin(qJ(5));
t754 = cos(qJ(5));
t650 = t652 * t754 - t653 * t751;
t688 = t712 * t754 - t713 * t751;
t665 = qJD(5) * t688 + t698 * t751 + t699 * t754;
t689 = t712 * t751 + t713 * t754;
t671 = -mrSges(6,1) * t688 + mrSges(6,2) * t689;
t744 = qJD(4) + qJD(5);
t683 = -mrSges(6,2) * t744 + mrSges(6,3) * t688;
t741 = qJDD(4) + qJDD(5);
t647 = m(6) * t650 + mrSges(6,1) * t741 - mrSges(6,3) * t665 - t671 * t689 + t683 * t744;
t651 = t652 * t751 + t653 * t754;
t664 = -qJD(5) * t689 + t698 * t754 - t699 * t751;
t684 = mrSges(6,1) * t744 - mrSges(6,3) * t689;
t648 = m(6) * t651 - mrSges(6,2) * t741 + mrSges(6,3) * t664 + t671 * t688 - t684 * t744;
t636 = t647 * t754 + t648 * t751;
t692 = -mrSges(5,1) * t712 + mrSges(5,2) * t713;
t702 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t712;
t633 = m(5) * t657 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t699 + qJD(4) * t702 - t692 * t713 + t636;
t703 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t713;
t769 = -t647 * t751 + t648 * t754;
t634 = m(5) * t658 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t698 - qJD(4) * t703 + t692 * t712 + t769;
t631 = t755 * t633 + t752 * t634;
t763 = m(4) * t680 + t631;
t628 = m(3) * t700 + ((-mrSges(4,2) - mrSges(3,3)) * qJDD(1) + (-t716 - t717) * qJD(1)) * t749 - t763;
t770 = -t752 * t633 + t634 * t755;
t765 = m(4) * t682 + mrSges(4,2) * t775 + t716 * t778 + t770;
t629 = m(3) * t701 + (qJDD(1) * mrSges(3,3) + qJD(1) * t717) * t750 + t765;
t771 = -t628 * t749 + t629 * t750;
t620 = m(2) * t722 - mrSges(2,1) * t757 - qJDD(1) * mrSges(2,2) + t771;
t679 = pkin(3) * t775 + (-t742 - t743) * t757 * pkin(6) - t697;
t655 = -t698 * pkin(4) - t711 * pkin(7) + t713 * t704 + t679;
t768 = m(6) * t655 - mrSges(6,1) * t664 + mrSges(6,2) * t665 - t683 * t688 + t684 * t689;
t761 = -m(5) * t679 + t698 * mrSges(5,1) - t699 * mrSges(5,2) + t702 * t712 - t703 * t713 - t768;
t643 = m(4) * t697 - mrSges(4,1) * t775 - mrSges(4,2) * t791 - mrSges(4,3) * t776 + t761;
t758 = -m(3) * t710 + mrSges(3,1) * t775 + mrSges(3,3) * t791 - t643;
t640 = t758 + (mrSges(2,1) - t786) * qJDD(1) - t757 * mrSges(2,2) + m(2) * t721;
t783 = t620 * t753 + t640 * t756;
t622 = t628 * t750 + t629 * t749;
t781 = t794 * qJD(1);
t772 = t620 * t756 - t640 * t753;
t666 = Ifges(6,5) * t689 + Ifges(6,6) * t688 + Ifges(6,3) * t744;
t668 = Ifges(6,1) * t689 + Ifges(6,4) * t688 + Ifges(6,5) * t744;
t637 = -mrSges(6,1) * t655 + mrSges(6,3) * t651 + Ifges(6,4) * t665 + Ifges(6,2) * t664 + Ifges(6,6) * t741 - t666 * t689 + t668 * t744;
t667 = Ifges(6,4) * t689 + Ifges(6,2) * t688 + Ifges(6,6) * t744;
t638 = mrSges(6,2) * t655 - mrSges(6,3) * t650 + Ifges(6,1) * t665 + Ifges(6,4) * t664 + Ifges(6,5) * t741 + t666 * t688 - t667 * t744;
t685 = Ifges(5,5) * t713 + Ifges(5,6) * t712 + Ifges(5,3) * qJD(4);
t687 = Ifges(5,1) * t713 + Ifges(5,4) * t712 + Ifges(5,5) * qJD(4);
t623 = -mrSges(5,1) * t679 + mrSges(5,3) * t658 + Ifges(5,4) * t699 + Ifges(5,2) * t698 + Ifges(5,6) * qJDD(4) - pkin(4) * t768 + pkin(7) * t769 + qJD(4) * t687 + t754 * t637 + t751 * t638 - t713 * t685;
t686 = Ifges(5,4) * t713 + Ifges(5,2) * t712 + Ifges(5,6) * qJD(4);
t624 = mrSges(5,2) * t679 - mrSges(5,3) * t657 + Ifges(5,1) * t699 + Ifges(5,4) * t698 + Ifges(5,5) * qJDD(4) - pkin(7) * t636 - qJD(4) * t686 - t637 * t751 + t638 * t754 + t685 * t712;
t615 = -mrSges(3,1) * t710 + mrSges(3,3) * t701 - mrSges(4,1) * t697 + mrSges(4,2) * t682 - t752 * t624 - t755 * t623 - pkin(3) * t761 - pkin(6) * t770 - pkin(2) * t643 - t781 * t779 + (t750 * t787 + t790) * qJDD(1);
t617 = mrSges(3,2) * t710 + mrSges(4,2) * t680 - mrSges(3,3) * t700 - mrSges(4,3) * t697 - pkin(6) * t631 - qJ(3) * t643 - t752 * t623 + t755 * t624 + t781 * t778 + (t749 * t789 + t750 * t788) * qJDD(1);
t642 = mrSges(3,2) * t776 - t758;
t762 = mrSges(2,1) * t721 - mrSges(2,2) * t722 + Ifges(2,3) * qJDD(1) - pkin(1) * t642 + qJ(2) * t771 + t615 * t750 + t617 * t749;
t760 = mrSges(6,1) * t650 - mrSges(6,2) * t651 + Ifges(6,5) * t665 + Ifges(6,6) * t664 + Ifges(6,3) * t741 + t667 * t689 - t688 * t668;
t759 = mrSges(5,1) * t657 - mrSges(5,2) * t658 + Ifges(5,5) * t699 + Ifges(5,6) * t698 + Ifges(5,3) * qJDD(4) + pkin(4) * t636 + t686 * t713 - t712 * t687 + t760;
t630 = (qJDD(1) * mrSges(4,2) + qJD(1) * t716) * t749 + t763;
t613 = t759 + mrSges(2,1) * g(3) - pkin(1) * t622 + (Ifges(2,6) - t794) * qJDD(1) - qJ(3) * t765 + mrSges(2,3) * t722 + mrSges(3,2) * t701 - mrSges(3,1) * t700 + mrSges(4,1) * t680 - mrSges(4,3) * t682 + pkin(2) * t630 + pkin(3) * t631 + (t788 * t743 + (-t790 + (-t787 + t789) * t750) * t749 + Ifges(2,5)) * t757;
t612 = -mrSges(2,2) * g(3) - mrSges(2,3) * t721 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t757 - qJ(2) * t622 - t615 * t749 + t617 * t750;
t1 = [-m(1) * g(1) + t772; -m(1) * g(2) + t783; (-m(1) - m(2)) * g(3) + t622; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t783 + t612 * t756 - t613 * t753; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t772 + t753 * t612 + t756 * t613; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t762; t762; t642; t630; t759; t760;];
tauJB = t1;
