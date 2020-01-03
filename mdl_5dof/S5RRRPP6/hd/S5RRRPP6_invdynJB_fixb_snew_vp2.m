% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:00:27
% EndTime: 2019-12-31 21:00:35
% DurationCPUTime: 5.03s
% Computational Cost: add. (54958->288), mult. (112681->352), div. (0->0), fcn. (73998->8), ass. (0->116)
t761 = Ifges(5,1) + Ifges(6,1);
t754 = Ifges(5,4) - Ifges(6,5);
t753 = Ifges(5,5) + Ifges(6,4);
t760 = -Ifges(5,2) - Ifges(6,3);
t759 = -Ifges(6,2) - Ifges(5,3);
t752 = Ifges(5,6) - Ifges(6,6);
t723 = sin(qJ(1));
t726 = cos(qJ(1));
t714 = -t726 * g(1) - t723 * g(2);
t728 = qJD(1) ^ 2;
t697 = -t728 * pkin(1) + qJDD(1) * pkin(6) + t714;
t722 = sin(qJ(2));
t725 = cos(qJ(2));
t687 = -t725 * g(3) - t722 * t697;
t707 = (-pkin(2) * t725 - pkin(7) * t722) * qJD(1);
t727 = qJD(2) ^ 2;
t745 = qJD(1) * t722;
t665 = -qJDD(2) * pkin(2) - t727 * pkin(7) + t707 * t745 - t687;
t721 = sin(qJ(3));
t724 = cos(qJ(3));
t705 = t721 * qJD(2) + t724 * t745;
t743 = qJD(1) * qJD(2);
t739 = t725 * t743;
t708 = t722 * qJDD(1) + t739;
t679 = -t705 * qJD(3) + t724 * qJDD(2) - t721 * t708;
t744 = t725 * qJD(1);
t716 = qJD(3) - t744;
t685 = t716 * pkin(3) - t705 * qJ(4);
t704 = t724 * qJD(2) - t721 * t745;
t702 = t704 ^ 2;
t639 = -t679 * pkin(3) - t702 * qJ(4) + t705 * t685 + qJDD(4) + t665;
t680 = t704 * qJD(3) + t721 * qJDD(2) + t724 * t708;
t720 = sin(pkin(8));
t751 = cos(pkin(8));
t651 = -t751 * t679 + t720 * t680;
t652 = t720 * t679 + t751 * t680;
t681 = -t751 * t704 + t720 * t705;
t682 = t720 * t704 + t751 * t705;
t634 = -0.2e1 * qJD(5) * t682 + (t681 * t716 - t652) * qJ(5) + (t682 * t716 + t651) * pkin(4) + t639;
t669 = -t716 * mrSges(6,1) + t682 * mrSges(6,2);
t670 = -t681 * mrSges(6,2) + t716 * mrSges(6,3);
t627 = m(6) * t634 + t651 * mrSges(6,1) - t652 * mrSges(6,3) - t682 * t669 + t681 * t670;
t713 = t723 * g(1) - t726 * g(2);
t696 = -qJDD(1) * pkin(1) - t728 * pkin(6) - t713;
t740 = t722 * t743;
t709 = t725 * qJDD(1) - t740;
t661 = (-t708 - t739) * pkin(7) + (-t709 + t740) * pkin(2) + t696;
t688 = -t722 * g(3) + t725 * t697;
t666 = -t727 * pkin(2) + qJDD(2) * pkin(7) + t707 * t744 + t688;
t640 = t724 * t661 - t721 * t666;
t703 = qJDD(3) - t709;
t636 = (t704 * t716 - t680) * qJ(4) + (t704 * t705 + t703) * pkin(3) + t640;
t641 = t721 * t661 + t724 * t666;
t638 = -t702 * pkin(3) + t679 * qJ(4) - t716 * t685 + t641;
t756 = -2 * qJD(4);
t632 = t720 * t636 + t751 * t638 + t681 * t756;
t656 = t681 * pkin(4) - t682 * qJ(5);
t715 = t716 ^ 2;
t629 = -t715 * pkin(4) + t703 * qJ(5) + 0.2e1 * qJD(5) * t716 - t681 * t656 + t632;
t747 = -t754 * t681 + t761 * t682 + t753 * t716;
t749 = t752 * t681 - t753 * t682 + t759 * t716;
t613 = -mrSges(5,1) * t639 - mrSges(6,1) * t634 + mrSges(6,2) * t629 + mrSges(5,3) * t632 - pkin(4) * t627 + t760 * t651 + t754 * t652 + t749 * t682 + t752 * t703 + t747 * t716;
t733 = t751 * t636 - t720 * t638;
t630 = -t703 * pkin(4) - t715 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t656) * t682 - t733;
t631 = t682 * t756 + t733;
t748 = t760 * t681 + t754 * t682 + t752 * t716;
t614 = mrSges(5,2) * t639 + mrSges(6,2) * t630 - mrSges(5,3) * t631 - mrSges(6,3) * t634 - qJ(5) * t627 - t754 * t651 + t761 * t652 + t749 * t681 + t753 * t703 - t748 * t716;
t667 = -t716 * mrSges(5,2) - t681 * mrSges(5,3);
t668 = t716 * mrSges(5,1) - t682 * mrSges(5,3);
t624 = m(5) * t639 + t651 * mrSges(5,1) + t652 * mrSges(5,2) + t681 * t667 + t682 * t668 + t627;
t672 = Ifges(4,5) * t705 + Ifges(4,6) * t704 + Ifges(4,3) * t716;
t674 = Ifges(4,1) * t705 + Ifges(4,4) * t704 + Ifges(4,5) * t716;
t741 = m(6) * t629 + t703 * mrSges(6,3) + t716 * t669;
t657 = t681 * mrSges(6,1) - t682 * mrSges(6,3);
t746 = -t681 * mrSges(5,1) - t682 * mrSges(5,2) - t657;
t755 = -mrSges(5,3) - mrSges(6,2);
t618 = m(5) * t632 - t703 * mrSges(5,2) + t755 * t651 - t716 * t668 + t746 * t681 + t741;
t735 = -m(6) * t630 + t703 * mrSges(6,1) + t716 * t670;
t620 = m(5) * t631 + t703 * mrSges(5,1) + t755 * t652 + t716 * t667 + t746 * t682 + t735;
t736 = t751 * t618 - t720 * t620;
t595 = -mrSges(4,1) * t665 + mrSges(4,3) * t641 + Ifges(4,4) * t680 + Ifges(4,2) * t679 + Ifges(4,6) * t703 - pkin(3) * t624 + qJ(4) * t736 + t751 * t613 + t720 * t614 - t705 * t672 + t716 * t674;
t615 = t720 * t618 + t751 * t620;
t673 = Ifges(4,4) * t705 + Ifges(4,2) * t704 + Ifges(4,6) * t716;
t596 = mrSges(4,2) * t665 - mrSges(4,3) * t640 + Ifges(4,1) * t680 + Ifges(4,4) * t679 + Ifges(4,5) * t703 - qJ(4) * t615 - t720 * t613 + t751 * t614 + t704 * t672 - t716 * t673;
t683 = -t704 * mrSges(4,1) + t705 * mrSges(4,2);
t684 = -t716 * mrSges(4,2) + t704 * mrSges(4,3);
t611 = m(4) * t640 + t703 * mrSges(4,1) - t680 * mrSges(4,3) - t705 * t683 + t716 * t684 + t615;
t686 = t716 * mrSges(4,1) - t705 * mrSges(4,3);
t612 = m(4) * t641 - t703 * mrSges(4,2) + t679 * mrSges(4,3) + t704 * t683 - t716 * t686 + t736;
t609 = -t721 * t611 + t724 * t612;
t623 = -m(4) * t665 + t679 * mrSges(4,1) - t680 * mrSges(4,2) + t704 * t684 - t705 * t686 - t624;
t694 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t722 + Ifges(3,2) * t725) * qJD(1);
t695 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t722 + Ifges(3,4) * t725) * qJD(1);
t758 = mrSges(3,1) * t687 - mrSges(3,2) * t688 + Ifges(3,5) * t708 + Ifges(3,6) * t709 + Ifges(3,3) * qJDD(2) + pkin(2) * t623 + pkin(7) * t609 + t724 * t595 + t721 * t596 + (t722 * t694 - t725 * t695) * qJD(1);
t626 = t652 * mrSges(6,2) + t682 * t657 - t735;
t757 = -t752 * t651 + t753 * t652 + t747 * t681 + t748 * t682 + (Ifges(4,3) - t759) * t703 + mrSges(4,1) * t640 + mrSges(5,1) * t631 - mrSges(6,1) * t630 - mrSges(4,2) * t641 - mrSges(5,2) * t632 + mrSges(6,3) * t629 + Ifges(4,5) * t680 + Ifges(4,6) * t679 + pkin(3) * t615 - pkin(4) * t626 + qJ(5) * (-t651 * mrSges(6,2) - t681 * t657 + t741) + t705 * t673 - t704 * t674;
t706 = (-mrSges(3,1) * t725 + mrSges(3,2) * t722) * qJD(1);
t711 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t745;
t607 = m(3) * t688 - qJDD(2) * mrSges(3,2) + t709 * mrSges(3,3) - qJD(2) * t711 + t706 * t744 + t609;
t712 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t744;
t622 = m(3) * t687 + qJDD(2) * mrSges(3,1) - t708 * mrSges(3,3) + qJD(2) * t712 - t706 * t745 + t623;
t737 = t725 * t607 - t722 * t622;
t599 = m(2) * t714 - t728 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t737;
t608 = t724 * t611 + t721 * t612;
t731 = -m(3) * t696 + t709 * mrSges(3,1) - t708 * mrSges(3,2) - t711 * t745 + t712 * t744 - t608;
t603 = m(2) * t713 + qJDD(1) * mrSges(2,1) - t728 * mrSges(2,2) + t731;
t750 = t723 * t599 + t726 * t603;
t601 = t722 * t607 + t725 * t622;
t738 = t726 * t599 - t723 * t603;
t693 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t722 + Ifges(3,6) * t725) * qJD(1);
t592 = mrSges(3,2) * t696 - mrSges(3,3) * t687 + Ifges(3,1) * t708 + Ifges(3,4) * t709 + Ifges(3,5) * qJDD(2) - pkin(7) * t608 - qJD(2) * t694 - t721 * t595 + t724 * t596 + t693 * t744;
t594 = -mrSges(3,1) * t696 + mrSges(3,3) * t688 + Ifges(3,4) * t708 + Ifges(3,2) * t709 + Ifges(3,6) * qJDD(2) - pkin(2) * t608 + qJD(2) * t695 - t693 * t745 - t757;
t732 = mrSges(2,1) * t713 - mrSges(2,2) * t714 + Ifges(2,3) * qJDD(1) + pkin(1) * t731 + pkin(6) * t737 + t722 * t592 + t725 * t594;
t590 = mrSges(2,1) * g(3) + mrSges(2,3) * t714 + t728 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t601 - t758;
t589 = -mrSges(2,2) * g(3) - mrSges(2,3) * t713 + Ifges(2,5) * qJDD(1) - t728 * Ifges(2,6) - pkin(6) * t601 + t725 * t592 - t722 * t594;
t1 = [-m(1) * g(1) + t738; -m(1) * g(2) + t750; (-m(1) - m(2)) * g(3) + t601; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t750 + t726 * t589 - t723 * t590; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t738 + t723 * t589 + t726 * t590; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t732; t732; t758; t757; t624; t626;];
tauJB = t1;
