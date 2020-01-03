% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPP8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:46
% EndTime: 2019-12-31 21:07:49
% DurationCPUTime: 2.33s
% Computational Cost: add. (20169->267), mult. (39166->307), div. (0->0), fcn. (22911->6), ass. (0->108)
t758 = Ifges(4,1) + Ifges(5,2) + Ifges(6,3);
t741 = Ifges(4,4) + Ifges(5,6) - Ifges(6,6);
t740 = Ifges(4,5) - Ifges(5,4) + Ifges(6,5);
t757 = -Ifges(4,2) - Ifges(5,3) - Ifges(6,2);
t739 = -Ifges(4,6) + Ifges(5,5) + Ifges(6,4);
t756 = Ifges(4,3) + Ifges(5,1) + Ifges(6,1);
t714 = sin(qJ(1));
t716 = cos(qJ(1));
t702 = t714 * g(1) - t716 * g(2);
t718 = qJD(1) ^ 2;
t684 = -qJDD(1) * pkin(1) - t718 * pkin(6) - t702;
t713 = sin(qJ(2));
t715 = cos(qJ(2));
t742 = qJD(1) * qJD(2);
t732 = t715 * t742;
t697 = t713 * qJDD(1) + t732;
t733 = t713 * t742;
t698 = t715 * qJDD(1) - t733;
t633 = (-t697 - t732) * pkin(7) + (-t698 + t733) * pkin(2) + t684;
t703 = -t716 * g(1) - t714 * g(2);
t685 = -t718 * pkin(1) + qJDD(1) * pkin(6) + t703;
t675 = -t713 * g(3) + t715 * t685;
t696 = (-pkin(2) * t715 - pkin(7) * t713) * qJD(1);
t717 = qJD(2) ^ 2;
t743 = t715 * qJD(1);
t637 = -t717 * pkin(2) + qJDD(2) * pkin(7) + t696 * t743 + t675;
t712 = sin(qJ(3));
t750 = cos(qJ(3));
t630 = t750 * t633 - t712 * t637;
t744 = qJD(1) * t713;
t693 = -t750 * qJD(2) + t712 * t744;
t694 = t712 * qJD(2) + t750 * t744;
t664 = t693 * pkin(3) - t694 * qJ(4);
t692 = qJDD(3) - t698;
t705 = -qJD(3) + t743;
t704 = t705 ^ 2;
t629 = -t692 * pkin(3) - t704 * qJ(4) + t694 * t664 + qJDD(4) - t630;
t659 = -t693 * qJD(3) + t712 * qJDD(2) + t750 * t697;
t666 = -t693 * mrSges(5,2) - t694 * mrSges(5,3);
t755 = -m(5) * t629 - t659 * mrSges(5,1) - t694 * t666;
t673 = t694 * mrSges(5,1) - t705 * mrSges(5,2);
t658 = t694 * qJD(3) - t750 * qJDD(2) + t712 * t697;
t674 = -t715 * g(3) - t713 * t685;
t636 = -qJDD(2) * pkin(2) - t717 * pkin(7) + t696 * t744 - t674;
t747 = t693 * t705;
t752 = -2 * qJD(4);
t722 = (-t659 - t747) * qJ(4) + t636 + (-t705 * pkin(3) + t752) * t694;
t628 = t658 * pkin(3) + t722;
t671 = t693 * mrSges(5,1) + t705 * mrSges(5,3);
t669 = t694 * pkin(4) + t705 * qJ(5);
t691 = t693 ^ 2;
t751 = 2 * qJD(5);
t624 = -t691 * pkin(4) + t693 * t751 - t694 * t669 + (pkin(3) + qJ(5)) * t658 + t722;
t670 = t694 * mrSges(6,1) + t705 * mrSges(6,3);
t672 = -t693 * mrSges(6,1) - t705 * mrSges(6,2);
t728 = m(6) * t624 - t659 * mrSges(6,2) + t658 * mrSges(6,3) - t694 * t670 + t693 * t672;
t724 = -m(5) * t628 + t658 * mrSges(5,2) + t693 * t671 - t728;
t617 = -t659 * mrSges(5,3) - t694 * t673 - t724;
t663 = -t694 * mrSges(6,2) + t693 * mrSges(6,3);
t631 = t712 * t633 + t750 * t637;
t723 = -t704 * pkin(3) + t692 * qJ(4) - t693 * t664 + t631;
t626 = -t658 * pkin(4) - t691 * qJ(5) + qJDD(5) + (t752 - t669) * t705 + t723;
t737 = m(6) * t626 + t692 * mrSges(6,2) - t705 * t670;
t621 = -t658 * mrSges(6,1) - t693 * t663 + t737;
t627 = 0.2e1 * qJD(4) * t705 - t723;
t734 = -t741 * t693 + t758 * t694 - t740 * t705;
t736 = -t739 * t693 - t740 * t694 + t756 * t705;
t596 = -mrSges(4,1) * t636 - mrSges(5,1) * t627 + mrSges(6,1) * t626 + mrSges(5,2) * t628 + mrSges(4,3) * t631 - mrSges(6,3) * t624 - pkin(3) * t617 + pkin(4) * t621 - qJ(5) * t728 + t757 * t658 + t741 * t659 - t739 * t692 + t736 * t694 - t734 * t705;
t622 = t705 * t751 + (t693 * t694 - t692) * qJ(5) + (t659 - t747) * pkin(4) + t629;
t729 = -m(6) * t622 + t692 * mrSges(6,3) - t705 * t672;
t620 = t659 * mrSges(6,1) + t694 * t663 - t729;
t735 = t757 * t693 + t741 * t694 + t739 * t705;
t602 = mrSges(5,1) * t629 + mrSges(6,1) * t622 + mrSges(4,2) * t636 - mrSges(6,2) * t624 - mrSges(4,3) * t630 - mrSges(5,3) * t628 + pkin(4) * t620 - qJ(4) * t617 - t741 * t658 + t758 * t659 + t740 * t692 + t736 * t693 + t735 * t705;
t665 = t693 * mrSges(4,1) + t694 * mrSges(4,2);
t667 = t705 * mrSges(4,2) - t693 * mrSges(4,3);
t748 = -mrSges(6,1) - mrSges(4,3);
t613 = m(4) * t630 + (-t667 + t671) * t705 + (-t663 - t665) * t694 + (mrSges(4,1) - mrSges(5,2)) * t692 + t748 * t659 + t729 + t755;
t668 = -t705 * mrSges(4,1) - t694 * mrSges(4,3);
t726 = -m(5) * t627 + t692 * mrSges(5,3) - t705 * t673 + t737;
t745 = -t663 - t666;
t616 = m(4) * t631 - t692 * mrSges(4,2) + t705 * t668 + (-t665 + t745) * t693 + (-mrSges(5,1) + t748) * t658 + t726;
t610 = -t712 * t613 + t750 * t616;
t614 = -m(4) * t636 - t658 * mrSges(4,1) - t693 * t667 + (-t668 + t673) * t694 + (-mrSges(4,2) + mrSges(5,3)) * t659 + t724;
t682 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t713 + Ifges(3,2) * t715) * qJD(1);
t683 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t713 + Ifges(3,4) * t715) * qJD(1);
t754 = mrSges(3,1) * t674 - mrSges(3,2) * t675 + Ifges(3,5) * t697 + Ifges(3,6) * t698 + Ifges(3,3) * qJDD(2) + pkin(2) * t614 + pkin(7) * t610 + (t682 * t713 - t683 * t715) * qJD(1) + t750 * t596 + t712 * t602;
t618 = t692 * mrSges(5,2) - t705 * t671 + t620 - t755;
t753 = t739 * t658 + t740 * t659 + t756 * t692 + t734 * t693 + t735 * t694 + mrSges(4,1) * t630 - mrSges(4,2) * t631 + mrSges(5,2) * t629 + mrSges(6,2) * t626 - mrSges(5,3) * t627 - mrSges(6,3) * t622 - pkin(3) * t618 + qJ(4) * (t745 * t693 + (-mrSges(5,1) - mrSges(6,1)) * t658 + t726) - qJ(5) * t620;
t695 = (-mrSges(3,1) * t715 + mrSges(3,2) * t713) * qJD(1);
t700 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t744;
t608 = m(3) * t675 - qJDD(2) * mrSges(3,2) + t698 * mrSges(3,3) - qJD(2) * t700 + t695 * t743 + t610;
t701 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t743;
t612 = m(3) * t674 + qJDD(2) * mrSges(3,1) - t697 * mrSges(3,3) + qJD(2) * t701 - t695 * t744 + t614;
t730 = t715 * t608 - t713 * t612;
t599 = m(2) * t703 - t718 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t730;
t609 = t750 * t613 + t712 * t616;
t721 = -m(3) * t684 + t698 * mrSges(3,1) - t697 * mrSges(3,2) - t700 * t744 + t701 * t743 - t609;
t604 = m(2) * t702 + qJDD(1) * mrSges(2,1) - t718 * mrSges(2,2) + t721;
t746 = t714 * t599 + t716 * t604;
t601 = t713 * t608 + t715 * t612;
t731 = t716 * t599 - t714 * t604;
t681 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t713 + Ifges(3,6) * t715) * qJD(1);
t593 = mrSges(3,2) * t684 - mrSges(3,3) * t674 + Ifges(3,1) * t697 + Ifges(3,4) * t698 + Ifges(3,5) * qJDD(2) - pkin(7) * t609 - qJD(2) * t682 - t712 * t596 + t750 * t602 + t681 * t743;
t595 = -mrSges(3,1) * t684 + mrSges(3,3) * t675 + Ifges(3,4) * t697 + Ifges(3,2) * t698 + Ifges(3,6) * qJDD(2) - pkin(2) * t609 + qJD(2) * t683 - t681 * t744 - t753;
t725 = mrSges(2,1) * t702 - mrSges(2,2) * t703 + Ifges(2,3) * qJDD(1) + pkin(1) * t721 + pkin(6) * t730 + t713 * t593 + t715 * t595;
t591 = mrSges(2,1) * g(3) + mrSges(2,3) * t703 + t718 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t601 - t754;
t590 = -mrSges(2,2) * g(3) - mrSges(2,3) * t702 + Ifges(2,5) * qJDD(1) - t718 * Ifges(2,6) - pkin(6) * t601 + t715 * t593 - t713 * t595;
t1 = [-m(1) * g(1) + t731; -m(1) * g(2) + t746; (-m(1) - m(2)) * g(3) + t601; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t746 + t716 * t590 - t714 * t591; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t731 + t714 * t590 + t716 * t591; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t725; t725; t754; t753; t618; t621;];
tauJB = t1;
