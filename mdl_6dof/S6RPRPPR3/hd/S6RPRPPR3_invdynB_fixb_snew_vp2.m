% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-05-05 16:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:41:42
% EndTime: 2019-05-05 16:41:47
% DurationCPUTime: 2.81s
% Computational Cost: add. (23434->307), mult. (46422->358), div. (0->0), fcn. (22077->8), ass. (0->124)
t767 = Ifges(4,1) + Ifges(5,1) + Ifges(6,2);
t746 = Ifges(4,4) - Ifges(5,5) + Ifges(6,4);
t745 = Ifges(4,5) + Ifges(5,4) + Ifges(6,6);
t766 = Ifges(4,2) + Ifges(5,3) + Ifges(6,1);
t765 = -Ifges(5,2) - Ifges(4,3) - Ifges(6,3);
t744 = Ifges(4,6) - Ifges(5,6) + Ifges(6,5);
t718 = sin(qJ(1));
t721 = cos(qJ(1));
t694 = t718 * g(1) - t721 * g(2);
t675 = qJDD(1) * pkin(1) + t694;
t695 = -t721 * g(1) - t718 * g(2);
t723 = qJD(1) ^ 2;
t680 = -t723 * pkin(1) + t695;
t713 = sin(pkin(9));
t714 = cos(pkin(9));
t638 = t713 * t675 + t714 * t680;
t629 = -t723 * pkin(2) + qJDD(1) * pkin(7) + t638;
t711 = -g(3) + qJDD(2);
t717 = sin(qJ(3));
t720 = cos(qJ(3));
t625 = t720 * t629 + t717 * t711;
t676 = (-pkin(3) * t720 - qJ(4) * t717) * qJD(1);
t751 = qJD(1) * t720;
t764 = qJDD(3) * qJ(4) + t676 * t751 + t625;
t748 = qJD(1) * qJD(3);
t738 = t720 * t748;
t682 = t717 * qJDD(1) + t738;
t747 = qJD(1) * qJD(5);
t763 = -0.2e1 * t717 * t747 + (-t682 + t738) * qJ(5);
t762 = -2 * qJD(4);
t761 = 2 * qJD(4);
t760 = -pkin(3) - pkin(8);
t759 = pkin(4) + pkin(8);
t722 = qJD(3) ^ 2;
t758 = t722 * pkin(3);
t757 = mrSges(4,3) + mrSges(5,2);
t755 = t720 ^ 2 * t723;
t754 = t720 * t723;
t624 = -t717 * t629 + t720 * t711;
t677 = (-mrSges(5,1) * t720 - mrSges(5,3) * t717) * qJD(1);
t678 = (-mrSges(4,1) * t720 + mrSges(4,2) * t717) * qJD(1);
t691 = -qJD(3) * mrSges(6,1) + mrSges(6,3) * t751;
t692 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t751;
t749 = t717 * qJD(1);
t733 = t676 * t749 + qJDD(4) - t624;
t623 = -qJDD(3) * pkin(3) - t722 * qJ(4) + t733;
t693 = mrSges(5,2) * t751 + qJD(3) * mrSges(5,3);
t619 = (-t717 * t754 - qJDD(3)) * pkin(4) + t623 + t763;
t679 = (mrSges(6,1) * t717 - mrSges(6,2) * t720) * qJD(1);
t739 = t717 * t748;
t683 = t720 * qJDD(1) - t739;
t687 = -qJD(3) * pkin(4) - qJ(5) * t749;
t637 = t714 * t675 - t713 * t680;
t628 = -qJDD(1) * pkin(2) - t723 * pkin(7) - t637;
t731 = -t683 * pkin(3) + t628 + (-t682 - t738) * qJ(4);
t726 = -qJ(5) * t755 + qJDD(5) - t731 + (t687 + t761) * t749;
t612 = t726 + (pkin(5) * t720 + t717 * t760) * t748 + t682 * pkin(5) + t759 * t683;
t681 = (pkin(5) * t717 + pkin(8) * t720) * qJD(1);
t615 = (-pkin(5) - qJ(4)) * t722 + (-pkin(4) * t754 - qJD(1) * t681) * t717 + (-pkin(3) - t759) * qJDD(3) + t733 + t763;
t716 = sin(qJ(6));
t719 = cos(qJ(6));
t610 = t719 * t612 - t716 * t615;
t673 = -t719 * qJD(3) + t716 * t751;
t636 = t673 * qJD(6) - t716 * qJDD(3) - t719 * t683;
t674 = -t716 * qJD(3) - t719 * t751;
t639 = -t673 * mrSges(7,1) + t674 * mrSges(7,2);
t698 = qJD(6) + t749;
t640 = -t698 * mrSges(7,2) + t673 * mrSges(7,3);
t671 = qJDD(6) + t682;
t608 = m(7) * t610 + t671 * mrSges(7,1) - t636 * mrSges(7,3) - t674 * t639 + t698 * t640;
t611 = t716 * t612 + t719 * t615;
t635 = -t674 * qJD(6) - t719 * qJDD(3) + t716 * t683;
t641 = t698 * mrSges(7,1) - t674 * mrSges(7,3);
t609 = m(7) * t611 - t671 * mrSges(7,2) + t635 * mrSges(7,3) + t673 * t639 - t698 * t641;
t752 = -t716 * t608 + t719 * t609;
t734 = -m(6) * t619 + t679 * t749 - t752;
t729 = -m(5) * t623 + qJDD(3) * mrSges(5,1) + qJD(3) * t693 + t734;
t595 = m(4) * t624 + (mrSges(4,1) - mrSges(6,2)) * qJDD(3) + (-t691 + t692) * qJD(3) + (-t677 - t678) * t749 + (mrSges(6,3) - t757) * t682 + t729;
t689 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t749;
t704 = qJD(3) * t761;
t622 = t704 - t758 + t764;
t690 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t749;
t728 = pkin(4) * t755 + t683 * qJ(5) - t764;
t618 = 0.2e1 * t720 * t747 + t758 + (t762 - t687) * qJD(3) + t728;
t688 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t749;
t614 = qJDD(3) * pkin(5) + qJD(3) * t687 + t704 + t760 * t722 + (-0.2e1 * qJD(5) - t681) * t751 - t728;
t730 = -m(7) * t614 + t635 * mrSges(7,1) - t636 * mrSges(7,2) + t673 * t640 - t674 * t641;
t727 = -m(6) * t618 + qJDD(3) * mrSges(6,1) - t683 * mrSges(6,3) + qJD(3) * t688 - t730;
t725 = m(5) * t622 + qJDD(3) * mrSges(5,3) + qJD(3) * t690 + t677 * t751 + t727;
t601 = t725 + (t678 - t679) * t751 - qJDD(3) * mrSges(4,2) - qJD(3) * t689 + m(4) * t625 + t757 * t683;
t735 = -t717 * t595 + t720 * t601;
t590 = m(3) * t638 - t723 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t735;
t620 = (pkin(3) * qJD(3) + t762) * t749 + t731;
t598 = t719 * t608 + t716 * t609;
t617 = -pkin(3) * t739 + t683 * pkin(4) + t726;
t732 = -m(6) * t617 - t682 * mrSges(6,1) + t683 * mrSges(6,2) - t688 * t749 + t691 * t751 - t598;
t596 = m(5) * t620 - t683 * mrSges(5,1) - t682 * mrSges(5,3) - t690 * t749 - t693 * t751 + t732;
t724 = -m(4) * t628 + t683 * mrSges(4,1) - t682 * mrSges(4,2) - t689 * t749 + t692 * t751 - t596;
t593 = m(3) * t637 + qJDD(1) * mrSges(3,1) - t723 * mrSges(3,2) + t724;
t587 = t713 * t590 + t714 * t593;
t585 = m(2) * t694 + qJDD(1) * mrSges(2,1) - t723 * mrSges(2,2) + t587;
t736 = t714 * t590 - t713 * t593;
t586 = m(2) * t695 - t723 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t736;
t753 = t721 * t585 + t718 * t586;
t591 = t720 * t595 + t717 * t601;
t750 = qJD(3) * t691;
t743 = t765 * qJD(3) + (-t717 * t745 - t720 * t744) * qJD(1);
t742 = -t744 * qJD(3) + (-t717 * t746 - t720 * t766) * qJD(1);
t741 = t745 * qJD(3) + (t717 * t767 + t746 * t720) * qJD(1);
t740 = m(3) * t711 + t591;
t737 = -t718 * t585 + t721 * t586;
t632 = Ifges(7,1) * t674 + Ifges(7,4) * t673 + Ifges(7,5) * t698;
t631 = Ifges(7,4) * t674 + Ifges(7,2) * t673 + Ifges(7,6) * t698;
t630 = Ifges(7,5) * t674 + Ifges(7,6) * t673 + Ifges(7,3) * t698;
t603 = mrSges(7,2) * t614 - mrSges(7,3) * t610 + Ifges(7,1) * t636 + Ifges(7,4) * t635 + Ifges(7,5) * t671 + t673 * t630 - t698 * t631;
t602 = -mrSges(7,1) * t614 + mrSges(7,3) * t611 + Ifges(7,4) * t636 + Ifges(7,2) * t635 + Ifges(7,6) * t671 - t674 * t630 + t698 * t632;
t597 = qJDD(3) * mrSges(6,2) - t682 * mrSges(6,3) - t734 + t750;
t581 = t674 * t631 + Ifges(7,3) * t671 - t673 * t632 + Ifges(7,5) * t636 + mrSges(5,2) * t623 - mrSges(4,3) * t624 + mrSges(4,2) * t628 + Ifges(7,6) * t635 + mrSges(6,1) * t617 - mrSges(6,3) * t619 - mrSges(5,3) * t620 + mrSges(7,1) * t610 - mrSges(7,2) * t611 + pkin(5) * t598 - qJ(4) * t596 - qJ(5) * t597 + t746 * t683 + t767 * t682 + t745 * qJDD(3) + t742 * qJD(3) - t743 * t751;
t580 = -qJ(5) * t727 + t716 * t602 - t719 * t603 - pkin(4) * t732 + mrSges(4,3) * t625 - mrSges(4,1) * t628 - mrSges(6,2) * t617 + mrSges(6,3) * t618 - mrSges(5,1) * t620 + mrSges(5,2) * t622 + pkin(8) * t598 - pkin(3) * t596 + t766 * t683 + t746 * t682 + t744 * qJDD(3) + t741 * qJD(3) + (qJ(5) * t679 * t720 + t717 * t743) * qJD(1);
t579 = Ifges(3,6) * qJDD(1) + t723 * Ifges(3,5) + t716 * t603 + t719 * t602 - mrSges(3,1) * t711 + pkin(5) * t730 + mrSges(3,3) * t638 + mrSges(5,1) * t623 - mrSges(4,1) * t624 + mrSges(4,2) * t625 + mrSges(6,1) * t618 - mrSges(6,2) * t619 - mrSges(5,3) * t622 + pkin(8) * t752 + pkin(4) * t597 - pkin(2) * t591 - qJ(4) * t725 - pkin(3) * (t729 - t750) + (-qJ(4) * mrSges(5,2) - t744) * t683 + (pkin(3) * mrSges(6,2) + t765) * qJDD(3) + (-pkin(3) * (-mrSges(5,2) + mrSges(6,3)) - t745) * t682 + ((qJ(4) * t679 + t741) * t720 + (pkin(3) * t677 + t742) * t717) * qJD(1);
t578 = mrSges(3,2) * t711 - mrSges(3,3) * t637 + Ifges(3,5) * qJDD(1) - t723 * Ifges(3,6) - pkin(7) * t591 - t717 * t580 + t720 * t581;
t577 = -mrSges(2,2) * g(3) - mrSges(2,3) * t694 + Ifges(2,5) * qJDD(1) - t723 * Ifges(2,6) - qJ(2) * t587 + t714 * t578 - t713 * t579;
t576 = mrSges(2,1) * g(3) + mrSges(2,3) * t695 + t723 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t740 + qJ(2) * t736 + t713 * t578 + t714 * t579;
t1 = [-m(1) * g(1) + t737; -m(1) * g(2) + t753; (-m(1) - m(2)) * g(3) + t740; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t753 - t718 * t576 + t721 * t577; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t737 + t721 * t576 + t718 * t577; pkin(1) * t587 + mrSges(2,1) * t694 - mrSges(2,2) * t695 + t720 * t580 + pkin(2) * t724 + pkin(7) * t735 - mrSges(3,2) * t638 + t717 * t581 + mrSges(3,1) * t637 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
