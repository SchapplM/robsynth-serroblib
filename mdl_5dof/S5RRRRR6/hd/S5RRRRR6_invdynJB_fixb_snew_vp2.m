% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:28
% EndTime: 2019-12-05 19:00:34
% DurationCPUTime: 6.34s
% Computational Cost: add. (117086->275), mult. (150157->346), div. (0->0), fcn. (97178->10), ass. (0->115)
t680 = qJDD(1) + qJDD(2);
t687 = sin(qJ(3));
t692 = cos(qJ(3));
t682 = qJD(1) + qJD(2);
t712 = qJD(3) * t682;
t657 = t687 * t680 + t692 * t712;
t689 = sin(qJ(1));
t694 = cos(qJ(1));
t669 = t694 * g(2) + t689 * g(3);
t662 = qJDD(1) * pkin(1) + t669;
t668 = t689 * g(2) - t694 * g(3);
t695 = qJD(1) ^ 2;
t663 = -t695 * pkin(1) + t668;
t688 = sin(qJ(2));
t693 = cos(qJ(2));
t643 = t688 * t662 + t693 * t663;
t678 = t682 ^ 2;
t640 = -t678 * pkin(2) + t680 * pkin(7) + t643;
t713 = t687 * t640;
t716 = pkin(3) * t678;
t618 = qJDD(3) * pkin(3) - t657 * pkin(8) - t713 + (pkin(8) * t712 + t687 * t716 - g(1)) * t692;
t630 = -t687 * g(1) + t692 * t640;
t658 = t692 * t680 - t687 * t712;
t715 = t682 * t687;
t666 = qJD(3) * pkin(3) - pkin(8) * t715;
t684 = t692 ^ 2;
t619 = t658 * pkin(8) - qJD(3) * t666 - t684 * t716 + t630;
t686 = sin(qJ(4));
t691 = cos(qJ(4));
t600 = t691 * t618 - t686 * t619;
t651 = (-t686 * t687 + t691 * t692) * t682;
t626 = t651 * qJD(4) + t691 * t657 + t686 * t658;
t652 = (t686 * t692 + t687 * t691) * t682;
t679 = qJDD(3) + qJDD(4);
t681 = qJD(3) + qJD(4);
t595 = (t651 * t681 - t626) * pkin(9) + (t651 * t652 + t679) * pkin(4) + t600;
t601 = t686 * t618 + t691 * t619;
t625 = -t652 * qJD(4) - t686 * t657 + t691 * t658;
t646 = t681 * pkin(4) - t652 * pkin(9);
t647 = t651 ^ 2;
t596 = -t647 * pkin(4) + t625 * pkin(9) - t681 * t646 + t601;
t685 = sin(qJ(5));
t690 = cos(qJ(5));
t593 = t690 * t595 - t685 * t596;
t635 = t690 * t651 - t685 * t652;
t607 = t635 * qJD(5) + t685 * t625 + t690 * t626;
t636 = t685 * t651 + t690 * t652;
t614 = -t635 * mrSges(6,1) + t636 * mrSges(6,2);
t674 = qJD(5) + t681;
t627 = -t674 * mrSges(6,2) + t635 * mrSges(6,3);
t673 = qJDD(5) + t679;
t590 = m(6) * t593 + t673 * mrSges(6,1) - t607 * mrSges(6,3) - t636 * t614 + t674 * t627;
t594 = t685 * t595 + t690 * t596;
t606 = -t636 * qJD(5) + t690 * t625 - t685 * t626;
t628 = t674 * mrSges(6,1) - t636 * mrSges(6,3);
t591 = m(6) * t594 - t673 * mrSges(6,2) + t606 * mrSges(6,3) + t635 * t614 - t674 * t628;
t581 = t690 * t590 + t685 * t591;
t638 = -t651 * mrSges(5,1) + t652 * mrSges(5,2);
t644 = -t681 * mrSges(5,2) + t651 * mrSges(5,3);
t578 = m(5) * t600 + t679 * mrSges(5,1) - t626 * mrSges(5,3) - t652 * t638 + t681 * t644 + t581;
t645 = t681 * mrSges(5,1) - t652 * mrSges(5,3);
t707 = -t685 * t590 + t690 * t591;
t579 = m(5) * t601 - t679 * mrSges(5,2) + t625 * mrSges(5,3) + t651 * t638 - t681 * t645 + t707;
t574 = t691 * t578 + t686 * t579;
t629 = -t692 * g(1) - t713;
t649 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t687 + Ifges(4,2) * t692) * t682;
t650 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t687 + Ifges(4,4) * t692) * t682;
t632 = Ifges(5,4) * t652 + Ifges(5,2) * t651 + Ifges(5,6) * t681;
t633 = Ifges(5,1) * t652 + Ifges(5,4) * t651 + Ifges(5,5) * t681;
t610 = Ifges(6,4) * t636 + Ifges(6,2) * t635 + Ifges(6,6) * t674;
t611 = Ifges(6,1) * t636 + Ifges(6,4) * t635 + Ifges(6,5) * t674;
t701 = -mrSges(6,1) * t593 + mrSges(6,2) * t594 - Ifges(6,5) * t607 - Ifges(6,6) * t606 - Ifges(6,3) * t673 - t636 * t610 + t635 * t611;
t698 = -mrSges(5,1) * t600 + mrSges(5,2) * t601 - Ifges(5,5) * t626 - Ifges(5,6) * t625 - Ifges(5,3) * t679 - pkin(4) * t581 - t652 * t632 + t651 * t633 + t701;
t717 = mrSges(4,1) * t629 - mrSges(4,2) * t630 + Ifges(4,5) * t657 + Ifges(4,6) * t658 + Ifges(4,3) * qJDD(3) + pkin(3) * t574 + (t687 * t649 - t692 * t650) * t682 - t698;
t714 = t682 * t692;
t656 = (-mrSges(4,1) * t692 + mrSges(4,2) * t687) * t682;
t665 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t714;
t572 = m(4) * t629 + qJDD(3) * mrSges(4,1) - t657 * mrSges(4,3) + qJD(3) * t665 - t656 * t715 + t574;
t664 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t715;
t708 = -t686 * t578 + t691 * t579;
t573 = m(4) * t630 - qJDD(3) * mrSges(4,2) + t658 * mrSges(4,3) - qJD(3) * t664 + t656 * t714 + t708;
t709 = -t687 * t572 + t692 * t573;
t564 = m(3) * t643 - t678 * mrSges(3,1) - t680 * mrSges(3,2) + t709;
t642 = t693 * t662 - t688 * t663;
t703 = -t680 * pkin(2) - t642;
t639 = -t678 * pkin(7) + t703;
t620 = -t658 * pkin(3) + t666 * t715 + (-pkin(8) * t684 - pkin(7)) * t678 + t703;
t598 = -t625 * pkin(4) - t647 * pkin(9) + t652 * t646 + t620;
t706 = m(6) * t598 - t606 * mrSges(6,1) + t607 * mrSges(6,2) - t635 * t627 + t636 * t628;
t699 = m(5) * t620 - t625 * mrSges(5,1) + t626 * mrSges(5,2) - t651 * t644 + t652 * t645 + t706;
t697 = -m(4) * t639 + t658 * mrSges(4,1) - t657 * mrSges(4,2) - t664 * t715 + t665 * t714 - t699;
t585 = m(3) * t642 + t680 * mrSges(3,1) - t678 * mrSges(3,2) + t697;
t561 = t688 * t564 + t693 * t585;
t566 = t692 * t572 + t687 * t573;
t710 = t693 * t564 - t688 * t585;
t558 = m(2) * t668 - t695 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t710;
t559 = m(2) * t669 + qJDD(1) * mrSges(2,1) - t695 * mrSges(2,2) + t561;
t711 = t694 * t558 - t689 * t559;
t705 = -t689 * t558 - t694 * t559;
t609 = Ifges(6,5) * t636 + Ifges(6,6) * t635 + Ifges(6,3) * t674;
t582 = -mrSges(6,1) * t598 + mrSges(6,3) * t594 + Ifges(6,4) * t607 + Ifges(6,2) * t606 + Ifges(6,6) * t673 - t636 * t609 + t674 * t611;
t583 = mrSges(6,2) * t598 - mrSges(6,3) * t593 + Ifges(6,1) * t607 + Ifges(6,4) * t606 + Ifges(6,5) * t673 + t635 * t609 - t674 * t610;
t631 = Ifges(5,5) * t652 + Ifges(5,6) * t651 + Ifges(5,3) * t681;
t567 = -mrSges(5,1) * t620 + mrSges(5,3) * t601 + Ifges(5,4) * t626 + Ifges(5,2) * t625 + Ifges(5,6) * t679 - pkin(4) * t706 + pkin(9) * t707 + t690 * t582 + t685 * t583 - t652 * t631 + t681 * t633;
t568 = mrSges(5,2) * t620 - mrSges(5,3) * t600 + Ifges(5,1) * t626 + Ifges(5,4) * t625 + Ifges(5,5) * t679 - pkin(9) * t581 - t685 * t582 + t690 * t583 + t651 * t631 - t681 * t632;
t648 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t687 + Ifges(4,6) * t692) * t682;
t554 = -mrSges(4,1) * t639 + mrSges(4,3) * t630 + Ifges(4,4) * t657 + Ifges(4,2) * t658 + Ifges(4,6) * qJDD(3) - pkin(3) * t699 + pkin(8) * t708 + qJD(3) * t650 + t691 * t567 + t686 * t568 - t648 * t715;
t556 = mrSges(4,2) * t639 - mrSges(4,3) * t629 + Ifges(4,1) * t657 + Ifges(4,4) * t658 + Ifges(4,5) * qJDD(3) - pkin(8) * t574 - qJD(3) * t649 - t686 * t567 + t691 * t568 + t648 * t714;
t702 = mrSges(3,1) * t642 - mrSges(3,2) * t643 + Ifges(3,3) * t680 + pkin(2) * t697 + pkin(7) * t709 + t692 * t554 + t687 * t556;
t700 = mrSges(2,1) * t669 - mrSges(2,2) * t668 + Ifges(2,3) * qJDD(1) + pkin(1) * t561 + t702;
t552 = mrSges(3,1) * g(1) + mrSges(3,3) * t643 + t678 * Ifges(3,5) + Ifges(3,6) * t680 - pkin(2) * t566 - t717;
t551 = -mrSges(3,2) * g(1) - mrSges(3,3) * t642 + Ifges(3,5) * t680 - t678 * Ifges(3,6) - pkin(7) * t566 - t687 * t554 + t692 * t556;
t550 = -mrSges(2,2) * g(1) - mrSges(2,3) * t669 + Ifges(2,5) * qJDD(1) - t695 * Ifges(2,6) - pkin(6) * t561 + t693 * t551 - t688 * t552;
t549 = Ifges(2,6) * qJDD(1) + t695 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t668 + t688 * t551 + t693 * t552 - pkin(1) * (-m(3) * g(1) + t566) + pkin(6) * t710;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t566; -m(1) * g(2) + t705; -m(1) * g(3) + t711; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t700; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t711 - t694 * t549 - t689 * t550; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t705 - t689 * t549 + t694 * t550; t700; t702; t717; -t698; -t701;];
tauJB = t1;
