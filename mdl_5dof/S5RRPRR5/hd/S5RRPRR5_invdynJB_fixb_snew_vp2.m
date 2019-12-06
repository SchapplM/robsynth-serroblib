% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:33:57
% EndTime: 2019-12-05 18:34:03
% DurationCPUTime: 5.53s
% Computational Cost: add. (94999->251), mult. (131436->313), div. (0->0), fcn. (88562->10), ass. (0->114)
t673 = qJD(1) + qJD(2);
t667 = t673 ^ 2;
t676 = cos(pkin(9));
t713 = pkin(3) * t676;
t675 = sin(pkin(9));
t712 = mrSges(4,2) * t675;
t671 = t676 ^ 2;
t711 = t667 * t671;
t669 = qJDD(1) + qJDD(2);
t710 = t669 * t676;
t698 = Ifges(4,5) * t675 + Ifges(4,6) * t676;
t709 = t667 * t698;
t680 = sin(qJ(1));
t684 = cos(qJ(1));
t658 = t684 * g(2) + t680 * g(3);
t652 = qJDD(1) * pkin(1) + t658;
t657 = t680 * g(2) - t684 * g(3);
t685 = qJD(1) ^ 2;
t653 = -t685 * pkin(1) + t657;
t679 = sin(qJ(2));
t683 = cos(qJ(2));
t640 = t679 * t652 + t683 * t653;
t637 = -t667 * pkin(2) + t669 * qJ(3) + t640;
t708 = qJD(3) * t673;
t706 = -t676 * g(1) - 0.2e1 * t675 * t708;
t618 = (-pkin(7) * t669 + t667 * t713 - t637) * t675 + t706;
t622 = -t675 * g(1) + (t637 + 0.2e1 * t708) * t676;
t619 = -pkin(3) * t711 + pkin(7) * t710 + t622;
t678 = sin(qJ(4));
t682 = cos(qJ(4));
t600 = t682 * t618 - t678 * t619;
t693 = t675 * t682 + t676 * t678;
t692 = -t675 * t678 + t676 * t682;
t645 = t692 * t673;
t707 = t645 * qJD(4);
t636 = t693 * t669 + t707;
t646 = t693 * t673;
t596 = (-t636 + t707) * pkin(8) + (t645 * t646 + qJDD(4)) * pkin(4) + t600;
t601 = t678 * t618 + t682 * t619;
t635 = -t646 * qJD(4) + t692 * t669;
t643 = qJD(4) * pkin(4) - t646 * pkin(8);
t644 = t645 ^ 2;
t597 = -t644 * pkin(4) + t635 * pkin(8) - qJD(4) * t643 + t601;
t677 = sin(qJ(5));
t681 = cos(qJ(5));
t594 = t681 * t596 - t677 * t597;
t628 = t681 * t645 - t677 * t646;
t608 = t628 * qJD(5) + t677 * t635 + t681 * t636;
t629 = t677 * t645 + t681 * t646;
t614 = -t628 * mrSges(6,1) + t629 * mrSges(6,2);
t672 = qJD(4) + qJD(5);
t623 = -t672 * mrSges(6,2) + t628 * mrSges(6,3);
t668 = qJDD(4) + qJDD(5);
t591 = m(6) * t594 + t668 * mrSges(6,1) - t608 * mrSges(6,3) - t629 * t614 + t672 * t623;
t595 = t677 * t596 + t681 * t597;
t607 = -t629 * qJD(5) + t681 * t635 - t677 * t636;
t624 = t672 * mrSges(6,1) - t629 * mrSges(6,3);
t592 = m(6) * t595 - t668 * mrSges(6,2) + t607 * mrSges(6,3) + t628 * t614 - t672 * t624;
t581 = t681 * t591 + t677 * t592;
t633 = -t645 * mrSges(5,1) + t646 * mrSges(5,2);
t641 = -qJD(4) * mrSges(5,2) + t645 * mrSges(5,3);
t579 = m(5) * t600 + qJDD(4) * mrSges(5,1) - t636 * mrSges(5,3) + qJD(4) * t641 - t646 * t633 + t581;
t642 = qJD(4) * mrSges(5,1) - t646 * mrSges(5,3);
t701 = -t677 * t591 + t681 * t592;
t580 = m(5) * t601 - qJDD(4) * mrSges(5,2) + t635 * mrSges(5,3) - qJD(4) * t642 + t645 * t633 + t701;
t575 = t682 * t579 + t678 * t580;
t621 = -t675 * t637 + t706;
t695 = mrSges(4,3) * t669 + (-mrSges(4,1) * t676 + t712) * t667;
t573 = m(4) * t621 - t695 * t675 + t575;
t702 = -t678 * t579 + t682 * t580;
t574 = m(4) * t622 + t695 * t676 + t702;
t703 = -t675 * t573 + t676 * t574;
t565 = m(3) * t640 - t667 * mrSges(3,1) - t669 * mrSges(3,2) + t703;
t639 = t683 * t652 - t679 * t653;
t697 = qJDD(3) - t639;
t634 = -t669 * pkin(2) - t667 * qJ(3) + t697;
t670 = t675 ^ 2;
t620 = (-pkin(2) - t713) * t669 + (-qJ(3) + (-t670 - t671) * pkin(7)) * t667 + t697;
t599 = -t635 * pkin(4) - t644 * pkin(8) + t646 * t643 + t620;
t696 = m(6) * t599 - t607 * mrSges(6,1) + t608 * mrSges(6,2) - t628 * t623 + t629 * t624;
t688 = m(5) * t620 - t635 * mrSges(5,1) + t636 * mrSges(5,2) - t645 * t641 + t646 * t642 + t696;
t687 = -m(4) * t634 + mrSges(4,1) * t710 - t688 + (t667 * t670 + t711) * mrSges(4,3);
t585 = t687 - t667 * mrSges(3,2) + m(3) * t639 + (mrSges(3,1) - t712) * t669;
t562 = t679 * t565 + t683 * t585;
t567 = t676 * t573 + t675 * t574;
t704 = t683 * t565 - t679 * t585;
t559 = m(2) * t657 - t685 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t704;
t560 = m(2) * t658 + qJDD(1) * mrSges(2,1) - t685 * mrSges(2,2) + t562;
t705 = t684 * t559 - t680 * t560;
t700 = Ifges(4,1) * t675 + Ifges(4,4) * t676;
t699 = Ifges(4,4) * t675 + Ifges(4,2) * t676;
t694 = -t680 * t559 - t684 * t560;
t609 = Ifges(6,5) * t629 + Ifges(6,6) * t628 + Ifges(6,3) * t672;
t611 = Ifges(6,1) * t629 + Ifges(6,4) * t628 + Ifges(6,5) * t672;
t582 = -mrSges(6,1) * t599 + mrSges(6,3) * t595 + Ifges(6,4) * t608 + Ifges(6,2) * t607 + Ifges(6,6) * t668 - t629 * t609 + t672 * t611;
t610 = Ifges(6,4) * t629 + Ifges(6,2) * t628 + Ifges(6,6) * t672;
t583 = mrSges(6,2) * t599 - mrSges(6,3) * t594 + Ifges(6,1) * t608 + Ifges(6,4) * t607 + Ifges(6,5) * t668 + t628 * t609 - t672 * t610;
t625 = Ifges(5,5) * t646 + Ifges(5,6) * t645 + Ifges(5,3) * qJD(4);
t627 = Ifges(5,1) * t646 + Ifges(5,4) * t645 + Ifges(5,5) * qJD(4);
t568 = -mrSges(5,1) * t620 + mrSges(5,3) * t601 + Ifges(5,4) * t636 + Ifges(5,2) * t635 + Ifges(5,6) * qJDD(4) - pkin(4) * t696 + pkin(8) * t701 + qJD(4) * t627 + t681 * t582 + t677 * t583 - t646 * t625;
t626 = Ifges(5,4) * t646 + Ifges(5,2) * t645 + Ifges(5,6) * qJD(4);
t569 = mrSges(5,2) * t620 - mrSges(5,3) * t600 + Ifges(5,1) * t636 + Ifges(5,4) * t635 + Ifges(5,5) * qJDD(4) - pkin(8) * t581 - qJD(4) * t626 - t677 * t582 + t681 * t583 + t645 * t625;
t555 = -mrSges(4,1) * t634 + mrSges(4,3) * t622 - pkin(3) * t688 + pkin(7) * t702 + t682 * t568 + t678 * t569 + t699 * t669 - t675 * t709;
t557 = mrSges(4,2) * t634 - mrSges(4,3) * t621 - pkin(7) * t575 - t678 * t568 + t682 * t569 + t700 * t669 + t676 * t709;
t587 = t669 * t712 - t687;
t691 = mrSges(3,1) * t639 - mrSges(3,2) * t640 + Ifges(3,3) * t669 - pkin(2) * t587 + qJ(3) * t703 + t676 * t555 + t675 * t557;
t690 = -mrSges(6,1) * t594 + mrSges(6,2) * t595 - Ifges(6,5) * t608 - Ifges(6,6) * t607 - Ifges(6,3) * t668 - t629 * t610 + t628 * t611;
t689 = mrSges(2,1) * t658 - mrSges(2,2) * t657 + Ifges(2,3) * qJDD(1) + pkin(1) * t562 + t691;
t686 = mrSges(5,1) * t600 - mrSges(5,2) * t601 + Ifges(5,5) * t636 + Ifges(5,6) * t635 + Ifges(5,3) * qJDD(4) + pkin(4) * t581 + t646 * t626 - t645 * t627 - t690;
t553 = (Ifges(3,6) - t698) * t669 - t686 + mrSges(3,1) * g(1) + mrSges(3,3) * t640 - mrSges(4,1) * t621 + mrSges(4,2) * t622 - pkin(2) * t567 - pkin(3) * t575 + (-t675 * t699 + t676 * t700 + Ifges(3,5)) * t667;
t552 = -mrSges(3,2) * g(1) - mrSges(3,3) * t639 + Ifges(3,5) * t669 - t667 * Ifges(3,6) - qJ(3) * t567 - t675 * t555 + t676 * t557;
t551 = -mrSges(2,2) * g(1) - mrSges(2,3) * t658 + Ifges(2,5) * qJDD(1) - t685 * Ifges(2,6) - pkin(6) * t562 + t683 * t552 - t679 * t553;
t550 = Ifges(2,6) * qJDD(1) + t685 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t657 + t679 * t552 + t683 * t553 - pkin(1) * (-m(3) * g(1) + t567) + pkin(6) * t704;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t567; -m(1) * g(2) + t694; -m(1) * g(3) + t705; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t689; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t705 - t684 * t550 - t680 * t551; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t694 - t680 * t550 + t684 * t551; t689; t691; t587; t686; -t690;];
tauJB = t1;
