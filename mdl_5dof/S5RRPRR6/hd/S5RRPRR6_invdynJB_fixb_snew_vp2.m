% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR6
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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:36:01
% EndTime: 2019-12-05 18:36:05
% DurationCPUTime: 4.63s
% Computational Cost: add. (67358->251), mult. (92265->330), div. (0->0), fcn. (56620->10), ass. (0->117)
t674 = sin(qJ(1));
t678 = cos(qJ(1));
t654 = t678 * g(2) + t674 * g(3);
t644 = qJDD(1) * pkin(1) + t654;
t653 = t674 * g(2) - t678 * g(3);
t679 = qJD(1) ^ 2;
t645 = -t679 * pkin(1) + t653;
t673 = sin(qJ(2));
t677 = cos(qJ(2));
t626 = t673 * t644 + t677 * t645;
t666 = (qJD(1) + qJD(2));
t663 = t666 ^ 2;
t664 = qJDD(1) + qJDD(2);
t713 = -t663 * pkin(2) + t664 * qJ(3) + (2 * qJD(3) * t666) + t626;
t669 = sin(pkin(9));
t670 = cos(pkin(9));
t612 = -t670 * g(1) - t713 * t669;
t712 = mrSges(4,2) * t669;
t711 = mrSges(4,3) * t664;
t710 = t663 * t669 ^ 2;
t709 = t669 * t666;
t672 = sin(qJ(4));
t708 = t669 * t672;
t676 = cos(qJ(4));
t707 = t669 * t676;
t706 = t670 * t664;
t705 = t670 * t666;
t613 = -t669 * g(1) + t713 * t670;
t637 = (-mrSges(4,1) * t670 + t712) * t666;
t692 = -pkin(3) * t670 - pkin(7) * t669;
t639 = t692 * t666;
t599 = t639 * t705 + t613;
t625 = t677 * t644 - t673 * t645;
t687 = -t663 * qJ(3) + qJDD(3) - t625;
t611 = (-pkin(2) + t692) * t664 + t687;
t610 = t676 * t611;
t702 = qJD(4) * t666;
t636 = (t664 * t676 - t672 * t702) * t669;
t649 = qJDD(4) - t706;
t650 = qJD(4) - t705;
t591 = t649 * pkin(4) - t636 * pkin(8) + t610 + (-pkin(4) * t676 * t710 - pkin(8) * t650 * t709 - t599) * t672;
t594 = t676 * t599 + t672 * t611;
t699 = t666 * t707;
t634 = t650 * pkin(4) - pkin(8) * t699;
t635 = (-t664 * t672 - t676 * t702) * t669;
t701 = t672 ^ 2 * t710;
t592 = -pkin(4) * t701 + t635 * pkin(8) - t650 * t634 + t594;
t671 = sin(qJ(5));
t675 = cos(qJ(5));
t589 = t675 * t591 - t671 * t592;
t627 = (-t671 * t676 - t672 * t675) * t709;
t604 = t627 * qJD(5) + t671 * t635 + t675 * t636;
t628 = (-t671 * t672 + t675 * t676) * t709;
t614 = -t627 * mrSges(6,1) + t628 * mrSges(6,2);
t648 = qJD(5) + t650;
t619 = -t648 * mrSges(6,2) + t627 * mrSges(6,3);
t646 = qJDD(5) + t649;
t586 = m(6) * t589 + t646 * mrSges(6,1) - t604 * mrSges(6,3) - t628 * t614 + t648 * t619;
t590 = t671 * t591 + t675 * t592;
t603 = -t628 * qJD(5) + t675 * t635 - t671 * t636;
t620 = t648 * mrSges(6,1) - t628 * mrSges(6,3);
t587 = m(6) * t590 - t646 * mrSges(6,2) + t603 * mrSges(6,3) + t627 * t614 - t648 * t620;
t578 = t675 * t586 + t671 * t587;
t593 = -t672 * t599 + t610;
t700 = t666 * t708;
t631 = -t650 * mrSges(5,2) - mrSges(5,3) * t700;
t633 = (mrSges(5,1) * t672 + mrSges(5,2) * t676) * t709;
t576 = m(5) * t593 + t649 * mrSges(5,1) - t636 * mrSges(5,3) + t650 * t631 - t633 * t699 + t578;
t632 = t650 * mrSges(5,1) - mrSges(5,3) * t699;
t693 = -t671 * t586 + t675 * t587;
t577 = m(5) * t594 - t649 * mrSges(5,2) + t635 * mrSges(5,3) - t650 * t632 - t633 * t700 + t693;
t694 = -t672 * t576 + t676 * t577;
t571 = m(4) * t613 + (t637 * t666 + t711) * t670 + t694;
t598 = t639 * t709 - t612;
t595 = -t635 * pkin(4) - pkin(8) * t701 + t634 * t699 + t598;
t685 = m(6) * t595 - t603 * mrSges(6,1) + t604 * mrSges(6,2) - t627 * t619 + t628 * t620;
t681 = -m(5) * t598 + t635 * mrSges(5,1) - t636 * mrSges(5,2) - t685;
t582 = m(4) * t612 + (-t711 + (-t631 * t672 - t632 * t676 - t637) * t666) * t669 + t681;
t695 = t670 * t571 - t669 * t582;
t562 = m(3) * t626 - t663 * mrSges(3,1) - t664 * mrSges(3,2) + t695;
t574 = t676 * t576 + t672 * t577;
t617 = -t664 * pkin(2) + t687;
t683 = -m(4) * t617 + mrSges(4,1) * t706 - t574 + (t663 * t670 ^ 2 + t710) * mrSges(4,3);
t568 = m(3) * t625 - t663 * mrSges(3,2) + (mrSges(3,1) - t712) * t664 + t683;
t557 = t673 * t562 + t677 * t568;
t565 = t669 * t571 + t670 * t582;
t696 = t677 * t562 - t673 * t568;
t554 = m(2) * t653 - t679 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t696;
t555 = m(2) * t654 + qJDD(1) * mrSges(2,1) - t679 * mrSges(2,2) + t557;
t697 = t678 * t554 - t674 * t555;
t691 = Ifges(4,1) * t669 + Ifges(4,4) * t670;
t690 = Ifges(4,5) * t669 + Ifges(4,6) * t670;
t689 = -t674 * t554 - t678 * t555;
t622 = Ifges(5,6) * t650 + (Ifges(5,4) * t676 - Ifges(5,2) * t672) * t709;
t623 = Ifges(5,5) * t650 + (Ifges(5,1) * t676 - Ifges(5,4) * t672) * t709;
t688 = t622 * t676 + t623 * t672;
t605 = Ifges(6,5) * t628 + Ifges(6,6) * t627 + Ifges(6,3) * t648;
t607 = Ifges(6,1) * t628 + Ifges(6,4) * t627 + Ifges(6,5) * t648;
t579 = -mrSges(6,1) * t595 + mrSges(6,3) * t590 + Ifges(6,4) * t604 + Ifges(6,2) * t603 + Ifges(6,6) * t646 - t628 * t605 + t648 * t607;
t606 = Ifges(6,4) * t628 + Ifges(6,2) * t627 + Ifges(6,6) * t648;
t580 = mrSges(6,2) * t595 - mrSges(6,3) * t589 + Ifges(6,1) * t604 + Ifges(6,4) * t603 + Ifges(6,5) * t646 + t627 * t605 - t648 * t606;
t621 = Ifges(5,3) * t650 + (Ifges(5,5) * t676 - Ifges(5,6) * t672) * t709;
t563 = -mrSges(5,1) * t598 + mrSges(5,3) * t594 + Ifges(5,4) * t636 + Ifges(5,2) * t635 + Ifges(5,6) * t649 - pkin(4) * t685 + pkin(8) * t693 + t675 * t579 + t671 * t580 - t621 * t699 + t650 * t623;
t566 = mrSges(5,2) * t598 - mrSges(5,3) * t593 + Ifges(5,1) * t636 + Ifges(5,4) * t635 + Ifges(5,5) * t649 - pkin(8) * t578 - t671 * t579 + t675 * t580 - t621 * t700 - t650 * t622;
t638 = t690 * t666;
t552 = mrSges(4,2) * t617 - mrSges(4,3) * t612 - pkin(7) * t574 - t672 * t563 + t676 * t566 + t638 * t705 + t691 * t664;
t684 = -mrSges(6,1) * t589 + mrSges(6,2) * t590 - Ifges(6,5) * t604 - Ifges(6,6) * t603 - Ifges(6,3) * t646 - t628 * t606 + t627 * t607;
t680 = mrSges(5,1) * t593 - mrSges(5,2) * t594 + Ifges(5,5) * t636 + Ifges(5,6) * t635 + Ifges(5,3) * t649 + pkin(4) * t578 - t684;
t559 = -t680 - pkin(3) * t574 + mrSges(4,3) * t613 + Ifges(4,2) * t706 - mrSges(4,1) * t617 + (Ifges(4,4) * t664 + (-t638 - t688) * t666) * t669;
t573 = t664 * t712 - t683;
t686 = mrSges(3,1) * t625 - mrSges(3,2) * t626 + Ifges(3,3) * t664 - pkin(2) * t573 + qJ(3) * t695 + t669 * t552 + t670 * t559;
t682 = mrSges(2,1) * t654 - mrSges(2,2) * t653 + Ifges(2,3) * qJDD(1) + pkin(1) * t557 + t686;
t550 = t663 * Ifges(3,5) + mrSges(3,1) * g(1) + mrSges(3,3) * t626 - mrSges(4,1) * t612 + mrSges(4,2) * t613 - t672 * t566 - t676 * t563 - pkin(3) * t681 - pkin(7) * t694 - pkin(2) * t565 + (Ifges(3,6) - t690) * t664 + (-pkin(3) * (-t631 * t708 - t632 * t707) + (-t669 * (Ifges(4,4) * t669 + Ifges(4,2) * t670) + t670 * t691) * t666) * t666;
t549 = -mrSges(3,2) * g(1) - mrSges(3,3) * t625 + Ifges(3,5) * t664 - t663 * Ifges(3,6) - qJ(3) * t565 + t670 * t552 - t669 * t559;
t548 = -mrSges(2,2) * g(1) - mrSges(2,3) * t654 + Ifges(2,5) * qJDD(1) - t679 * Ifges(2,6) - pkin(6) * t557 + t677 * t549 - t673 * t550;
t547 = Ifges(2,6) * qJDD(1) + t679 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t653 + t673 * t549 + t677 * t550 - pkin(1) * (-m(3) * g(1) + t565) + pkin(6) * t696;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t565; -m(1) * g(2) + t689; -m(1) * g(3) + t697; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t682; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t697 - t678 * t547 - t674 * t548; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t689 - t674 * t547 + t678 * t548; t682; t686; t573; t688 * t709 + t680; -t684;];
tauJB = t1;
