% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:32
% EndTime: 2019-12-05 17:59:35
% DurationCPUTime: 2.11s
% Computational Cost: add. (16202->245), mult. (32220->284), div. (0->0), fcn. (18535->6), ass. (0->99)
t706 = Ifges(5,4) + Ifges(6,4);
t717 = Ifges(5,2) + Ifges(6,2);
t714 = Ifges(5,6) + Ifges(6,6);
t716 = Ifges(5,1) + Ifges(6,1);
t715 = Ifges(5,5) + Ifges(6,5);
t713 = Ifges(5,3) + Ifges(6,3);
t671 = sin(qJ(4));
t672 = sin(qJ(3));
t674 = cos(qJ(4));
t675 = cos(qJ(3));
t638 = (-t671 * t675 - t672 * t674) * qJD(1);
t639 = (-t671 * t672 + t674 * t675) * qJD(1);
t662 = qJD(3) + qJD(4);
t712 = t717 * t638 + t706 * t639 + t714 * t662;
t677 = qJD(1) ^ 2;
t673 = sin(qJ(1));
t676 = cos(qJ(1));
t654 = -t676 * g(1) - t673 * g(2);
t685 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t654;
t709 = -pkin(1) - pkin(6);
t627 = t709 * t677 + t685;
t696 = qJD(1) * qJD(3);
t647 = -t672 * qJDD(1) - t675 * t696;
t691 = t672 * t696;
t648 = t675 * qJDD(1) - t691;
t698 = qJD(1) * t672;
t650 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t698;
t697 = qJD(1) * t675;
t651 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t697;
t652 = qJD(3) * pkin(3) - pkin(7) * t697;
t668 = t672 ^ 2;
t594 = -t647 * pkin(3) + t652 * t697 + (-pkin(7) * t668 + t709) * t677 + t685;
t603 = -t639 * qJD(4) + t674 * t647 - t671 * t648;
t604 = t638 * qJD(4) + t671 * t647 + t674 * t648;
t622 = -t662 * mrSges(6,2) + t638 * mrSges(6,3);
t623 = -t662 * mrSges(5,2) + t638 * mrSges(5,3);
t626 = t662 * mrSges(5,1) - t639 * mrSges(5,3);
t624 = t662 * pkin(4) - t639 * qJ(5);
t634 = t638 ^ 2;
t582 = -t603 * pkin(4) - t634 * qJ(5) + t639 * t624 + qJDD(5) + t594;
t625 = t662 * mrSges(6,1) - t639 * mrSges(6,3);
t692 = m(6) * t582 + t604 * mrSges(6,2) + t639 * t625;
t682 = m(5) * t594 + t604 * mrSges(5,2) - (mrSges(5,1) + mrSges(6,1)) * t603 + t639 * t626 - (t622 + t623) * t638 + t692;
t711 = -m(4) * t627 + t647 * mrSges(4,1) - t648 * mrSges(4,2) - t650 * t698 - t651 * t697 - t682;
t708 = mrSges(2,1) - mrSges(3,2);
t705 = Ifges(2,5) - Ifges(3,4);
t704 = -Ifges(2,6) + Ifges(3,5);
t653 = t673 * g(1) - t676 * g(2);
t684 = -t677 * qJ(2) + qJDD(2) - t653;
t628 = t709 * qJDD(1) + t684;
t616 = t672 * g(3) + t675 * t628;
t589 = (-t648 - t691) * pkin(7) + (-t672 * t675 * t677 + qJDD(3)) * pkin(3) + t616;
t617 = -t675 * g(3) + t672 * t628;
t590 = -t668 * t677 * pkin(3) + t647 * pkin(7) - qJD(3) * t652 + t617;
t584 = t674 * t589 - t671 * t590;
t613 = -t638 * mrSges(6,1) + t639 * mrSges(6,2);
t614 = -t638 * mrSges(5,1) + t639 * mrSges(5,2);
t661 = qJDD(3) + qJDD(4);
t578 = -0.2e1 * qJD(5) * t639 + (t638 * t662 - t604) * qJ(5) + (t638 * t639 + t661) * pkin(4) + t584;
t694 = m(6) * t578 + t661 * mrSges(6,1) + t662 * t622;
t569 = m(5) * t584 + t661 * mrSges(5,1) + t662 * t623 + (-t613 - t614) * t639 + (-mrSges(5,3) - mrSges(6,3)) * t604 + t694;
t585 = t671 * t589 + t674 * t590;
t580 = -t634 * pkin(4) + t603 * qJ(5) + 0.2e1 * qJD(5) * t638 - t662 * t624 + t585;
t693 = m(6) * t580 + t603 * mrSges(6,3) + t638 * t613;
t572 = m(5) * t585 + t603 * mrSges(5,3) + t638 * t614 + (-t625 - t626) * t662 + (-mrSges(5,2) - mrSges(6,2)) * t661 + t693;
t563 = t674 * t569 + t671 * t572;
t646 = (mrSges(4,1) * t672 + mrSges(4,2) * t675) * qJD(1);
t560 = m(4) * t616 + qJDD(3) * mrSges(4,1) - t648 * mrSges(4,3) + qJD(3) * t650 - t646 * t697 + t563;
t688 = -t671 * t569 + t674 * t572;
t561 = m(4) * t617 - qJDD(3) * mrSges(4,2) + t647 * mrSges(4,3) - qJD(3) * t651 - t646 * t698 + t688;
t556 = t675 * t560 + t672 * t561;
t633 = -qJDD(1) * pkin(1) + t684;
t683 = -m(3) * t633 + t677 * mrSges(3,3) - t556;
t552 = m(2) * t653 - t677 * mrSges(2,2) + t708 * qJDD(1) + t683;
t631 = t677 * pkin(1) - t685;
t679 = -m(3) * t631 + t677 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t711;
t566 = m(2) * t654 - t677 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t679;
t702 = t676 * t552 + t673 * t566;
t701 = -t714 * t638 - t715 * t639 - t713 * t662;
t700 = -t706 * t638 - t716 * t639 - t715 * t662;
t690 = -t673 * t552 + t676 * t566;
t689 = -t672 * t560 + t675 * t561;
t575 = -t603 * mrSges(6,1) - t638 * t622 + t692;
t557 = -mrSges(5,1) * t594 + mrSges(5,3) * t585 - mrSges(6,1) * t582 + mrSges(6,3) * t580 - pkin(4) * t575 + qJ(5) * t693 + (-qJ(5) * t625 - t700) * t662 + (-qJ(5) * mrSges(6,2) + t714) * t661 + t701 * t639 + t706 * t604 + t717 * t603;
t574 = -t604 * mrSges(6,3) - t639 * t613 + t694;
t558 = mrSges(5,2) * t594 + mrSges(6,2) * t582 - mrSges(5,3) * t584 - mrSges(6,3) * t578 - qJ(5) * t574 + t706 * t603 + t716 * t604 - t701 * t638 + t715 * t661 - t712 * t662;
t635 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t675 - Ifges(4,6) * t672) * qJD(1);
t637 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t675 - Ifges(4,4) * t672) * qJD(1);
t548 = -mrSges(4,1) * t627 + mrSges(4,3) * t617 + Ifges(4,4) * t648 + Ifges(4,2) * t647 + Ifges(4,6) * qJDD(3) - pkin(3) * t682 + pkin(7) * t688 + qJD(3) * t637 + t674 * t557 + t671 * t558 - t635 * t697;
t636 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t675 - Ifges(4,2) * t672) * qJD(1);
t550 = mrSges(4,2) * t627 - mrSges(4,3) * t616 + Ifges(4,1) * t648 + Ifges(4,4) * t647 + Ifges(4,5) * qJDD(3) - pkin(7) * t563 - qJD(3) * t636 - t671 * t557 + t674 * t558 - t635 * t698;
t554 = qJDD(1) * mrSges(3,2) - t683;
t681 = mrSges(2,1) * t653 - mrSges(2,2) * t654 + mrSges(3,2) * t633 - mrSges(3,3) * t631 - pkin(1) * t554 - pkin(6) * t556 + qJ(2) * t679 - t672 * t548 + t675 * t550 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t680 = mrSges(5,1) * t584 + mrSges(6,1) * t578 - mrSges(5,2) * t585 - mrSges(6,2) * t580 + pkin(4) * t574 + t714 * t603 + t715 * t604 + t700 * t638 + t712 * t639 + t713 * t661;
t678 = mrSges(4,1) * t616 - mrSges(4,2) * t617 + Ifges(4,5) * t648 + Ifges(4,6) * t647 + Ifges(4,3) * qJDD(3) + pkin(3) * t563 + t636 * t697 + t637 * t698 + t680;
t555 = -m(3) * g(3) + t689;
t547 = t678 + t705 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t704 * t677 - mrSges(2,3) * t653 + mrSges(3,1) * t633 + pkin(2) * t556 - qJ(2) * t555;
t546 = -mrSges(3,1) * t631 + mrSges(2,3) * t654 - pkin(1) * t555 - pkin(2) * t711 - pkin(6) * t689 + t708 * g(3) - t704 * qJDD(1) - t675 * t548 - t672 * t550 + t705 * t677;
t1 = [-m(1) * g(1) + t690; -m(1) * g(2) + t702; (-m(1) - m(2) - m(3)) * g(3) + t689; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t702 - t673 * t546 + t676 * t547; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t690 + t676 * t546 + t673 * t547; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t681; t681; t554; t678; t680; t575;];
tauJB = t1;
