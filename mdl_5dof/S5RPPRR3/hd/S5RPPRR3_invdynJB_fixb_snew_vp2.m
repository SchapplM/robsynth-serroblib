% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:28:07
% EndTime: 2020-01-03 11:28:12
% DurationCPUTime: 5.48s
% Computational Cost: add. (55727->249), mult. (123458->311), div. (0->0), fcn. (83128->10), ass. (0->112)
t678 = qJD(1) ^ 2;
t670 = cos(pkin(9));
t705 = pkin(3) * t670;
t668 = sin(pkin(9));
t704 = mrSges(4,2) * t668;
t663 = t670 ^ 2;
t703 = t663 * t678;
t674 = sin(qJ(1));
t677 = cos(qJ(1));
t649 = -t674 * g(2) + t677 * g(3);
t650 = -t677 * g(2) - t674 * g(3);
t646 = qJDD(1) * pkin(1) + t650;
t647 = -t678 * pkin(1) + t649;
t669 = sin(pkin(8));
t671 = cos(pkin(8));
t634 = t669 * t646 + t671 * t647;
t626 = -t678 * pkin(2) + qJDD(1) * qJ(3) + t634;
t667 = -g(1) + qJDD(2);
t698 = qJD(1) * qJD(3);
t701 = t670 * t667 - 0.2e1 * t668 * t698;
t612 = (-pkin(6) * qJDD(1) + t678 * t705 - t626) * t668 + t701;
t616 = t668 * t667 + (t626 + 0.2e1 * t698) * t670;
t697 = qJDD(1) * t670;
t613 = -pkin(3) * t703 + pkin(6) * t697 + t616;
t673 = sin(qJ(4));
t676 = cos(qJ(4));
t594 = t676 * t612 - t673 * t613;
t686 = t668 * t676 + t670 * t673;
t685 = -t668 * t673 + t670 * t676;
t639 = t685 * qJD(1);
t699 = t639 * qJD(4);
t631 = t686 * qJDD(1) + t699;
t640 = t686 * qJD(1);
t590 = (-t631 + t699) * pkin(7) + (t639 * t640 + qJDD(4)) * pkin(4) + t594;
t595 = t673 * t612 + t676 * t613;
t630 = -t640 * qJD(4) + t685 * qJDD(1);
t637 = qJD(4) * pkin(4) - t640 * pkin(7);
t638 = t639 ^ 2;
t591 = -t638 * pkin(4) + t630 * pkin(7) - qJD(4) * t637 + t595;
t672 = sin(qJ(5));
t675 = cos(qJ(5));
t588 = t675 * t590 - t672 * t591;
t624 = t675 * t639 - t672 * t640;
t602 = t624 * qJD(5) + t672 * t630 + t675 * t631;
t625 = t672 * t639 + t675 * t640;
t608 = -t624 * mrSges(6,1) + t625 * mrSges(6,2);
t664 = qJD(4) + qJD(5);
t617 = -t664 * mrSges(6,2) + t624 * mrSges(6,3);
t661 = qJDD(4) + qJDD(5);
t585 = m(6) * t588 + t661 * mrSges(6,1) - t602 * mrSges(6,3) - t625 * t608 + t664 * t617;
t589 = t672 * t590 + t675 * t591;
t601 = -t625 * qJD(5) + t675 * t630 - t672 * t631;
t618 = t664 * mrSges(6,1) - t625 * mrSges(6,3);
t586 = m(6) * t589 - t661 * mrSges(6,2) + t601 * mrSges(6,3) + t624 * t608 - t664 * t618;
t575 = t675 * t585 + t672 * t586;
t628 = -t639 * mrSges(5,1) + t640 * mrSges(5,2);
t635 = -qJD(4) * mrSges(5,2) + t639 * mrSges(5,3);
t573 = m(5) * t594 + qJDD(4) * mrSges(5,1) - t631 * mrSges(5,3) + qJD(4) * t635 - t640 * t628 + t575;
t636 = qJD(4) * mrSges(5,1) - t640 * mrSges(5,3);
t692 = -t672 * t585 + t675 * t586;
t574 = m(5) * t595 - qJDD(4) * mrSges(5,2) + t630 * mrSges(5,3) - qJD(4) * t636 + t639 * t628 + t692;
t569 = t676 * t573 + t673 * t574;
t615 = -t668 * t626 + t701;
t684 = mrSges(4,3) * qJDD(1) + t678 * (-mrSges(4,1) * t670 + t704);
t567 = m(4) * t615 - t684 * t668 + t569;
t693 = -t673 * t573 + t676 * t574;
t568 = m(4) * t616 + t684 * t670 + t693;
t694 = -t668 * t567 + t670 * t568;
t558 = m(3) * t634 - t678 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t694;
t633 = t671 * t646 - t669 * t647;
t688 = qJDD(3) - t633;
t620 = -qJDD(1) * pkin(2) - t678 * qJ(3) + t688;
t662 = t668 ^ 2;
t614 = (-pkin(2) - t705) * qJDD(1) + (-qJ(3) + (-t662 - t663) * pkin(6)) * t678 + t688;
t593 = -t630 * pkin(4) - t638 * pkin(7) + t640 * t637 + t614;
t687 = m(6) * t593 - t601 * mrSges(6,1) + t602 * mrSges(6,2) - t624 * t617 + t625 * t618;
t681 = m(5) * t614 - t630 * mrSges(5,1) + t631 * mrSges(5,2) - t639 * t635 + t640 * t636 + t687;
t680 = -m(4) * t620 + mrSges(4,1) * t697 - t681 + (t662 * t678 + t703) * mrSges(4,3);
t579 = t680 + (mrSges(3,1) - t704) * qJDD(1) - t678 * mrSges(3,2) + m(3) * t633;
t695 = t671 * t558 - t669 * t579;
t552 = m(2) * t649 - t678 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t695;
t555 = t669 * t558 + t671 * t579;
t553 = m(2) * t650 + qJDD(1) * mrSges(2,1) - t678 * mrSges(2,2) + t555;
t702 = t674 * t552 + t677 * t553;
t561 = t670 * t567 + t668 * t568;
t689 = Ifges(4,5) * t668 + Ifges(4,6) * t670;
t700 = t678 * t689;
t559 = m(3) * t667 + t561;
t696 = -t677 * t552 + t674 * t553;
t691 = Ifges(4,1) * t668 + Ifges(4,4) * t670;
t690 = Ifges(4,4) * t668 + Ifges(4,2) * t670;
t604 = Ifges(6,4) * t625 + Ifges(6,2) * t624 + Ifges(6,6) * t664;
t605 = Ifges(6,1) * t625 + Ifges(6,4) * t624 + Ifges(6,5) * t664;
t683 = -mrSges(6,1) * t588 + mrSges(6,2) * t589 - Ifges(6,5) * t602 - Ifges(6,6) * t601 - Ifges(6,3) * t661 - t625 * t604 + t624 * t605;
t603 = Ifges(6,5) * t625 + Ifges(6,6) * t624 + Ifges(6,3) * t664;
t576 = -mrSges(6,1) * t593 + mrSges(6,3) * t589 + Ifges(6,4) * t602 + Ifges(6,2) * t601 + Ifges(6,6) * t661 - t625 * t603 + t664 * t605;
t577 = mrSges(6,2) * t593 - mrSges(6,3) * t588 + Ifges(6,1) * t602 + Ifges(6,4) * t601 + Ifges(6,5) * t661 + t624 * t603 - t664 * t604;
t621 = Ifges(5,5) * t640 + Ifges(5,6) * t639 + Ifges(5,3) * qJD(4);
t623 = Ifges(5,1) * t640 + Ifges(5,4) * t639 + Ifges(5,5) * qJD(4);
t562 = -mrSges(5,1) * t614 + mrSges(5,3) * t595 + Ifges(5,4) * t631 + Ifges(5,2) * t630 + Ifges(5,6) * qJDD(4) - pkin(4) * t687 + pkin(7) * t692 + qJD(4) * t623 + t675 * t576 + t672 * t577 - t640 * t621;
t622 = Ifges(5,4) * t640 + Ifges(5,2) * t639 + Ifges(5,6) * qJD(4);
t563 = mrSges(5,2) * t614 - mrSges(5,3) * t594 + Ifges(5,1) * t631 + Ifges(5,4) * t630 + Ifges(5,5) * qJDD(4) - pkin(7) * t575 - qJD(4) * t622 - t672 * t576 + t675 * t577 + t639 * t621;
t546 = -mrSges(4,1) * t620 + mrSges(4,3) * t616 - pkin(3) * t681 + pkin(6) * t693 + t690 * qJDD(1) + t676 * t562 + t673 * t563 - t668 * t700;
t548 = mrSges(4,2) * t620 - mrSges(4,3) * t615 - pkin(6) * t569 + t691 * qJDD(1) - t673 * t562 + t676 * t563 + t670 * t700;
t581 = qJDD(1) * t704 - t680;
t682 = mrSges(2,1) * t650 + mrSges(3,1) * t633 - mrSges(2,2) * t649 - mrSges(3,2) * t634 + pkin(1) * t555 - pkin(2) * t581 + qJ(3) * t694 + t670 * t546 + t668 * t548 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t679 = mrSges(5,1) * t594 - mrSges(5,2) * t595 + Ifges(5,5) * t631 + Ifges(5,6) * t630 + Ifges(5,3) * qJDD(4) + pkin(4) * t575 + t640 * t622 - t639 * t623 - t683;
t544 = -t679 + (Ifges(3,6) - t689) * qJDD(1) - mrSges(3,1) * t667 + mrSges(3,3) * t634 - mrSges(4,1) * t615 + mrSges(4,2) * t616 - pkin(3) * t569 - pkin(2) * t561 + (-t668 * t690 + t670 * t691 + Ifges(3,5)) * t678;
t543 = mrSges(3,2) * t667 - mrSges(3,3) * t633 + Ifges(3,5) * qJDD(1) - t678 * Ifges(3,6) - qJ(3) * t561 - t668 * t546 + t670 * t548;
t542 = -mrSges(2,2) * g(1) - mrSges(2,3) * t650 + Ifges(2,5) * qJDD(1) - t678 * Ifges(2,6) - qJ(2) * t555 + t671 * t543 - t669 * t544;
t541 = mrSges(2,1) * g(1) + mrSges(2,3) * t649 + t678 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t559 + qJ(2) * t695 + t669 * t543 + t671 * t544;
t1 = [(-m(1) - m(2)) * g(1) + t559; -m(1) * g(2) + t702; -m(1) * g(3) + t696; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t682; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t696 + t677 * t541 + t674 * t542; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t702 + t674 * t541 - t677 * t542; t682; t559; t581; t679; -t683;];
tauJB = t1;
