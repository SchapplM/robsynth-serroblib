% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:38
% EndTime: 2019-12-05 17:53:43
% DurationCPUTime: 5.30s
% Computational Cost: add. (62783->270), mult. (135182->344), div. (0->0), fcn. (85870->10), ass. (0->109)
t675 = sin(qJ(1));
t678 = cos(qJ(1));
t655 = t678 * g(2) + t675 * g(3);
t645 = qJDD(1) * pkin(1) + t655;
t654 = t675 * g(2) - t678 * g(3);
t679 = qJD(1) ^ 2;
t647 = -t679 * pkin(1) + t654;
t670 = sin(pkin(8));
t672 = cos(pkin(8));
t625 = t670 * t645 + t672 * t647;
t619 = -t679 * pkin(2) + qJDD(1) * pkin(6) + t625;
t668 = -g(1) + qJDD(2);
t674 = sin(qJ(3));
t677 = cos(qJ(3));
t608 = -t674 * t619 + t677 * t668;
t694 = qJD(1) * qJD(3);
t693 = t677 * t694;
t648 = t674 * qJDD(1) + t693;
t605 = (-t648 + t693) * qJ(4) + (t674 * t677 * t679 + qJDD(3)) * pkin(3) + t608;
t609 = t677 * t619 + t674 * t668;
t649 = t677 * qJDD(1) - t674 * t694;
t696 = qJD(1) * t674;
t651 = qJD(3) * pkin(3) - qJ(4) * t696;
t667 = t677 ^ 2;
t606 = -t667 * t679 * pkin(3) + t649 * qJ(4) - qJD(3) * t651 + t609;
t669 = sin(pkin(9));
t671 = cos(pkin(9));
t635 = (t669 * t677 + t671 * t674) * qJD(1);
t585 = -0.2e1 * qJD(4) * t635 + t671 * t605 - t669 * t606;
t627 = t671 * t648 + t669 * t649;
t634 = (-t669 * t674 + t671 * t677) * qJD(1);
t583 = (qJD(3) * t634 - t627) * pkin(7) + (t634 * t635 + qJDD(3)) * pkin(4) + t585;
t586 = 0.2e1 * qJD(4) * t634 + t669 * t605 + t671 * t606;
t626 = -t669 * t648 + t671 * t649;
t630 = qJD(3) * pkin(4) - t635 * pkin(7);
t633 = t634 ^ 2;
t584 = -t633 * pkin(4) + t626 * pkin(7) - qJD(3) * t630 + t586;
t673 = sin(qJ(5));
t676 = cos(qJ(5));
t581 = t676 * t583 - t673 * t584;
t616 = t676 * t634 - t673 * t635;
t595 = t616 * qJD(5) + t673 * t626 + t676 * t627;
t617 = t673 * t634 + t676 * t635;
t604 = -t616 * mrSges(6,1) + t617 * mrSges(6,2);
t664 = qJD(3) + qJD(5);
t610 = -t664 * mrSges(6,2) + t616 * mrSges(6,3);
t663 = qJDD(3) + qJDD(5);
t577 = m(6) * t581 + t663 * mrSges(6,1) - t595 * mrSges(6,3) - t617 * t604 + t664 * t610;
t582 = t673 * t583 + t676 * t584;
t594 = -t617 * qJD(5) + t676 * t626 - t673 * t627;
t611 = t664 * mrSges(6,1) - t617 * mrSges(6,3);
t578 = m(6) * t582 - t663 * mrSges(6,2) + t594 * mrSges(6,3) + t616 * t604 - t664 * t611;
t568 = t676 * t577 + t673 * t578;
t621 = -t634 * mrSges(5,1) + t635 * mrSges(5,2);
t628 = -qJD(3) * mrSges(5,2) + t634 * mrSges(5,3);
t566 = m(5) * t585 + qJDD(3) * mrSges(5,1) - t627 * mrSges(5,3) + qJD(3) * t628 - t635 * t621 + t568;
t629 = qJD(3) * mrSges(5,1) - t635 * mrSges(5,3);
t688 = -t673 * t577 + t676 * t578;
t567 = m(5) * t586 - qJDD(3) * mrSges(5,2) + t626 * mrSges(5,3) - qJD(3) * t629 + t634 * t621 + t688;
t562 = t671 * t566 + t669 * t567;
t614 = Ifges(5,4) * t635 + Ifges(5,2) * t634 + Ifges(5,6) * qJD(3);
t615 = Ifges(5,1) * t635 + Ifges(5,4) * t634 + Ifges(5,5) * qJD(3);
t640 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t674 + Ifges(4,2) * t677) * qJD(1);
t641 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t674 + Ifges(4,4) * t677) * qJD(1);
t597 = Ifges(6,4) * t617 + Ifges(6,2) * t616 + Ifges(6,6) * t664;
t598 = Ifges(6,1) * t617 + Ifges(6,4) * t616 + Ifges(6,5) * t664;
t683 = -mrSges(6,1) * t581 + mrSges(6,2) * t582 - Ifges(6,5) * t595 - Ifges(6,6) * t594 - Ifges(6,3) * t663 - t617 * t597 + t616 * t598;
t698 = mrSges(4,1) * t608 + mrSges(5,1) * t585 - mrSges(4,2) * t609 - mrSges(5,2) * t586 + Ifges(4,5) * t648 + Ifges(5,5) * t627 + Ifges(4,6) * t649 + Ifges(5,6) * t626 + pkin(3) * t562 + pkin(4) * t568 + (t674 * t640 - t677 * t641) * qJD(1) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t635 * t614 - t634 * t615 - t683;
t646 = (-mrSges(4,1) * t677 + mrSges(4,2) * t674) * qJD(1);
t695 = qJD(1) * t677;
t653 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t695;
t560 = m(4) * t608 + qJDD(3) * mrSges(4,1) - t648 * mrSges(4,3) + qJD(3) * t653 - t646 * t696 + t562;
t652 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t696;
t689 = -t669 * t566 + t671 * t567;
t561 = m(4) * t609 - qJDD(3) * mrSges(4,2) + t649 * mrSges(4,3) - qJD(3) * t652 + t646 * t695 + t689;
t690 = -t674 * t560 + t677 * t561;
t551 = m(3) * t625 - t679 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t690;
t624 = t672 * t645 - t670 * t647;
t684 = -qJDD(1) * pkin(2) - t624;
t607 = -t649 * pkin(3) + qJDD(4) + t651 * t696 + (-qJ(4) * t667 - pkin(6)) * t679 + t684;
t588 = -t626 * pkin(4) - t633 * pkin(7) + t635 * t630 + t607;
t687 = m(6) * t588 - t594 * mrSges(6,1) + t595 * mrSges(6,2) - t616 * t610 + t617 * t611;
t579 = m(5) * t607 - t626 * mrSges(5,1) + t627 * mrSges(5,2) - t634 * t628 + t635 * t629 + t687;
t618 = -t679 * pkin(6) + t684;
t681 = -m(4) * t618 + t649 * mrSges(4,1) - t648 * mrSges(4,2) - t652 * t696 + t653 * t695 - t579;
t572 = m(3) * t624 + qJDD(1) * mrSges(3,1) - t679 * mrSges(3,2) + t681;
t548 = t670 * t551 + t672 * t572;
t554 = t677 * t560 + t674 * t561;
t552 = m(3) * t668 + t554;
t691 = t672 * t551 - t670 * t572;
t545 = m(2) * t654 - t679 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t691;
t546 = m(2) * t655 + qJDD(1) * mrSges(2,1) - t679 * mrSges(2,2) + t548;
t692 = t678 * t545 - t675 * t546;
t686 = -t675 * t545 - t678 * t546;
t596 = Ifges(6,5) * t617 + Ifges(6,6) * t616 + Ifges(6,3) * t664;
t569 = -mrSges(6,1) * t588 + mrSges(6,3) * t582 + Ifges(6,4) * t595 + Ifges(6,2) * t594 + Ifges(6,6) * t663 - t617 * t596 + t664 * t598;
t570 = mrSges(6,2) * t588 - mrSges(6,3) * t581 + Ifges(6,1) * t595 + Ifges(6,4) * t594 + Ifges(6,5) * t663 + t616 * t596 - t664 * t597;
t613 = Ifges(5,5) * t635 + Ifges(5,6) * t634 + Ifges(5,3) * qJD(3);
t555 = -mrSges(5,1) * t607 + mrSges(5,3) * t586 + Ifges(5,4) * t627 + Ifges(5,2) * t626 + Ifges(5,6) * qJDD(3) - pkin(4) * t687 + pkin(7) * t688 + qJD(3) * t615 + t676 * t569 + t673 * t570 - t635 * t613;
t556 = mrSges(5,2) * t607 - mrSges(5,3) * t585 + Ifges(5,1) * t627 + Ifges(5,4) * t626 + Ifges(5,5) * qJDD(3) - pkin(7) * t568 - qJD(3) * t614 - t673 * t569 + t676 * t570 + t634 * t613;
t639 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t674 + Ifges(4,6) * t677) * qJD(1);
t541 = -mrSges(4,1) * t618 + mrSges(4,3) * t609 + Ifges(4,4) * t648 + Ifges(4,2) * t649 + Ifges(4,6) * qJDD(3) - pkin(3) * t579 + qJ(4) * t689 + qJD(3) * t641 + t671 * t555 + t669 * t556 - t639 * t696;
t543 = mrSges(4,2) * t618 - mrSges(4,3) * t608 + Ifges(4,1) * t648 + Ifges(4,4) * t649 + Ifges(4,5) * qJDD(3) - qJ(4) * t562 - qJD(3) * t640 - t669 * t555 + t671 * t556 + t639 * t695;
t682 = mrSges(2,1) * t655 + mrSges(3,1) * t624 - mrSges(2,2) * t654 - mrSges(3,2) * t625 + pkin(1) * t548 + pkin(2) * t681 + pkin(6) * t690 + t677 * t541 + t674 * t543 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t539 = -mrSges(3,1) * t668 + mrSges(3,3) * t625 + t679 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t554 - t698;
t538 = mrSges(3,2) * t668 - mrSges(3,3) * t624 + Ifges(3,5) * qJDD(1) - t679 * Ifges(3,6) - pkin(6) * t554 - t674 * t541 + t677 * t543;
t537 = -mrSges(2,2) * g(1) - mrSges(2,3) * t655 + Ifges(2,5) * qJDD(1) - t679 * Ifges(2,6) - qJ(2) * t548 + t672 * t538 - t670 * t539;
t536 = mrSges(2,1) * g(1) + mrSges(2,3) * t654 + t679 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t552 + qJ(2) * t691 + t670 * t538 + t672 * t539;
t1 = [(-m(1) - m(2)) * g(1) + t552; -m(1) * g(2) + t686; -m(1) * g(3) + t692; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t682; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t692 - t678 * t536 - t675 * t537; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t686 - t675 * t536 + t678 * t537; t682; t552; t698; t579; -t683;];
tauJB = t1;
