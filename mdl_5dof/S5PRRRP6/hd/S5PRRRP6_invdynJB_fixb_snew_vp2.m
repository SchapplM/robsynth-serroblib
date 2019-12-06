% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRRP6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP6_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP6_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP6_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:51:03
% EndTime: 2019-12-05 16:51:07
% DurationCPUTime: 2.66s
% Computational Cost: add. (22959->237), mult. (45457->291), div. (0->0), fcn. (28051->8), ass. (0->100)
t673 = Ifges(5,1) + Ifges(6,1);
t663 = Ifges(5,4) - Ifges(6,5);
t671 = Ifges(6,4) + Ifges(5,5);
t672 = Ifges(5,2) + Ifges(6,3);
t669 = Ifges(5,6) - Ifges(6,6);
t670 = -Ifges(6,2) - Ifges(5,3);
t637 = sin(qJ(4));
t640 = cos(qJ(3));
t657 = qJD(2) * t640;
t638 = sin(qJ(3));
t658 = qJD(2) * t638;
t665 = cos(qJ(4));
t605 = t637 * t658 - t657 * t665;
t606 = (t637 * t640 + t638 * t665) * qJD(2);
t633 = qJD(3) + qJD(4);
t668 = t605 * t672 - t606 * t663 - t633 * t669;
t667 = -t605 * t663 + t606 * t673 + t633 * t671;
t636 = sin(pkin(8));
t662 = cos(pkin(8));
t619 = -g(1) * t662 - g(2) * t636;
t635 = -g(3) + qJDD(1);
t639 = sin(qJ(2));
t641 = cos(qJ(2));
t601 = t641 * t619 + t639 * t635;
t642 = qJD(2) ^ 2;
t595 = -pkin(2) * t642 + qJDD(2) * pkin(6) + t601;
t618 = g(1) * t636 - g(2) * t662;
t577 = -t638 * t595 - t640 * t618;
t656 = qJD(2) * qJD(3);
t653 = t640 * t656;
t616 = qJDD(2) * t638 + t653;
t563 = (-t616 + t653) * pkin(7) + (t638 * t640 * t642 + qJDD(3)) * pkin(3) + t577;
t578 = t640 * t595 - t638 * t618;
t617 = qJDD(2) * t640 - t638 * t656;
t622 = qJD(3) * pkin(3) - pkin(7) * t658;
t634 = t640 ^ 2;
t564 = -pkin(3) * t634 * t642 + pkin(7) * t617 - qJD(3) * t622 + t578;
t560 = t637 * t563 + t665 * t564;
t575 = qJD(4) * t606 + t616 * t637 - t617 * t665;
t597 = mrSges(5,1) * t633 - mrSges(5,3) * t606;
t632 = qJDD(3) + qJDD(4);
t588 = pkin(4) * t605 - qJ(5) * t606;
t631 = t633 ^ 2;
t554 = -pkin(4) * t631 + qJ(5) * t632 + 0.2e1 * qJD(5) * t633 - t588 * t605 + t560;
t598 = -mrSges(6,1) * t633 + mrSges(6,2) * t606;
t655 = m(6) * t554 + t632 * mrSges(6,3) + t633 * t598;
t589 = mrSges(6,1) * t605 - mrSges(6,3) * t606;
t659 = -mrSges(5,1) * t605 - mrSges(5,2) * t606 - t589;
t664 = -mrSges(5,3) - mrSges(6,2);
t544 = m(5) * t560 - t632 * mrSges(5,2) + t575 * t664 - t633 * t597 + t605 * t659 + t655;
t559 = t563 * t665 - t564 * t637;
t576 = -qJD(4) * t605 + t616 * t665 + t617 * t637;
t596 = -mrSges(5,2) * t633 - mrSges(5,3) * t605;
t555 = -pkin(4) * t632 - qJ(5) * t631 + t588 * t606 + qJDD(5) - t559;
t599 = -mrSges(6,2) * t605 + mrSges(6,3) * t633;
t649 = -m(6) * t555 + t632 * mrSges(6,1) + t633 * t599;
t546 = m(5) * t559 + t632 * mrSges(5,1) + t576 * t664 + t633 * t596 + t606 * t659 + t649;
t538 = t544 * t637 + t546 * t665;
t603 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t638 + Ifges(4,2) * t640) * qJD(2);
t604 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t638 + Ifges(4,4) * t640) * qJD(2);
t551 = t576 * mrSges(6,2) + t606 * t589 - t649;
t644 = -mrSges(5,1) * t559 + mrSges(6,1) * t555 + mrSges(5,2) * t560 - mrSges(6,3) * t554 + pkin(4) * t551 - qJ(5) * t655 + t670 * t632 + t668 * t606 + (qJ(5) * t589 - t667) * t605 - t671 * t576 + (mrSges(6,2) * qJ(5) + t669) * t575;
t666 = mrSges(4,1) * t577 - mrSges(4,2) * t578 + Ifges(4,5) * t616 + Ifges(4,6) * t617 + Ifges(4,3) * qJDD(3) + pkin(3) * t538 + (t603 * t638 - t604 * t640) * qJD(2) - t644;
t615 = (-mrSges(4,1) * t640 + mrSges(4,2) * t638) * qJD(2);
t621 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t657;
t536 = m(4) * t577 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t616 + qJD(3) * t621 - t615 * t658 + t538;
t620 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t658;
t650 = t544 * t665 - t546 * t637;
t537 = m(4) * t578 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t617 - qJD(3) * t620 + t615 * t657 + t650;
t532 = -t536 * t638 + t537 * t640;
t528 = m(3) * t601 - mrSges(3,1) * t642 - qJDD(2) * mrSges(3,2) + t532;
t600 = -t639 * t619 + t641 * t635;
t647 = -qJDD(2) * pkin(2) - t600;
t594 = -t642 * pkin(6) + t647;
t565 = -t617 * pkin(3) + t622 * t658 + (-pkin(7) * t634 - pkin(6)) * t642 + t647;
t557 = -0.2e1 * qJD(5) * t606 + (t605 * t633 - t576) * qJ(5) + (t606 * t633 + t575) * pkin(4) + t565;
t547 = m(6) * t557 + t575 * mrSges(6,1) - t576 * mrSges(6,3) - t606 * t598 + t605 * t599;
t646 = m(5) * t565 + t575 * mrSges(5,1) + t576 * mrSges(5,2) + t605 * t596 + t597 * t606 + t547;
t541 = -m(4) * t594 + t617 * mrSges(4,1) - mrSges(4,2) * t616 - t620 * t658 + t621 * t657 - t646;
t540 = m(3) * t600 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t642 + t541;
t651 = t528 * t641 - t540 * t639;
t524 = m(2) * t619 + t651;
t531 = t640 * t536 + t638 * t537;
t530 = (m(2) + m(3)) * t618 - t531;
t661 = t524 * t636 + t530 * t662;
t525 = t528 * t639 + t540 * t641;
t660 = t605 * t669 - t606 * t671 + t633 * t670;
t654 = m(2) * t635 + t525;
t652 = t524 * t662 - t530 * t636;
t533 = -mrSges(5,1) * t565 - mrSges(6,1) * t557 + mrSges(6,2) * t554 + mrSges(5,3) * t560 - pkin(4) * t547 - t575 * t672 + t663 * t576 + t660 * t606 + t669 * t632 + t667 * t633;
t534 = mrSges(5,2) * t565 + mrSges(6,2) * t555 - mrSges(5,3) * t559 - mrSges(6,3) * t557 - qJ(5) * t547 - t663 * t575 + t576 * t673 + t660 * t605 + t671 * t632 + t668 * t633;
t602 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t638 + Ifges(4,6) * t640) * qJD(2);
t520 = -mrSges(4,1) * t594 + mrSges(4,3) * t578 + Ifges(4,4) * t616 + Ifges(4,2) * t617 + Ifges(4,6) * qJDD(3) - pkin(3) * t646 + pkin(7) * t650 + qJD(3) * t604 + t533 * t665 + t637 * t534 - t602 * t658;
t521 = mrSges(4,2) * t594 - mrSges(4,3) * t577 + Ifges(4,1) * t616 + Ifges(4,4) * t617 + Ifges(4,5) * qJDD(3) - pkin(7) * t538 - qJD(3) * t603 - t533 * t637 + t534 * t665 + t602 * t657;
t645 = mrSges(3,1) * t600 - mrSges(3,2) * t601 + Ifges(3,3) * qJDD(2) + pkin(2) * t541 + pkin(6) * t532 + t520 * t640 + t521 * t638;
t519 = mrSges(3,1) * t618 + mrSges(3,3) * t601 + t642 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t531 - t666;
t518 = -mrSges(3,2) * t618 - mrSges(3,3) * t600 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t642 - pkin(6) * t531 - t520 * t638 + t521 * t640;
t517 = -mrSges(2,1) * t635 + mrSges(2,3) * t619 - pkin(1) * t525 - t645;
t516 = mrSges(2,2) * t635 - mrSges(2,3) * t618 - pkin(5) * t525 + t518 * t641 - t519 * t639;
t1 = [-m(1) * g(1) + t652; -m(1) * g(2) + t661; -m(1) * g(3) + t654; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t661 + t516 * t662 - t517 * t636; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t652 + t636 * t516 + t517 * t662; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t618 - mrSges(2,2) * t619 + t639 * t518 + t641 * t519 + pkin(1) * (m(3) * t618 - t531) + pkin(5) * t651; t654; t645; t666; -t644; t551;];
tauJB = t1;
