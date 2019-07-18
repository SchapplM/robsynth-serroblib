% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR3
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_invdynJB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:17:45
% EndTime: 2019-07-18 17:17:52
% DurationCPUTime: 5.60s
% Computational Cost: add. (67429->292), mult. (136740->372), div. (0->0), fcn. (102620->10), ass. (0->115)
t668 = sin(qJ(1));
t673 = cos(qJ(1));
t655 = -g(1) * t673 - g(2) * t668;
t667 = sin(qJ(2));
t672 = cos(qJ(2));
t638 = -t672 * g(3) - t667 * t655;
t674 = qJD(1) ^ 2;
t630 = (t667 * t672 * t674 + qJDD(2)) * pkin(1) + t638;
t639 = -t667 * g(3) + t672 * t655;
t633 = (-t672 ^ 2 * t674 - qJD(2) ^ 2) * pkin(1) + t639;
t666 = sin(qJ(3));
t671 = cos(qJ(3));
t606 = t666 * t630 + t671 * t633;
t645 = (t666 * t672 + t667 * t671) * qJD(1);
t687 = qJD(1) * qJD(2);
t649 = qJDD(1) * t667 + t672 * t687;
t686 = t667 * t687;
t650 = qJDD(1) * t672 - t686;
t617 = -t645 * qJD(3) - t666 * t649 + t650 * t671;
t688 = qJD(1) * t672;
t689 = qJD(1) * t667;
t644 = -t666 * t689 + t671 * t688;
t624 = -mrSges(4,1) * t644 + mrSges(4,2) * t645;
t662 = qJD(2) + qJD(3);
t635 = mrSges(4,1) * t662 - mrSges(4,3) * t645;
t661 = qJDD(2) + qJDD(3);
t618 = qJD(3) * t644 + t649 * t671 + t650 * t666;
t654 = t668 * g(1) - t673 * g(2);
t627 = -t654 + (-t650 + t686) * pkin(1);
t586 = (-t644 * t662 - t618) * pkin(5) + (t645 * t662 - t617) * pkin(2) + t627;
t625 = -pkin(2) * t644 - pkin(5) * t645;
t660 = t662 ^ 2;
t594 = -pkin(2) * t660 + pkin(5) * t661 + t625 * t644 + t606;
t665 = sin(qJ(4));
t670 = cos(qJ(4));
t575 = t670 * t586 - t665 * t594;
t616 = qJDD(4) - t617;
t631 = -t645 * t665 + t662 * t670;
t632 = t645 * t670 + t662 * t665;
t573 = (t631 * t632 + t616) * pkin(3) + t575;
t576 = t665 * t586 + t670 * t594;
t640 = qJD(4) - t644;
t574 = (-t631 ^ 2 - t640 ^ 2) * pkin(3) + t576;
t664 = sin(qJ(5));
t669 = cos(qJ(5));
t571 = t573 * t669 - t574 * t664;
t596 = -qJD(4) * t632 - t618 * t665 + t661 * t670;
t597 = qJD(4) * t631 + t618 * t670 + t661 * t665;
t607 = t631 * t669 - t632 * t664;
t582 = qJD(5) * t607 + t596 * t664 + t597 * t669;
t608 = t631 * t664 + t632 * t669;
t592 = -mrSges(6,1) * t607 + mrSges(6,2) * t608;
t636 = qJD(5) + t640;
t598 = -mrSges(6,2) * t636 + mrSges(6,3) * t607;
t612 = qJDD(5) + t616;
t568 = m(6) * t571 + mrSges(6,1) * t612 - mrSges(6,3) * t582 - t592 * t608 + t598 * t636;
t572 = t573 * t664 + t574 * t669;
t581 = -qJD(5) * t608 + t596 * t669 - t597 * t664;
t599 = mrSges(6,1) * t636 - mrSges(6,3) * t608;
t569 = m(6) * t572 - mrSges(6,2) * t612 + mrSges(6,3) * t581 + t592 * t607 - t599 * t636;
t559 = t669 * t568 + t664 * t569;
t609 = -mrSges(5,1) * t631 + mrSges(5,2) * t632;
t619 = -mrSges(5,2) * t640 + mrSges(5,3) * t631;
t557 = m(5) * t575 + mrSges(5,1) * t616 - mrSges(5,3) * t597 - t609 * t632 + t619 * t640 + t559;
t620 = mrSges(5,1) * t640 - mrSges(5,3) * t632;
t558 = m(5) * t576 - mrSges(5,2) * t616 + mrSges(5,3) * t596 - t568 * t664 + t569 * t669 + t609 * t631 - t620 * t640;
t684 = -t557 * t665 + t670 * t558;
t546 = m(4) * t606 - mrSges(4,2) * t661 + mrSges(4,3) * t617 + t624 * t644 - t635 * t662 + t684;
t605 = t671 * t630 - t666 * t633;
t634 = -mrSges(4,2) * t662 + mrSges(4,3) * t644;
t593 = -t661 * pkin(2) - t660 * pkin(5) + t645 * t625 - t605;
t577 = (t632 * t640 - t596) * pkin(3) + t593;
t681 = m(6) * t577 - t581 * mrSges(6,1) + mrSges(6,2) * t582 - t607 * t598 + t599 * t608;
t677 = -m(5) * t593 + t596 * mrSges(5,1) - mrSges(5,2) * t597 + t631 * t619 - t620 * t632 - t681;
t563 = m(4) * t605 + mrSges(4,1) * t661 - mrSges(4,3) * t618 - t624 * t645 + t634 * t662 + t677;
t542 = t666 * t546 + t671 * t563;
t642 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t667 + Ifges(3,2) * t672) * qJD(1);
t643 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t667 + Ifges(3,4) * t672) * qJD(1);
t587 = Ifges(6,5) * t608 + Ifges(6,6) * t607 + Ifges(6,3) * t636;
t589 = Ifges(6,1) * t608 + Ifges(6,4) * t607 + Ifges(6,5) * t636;
t560 = -mrSges(6,1) * t577 + mrSges(6,3) * t572 + Ifges(6,4) * t582 + Ifges(6,2) * t581 + Ifges(6,6) * t612 - t587 * t608 + t589 * t636;
t588 = Ifges(6,4) * t608 + Ifges(6,2) * t607 + Ifges(6,6) * t636;
t561 = mrSges(6,2) * t577 - mrSges(6,3) * t571 + Ifges(6,1) * t582 + Ifges(6,4) * t581 + Ifges(6,5) * t612 + t587 * t607 - t588 * t636;
t600 = Ifges(5,5) * t632 + Ifges(5,6) * t631 + Ifges(5,3) * t640;
t602 = Ifges(5,1) * t632 + Ifges(5,4) * t631 + Ifges(5,5) * t640;
t551 = -mrSges(5,1) * t593 + mrSges(5,3) * t576 + Ifges(5,4) * t597 + Ifges(5,2) * t596 + Ifges(5,6) * t616 - pkin(3) * t681 + t669 * t560 + t664 * t561 - t632 * t600 + t640 * t602;
t601 = Ifges(5,4) * t632 + Ifges(5,2) * t631 + Ifges(5,6) * t640;
t553 = mrSges(5,2) * t593 - mrSges(5,3) * t575 + Ifges(5,1) * t597 + Ifges(5,4) * t596 + Ifges(5,5) * t616 - t560 * t664 + t561 * t669 + t600 * t631 - t601 * t640;
t622 = Ifges(4,4) * t645 + Ifges(4,2) * t644 + Ifges(4,6) * t662;
t623 = Ifges(4,1) * t645 + Ifges(4,4) * t644 + Ifges(4,5) * t662;
t679 = -mrSges(4,1) * t605 + mrSges(4,2) * t606 - Ifges(4,5) * t618 - Ifges(4,6) * t617 - Ifges(4,3) * t661 - pkin(2) * t677 - pkin(5) * t684 - t670 * t551 - t665 * t553 - t645 * t622 + t644 * t623;
t691 = mrSges(3,1) * t638 - mrSges(3,2) * t639 + Ifges(3,5) * t649 + Ifges(3,6) * t650 + Ifges(3,3) * qJDD(2) + pkin(1) * t542 + (t642 * t667 - t643 * t672) * qJD(1) - t679;
t648 = (-mrSges(3,1) * t672 + mrSges(3,2) * t667) * qJD(1);
t653 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t688;
t540 = m(3) * t638 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t649 + qJD(2) * t653 - t648 * t689 + t542;
t652 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t689;
t541 = m(3) * t639 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t650 - qJD(2) * t652 + t546 * t671 - t563 * t666 + t648 * t688;
t538 = m(2) * t655 - mrSges(2,1) * t674 - qJDD(1) * mrSges(2,2) - t540 * t667 + t541 * t672;
t548 = t670 * t557 + t665 * t558;
t678 = m(4) * t627 - t617 * mrSges(4,1) + t618 * mrSges(4,2) - t644 * t634 + t645 * t635 + t548;
t544 = qJDD(1) * mrSges(2,1) + t650 * mrSges(3,1) - t674 * mrSges(2,2) - t649 * mrSges(3,2) + (m(2) + m(3)) * t654 + (-t652 * t667 + t653 * t672) * qJD(1) - t678;
t690 = t668 * t538 + t673 * t544;
t685 = t673 * t538 - t544 * t668;
t621 = Ifges(4,5) * t645 + Ifges(4,6) * t644 + Ifges(4,3) * t662;
t535 = mrSges(4,2) * t627 - mrSges(4,3) * t605 + Ifges(4,1) * t618 + Ifges(4,4) * t617 + Ifges(4,5) * t661 - pkin(5) * t548 - t551 * t665 + t553 * t670 + t621 * t644 - t622 * t662;
t680 = -mrSges(6,1) * t571 + mrSges(6,2) * t572 - Ifges(6,5) * t582 - Ifges(6,6) * t581 - Ifges(6,3) * t612 - t608 * t588 + t607 * t589;
t675 = mrSges(5,1) * t575 - mrSges(5,2) * t576 + Ifges(5,5) * t597 + Ifges(5,6) * t596 + Ifges(5,3) * t616 + pkin(3) * t559 + t632 * t601 - t631 * t602 - t680;
t539 = -mrSges(4,1) * t627 + mrSges(4,3) * t606 + Ifges(4,4) * t618 + Ifges(4,2) * t617 + Ifges(4,6) * t661 - pkin(2) * t548 - t645 * t621 + t662 * t623 - t675;
t641 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t667 + Ifges(3,6) * t672) * qJD(1);
t531 = mrSges(3,1) * t654 + mrSges(3,3) * t639 + Ifges(3,4) * t649 + Ifges(3,2) * t650 + Ifges(3,6) * qJDD(2) - pkin(1) * t678 + qJD(2) * t643 + t666 * t535 + t671 * t539 - t641 * t689;
t533 = -mrSges(3,2) * t654 - mrSges(3,3) * t638 + Ifges(3,1) * t649 + Ifges(3,4) * t650 + Ifges(3,5) * qJDD(2) - qJD(2) * t642 + t535 * t671 - t539 * t666 + t641 * t688;
t682 = mrSges(2,1) * t654 - mrSges(2,2) * t655 + Ifges(2,3) * qJDD(1) + t672 * t531 + t667 * t533;
t534 = mrSges(2,1) * g(3) + mrSges(2,3) * t655 + t674 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - t691;
t529 = -mrSges(2,2) * g(3) - mrSges(2,3) * t654 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t674 - t531 * t667 + t533 * t672;
t1 = [-m(1) * g(1) + t685; -m(1) * g(2) + t690; t672 * t540 + t667 * t541 + (-m(1) - m(2)) * g(3); -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t690 + t673 * t529 - t668 * t534; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t685 + t668 * t529 + t673 * t534; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t682; t682; t691; -t679; t675; -t680;];
tauJB  = t1;
