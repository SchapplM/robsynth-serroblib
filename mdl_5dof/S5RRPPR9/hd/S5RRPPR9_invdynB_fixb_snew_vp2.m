% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPPR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPPR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR9_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR9_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR9_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:51
% EndTime: 2019-12-31 19:40:53
% DurationCPUTime: 1.63s
% Computational Cost: add. (9507->276), mult. (19682->323), div. (0->0), fcn. (9106->6), ass. (0->110)
t682 = Ifges(3,1) + Ifges(4,1) + Ifges(5,2);
t661 = Ifges(3,4) - Ifges(4,5) + Ifges(5,4);
t660 = Ifges(3,5) + Ifges(4,4) + Ifges(5,6);
t681 = Ifges(3,2) + Ifges(5,1) + Ifges(4,3);
t659 = Ifges(3,6) - Ifges(4,6) + Ifges(5,5);
t680 = -Ifges(3,3) - Ifges(4,2) - Ifges(5,3);
t635 = sin(qJ(1));
t638 = cos(qJ(1));
t615 = -t638 * g(1) - t635 * g(2);
t640 = qJD(1) ^ 2;
t580 = -t640 * pkin(1) + qJDD(1) * pkin(6) + t615;
t634 = sin(qJ(2));
t637 = cos(qJ(2));
t563 = -t634 * g(3) + t637 * t580;
t597 = (-pkin(2) * t637 - qJ(3) * t634) * qJD(1);
t666 = qJD(1) * t637;
t679 = qJDD(2) * qJ(3) + t597 * t666 + t563;
t663 = qJD(1) * qJD(2);
t654 = t637 * t663;
t602 = t634 * qJDD(1) + t654;
t662 = qJD(1) * qJD(4);
t678 = -0.2e1 * t634 * t662 + (-t602 + t654) * qJ(4);
t677 = -2 * qJD(3);
t676 = 2 * qJD(3);
t675 = -pkin(2) - pkin(7);
t674 = pkin(3) + pkin(7);
t639 = qJD(2) ^ 2;
t673 = t639 * pkin(2);
t672 = mrSges(3,3) + mrSges(4,2);
t670 = t637 ^ 2 * t640;
t669 = t637 * t640;
t562 = -t637 * g(3) - t634 * t580;
t598 = (-mrSges(4,1) * t637 - mrSges(4,3) * t634) * qJD(1);
t599 = (-mrSges(3,1) * t637 + mrSges(3,2) * t634) * qJD(1);
t611 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t666;
t612 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t666;
t664 = t634 * qJD(1);
t651 = t597 * t664 + qJDD(3) - t562;
t551 = -qJDD(2) * pkin(2) - t639 * qJ(3) + t651;
t613 = mrSges(4,2) * t666 + qJD(2) * mrSges(4,3);
t547 = (-t634 * t669 - qJDD(2)) * pkin(3) + t551 + t678;
t600 = (mrSges(5,1) * t634 - mrSges(5,2) * t637) * qJD(1);
t655 = t634 * t663;
t603 = t637 * qJDD(1) - t655;
t607 = -qJD(2) * pkin(3) - qJ(4) * t664;
t614 = t635 * g(1) - t638 * g(2);
t579 = -qJDD(1) * pkin(1) - t640 * pkin(6) - t614;
t648 = -t603 * pkin(2) + t579 + (-t602 - t654) * qJ(3);
t643 = -qJ(4) * t670 + qJDD(4) - t648 + (t607 + t676) * t664;
t540 = t674 * t603 + t643 + (pkin(4) * t637 + t675 * t634) * t663 + t602 * pkin(4);
t601 = (pkin(4) * t634 + pkin(7) * t637) * qJD(1);
t543 = (-pkin(4) - qJ(3)) * t639 + (-pkin(3) * t669 - qJD(1) * t601) * t634 + (-pkin(2) - t674) * qJDD(2) + t651 + t678;
t633 = sin(qJ(5));
t636 = cos(qJ(5));
t538 = t636 * t540 - t633 * t543;
t595 = -t636 * qJD(2) + t633 * t666;
t558 = t595 * qJD(5) - t633 * qJDD(2) - t636 * t603;
t596 = -t633 * qJD(2) - t636 * t666;
t559 = -t595 * mrSges(6,1) + t596 * mrSges(6,2);
t618 = qJD(5) + t664;
t560 = -t618 * mrSges(6,2) + t595 * mrSges(6,3);
t593 = qJDD(5) + t602;
t536 = m(6) * t538 + t593 * mrSges(6,1) - t558 * mrSges(6,3) - t596 * t559 + t618 * t560;
t539 = t633 * t540 + t636 * t543;
t557 = -t596 * qJD(5) - t636 * qJDD(2) + t633 * t603;
t561 = t618 * mrSges(6,1) - t596 * mrSges(6,3);
t537 = m(6) * t539 - t593 * mrSges(6,2) + t557 * mrSges(6,3) + t595 * t559 - t618 * t561;
t667 = -t633 * t536 + t636 * t537;
t650 = -m(5) * t547 + t600 * t664 - t667;
t646 = -m(4) * t551 + qJDD(2) * mrSges(4,1) + qJD(2) * t613 + t650;
t523 = m(3) * t562 + (mrSges(3,1) - mrSges(5,2)) * qJDD(2) + (-t611 + t612) * qJD(2) + (-t598 - t599) * t664 + (mrSges(5,3) - t672) * t602 + t646;
t609 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t664;
t622 = qJD(2) * t676;
t550 = t622 - t673 + t679;
t610 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t664;
t645 = pkin(3) * t670 + t603 * qJ(4) - t679;
t546 = 0.2e1 * t637 * t662 + t673 + (t677 - t607) * qJD(2) + t645;
t608 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t664;
t542 = qJDD(2) * pkin(4) + qJD(2) * t607 + t622 + t675 * t639 + (-0.2e1 * qJD(4) - t601) * t666 - t645;
t647 = -m(6) * t542 + t557 * mrSges(6,1) - t558 * mrSges(6,2) + t595 * t560 - t596 * t561;
t644 = -m(5) * t546 + qJDD(2) * mrSges(5,1) - t603 * mrSges(5,3) + qJD(2) * t608 - t647;
t642 = m(4) * t550 + qJDD(2) * mrSges(4,3) + qJD(2) * t610 + t598 * t666 + t644;
t529 = t672 * t603 + t642 + (t599 - t600) * t666 - qJDD(2) * mrSges(3,2) - qJD(2) * t609 + m(3) * t563;
t652 = -t634 * t523 + t637 * t529;
t518 = m(2) * t615 - t640 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t652;
t548 = (pkin(2) * qJD(2) + t677) * t664 + t648;
t526 = t636 * t536 + t633 * t537;
t545 = -pkin(2) * t655 + t603 * pkin(3) + t643;
t649 = -m(5) * t545 - t602 * mrSges(5,1) + t603 * mrSges(5,2) - t608 * t664 + t611 * t666 - t526;
t524 = m(4) * t548 - t603 * mrSges(4,1) - t602 * mrSges(4,3) - t610 * t664 - t613 * t666 + t649;
t641 = -m(3) * t579 + t603 * mrSges(3,1) - t602 * mrSges(3,2) - t609 * t664 + t612 * t666 - t524;
t521 = m(2) * t614 + qJDD(1) * mrSges(2,1) - t640 * mrSges(2,2) + t641;
t668 = t635 * t518 + t638 * t521;
t519 = t637 * t523 + t634 * t529;
t665 = qJD(2) * t611;
t658 = t680 * qJD(2) + (-t660 * t634 - t659 * t637) * qJD(1);
t657 = -t659 * qJD(2) + (-t661 * t634 - t681 * t637) * qJD(1);
t656 = t660 * qJD(2) + (t682 * t634 + t661 * t637) * qJD(1);
t653 = t638 * t518 - t635 * t521;
t554 = Ifges(6,1) * t596 + Ifges(6,4) * t595 + Ifges(6,5) * t618;
t553 = Ifges(6,4) * t596 + Ifges(6,2) * t595 + Ifges(6,6) * t618;
t552 = Ifges(6,5) * t596 + Ifges(6,6) * t595 + Ifges(6,3) * t618;
t531 = mrSges(6,2) * t542 - mrSges(6,3) * t538 + Ifges(6,1) * t558 + Ifges(6,4) * t557 + Ifges(6,5) * t593 + t595 * t552 - t618 * t553;
t530 = -mrSges(6,1) * t542 + mrSges(6,3) * t539 + Ifges(6,4) * t558 + Ifges(6,2) * t557 + Ifges(6,6) * t593 - t596 * t552 + t618 * t554;
t525 = qJDD(2) * mrSges(5,2) - t602 * mrSges(5,3) - t650 + t665;
t515 = Ifges(6,3) * t593 - t595 * t554 + t596 * t553 + mrSges(3,2) * t579 + mrSges(4,2) * t551 + Ifges(6,6) * t557 + Ifges(6,5) * t558 - mrSges(3,3) * t562 + mrSges(5,1) * t545 - mrSges(5,3) * t547 - mrSges(4,3) * t548 + mrSges(6,1) * t538 - mrSges(6,2) * t539 - qJ(4) * t525 + pkin(4) * t526 - qJ(3) * t524 + t661 * t603 + t682 * t602 + t660 * qJDD(2) + t657 * qJD(2) - t658 * t666;
t514 = -t636 * t531 - qJ(4) * t644 + t633 * t530 - pkin(3) * t649 - mrSges(3,1) * t579 + mrSges(3,3) * t563 - mrSges(5,2) * t545 + mrSges(5,3) * t546 - mrSges(4,1) * t548 + mrSges(4,2) * t550 + pkin(7) * t526 - pkin(2) * t524 + t681 * t603 + t661 * t602 + t659 * qJDD(2) + t656 * qJD(2) + (qJ(4) * t600 * t637 + t658 * t634) * qJD(1);
t513 = Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) - qJ(3) * t642 - pkin(2) * (t646 - t665) + t636 * t530 + t640 * Ifges(2,5) + t633 * t531 + mrSges(2,3) * t615 + pkin(4) * t647 - mrSges(4,3) * t550 + mrSges(4,1) * t551 - mrSges(3,1) * t562 + mrSges(3,2) * t563 + mrSges(5,1) * t546 - mrSges(5,2) * t547 + pkin(7) * t667 + pkin(3) * t525 - pkin(1) * t519 + (-qJ(3) * mrSges(4,2) - t659) * t603 + (pkin(2) * mrSges(5,2) + t680) * qJDD(2) + (-pkin(2) * (-mrSges(4,2) + mrSges(5,3)) - t660) * t602 + ((qJ(3) * t600 + t656) * t637 + (pkin(2) * t598 + t657) * t634) * qJD(1);
t512 = -mrSges(2,2) * g(3) - mrSges(2,3) * t614 + Ifges(2,5) * qJDD(1) - t640 * Ifges(2,6) - pkin(6) * t519 - t634 * t514 + t637 * t515;
t1 = [-m(1) * g(1) + t653; -m(1) * g(2) + t668; (-m(1) - m(2)) * g(3) + t519; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t668 + t638 * t512 - t635 * t513; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t653 + t635 * t512 + t638 * t513; -mrSges(1,1) * g(2) + mrSges(2,1) * t614 + mrSges(1,2) * g(1) - mrSges(2,2) * t615 + Ifges(2,3) * qJDD(1) + pkin(1) * t641 + pkin(6) * t652 + t637 * t514 + t634 * t515;];
tauB = t1;
