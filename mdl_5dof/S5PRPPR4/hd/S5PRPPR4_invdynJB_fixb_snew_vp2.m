% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPPR4
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:53
% EndTime: 2019-12-31 17:36:54
% DurationCPUTime: 1.58s
% Computational Cost: add. (12293->202), mult. (26650->251), div. (0->0), fcn. (16160->8), ass. (0->97)
t595 = sin(pkin(7));
t597 = cos(pkin(7));
t574 = t595 * g(1) - t597 * g(2);
t575 = -t597 * g(1) - t595 * g(2);
t599 = sin(qJ(2));
t601 = cos(qJ(2));
t557 = t599 * t574 + t601 * t575;
t602 = qJD(2) ^ 2;
t642 = -t602 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(2) * qJD(3)) + t557;
t594 = sin(pkin(8));
t596 = cos(pkin(8));
t641 = (Ifges(4,6) - Ifges(5,6)) * t596 + (Ifges(5,4) + Ifges(4,5)) * t594;
t589 = t594 ^ 2;
t590 = t596 ^ 2;
t630 = t590 * t602;
t638 = t589 * t602 + t630;
t635 = Ifges(4,4) - Ifges(5,5);
t637 = t635 * t594;
t593 = -g(3) + qJDD(1);
t537 = t596 * t593 - t642 * t594;
t636 = Ifges(4,1) + Ifges(5,1);
t634 = Ifges(4,2) + Ifges(5,3);
t633 = mrSges(4,2) * t594;
t632 = qJ(4) * t594;
t629 = t602 * qJ(3);
t570 = (-mrSges(5,1) * t596 - mrSges(5,3) * t594) * qJD(2);
t571 = (-mrSges(4,1) * t596 + t633) * qJD(2);
t612 = -pkin(3) * t596 - t632;
t569 = t612 * qJD(2);
t626 = t594 * qJD(2);
t532 = t569 * t626 + qJDD(4) - t537;
t528 = (-pkin(4) * t596 * t602 - pkin(6) * qJDD(2)) * t594 + t532;
t538 = t594 * t593 + t596 * t642;
t625 = t596 * qJD(2);
t534 = t569 * t625 + t538;
t622 = qJDD(2) * t596;
t529 = -pkin(4) * t630 - pkin(6) * t622 + t534;
t598 = sin(qJ(5));
t600 = cos(qJ(5));
t526 = t600 * t528 - t598 * t529;
t609 = -t594 * t598 - t596 * t600;
t563 = t609 * qJD(2);
t610 = t594 * t600 - t596 * t598;
t564 = t610 * qJD(2);
t546 = -t563 * mrSges(6,1) + t564 * mrSges(6,2);
t553 = t563 * qJD(5) + t610 * qJDD(2);
t558 = -qJD(5) * mrSges(6,2) + t563 * mrSges(6,3);
t523 = m(6) * t526 + qJDD(5) * mrSges(6,1) - t553 * mrSges(6,3) + qJD(5) * t558 - t564 * t546;
t527 = t598 * t528 + t600 * t529;
t552 = -t564 * qJD(5) + t609 * qJDD(2);
t559 = qJD(5) * mrSges(6,1) - t564 * mrSges(6,3);
t524 = m(6) * t527 - qJDD(5) * mrSges(6,2) + t552 * mrSges(6,3) - qJD(5) * t559 + t563 * t546;
t514 = t600 * t523 + t598 * t524;
t606 = m(5) * t532 + t514;
t511 = m(4) * t537 + ((-mrSges(5,2) - mrSges(4,3)) * qJDD(2) + (-t570 - t571) * qJD(2)) * t594 - t606;
t615 = -t598 * t523 + t600 * t524;
t608 = m(5) * t534 + mrSges(5,2) * t622 + t570 * t625 + t615;
t512 = m(4) * t538 + (qJDD(2) * mrSges(4,3) + qJD(2) * t571) * t596 + t608;
t616 = -t594 * t511 + t596 * t512;
t505 = m(3) * t557 - t602 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t616;
t556 = t601 * t574 - t599 * t575;
t619 = -qJDD(3) + t556;
t607 = -0.2e1 * qJD(4) * t626 - t619;
t536 = -t629 + (-pkin(2) + t612) * qJDD(2) + t607;
t531 = (qJ(3) + (-t589 - t590) * pkin(6)) * t602 + (t632 + pkin(2) + (pkin(3) + pkin(4)) * t596) * qJDD(2) - t607;
t611 = -m(6) * t531 + t552 * mrSges(6,1) - t553 * mrSges(6,2) + t563 * t558 - t564 * t559;
t623 = qJDD(2) * t594;
t521 = m(5) * t536 - mrSges(5,1) * t622 - t638 * mrSges(5,2) - mrSges(5,3) * t623 + t611;
t551 = -qJDD(2) * pkin(2) - t619 - t629;
t603 = -m(4) * t551 + mrSges(4,1) * t622 + t638 * mrSges(4,3) - t521;
t518 = t603 + (mrSges(3,1) - t633) * qJDD(2) - t602 * mrSges(3,2) + m(3) * t556;
t502 = t599 * t505 + t601 * t518;
t500 = m(2) * t574 + t502;
t617 = t601 * t505 - t599 * t518;
t501 = m(2) * t575 + t617;
t628 = t597 * t500 + t595 * t501;
t507 = t596 * t511 + t594 * t512;
t627 = t641 * qJD(2);
t621 = m(3) * t593 + t507;
t618 = -t595 * t500 + t597 * t501;
t613 = m(2) * t593 + t621;
t539 = Ifges(6,5) * t564 + Ifges(6,6) * t563 + Ifges(6,3) * qJD(5);
t541 = Ifges(6,1) * t564 + Ifges(6,4) * t563 + Ifges(6,5) * qJD(5);
t515 = -mrSges(6,1) * t531 + mrSges(6,3) * t527 + Ifges(6,4) * t553 + Ifges(6,2) * t552 + Ifges(6,6) * qJDD(5) + qJD(5) * t541 - t564 * t539;
t540 = Ifges(6,4) * t564 + Ifges(6,2) * t563 + Ifges(6,6) * qJD(5);
t516 = mrSges(6,2) * t531 - mrSges(6,3) * t526 + Ifges(6,1) * t553 + Ifges(6,4) * t552 + Ifges(6,5) * qJDD(5) - qJD(5) * t540 + t563 * t539;
t494 = -mrSges(4,1) * t551 + mrSges(4,3) * t538 - mrSges(5,1) * t536 + mrSges(5,2) * t534 - t598 * t516 - t600 * t515 - pkin(4) * t611 - pkin(6) * t615 - pkin(3) * t521 - t627 * t626 + (t634 * t596 + t637) * qJDD(2);
t496 = mrSges(4,2) * t551 + mrSges(5,2) * t532 - mrSges(4,3) * t537 - mrSges(5,3) * t536 - pkin(6) * t514 - qJ(4) * t521 - t598 * t515 + t600 * t516 + t627 * t625 + (t636 * t594 + t635 * t596) * qJDD(2);
t520 = mrSges(4,2) * t623 - t603;
t605 = mrSges(3,1) * t556 - mrSges(3,2) * t557 + Ifges(3,3) * qJDD(2) - pkin(2) * t520 + qJ(3) * t616 + t596 * t494 + t594 * t496;
t604 = mrSges(6,1) * t526 - mrSges(6,2) * t527 + Ifges(6,5) * t553 + Ifges(6,6) * t552 + Ifges(6,3) * qJDD(5) + t564 * t540 - t563 * t541;
t513 = (qJDD(2) * mrSges(5,2) + qJD(2) * t570) * t594 + t606;
t492 = (Ifges(3,6) - t641) * qJDD(2) - qJ(4) * t608 - mrSges(3,1) * t593 + mrSges(3,3) * t557 - mrSges(5,3) * t534 - mrSges(4,1) * t537 + mrSges(4,2) * t538 + mrSges(5,1) * t532 + pkin(3) * t513 + pkin(4) * t514 - pkin(2) * t507 + t604 + (t635 * t590 + (-t637 + (-t634 + t636) * t596) * t594 + Ifges(3,5)) * t602;
t491 = mrSges(3,2) * t593 - mrSges(3,3) * t556 + Ifges(3,5) * qJDD(2) - t602 * Ifges(3,6) - qJ(3) * t507 - t594 * t494 + t596 * t496;
t490 = mrSges(2,2) * t593 - mrSges(2,3) * t574 - pkin(5) * t502 + t601 * t491 - t599 * t492;
t489 = -mrSges(2,1) * t593 + mrSges(2,3) * t575 - pkin(1) * t621 + pkin(5) * t617 + t599 * t491 + t601 * t492;
t1 = [-m(1) * g(1) + t618; -m(1) * g(2) + t628; -m(1) * g(3) + t613; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t628 - t595 * t489 + t597 * t490; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t618 + t597 * t489 + t595 * t490; -mrSges(1,1) * g(2) + mrSges(2,1) * t574 + mrSges(1,2) * g(1) - mrSges(2,2) * t575 + pkin(1) * t502 + t605; t613; t605; t520; t513; t604;];
tauJB = t1;
