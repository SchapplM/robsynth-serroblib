% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:48
% EndTime: 2019-12-31 20:57:51
% DurationCPUTime: 2.15s
% Computational Cost: add. (16947->269), mult. (35448->316), div. (0->0), fcn. (21482->6), ass. (0->102)
t598 = sin(qJ(3));
t601 = cos(qJ(2));
t625 = qJD(1) * t601;
t599 = sin(qJ(2));
t626 = qJD(1) * t599;
t634 = cos(qJ(3));
t571 = t598 * t626 - t634 * t625;
t624 = qJD(1) * qJD(2);
t580 = t599 * qJDD(1) + t601 * t624;
t581 = t601 * qJDD(1) - t599 * t624;
t535 = -t571 * qJD(3) + t634 * t580 + t598 * t581;
t584 = qJD(2) * pkin(2) - pkin(7) * t626;
t597 = t601 ^ 2;
t603 = qJD(1) ^ 2;
t600 = sin(qJ(1));
t602 = cos(qJ(1));
t585 = t600 * g(1) - t602 * g(2);
t610 = qJDD(1) * pkin(1) + t585;
t536 = -t581 * pkin(2) + t584 * t626 - (pkin(7) * t597 + pkin(6)) * t603 - t610;
t596 = qJD(2) + qJD(3);
t631 = t571 * t596;
t641 = t536 + (-t535 + t631) * qJ(4);
t640 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t623 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t622 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t639 = Ifges(4,2) + Ifges(5,3) + Ifges(6,2);
t638 = -Ifges(5,2) - Ifges(4,3) - Ifges(6,3);
t621 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t572 = (t598 * t601 + t634 * t599) * qJD(1);
t636 = -0.2e1 * t572;
t635 = 2 * qJD(4);
t633 = pkin(2) * t603;
t632 = -mrSges(4,3) - mrSges(5,2);
t559 = t596 * mrSges(6,2) + t571 * mrSges(6,3);
t630 = t596 * t559;
t586 = -t602 * g(1) - t600 * g(2);
t574 = -t603 * pkin(1) + qJDD(1) * pkin(6) + t586;
t629 = t599 * t574;
t521 = qJDD(2) * pkin(2) - t580 * pkin(7) - t629 + (pkin(7) * t624 + t599 * t633 - g(3)) * t601;
t558 = -t599 * g(3) + t601 * t574;
t522 = t581 * pkin(7) - qJD(2) * t584 - t597 * t633 + t558;
t516 = t634 * t521 - t598 * t522;
t560 = -t596 * mrSges(4,2) - t571 * mrSges(4,3);
t595 = qJDD(2) + qJDD(3);
t551 = t571 * pkin(3) - t572 * qJ(4);
t594 = t596 ^ 2;
t515 = -t595 * pkin(3) - t594 * qJ(4) + t572 * t551 + qJDD(4) - t516;
t565 = -t571 * mrSges(5,2) + t596 * mrSges(5,3);
t508 = qJD(5) * t636 + (-t535 - t631) * qJ(5) + (t571 * t572 - t595) * pkin(4) + t515;
t553 = -t571 * mrSges(6,1) + t572 * mrSges(6,2);
t611 = -m(6) * t508 + t535 * mrSges(6,3) + t572 * t553;
t606 = -m(5) * t515 + t595 * mrSges(5,1) + t596 * t565 + t611;
t552 = t571 * mrSges(5,1) - t572 * mrSges(5,3);
t627 = -t571 * mrSges(4,1) - t572 * mrSges(4,2) - t552;
t500 = m(4) * t516 + (t559 + t560) * t596 + (mrSges(4,1) + mrSges(6,1)) * t595 + t627 * t572 + t632 * t535 + t606;
t517 = t598 * t521 + t634 * t522;
t534 = t572 * qJD(3) + t598 * t580 - t634 * t581;
t562 = -t596 * mrSges(6,1) - t572 * mrSges(6,3);
t563 = t596 * mrSges(4,1) - t572 * mrSges(4,3);
t514 = -t594 * pkin(3) + t595 * qJ(4) - t571 * t551 + t596 * t635 + t517;
t564 = -t596 * mrSges(5,1) + t572 * mrSges(5,2);
t561 = -t596 * pkin(4) - t572 * qJ(5);
t567 = t571 ^ 2;
t510 = -t567 * pkin(4) + t534 * qJ(5) + 0.2e1 * qJD(5) * t571 + t596 * t561 + t514;
t620 = m(6) * t510 + t534 * mrSges(6,3) + t571 * t553;
t608 = m(5) * t514 + t595 * mrSges(5,3) + t596 * t564 + t620;
t503 = m(4) * t517 + (t562 - t563) * t596 + (-mrSges(4,2) + mrSges(6,2)) * t595 + t627 * t571 + t632 * t534 + t608;
t496 = t634 * t500 + t598 * t503;
t557 = -t601 * g(3) - t629;
t579 = (-mrSges(3,1) * t601 + mrSges(3,2) * t599) * qJD(1);
t583 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t625;
t494 = m(3) * t557 + qJDD(2) * mrSges(3,1) - t580 * mrSges(3,3) + qJD(2) * t583 - t579 * t626 + t496;
t582 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t626;
t613 = -t598 * t500 + t634 * t503;
t495 = m(3) * t558 - qJDD(2) * mrSges(3,2) + t581 * mrSges(3,3) - qJD(2) * t582 + t579 * t625 + t613;
t614 = -t599 * t494 + t601 * t495;
t487 = m(2) * t586 - t603 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t614;
t573 = -t603 * pkin(6) - t610;
t512 = qJD(4) * t636 + (t572 * t596 + t534) * pkin(3) + t641;
t507 = -t567 * qJ(5) + qJDD(5) + (-pkin(3) - pkin(4)) * t534 + (-pkin(3) * t596 + t561 + t635) * t572 - t641;
t609 = -m(6) * t507 + t534 * mrSges(6,1) - t535 * mrSges(6,2) + t571 * t559 - t572 * t562;
t504 = m(5) * t512 + t534 * mrSges(5,1) - t535 * mrSges(5,3) - t572 * t564 + t571 * t565 + t609;
t605 = m(4) * t536 + t534 * mrSges(4,1) + t535 * mrSges(4,2) + t571 * t560 + t572 * t563 + t504;
t604 = -m(3) * t573 + t581 * mrSges(3,1) - t580 * mrSges(3,2) - t582 * t626 + t583 * t625 - t605;
t498 = m(2) * t585 + qJDD(1) * mrSges(2,1) - t603 * mrSges(2,2) + t604;
t628 = t600 * t487 + t602 * t498;
t488 = t601 * t494 + t599 * t495;
t619 = t621 * t571 - t622 * t572 + t638 * t596;
t618 = t639 * t571 - t623 * t572 - t621 * t596;
t617 = t623 * t571 - t640 * t572 - t622 * t596;
t615 = t602 * t487 - t600 * t498;
t570 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t599 + Ifges(3,4) * t601) * qJD(1);
t569 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t599 + Ifges(3,2) * t601) * qJD(1);
t568 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t599 + Ifges(3,6) * t601) * qJD(1);
t505 = -t595 * mrSges(6,1) - t611 - t630;
t490 = mrSges(4,2) * t536 + mrSges(5,2) * t515 + mrSges(6,2) * t507 - mrSges(4,3) * t516 - mrSges(5,3) * t512 - mrSges(6,3) * t508 - qJ(4) * t504 - qJ(5) * t505 - t623 * t534 + t640 * t535 + t619 * t571 + t622 * t595 + t618 * t596;
t489 = -mrSges(4,1) * t536 + mrSges(4,3) * t517 - mrSges(5,1) * t512 + mrSges(5,2) * t514 + mrSges(6,1) * t507 - mrSges(6,3) * t510 - pkin(4) * t609 - qJ(5) * t620 - pkin(3) * t504 + (-qJ(5) * t562 - t617) * t596 + (-qJ(5) * mrSges(6,2) + t621) * t595 + t619 * t572 + t623 * t535 - t639 * t534;
t484 = mrSges(3,2) * t573 - mrSges(3,3) * t557 + Ifges(3,1) * t580 + Ifges(3,4) * t581 + Ifges(3,5) * qJDD(2) - pkin(7) * t496 - qJD(2) * t569 - t598 * t489 + t634 * t490 + t568 * t625;
t483 = -mrSges(3,1) * t573 + mrSges(3,3) * t558 + Ifges(3,4) * t580 + Ifges(3,2) * t581 + Ifges(3,6) * qJDD(2) - pkin(2) * t605 + pkin(7) * t613 + qJD(2) * t570 + t634 * t489 + t598 * t490 - t568 * t626;
t482 = (-pkin(3) * mrSges(6,1) - qJ(4) * mrSges(6,2) + t638) * t595 + Ifges(2,6) * qJDD(1) - Ifges(3,3) * qJDD(2) + (-t599 * t569 + t601 * t570) * qJD(1) + mrSges(2,1) * g(3) - pkin(3) * (t606 + t630) + (qJ(4) * mrSges(5,2) + t621) * t534 + (pkin(3) * mrSges(5,2) - t622) * t535 + (qJ(4) * t552 + t617) * t571 + (pkin(3) * t552 + t618) * t572 - qJ(4) * (t596 * t562 + t608) + t603 * Ifges(2,5) - Ifges(3,5) * t580 - Ifges(3,6) * t581 + mrSges(2,3) * t586 - mrSges(3,1) * t557 + mrSges(3,2) * t558 + mrSges(5,1) * t515 - mrSges(4,1) * t516 + mrSges(4,2) * t517 - mrSges(6,2) * t510 - mrSges(5,3) * t514 + pkin(4) * t505 + mrSges(6,1) * t508 - pkin(2) * t496 - pkin(1) * t488;
t481 = -mrSges(2,2) * g(3) - mrSges(2,3) * t585 + Ifges(2,5) * qJDD(1) - t603 * Ifges(2,6) - pkin(6) * t488 - t599 * t483 + t601 * t484;
t1 = [-m(1) * g(1) + t615; -m(1) * g(2) + t628; (-m(1) - m(2)) * g(3) + t488; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t628 + t602 * t481 - t600 * t482; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t615 + t600 * t481 + t602 * t482; -mrSges(1,1) * g(2) + mrSges(2,1) * t585 + mrSges(1,2) * g(1) - mrSges(2,2) * t586 + Ifges(2,3) * qJDD(1) + pkin(1) * t604 + pkin(6) * t614 + t601 * t483 + t599 * t484;];
tauB = t1;
