% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:18
% EndTime: 2020-01-03 11:45:20
% DurationCPUTime: 1.83s
% Computational Cost: add. (21678->207), mult. (29571->247), div. (0->0), fcn. (14363->8), ass. (0->94)
t624 = Ifges(5,1) + Ifges(6,1);
t616 = Ifges(5,4) + Ifges(6,4);
t615 = Ifges(5,5) + Ifges(6,5);
t623 = Ifges(5,2) + Ifges(6,2);
t614 = Ifges(5,6) + Ifges(6,6);
t622 = Ifges(5,3) + Ifges(6,3);
t575 = qJD(1) + qJD(3);
t583 = sin(qJ(4));
t586 = cos(qJ(4));
t549 = (-mrSges(6,1) * t586 + mrSges(6,2) * t583) * t575;
t574 = qJDD(1) + qJDD(3);
t605 = qJD(4) * t575;
t600 = t586 * t605;
t551 = t583 * t574 + t600;
t585 = sin(qJ(1));
t588 = cos(qJ(1));
t566 = -t588 * g(2) - t585 * g(3);
t557 = qJDD(1) * pkin(1) + t566;
t565 = -t585 * g(2) + t588 * g(3);
t589 = qJD(1) ^ 2;
t558 = -t589 * pkin(1) + t565;
t581 = sin(pkin(8));
t582 = cos(pkin(8));
t535 = t582 * t557 - t581 * t558;
t532 = qJDD(1) * pkin(2) + t535;
t536 = t581 * t557 + t582 * t558;
t533 = -t589 * pkin(2) + t536;
t584 = sin(qJ(3));
t587 = cos(qJ(3));
t528 = t584 * t532 + t587 * t533;
t573 = t575 ^ 2;
t525 = -t573 * pkin(3) + t574 * pkin(7) + t528;
t580 = -g(1) + qJDD(2);
t568 = t586 * t580;
t604 = qJD(5) * t575;
t618 = pkin(4) * t573;
t518 = qJDD(4) * pkin(4) + t568 + (-t551 + t600) * qJ(5) + (t586 * t618 - t525 - 0.2e1 * t604) * t583;
t611 = t575 * t586;
t562 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t611;
t603 = m(6) * t518 + qJDD(4) * mrSges(6,1) + qJD(4) * t562;
t612 = t575 * t583;
t515 = -t551 * mrSges(6,3) - t549 * t612 + t603;
t522 = t586 * t525 + t583 * t580;
t552 = t586 * t574 - t583 * t605;
t559 = qJD(4) * pkin(4) - qJ(5) * t612;
t579 = t586 ^ 2;
t519 = t552 * qJ(5) - qJD(4) * t559 - t579 * t618 + 0.2e1 * t586 * t604 + t522;
t521 = -t583 * t525 + t568;
t607 = (t624 * t583 + t616 * t586) * t575 + t615 * qJD(4);
t608 = (t616 * t583 + t623 * t586) * t575 + t614 * qJD(4);
t621 = mrSges(5,1) * t521 + mrSges(6,1) * t518 - mrSges(5,2) * t522 - mrSges(6,2) * t519 + pkin(4) * t515 + t622 * qJDD(4) + t615 * t551 + t614 * t552 + (t608 * t583 - t607 * t586) * t575;
t617 = -mrSges(5,2) - mrSges(6,2);
t550 = (-mrSges(5,1) * t586 + mrSges(5,2) * t583) * t575;
t563 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t611;
t512 = m(5) * t521 + qJDD(4) * mrSges(5,1) + qJD(4) * t563 + (-t549 - t550) * t612 + (-mrSges(5,3) - mrSges(6,3)) * t551 + t603;
t602 = m(6) * t519 + t552 * mrSges(6,3) + t549 * t611;
t560 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t612;
t606 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t612 - t560;
t513 = m(5) * t522 + t552 * mrSges(5,3) + t606 * qJD(4) + t617 * qJDD(4) + t550 * t611 + t602;
t596 = -t583 * t512 + t586 * t513;
t502 = m(4) * t528 - t573 * mrSges(4,1) - t574 * mrSges(4,2) + t596;
t527 = t587 * t532 - t584 * t533;
t594 = -t574 * pkin(3) - t527;
t524 = -t573 * pkin(7) + t594;
t520 = t559 * t612 - t552 * pkin(4) + qJDD(5) + (-qJ(5) * t579 - pkin(7)) * t573 + t594;
t595 = -m(6) * t520 + t552 * mrSges(6,1) + t562 * t611;
t591 = -m(5) * t524 + t552 * mrSges(5,1) + t617 * t551 + t563 * t611 + t606 * t612 + t595;
t507 = m(4) * t527 + t574 * mrSges(4,1) - t573 * mrSges(4,2) + t591;
t495 = t584 * t502 + t587 * t507;
t492 = m(3) * t535 + qJDD(1) * mrSges(3,1) - t589 * mrSges(3,2) + t495;
t597 = t587 * t502 - t584 * t507;
t493 = m(3) * t536 - t589 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t597;
t598 = -t581 * t492 + t582 * t493;
t484 = m(2) * t565 - t589 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t598;
t487 = t582 * t492 + t581 * t493;
t485 = m(2) * t566 + qJDD(1) * mrSges(2,1) - t589 * mrSges(2,2) + t487;
t610 = t585 * t484 + t588 * t485;
t505 = t586 * t512 + t583 * t513;
t609 = (t615 * t583 + t614 * t586) * t575 + t622 * qJD(4);
t601 = m(4) * t580 + t505;
t599 = -t588 * t484 + t585 * t485;
t503 = m(3) * t580 + t601;
t514 = t551 * mrSges(6,2) + t560 * t612 - t595;
t497 = -mrSges(5,1) * t524 + mrSges(5,3) * t522 - mrSges(6,1) * t520 + mrSges(6,3) * t519 - pkin(4) * t514 + qJ(5) * t602 - t609 * t612 + t623 * t552 + t616 * t551 + (-qJ(5) * mrSges(6,2) + t614) * qJDD(4) + (-qJ(5) * t560 + t607) * qJD(4);
t499 = mrSges(5,2) * t524 + mrSges(6,2) * t520 - mrSges(5,3) * t521 - mrSges(6,3) * t518 - qJ(5) * t515 - t608 * qJD(4) + t615 * qJDD(4) + t624 * t551 + t616 * t552 + t609 * t611;
t593 = mrSges(4,1) * t527 - mrSges(4,2) * t528 + Ifges(4,3) * t574 + pkin(3) * t591 + pkin(7) * t596 + t586 * t497 + t583 * t499;
t590 = mrSges(2,1) * t566 + mrSges(3,1) * t535 - mrSges(2,2) * t565 - mrSges(3,2) * t536 + pkin(1) * t487 + pkin(2) * t495 + t593 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t488 = -mrSges(4,1) * t580 + mrSges(4,3) * t528 + t573 * Ifges(4,5) + Ifges(4,6) * t574 - pkin(3) * t505 - t621;
t480 = mrSges(4,2) * t580 - mrSges(4,3) * t527 + Ifges(4,5) * t574 - t573 * Ifges(4,6) - pkin(7) * t505 - t583 * t497 + t586 * t499;
t479 = mrSges(3,2) * t580 - mrSges(3,3) * t535 + Ifges(3,5) * qJDD(1) - t589 * Ifges(3,6) - pkin(6) * t495 + t587 * t480 - t584 * t488;
t478 = -mrSges(3,1) * t580 + mrSges(3,3) * t536 + t589 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t601 + pkin(6) * t597 + t584 * t480 + t587 * t488;
t477 = -mrSges(2,2) * g(1) - mrSges(2,3) * t566 + Ifges(2,5) * qJDD(1) - t589 * Ifges(2,6) - qJ(2) * t487 - t581 * t478 + t582 * t479;
t476 = mrSges(2,1) * g(1) + mrSges(2,3) * t565 + t589 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t503 + qJ(2) * t598 + t582 * t478 + t581 * t479;
t1 = [(-m(1) - m(2)) * g(1) + t503; -m(1) * g(2) + t610; -m(1) * g(3) + t599; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t590; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t599 + t588 * t476 + t585 * t477; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t610 + t585 * t476 - t588 * t477; t590; t503; t593; t621; t514;];
tauJB = t1;
