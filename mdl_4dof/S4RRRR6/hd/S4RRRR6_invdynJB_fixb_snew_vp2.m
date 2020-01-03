% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR6_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR6_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:44
% EndTime: 2019-12-31 17:29:49
% DurationCPUTime: 4.62s
% Computational Cost: add. (52596->253), mult. (113868->338), div. (0->0), fcn. (84383->10), ass. (0->111)
t583 = sin(pkin(4));
t587 = sin(qJ(2));
t591 = cos(qJ(2));
t603 = qJD(1) * qJD(2);
t571 = (-qJDD(1) * t591 + t587 * t603) * t583;
t613 = pkin(6) * t583;
t584 = cos(pkin(4));
t612 = t584 * g(3);
t611 = t583 * t587;
t610 = t583 * t591;
t609 = t584 * t587;
t608 = t584 * t591;
t588 = sin(qJ(1));
t592 = cos(qJ(1));
t575 = t588 * g(1) - t592 * g(2);
t593 = qJD(1) ^ 2;
t566 = qJDD(1) * pkin(1) + t593 * t613 + t575;
t576 = -t592 * g(1) - t588 * g(2);
t567 = -t593 * pkin(1) + qJDD(1) * t613 + t576;
t606 = t566 * t609 + t591 * t567;
t545 = -g(3) * t611 + t606;
t579 = t584 * qJD(1) + qJD(2);
t605 = qJD(1) * t583;
t602 = t587 * t605;
t564 = t579 * mrSges(3,1) - mrSges(3,3) * t602;
t568 = (-mrSges(3,1) * t591 + mrSges(3,2) * t587) * t605;
t578 = t584 * qJDD(1) + qJDD(2);
t569 = (-pkin(2) * t591 - pkin(7) * t587) * t605;
t577 = t579 ^ 2;
t604 = qJD(1) * t591;
t531 = -t577 * pkin(2) + t578 * pkin(7) + (-g(3) * t587 + t569 * t604) * t583 + t606;
t570 = (qJDD(1) * t587 + t591 * t603) * t583;
t532 = t571 * pkin(2) - t570 * pkin(7) - t612 + (-t566 + (pkin(2) * t587 - pkin(7) * t591) * t579 * qJD(1)) * t583;
t586 = sin(qJ(3));
t590 = cos(qJ(3));
t520 = t590 * t531 + t586 * t532;
t559 = t590 * t579 - t586 * t602;
t560 = t586 * t579 + t590 * t602;
t547 = -t559 * pkin(3) - t560 * pkin(8);
t563 = qJDD(3) + t571;
t601 = t583 * t604;
t574 = qJD(3) - t601;
t573 = t574 ^ 2;
t517 = -t573 * pkin(3) + t563 * pkin(8) + t559 * t547 + t520;
t544 = -g(3) * t610 + t566 * t608 - t587 * t567;
t530 = -t578 * pkin(2) - t577 * pkin(7) + t569 * t602 - t544;
t542 = -t560 * qJD(3) - t586 * t570 + t590 * t578;
t543 = t559 * qJD(3) + t590 * t570 + t586 * t578;
t518 = (-t559 * t574 - t543) * pkin(8) + (t560 * t574 - t542) * pkin(3) + t530;
t585 = sin(qJ(4));
t589 = cos(qJ(4));
t514 = -t585 * t517 + t589 * t518;
t548 = -t585 * t560 + t589 * t574;
t523 = t548 * qJD(4) + t589 * t543 + t585 * t563;
t549 = t589 * t560 + t585 * t574;
t533 = -t548 * mrSges(5,1) + t549 * mrSges(5,2);
t558 = qJD(4) - t559;
t534 = -t558 * mrSges(5,2) + t548 * mrSges(5,3);
t540 = qJDD(4) - t542;
t511 = m(5) * t514 + t540 * mrSges(5,1) - t523 * mrSges(5,3) - t549 * t533 + t558 * t534;
t515 = t589 * t517 + t585 * t518;
t522 = -t549 * qJD(4) - t585 * t543 + t589 * t563;
t535 = t558 * mrSges(5,1) - t549 * mrSges(5,3);
t512 = m(5) * t515 - t540 * mrSges(5,2) + t522 * mrSges(5,3) + t548 * t533 - t558 * t535;
t505 = -t585 * t511 + t589 * t512;
t546 = -t559 * mrSges(4,1) + t560 * mrSges(4,2);
t551 = t574 * mrSges(4,1) - t560 * mrSges(4,3);
t503 = m(4) * t520 - t563 * mrSges(4,2) + t542 * mrSges(4,3) + t559 * t546 - t574 * t551 + t505;
t519 = -t586 * t531 + t590 * t532;
t516 = -t563 * pkin(3) - t573 * pkin(8) + t560 * t547 - t519;
t513 = -m(5) * t516 + t522 * mrSges(5,1) - t523 * mrSges(5,2) + t548 * t534 - t549 * t535;
t550 = -t574 * mrSges(4,2) + t559 * mrSges(4,3);
t509 = m(4) * t519 + t563 * mrSges(4,1) - t543 * mrSges(4,3) - t560 * t546 + t574 * t550 + t513;
t599 = t590 * t503 - t586 * t509;
t494 = m(3) * t545 - t578 * mrSges(3,2) - t571 * mrSges(3,3) - t579 * t564 + t568 * t601 + t599;
t497 = t586 * t503 + t590 * t509;
t555 = -t583 * t566 - t612;
t565 = -t579 * mrSges(3,2) + mrSges(3,3) * t601;
t496 = m(3) * t555 + t571 * mrSges(3,1) + t570 * mrSges(3,2) + (t564 * t587 - t565 * t591) * t605 + t497;
t504 = t589 * t511 + t585 * t512;
t596 = -m(4) * t530 + t542 * mrSges(4,1) - t543 * mrSges(4,2) + t559 * t550 - t560 * t551 - t504;
t500 = m(3) * t544 + t578 * mrSges(3,1) - t570 * mrSges(3,3) + t579 * t565 - t568 * t602 + t596;
t483 = t494 * t609 - t583 * t496 + t500 * t608;
t480 = m(2) * t575 + qJDD(1) * mrSges(2,1) - t593 * mrSges(2,2) + t483;
t488 = t591 * t494 - t587 * t500;
t486 = m(2) * t576 - t593 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t488;
t607 = t592 * t480 + t588 * t486;
t482 = t494 * t611 + t584 * t496 + t500 * t610;
t600 = -t588 * t480 + t592 * t486;
t524 = Ifges(5,5) * t549 + Ifges(5,6) * t548 + Ifges(5,3) * t558;
t526 = Ifges(5,1) * t549 + Ifges(5,4) * t548 + Ifges(5,5) * t558;
t506 = -mrSges(5,1) * t516 + mrSges(5,3) * t515 + Ifges(5,4) * t523 + Ifges(5,2) * t522 + Ifges(5,6) * t540 - t549 * t524 + t558 * t526;
t525 = Ifges(5,4) * t549 + Ifges(5,2) * t548 + Ifges(5,6) * t558;
t507 = mrSges(5,2) * t516 - mrSges(5,3) * t514 + Ifges(5,1) * t523 + Ifges(5,4) * t522 + Ifges(5,5) * t540 + t548 * t524 - t558 * t525;
t536 = Ifges(4,5) * t560 + Ifges(4,6) * t559 + Ifges(4,3) * t574;
t537 = Ifges(4,4) * t560 + Ifges(4,2) * t559 + Ifges(4,6) * t574;
t489 = mrSges(4,2) * t530 - mrSges(4,3) * t519 + Ifges(4,1) * t543 + Ifges(4,4) * t542 + Ifges(4,5) * t563 - pkin(8) * t504 - t585 * t506 + t589 * t507 + t559 * t536 - t574 * t537;
t538 = Ifges(4,1) * t560 + Ifges(4,4) * t559 + Ifges(4,5) * t574;
t595 = mrSges(5,1) * t514 - mrSges(5,2) * t515 + Ifges(5,5) * t523 + Ifges(5,6) * t522 + Ifges(5,3) * t540 + t549 * t525 - t548 * t526;
t490 = -mrSges(4,1) * t530 + mrSges(4,3) * t520 + Ifges(4,4) * t543 + Ifges(4,2) * t542 + Ifges(4,6) * t563 - pkin(3) * t504 - t560 * t536 + t574 * t538 - t595;
t553 = Ifges(3,6) * t579 + (Ifges(3,4) * t587 + Ifges(3,2) * t591) * t605;
t554 = Ifges(3,5) * t579 + (Ifges(3,1) * t587 + Ifges(3,4) * t591) * t605;
t474 = Ifges(3,5) * t570 - Ifges(3,6) * t571 + Ifges(3,3) * t578 + mrSges(3,1) * t544 - mrSges(3,2) * t545 + t586 * t489 + t590 * t490 + pkin(2) * t596 + pkin(7) * t599 + (t553 * t587 - t554 * t591) * t605;
t552 = Ifges(3,3) * t579 + (Ifges(3,5) * t587 + Ifges(3,6) * t591) * t605;
t476 = mrSges(3,2) * t555 - mrSges(3,3) * t544 + Ifges(3,1) * t570 - Ifges(3,4) * t571 + Ifges(3,5) * t578 - pkin(7) * t497 + t590 * t489 - t586 * t490 + t552 * t601 - t579 * t553;
t594 = mrSges(4,1) * t519 - mrSges(4,2) * t520 + Ifges(4,5) * t543 + Ifges(4,6) * t542 + Ifges(4,3) * t563 + pkin(3) * t513 + pkin(8) * t505 + t589 * t506 + t585 * t507 + t560 * t537 - t559 * t538;
t478 = -mrSges(3,1) * t555 + mrSges(3,3) * t545 + Ifges(3,4) * t570 - Ifges(3,2) * t571 + Ifges(3,6) * t578 - pkin(2) * t497 - t552 * t602 + t579 * t554 - t594;
t597 = mrSges(2,1) * t575 - mrSges(2,2) * t576 + Ifges(2,3) * qJDD(1) + pkin(1) * t483 + t584 * t474 + t476 * t611 + t478 * t610 + t488 * t613;
t472 = -mrSges(2,2) * g(3) - mrSges(2,3) * t575 + Ifges(2,5) * qJDD(1) - t593 * Ifges(2,6) + t591 * t476 - t587 * t478 + (-t482 * t583 - t483 * t584) * pkin(6);
t471 = mrSges(2,1) * g(3) + mrSges(2,3) * t576 + t593 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t482 - t583 * t474 + (pkin(6) * t488 + t476 * t587 + t478 * t591) * t584;
t1 = [-m(1) * g(1) + t600; -m(1) * g(2) + t607; (-m(1) - m(2)) * g(3) + t482; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t607 - t588 * t471 + t592 * t472; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t600 + t592 * t471 + t588 * t472; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t597; t597; t474; t594; t595;];
tauJB = t1;
