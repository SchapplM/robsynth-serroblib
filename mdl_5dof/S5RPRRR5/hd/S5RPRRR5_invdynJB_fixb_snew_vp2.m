% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% m [6x1]
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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:48:53
% EndTime: 2022-01-20 09:48:57
% DurationCPUTime: 3.53s
% Computational Cost: add. (49986->229), mult. (67543->290), div. (0->0), fcn. (38166->10), ass. (0->100)
t587 = sin(qJ(1));
t591 = cos(qJ(1));
t564 = t587 * g(1) - t591 * g(2);
t558 = qJDD(1) * pkin(1) + t564;
t565 = -t591 * g(1) - t587 * g(2);
t592 = qJD(1) ^ 2;
t559 = -t592 * pkin(1) + t565;
t582 = sin(pkin(9));
t583 = cos(pkin(9));
t541 = t583 * t558 - t582 * t559;
t538 = qJDD(1) * pkin(2) + t541;
t542 = t582 * t558 + t583 * t559;
t539 = -t592 * pkin(2) + t542;
t586 = sin(qJ(3));
t590 = cos(qJ(3));
t523 = t586 * t538 + t590 * t539;
t577 = qJD(1) + qJD(3);
t573 = t577 ^ 2;
t575 = qJDD(1) + qJDD(3);
t520 = -t573 * pkin(3) + t575 * pkin(7) + t523;
t581 = -g(3) + qJDD(2);
t585 = sin(qJ(4));
t589 = cos(qJ(4));
t516 = -t585 * t520 + t589 * t581;
t608 = qJD(4) * t577;
t606 = t589 * t608;
t553 = t585 * t575 + t606;
t513 = (-t553 + t606) * pkin(8) + (t573 * t585 * t589 + qJDD(4)) * pkin(4) + t516;
t517 = t589 * t520 + t585 * t581;
t554 = t589 * t575 - t585 * t608;
t611 = t577 * t585;
t562 = qJD(4) * pkin(4) - pkin(8) * t611;
t580 = t589 ^ 2;
t514 = -t580 * t573 * pkin(4) + t554 * pkin(8) - qJD(4) * t562 + t517;
t584 = sin(qJ(5));
t588 = cos(qJ(5));
t511 = t588 * t513 - t584 * t514;
t548 = (-t584 * t585 + t588 * t589) * t577;
t529 = t548 * qJD(5) + t588 * t553 + t584 * t554;
t549 = (t584 * t589 + t585 * t588) * t577;
t534 = -t548 * mrSges(6,1) + t549 * mrSges(6,2);
t576 = qJD(4) + qJD(5);
t543 = -t576 * mrSges(6,2) + t548 * mrSges(6,3);
t574 = qJDD(4) + qJDD(5);
t508 = m(6) * t511 + t574 * mrSges(6,1) - t529 * mrSges(6,3) - t549 * t534 + t576 * t543;
t512 = t584 * t513 + t588 * t514;
t528 = -t549 * qJD(5) - t584 * t553 + t588 * t554;
t544 = t576 * mrSges(6,1) - t549 * mrSges(6,3);
t509 = m(6) * t512 - t574 * mrSges(6,2) + t528 * mrSges(6,3) + t548 * t534 - t576 * t544;
t499 = t588 * t508 + t584 * t509;
t546 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t585 + Ifges(5,2) * t589) * t577;
t547 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t585 + Ifges(5,4) * t589) * t577;
t531 = Ifges(6,4) * t549 + Ifges(6,2) * t548 + Ifges(6,6) * t576;
t532 = Ifges(6,1) * t549 + Ifges(6,4) * t548 + Ifges(6,5) * t576;
t596 = -mrSges(6,1) * t511 + mrSges(6,2) * t512 - Ifges(6,5) * t529 - Ifges(6,6) * t528 - Ifges(6,3) * t574 - t549 * t531 + t548 * t532;
t612 = mrSges(5,1) * t516 - mrSges(5,2) * t517 + Ifges(5,5) * t553 + Ifges(5,6) * t554 + Ifges(5,3) * qJDD(4) + pkin(4) * t499 + (t585 * t546 - t589 * t547) * t577 - t596;
t610 = t577 * t589;
t552 = (-mrSges(5,1) * t589 + mrSges(5,2) * t585) * t577;
t561 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t610;
t497 = m(5) * t516 + qJDD(4) * mrSges(5,1) - t553 * mrSges(5,3) + qJD(4) * t561 - t552 * t611 + t499;
t560 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t611;
t601 = -t584 * t508 + t588 * t509;
t498 = m(5) * t517 - qJDD(4) * mrSges(5,2) + t554 * mrSges(5,3) - qJD(4) * t560 + t552 * t610 + t601;
t602 = -t585 * t497 + t589 * t498;
t490 = m(4) * t523 - t573 * mrSges(4,1) - t575 * mrSges(4,2) + t602;
t522 = t590 * t538 - t586 * t539;
t599 = -t575 * pkin(3) - t522;
t519 = -t573 * pkin(7) + t599;
t515 = t562 * t611 - t554 * pkin(4) + (-pkin(8) * t580 - pkin(7)) * t573 + t599;
t597 = m(6) * t515 - t528 * mrSges(6,1) + t529 * mrSges(6,2) - t548 * t543 + t549 * t544;
t594 = -m(5) * t519 + t554 * mrSges(5,1) - t553 * mrSges(5,2) - t560 * t611 + t561 * t610 - t597;
t503 = m(4) * t522 + t575 * mrSges(4,1) - t573 * mrSges(4,2) + t594;
t485 = t586 * t490 + t590 * t503;
t482 = m(3) * t541 + qJDD(1) * mrSges(3,1) - t592 * mrSges(3,2) + t485;
t603 = t590 * t490 - t586 * t503;
t483 = m(3) * t542 - t592 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t603;
t475 = t583 * t482 + t582 * t483;
t472 = m(2) * t564 + qJDD(1) * mrSges(2,1) - t592 * mrSges(2,2) + t475;
t604 = -t582 * t482 + t583 * t483;
t473 = m(2) * t565 - t592 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t604;
t609 = t591 * t472 + t587 * t473;
t493 = t589 * t497 + t585 * t498;
t607 = m(4) * t581 + t493;
t605 = -t587 * t472 + t591 * t473;
t491 = m(3) * t581 + t607;
t530 = Ifges(6,5) * t549 + Ifges(6,6) * t548 + Ifges(6,3) * t576;
t500 = -mrSges(6,1) * t515 + mrSges(6,3) * t512 + Ifges(6,4) * t529 + Ifges(6,2) * t528 + Ifges(6,6) * t574 - t549 * t530 + t576 * t532;
t501 = mrSges(6,2) * t515 - mrSges(6,3) * t511 + Ifges(6,1) * t529 + Ifges(6,4) * t528 + Ifges(6,5) * t574 + t548 * t530 - t576 * t531;
t545 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t585 + Ifges(5,6) * t589) * t577;
t478 = -mrSges(5,1) * t519 + mrSges(5,3) * t517 + Ifges(5,4) * t553 + Ifges(5,2) * t554 + Ifges(5,6) * qJDD(4) - pkin(4) * t597 + pkin(8) * t601 + qJD(4) * t547 + t588 * t500 + t584 * t501 - t545 * t611;
t487 = mrSges(5,2) * t519 - mrSges(5,3) * t516 + Ifges(5,1) * t553 + Ifges(5,4) * t554 + Ifges(5,5) * qJDD(4) - pkin(8) * t499 - qJD(4) * t546 - t584 * t500 + t588 * t501 + t545 * t610;
t598 = mrSges(4,1) * t522 - mrSges(4,2) * t523 + Ifges(4,3) * t575 + pkin(3) * t594 + pkin(7) * t602 + t589 * t478 + t585 * t487;
t595 = mrSges(2,1) * t564 + mrSges(3,1) * t541 - mrSges(2,2) * t565 - mrSges(3,2) * t542 + pkin(1) * t475 + pkin(2) * t485 + t598 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t476 = -mrSges(4,1) * t581 + mrSges(4,3) * t523 + t573 * Ifges(4,5) + Ifges(4,6) * t575 - pkin(3) * t493 - t612;
t468 = mrSges(4,2) * t581 - mrSges(4,3) * t522 + Ifges(4,5) * t575 - t573 * Ifges(4,6) - pkin(7) * t493 - t585 * t478 + t589 * t487;
t467 = mrSges(3,2) * t581 - mrSges(3,3) * t541 + Ifges(3,5) * qJDD(1) - t592 * Ifges(3,6) - pkin(6) * t485 + t590 * t468 - t586 * t476;
t466 = -mrSges(3,1) * t581 + mrSges(3,3) * t542 + t592 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t607 + pkin(6) * t603 + t586 * t468 + t590 * t476;
t465 = -mrSges(2,2) * g(3) - mrSges(2,3) * t564 + Ifges(2,5) * qJDD(1) - t592 * Ifges(2,6) - qJ(2) * t475 - t582 * t466 + t583 * t467;
t464 = mrSges(2,1) * g(3) + mrSges(2,3) * t565 + t592 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t491 + qJ(2) * t604 + t583 * t466 + t582 * t467;
t1 = [-m(1) * g(1) + t605; -m(1) * g(2) + t609; (-m(1) - m(2)) * g(3) + t491; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t609 - t587 * t464 + t591 * t465; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t605 + t591 * t464 + t587 * t465; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t595; t595; t491; t598; t612; -t596;];
tauJB = t1;
