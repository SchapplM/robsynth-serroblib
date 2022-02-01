% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:46
% EndTime: 2022-01-20 10:19:49
% DurationCPUTime: 1.92s
% Computational Cost: add. (23252->208), mult. (29959->247), div. (0->0), fcn. (14555->8), ass. (0->93)
t621 = Ifges(5,1) + Ifges(6,1);
t613 = Ifges(5,4) + Ifges(6,4);
t612 = Ifges(5,5) + Ifges(6,5);
t620 = Ifges(5,2) + Ifges(6,2);
t611 = Ifges(5,6) + Ifges(6,6);
t619 = Ifges(5,3) + Ifges(6,3);
t574 = qJD(1) + qJD(2);
t581 = sin(qJ(4));
t584 = cos(qJ(4));
t548 = (-mrSges(6,1) * t584 + mrSges(6,2) * t581) * t574;
t573 = qJDD(1) + qJDD(2);
t602 = qJD(4) * t574;
t598 = t584 * t602;
t550 = t573 * t581 + t598;
t583 = sin(qJ(1));
t586 = cos(qJ(1));
t564 = t583 * g(1) - g(2) * t586;
t556 = qJDD(1) * pkin(1) + t564;
t565 = -g(1) * t586 - g(2) * t583;
t587 = qJD(1) ^ 2;
t557 = -pkin(1) * t587 + t565;
t582 = sin(qJ(2));
t585 = cos(qJ(2));
t534 = t585 * t556 - t557 * t582;
t531 = pkin(2) * t573 + t534;
t535 = t582 * t556 + t585 * t557;
t572 = t574 ^ 2;
t532 = -pkin(2) * t572 + t535;
t579 = sin(pkin(8));
t580 = cos(pkin(8));
t527 = t579 * t531 + t580 * t532;
t524 = -pkin(3) * t572 + pkin(7) * t573 + t527;
t578 = -g(3) + qJDD(3);
t567 = t584 * t578;
t601 = qJD(5) * t574;
t615 = pkin(4) * t572;
t517 = qJDD(4) * pkin(4) + t567 + (-t550 + t598) * qJ(5) + (t584 * t615 - t524 - 0.2e1 * t601) * t581;
t608 = t574 * t584;
t561 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t608;
t600 = m(6) * t517 + qJDD(4) * mrSges(6,1) + qJD(4) * t561;
t609 = t574 * t581;
t514 = -t550 * mrSges(6,3) - t548 * t609 + t600;
t521 = t524 * t584 + t581 * t578;
t551 = t573 * t584 - t581 * t602;
t558 = qJD(4) * pkin(4) - qJ(5) * t609;
t577 = t584 ^ 2;
t518 = qJ(5) * t551 - qJD(4) * t558 - t577 * t615 + 0.2e1 * t584 * t601 + t521;
t520 = -t581 * t524 + t567;
t604 = (t581 * t621 + t584 * t613) * t574 + t612 * qJD(4);
t605 = (t581 * t613 + t584 * t620) * t574 + t611 * qJD(4);
t618 = mrSges(5,1) * t520 + mrSges(6,1) * t517 - mrSges(5,2) * t521 - mrSges(6,2) * t518 + pkin(4) * t514 + t619 * qJDD(4) + t612 * t550 + t611 * t551 + (t581 * t605 - t584 * t604) * t574;
t614 = -mrSges(5,2) - mrSges(6,2);
t549 = (-mrSges(5,1) * t584 + mrSges(5,2) * t581) * t574;
t562 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t608;
t511 = m(5) * t520 + qJDD(4) * mrSges(5,1) + qJD(4) * t562 + (-t548 - t549) * t609 + (-mrSges(5,3) - mrSges(6,3)) * t550 + t600;
t599 = m(6) * t518 + t551 * mrSges(6,3) + t548 * t608;
t559 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t609;
t603 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t609 - t559;
t512 = m(5) * t521 + t551 * mrSges(5,3) + qJD(4) * t603 + qJDD(4) * t614 + t549 * t608 + t599;
t594 = -t511 * t581 + t512 * t584;
t501 = m(4) * t527 - mrSges(4,1) * t572 - mrSges(4,2) * t573 + t594;
t526 = t580 * t531 - t579 * t532;
t592 = -t573 * pkin(3) - t526;
t523 = -pkin(7) * t572 + t592;
t519 = t558 * t609 - t551 * pkin(4) + qJDD(5) + (-qJ(5) * t577 - pkin(7)) * t572 + t592;
t593 = -m(6) * t519 + t551 * mrSges(6,1) + t561 * t608;
t589 = -m(5) * t523 + t551 * mrSges(5,1) + t550 * t614 + t562 * t608 + t603 * t609 + t593;
t506 = m(4) * t526 + t573 * mrSges(4,1) - t572 * mrSges(4,2) + t589;
t494 = t501 * t579 + t506 * t580;
t491 = m(3) * t534 + mrSges(3,1) * t573 - mrSges(3,2) * t572 + t494;
t595 = t501 * t580 - t506 * t579;
t492 = m(3) * t535 - mrSges(3,1) * t572 - mrSges(3,2) * t573 + t595;
t486 = t491 * t585 + t492 * t582;
t483 = m(2) * t564 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t587 + t486;
t596 = -t491 * t582 + t492 * t585;
t484 = m(2) * t565 - mrSges(2,1) * t587 - qJDD(1) * mrSges(2,2) + t596;
t607 = t483 * t586 + t484 * t583;
t504 = t511 * t584 + t512 * t581;
t606 = (t581 * t612 + t584 * t611) * t574 + t619 * qJD(4);
t502 = m(4) * t578 + t504;
t597 = -t483 * t583 + t484 * t586;
t513 = t550 * mrSges(6,2) + t559 * t609 - t593;
t496 = -mrSges(5,1) * t523 + mrSges(5,3) * t521 - mrSges(6,1) * t519 + mrSges(6,3) * t518 - pkin(4) * t513 + qJ(5) * t599 - t606 * t609 + t620 * t551 + t613 * t550 + (-mrSges(6,2) * qJ(5) + t611) * qJDD(4) + (-qJ(5) * t559 + t604) * qJD(4);
t498 = mrSges(5,2) * t523 + mrSges(6,2) * t519 - mrSges(5,3) * t520 - mrSges(6,3) * t517 - qJ(5) * t514 - t605 * qJD(4) + t612 * qJDD(4) + t550 * t621 + t613 * t551 + t606 * t608;
t590 = mrSges(3,1) * t534 + mrSges(4,1) * t526 - mrSges(3,2) * t535 - mrSges(4,2) * t527 + pkin(2) * t494 + pkin(3) * t589 + pkin(7) * t594 + t584 * t496 + t581 * t498 + (Ifges(4,3) + Ifges(3,3)) * t573;
t588 = mrSges(2,1) * t564 - mrSges(2,2) * t565 + Ifges(2,3) * qJDD(1) + pkin(1) * t486 + t590;
t487 = -mrSges(4,1) * t578 + mrSges(4,3) * t527 + t572 * Ifges(4,5) + Ifges(4,6) * t573 - pkin(3) * t504 - t618;
t479 = mrSges(4,2) * t578 - mrSges(4,3) * t526 + Ifges(4,5) * t573 - Ifges(4,6) * t572 - pkin(7) * t504 - t496 * t581 + t498 * t584;
t478 = -mrSges(3,2) * g(3) - mrSges(3,3) * t534 + Ifges(3,5) * t573 - Ifges(3,6) * t572 - qJ(3) * t494 + t479 * t580 - t487 * t579;
t477 = mrSges(3,1) * g(3) + mrSges(3,3) * t535 + t572 * Ifges(3,5) + Ifges(3,6) * t573 - pkin(2) * t502 + qJ(3) * t595 + t579 * t479 + t580 * t487;
t476 = -mrSges(2,2) * g(3) - mrSges(2,3) * t564 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t587 - pkin(6) * t486 - t477 * t582 + t478 * t585;
t475 = Ifges(2,6) * qJDD(1) + t587 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t565 + t582 * t478 + t585 * t477 - pkin(1) * (-m(3) * g(3) + t502) + pkin(6) * t596;
t1 = [-m(1) * g(1) + t597; -m(1) * g(2) + t607; (-m(1) - m(2) - m(3)) * g(3) + t502; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t607 - t475 * t583 + t476 * t586; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t597 + t586 * t475 + t583 * t476; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t588; t588; t590; t502; t618; t513;];
tauJB = t1;
