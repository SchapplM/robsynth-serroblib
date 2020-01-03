% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRRP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP5_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP5_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP5_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:46
% EndTime: 2019-12-31 17:16:48
% DurationCPUTime: 1.61s
% Computational Cost: add. (11342->217), mult. (23510->266), div. (0->0), fcn. (13803->6), ass. (0->89)
t601 = Ifges(4,1) + Ifges(5,1);
t590 = Ifges(4,4) - Ifges(5,5);
t599 = Ifges(5,4) + Ifges(4,5);
t600 = Ifges(4,2) + Ifges(5,3);
t597 = Ifges(4,6) - Ifges(5,6);
t598 = -Ifges(5,2) - Ifges(4,3);
t565 = sin(qJ(3));
t568 = cos(qJ(2));
t584 = qJD(1) * t568;
t566 = sin(qJ(2));
t585 = qJD(1) * t566;
t593 = cos(qJ(3));
t537 = t565 * t585 - t593 * t584;
t538 = (t565 * t568 + t593 * t566) * qJD(1);
t562 = qJD(2) + qJD(3);
t596 = t600 * t537 - t590 * t538 - t597 * t562;
t595 = -t590 * t537 + t601 * t538 + t599 * t562;
t583 = qJD(1) * qJD(2);
t546 = t566 * qJDD(1) + t568 * t583;
t567 = sin(qJ(1));
t569 = cos(qJ(1));
t553 = -t569 * g(1) - t567 * g(2);
t570 = qJD(1) ^ 2;
t540 = -t570 * pkin(1) + qJDD(1) * pkin(5) + t553;
t589 = t566 * t540;
t592 = pkin(2) * t570;
t501 = qJDD(2) * pkin(2) - t546 * pkin(6) - t589 + (pkin(6) * t583 + t566 * t592 - g(3)) * t568;
t528 = -t566 * g(3) + t568 * t540;
t547 = t568 * qJDD(1) - t566 * t583;
t551 = qJD(2) * pkin(2) - pkin(6) * t585;
t564 = t568 ^ 2;
t502 = t547 * pkin(6) - qJD(2) * t551 - t564 * t592 + t528;
t496 = t565 * t501 + t593 * t502;
t510 = t538 * qJD(3) + t565 * t546 - t593 * t547;
t530 = t562 * mrSges(4,1) - t538 * mrSges(4,3);
t561 = qJDD(2) + qJDD(3);
t522 = t537 * pkin(3) - t538 * qJ(4);
t560 = t562 ^ 2;
t492 = -t560 * pkin(3) + t561 * qJ(4) + 0.2e1 * qJD(4) * t562 - t537 * t522 + t496;
t531 = -t562 * mrSges(5,1) + t538 * mrSges(5,2);
t582 = m(5) * t492 + t561 * mrSges(5,3) + t562 * t531;
t523 = t537 * mrSges(5,1) - t538 * mrSges(5,3);
t586 = -t537 * mrSges(4,1) - t538 * mrSges(4,2) - t523;
t591 = -mrSges(4,3) - mrSges(5,2);
t480 = m(4) * t496 - t561 * mrSges(4,2) + t591 * t510 - t562 * t530 + t586 * t537 + t582;
t495 = t593 * t501 - t565 * t502;
t511 = -t537 * qJD(3) + t593 * t546 + t565 * t547;
t529 = -t562 * mrSges(4,2) - t537 * mrSges(4,3);
t493 = -t561 * pkin(3) - t560 * qJ(4) + t538 * t522 + qJDD(4) - t495;
t532 = -t537 * mrSges(5,2) + t562 * mrSges(5,3);
t578 = -m(5) * t493 + t561 * mrSges(5,1) + t562 * t532;
t482 = m(4) * t495 + t561 * mrSges(4,1) + t591 * t511 + t562 * t529 + t586 * t538 + t578;
t474 = t565 * t480 + t593 * t482;
t527 = -t568 * g(3) - t589;
t535 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t566 + Ifges(3,2) * t568) * qJD(1);
t536 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t566 + Ifges(3,4) * t568) * qJD(1);
t487 = t511 * mrSges(5,2) + t538 * t523 - t578;
t573 = -mrSges(4,1) * t495 + mrSges(5,1) * t493 + mrSges(4,2) * t496 - mrSges(5,3) * t492 + pkin(3) * t487 - qJ(4) * t582 + t598 * t561 + t596 * t538 + (qJ(4) * t523 - t595) * t537 - t599 * t511 + (qJ(4) * mrSges(5,2) + t597) * t510;
t594 = mrSges(3,1) * t527 - mrSges(3,2) * t528 + Ifges(3,5) * t546 + Ifges(3,6) * t547 + Ifges(3,3) * qJDD(2) + pkin(2) * t474 + (t566 * t535 - t568 * t536) * qJD(1) - t573;
t545 = (-mrSges(3,1) * t568 + mrSges(3,2) * t566) * qJD(1);
t550 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t584;
t472 = m(3) * t527 + qJDD(2) * mrSges(3,1) - t546 * mrSges(3,3) + qJD(2) * t550 - t545 * t585 + t474;
t549 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t585;
t579 = t593 * t480 - t565 * t482;
t473 = m(3) * t528 - qJDD(2) * mrSges(3,2) + t547 * mrSges(3,3) - qJD(2) * t549 + t545 * t584 + t579;
t580 = -t566 * t472 + t568 * t473;
t464 = m(2) * t553 - t570 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t580;
t552 = t567 * g(1) - t569 * g(2);
t576 = -qJDD(1) * pkin(1) - t552;
t539 = -t570 * pkin(5) + t576;
t512 = -t547 * pkin(2) + t551 * t585 + (-pkin(6) * t564 - pkin(5)) * t570 + t576;
t489 = -0.2e1 * qJD(4) * t538 + (t537 * t562 - t511) * qJ(4) + (t538 * t562 + t510) * pkin(3) + t512;
t483 = m(5) * t489 + t510 * mrSges(5,1) - t511 * mrSges(5,3) - t538 * t531 + t537 * t532;
t574 = m(4) * t512 + t510 * mrSges(4,1) + t511 * mrSges(4,2) + t537 * t529 + t538 * t530 + t483;
t572 = -m(3) * t539 + t547 * mrSges(3,1) - t546 * mrSges(3,2) - t549 * t585 + t550 * t584 - t574;
t476 = m(2) * t552 + qJDD(1) * mrSges(2,1) - t570 * mrSges(2,2) + t572;
t588 = t567 * t464 + t569 * t476;
t466 = t568 * t472 + t566 * t473;
t587 = t597 * t537 - t599 * t538 + t598 * t562;
t581 = t569 * t464 - t567 * t476;
t467 = -mrSges(4,1) * t512 - mrSges(5,1) * t489 + mrSges(5,2) * t492 + mrSges(4,3) * t496 - pkin(3) * t483 - t600 * t510 + t590 * t511 + t587 * t538 + t597 * t561 + t595 * t562;
t468 = mrSges(4,2) * t512 + mrSges(5,2) * t493 - mrSges(4,3) * t495 - mrSges(5,3) * t489 - qJ(4) * t483 - t590 * t510 + t601 * t511 + t587 * t537 + t599 * t561 + t596 * t562;
t534 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t566 + Ifges(3,6) * t568) * qJD(1);
t459 = -mrSges(3,1) * t539 + mrSges(3,3) * t528 + Ifges(3,4) * t546 + Ifges(3,2) * t547 + Ifges(3,6) * qJDD(2) - pkin(2) * t574 + pkin(6) * t579 + qJD(2) * t536 + t593 * t467 + t565 * t468 - t534 * t585;
t461 = mrSges(3,2) * t539 - mrSges(3,3) * t527 + Ifges(3,1) * t546 + Ifges(3,4) * t547 + Ifges(3,5) * qJDD(2) - pkin(6) * t474 - qJD(2) * t535 - t565 * t467 + t593 * t468 + t534 * t584;
t575 = mrSges(2,1) * t552 - mrSges(2,2) * t553 + Ifges(2,3) * qJDD(1) + pkin(1) * t572 + pkin(5) * t580 + t568 * t459 + t566 * t461;
t457 = mrSges(2,1) * g(3) + mrSges(2,3) * t553 + t570 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t466 - t594;
t456 = -mrSges(2,2) * g(3) - mrSges(2,3) * t552 + Ifges(2,5) * qJDD(1) - t570 * Ifges(2,6) - pkin(5) * t466 - t566 * t459 + t568 * t461;
t1 = [-m(1) * g(1) + t581; -m(1) * g(2) + t588; (-m(1) - m(2)) * g(3) + t466; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t588 + t569 * t456 - t567 * t457; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t581 + t567 * t456 + t569 * t457; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t575; t575; t594; -t573; t487;];
tauJB = t1;
