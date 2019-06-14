% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:36:34
% EndTime: 2019-05-05 21:36:38
% DurationCPUTime: 2.32s
% Computational Cost: add. (16749->277), mult. (39224->329), div. (0->0), fcn. (28162->8), ass. (0->116)
t578 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t556 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t555 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t577 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t554 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t576 = Ifges(5,3) + Ifges(6,2) + Ifges(7,3);
t531 = qJD(1) ^ 2;
t528 = sin(qJ(1));
t529 = cos(qJ(1));
t541 = -g(1) * t529 - g(2) * t528;
t512 = -pkin(1) * t531 + qJDD(1) * qJ(2) + t541;
t524 = sin(pkin(9));
t525 = cos(pkin(9));
t557 = qJD(1) * qJD(2);
t546 = -t525 * g(3) - 0.2e1 * t524 * t557;
t565 = pkin(7) * qJDD(1);
t567 = pkin(2) * t531;
t484 = (t525 * t567 - t512 - t565) * t524 + t546;
t499 = -g(3) * t524 + (t512 + 0.2e1 * t557) * t525;
t523 = t525 ^ 2;
t485 = -t523 * t567 + t525 * t565 + t499;
t527 = sin(qJ(3));
t570 = cos(qJ(3));
t440 = t527 * t484 + t570 * t485;
t548 = t525 * t570;
t560 = qJD(1) * t524;
t510 = qJD(1) * t548 - t527 * t560;
t536 = t524 * t570 + t525 * t527;
t511 = t536 * qJD(1);
t494 = -pkin(3) * t510 - pkin(8) * t511;
t530 = qJD(3) ^ 2;
t436 = -pkin(3) * t530 + qJDD(3) * pkin(8) + t494 * t510 + t440;
t547 = t528 * g(1) - t529 * g(2);
t540 = qJDD(2) - t547;
t561 = -t524 ^ 2 - t523;
t495 = (-pkin(2) * t525 - pkin(1)) * qJDD(1) + (pkin(7) * t561 - qJ(2)) * t531 + t540;
t559 = qJD(3) * t511;
t496 = -t559 + (-t524 * t527 + t548) * qJDD(1);
t558 = t510 * qJD(3);
t497 = qJDD(1) * t536 + t558;
t438 = (-t497 - t558) * pkin(8) + (-t496 + t559) * pkin(3) + t495;
t526 = sin(qJ(4));
t569 = cos(qJ(4));
t431 = -t526 * t436 + t438 * t569;
t501 = -qJD(3) * t569 + t526 * t511;
t502 = t526 * qJD(3) + t511 * t569;
t468 = pkin(4) * t501 - qJ(5) * t502;
t493 = qJDD(4) - t496;
t508 = qJD(4) - t510;
t507 = t508 ^ 2;
t429 = -t493 * pkin(4) - t507 * qJ(5) + t502 * t468 + qJDD(5) - t431;
t475 = -mrSges(6,2) * t501 + mrSges(6,3) * t508;
t575 = -m(6) * t429 + t493 * mrSges(6,1) + t508 * t475;
t465 = -t501 * qJD(4) + t526 * qJDD(3) + t497 * t569;
t439 = t570 * t484 - t527 * t485;
t535 = qJDD(3) * pkin(3) + t530 * pkin(8) - t511 * t494 + t439;
t564 = t501 * t508;
t574 = (-t465 + t564) * qJ(5) - t535;
t476 = mrSges(7,2) * t508 + mrSges(7,3) * t501;
t572 = -0.2e1 * t502;
t422 = qJD(6) * t572 + (-t465 - t564) * qJ(6) + (t501 * t502 - t493) * pkin(5) + t429;
t470 = -mrSges(7,1) * t501 + mrSges(7,2) * t502;
t542 = -m(7) * t422 + t465 * mrSges(7,3) + t502 * t470;
t420 = -t493 * mrSges(7,1) - t508 * t476 - t542;
t469 = mrSges(6,1) * t501 - mrSges(6,3) * t502;
t417 = t465 * mrSges(6,2) + t502 * t469 + t420 - t575;
t432 = t569 * t436 + t526 * t438;
t571 = 2 * qJD(5);
t428 = -pkin(4) * t507 + t493 * qJ(5) - t501 * t468 + t508 * t571 + t432;
t464 = t502 * qJD(4) - qJDD(3) * t569 + t526 * t497;
t478 = -pkin(5) * t508 - qJ(6) * t502;
t500 = t501 ^ 2;
t424 = -pkin(5) * t500 + qJ(6) * t464 + 0.2e1 * qJD(6) * t501 + t478 * t508 + t428;
t479 = -mrSges(7,1) * t508 - mrSges(7,3) * t502;
t481 = -mrSges(6,1) * t508 + mrSges(6,2) * t502;
t552 = m(7) * t424 + t464 * mrSges(7,3) + t501 * t470;
t537 = m(6) * t428 + t493 * mrSges(6,3) + t508 * t481 + t552;
t549 = -t556 * t501 + t502 * t578 + t555 * t508;
t550 = t501 * t577 + t502 * t556 + t508 * t554;
t573 = -t464 * t554 + t465 * t555 + t576 * t493 + t501 * t549 + t502 * t550 + mrSges(5,1) * t431 - mrSges(6,1) * t429 - mrSges(7,1) * t422 - mrSges(5,2) * t432 + mrSges(7,2) * t424 + mrSges(6,3) * t428 - pkin(4) * t417 - pkin(5) * t420 + qJ(5) * (-t464 * mrSges(6,2) + t493 * mrSges(7,2) - t501 * t469 + t508 * t479 + t537);
t566 = -mrSges(5,3) - mrSges(6,2);
t492 = -mrSges(4,1) * t510 + mrSges(4,2) * t511;
t504 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t511;
t477 = -mrSges(5,2) * t508 - mrSges(5,3) * t501;
t562 = -mrSges(5,1) * t501 - mrSges(5,2) * t502 - t469;
t413 = m(5) * t431 + (t476 + t477) * t508 + t562 * t502 + (mrSges(5,1) + mrSges(7,1)) * t493 + t566 * t465 + t542 + t575;
t480 = mrSges(5,1) * t508 - mrSges(5,3) * t502;
t416 = m(5) * t432 + (t479 - t480) * t508 + t562 * t501 + (-mrSges(5,2) + mrSges(7,2)) * t493 + t566 * t464 + t537;
t544 = -t413 * t526 + t569 * t416;
t408 = m(4) * t440 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t496 - qJD(3) * t504 + t492 * t510 + t544;
t503 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t510;
t426 = -t500 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t464 + (-pkin(4) * t508 + t478 + t571) * t502 - t574;
t421 = m(7) * t426 - t464 * mrSges(7,1) + t465 * mrSges(7,2) - t501 * t476 + t502 * t479;
t430 = qJD(5) * t572 + (t502 * t508 + t464) * pkin(4) + t574;
t419 = m(6) * t430 + mrSges(6,1) * t464 - t465 * mrSges(6,3) + t475 * t501 - t502 * t481 - t421;
t532 = m(5) * t535 - t464 * mrSges(5,1) - mrSges(5,2) * t465 - t501 * t477 - t480 * t502 - t419;
t411 = m(4) * t439 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t497 + qJD(3) * t503 - t492 * t511 + t532;
t563 = t527 * t408 + t570 * t411;
t409 = t569 * t413 + t526 * t416;
t551 = t501 * t554 - t502 * t555 - t508 * t576;
t545 = t570 * t408 - t527 * t411;
t539 = -mrSges(3,1) * t525 + mrSges(3,2) * t524;
t538 = mrSges(3,3) * qJDD(1) + t531 * t539;
t534 = m(4) * t495 - t496 * mrSges(4,1) + t497 * mrSges(4,2) - t510 * t503 + t511 * t504 + t409;
t514 = (Ifges(3,5) * t524 + Ifges(3,6) * t525) * qJD(1);
t509 = -qJDD(1) * pkin(1) - t531 * qJ(2) + t540;
t498 = -t524 * t512 + t546;
t488 = Ifges(4,1) * t511 + Ifges(4,4) * t510 + Ifges(4,5) * qJD(3);
t487 = Ifges(4,4) * t511 + Ifges(4,2) * t510 + Ifges(4,6) * qJD(3);
t486 = Ifges(4,5) * t511 + Ifges(4,6) * t510 + Ifges(4,3) * qJD(3);
t405 = mrSges(3,3) * t531 * t561 + m(3) * t509 + qJDD(1) * t539 + t534;
t404 = -mrSges(5,2) * t535 + mrSges(6,2) * t429 + mrSges(7,2) * t426 - mrSges(5,3) * t431 - mrSges(6,3) * t430 - mrSges(7,3) * t422 - qJ(5) * t419 - qJ(6) * t420 - t556 * t464 + t465 * t578 + t555 * t493 + t551 * t501 - t550 * t508;
t403 = mrSges(5,1) * t535 + mrSges(5,3) * t432 - mrSges(6,1) * t430 + mrSges(6,2) * t428 + mrSges(7,1) * t426 - mrSges(7,3) * t424 + pkin(5) * t421 - qJ(6) * t552 - pkin(4) * t419 + (-qJ(6) * t479 + t549) * t508 + t551 * t502 + (-mrSges(7,2) * qJ(6) + t554) * t493 + t556 * t465 + t577 * t464;
t402 = -mrSges(4,1) * t495 + mrSges(4,3) * t440 + Ifges(4,4) * t497 + Ifges(4,2) * t496 + Ifges(4,6) * qJDD(3) - pkin(3) * t409 + qJD(3) * t488 - t511 * t486 - t573;
t401 = mrSges(4,2) * t495 - mrSges(4,3) * t439 + Ifges(4,1) * t497 + Ifges(4,4) * t496 + Ifges(4,5) * qJDD(3) - pkin(8) * t409 - qJD(3) * t487 - t526 * t403 + t404 * t569 + t510 * t486;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t547 - mrSges(2,2) * t541 + t524 * (t525 * qJD(1) * t514 + mrSges(3,2) * t509 - mrSges(3,3) * t498 + t570 * t401 - t527 * t402 - pkin(7) * t563 + (Ifges(3,1) * t524 + Ifges(3,4) * t525) * qJDD(1)) + t525 * (-t514 * t560 - mrSges(3,1) * t509 + mrSges(3,3) * t499 + t527 * t401 + t570 * t402 - pkin(2) * t534 + pkin(7) * t545 + (Ifges(3,4) * t524 + Ifges(3,2) * t525) * qJDD(1)) - pkin(1) * t405 + qJ(2) * ((m(3) * t499 + t525 * t538 + t545) * t525 + (-m(3) * t498 + t524 * t538 - t563) * t524); t405; mrSges(4,1) * t439 - mrSges(4,2) * t440 + Ifges(4,5) * t497 + Ifges(4,6) * t496 + Ifges(4,3) * qJDD(3) + pkin(3) * t532 + pkin(8) * t544 + t403 * t569 + t526 * t404 + t511 * t487 - t510 * t488; t573; t417; t421;];
tauJ  = t1;
