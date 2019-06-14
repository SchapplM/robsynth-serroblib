% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPR11
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-05-06 16:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR11_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR11_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:00:20
% EndTime: 2019-05-06 16:00:30
% DurationCPUTime: 5.70s
% Computational Cost: add. (53814->334), mult. (114261->408), div. (0->0), fcn. (72821->10), ass. (0->133)
t579 = -2 * qJD(3);
t578 = Ifges(3,1) + Ifges(4,2);
t573 = Ifges(3,4) + Ifges(4,6);
t572 = Ifges(3,5) - Ifges(4,4);
t577 = Ifges(3,2) + Ifges(4,3);
t571 = Ifges(3,6) - Ifges(4,5);
t576 = (Ifges(3,3) + Ifges(4,1));
t546 = qJD(1) ^ 2;
t540 = sin(qJ(1));
t544 = cos(qJ(1));
t557 = -g(1) * t544 - g(2) * t540;
t502 = -pkin(1) * t546 + qJDD(1) * pkin(7) + t557;
t539 = sin(qJ(2));
t543 = cos(qJ(2));
t488 = -g(3) * t539 + t543 * t502;
t511 = (-pkin(2) * t543 - qJ(3) * t539) * qJD(1);
t545 = qJD(2) ^ 2;
t566 = qJD(1) * t543;
t463 = pkin(2) * t545 - qJDD(2) * qJ(3) + (qJD(2) * t579) - t511 * t566 - t488;
t575 = pkin(7) * t546;
t574 = mrSges(3,1) - mrSges(4,2);
t565 = qJD(1) * qJD(2);
t562 = t539 * t565;
t515 = qJDD(1) * t543 - t562;
t528 = t539 * qJD(1);
t522 = pkin(3) * t528 - (qJD(2) * pkin(8));
t534 = t543 ^ 2;
t563 = t543 * t565;
t514 = qJDD(1) * t539 + t563;
t561 = g(1) * t540 - t544 * g(2);
t555 = -qJDD(1) * pkin(1) - t561;
t550 = pkin(2) * t562 + t528 * t579 + (-t514 - t563) * qJ(3) + t555;
t441 = -t522 * t528 + (-pkin(3) * t534 - pkin(7)) * t546 + (-pkin(2) - pkin(8)) * t515 + t550;
t487 = -t543 * g(3) - t539 * t502;
t464 = -qJDD(2) * pkin(2) - qJ(3) * t545 + t511 * t528 + qJDD(3) - t487;
t451 = (-t539 * t543 * t546 - qJDD(2)) * pkin(8) + (t514 - t563) * pkin(3) + t464;
t538 = sin(qJ(4));
t542 = cos(qJ(4));
t430 = -t441 * t538 + t542 * t451;
t509 = -qJD(2) * t538 - t542 * t566;
t479 = qJD(4) * t509 + qJDD(2) * t542 - t515 * t538;
t508 = qJDD(4) + t514;
t510 = qJD(2) * t542 - t538 * t566;
t525 = t528 + qJD(4);
t420 = (t509 * t525 - t479) * qJ(5) + (t509 * t510 + t508) * pkin(4) + t430;
t431 = t542 * t441 + t538 * t451;
t478 = -qJD(4) * t510 - qJDD(2) * t538 - t515 * t542;
t485 = pkin(4) * t525 - qJ(5) * t510;
t507 = t509 ^ 2;
t422 = -pkin(4) * t507 + qJ(5) * t478 - t485 * t525 + t431;
t535 = sin(pkin(10));
t536 = cos(pkin(10));
t482 = t509 * t535 + t510 * t536;
t414 = -0.2e1 * qJD(5) * t482 + t536 * t420 - t422 * t535;
t456 = t478 * t535 + t479 * t536;
t481 = t509 * t536 - t510 * t535;
t411 = (t481 * t525 - t456) * pkin(9) + (t481 * t482 + t508) * pkin(5) + t414;
t415 = 0.2e1 * qJD(5) * t481 + t535 * t420 + t536 * t422;
t455 = t478 * t536 - t479 * t535;
t467 = pkin(5) * t525 - pkin(9) * t482;
t480 = t481 ^ 2;
t412 = -pkin(5) * t480 + pkin(9) * t455 - t467 * t525 + t415;
t537 = sin(qJ(6));
t541 = cos(qJ(6));
t409 = t411 * t541 - t412 * t537;
t458 = t481 * t541 - t482 * t537;
t427 = qJD(6) * t458 + t455 * t537 + t456 * t541;
t459 = t481 * t537 + t482 * t541;
t438 = -mrSges(7,1) * t458 + mrSges(7,2) * t459;
t523 = qJD(6) + t525;
t442 = -mrSges(7,2) * t523 + mrSges(7,3) * t458;
t505 = qJDD(6) + t508;
t405 = m(7) * t409 + mrSges(7,1) * t505 - mrSges(7,3) * t427 - t438 * t459 + t442 * t523;
t410 = t411 * t537 + t412 * t541;
t426 = -qJD(6) * t459 + t455 * t541 - t456 * t537;
t443 = mrSges(7,1) * t523 - mrSges(7,3) * t459;
t406 = m(7) * t410 - mrSges(7,2) * t505 + mrSges(7,3) * t426 + t438 * t458 - t443 * t523;
t399 = t541 * t405 + t537 * t406;
t460 = -mrSges(6,1) * t481 + mrSges(6,2) * t482;
t465 = -mrSges(6,2) * t525 + mrSges(6,3) * t481;
t396 = m(6) * t414 + mrSges(6,1) * t508 - mrSges(6,3) * t456 - t460 * t482 + t465 * t525 + t399;
t466 = mrSges(6,1) * t525 - mrSges(6,3) * t482;
t558 = -t405 * t537 + t541 * t406;
t397 = m(6) * t415 - mrSges(6,2) * t508 + mrSges(6,3) * t455 + t460 * t481 - t466 * t525 + t558;
t392 = t536 * t396 + t535 * t397;
t570 = (t576 * qJD(2)) + (t539 * t572 + t543 * t571) * qJD(1);
t569 = t571 * qJD(2) + (t539 * t573 + t543 * t577) * qJD(1);
t568 = t572 * qJD(2) + (t539 * t578 + t543 * t573) * qJD(1);
t520 = -mrSges(4,1) * t566 - (qJD(2) * mrSges(4,3));
t567 = -(qJD(2) * mrSges(3,2)) + mrSges(3,3) * t566 - t520;
t483 = -mrSges(5,1) * t509 + mrSges(5,2) * t510;
t484 = -mrSges(5,2) * t525 + mrSges(5,3) * t509;
t389 = m(5) * t430 + mrSges(5,1) * t508 - mrSges(5,3) * t479 - t483 * t510 + t484 * t525 + t392;
t486 = mrSges(5,1) * t525 - mrSges(5,3) * t510;
t559 = -t396 * t535 + t536 * t397;
t390 = m(5) * t431 - mrSges(5,2) * t508 + mrSges(5,3) * t478 + t483 * t509 - t486 * t525 + t559;
t560 = -t389 * t538 + t542 * t390;
t385 = t542 * t389 + t538 * t390;
t461 = -pkin(2) * t515 + t550 - t575;
t556 = m(4) * t461 + t560;
t450 = -pkin(8) * t534 * t546 + pkin(3) * t515 + qJD(2) * t522 - t463;
t433 = -pkin(4) * t478 - qJ(5) * t507 + t510 * t485 + qJDD(5) + t450;
t417 = -pkin(5) * t455 - pkin(9) * t480 + t467 * t482 + t433;
t553 = m(7) * t417 - t426 * mrSges(7,1) + t427 * mrSges(7,2) - t458 * t442 + t459 * t443;
t552 = m(4) * t464 + t514 * mrSges(4,1) + t385;
t435 = Ifges(7,4) * t459 + Ifges(7,2) * t458 + Ifges(7,6) * t523;
t436 = Ifges(7,1) * t459 + Ifges(7,4) * t458 + Ifges(7,5) * t523;
t551 = mrSges(7,1) * t409 - mrSges(7,2) * t410 + Ifges(7,5) * t427 + Ifges(7,6) * t426 + Ifges(7,3) * t505 + t459 * t435 - t458 * t436;
t407 = m(6) * t433 - t455 * mrSges(6,1) + t456 * mrSges(6,2) - t481 * t465 + t482 * t466 + t553;
t549 = -m(5) * t450 + t478 * mrSges(5,1) - t479 * mrSges(5,2) + t509 * t484 - t510 * t486 - t407;
t512 = (mrSges(4,2) * t543 - mrSges(4,3) * t539) * qJD(1);
t521 = mrSges(4,1) * t528 + qJD(2) * mrSges(4,2);
t548 = -m(4) * t463 + qJDD(2) * mrSges(4,3) + qJD(2) * t521 + t512 * t566 - t549;
t453 = Ifges(6,4) * t482 + Ifges(6,2) * t481 + Ifges(6,6) * t525;
t454 = Ifges(6,1) * t482 + Ifges(6,4) * t481 + Ifges(6,5) * t525;
t469 = Ifges(5,4) * t510 + Ifges(5,2) * t509 + Ifges(5,6) * t525;
t470 = Ifges(5,1) * t510 + Ifges(5,4) * t509 + Ifges(5,5) * t525;
t547 = mrSges(5,1) * t430 + mrSges(6,1) * t414 - mrSges(5,2) * t431 - mrSges(6,2) * t415 + Ifges(5,5) * t479 + Ifges(6,5) * t456 + Ifges(5,6) * t478 + Ifges(6,6) * t455 + pkin(4) * t392 + pkin(5) * t399 + t482 * t453 - t481 * t454 + t510 * t469 - t509 * t470 + t551 + (Ifges(6,3) + Ifges(5,3)) * t508;
t518 = (qJD(2) * mrSges(3,1)) - mrSges(3,3) * t528;
t513 = (-mrSges(3,1) * t543 + mrSges(3,2) * t539) * qJD(1);
t501 = t555 - t575;
t468 = Ifges(5,5) * t510 + Ifges(5,6) * t509 + Ifges(5,3) * t525;
t452 = Ifges(6,5) * t482 + Ifges(6,6) * t481 + Ifges(6,3) * t525;
t434 = Ifges(7,5) * t459 + Ifges(7,6) * t458 + Ifges(7,3) * t523;
t401 = mrSges(7,2) * t417 - mrSges(7,3) * t409 + Ifges(7,1) * t427 + Ifges(7,4) * t426 + Ifges(7,5) * t505 + t434 * t458 - t435 * t523;
t400 = -mrSges(7,1) * t417 + mrSges(7,3) * t410 + Ifges(7,4) * t427 + Ifges(7,2) * t426 + Ifges(7,6) * t505 - t434 * t459 + t436 * t523;
t387 = mrSges(6,2) * t433 - mrSges(6,3) * t414 + Ifges(6,1) * t456 + Ifges(6,4) * t455 + Ifges(6,5) * t508 - pkin(9) * t399 - t400 * t537 + t401 * t541 + t452 * t481 - t453 * t525;
t386 = -mrSges(6,1) * t433 + mrSges(6,3) * t415 + Ifges(6,4) * t456 + Ifges(6,2) * t455 + Ifges(6,6) * t508 - pkin(5) * t553 + pkin(9) * t558 + t541 * t400 + t537 * t401 - t482 * t452 + t525 * t454;
t384 = qJDD(2) * mrSges(4,2) + qJD(2) * t520 + t512 * t528 + t552;
t383 = mrSges(4,2) * t515 - mrSges(4,3) * t514 + (t520 * t543 - t521 * t539) * qJD(1) + t556;
t382 = mrSges(5,2) * t450 - mrSges(5,3) * t430 + Ifges(5,1) * t479 + Ifges(5,4) * t478 + Ifges(5,5) * t508 - qJ(5) * t392 - t386 * t535 + t387 * t536 + t468 * t509 - t469 * t525;
t381 = -mrSges(5,1) * t450 + mrSges(5,3) * t431 + Ifges(5,4) * t479 + Ifges(5,2) * t478 + Ifges(5,6) * t508 - pkin(4) * t407 + qJ(5) * t559 + t536 * t386 + t535 * t387 - t510 * t468 + t525 * t470;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t561 - mrSges(2,2) * t557 + t539 * (mrSges(4,1) * t464 + mrSges(3,2) * t501 - mrSges(3,3) * t487 - mrSges(4,3) * t461 + pkin(3) * t385 - qJ(3) * t383 - t569 * qJD(2) + t572 * qJDD(2) + t578 * t514 + t573 * t515 + t570 * t566 + t547) + t543 * (-mrSges(3,1) * t501 - mrSges(4,1) * t463 + mrSges(4,2) * t461 + mrSges(3,3) * t488 - pkin(2) * t383 - pkin(3) * t549 - pkin(8) * t560 + t568 * qJD(2) + t571 * qJDD(2) - t542 * t381 - t538 * t382 + t573 * t514 + t577 * t515 - t570 * t528) + pkin(1) * (-m(3) * t501 + t574 * t515 + (-mrSges(3,2) + mrSges(4,3)) * t514 + (t567 * t543 + (-t518 + t521) * t539) * qJD(1) - t556) + pkin(7) * (t543 * (-qJDD(2) * mrSges(3,2) + t513 * t566 + t548 + (mrSges(3,3) + mrSges(4,1)) * t515 + m(3) * t488 - qJD(2) * t518) + (-m(3) * t487 + t514 * mrSges(3,3) - t574 * qJDD(2) - t567 * qJD(2) + (t512 + t513) * t528 + t552) * t539); mrSges(3,1) * t487 - mrSges(3,2) * t488 + mrSges(4,2) * t464 - mrSges(4,3) * t463 + t542 * t382 - t538 * t381 - pkin(8) * t385 - pkin(2) * t384 + qJ(3) * t548 + (mrSges(4,1) * qJ(3) + t571) * t515 + t572 * t514 + t576 * qJDD(2) + (t539 * t569 - t543 * t568) * qJD(1); t384; t547; t407; t551;];
tauJ  = t1;
