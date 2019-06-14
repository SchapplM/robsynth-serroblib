% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPR7
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
% Datum: 2019-05-06 14:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:43:38
% EndTime: 2019-05-06 14:43:45
% DurationCPUTime: 4.42s
% Computational Cost: add. (38728->334), mult. (85988->410), div. (0->0), fcn. (56964->10), ass. (0->133)
t586 = Ifges(3,1) + Ifges(4,1);
t579 = Ifges(3,4) - Ifges(4,5);
t578 = Ifges(3,5) + Ifges(4,4);
t585 = Ifges(3,2) + Ifges(4,3);
t577 = Ifges(3,6) - Ifges(4,6);
t584 = (Ifges(3,3) + Ifges(4,2));
t554 = qJD(1) ^ 2;
t548 = sin(qJ(1));
t552 = cos(qJ(1));
t563 = -g(1) * t552 - g(2) * t548;
t511 = -pkin(1) * t554 + qJDD(1) * pkin(7) + t563;
t547 = sin(qJ(2));
t551 = cos(qJ(2));
t492 = -g(3) * t547 + t551 * t511;
t515 = (-t551 * pkin(2) - t547 * qJ(3)) * qJD(1);
t553 = qJD(2) ^ 2;
t569 = qJD(1) * t551;
t582 = 2 * qJD(3);
t472 = -pkin(2) * t553 + qJDD(2) * qJ(3) + (qJD(2) * t582) + t515 * t569 + t492;
t568 = qJD(1) * qJD(2);
t566 = t547 * t568;
t519 = qJDD(1) * t551 - t566;
t570 = qJD(1) * t547;
t526 = -(qJD(2) * pkin(3)) - pkin(8) * t570;
t575 = t551 ^ 2 * t554;
t468 = -pkin(3) * t575 - pkin(8) * t519 + qJD(2) * t526 + t472;
t491 = -t551 * g(3) - t547 * t511;
t476 = -qJDD(2) * pkin(2) - qJ(3) * t553 + t515 * t570 + qJDD(3) - t491;
t567 = t551 * t568;
t518 = qJDD(1) * t547 + t567;
t469 = (-t518 + t567) * pkin(8) + (-t547 * t551 * t554 - qJDD(2)) * pkin(3) + t476;
t546 = sin(qJ(4));
t550 = cos(qJ(4));
t441 = -t468 * t546 + t550 * t469;
t508 = (-t547 * t546 - t551 * t550) * qJD(1);
t478 = qJD(4) * t508 + t518 * t550 - t519 * t546;
t509 = (-t551 * t546 + t547 * t550) * qJD(1);
t535 = -qJDD(2) + qJDD(4);
t536 = -qJD(2) + qJD(4);
t435 = (t508 * t536 - t478) * qJ(5) + (t508 * t509 + t535) * pkin(4) + t441;
t442 = t550 * t468 + t546 * t469;
t477 = -qJD(4) * t509 - t518 * t546 - t519 * t550;
t494 = pkin(4) * t536 - qJ(5) * t509;
t501 = t508 ^ 2;
t437 = -pkin(4) * t501 + qJ(5) * t477 - t494 * t536 + t442;
t542 = sin(pkin(10));
t543 = cos(pkin(10));
t488 = t508 * t543 - t509 * t542;
t581 = 2 * qJD(5);
t432 = t542 * t435 + t543 * t437 + t488 * t581;
t489 = t508 * t542 + t509 * t543;
t467 = -pkin(5) * t488 - pkin(9) * t489;
t534 = t536 ^ 2;
t430 = -pkin(5) * t534 + pkin(9) * t535 + t467 * t488 + t432;
t571 = t548 * g(1) - t552 * g(2);
t510 = -qJDD(1) * pkin(1) - t554 * pkin(7) - t571;
t561 = -t519 * pkin(2) + t510 + (-t518 - t567) * qJ(3);
t455 = -pkin(2) * t566 + t519 * pkin(3) - pkin(8) * t575 - t561 + (t526 + t582) * t570;
t439 = -t477 * pkin(4) - t501 * qJ(5) + t509 * t494 + qJDD(5) + t455;
t452 = t477 * t543 - t478 * t542;
t453 = t477 * t542 + t478 * t543;
t433 = (-t488 * t536 - t453) * pkin(9) + (t489 * t536 - t452) * pkin(5) + t439;
t545 = sin(qJ(6));
t549 = cos(qJ(6));
t427 = -t430 * t545 + t433 * t549;
t474 = -t489 * t545 + t536 * t549;
t444 = qJD(6) * t474 + t453 * t549 + t535 * t545;
t451 = qJDD(6) - t452;
t475 = t489 * t549 + t536 * t545;
t454 = -mrSges(7,1) * t474 + mrSges(7,2) * t475;
t482 = qJD(6) - t488;
t456 = -mrSges(7,2) * t482 + mrSges(7,3) * t474;
t424 = m(7) * t427 + mrSges(7,1) * t451 - mrSges(7,3) * t444 - t454 * t475 + t456 * t482;
t428 = t430 * t549 + t433 * t545;
t443 = -qJD(6) * t475 - t453 * t545 + t535 * t549;
t457 = mrSges(7,1) * t482 - mrSges(7,3) * t475;
t425 = m(7) * t428 - mrSges(7,2) * t451 + mrSges(7,3) * t443 + t454 * t474 - t457 * t482;
t416 = -t424 * t545 + t549 * t425;
t466 = -mrSges(6,1) * t488 + mrSges(6,2) * t489;
t480 = mrSges(6,1) * t536 - mrSges(6,3) * t489;
t413 = m(6) * t432 - mrSges(6,2) * t535 + mrSges(6,3) * t452 + t466 * t488 - t480 * t536 + t416;
t562 = -t435 * t543 + t437 * t542;
t429 = -pkin(5) * t535 - pkin(9) * t534 + (t581 + t467) * t489 + t562;
t426 = -m(7) * t429 + t443 * mrSges(7,1) - mrSges(7,2) * t444 + t474 * t456 - t457 * t475;
t431 = -0.2e1 * qJD(5) * t489 - t562;
t479 = -mrSges(6,2) * t536 + mrSges(6,3) * t488;
t420 = m(6) * t431 + mrSges(6,1) * t535 - mrSges(6,3) * t453 - t466 * t489 + t479 * t536 + t426;
t409 = t542 * t413 + t543 * t420;
t445 = Ifges(7,5) * t475 + Ifges(7,6) * t474 + Ifges(7,3) * t482;
t447 = Ifges(7,1) * t475 + Ifges(7,4) * t474 + Ifges(7,5) * t482;
t417 = -mrSges(7,1) * t429 + mrSges(7,3) * t428 + Ifges(7,4) * t444 + Ifges(7,2) * t443 + Ifges(7,6) * t451 - t445 * t475 + t447 * t482;
t446 = Ifges(7,4) * t475 + Ifges(7,2) * t474 + Ifges(7,6) * t482;
t418 = mrSges(7,2) * t429 - mrSges(7,3) * t427 + Ifges(7,1) * t444 + Ifges(7,4) * t443 + Ifges(7,5) * t451 + t445 * t474 - t446 * t482;
t459 = Ifges(6,4) * t489 + Ifges(6,2) * t488 + Ifges(6,6) * t536;
t460 = Ifges(6,1) * t489 + Ifges(6,4) * t488 + Ifges(6,5) * t536;
t484 = Ifges(5,4) * t509 + Ifges(5,2) * t508 + Ifges(5,6) * t536;
t485 = Ifges(5,1) * t509 + Ifges(5,4) * t508 + Ifges(5,5) * t536;
t583 = Ifges(5,5) * t478 + Ifges(5,6) * t477 + t509 * t484 - t508 * t485 + mrSges(5,1) * t441 - mrSges(5,2) * t442 + Ifges(6,5) * t453 + Ifges(6,6) * t452 + t489 * t459 - t488 * t460 + mrSges(6,1) * t431 - mrSges(6,2) * t432 + t545 * t418 + t549 * t417 + pkin(5) * t426 + pkin(9) * t416 + pkin(4) * t409 + (Ifges(5,3) + Ifges(6,3)) * t535;
t580 = mrSges(3,3) + mrSges(4,2);
t415 = t549 * t424 + t545 * t425;
t574 = (t584 * qJD(2)) + (t547 * t578 + t551 * t577) * qJD(1);
t573 = -t577 * qJD(2) + (-t547 * t579 - t585 * t551) * qJD(1);
t572 = t578 * qJD(2) + (t547 * t586 + t551 * t579) * qJD(1);
t490 = -mrSges(5,1) * t508 + mrSges(5,2) * t509;
t493 = -mrSges(5,2) * t536 + mrSges(5,3) * t508;
t407 = m(5) * t441 + mrSges(5,1) * t535 - mrSges(5,3) * t478 - t490 * t509 + t493 * t536 + t409;
t495 = mrSges(5,1) * t536 - mrSges(5,3) * t509;
t564 = t543 * t413 - t420 * t542;
t408 = m(5) * t442 - mrSges(5,2) * t535 + mrSges(5,3) * t477 + t490 * t508 - t495 * t536 + t564;
t565 = -t546 * t407 + t550 * t408;
t403 = t550 * t407 + t546 * t408;
t414 = m(6) * t439 - t452 * mrSges(6,1) + t453 * mrSges(6,2) - t488 * t479 + t489 * t480 + t415;
t516 = (-t551 * mrSges(4,1) - t547 * mrSges(4,3)) * qJD(1);
t523 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t570;
t560 = m(4) * t472 + qJDD(2) * mrSges(4,3) + qJD(2) * t523 + t516 * t569 + t565;
t525 = mrSges(4,2) * t569 + qJD(2) * mrSges(4,3);
t559 = m(4) * t476 - qJDD(2) * mrSges(4,1) - qJD(2) * t525 + t403;
t558 = m(5) * t455 - t477 * mrSges(5,1) + t478 * mrSges(5,2) - t508 * t493 + t509 * t495 + t414;
t557 = mrSges(7,1) * t427 - mrSges(7,2) * t428 + Ifges(7,5) * t444 + Ifges(7,6) * t443 + Ifges(7,3) * t451 + t446 * t475 - t474 * t447;
t470 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t570 + t561;
t556 = m(4) * t470 - t558;
t524 = -(qJD(2) * mrSges(3,2)) + mrSges(3,3) * t569;
t522 = (qJD(2) * mrSges(3,1)) - mrSges(3,3) * t570;
t517 = (-t551 * mrSges(3,1) + t547 * mrSges(3,2)) * qJD(1);
t483 = Ifges(5,5) * t509 + Ifges(5,6) * t508 + Ifges(5,3) * t536;
t458 = Ifges(6,5) * t489 + Ifges(6,6) * t488 + Ifges(6,3) * t536;
t410 = (-t523 * t547 - t525 * t551) * qJD(1) + t556 - t519 * mrSges(4,1) - t518 * mrSges(4,3);
t405 = -mrSges(6,1) * t439 + mrSges(6,3) * t432 + Ifges(6,4) * t453 + Ifges(6,2) * t452 + Ifges(6,6) * t535 - pkin(5) * t415 - t458 * t489 + t460 * t536 - t557;
t404 = mrSges(6,2) * t439 - mrSges(6,3) * t431 + Ifges(6,1) * t453 + Ifges(6,4) * t452 + Ifges(6,5) * t535 - pkin(9) * t415 - t417 * t545 + t418 * t549 + t458 * t488 - t459 * t536;
t402 = t518 * mrSges(4,2) + t516 * t570 + t559;
t401 = mrSges(5,2) * t455 - mrSges(5,3) * t441 + Ifges(5,1) * t478 + Ifges(5,4) * t477 + Ifges(5,5) * t535 - qJ(5) * t409 + t404 * t543 - t405 * t542 + t483 * t508 - t484 * t536;
t400 = -mrSges(5,1) * t455 + mrSges(5,3) * t442 + Ifges(5,4) * t478 + Ifges(5,2) * t477 + Ifges(5,6) * t535 - pkin(4) * t414 + qJ(5) * t564 + t542 * t404 + t543 * t405 - t509 * t483 + t536 * t485;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t571 - mrSges(2,2) * t563 + t547 * (mrSges(3,2) * t510 + mrSges(4,2) * t476 - mrSges(3,3) * t491 - mrSges(4,3) * t470 - pkin(8) * t403 - qJ(3) * t410 + t573 * qJD(2) + t578 * qJDD(2) - t546 * t400 + t550 * t401 + t586 * t518 + t579 * t519 + t574 * t569) + t551 * (-mrSges(3,1) * t510 - mrSges(4,1) * t470 + mrSges(4,2) * t472 + mrSges(3,3) * t492 - pkin(2) * t410 + pkin(3) * t558 - pkin(8) * t565 + t572 * qJD(2) + t577 * qJDD(2) - t550 * t400 - t546 * t401 + t579 * t518 + t585 * t519 - t574 * t570) + pkin(1) * (((t524 + t525) * t551 + (-t522 + t523) * t547) * qJD(1) + (mrSges(3,1) + mrSges(4,1)) * t519 + (-mrSges(3,2) + mrSges(4,3)) * t518 - t556 - m(3) * t510) + pkin(7) * (t551 * (m(3) * t492 - qJDD(2) * mrSges(3,2) - qJD(2) * t522 + t517 * t569 + t580 * t519 + t560) + (-m(3) * t491 - qJDD(2) * mrSges(3,1) - qJD(2) * t524 + t580 * t518 + (t516 + t517) * t570 + t559) * t547); (-t573 * t547 - t572 * t551) * qJD(1) + qJ(3) * t560 + (qJ(3) * mrSges(4,2) + t577) * t519 + t578 * t518 + t584 * qJDD(2) + mrSges(3,1) * t491 - mrSges(3,2) * t492 - mrSges(4,1) * t476 + mrSges(4,3) * t472 - pkin(3) * t403 - pkin(2) * t402 - t583; t402; t583; t414; t557;];
tauJ  = t1;
