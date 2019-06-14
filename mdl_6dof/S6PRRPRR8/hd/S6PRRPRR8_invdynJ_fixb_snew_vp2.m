% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPRR8
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 06:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:16:11
% EndTime: 2019-05-05 06:16:14
% DurationCPUTime: 3.15s
% Computational Cost: add. (24557->284), mult. (52728->359), div. (0->0), fcn. (38990->14), ass. (0->132)
t523 = sin(qJ(3));
t527 = cos(qJ(3));
t562 = Ifges(4,4) + Ifges(5,6);
t573 = t523 * t562 + t527 * (Ifges(4,2) + Ifges(5,3));
t572 = t523 * (Ifges(4,1) + Ifges(5,2)) + t527 * t562;
t571 = -2 * qJD(4);
t561 = Ifges(4,5) - Ifges(5,4);
t560 = Ifges(4,6) - Ifges(5,5);
t515 = sin(pkin(12));
t518 = cos(pkin(12));
t501 = g(1) * t515 - g(2) * t518;
t514 = -g(3) + qJDD(1);
t517 = sin(pkin(6));
t520 = cos(pkin(6));
t570 = t501 * t520 + t514 * t517;
t479 = -t501 * t517 + t514 * t520;
t516 = sin(pkin(7));
t557 = t479 * t516;
t502 = -g(1) * t518 - g(2) * t515;
t524 = sin(qJ(2));
t528 = cos(qJ(2));
t456 = -t502 * t524 + t570 * t528;
t529 = qJD(2) ^ 2;
t452 = pkin(9) * t516 * t529 + qJDD(2) * pkin(2) + t456;
t519 = cos(pkin(7));
t558 = t452 * t519;
t567 = (t557 + t558) * t527;
t511 = qJD(2) * t519 + qJD(3);
t549 = qJD(2) * t516;
t544 = t523 * t549;
t566 = (pkin(3) * t511 + t571) * t544;
t457 = t528 * t502 + t570 * t524;
t546 = qJDD(2) * t516;
t453 = -pkin(2) * t529 + pkin(9) * t546 + t457;
t553 = t519 * t523;
t426 = t452 * t553 + t527 * t453 + t523 * t557;
t488 = (-pkin(3) * t527 - qJ(4) * t523) * t549;
t509 = t511 ^ 2;
t510 = qJDD(2) * t519 + qJDD(3);
t548 = qJD(2) * t527;
t543 = t516 * t548;
t421 = t509 * pkin(3) - t510 * qJ(4) - t488 * t543 + t511 * t571 - t426;
t565 = (t560 * t511 + t549 * t573) * t523 - t527 * (t561 * t511 + t549 * t572);
t564 = -pkin(3) - pkin(10);
t563 = mrSges(4,1) - mrSges(5,2);
t450 = t523 * t453;
t425 = -t450 + t567;
t485 = -mrSges(4,2) * t511 + mrSges(4,3) * t543;
t486 = -mrSges(5,1) * t543 - mrSges(5,3) * t511;
t489 = (mrSges(5,2) * t527 - mrSges(5,3) * t523) * t549;
t490 = (-mrSges(4,1) * t527 + mrSges(4,2) * t523) * t549;
t492 = (qJD(3) * t548 + qJDD(2) * t523) * t516;
t538 = -t509 * qJ(4) + t488 * t544 + qJDD(4) + t450;
t555 = t516 ^ 2 * t529;
t418 = t492 * pkin(4) + t564 * t510 + (-pkin(10) * t523 * t555 - t558 + (-pkin(4) * qJD(2) * t511 - t479) * t516) * t527 + t538;
t475 = t519 * t479;
t491 = pkin(4) * t544 - pkin(10) * t511;
t493 = -qJD(3) * t544 + t527 * t546;
t545 = t527 ^ 2 * t555;
t420 = -pkin(4) * t545 - t492 * qJ(4) + t475 + t564 * t493 + (-t452 + (-qJ(4) * t511 * t527 - t491 * t523) * qJD(2)) * t516 + t566;
t522 = sin(qJ(5));
t526 = cos(qJ(5));
t413 = t522 * t418 + t526 * t420;
t478 = t511 * t526 - t522 * t543;
t447 = -qJD(5) * t478 - t493 * t526 - t510 * t522;
t477 = -t511 * t522 - t526 * t543;
t454 = -mrSges(6,1) * t477 + mrSges(6,2) * t478;
t500 = qJD(5) + t544;
t461 = mrSges(6,1) * t500 - mrSges(6,3) * t478;
t483 = qJDD(5) + t492;
t455 = -pkin(5) * t477 - pkin(11) * t478;
t498 = t500 ^ 2;
t410 = -pkin(5) * t498 + pkin(11) * t483 + t455 * t477 + t413;
t416 = t493 * pkin(4) - pkin(10) * t545 + t511 * t491 - t421;
t448 = qJD(5) * t477 - t493 * t522 + t510 * t526;
t414 = (-t477 * t500 - t448) * pkin(11) + (t478 * t500 - t447) * pkin(5) + t416;
t521 = sin(qJ(6));
t525 = cos(qJ(6));
t407 = -t410 * t521 + t414 * t525;
t458 = -t478 * t521 + t500 * t525;
t429 = qJD(6) * t458 + t448 * t525 + t483 * t521;
t459 = t478 * t525 + t500 * t521;
t434 = -mrSges(7,1) * t458 + mrSges(7,2) * t459;
t476 = qJD(6) - t477;
t437 = -mrSges(7,2) * t476 + mrSges(7,3) * t458;
t445 = qJDD(6) - t447;
t404 = m(7) * t407 + mrSges(7,1) * t445 - mrSges(7,3) * t429 - t434 * t459 + t437 * t476;
t408 = t410 * t525 + t414 * t521;
t428 = -qJD(6) * t459 - t448 * t521 + t483 * t525;
t438 = mrSges(7,1) * t476 - mrSges(7,3) * t459;
t405 = m(7) * t408 - mrSges(7,2) * t445 + mrSges(7,3) * t428 + t434 * t458 - t438 * t476;
t540 = -t404 * t521 + t525 * t405;
t393 = m(6) * t413 - mrSges(6,2) * t483 + mrSges(6,3) * t447 + t454 * t477 - t461 * t500 + t540;
t412 = t418 * t526 - t420 * t522;
t460 = -mrSges(6,2) * t500 + mrSges(6,3) * t477;
t409 = -pkin(5) * t483 - pkin(11) * t498 + t455 * t478 - t412;
t534 = -m(7) * t409 + t428 * mrSges(7,1) - mrSges(7,2) * t429 + t458 * t437 - t438 * t459;
t400 = m(6) * t412 + mrSges(6,1) * t483 - mrSges(6,3) * t448 - t454 * t478 + t460 * t500 + t534;
t388 = t522 * t393 + t526 * t400;
t422 = -t510 * pkin(3) + t538 - t567;
t535 = -m(5) * t422 - t492 * mrSges(5,1) - t388;
t384 = m(4) * t425 - t492 * mrSges(4,3) + (t485 - t486) * t511 + t563 * t510 + (-t489 - t490) * t544 + t535;
t559 = t384 * t527;
t484 = mrSges(4,1) * t511 - mrSges(4,3) * t544;
t487 = mrSges(5,1) * t544 + mrSges(5,2) * t511;
t395 = t525 * t404 + t521 * t405;
t533 = -m(6) * t416 + t447 * mrSges(6,1) - t448 * mrSges(6,2) + t477 * t460 - t478 * t461 - t395;
t531 = -m(5) * t421 + t510 * mrSges(5,3) + t511 * t487 + t489 * t543 - t533;
t391 = t531 + t490 * t543 + (mrSges(4,3) + mrSges(5,1)) * t493 - t510 * mrSges(4,2) - t511 * t484 + m(4) * t426;
t552 = t391 * t553 + t519 * t559;
t542 = -t384 * t523 + t527 * t391;
t541 = t526 * t393 - t522 * t400;
t436 = -t516 * t452 + t475;
t424 = -t493 * pkin(3) + (-t511 * t543 - t492) * qJ(4) + t436 + t566;
t537 = m(5) * t424 - t492 * mrSges(5,3) + t486 * t543 + t541;
t430 = Ifges(7,5) * t459 + Ifges(7,6) * t458 + Ifges(7,3) * t476;
t432 = Ifges(7,1) * t459 + Ifges(7,4) * t458 + Ifges(7,5) * t476;
t398 = -mrSges(7,1) * t409 + mrSges(7,3) * t408 + Ifges(7,4) * t429 + Ifges(7,2) * t428 + Ifges(7,6) * t445 - t430 * t459 + t432 * t476;
t431 = Ifges(7,4) * t459 + Ifges(7,2) * t458 + Ifges(7,6) * t476;
t399 = mrSges(7,2) * t409 - mrSges(7,3) * t407 + Ifges(7,1) * t429 + Ifges(7,4) * t428 + Ifges(7,5) * t445 + t430 * t458 - t431 * t476;
t440 = Ifges(6,4) * t478 + Ifges(6,2) * t477 + Ifges(6,6) * t500;
t441 = Ifges(6,1) * t478 + Ifges(6,4) * t477 + Ifges(6,5) * t500;
t532 = mrSges(6,1) * t412 - mrSges(6,2) * t413 + Ifges(6,5) * t448 + Ifges(6,6) * t447 + Ifges(6,3) * t483 + pkin(5) * t534 + pkin(11) * t540 + t525 * t398 + t521 * t399 + t478 * t440 - t477 * t441;
t530 = mrSges(7,1) * t407 - mrSges(7,2) * t408 + Ifges(7,5) * t429 + Ifges(7,6) * t428 + Ifges(7,3) * t445 + t431 * t459 - t432 * t458;
t439 = Ifges(6,5) * t478 + Ifges(6,6) * t477 + Ifges(6,3) * t500;
t387 = t493 * mrSges(5,2) - t487 * t544 + t537;
t386 = t510 * mrSges(5,2) + t511 * t486 + t489 * t544 - t535;
t385 = m(4) * t436 + t492 * mrSges(4,2) - t563 * t493 + (-t485 * t527 + (t484 - t487) * t523) * t549 + t537;
t382 = -mrSges(6,1) * t416 + mrSges(6,3) * t413 + Ifges(6,4) * t448 + Ifges(6,2) * t447 + Ifges(6,6) * t483 - pkin(5) * t395 - t439 * t478 + t441 * t500 - t530;
t381 = mrSges(6,2) * t416 - mrSges(6,3) * t412 + Ifges(6,1) * t448 + Ifges(6,4) * t447 + Ifges(6,5) * t483 - pkin(11) * t395 - t398 * t521 + t399 * t525 + t439 * t477 - t440 * t500;
t380 = mrSges(4,1) * t425 - mrSges(4,2) * t426 + mrSges(5,2) * t422 - mrSges(5,3) * t421 + t526 * t381 - t522 * t382 - pkin(10) * t388 - pkin(3) * t386 + qJ(4) * t531 + (Ifges(4,3) + Ifges(5,1)) * t510 + (mrSges(5,1) * qJ(4) + t560) * t493 + t561 * t492 + t565 * t549;
t1 = [m(2) * t514 + t520 * (m(3) * t479 + t519 * t385 + (t391 * t523 + t559) * t516) + (t524 * (m(3) * t457 - mrSges(3,1) * t529 - qJDD(2) * mrSges(3,2) + t542) + t528 * (m(3) * t456 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t529 - t385 * t516 + t552)) * t517; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t456 - mrSges(3,2) * t457 + t519 * t380 + pkin(2) * t552 + (t523 * (mrSges(5,1) * t422 + mrSges(4,2) * t436 - mrSges(4,3) * t425 - mrSges(5,3) * t424 + pkin(4) * t388 - qJ(4) * t387 + t532) + t527 * (-mrSges(4,1) * t436 - mrSges(5,1) * t421 + mrSges(5,2) * t424 + mrSges(4,3) * t426 - pkin(3) * t387 - pkin(4) * t533 - pkin(10) * t541 - t522 * t381 - t526 * t382) - pkin(2) * t385 + pkin(9) * t542 - t565 * t511 + (t523 * t561 + t527 * t560) * t510 + t573 * t493 + t572 * t492) * t516; t380; t386; t532; t530;];
tauJ  = t1;
