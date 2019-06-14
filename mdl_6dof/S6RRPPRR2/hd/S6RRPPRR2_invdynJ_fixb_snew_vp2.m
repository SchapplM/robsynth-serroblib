% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-05-06 09:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:45:13
% EndTime: 2019-05-06 09:45:23
% DurationCPUTime: 8.72s
% Computational Cost: add. (119538->349), mult. (284040->445), div. (0->0), fcn. (207254->12), ass. (0->136)
t558 = -2 * qJD(3);
t531 = sin(qJ(2));
t535 = cos(qJ(2));
t551 = qJD(1) * qJD(2);
t518 = qJDD(1) * t531 + t535 * t551;
t538 = qJD(1) ^ 2;
t532 = sin(qJ(1));
t536 = cos(qJ(1));
t545 = -g(1) * t536 - g(2) * t532;
t515 = -pkin(1) * t538 + qJDD(1) * pkin(7) + t545;
t555 = t531 * t515;
t557 = pkin(2) * t538;
t473 = qJDD(2) * pkin(2) - t518 * qJ(3) - t555 + (qJ(3) * t551 + t531 * t557 - g(3)) * t535;
t500 = -g(3) * t531 + t535 * t515;
t519 = qJDD(1) * t535 - t531 * t551;
t554 = qJD(1) * t531;
t520 = qJD(2) * pkin(2) - qJ(3) * t554;
t525 = t535 ^ 2;
t475 = qJ(3) * t519 - qJD(2) * t520 - t525 * t557 + t500;
t527 = sin(pkin(10));
t556 = cos(pkin(10));
t509 = (t527 * t535 + t531 * t556) * qJD(1);
t454 = t473 * t556 - t527 * t475 + t509 * t558;
t553 = qJD(1) * t535;
t508 = t527 * t554 - t556 * t553;
t455 = t527 * t473 + t556 * t475 + t508 * t558;
t488 = mrSges(4,1) * t508 + mrSges(4,2) * t509;
t492 = t518 * t527 - t556 * t519;
t502 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t509;
t487 = pkin(3) * t508 - qJ(4) * t509;
t537 = qJD(2) ^ 2;
t441 = -pkin(3) * t537 + qJDD(2) * qJ(4) - t487 * t508 + t455;
t550 = t532 * g(1) - t536 * g(2);
t542 = -qJDD(1) * pkin(1) - t550;
t477 = -t519 * pkin(2) + qJDD(3) + t520 * t554 + (-qJ(3) * t525 - pkin(7)) * t538 + t542;
t493 = t518 * t556 + t527 * t519;
t444 = (qJD(2) * t508 - t493) * qJ(4) + (qJD(2) * t509 + t492) * pkin(3) + t477;
t526 = sin(pkin(11));
t528 = cos(pkin(11));
t498 = qJD(2) * t526 + t509 * t528;
t430 = -0.2e1 * qJD(4) * t498 - t526 * t441 + t528 * t444;
t483 = qJDD(2) * t526 + t493 * t528;
t497 = qJD(2) * t528 - t509 * t526;
t420 = (t497 * t508 - t483) * pkin(8) + (t497 * t498 + t492) * pkin(4) + t430;
t431 = 0.2e1 * qJD(4) * t497 + t528 * t441 + t526 * t444;
t480 = pkin(4) * t508 - pkin(8) * t498;
t482 = qJDD(2) * t528 - t493 * t526;
t496 = t497 ^ 2;
t422 = -pkin(4) * t496 + pkin(8) * t482 - t480 * t508 + t431;
t530 = sin(qJ(5));
t534 = cos(qJ(5));
t414 = t534 * t420 - t530 * t422;
t469 = t497 * t534 - t498 * t530;
t447 = qJD(5) * t469 + t482 * t530 + t483 * t534;
t470 = t497 * t530 + t498 * t534;
t491 = qJDD(5) + t492;
t507 = qJD(5) + t508;
t412 = (t469 * t507 - t447) * pkin(9) + (t469 * t470 + t491) * pkin(5) + t414;
t415 = t530 * t420 + t534 * t422;
t446 = -qJD(5) * t470 + t482 * t534 - t483 * t530;
t461 = pkin(5) * t507 - pkin(9) * t470;
t468 = t469 ^ 2;
t413 = -pkin(5) * t468 + pkin(9) * t446 - t461 * t507 + t415;
t529 = sin(qJ(6));
t533 = cos(qJ(6));
t410 = t412 * t533 - t413 * t529;
t456 = t469 * t533 - t470 * t529;
t428 = qJD(6) * t456 + t446 * t529 + t447 * t533;
t457 = t469 * t529 + t470 * t533;
t438 = -mrSges(7,1) * t456 + mrSges(7,2) * t457;
t503 = qJD(6) + t507;
t448 = -mrSges(7,2) * t503 + mrSges(7,3) * t456;
t490 = qJDD(6) + t491;
t404 = m(7) * t410 + mrSges(7,1) * t490 - mrSges(7,3) * t428 - t438 * t457 + t448 * t503;
t411 = t412 * t529 + t413 * t533;
t427 = -qJD(6) * t457 + t446 * t533 - t447 * t529;
t449 = mrSges(7,1) * t503 - mrSges(7,3) * t457;
t405 = m(7) * t411 - mrSges(7,2) * t490 + mrSges(7,3) * t427 + t438 * t456 - t449 * t503;
t398 = t533 * t404 + t529 * t405;
t458 = -mrSges(6,1) * t469 + mrSges(6,2) * t470;
t459 = -mrSges(6,2) * t507 + mrSges(6,3) * t469;
t396 = m(6) * t414 + mrSges(6,1) * t491 - mrSges(6,3) * t447 - t458 * t470 + t459 * t507 + t398;
t460 = mrSges(6,1) * t507 - mrSges(6,3) * t470;
t546 = -t404 * t529 + t533 * t405;
t397 = m(6) * t415 - mrSges(6,2) * t491 + mrSges(6,3) * t446 + t458 * t469 - t460 * t507 + t546;
t392 = t534 * t396 + t530 * t397;
t474 = -mrSges(5,1) * t497 + mrSges(5,2) * t498;
t478 = -mrSges(5,2) * t508 + mrSges(5,3) * t497;
t390 = m(5) * t430 + mrSges(5,1) * t492 - mrSges(5,3) * t483 - t474 * t498 + t478 * t508 + t392;
t479 = mrSges(5,1) * t508 - mrSges(5,3) * t498;
t547 = -t396 * t530 + t534 * t397;
t391 = m(5) * t431 - mrSges(5,2) * t492 + mrSges(5,3) * t482 + t474 * t497 - t479 * t508 + t547;
t548 = -t390 * t526 + t528 * t391;
t382 = m(4) * t455 - qJDD(2) * mrSges(4,2) - mrSges(4,3) * t492 - qJD(2) * t502 - t488 * t508 + t548;
t440 = -qJDD(2) * pkin(3) - t537 * qJ(4) + t509 * t487 + qJDD(4) - t454;
t432 = -t482 * pkin(4) - t496 * pkin(8) + t498 * t480 + t440;
t417 = -t446 * pkin(5) - t468 * pkin(9) + t470 * t461 + t432;
t543 = m(7) * t417 - t427 * mrSges(7,1) + t428 * mrSges(7,2) - t456 * t448 + t457 * t449;
t540 = m(6) * t432 - t446 * mrSges(6,1) + mrSges(6,2) * t447 - t469 * t459 + t460 * t470 + t543;
t408 = m(5) * t440 - t482 * mrSges(5,1) + mrSges(5,2) * t483 - t497 * t478 + t479 * t498 + t540;
t501 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t508;
t407 = m(4) * t454 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t493 + qJD(2) * t501 - t488 * t509 - t408;
t379 = t527 * t382 + t556 * t407;
t384 = t528 * t390 + t526 * t391;
t549 = t556 * t382 - t407 * t527;
t434 = Ifges(7,4) * t457 + Ifges(7,2) * t456 + Ifges(7,6) * t503;
t435 = Ifges(7,1) * t457 + Ifges(7,4) * t456 + Ifges(7,5) * t503;
t541 = -mrSges(7,1) * t410 + mrSges(7,2) * t411 - Ifges(7,5) * t428 - Ifges(7,6) * t427 - Ifges(7,3) * t490 - t457 * t434 + t456 * t435;
t383 = m(4) * t477 + t492 * mrSges(4,1) + t493 * mrSges(4,2) + t508 * t501 + t509 * t502 + t384;
t451 = Ifges(6,4) * t470 + Ifges(6,2) * t469 + Ifges(6,6) * t507;
t452 = Ifges(6,1) * t470 + Ifges(6,4) * t469 + Ifges(6,5) * t507;
t539 = mrSges(6,1) * t414 - mrSges(6,2) * t415 + Ifges(6,5) * t447 + Ifges(6,6) * t446 + Ifges(6,3) * t491 + pkin(5) * t398 + t470 * t451 - t469 * t452 - t541;
t522 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t553;
t521 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t554;
t517 = (-mrSges(3,1) * t535 + mrSges(3,2) * t531) * qJD(1);
t514 = -t538 * pkin(7) + t542;
t512 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t531 + Ifges(3,4) * t535) * qJD(1);
t511 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t531 + Ifges(3,2) * t535) * qJD(1);
t499 = -t535 * g(3) - t555;
t486 = Ifges(4,1) * t509 - Ifges(4,4) * t508 + Ifges(4,5) * qJD(2);
t485 = Ifges(4,4) * t509 - Ifges(4,2) * t508 + Ifges(4,6) * qJD(2);
t484 = Ifges(4,5) * t509 - Ifges(4,6) * t508 + Ifges(4,3) * qJD(2);
t464 = Ifges(5,1) * t498 + Ifges(5,4) * t497 + Ifges(5,5) * t508;
t463 = Ifges(5,4) * t498 + Ifges(5,2) * t497 + Ifges(5,6) * t508;
t462 = Ifges(5,5) * t498 + Ifges(5,6) * t497 + Ifges(5,3) * t508;
t450 = Ifges(6,5) * t470 + Ifges(6,6) * t469 + Ifges(6,3) * t507;
t433 = Ifges(7,5) * t457 + Ifges(7,6) * t456 + Ifges(7,3) * t503;
t400 = mrSges(7,2) * t417 - mrSges(7,3) * t410 + Ifges(7,1) * t428 + Ifges(7,4) * t427 + Ifges(7,5) * t490 + t433 * t456 - t434 * t503;
t399 = -mrSges(7,1) * t417 + mrSges(7,3) * t411 + Ifges(7,4) * t428 + Ifges(7,2) * t427 + Ifges(7,6) * t490 - t433 * t457 + t435 * t503;
t386 = mrSges(6,2) * t432 - mrSges(6,3) * t414 + Ifges(6,1) * t447 + Ifges(6,4) * t446 + Ifges(6,5) * t491 - pkin(9) * t398 - t399 * t529 + t400 * t533 + t450 * t469 - t451 * t507;
t385 = -mrSges(6,1) * t432 + mrSges(6,3) * t415 + Ifges(6,4) * t447 + Ifges(6,2) * t446 + Ifges(6,6) * t491 - pkin(5) * t543 + pkin(9) * t546 + t533 * t399 + t529 * t400 - t470 * t450 + t507 * t452;
t378 = mrSges(5,2) * t440 - mrSges(5,3) * t430 + Ifges(5,1) * t483 + Ifges(5,4) * t482 + Ifges(5,5) * t492 - pkin(8) * t392 - t385 * t530 + t386 * t534 + t462 * t497 - t463 * t508;
t377 = -mrSges(5,1) * t440 + mrSges(5,3) * t431 + Ifges(5,4) * t483 + Ifges(5,2) * t482 + Ifges(5,6) * t492 - pkin(4) * t540 + pkin(8) * t547 + t534 * t385 + t530 * t386 - t498 * t462 + t508 * t464;
t376 = -t509 * t484 + Ifges(4,4) * t493 + t497 * t464 - t498 * t463 - Ifges(5,6) * t482 - Ifges(5,5) * t483 + qJD(2) * t486 - mrSges(4,1) * t477 + mrSges(4,3) * t455 - mrSges(5,1) * t430 + mrSges(5,2) * t431 + (-Ifges(5,3) - Ifges(4,2)) * t492 - pkin(4) * t392 - t539 - pkin(3) * t384 + Ifges(4,6) * qJDD(2);
t375 = mrSges(4,2) * t477 - mrSges(4,3) * t454 + Ifges(4,1) * t493 - Ifges(4,4) * t492 + Ifges(4,5) * qJDD(2) - qJ(4) * t384 - qJD(2) * t485 - t377 * t526 + t378 * t528 - t484 * t508;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t550 - mrSges(2,2) * t545 + t531 * (mrSges(3,2) * t514 - mrSges(3,3) * t499 + Ifges(3,1) * t518 + Ifges(3,4) * t519 + Ifges(3,5) * qJDD(2) - qJ(3) * t379 - qJD(2) * t511 + t556 * t375 - t527 * t376) + t535 * (-mrSges(3,1) * t514 + mrSges(3,3) * t500 + Ifges(3,4) * t518 + Ifges(3,2) * t519 + Ifges(3,6) * qJDD(2) - pkin(2) * t383 + qJ(3) * t549 + qJD(2) * t512 + t527 * t375 + t556 * t376) + pkin(1) * (-m(3) * t514 + t519 * mrSges(3,1) - t518 * mrSges(3,2) + (-t521 * t531 + t522 * t535) * qJD(1) - t383) + pkin(7) * (t535 * (m(3) * t500 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t519 - qJD(2) * t521 + t517 * t553 + t549) - t531 * (m(3) * t499 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t518 + qJD(2) * t522 - t517 * t554 + t379)); Ifges(3,5) * t518 + Ifges(3,6) * t519 + mrSges(3,1) * t499 - mrSges(3,2) * t500 + Ifges(4,5) * t493 - Ifges(4,6) * t492 + t509 * t485 + t508 * t486 + mrSges(4,1) * t454 - mrSges(4,2) * t455 + t526 * t378 + t528 * t377 - pkin(3) * t408 + qJ(4) * t548 + pkin(2) * t379 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t531 * t511 - t535 * t512) * qJD(1); t383; t408; t539; -t541;];
tauJ  = t1;
