% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-05-05 21:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:30:45
% EndTime: 2019-05-05 21:30:51
% DurationCPUTime: 4.59s
% Computational Cost: add. (44501->295), mult. (106642->363), div. (0->0), fcn. (80064->10), ass. (0->126)
t549 = Ifges(6,1) + Ifges(7,1);
t542 = Ifges(6,4) - Ifges(7,5);
t541 = Ifges(6,5) + Ifges(7,4);
t548 = -Ifges(6,2) - Ifges(7,3);
t547 = -Ifges(7,2) - Ifges(6,3);
t540 = Ifges(6,6) - Ifges(7,6);
t509 = qJD(1) ^ 2;
t500 = sin(pkin(9));
t501 = cos(pkin(9));
t503 = sin(qJ(3));
t506 = cos(qJ(3));
t515 = t500 * t503 - t501 * t506;
t489 = t515 * qJD(1);
t516 = t500 * t506 + t501 * t503;
t490 = t516 * qJD(1);
t529 = t490 * qJD(3);
t475 = -t515 * qJDD(1) - t529;
t504 = sin(qJ(1));
t507 = cos(qJ(1));
t519 = -g(1) * t507 - g(2) * t504;
t491 = -pkin(1) * t509 + qJDD(1) * qJ(2) + t519;
t528 = qJD(1) * qJD(2);
t524 = -t501 * g(3) - 0.2e1 * t500 * t528;
t538 = pkin(7) * qJDD(1);
t544 = pkin(2) * t509;
t464 = (t501 * t544 - t491 - t538) * t500 + t524;
t478 = -g(3) * t500 + (t491 + 0.2e1 * t528) * t501;
t498 = t501 ^ 2;
t465 = -t498 * t544 + t501 * t538 + t478;
t439 = t503 * t464 + t506 * t465;
t473 = pkin(3) * t489 - pkin(8) * t490;
t508 = qJD(3) ^ 2;
t415 = -pkin(3) * t508 + qJDD(3) * pkin(8) - t473 * t489 + t439;
t525 = t504 * g(1) - t507 * g(2);
t518 = qJDD(2) - t525;
t532 = -t500 ^ 2 - t498;
t474 = (-pkin(2) * t501 - pkin(1)) * qJDD(1) + (t532 * pkin(7) - qJ(2)) * t509 + t518;
t530 = t489 * qJD(3);
t476 = t516 * qJDD(1) - t530;
t423 = (-t476 + t530) * pkin(8) + (-t475 + t529) * pkin(3) + t474;
t502 = sin(qJ(4));
t505 = cos(qJ(4));
t411 = -t502 * t415 + t505 * t423;
t481 = qJD(3) * t505 - t490 * t502;
t452 = qJD(4) * t481 + qJDD(3) * t502 + t476 * t505;
t472 = qJDD(4) - t475;
t482 = qJD(3) * t502 + t490 * t505;
t487 = qJD(4) + t489;
t407 = (t481 * t487 - t452) * qJ(5) + (t481 * t482 + t472) * pkin(4) + t411;
t412 = t505 * t415 + t502 * t423;
t451 = -qJD(4) * t482 + qJDD(3) * t505 - t476 * t502;
t460 = pkin(4) * t487 - qJ(5) * t482;
t480 = t481 ^ 2;
t409 = -pkin(4) * t480 + qJ(5) * t451 - t460 * t487 + t412;
t499 = sin(pkin(10));
t539 = cos(pkin(10));
t454 = -t539 * t481 + t499 * t482;
t545 = -2 * qJD(5);
t403 = t499 * t407 + t539 * t409 + t454 * t545;
t419 = -t539 * t451 + t499 * t452;
t455 = t499 * t481 + t539 * t482;
t442 = mrSges(6,1) * t487 - mrSges(6,3) * t455;
t433 = pkin(5) * t454 - qJ(6) * t455;
t486 = t487 ^ 2;
t400 = -pkin(5) * t486 + qJ(6) * t472 + 0.2e1 * qJD(6) * t487 - t433 * t454 + t403;
t443 = -mrSges(7,1) * t487 + mrSges(7,2) * t455;
t526 = m(7) * t400 + t472 * mrSges(7,3) + t487 * t443;
t434 = mrSges(7,1) * t454 - mrSges(7,3) * t455;
t533 = -mrSges(6,1) * t454 - mrSges(6,2) * t455 - t434;
t543 = -mrSges(6,3) - mrSges(7,2);
t390 = m(6) * t403 - t472 * mrSges(6,2) + t543 * t419 - t487 * t442 + t533 * t454 + t526;
t513 = t539 * t407 - t499 * t409;
t402 = t455 * t545 + t513;
t420 = t499 * t451 + t539 * t452;
t441 = -mrSges(6,2) * t487 - mrSges(6,3) * t454;
t401 = -t472 * pkin(5) - t486 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t433) * t455 - t513;
t440 = -mrSges(7,2) * t454 + mrSges(7,3) * t487;
t520 = -m(7) * t401 + t472 * mrSges(7,1) + t487 * t440;
t392 = m(6) * t402 + t472 * mrSges(6,1) + t543 * t420 + t487 * t441 + t533 * t455 + t520;
t385 = t499 * t390 + t539 * t392;
t397 = t420 * mrSges(7,2) + t455 * t434 - t520;
t446 = Ifges(5,4) * t482 + Ifges(5,2) * t481 + Ifges(5,6) * t487;
t447 = Ifges(5,1) * t482 + Ifges(5,4) * t481 + Ifges(5,5) * t487;
t534 = -t542 * t454 + t549 * t455 + t541 * t487;
t535 = t548 * t454 + t542 * t455 + t540 * t487;
t546 = -t540 * t419 + t541 * t420 + t534 * t454 + t535 * t455 - (-Ifges(5,3) + t547) * t472 + mrSges(5,1) * t411 + mrSges(6,1) * t402 - mrSges(7,1) * t401 - mrSges(5,2) * t412 - mrSges(6,2) * t403 + mrSges(7,3) * t400 + Ifges(5,5) * t452 + Ifges(5,6) * t451 + pkin(4) * t385 - pkin(5) * t397 + qJ(6) * (-t419 * mrSges(7,2) - t454 * t434 + t526) + t482 * t446 - t481 * t447;
t471 = mrSges(4,1) * t489 + mrSges(4,2) * t490;
t484 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t490;
t456 = -mrSges(5,1) * t481 + mrSges(5,2) * t482;
t459 = -mrSges(5,2) * t487 + mrSges(5,3) * t481;
t383 = m(5) * t411 + mrSges(5,1) * t472 - mrSges(5,3) * t452 - t456 * t482 + t459 * t487 + t385;
t461 = mrSges(5,1) * t487 - mrSges(5,3) * t482;
t521 = t539 * t390 - t392 * t499;
t384 = m(5) * t412 - mrSges(5,2) * t472 + mrSges(5,3) * t451 + t456 * t481 - t461 * t487 + t521;
t522 = -t383 * t502 + t505 * t384;
t378 = m(4) * t439 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t475 - qJD(3) * t484 - t471 * t489 + t522;
t438 = t506 * t464 - t503 * t465;
t483 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t489;
t414 = -qJDD(3) * pkin(3) - t508 * pkin(8) + t490 * t473 - t438;
t410 = -t451 * pkin(4) - t480 * qJ(5) + t482 * t460 + qJDD(5) + t414;
t405 = -0.2e1 * qJD(6) * t455 + (t454 * t487 - t420) * qJ(6) + (t455 * t487 + t419) * pkin(5) + t410;
t398 = m(7) * t405 + t419 * mrSges(7,1) - t420 * mrSges(7,3) + t454 * t440 - t455 * t443;
t395 = m(6) * t410 + t419 * mrSges(6,1) + mrSges(6,2) * t420 + t454 * t441 + t442 * t455 + t398;
t511 = -m(5) * t414 + t451 * mrSges(5,1) - mrSges(5,2) * t452 + t481 * t459 - t461 * t482 - t395;
t394 = m(4) * t438 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t476 + qJD(3) * t483 - t471 * t490 + t511;
t537 = t503 * t378 + t506 * t394;
t379 = t505 * t383 + t502 * t384;
t536 = t540 * t454 - t541 * t455 + t547 * t487;
t523 = t506 * t378 - t503 * t394;
t517 = -mrSges(3,1) * t501 + mrSges(3,2) * t500;
t514 = mrSges(3,3) * qJDD(1) + t509 * t517;
t512 = m(4) * t474 - t475 * mrSges(4,1) + t476 * mrSges(4,2) + t489 * t483 + t490 * t484 + t379;
t488 = -qJDD(1) * pkin(1) - t509 * qJ(2) + t518;
t477 = -t500 * t491 + t524;
t468 = Ifges(4,1) * t490 - Ifges(4,4) * t489 + Ifges(4,5) * qJD(3);
t467 = Ifges(4,4) * t490 - Ifges(4,2) * t489 + Ifges(4,6) * qJD(3);
t466 = Ifges(4,5) * t490 - Ifges(4,6) * t489 + Ifges(4,3) * qJD(3);
t445 = Ifges(5,5) * t482 + Ifges(5,6) * t481 + Ifges(5,3) * t487;
t387 = mrSges(6,2) * t410 + mrSges(7,2) * t401 - mrSges(6,3) * t402 - mrSges(7,3) * t405 - qJ(6) * t398 - t542 * t419 + t549 * t420 + t536 * t454 + t541 * t472 - t535 * t487;
t386 = -mrSges(6,1) * t410 - mrSges(7,1) * t405 + mrSges(7,2) * t400 + mrSges(6,3) * t403 - pkin(5) * t398 + t548 * t419 + t542 * t420 + t536 * t455 + t540 * t472 + t534 * t487;
t375 = t532 * t509 * mrSges(3,3) + m(3) * t488 + t517 * qJDD(1) + t512;
t374 = mrSges(5,2) * t414 - mrSges(5,3) * t411 + Ifges(5,1) * t452 + Ifges(5,4) * t451 + Ifges(5,5) * t472 - qJ(5) * t385 - t499 * t386 + t539 * t387 + t481 * t445 - t487 * t446;
t373 = -mrSges(5,1) * t414 + mrSges(5,3) * t412 + Ifges(5,4) * t452 + Ifges(5,2) * t451 + Ifges(5,6) * t472 - pkin(4) * t395 + qJ(5) * t521 + t539 * t386 + t499 * t387 - t482 * t445 + t487 * t447;
t372 = -mrSges(4,1) * t474 + mrSges(4,3) * t439 + Ifges(4,4) * t476 + Ifges(4,2) * t475 + Ifges(4,6) * qJDD(3) - pkin(3) * t379 + qJD(3) * t468 - t490 * t466 - t546;
t371 = mrSges(4,2) * t474 - mrSges(4,3) * t438 + Ifges(4,1) * t476 + Ifges(4,4) * t475 + Ifges(4,5) * qJDD(3) - pkin(8) * t379 - qJD(3) * t467 - t373 * t502 + t374 * t505 - t466 * t489;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t525 - mrSges(2,2) * t519 + t500 * (mrSges(3,2) * t488 - mrSges(3,3) * t477 + t506 * t371 - t503 * t372 - pkin(7) * t537 + (Ifges(3,1) * t500 + Ifges(3,4) * t501) * qJDD(1)) + t501 * (-mrSges(3,1) * t488 + mrSges(3,3) * t478 + t503 * t371 + t506 * t372 - pkin(2) * t512 + pkin(7) * t523 + (Ifges(3,4) * t500 + Ifges(3,2) * t501) * qJDD(1)) - pkin(1) * t375 + qJ(2) * ((m(3) * t478 + t514 * t501 + t523) * t501 + (-m(3) * t477 + t514 * t500 - t537) * t500); t375; mrSges(4,1) * t438 - mrSges(4,2) * t439 + Ifges(4,5) * t476 + Ifges(4,6) * t475 + Ifges(4,3) * qJDD(3) + pkin(3) * t511 + pkin(8) * t522 + t505 * t373 + t502 * t374 + t490 * t467 + t489 * t468; t546; t395; t397;];
tauJ  = t1;
