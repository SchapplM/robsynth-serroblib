% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-05-05 04:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:07:12
% EndTime: 2019-05-05 04:07:15
% DurationCPUTime: 1.70s
% Computational Cost: add. (6310->248), mult. (12305->290), div. (0->0), fcn. (7263->10), ass. (0->111)
t549 = Ifges(6,1) + Ifges(7,1);
t528 = Ifges(6,4) - Ifges(7,5);
t544 = Ifges(7,4) + Ifges(6,5);
t548 = Ifges(6,2) + Ifges(7,3);
t543 = Ifges(6,6) - Ifges(7,6);
t490 = sin(qJ(3));
t493 = cos(qJ(3));
t529 = Ifges(4,4) + Ifges(5,6);
t547 = t490 * t529 + t493 * (Ifges(4,2) + Ifges(5,3));
t546 = t490 * (Ifges(4,1) + Ifges(5,2)) + t493 * t529;
t545 = -2 * qJD(4);
t527 = Ifges(4,5) - Ifges(5,4);
t526 = Ifges(4,6) - Ifges(5,5);
t542 = Ifges(6,3) + Ifges(7,2);
t489 = sin(qJ(5));
t492 = cos(qJ(5));
t516 = qJD(2) * t493;
t458 = qJD(3) * t489 + t492 * t516;
t459 = qJD(3) * t492 - t489 * t516;
t517 = qJD(2) * t490;
t477 = qJD(5) + t517;
t541 = t548 * t458 - t528 * t459 - t543 * t477;
t540 = -t528 * t458 + t549 * t459 + t544 * t477;
t485 = sin(pkin(10));
t487 = cos(pkin(10));
t467 = g(1) * t485 - g(2) * t487;
t484 = -g(3) + qJDD(1);
t486 = sin(pkin(6));
t488 = cos(pkin(6));
t539 = t467 * t488 + t484 * t486;
t468 = -g(1) * t487 - g(2) * t485;
t491 = sin(qJ(2));
t494 = cos(qJ(2));
t407 = t494 * t468 + t539 * t491;
t496 = qJD(2) ^ 2;
t405 = -pkin(2) * t496 + qJDD(2) * pkin(8) + t407;
t434 = -t467 * t486 + t484 * t488;
t399 = t493 * t405 + t490 * t434;
t460 = (-t493 * pkin(3) - t490 * qJ(4)) * qJD(2);
t495 = qJD(3) ^ 2;
t395 = pkin(3) * t495 - qJDD(3) * qJ(4) + qJD(3) * t545 - t460 * t516 - t399;
t406 = -t491 * t468 + t494 * t539;
t501 = -qJDD(2) * pkin(2) - t406;
t533 = pkin(8) * t496;
t404 = t501 - t533;
t515 = qJD(2) * qJD(3);
t512 = t493 * t515;
t463 = qJDD(2) * t490 + t512;
t511 = t490 * t515;
t464 = qJDD(2) * t493 - t511;
t469 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t517;
t470 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t516;
t471 = -mrSges(5,1) * t516 - qJD(3) * mrSges(5,3);
t499 = pkin(3) * t511 + t517 * t545 + (-t463 - t512) * qJ(4) + t501;
t397 = -pkin(3) * t464 + t499 - t533;
t472 = mrSges(5,1) * t517 + qJD(3) * mrSges(5,2);
t402 = t490 * t405;
t507 = -qJ(4) * t495 + t460 * t517 + qJDD(4) + t402;
t532 = pkin(9) * t496;
t534 = -pkin(3) - pkin(9);
t392 = pkin(4) * t463 + t534 * qJDD(3) + (-pkin(4) * t515 - t490 * t532 - t434) * t493 + t507;
t474 = pkin(4) * t517 - qJD(3) * pkin(9);
t483 = t493 ^ 2;
t394 = -t474 * t517 + (-pkin(4) * t483 - pkin(8)) * t496 + t534 * t464 + t499;
t387 = t489 * t392 + t492 * t394;
t422 = qJD(5) * t459 + qJDD(3) * t489 + t492 * t464;
t432 = mrSges(6,1) * t477 - mrSges(6,3) * t459;
t455 = qJDD(5) + t463;
t426 = pkin(5) * t458 - qJ(6) * t459;
t475 = t477 ^ 2;
t382 = -pkin(5) * t475 + qJ(6) * t455 + 0.2e1 * qJD(6) * t477 - t426 * t458 + t387;
t433 = -mrSges(7,1) * t477 + mrSges(7,2) * t459;
t513 = m(7) * t382 + t455 * mrSges(7,3) + t477 * t433;
t427 = mrSges(7,1) * t458 - mrSges(7,3) * t459;
t520 = -mrSges(6,1) * t458 - mrSges(6,2) * t459 - t427;
t530 = -mrSges(6,3) - mrSges(7,2);
t373 = m(6) * t387 - mrSges(6,2) * t455 + t422 * t530 - t432 * t477 + t458 * t520 + t513;
t386 = t392 * t492 - t394 * t489;
t423 = -qJD(5) * t458 + qJDD(3) * t492 - t464 * t489;
t430 = -mrSges(6,2) * t477 - mrSges(6,3) * t458;
t383 = -pkin(5) * t455 - qJ(6) * t475 + t426 * t459 + qJDD(6) - t386;
t431 = -mrSges(7,2) * t458 + mrSges(7,3) * t477;
t509 = -m(7) * t383 + t455 * mrSges(7,1) + t477 * t431;
t375 = m(6) * t386 + mrSges(6,1) * t455 + t423 * t530 + t430 * t477 + t459 * t520 + t509;
t522 = t492 * t373 - t489 * t375;
t506 = m(5) * t397 + t464 * mrSges(5,2) - t472 * t517 + t522;
t536 = -m(4) * t404 + t464 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t463 + t470 * t516 + (-t469 * t490 - t471 * t493) * qJD(2) - t506;
t535 = t490 * (t547 * qJD(2) + t526 * qJD(3)) - t493 * (qJD(2) * t546 + t527 * qJD(3));
t525 = t434 * t493;
t521 = t543 * t458 - t544 * t459 - t542 * t477;
t398 = -t402 + t525;
t461 = (t493 * mrSges(5,2) - t490 * mrSges(5,3)) * qJD(2);
t462 = (-t493 * mrSges(4,1) + t490 * mrSges(4,2)) * qJD(2);
t369 = t373 * t489 + t375 * t492;
t396 = -qJDD(3) * pkin(3) + t507 - t525;
t503 = -m(5) * t396 - t463 * mrSges(5,1) - t369;
t364 = m(4) * t398 - mrSges(4,3) * t463 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t470 - t471) * qJD(3) + (-t461 - t462) * t517 + t503;
t391 = pkin(4) * t464 + qJD(3) * t474 - t483 * t532 - t395;
t388 = -0.2e1 * qJD(6) * t459 + (t458 * t477 - t423) * qJ(6) + (t459 * t477 + t422) * pkin(5) + t391;
t379 = m(7) * t388 + t422 * mrSges(7,1) - t423 * mrSges(7,3) + t458 * t431 - t459 * t433;
t500 = m(6) * t391 + t422 * mrSges(6,1) + t423 * mrSges(6,2) + t458 * t430 + t459 * t432 + t379;
t498 = -m(5) * t395 + qJDD(3) * mrSges(5,3) + qJD(3) * t472 + t461 * t516 + t500;
t371 = m(4) * t399 - qJD(3) * t469 + (mrSges(4,3) + mrSges(5,1)) * t464 - qJDD(3) * mrSges(4,2) + t498 + t462 * t516;
t510 = -t490 * t364 + t493 * t371;
t378 = mrSges(7,2) * t423 + t427 * t459 - t509;
t497 = mrSges(6,1) * t386 - mrSges(7,1) * t383 - mrSges(6,2) * t387 + mrSges(7,3) * t382 - pkin(5) * t378 + qJ(6) * t513 - t541 * t459 + (-qJ(6) * t427 + t540) * t458 + t542 * t455 + t544 * t423 + (-qJ(6) * mrSges(7,2) - t543) * t422;
t368 = mrSges(6,2) * t391 + mrSges(7,2) * t383 - mrSges(6,3) * t386 - mrSges(7,3) * t388 - qJ(6) * t379 - t528 * t422 + t549 * t423 + t544 * t455 + t521 * t458 + t541 * t477;
t367 = -mrSges(6,1) * t391 - mrSges(7,1) * t388 + mrSges(7,2) * t382 + mrSges(6,3) * t387 - pkin(5) * t379 - t548 * t422 + t528 * t423 + t543 * t455 + t521 * t459 + t540 * t477;
t366 = qJDD(3) * mrSges(5,2) + qJD(3) * t471 + t461 * t517 - t503;
t365 = -mrSges(5,3) * t463 + t471 * t516 + t506;
t1 = [m(2) * t484 + t488 * (m(3) * t434 + t364 * t493 + t371 * t490) + (t491 * (m(3) * t407 - mrSges(3,1) * t496 - qJDD(2) * mrSges(3,2) + t510) + t494 * (m(3) * t406 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t496 + t536)) * t486; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t406 - mrSges(3,2) * t407 + t490 * (mrSges(5,1) * t396 + mrSges(4,2) * t404 - mrSges(4,3) * t398 - mrSges(5,3) * t397 + pkin(4) * t369 - qJ(4) * t365 + t497) + t493 * (-mrSges(4,1) * t404 - mrSges(5,1) * t395 + mrSges(5,2) * t397 + mrSges(4,3) * t399 - pkin(3) * t365 + pkin(4) * t500 - pkin(9) * t522 - t492 * t367 - t489 * t368) + pkin(8) * t510 + t547 * t464 + (t490 * t527 + t493 * t526) * qJDD(3) - t535 * qJD(3) + t546 * t463 + t536 * pkin(2); mrSges(4,1) * t398 - mrSges(4,2) * t399 + mrSges(5,2) * t396 - mrSges(5,3) * t395 + t492 * t368 - t489 * t367 - pkin(9) * t369 - pkin(3) * t366 + qJ(4) * t498 + (qJ(4) * mrSges(5,1) + t526) * t464 + t527 * t463 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) + t535 * qJD(2); t366; t497; t378;];
tauJ  = t1;
