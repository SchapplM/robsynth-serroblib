% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 09:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:42:11
% EndTime: 2019-05-07 09:42:29
% DurationCPUTime: 10.49s
% Computational Cost: add. (147327->350), mult. (344798->446), div. (0->0), fcn. (262782->12), ass. (0->138)
t543 = qJD(1) ^ 2;
t562 = pkin(2) * t543;
t537 = sin(qJ(1));
t542 = cos(qJ(1));
t552 = -g(1) * t542 - g(2) * t537;
t513 = -pkin(1) * t543 + qJDD(1) * pkin(7) + t552;
t536 = sin(qJ(2));
t561 = t513 * t536;
t541 = cos(qJ(2));
t558 = qJD(1) * qJD(2);
t516 = qJDD(1) * t536 + t541 * t558;
t480 = qJDD(2) * pkin(2) - pkin(8) * t516 - t561 + (pkin(8) * t558 + t536 * t562 - g(3)) * t541;
t501 = -g(3) * t536 + t541 * t513;
t517 = qJDD(1) * t541 - t536 * t558;
t560 = qJD(1) * t536;
t520 = qJD(2) * pkin(2) - pkin(8) * t560;
t530 = t541 ^ 2;
t481 = pkin(8) * t517 - qJD(2) * t520 - t530 * t562 + t501;
t535 = sin(qJ(3));
t540 = cos(qJ(3));
t456 = t540 * t480 - t481 * t535;
t510 = (-t535 * t536 + t540 * t541) * qJD(1);
t485 = qJD(3) * t510 + t516 * t540 + t517 * t535;
t511 = (t535 * t541 + t536 * t540) * qJD(1);
t528 = qJDD(2) + qJDD(3);
t529 = qJD(2) + qJD(3);
t440 = (t510 * t529 - t485) * qJ(4) + (t510 * t511 + t528) * pkin(3) + t456;
t457 = t535 * t480 + t540 * t481;
t484 = -qJD(3) * t511 - t516 * t535 + t517 * t540;
t503 = pkin(3) * t529 - qJ(4) * t511;
t506 = t510 ^ 2;
t442 = -pkin(3) * t506 + qJ(4) * t484 - t503 * t529 + t457;
t531 = sin(pkin(11));
t532 = cos(pkin(11));
t498 = t510 * t531 + t511 * t532;
t420 = -0.2e1 * qJD(4) * t498 + t532 * t440 - t442 * t531;
t462 = t484 * t531 + t485 * t532;
t497 = t510 * t532 - t511 * t531;
t416 = (t497 * t529 - t462) * pkin(9) + (t497 * t498 + t528) * pkin(4) + t420;
t421 = 0.2e1 * qJD(4) * t497 + t531 * t440 + t532 * t442;
t461 = t484 * t532 - t485 * t531;
t489 = pkin(4) * t529 - pkin(9) * t498;
t495 = t497 ^ 2;
t418 = -pkin(4) * t495 + pkin(9) * t461 - t489 * t529 + t421;
t534 = sin(qJ(5));
t539 = cos(qJ(5));
t413 = t534 * t416 + t539 * t418;
t474 = t497 * t534 + t498 * t539;
t431 = -qJD(5) * t474 + t461 * t539 - t462 * t534;
t473 = t497 * t539 - t498 * t534;
t450 = -mrSges(6,1) * t473 + mrSges(6,2) * t474;
t526 = qJD(5) + t529;
t466 = mrSges(6,1) * t526 - mrSges(6,3) * t474;
t525 = qJDD(5) + t528;
t451 = -pkin(5) * t473 - pkin(10) * t474;
t524 = t526 ^ 2;
t410 = -pkin(5) * t524 + pkin(10) * t525 + t451 * t473 + t413;
t557 = g(1) * t537 - t542 * g(2);
t551 = -qJDD(1) * pkin(1) - t557;
t486 = -pkin(2) * t517 + t520 * t560 + (-pkin(8) * t530 - pkin(7)) * t543 + t551;
t453 = -pkin(3) * t484 - qJ(4) * t506 + t511 * t503 + qJDD(4) + t486;
t426 = -pkin(4) * t461 - pkin(9) * t495 + t498 * t489 + t453;
t432 = qJD(5) * t473 + t461 * t534 + t462 * t539;
t414 = (-t473 * t526 - t432) * pkin(10) + t426 + (t474 * t526 - t431) * pkin(5);
t533 = sin(qJ(6));
t538 = cos(qJ(6));
t407 = -t410 * t533 + t414 * t538;
t463 = -t474 * t533 + t526 * t538;
t424 = qJD(6) * t463 + t432 * t538 + t525 * t533;
t430 = qJDD(6) - t431;
t464 = t474 * t538 + t526 * t533;
t443 = -mrSges(7,1) * t463 + mrSges(7,2) * t464;
t470 = qJD(6) - t473;
t444 = -mrSges(7,2) * t470 + mrSges(7,3) * t463;
t404 = m(7) * t407 + mrSges(7,1) * t430 - mrSges(7,3) * t424 - t443 * t464 + t444 * t470;
t408 = t410 * t538 + t414 * t533;
t423 = -qJD(6) * t464 - t432 * t533 + t525 * t538;
t445 = mrSges(7,1) * t470 - mrSges(7,3) * t464;
t405 = m(7) * t408 - mrSges(7,2) * t430 + mrSges(7,3) * t423 + t443 * t463 - t445 * t470;
t553 = -t404 * t533 + t538 * t405;
t391 = m(6) * t413 - mrSges(6,2) * t525 + mrSges(6,3) * t431 + t450 * t473 - t466 * t526 + t553;
t412 = t416 * t539 - t418 * t534;
t465 = -mrSges(6,2) * t526 + mrSges(6,3) * t473;
t409 = -pkin(5) * t525 - pkin(10) * t524 + t451 * t474 - t412;
t550 = -m(7) * t409 + t423 * mrSges(7,1) - mrSges(7,2) * t424 + t463 * t444 - t445 * t464;
t400 = m(6) * t412 + mrSges(6,1) * t525 - mrSges(6,3) * t432 - t450 * t474 + t465 * t526 + t550;
t388 = t534 * t391 + t539 * t400;
t475 = -mrSges(5,1) * t497 + mrSges(5,2) * t498;
t487 = -mrSges(5,2) * t529 + mrSges(5,3) * t497;
t385 = m(5) * t420 + mrSges(5,1) * t528 - mrSges(5,3) * t462 - t475 * t498 + t487 * t529 + t388;
t488 = mrSges(5,1) * t529 - mrSges(5,3) * t498;
t554 = t539 * t391 - t400 * t534;
t386 = m(5) * t421 - mrSges(5,2) * t528 + mrSges(5,3) * t461 + t475 * t497 - t488 * t529 + t554;
t379 = t532 * t385 + t531 * t386;
t499 = -mrSges(4,1) * t510 + mrSges(4,2) * t511;
t502 = -mrSges(4,2) * t529 + mrSges(4,3) * t510;
t376 = m(4) * t456 + mrSges(4,1) * t528 - mrSges(4,3) * t485 - t499 * t511 + t502 * t529 + t379;
t504 = mrSges(4,1) * t529 - mrSges(4,3) * t511;
t555 = -t385 * t531 + t532 * t386;
t377 = m(4) * t457 - mrSges(4,2) * t528 + mrSges(4,3) * t484 + t499 * t510 - t504 * t529 + t555;
t371 = t540 * t376 + t535 * t377;
t394 = t538 * t404 + t533 * t405;
t559 = qJD(1) * t541;
t556 = -t376 * t535 + t540 * t377;
t549 = -m(6) * t426 + t431 * mrSges(6,1) - t432 * mrSges(6,2) + t473 * t465 - t474 * t466 - t394;
t433 = Ifges(7,5) * t464 + Ifges(7,6) * t463 + Ifges(7,3) * t470;
t435 = Ifges(7,1) * t464 + Ifges(7,4) * t463 + Ifges(7,5) * t470;
t397 = -mrSges(7,1) * t409 + mrSges(7,3) * t408 + Ifges(7,4) * t424 + Ifges(7,2) * t423 + Ifges(7,6) * t430 - t433 * t464 + t435 * t470;
t434 = Ifges(7,4) * t464 + Ifges(7,2) * t463 + Ifges(7,6) * t470;
t398 = mrSges(7,2) * t409 - mrSges(7,3) * t407 + Ifges(7,1) * t424 + Ifges(7,4) * t423 + Ifges(7,5) * t430 + t433 * t463 - t434 * t470;
t447 = Ifges(6,4) * t474 + Ifges(6,2) * t473 + Ifges(6,6) * t526;
t448 = Ifges(6,1) * t474 + Ifges(6,4) * t473 + Ifges(6,5) * t526;
t548 = mrSges(6,1) * t412 - mrSges(6,2) * t413 + Ifges(6,5) * t432 + Ifges(6,6) * t431 + Ifges(6,3) * t525 + pkin(5) * t550 + pkin(10) * t553 + t538 * t397 + t533 * t398 + t474 * t447 - t473 * t448;
t547 = mrSges(7,1) * t407 - mrSges(7,2) * t408 + Ifges(7,5) * t424 + Ifges(7,6) * t423 + Ifges(7,3) * t430 + t434 * t464 - t435 * t463;
t546 = -m(5) * t453 + t461 * mrSges(5,1) - t462 * mrSges(5,2) + t497 * t487 - t498 * t488 + t549;
t545 = m(4) * t486 - t484 * mrSges(4,1) + t485 * mrSges(4,2) - t510 * t502 + t511 * t504 - t546;
t468 = Ifges(5,4) * t498 + Ifges(5,2) * t497 + Ifges(5,6) * t529;
t469 = Ifges(5,1) * t498 + Ifges(5,4) * t497 + Ifges(5,5) * t529;
t493 = Ifges(4,4) * t511 + Ifges(4,2) * t510 + Ifges(4,6) * t529;
t494 = Ifges(4,1) * t511 + Ifges(4,4) * t510 + Ifges(4,5) * t529;
t544 = mrSges(4,1) * t456 + mrSges(5,1) * t420 - mrSges(4,2) * t457 - mrSges(5,2) * t421 + Ifges(4,5) * t485 + Ifges(5,5) * t462 + Ifges(4,6) * t484 + Ifges(5,6) * t461 + pkin(3) * t379 + pkin(4) * t388 + t498 * t468 - t497 * t469 + t511 * t493 - t510 * t494 + t548 + (Ifges(5,3) + Ifges(4,3)) * t528;
t519 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t559;
t518 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t560;
t515 = (-mrSges(3,1) * t541 + mrSges(3,2) * t536) * qJD(1);
t512 = -pkin(7) * t543 + t551;
t509 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t536 + Ifges(3,4) * t541) * qJD(1);
t508 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t536 + Ifges(3,2) * t541) * qJD(1);
t500 = -g(3) * t541 - t561;
t492 = Ifges(4,5) * t511 + Ifges(4,6) * t510 + Ifges(4,3) * t529;
t467 = Ifges(5,5) * t498 + Ifges(5,6) * t497 + Ifges(5,3) * t529;
t446 = Ifges(6,5) * t474 + Ifges(6,6) * t473 + Ifges(6,3) * t526;
t381 = -mrSges(6,1) * t426 + mrSges(6,3) * t413 + Ifges(6,4) * t432 + Ifges(6,2) * t431 + Ifges(6,6) * t525 - pkin(5) * t394 - t446 * t474 + t448 * t526 - t547;
t380 = mrSges(6,2) * t426 - mrSges(6,3) * t412 + Ifges(6,1) * t432 + Ifges(6,4) * t431 + Ifges(6,5) * t525 - pkin(10) * t394 - t397 * t533 + t398 * t538 + t446 * t473 - t447 * t526;
t372 = mrSges(5,2) * t453 - mrSges(5,3) * t420 + Ifges(5,1) * t462 + Ifges(5,4) * t461 + Ifges(5,5) * t528 - pkin(9) * t388 + t380 * t539 - t381 * t534 + t467 * t497 - t468 * t529;
t370 = -mrSges(5,1) * t453 + mrSges(5,3) * t421 + Ifges(5,4) * t462 + Ifges(5,2) * t461 + Ifges(5,6) * t528 + pkin(4) * t549 + pkin(9) * t554 + t534 * t380 + t539 * t381 - t498 * t467 + t529 * t469;
t369 = mrSges(4,2) * t486 - mrSges(4,3) * t456 + Ifges(4,1) * t485 + Ifges(4,4) * t484 + Ifges(4,5) * t528 - qJ(4) * t379 - t370 * t531 + t372 * t532 + t492 * t510 - t493 * t529;
t368 = -mrSges(4,1) * t486 + mrSges(4,3) * t457 + Ifges(4,4) * t485 + Ifges(4,2) * t484 + Ifges(4,6) * t528 + pkin(3) * t546 + qJ(4) * t555 + t532 * t370 + t531 * t372 - t511 * t492 + t529 * t494;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t557 - mrSges(2,2) * t552 + t536 * (mrSges(3,2) * t512 - mrSges(3,3) * t500 + Ifges(3,1) * t516 + Ifges(3,4) * t517 + Ifges(3,5) * qJDD(2) - pkin(8) * t371 - qJD(2) * t508 - t368 * t535 + t369 * t540) + t541 * (-mrSges(3,1) * t512 + mrSges(3,3) * t501 + Ifges(3,4) * t516 + Ifges(3,2) * t517 + Ifges(3,6) * qJDD(2) - pkin(2) * t545 + pkin(8) * t556 + qJD(2) * t509 + t540 * t368 + t535 * t369) + pkin(1) * ((-t518 * t536 + t519 * t541) * qJD(1) - t516 * mrSges(3,2) + t517 * mrSges(3,1) - m(3) * t512 - t545) + pkin(7) * (t541 * (m(3) * t501 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t517 - qJD(2) * t518 + t515 * t559 + t556) - t536 * (m(3) * t500 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t516 + qJD(2) * t519 - t515 * t560 + t371)); Ifges(3,3) * qJDD(2) + t544 + pkin(2) * t371 + Ifges(3,5) * t516 + Ifges(3,6) * t517 + mrSges(3,1) * t500 - mrSges(3,2) * t501 + (t536 * t508 - t541 * t509) * qJD(1); t544; -t546; t548; t547;];
tauJ  = t1;
