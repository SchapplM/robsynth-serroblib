% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 20:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:47:37
% EndTime: 2019-05-07 20:47:55
% DurationCPUTime: 13.11s
% Computational Cost: add. (208812->349), mult. (436805->440), div. (0->0), fcn. (324691->12), ass. (0->137)
t541 = qJD(1) ^ 2;
t534 = sin(qJ(1));
t539 = cos(qJ(1));
t553 = t534 * g(1) - t539 * g(2);
t509 = -qJDD(1) * pkin(1) - t541 * pkin(7) - t553;
t533 = sin(qJ(2));
t538 = cos(qJ(2));
t555 = qJD(1) * qJD(2);
t554 = t538 * t555;
t518 = qJDD(1) * t533 + t554;
t525 = t533 * t555;
t519 = qJDD(1) * t538 - t525;
t475 = (-t518 - t554) * pkin(8) + (-t519 + t525) * pkin(2) + t509;
t548 = -g(1) * t539 - g(2) * t534;
t510 = -pkin(1) * t541 + qJDD(1) * pkin(7) + t548;
t497 = -g(3) * t533 + t538 * t510;
t517 = (-pkin(2) * t538 - pkin(8) * t533) * qJD(1);
t540 = qJD(2) ^ 2;
t556 = qJD(1) * t538;
t478 = -pkin(2) * t540 + qJDD(2) * pkin(8) + t517 * t556 + t497;
t532 = sin(qJ(3));
t537 = cos(qJ(3));
t457 = t537 * t475 - t532 * t478;
t557 = qJD(1) * t533;
t514 = qJD(2) * t537 - t532 * t557;
t489 = qJD(3) * t514 + qJDD(2) * t532 + t518 * t537;
t513 = qJDD(3) - t519;
t515 = qJD(2) * t532 + t537 * t557;
t524 = qJD(3) - t556;
t437 = (t514 * t524 - t489) * pkin(9) + (t514 * t515 + t513) * pkin(3) + t457;
t458 = t532 * t475 + t537 * t478;
t488 = -qJD(3) * t515 + qJDD(2) * t537 - t518 * t532;
t498 = pkin(3) * t524 - pkin(9) * t515;
t512 = t514 ^ 2;
t439 = -pkin(3) * t512 + pkin(9) * t488 - t498 * t524 + t458;
t531 = sin(qJ(4));
t536 = cos(qJ(4));
t418 = t536 * t437 - t531 * t439;
t491 = t514 * t536 - t515 * t531;
t456 = qJD(4) * t491 + t488 * t531 + t489 * t536;
t492 = t514 * t531 + t515 * t536;
t511 = qJDD(4) + t513;
t523 = qJD(4) + t524;
t407 = (t491 * t523 - t456) * qJ(5) + (t491 * t492 + t511) * pkin(4) + t418;
t419 = t531 * t437 + t536 * t439;
t455 = -qJD(4) * t492 + t488 * t536 - t489 * t531;
t480 = pkin(4) * t523 - qJ(5) * t492;
t490 = t491 ^ 2;
t415 = -pkin(4) * t490 + qJ(5) * t455 - t480 * t523 + t419;
t528 = sin(pkin(11));
t529 = cos(pkin(11));
t471 = t491 * t528 + t492 * t529;
t401 = -0.2e1 * qJD(5) * t471 + t529 * t407 - t528 * t415;
t434 = t455 * t528 + t456 * t529;
t470 = t491 * t529 - t492 * t528;
t398 = (t470 * t523 - t434) * pkin(10) + (t470 * t471 + t511) * pkin(5) + t401;
t402 = 0.2e1 * qJD(5) * t470 + t528 * t407 + t529 * t415;
t433 = t455 * t529 - t456 * t528;
t461 = pkin(5) * t523 - pkin(10) * t471;
t469 = t470 ^ 2;
t399 = -pkin(5) * t469 + pkin(10) * t433 - t461 * t523 + t402;
t530 = sin(qJ(6));
t535 = cos(qJ(6));
t396 = t398 * t535 - t399 * t530;
t447 = t470 * t535 - t471 * t530;
t413 = qJD(6) * t447 + t433 * t530 + t434 * t535;
t448 = t470 * t530 + t471 * t535;
t427 = -mrSges(7,1) * t447 + mrSges(7,2) * t448;
t520 = qJD(6) + t523;
t440 = -mrSges(7,2) * t520 + mrSges(7,3) * t447;
t505 = qJDD(6) + t511;
t392 = m(7) * t396 + mrSges(7,1) * t505 - mrSges(7,3) * t413 - t427 * t448 + t440 * t520;
t397 = t398 * t530 + t399 * t535;
t412 = -qJD(6) * t448 + t433 * t535 - t434 * t530;
t441 = mrSges(7,1) * t520 - mrSges(7,3) * t448;
t393 = m(7) * t397 - mrSges(7,2) * t505 + mrSges(7,3) * t412 + t427 * t447 - t441 * t520;
t386 = t535 * t392 + t530 * t393;
t449 = -mrSges(6,1) * t470 + mrSges(6,2) * t471;
t459 = -mrSges(6,2) * t523 + mrSges(6,3) * t470;
t383 = m(6) * t401 + mrSges(6,1) * t511 - mrSges(6,3) * t434 - t449 * t471 + t459 * t523 + t386;
t460 = mrSges(6,1) * t523 - mrSges(6,3) * t471;
t549 = -t392 * t530 + t535 * t393;
t384 = m(6) * t402 - mrSges(6,2) * t511 + mrSges(6,3) * t433 + t449 * t470 - t460 * t523 + t549;
t379 = t529 * t383 + t528 * t384;
t472 = -mrSges(5,1) * t491 + mrSges(5,2) * t492;
t479 = -mrSges(5,2) * t523 + mrSges(5,3) * t491;
t376 = m(5) * t418 + mrSges(5,1) * t511 - mrSges(5,3) * t456 - t472 * t492 + t479 * t523 + t379;
t481 = mrSges(5,1) * t523 - mrSges(5,3) * t492;
t550 = -t383 * t528 + t529 * t384;
t377 = m(5) * t419 - mrSges(5,2) * t511 + mrSges(5,3) * t455 + t472 * t491 - t481 * t523 + t550;
t370 = t536 * t376 + t531 * t377;
t496 = -t538 * g(3) - t533 * t510;
t493 = -mrSges(4,1) * t514 + mrSges(4,2) * t515;
t494 = -mrSges(4,2) * t524 + mrSges(4,3) * t514;
t368 = m(4) * t457 + mrSges(4,1) * t513 - mrSges(4,3) * t489 - t493 * t515 + t494 * t524 + t370;
t495 = mrSges(4,1) * t524 - mrSges(4,3) * t515;
t551 = -t376 * t531 + t536 * t377;
t369 = m(4) * t458 - mrSges(4,2) * t513 + mrSges(4,3) * t488 + t493 * t514 - t495 * t524 + t551;
t552 = -t368 * t532 + t537 * t369;
t477 = -qJDD(2) * pkin(2) - pkin(8) * t540 + t517 * t557 - t496;
t450 = -pkin(3) * t488 - pkin(9) * t512 + t515 * t498 + t477;
t421 = -pkin(4) * t455 - qJ(5) * t490 + t492 * t480 + qJDD(5) + t450;
t404 = -pkin(5) * t433 - pkin(10) * t469 + t461 * t471 + t421;
t547 = m(7) * t404 - t412 * mrSges(7,1) + t413 * mrSges(7,2) - t447 * t440 + t448 * t441;
t364 = t537 * t368 + t532 * t369;
t394 = m(6) * t421 - t433 * mrSges(6,1) + t434 * mrSges(6,2) - t470 * t459 + t471 * t460 + t547;
t423 = Ifges(7,4) * t448 + Ifges(7,2) * t447 + Ifges(7,6) * t520;
t424 = Ifges(7,1) * t448 + Ifges(7,4) * t447 + Ifges(7,5) * t520;
t546 = -mrSges(7,1) * t396 + mrSges(7,2) * t397 - Ifges(7,5) * t413 - Ifges(7,6) * t412 - Ifges(7,3) * t505 - t448 * t423 + t447 * t424;
t545 = m(5) * t450 - t455 * mrSges(5,1) + t456 * mrSges(5,2) - t491 * t479 + t492 * t481 + t394;
t544 = -m(4) * t477 + t488 * mrSges(4,1) - t489 * mrSges(4,2) + t514 * t494 - t515 * t495 - t545;
t443 = Ifges(6,4) * t471 + Ifges(6,2) * t470 + Ifges(6,6) * t523;
t444 = Ifges(6,1) * t471 + Ifges(6,4) * t470 + Ifges(6,5) * t523;
t463 = Ifges(5,4) * t492 + Ifges(5,2) * t491 + Ifges(5,6) * t523;
t464 = Ifges(5,1) * t492 + Ifges(5,4) * t491 + Ifges(5,5) * t523;
t543 = -mrSges(5,1) * t418 - mrSges(6,1) * t401 + mrSges(5,2) * t419 + mrSges(6,2) * t402 - Ifges(5,5) * t456 - Ifges(6,5) * t434 - Ifges(5,6) * t455 - Ifges(6,6) * t433 - pkin(4) * t379 - pkin(5) * t386 - t471 * t443 + t470 * t444 - t492 * t463 + t491 * t464 + t546 + (-Ifges(6,3) - Ifges(5,3)) * t511;
t483 = Ifges(4,4) * t515 + Ifges(4,2) * t514 + Ifges(4,6) * t524;
t484 = Ifges(4,1) * t515 + Ifges(4,4) * t514 + Ifges(4,5) * t524;
t542 = mrSges(4,1) * t457 - mrSges(4,2) * t458 + Ifges(4,5) * t489 + Ifges(4,6) * t488 + Ifges(4,3) * t513 + pkin(3) * t370 + t515 * t483 - t514 * t484 - t543;
t522 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t556;
t521 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t557;
t516 = (-mrSges(3,1) * t538 + mrSges(3,2) * t533) * qJD(1);
t508 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t533 + Ifges(3,4) * t538) * qJD(1);
t507 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t533 + Ifges(3,2) * t538) * qJD(1);
t482 = Ifges(4,5) * t515 + Ifges(4,6) * t514 + Ifges(4,3) * t524;
t462 = Ifges(5,5) * t492 + Ifges(5,6) * t491 + Ifges(5,3) * t523;
t442 = Ifges(6,5) * t471 + Ifges(6,6) * t470 + Ifges(6,3) * t523;
t422 = Ifges(7,5) * t448 + Ifges(7,6) * t447 + Ifges(7,3) * t520;
t388 = mrSges(7,2) * t404 - mrSges(7,3) * t396 + Ifges(7,1) * t413 + Ifges(7,4) * t412 + Ifges(7,5) * t505 + t422 * t447 - t423 * t520;
t387 = -mrSges(7,1) * t404 + mrSges(7,3) * t397 + Ifges(7,4) * t413 + Ifges(7,2) * t412 + Ifges(7,6) * t505 - t422 * t448 + t424 * t520;
t372 = mrSges(6,2) * t421 - mrSges(6,3) * t401 + Ifges(6,1) * t434 + Ifges(6,4) * t433 + Ifges(6,5) * t511 - pkin(10) * t386 - t387 * t530 + t388 * t535 + t442 * t470 - t443 * t523;
t371 = -mrSges(6,1) * t421 + mrSges(6,3) * t402 + Ifges(6,4) * t434 + Ifges(6,2) * t433 + Ifges(6,6) * t511 - pkin(5) * t547 + pkin(10) * t549 + t535 * t387 + t530 * t388 - t471 * t442 + t523 * t444;
t366 = mrSges(5,2) * t450 - mrSges(5,3) * t418 + Ifges(5,1) * t456 + Ifges(5,4) * t455 + Ifges(5,5) * t511 - qJ(5) * t379 - t371 * t528 + t372 * t529 + t462 * t491 - t463 * t523;
t365 = -mrSges(5,1) * t450 + mrSges(5,3) * t419 + Ifges(5,4) * t456 + Ifges(5,2) * t455 + Ifges(5,6) * t511 - pkin(4) * t394 + qJ(5) * t550 + t529 * t371 + t528 * t372 - t492 * t462 + t523 * t464;
t363 = mrSges(4,2) * t477 - mrSges(4,3) * t457 + Ifges(4,1) * t489 + Ifges(4,4) * t488 + Ifges(4,5) * t513 - pkin(9) * t370 - t365 * t531 + t366 * t536 + t482 * t514 - t483 * t524;
t362 = -mrSges(4,1) * t477 + mrSges(4,3) * t458 + Ifges(4,4) * t489 + Ifges(4,2) * t488 + Ifges(4,6) * t513 - pkin(3) * t545 + pkin(9) * t551 + t536 * t365 + t531 * t366 - t515 * t482 + t524 * t484;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t553 - mrSges(2,2) * t548 + t533 * (mrSges(3,2) * t509 - mrSges(3,3) * t496 + Ifges(3,1) * t518 + Ifges(3,4) * t519 + Ifges(3,5) * qJDD(2) - pkin(8) * t364 - qJD(2) * t507 - t532 * t362 + t537 * t363) + t538 * (-mrSges(3,1) * t509 + mrSges(3,3) * t497 + Ifges(3,4) * t518 + Ifges(3,2) * t519 + Ifges(3,6) * qJDD(2) - pkin(2) * t364 + qJD(2) * t508 - t542) + pkin(1) * (-m(3) * t509 + t519 * mrSges(3,1) - t518 * mrSges(3,2) + (-t521 * t533 + t522 * t538) * qJD(1) - t364) + pkin(7) * (t538 * (m(3) * t497 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t519 - qJD(2) * t521 + t516 * t556 + t552) - t533 * (m(3) * t496 + qJDD(2) * mrSges(3,1) - t518 * mrSges(3,3) + qJD(2) * t522 - t516 * t557 + t544)); Ifges(3,5) * t518 + Ifges(3,6) * t519 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t496 - mrSges(3,2) * t497 + t532 * t363 + t537 * t362 + pkin(2) * t544 + pkin(8) * t552 + (t533 * t507 - t538 * t508) * qJD(1); t542; -t543; t394; -t546;];
tauJ  = t1;
