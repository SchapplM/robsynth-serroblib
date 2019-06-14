% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-05-07 05:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:02:00
% EndTime: 2019-05-07 05:02:16
% DurationCPUTime: 15.45s
% Computational Cost: add. (241359->360), mult. (537750->472), div. (0->0), fcn. (435382->14), ass. (0->150)
t585 = -2 * qJD(4);
t546 = sin(pkin(6));
t551 = sin(qJ(2));
t555 = cos(qJ(2));
t572 = qJD(1) * qJD(2);
t536 = (-qJDD(1) * t555 + t551 * t572) * t546;
t575 = qJD(1) * t546;
t534 = (-t555 * pkin(2) - t551 * pkin(9)) * t575;
t548 = cos(pkin(6));
t541 = qJD(1) * t548 + qJD(2);
t539 = t541 ^ 2;
t540 = qJDD(1) * t548 + qJDD(2);
t574 = qJD(1) * t555;
t557 = qJD(1) ^ 2;
t552 = sin(qJ(1));
t556 = cos(qJ(1));
t565 = -g(1) * t556 - g(2) * t552;
t583 = pkin(8) * t546;
t532 = -pkin(1) * t557 + qJDD(1) * t583 + t565;
t569 = t552 * g(1) - g(2) * t556;
t531 = qJDD(1) * pkin(1) + t557 * t583 + t569;
t579 = t531 * t548;
t576 = t555 * t532 + t551 * t579;
t487 = -pkin(2) * t539 + pkin(9) * t540 + (-g(3) * t551 + t534 * t574) * t546 + t576;
t535 = (qJDD(1) * t551 + t555 * t572) * t546;
t582 = g(3) * t548;
t488 = pkin(2) * t536 - pkin(9) * t535 - t582 + (-t531 + (pkin(2) * t551 - pkin(9) * t555) * t541 * qJD(1)) * t546;
t550 = sin(qJ(3));
t554 = cos(qJ(3));
t458 = -t487 * t550 + t554 * t488;
t571 = t551 * t575;
t524 = t541 * t554 - t550 * t571;
t504 = qJD(3) * t524 + t535 * t554 + t540 * t550;
t525 = t541 * t550 + t554 * t571;
t528 = qJDD(3) + t536;
t570 = t546 * t574;
t538 = qJD(3) - t570;
t446 = (t524 * t538 - t504) * qJ(4) + (t524 * t525 + t528) * pkin(3) + t458;
t459 = t554 * t487 + t550 * t488;
t503 = -qJD(3) * t525 - t535 * t550 + t540 * t554;
t514 = pkin(3) * t538 - qJ(4) * t525;
t523 = t524 ^ 2;
t453 = -pkin(3) * t523 + qJ(4) * t503 - t514 * t538 + t459;
t545 = sin(pkin(11));
t580 = cos(pkin(11));
t511 = t545 * t524 + t580 * t525;
t437 = t580 * t446 - t545 * t453 + t511 * t585;
t510 = -t580 * t524 + t525 * t545;
t438 = t545 * t446 + t580 * t453 + t510 * t585;
t481 = pkin(4) * t510 - qJ(5) * t511;
t537 = t538 ^ 2;
t436 = -pkin(4) * t537 + qJ(5) * t528 - t481 * t510 + t438;
t577 = t546 * t555;
t505 = -g(3) * t577 - t551 * t532 + t555 * t579;
t486 = -pkin(2) * t540 - pkin(9) * t539 + t534 * t571 - t505;
t455 = -pkin(3) * t503 - qJ(4) * t523 + t525 * t514 + qJDD(4) + t486;
t475 = -t580 * t503 + t504 * t545;
t476 = t545 * t503 + t580 * t504;
t441 = (t510 * t538 - t476) * qJ(5) + (t511 * t538 + t475) * pkin(4) + t455;
t544 = sin(pkin(12));
t547 = cos(pkin(12));
t493 = t511 * t547 + t538 * t544;
t431 = -0.2e1 * qJD(5) * t493 - t436 * t544 + t547 * t441;
t467 = t476 * t547 + t528 * t544;
t492 = -t511 * t544 + t538 * t547;
t429 = (t492 * t510 - t467) * pkin(10) + (t492 * t493 + t475) * pkin(5) + t431;
t432 = 0.2e1 * qJD(5) * t492 + t547 * t436 + t544 * t441;
t466 = -t476 * t544 + t528 * t547;
t472 = pkin(5) * t510 - pkin(10) * t493;
t491 = t492 ^ 2;
t430 = -pkin(5) * t491 + pkin(10) * t466 - t472 * t510 + t432;
t549 = sin(qJ(6));
t553 = cos(qJ(6));
t427 = t429 * t553 - t430 * t549;
t468 = t492 * t553 - t493 * t549;
t444 = qJD(6) * t468 + t466 * t549 + t467 * t553;
t469 = t492 * t549 + t493 * t553;
t454 = -mrSges(7,1) * t468 + mrSges(7,2) * t469;
t509 = qJD(6) + t510;
t456 = -mrSges(7,2) * t509 + mrSges(7,3) * t468;
t474 = qJDD(6) + t475;
t424 = m(7) * t427 + mrSges(7,1) * t474 - mrSges(7,3) * t444 - t454 * t469 + t456 * t509;
t428 = t429 * t549 + t430 * t553;
t443 = -qJD(6) * t469 + t466 * t553 - t467 * t549;
t457 = mrSges(7,1) * t509 - mrSges(7,3) * t469;
t425 = m(7) * t428 - mrSges(7,2) * t474 + mrSges(7,3) * t443 + t454 * t468 - t457 * t509;
t416 = t553 * t424 + t549 * t425;
t470 = -mrSges(6,1) * t492 + mrSges(6,2) * t493;
t564 = -mrSges(6,2) * t510 + mrSges(6,3) * t492;
t414 = m(6) * t431 + t475 * mrSges(6,1) - t467 * mrSges(6,3) - t493 * t470 + t510 * t564 + t416;
t471 = mrSges(6,1) * t510 - mrSges(6,3) * t493;
t566 = -t424 * t549 + t553 * t425;
t415 = m(6) * t432 - mrSges(6,2) * t475 + mrSges(6,3) * t466 + t470 * t492 - t471 * t510 + t566;
t410 = -t414 * t544 + t547 * t415;
t482 = mrSges(5,1) * t510 + mrSges(5,2) * t511;
t495 = mrSges(5,1) * t538 - mrSges(5,3) * t511;
t407 = m(5) * t438 - mrSges(5,2) * t528 - mrSges(5,3) * t475 - t482 * t510 - t495 * t538 + t410;
t435 = -t528 * pkin(4) - t537 * qJ(5) + t511 * t481 + qJDD(5) - t437;
t433 = -t466 * pkin(5) - t491 * pkin(10) + t493 * t472 + t435;
t561 = m(7) * t433 - t443 * mrSges(7,1) + mrSges(7,2) * t444 - t468 * t456 + t457 * t469;
t426 = m(6) * t435 - t466 * mrSges(6,1) + mrSges(6,2) * t467 + t471 * t493 - t492 * t564 + t561;
t494 = -mrSges(5,2) * t538 - mrSges(5,3) * t510;
t420 = m(5) * t437 + mrSges(5,1) * t528 - mrSges(5,3) * t476 - t482 * t511 + t494 * t538 - t426;
t401 = t545 * t407 + t580 * t420;
t447 = Ifges(7,5) * t469 + Ifges(7,6) * t468 + Ifges(7,3) * t509;
t449 = Ifges(7,1) * t469 + Ifges(7,4) * t468 + Ifges(7,5) * t509;
t417 = -mrSges(7,1) * t433 + mrSges(7,3) * t428 + Ifges(7,4) * t444 + Ifges(7,2) * t443 + Ifges(7,6) * t474 - t447 * t469 + t449 * t509;
t448 = Ifges(7,4) * t469 + Ifges(7,2) * t468 + Ifges(7,6) * t509;
t418 = mrSges(7,2) * t433 - mrSges(7,3) * t427 + Ifges(7,1) * t444 + Ifges(7,4) * t443 + Ifges(7,5) * t474 + t447 * t468 - t448 * t509;
t460 = Ifges(6,5) * t493 + Ifges(6,6) * t492 + Ifges(6,3) * t510;
t462 = Ifges(6,1) * t493 + Ifges(6,4) * t492 + Ifges(6,5) * t510;
t402 = -mrSges(6,1) * t435 + mrSges(6,3) * t432 + Ifges(6,4) * t467 + Ifges(6,2) * t466 + Ifges(6,6) * t475 - pkin(5) * t561 + pkin(10) * t566 + t553 * t417 + t549 * t418 - t493 * t460 + t510 * t462;
t461 = Ifges(6,4) * t493 + Ifges(6,2) * t492 + Ifges(6,6) * t510;
t403 = mrSges(6,2) * t435 - mrSges(6,3) * t431 + Ifges(6,1) * t467 + Ifges(6,4) * t466 + Ifges(6,5) * t475 - pkin(10) * t416 - t417 * t549 + t418 * t553 + t460 * t492 - t461 * t510;
t478 = Ifges(5,4) * t511 - Ifges(5,2) * t510 + Ifges(5,6) * t538;
t479 = Ifges(5,1) * t511 - Ifges(5,4) * t510 + Ifges(5,5) * t538;
t498 = Ifges(4,4) * t525 + Ifges(4,2) * t524 + Ifges(4,6) * t538;
t499 = Ifges(4,1) * t525 + Ifges(4,4) * t524 + Ifges(4,5) * t538;
t584 = Ifges(4,5) * t504 + Ifges(4,6) * t503 + t525 * t498 - t524 * t499 + mrSges(4,1) * t458 - mrSges(4,2) * t459 + Ifges(5,5) * t476 - Ifges(5,6) * t475 + t511 * t478 + t510 * t479 + mrSges(5,1) * t437 - mrSges(5,2) * t438 + t544 * t403 + t547 * t402 - pkin(4) * t426 + qJ(5) * t410 + pkin(3) * t401 + (Ifges(4,3) + Ifges(5,3)) * t528;
t578 = t546 * t551;
t512 = -mrSges(4,1) * t524 + mrSges(4,2) * t525;
t513 = -mrSges(4,2) * t538 + mrSges(4,3) * t524;
t399 = m(4) * t458 + mrSges(4,1) * t528 - mrSges(4,3) * t504 - t512 * t525 + t513 * t538 + t401;
t515 = mrSges(4,1) * t538 - mrSges(4,3) * t525;
t567 = t580 * t407 - t420 * t545;
t400 = m(4) * t459 - mrSges(4,2) * t528 + mrSges(4,3) * t503 + t512 * t524 - t515 * t538 + t567;
t394 = t554 * t399 + t550 * t400;
t409 = t547 * t414 + t544 * t415;
t568 = -t399 * t550 + t554 * t400;
t408 = m(5) * t455 + t475 * mrSges(5,1) + mrSges(5,2) * t476 + t510 * t494 + t495 * t511 + t409;
t560 = mrSges(7,1) * t427 - mrSges(7,2) * t428 + Ifges(7,5) * t444 + Ifges(7,6) * t443 + Ifges(7,3) * t474 + t469 * t448 - t468 * t449;
t559 = -m(4) * t486 + t503 * mrSges(4,1) - mrSges(4,2) * t504 + t524 * t513 - t515 * t525 - t408;
t533 = (-t555 * mrSges(3,1) + t551 * mrSges(3,2)) * t575;
t530 = -mrSges(3,2) * t541 + mrSges(3,3) * t570;
t529 = mrSges(3,1) * t541 - mrSges(3,3) * t571;
t519 = -t531 * t546 - t582;
t518 = Ifges(3,5) * t541 + (t551 * Ifges(3,1) + t555 * Ifges(3,4)) * t575;
t517 = Ifges(3,6) * t541 + (t551 * Ifges(3,4) + t555 * Ifges(3,2)) * t575;
t516 = Ifges(3,3) * t541 + (t551 * Ifges(3,5) + t555 * Ifges(3,6)) * t575;
t506 = -g(3) * t578 + t576;
t497 = Ifges(4,5) * t525 + Ifges(4,6) * t524 + Ifges(4,3) * t538;
t477 = Ifges(5,5) * t511 - Ifges(5,6) * t510 + Ifges(5,3) * t538;
t404 = m(3) * t505 + mrSges(3,1) * t540 - mrSges(3,3) * t535 + t530 * t541 - t533 * t571 + t559;
t395 = -t560 + t538 * t479 + Ifges(5,6) * t528 - t511 * t477 - t493 * t461 + t492 * t462 + Ifges(5,4) * t476 - Ifges(6,6) * t466 - Ifges(6,5) * t467 - mrSges(5,1) * t455 + mrSges(5,3) * t438 - mrSges(6,1) * t431 + mrSges(6,2) * t432 - pkin(5) * t416 + (-Ifges(5,2) - Ifges(6,3)) * t475 - pkin(4) * t409;
t393 = m(3) * t506 - mrSges(3,2) * t540 - mrSges(3,3) * t536 - t529 * t541 + t533 * t570 + t568;
t392 = mrSges(5,2) * t455 - mrSges(5,3) * t437 + Ifges(5,1) * t476 - Ifges(5,4) * t475 + Ifges(5,5) * t528 - qJ(5) * t409 - t402 * t544 + t403 * t547 - t477 * t510 - t478 * t538;
t391 = mrSges(4,2) * t486 - mrSges(4,3) * t458 + Ifges(4,1) * t504 + Ifges(4,4) * t503 + Ifges(4,5) * t528 - qJ(4) * t401 + t580 * t392 - t545 * t395 + t524 * t497 - t538 * t498;
t390 = -mrSges(4,1) * t486 + mrSges(4,3) * t459 + Ifges(4,4) * t504 + Ifges(4,2) * t503 + Ifges(4,6) * t528 - pkin(3) * t408 + qJ(4) * t567 + t545 * t392 + t580 * t395 - t525 * t497 + t538 * t499;
t389 = Ifges(3,5) * t535 - Ifges(3,6) * t536 + Ifges(3,3) * t540 + mrSges(3,1) * t505 - mrSges(3,2) * t506 + t550 * t391 + t554 * t390 + pkin(2) * t559 + pkin(9) * t568 + (t551 * t517 - t555 * t518) * t575;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t569 - mrSges(2,2) * t565 + (mrSges(3,2) * t519 - mrSges(3,3) * t505 + Ifges(3,1) * t535 - Ifges(3,4) * t536 + Ifges(3,5) * t540 - pkin(9) * t394 - t390 * t550 + t391 * t554 + t516 * t570 - t517 * t541) * t578 + (-mrSges(3,1) * t519 + mrSges(3,3) * t506 + Ifges(3,4) * t535 - Ifges(3,2) * t536 + Ifges(3,6) * t540 - pkin(2) * t394 - t516 * t571 + t541 * t518 - t584) * t577 + t548 * t389 + pkin(1) * ((t551 * t393 + t555 * t404) * t548 + (-m(3) * t519 - t536 * mrSges(3,1) - t535 * mrSges(3,2) + (-t529 * t551 + t530 * t555) * t575 - t394) * t546) + (t393 * t555 - t404 * t551) * t583; t389; t584; t408; t426; t560;];
tauJ  = t1;
