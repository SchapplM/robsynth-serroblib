% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 14:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:57:38
% EndTime: 2019-05-06 13:57:51
% DurationCPUTime: 12.27s
% Computational Cost: add. (170633->360), mult. (452149->472), div. (0->0), fcn. (360544->14), ass. (0->149)
t580 = -2 * qJD(3);
t541 = sin(pkin(11));
t544 = cos(pkin(11));
t548 = sin(qJ(2));
t551 = cos(qJ(2));
t542 = sin(pkin(6));
t573 = qJD(1) * t542;
t521 = (t548 * t541 - t551 * t544) * t573;
t571 = qJD(1) * qJD(2);
t530 = (qJDD(1) * t548 + t551 * t571) * t542;
t545 = cos(pkin(6));
t535 = qJDD(1) * t545 + qJDD(2);
t536 = qJD(1) * t545 + qJD(2);
t553 = qJD(1) ^ 2;
t549 = sin(qJ(1));
t552 = cos(qJ(1));
t561 = -g(1) * t552 - g(2) * t549;
t578 = pkin(8) * t542;
t528 = -pkin(1) * t553 + qJDD(1) * t578 + t561;
t567 = t549 * g(1) - g(2) * t552;
t527 = qJDD(1) * pkin(1) + t553 * t578 + t567;
t577 = t527 * t545;
t562 = -t528 * t548 + t551 * t577;
t576 = t542 ^ 2 * t553;
t469 = pkin(2) * t535 - qJ(3) * t530 + (pkin(2) * t548 * t576 + (qJ(3) * qJD(1) * t536 - g(3)) * t542) * t551 + t562;
t575 = t542 * t548;
t498 = -g(3) * t575 + t551 * t528 + t548 * t577;
t569 = t548 * t573;
t524 = pkin(2) * t536 - qJ(3) * t569;
t531 = (qJDD(1) * t551 - t548 * t571) * t542;
t570 = t551 ^ 2 * t576;
t472 = -pkin(2) * t570 + qJ(3) * t531 - t524 * t536 + t498;
t522 = (t551 * t541 + t548 * t544) * t573;
t449 = t469 * t544 - t541 * t472 + t522 * t580;
t579 = cos(qJ(4));
t574 = t542 * t551;
t450 = t541 * t469 + t544 * t472 + t521 * t580;
t499 = mrSges(4,1) * t521 + mrSges(4,2) * t522;
t503 = -t530 * t541 + t531 * t544;
t509 = mrSges(4,1) * t536 - mrSges(4,3) * t522;
t500 = pkin(3) * t521 - pkin(9) * t522;
t534 = t536 ^ 2;
t443 = -pkin(3) * t534 + pkin(9) * t535 - t500 * t521 + t450;
t513 = -g(3) * t545 - t527 * t542;
t482 = -pkin(2) * t531 - qJ(3) * t570 + t524 * t569 + qJDD(3) + t513;
t504 = t530 * t544 + t531 * t541;
t452 = (t521 * t536 - t504) * pkin(9) + (t522 * t536 - t503) * pkin(3) + t482;
t547 = sin(qJ(4));
t436 = t579 * t443 + t547 * t452;
t506 = t522 * t547 - t579 * t536;
t507 = t522 * t579 + t547 * t536;
t483 = pkin(4) * t506 - qJ(5) * t507;
t502 = qJDD(4) - t503;
t520 = qJD(4) + t521;
t519 = t520 ^ 2;
t431 = -pkin(4) * t519 + qJ(5) * t502 - t483 * t506 + t436;
t442 = -pkin(3) * t535 - pkin(9) * t534 + t522 * t500 - t449;
t479 = qJD(4) * t507 + t504 * t547 - t579 * t535;
t480 = -t506 * qJD(4) + t504 * t579 + t547 * t535;
t434 = (t506 * t520 - t480) * qJ(5) + (t507 * t520 + t479) * pkin(4) + t442;
t540 = sin(pkin(12));
t543 = cos(pkin(12));
t490 = t507 * t543 + t520 * t540;
t426 = -0.2e1 * qJD(5) * t490 - t431 * t540 + t543 * t434;
t461 = t480 * t543 + t502 * t540;
t489 = -t507 * t540 + t520 * t543;
t424 = (t489 * t506 - t461) * pkin(10) + (t489 * t490 + t479) * pkin(5) + t426;
t427 = 0.2e1 * qJD(5) * t489 + t543 * t431 + t540 * t434;
t460 = -t480 * t540 + t502 * t543;
t471 = pkin(5) * t506 - pkin(10) * t490;
t488 = t489 ^ 2;
t425 = -pkin(5) * t488 + pkin(10) * t460 - t471 * t506 + t427;
t546 = sin(qJ(6));
t550 = cos(qJ(6));
t422 = t424 * t550 - t425 * t546;
t462 = t489 * t550 - t490 * t546;
t439 = qJD(6) * t462 + t460 * t546 + t461 * t550;
t463 = t489 * t546 + t490 * t550;
t448 = -mrSges(7,1) * t462 + mrSges(7,2) * t463;
t505 = qJD(6) + t506;
t453 = -mrSges(7,2) * t505 + mrSges(7,3) * t462;
t478 = qJDD(6) + t479;
t419 = m(7) * t422 + mrSges(7,1) * t478 - mrSges(7,3) * t439 - t448 * t463 + t453 * t505;
t423 = t424 * t546 + t425 * t550;
t438 = -qJD(6) * t463 + t460 * t550 - t461 * t546;
t454 = mrSges(7,1) * t505 - mrSges(7,3) * t463;
t420 = m(7) * t423 - mrSges(7,2) * t478 + mrSges(7,3) * t438 + t448 * t462 - t454 * t505;
t411 = t550 * t419 + t546 * t420;
t464 = -mrSges(6,1) * t489 + mrSges(6,2) * t490;
t560 = -mrSges(6,2) * t506 + mrSges(6,3) * t489;
t409 = m(6) * t426 + t479 * mrSges(6,1) - t461 * mrSges(6,3) - t490 * t464 + t506 * t560 + t411;
t470 = mrSges(6,1) * t506 - mrSges(6,3) * t490;
t564 = -t419 * t546 + t550 * t420;
t410 = m(6) * t427 - mrSges(6,2) * t479 + mrSges(6,3) * t460 + t464 * t489 - t470 * t506 + t564;
t407 = -t409 * t540 + t543 * t410;
t484 = mrSges(5,1) * t506 + mrSges(5,2) * t507;
t492 = mrSges(5,1) * t520 - mrSges(5,3) * t507;
t405 = m(5) * t436 - mrSges(5,2) * t502 - mrSges(5,3) * t479 - t484 * t506 - t492 * t520 + t407;
t435 = -t547 * t443 + t452 * t579;
t430 = -t502 * pkin(4) - t519 * qJ(5) + t507 * t483 + qJDD(5) - t435;
t428 = -t460 * pkin(5) - t488 * pkin(10) + t490 * t471 + t430;
t558 = m(7) * t428 - t438 * mrSges(7,1) + mrSges(7,2) * t439 - t462 * t453 + t454 * t463;
t421 = m(6) * t430 - t460 * mrSges(6,1) + mrSges(6,2) * t461 + t470 * t490 - t489 * t560 + t558;
t491 = -mrSges(5,2) * t520 - mrSges(5,3) * t506;
t418 = m(5) * t435 + mrSges(5,1) * t502 - mrSges(5,3) * t480 - t484 * t507 + t491 * t520 - t421;
t565 = t579 * t405 - t418 * t547;
t396 = m(4) * t450 - mrSges(4,2) * t535 + mrSges(4,3) * t503 - t499 * t521 - t509 * t536 + t565;
t508 = -mrSges(4,2) * t536 - mrSges(4,3) * t521;
t406 = t409 * t543 + t410 * t540;
t556 = -m(5) * t442 - t479 * mrSges(5,1) - mrSges(5,2) * t480 - t506 * t491 - t492 * t507 - t406;
t402 = m(4) * t449 + mrSges(4,1) * t535 - mrSges(4,3) * t504 - t499 * t522 + t508 * t536 + t556;
t392 = t541 * t396 + t544 * t402;
t398 = t547 * t405 + t579 * t418;
t568 = t551 * t573;
t566 = t544 * t396 - t402 * t541;
t557 = -m(4) * t482 + t503 * mrSges(4,1) - t504 * mrSges(4,2) - t521 * t508 - t522 * t509 - t398;
t445 = Ifges(7,4) * t463 + Ifges(7,2) * t462 + Ifges(7,6) * t505;
t446 = Ifges(7,1) * t463 + Ifges(7,4) * t462 + Ifges(7,5) * t505;
t555 = mrSges(7,1) * t422 - mrSges(7,2) * t423 + Ifges(7,5) * t439 + Ifges(7,6) * t438 + Ifges(7,3) * t478 + t463 * t445 - t462 * t446;
t444 = Ifges(7,5) * t463 + Ifges(7,6) * t462 + Ifges(7,3) * t505;
t412 = -mrSges(7,1) * t428 + mrSges(7,3) * t423 + Ifges(7,4) * t439 + Ifges(7,2) * t438 + Ifges(7,6) * t478 - t444 * t463 + t446 * t505;
t413 = mrSges(7,2) * t428 - mrSges(7,3) * t422 + Ifges(7,1) * t439 + Ifges(7,4) * t438 + Ifges(7,5) * t478 + t444 * t462 - t445 * t505;
t455 = Ifges(6,5) * t490 + Ifges(6,6) * t489 + Ifges(6,3) * t506;
t457 = Ifges(6,1) * t490 + Ifges(6,4) * t489 + Ifges(6,5) * t506;
t399 = -mrSges(6,1) * t430 + mrSges(6,3) * t427 + Ifges(6,4) * t461 + Ifges(6,2) * t460 + Ifges(6,6) * t479 - pkin(5) * t558 + pkin(10) * t564 + t550 * t412 + t546 * t413 - t490 * t455 + t506 * t457;
t456 = Ifges(6,4) * t490 + Ifges(6,2) * t489 + Ifges(6,6) * t506;
t400 = mrSges(6,2) * t430 - mrSges(6,3) * t426 + Ifges(6,1) * t461 + Ifges(6,4) * t460 + Ifges(6,5) * t479 - pkin(10) * t411 - t412 * t546 + t413 * t550 + t455 * t489 - t456 * t506;
t474 = Ifges(5,4) * t507 - Ifges(5,2) * t506 + Ifges(5,6) * t520;
t475 = Ifges(5,1) * t507 - Ifges(5,4) * t506 + Ifges(5,5) * t520;
t554 = mrSges(5,1) * t435 - mrSges(5,2) * t436 + Ifges(5,5) * t480 - Ifges(5,6) * t479 + Ifges(5,3) * t502 - pkin(4) * t421 + qJ(5) * t407 + t543 * t399 + t540 * t400 + t507 * t474 + t506 * t475;
t529 = (-t551 * mrSges(3,1) + t548 * mrSges(3,2)) * t573;
t526 = -mrSges(3,2) * t536 + mrSges(3,3) * t568;
t525 = mrSges(3,1) * t536 - mrSges(3,3) * t569;
t512 = Ifges(3,5) * t536 + (t548 * Ifges(3,1) + t551 * Ifges(3,4)) * t573;
t511 = Ifges(3,6) * t536 + (t548 * Ifges(3,4) + t551 * Ifges(3,2)) * t573;
t510 = Ifges(3,3) * t536 + (t548 * Ifges(3,5) + t551 * Ifges(3,6)) * t573;
t497 = -g(3) * t574 + t562;
t495 = Ifges(4,1) * t522 - Ifges(4,4) * t521 + Ifges(4,5) * t536;
t494 = Ifges(4,4) * t522 - Ifges(4,2) * t521 + Ifges(4,6) * t536;
t493 = Ifges(4,5) * t522 - Ifges(4,6) * t521 + Ifges(4,3) * t536;
t473 = Ifges(5,5) * t507 - Ifges(5,6) * t506 + Ifges(5,3) * t520;
t393 = -pkin(4) * t406 - t555 + t520 * t475 - t507 * t473 + Ifges(5,6) * t502 + t489 * t457 - t490 * t456 + Ifges(5,4) * t480 - Ifges(6,6) * t460 - Ifges(6,5) * t461 - mrSges(5,1) * t442 + mrSges(5,3) * t436 + mrSges(6,2) * t427 - mrSges(6,1) * t426 + (-Ifges(5,2) - Ifges(6,3)) * t479 - pkin(5) * t411;
t391 = m(3) * t498 - mrSges(3,2) * t535 + mrSges(3,3) * t531 - t525 * t536 + t529 * t568 + t566;
t390 = m(3) * t497 + mrSges(3,1) * t535 - mrSges(3,3) * t530 + t526 * t536 - t529 * t569 + t392;
t389 = mrSges(5,2) * t442 - mrSges(5,3) * t435 + Ifges(5,1) * t480 - Ifges(5,4) * t479 + Ifges(5,5) * t502 - qJ(5) * t406 - t399 * t540 + t400 * t543 - t473 * t506 - t474 * t520;
t388 = -mrSges(4,1) * t482 + mrSges(4,3) * t450 + Ifges(4,4) * t504 + Ifges(4,2) * t503 + Ifges(4,6) * t535 - pkin(3) * t398 - t522 * t493 + t536 * t495 - t554;
t387 = mrSges(4,2) * t482 - mrSges(4,3) * t449 + Ifges(4,1) * t504 + Ifges(4,4) * t503 + Ifges(4,5) * t535 - pkin(9) * t398 + t389 * t579 - t547 * t393 - t521 * t493 - t536 * t494;
t386 = Ifges(3,5) * t530 + Ifges(3,6) * t531 + mrSges(3,1) * t497 - mrSges(3,2) * t498 + Ifges(4,5) * t504 + Ifges(4,6) * t503 + t522 * t494 + t521 * t495 + mrSges(4,1) * t449 - mrSges(4,2) * t450 + t547 * t389 + t579 * t393 + pkin(3) * t556 + pkin(9) * t565 + pkin(2) * t392 + (Ifges(3,3) + Ifges(4,3)) * t535 + (t548 * t511 - t551 * t512) * t573;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t567 - mrSges(2,2) * t561 + (mrSges(3,2) * t513 - mrSges(3,3) * t497 + Ifges(3,1) * t530 + Ifges(3,4) * t531 + Ifges(3,5) * t535 - qJ(3) * t392 + t387 * t544 - t388 * t541 + t510 * t568 - t536 * t511) * t575 + (-mrSges(3,1) * t513 + mrSges(3,3) * t498 + Ifges(3,4) * t530 + Ifges(3,2) * t531 + Ifges(3,6) * t535 + pkin(2) * t557 + qJ(3) * t566 + t541 * t387 + t544 * t388 - t510 * t569 + t536 * t512) * t574 + t545 * t386 + pkin(1) * ((t390 * t551 + t391 * t548) * t545 + (-m(3) * t513 + t531 * mrSges(3,1) - t530 * mrSges(3,2) + (-t525 * t548 + t526 * t551) * t573 + t557) * t542) + (-t390 * t548 + t391 * t551) * t578; t386; -t557; t554; t421; t555;];
tauJ  = t1;
