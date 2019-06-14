% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:59:42
% EndTime: 2019-05-06 17:59:48
% DurationCPUTime: 6.23s
% Computational Cost: add. (76758->337), mult. (200464->429), div. (0->0), fcn. (157557->12), ass. (0->142)
t581 = -2 * qJD(3);
t580 = Ifges(6,1) + Ifges(7,1);
t573 = Ifges(6,4) - Ifges(7,5);
t572 = -Ifges(6,5) - Ifges(7,4);
t579 = Ifges(6,2) + Ifges(7,3);
t571 = Ifges(6,6) - Ifges(7,6);
t578 = -Ifges(6,3) - Ifges(7,2);
t531 = sin(pkin(11));
t533 = cos(pkin(11));
t537 = sin(qJ(2));
t540 = cos(qJ(2));
t532 = sin(pkin(6));
t561 = qJD(1) * t532;
t514 = (t531 * t537 - t533 * t540) * t561;
t559 = qJD(1) * qJD(2);
t523 = (qJDD(1) * t537 + t540 * t559) * t532;
t534 = cos(pkin(6));
t526 = qJDD(1) * t534 + qJDD(2);
t527 = qJD(1) * t534 + qJD(2);
t542 = qJD(1) ^ 2;
t538 = sin(qJ(1));
t541 = cos(qJ(1));
t548 = -g(1) * t541 - g(2) * t538;
t575 = pkin(8) * t532;
t521 = -pkin(1) * t542 + qJDD(1) * t575 + t548;
t554 = t538 * g(1) - g(2) * t541;
t520 = qJDD(1) * pkin(1) + t542 * t575 + t554;
t569 = t520 * t534;
t550 = -t537 * t521 + t540 * t569;
t568 = t532 ^ 2 * t542;
t458 = t526 * pkin(2) - t523 * qJ(3) + (pkin(2) * t537 * t568 + (qJ(3) * qJD(1) * t527 - g(3)) * t532) * t540 + t550;
t567 = t532 * t537;
t488 = -g(3) * t567 + t540 * t521 + t537 * t569;
t556 = t537 * t561;
t517 = pkin(2) * t527 - qJ(3) * t556;
t524 = (qJDD(1) * t540 - t537 * t559) * t532;
t558 = t540 ^ 2 * t568;
t463 = -pkin(2) * t558 + qJ(3) * t524 - t517 * t527 + t488;
t515 = (t531 * t540 + t533 * t537) * t561;
t433 = t458 * t533 - t531 * t463 + t515 * t581;
t495 = t523 * t533 + t524 * t531;
t536 = sin(qJ(4));
t539 = cos(qJ(4));
t498 = -t515 * t536 + t527 * t539;
t473 = qJD(4) * t498 + t495 * t539 + t526 * t536;
t499 = t515 * t539 + t527 * t536;
t513 = qJD(4) + t514;
t535 = sin(qJ(5));
t576 = cos(qJ(5));
t479 = t499 * t535 - t513 * t576;
t494 = -t523 * t531 + t524 * t533;
t493 = qJDD(4) - t494;
t440 = -t479 * qJD(5) + t473 * t576 + t535 * t493;
t480 = t499 * t576 + t535 * t513;
t452 = mrSges(7,1) * t479 - mrSges(7,3) * t480;
t434 = t531 * t458 + t533 * t463 + t514 * t581;
t491 = pkin(3) * t514 - pkin(9) * t515;
t525 = t527 ^ 2;
t432 = -pkin(3) * t525 + pkin(9) * t526 - t491 * t514 + t434;
t505 = -g(3) * t534 - t520 * t532;
t475 = -pkin(2) * t524 - qJ(3) * t558 + t517 * t556 + qJDD(3) + t505;
t436 = (t514 * t527 - t495) * pkin(9) + (t515 * t527 - t494) * pkin(3) + t475;
t428 = t539 * t432 + t536 * t436;
t477 = -pkin(4) * t498 - pkin(10) * t499;
t512 = t513 ^ 2;
t424 = -pkin(4) * t512 + pkin(10) * t493 + t477 * t498 + t428;
t431 = -pkin(3) * t526 - pkin(9) * t525 + t515 * t491 - t433;
t472 = -qJD(4) * t499 - t495 * t536 + t526 * t539;
t426 = (-t498 * t513 - t473) * pkin(10) + (t499 * t513 - t472) * pkin(4) + t431;
t420 = -t535 * t424 + t426 * t576;
t451 = pkin(5) * t479 - qJ(6) * t480;
t471 = qJDD(5) - t472;
t497 = qJD(5) - t498;
t496 = t497 ^ 2;
t418 = -t471 * pkin(5) - t496 * qJ(6) + t480 * t451 + qJDD(6) - t420;
t459 = -mrSges(7,2) * t479 + mrSges(7,3) * t497;
t549 = -m(7) * t418 + t471 * mrSges(7,1) + t497 * t459;
t414 = t440 * mrSges(7,2) + t480 * t452 - t549;
t421 = t424 * t576 + t535 * t426;
t417 = -pkin(5) * t496 + qJ(6) * t471 + 0.2e1 * qJD(6) * t497 - t451 * t479 + t421;
t439 = qJD(5) * t480 + t473 * t535 - t493 * t576;
t462 = -mrSges(7,1) * t497 + mrSges(7,2) * t480;
t557 = m(7) * t417 + t471 * mrSges(7,3) + t497 * t462;
t563 = t479 * t573 - t480 * t580 + t497 * t572;
t564 = t479 * t579 - t480 * t573 - t497 * t571;
t577 = -t439 * t571 - t440 * t572 - t578 * t471 - t479 * t563 - t480 * t564 + mrSges(6,1) * t420 - mrSges(7,1) * t418 - mrSges(6,2) * t421 + mrSges(7,3) * t417 - pkin(5) * t414 + qJ(6) * (-t439 * mrSges(7,2) - t479 * t452 + t557);
t574 = -mrSges(6,3) - mrSges(7,2);
t566 = t532 * t540;
t489 = mrSges(4,1) * t514 + mrSges(4,2) * t515;
t501 = mrSges(4,1) * t527 - mrSges(4,3) * t515;
t461 = mrSges(6,1) * t497 - mrSges(6,3) * t480;
t562 = -mrSges(6,1) * t479 - mrSges(6,2) * t480 - t452;
t410 = m(6) * t421 - t471 * mrSges(6,2) + t439 * t574 - t497 * t461 + t479 * t562 + t557;
t460 = -mrSges(6,2) * t497 - mrSges(6,3) * t479;
t411 = m(6) * t420 + t471 * mrSges(6,1) + t440 * t574 + t497 * t460 + t480 * t562 + t549;
t406 = t410 * t576 - t411 * t535;
t476 = -mrSges(5,1) * t498 + mrSges(5,2) * t499;
t482 = mrSges(5,1) * t513 - mrSges(5,3) * t499;
t402 = m(5) * t428 - mrSges(5,2) * t493 + mrSges(5,3) * t472 + t476 * t498 - t482 * t513 + t406;
t427 = -t536 * t432 + t436 * t539;
t423 = -pkin(4) * t493 - pkin(10) * t512 + t499 * t477 - t427;
t419 = -0.2e1 * qJD(6) * t480 + (t479 * t497 - t440) * qJ(6) + (t480 * t497 + t439) * pkin(5) + t423;
t415 = m(7) * t419 + mrSges(7,1) * t439 - t440 * mrSges(7,3) + t459 * t479 - t480 * t462;
t412 = -m(6) * t423 - t439 * mrSges(6,1) - mrSges(6,2) * t440 - t479 * t460 - t461 * t480 - t415;
t481 = -mrSges(5,2) * t513 + mrSges(5,3) * t498;
t408 = m(5) * t427 + mrSges(5,1) * t493 - mrSges(5,3) * t473 - t476 * t499 + t481 * t513 + t412;
t552 = t539 * t402 - t408 * t536;
t395 = m(4) * t434 - mrSges(4,2) * t526 + mrSges(4,3) * t494 - t489 * t514 - t501 * t527 + t552;
t500 = -mrSges(4,2) * t527 - mrSges(4,3) * t514;
t405 = t535 * t410 + t411 * t576;
t544 = -m(5) * t431 + t472 * mrSges(5,1) - t473 * mrSges(5,2) + t498 * t481 - t499 * t482 - t405;
t399 = m(4) * t433 + t526 * mrSges(4,1) - t495 * mrSges(4,3) - t515 * t489 + t527 * t500 + t544;
t391 = t531 * t395 + t533 * t399;
t397 = t536 * t402 + t539 * t408;
t565 = t479 * t571 + t480 * t572 + t497 * t578;
t555 = t540 * t561;
t553 = t533 * t395 - t399 * t531;
t546 = -m(4) * t475 + t494 * mrSges(4,1) - t495 * mrSges(4,2) - t514 * t500 - t515 * t501 - t397;
t403 = -mrSges(6,1) * t423 - mrSges(7,1) * t419 + mrSges(7,2) * t417 + mrSges(6,3) * t421 - pkin(5) * t415 - t439 * t579 + t573 * t440 + t571 * t471 + t565 * t480 - t563 * t497;
t404 = mrSges(6,2) * t423 + mrSges(7,2) * t418 - mrSges(6,3) * t420 - mrSges(7,3) * t419 - qJ(6) * t415 - t573 * t439 + t440 * t580 - t572 * t471 + t565 * t479 + t564 * t497;
t465 = Ifges(5,4) * t499 + Ifges(5,2) * t498 + Ifges(5,6) * t513;
t466 = Ifges(5,1) * t499 + Ifges(5,4) * t498 + Ifges(5,5) * t513;
t543 = mrSges(5,1) * t427 - mrSges(5,2) * t428 + Ifges(5,5) * t473 + Ifges(5,6) * t472 + Ifges(5,3) * t493 + pkin(4) * t412 + pkin(10) * t406 + t403 * t576 + t535 * t404 + t499 * t465 - t498 * t466;
t522 = (-mrSges(3,1) * t540 + mrSges(3,2) * t537) * t561;
t519 = -mrSges(3,2) * t527 + mrSges(3,3) * t555;
t518 = mrSges(3,1) * t527 - mrSges(3,3) * t556;
t504 = Ifges(3,5) * t527 + (Ifges(3,1) * t537 + Ifges(3,4) * t540) * t561;
t503 = Ifges(3,6) * t527 + (Ifges(3,4) * t537 + Ifges(3,2) * t540) * t561;
t502 = Ifges(3,3) * t527 + (Ifges(3,5) * t537 + Ifges(3,6) * t540) * t561;
t487 = -g(3) * t566 + t550;
t485 = Ifges(4,1) * t515 - Ifges(4,4) * t514 + Ifges(4,5) * t527;
t484 = Ifges(4,4) * t515 - Ifges(4,2) * t514 + Ifges(4,6) * t527;
t483 = Ifges(4,5) * t515 - Ifges(4,6) * t514 + Ifges(4,3) * t527;
t464 = Ifges(5,5) * t499 + Ifges(5,6) * t498 + Ifges(5,3) * t513;
t392 = -mrSges(5,1) * t431 + mrSges(5,3) * t428 + Ifges(5,4) * t473 + Ifges(5,2) * t472 + Ifges(5,6) * t493 - pkin(4) * t405 - t499 * t464 + t513 * t466 - t577;
t390 = m(3) * t488 - mrSges(3,2) * t526 + mrSges(3,3) * t524 - t518 * t527 + t522 * t555 + t553;
t389 = m(3) * t487 + mrSges(3,1) * t526 - mrSges(3,3) * t523 + t519 * t527 - t522 * t556 + t391;
t388 = mrSges(5,2) * t431 - mrSges(5,3) * t427 + Ifges(5,1) * t473 + Ifges(5,4) * t472 + Ifges(5,5) * t493 - pkin(10) * t405 - t535 * t403 + t404 * t576 + t498 * t464 - t513 * t465;
t387 = -mrSges(4,1) * t475 + mrSges(4,3) * t434 + Ifges(4,4) * t495 + Ifges(4,2) * t494 + Ifges(4,6) * t526 - pkin(3) * t397 - t515 * t483 + t527 * t485 - t543;
t386 = mrSges(4,2) * t475 - mrSges(4,3) * t433 + Ifges(4,1) * t495 + Ifges(4,4) * t494 + Ifges(4,5) * t526 - pkin(9) * t397 + t388 * t539 - t392 * t536 - t483 * t514 - t484 * t527;
t385 = Ifges(3,5) * t523 + Ifges(3,6) * t524 + mrSges(3,1) * t487 - mrSges(3,2) * t488 + Ifges(4,5) * t495 + Ifges(4,6) * t494 + t515 * t484 + t514 * t485 + mrSges(4,1) * t433 - mrSges(4,2) * t434 + t536 * t388 + t539 * t392 + pkin(3) * t544 + pkin(9) * t552 + pkin(2) * t391 + (Ifges(3,3) + Ifges(4,3)) * t526 + (t503 * t537 - t504 * t540) * t561;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t554 - mrSges(2,2) * t548 + (mrSges(3,2) * t505 - mrSges(3,3) * t487 + Ifges(3,1) * t523 + Ifges(3,4) * t524 + Ifges(3,5) * t526 - qJ(3) * t391 + t386 * t533 - t387 * t531 + t502 * t555 - t503 * t527) * t567 + (-mrSges(3,1) * t505 + mrSges(3,3) * t488 + Ifges(3,4) * t523 + Ifges(3,2) * t524 + Ifges(3,6) * t526 + pkin(2) * t546 + qJ(3) * t553 + t531 * t386 + t533 * t387 - t502 * t556 + t527 * t504) * t566 + t534 * t385 + pkin(1) * ((t389 * t540 + t390 * t537) * t534 + (-m(3) * t505 + t524 * mrSges(3,1) - t523 * mrSges(3,2) + (-t518 * t537 + t519 * t540) * t561 + t546) * t532) + (-t389 * t537 + t390 * t540) * t575; t385; -t546; t543; t577; t414;];
tauJ  = t1;
