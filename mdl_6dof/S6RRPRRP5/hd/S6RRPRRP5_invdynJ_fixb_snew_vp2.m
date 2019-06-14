% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRP5
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
% Datum: 2019-05-06 17:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:46:50
% EndTime: 2019-05-06 17:46:59
% DurationCPUTime: 6.42s
% Computational Cost: add. (78723->339), mult. (205946->429), div. (0->0), fcn. (162206->12), ass. (0->143)
t583 = -2 * qJD(3);
t582 = Ifges(6,1) + Ifges(7,1);
t576 = Ifges(6,4) + Ifges(7,4);
t575 = Ifges(6,5) + Ifges(7,5);
t581 = Ifges(6,2) + Ifges(7,2);
t574 = Ifges(6,6) + Ifges(7,6);
t580 = Ifges(6,3) + Ifges(7,3);
t532 = sin(pkin(11));
t534 = cos(pkin(11));
t538 = sin(qJ(2));
t542 = cos(qJ(2));
t533 = sin(pkin(6));
t564 = qJD(1) * t533;
t515 = (t538 * t532 - t542 * t534) * t564;
t562 = qJD(1) * qJD(2);
t524 = (qJDD(1) * t538 + t542 * t562) * t533;
t535 = cos(pkin(6));
t527 = qJDD(1) * t535 + qJDD(2);
t528 = qJD(1) * t535 + qJD(2);
t544 = qJD(1) ^ 2;
t539 = sin(qJ(1));
t543 = cos(qJ(1));
t550 = -g(1) * t543 - g(2) * t539;
t578 = t533 * pkin(8);
t522 = -pkin(1) * t544 + qJDD(1) * t578 + t550;
t556 = t539 * g(1) - g(2) * t543;
t521 = qJDD(1) * pkin(1) + t544 * t578 + t556;
t572 = t521 * t535;
t552 = -t522 * t538 + t542 * t572;
t571 = t533 ^ 2 * t544;
t461 = pkin(2) * t527 - qJ(3) * t524 + (pkin(2) * t538 * t571 + (qJ(3) * qJD(1) * t528 - g(3)) * t533) * t542 + t552;
t570 = t533 * t538;
t492 = -g(3) * t570 + t542 * t522 + t538 * t572;
t558 = t538 * t564;
t518 = pkin(2) * t528 - qJ(3) * t558;
t525 = (qJDD(1) * t542 - t538 * t562) * t533;
t561 = t542 ^ 2 * t571;
t467 = -pkin(2) * t561 + qJ(3) * t525 - t518 * t528 + t492;
t516 = (t542 * t532 + t538 * t534) * t564;
t436 = t461 * t534 - t532 * t467 + t516 * t583;
t498 = t524 * t534 + t525 * t532;
t537 = sin(qJ(4));
t541 = cos(qJ(4));
t500 = -t516 * t537 + t528 * t541;
t476 = qJD(4) * t500 + t498 * t541 + t527 * t537;
t501 = t516 * t541 + t528 * t537;
t514 = qJD(4) + t515;
t536 = sin(qJ(5));
t540 = cos(qJ(5));
t483 = -t501 * t536 + t514 * t540;
t497 = -t524 * t532 + t525 * t534;
t496 = qJDD(4) - t497;
t444 = qJD(5) * t483 + t476 * t540 + t496 * t536;
t484 = t501 * t540 + t514 * t536;
t455 = -mrSges(7,1) * t483 + mrSges(7,2) * t484;
t437 = t532 * t461 + t534 * t467 + t515 * t583;
t494 = pkin(3) * t515 - pkin(9) * t516;
t526 = t528 ^ 2;
t435 = -pkin(3) * t526 + pkin(9) * t527 - t494 * t515 + t437;
t507 = -g(3) * t535 - t521 * t533;
t478 = -pkin(2) * t525 - qJ(3) * t561 + t518 * t558 + qJDD(3) + t507;
t439 = (t515 * t528 - t498) * pkin(9) + (t516 * t528 - t497) * pkin(3) + t478;
t431 = t541 * t435 + t537 * t439;
t480 = -pkin(4) * t500 - pkin(10) * t501;
t513 = t514 ^ 2;
t426 = -pkin(4) * t513 + pkin(10) * t496 + t480 * t500 + t431;
t434 = -pkin(3) * t527 - pkin(9) * t526 + t516 * t494 - t436;
t475 = -qJD(4) * t501 - t498 * t537 + t527 * t541;
t429 = (-t500 * t514 - t476) * pkin(10) + (t501 * t514 - t475) * pkin(4) + t434;
t421 = -t426 * t536 + t540 * t429;
t474 = qJDD(5) - t475;
t499 = qJD(5) - t500;
t417 = -0.2e1 * qJD(6) * t484 + (t483 * t499 - t444) * qJ(6) + (t483 * t484 + t474) * pkin(5) + t421;
t462 = -mrSges(7,2) * t499 + mrSges(7,3) * t483;
t560 = m(7) * t417 + t474 * mrSges(7,1) + t499 * t462;
t415 = -mrSges(7,3) * t444 - t455 * t484 + t560;
t422 = t540 * t426 + t536 * t429;
t443 = -qJD(5) * t484 - t476 * t536 + t496 * t540;
t464 = pkin(5) * t499 - qJ(6) * t484;
t482 = t483 ^ 2;
t420 = -pkin(5) * t482 + qJ(6) * t443 + 0.2e1 * qJD(6) * t483 - t464 * t499 + t422;
t566 = t483 * t576 + t484 * t582 + t499 * t575;
t567 = -t483 * t581 - t484 * t576 - t499 * t574;
t579 = mrSges(6,1) * t421 + mrSges(7,1) * t417 - mrSges(6,2) * t422 - mrSges(7,2) * t420 + pkin(5) * t415 + t443 * t574 + t444 * t575 + t474 * t580 - t483 * t566 - t484 * t567;
t577 = -mrSges(6,2) - mrSges(7,2);
t569 = t533 * t542;
t493 = mrSges(4,1) * t515 + mrSges(4,2) * t516;
t503 = mrSges(4,1) * t528 - mrSges(4,3) * t516;
t456 = -mrSges(6,1) * t483 + mrSges(6,2) * t484;
t463 = -mrSges(6,2) * t499 + mrSges(6,3) * t483;
t409 = m(6) * t421 + mrSges(6,1) * t474 + t463 * t499 + (-t455 - t456) * t484 + (-mrSges(6,3) - mrSges(7,3)) * t444 + t560;
t559 = m(7) * t420 + t443 * mrSges(7,3) + t483 * t455;
t465 = mrSges(7,1) * t499 - mrSges(7,3) * t484;
t565 = -mrSges(6,1) * t499 + mrSges(6,3) * t484 - t465;
t411 = m(6) * t422 + mrSges(6,3) * t443 + t456 * t483 + t474 * t577 + t499 * t565 + t559;
t408 = -t409 * t536 + t540 * t411;
t479 = -mrSges(5,1) * t500 + mrSges(5,2) * t501;
t486 = mrSges(5,1) * t514 - mrSges(5,3) * t501;
t405 = m(5) * t431 - mrSges(5,2) * t496 + mrSges(5,3) * t475 + t479 * t500 - t486 * t514 + t408;
t430 = -t537 * t435 + t439 * t541;
t425 = -pkin(4) * t496 - pkin(10) * t513 + t501 * t480 - t430;
t423 = -pkin(5) * t443 - qJ(6) * t482 + t464 * t484 + qJDD(6) + t425;
t551 = -m(7) * t423 + t443 * mrSges(7,1) + t483 * t462;
t414 = -m(6) * t425 + t443 * mrSges(6,1) + t444 * t577 + t483 * t463 + t484 * t565 + t551;
t485 = -mrSges(5,2) * t514 + mrSges(5,3) * t500;
t413 = m(5) * t430 + mrSges(5,1) * t496 - mrSges(5,3) * t476 - t479 * t501 + t485 * t514 + t414;
t554 = t541 * t405 - t413 * t537;
t397 = m(4) * t437 - mrSges(4,2) * t527 + mrSges(4,3) * t497 - t493 * t515 - t503 * t528 + t554;
t502 = -mrSges(4,2) * t528 - mrSges(4,3) * t515;
t407 = t409 * t540 + t411 * t536;
t546 = -m(5) * t434 + t475 * mrSges(5,1) - mrSges(5,2) * t476 + t500 * t485 - t486 * t501 - t407;
t402 = m(4) * t436 + mrSges(4,1) * t527 - mrSges(4,3) * t498 - t493 * t516 + t502 * t528 + t546;
t393 = t532 * t397 + t534 * t402;
t399 = t537 * t405 + t541 * t413;
t568 = -t483 * t574 - t484 * t575 - t499 * t580;
t557 = t542 * t564;
t555 = t534 * t397 - t402 * t532;
t548 = -m(4) * t478 + t497 * mrSges(4,1) - t498 * mrSges(4,2) - t515 * t502 - t516 * t503 - t399;
t418 = mrSges(7,2) * t444 + t465 * t484 - t551;
t400 = -mrSges(6,1) * t425 + mrSges(6,3) * t422 - mrSges(7,1) * t423 + mrSges(7,3) * t420 - pkin(5) * t418 + qJ(6) * t559 + (-qJ(6) * t465 + t566) * t499 + t568 * t484 + (-qJ(6) * mrSges(7,2) + t574) * t474 + t576 * t444 + t581 * t443;
t406 = mrSges(6,2) * t425 + mrSges(7,2) * t423 - mrSges(6,3) * t421 - mrSges(7,3) * t417 - qJ(6) * t415 + t576 * t443 + t444 * t582 + t575 * t474 - t568 * t483 + t567 * t499;
t469 = Ifges(5,4) * t501 + Ifges(5,2) * t500 + Ifges(5,6) * t514;
t470 = Ifges(5,1) * t501 + Ifges(5,4) * t500 + Ifges(5,5) * t514;
t545 = mrSges(5,1) * t430 - mrSges(5,2) * t431 + Ifges(5,5) * t476 + Ifges(5,6) * t475 + Ifges(5,3) * t496 + pkin(4) * t414 + pkin(10) * t408 + t540 * t400 + t536 * t406 + t501 * t469 - t500 * t470;
t523 = (-t542 * mrSges(3,1) + t538 * mrSges(3,2)) * t564;
t520 = -mrSges(3,2) * t528 + mrSges(3,3) * t557;
t519 = mrSges(3,1) * t528 - mrSges(3,3) * t558;
t506 = Ifges(3,5) * t528 + (t538 * Ifges(3,1) + t542 * Ifges(3,4)) * t564;
t505 = Ifges(3,6) * t528 + (t538 * Ifges(3,4) + t542 * Ifges(3,2)) * t564;
t504 = Ifges(3,3) * t528 + (t538 * Ifges(3,5) + t542 * Ifges(3,6)) * t564;
t491 = -g(3) * t569 + t552;
t489 = Ifges(4,1) * t516 - Ifges(4,4) * t515 + Ifges(4,5) * t528;
t488 = Ifges(4,4) * t516 - Ifges(4,2) * t515 + Ifges(4,6) * t528;
t487 = Ifges(4,5) * t516 - Ifges(4,6) * t515 + Ifges(4,3) * t528;
t468 = Ifges(5,5) * t501 + Ifges(5,6) * t500 + Ifges(5,3) * t514;
t394 = -mrSges(5,1) * t434 + mrSges(5,3) * t431 + Ifges(5,4) * t476 + Ifges(5,2) * t475 + Ifges(5,6) * t496 - pkin(4) * t407 - t501 * t468 + t514 * t470 - t579;
t392 = m(3) * t492 - mrSges(3,2) * t527 + mrSges(3,3) * t525 - t519 * t528 + t523 * t557 + t555;
t391 = m(3) * t491 + mrSges(3,1) * t527 - mrSges(3,3) * t524 + t520 * t528 - t523 * t558 + t393;
t390 = mrSges(5,2) * t434 - mrSges(5,3) * t430 + Ifges(5,1) * t476 + Ifges(5,4) * t475 + Ifges(5,5) * t496 - pkin(10) * t407 - t400 * t536 + t406 * t540 + t468 * t500 - t469 * t514;
t389 = -mrSges(4,1) * t478 + mrSges(4,3) * t437 + Ifges(4,4) * t498 + Ifges(4,2) * t497 + Ifges(4,6) * t527 - pkin(3) * t399 - t516 * t487 + t528 * t489 - t545;
t388 = mrSges(4,2) * t478 - mrSges(4,3) * t436 + Ifges(4,1) * t498 + Ifges(4,4) * t497 + Ifges(4,5) * t527 - pkin(9) * t399 + t390 * t541 - t394 * t537 - t487 * t515 - t488 * t528;
t387 = Ifges(3,5) * t524 + Ifges(3,6) * t525 + mrSges(3,1) * t491 - mrSges(3,2) * t492 + Ifges(4,5) * t498 + Ifges(4,6) * t497 + t516 * t488 + t515 * t489 + mrSges(4,1) * t436 - mrSges(4,2) * t437 + t537 * t390 + t541 * t394 + pkin(3) * t546 + pkin(9) * t554 + pkin(2) * t393 + (Ifges(3,3) + Ifges(4,3)) * t527 + (t538 * t505 - t542 * t506) * t564;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t556 - mrSges(2,2) * t550 + (mrSges(3,2) * t507 - mrSges(3,3) * t491 + Ifges(3,1) * t524 + Ifges(3,4) * t525 + Ifges(3,5) * t527 - qJ(3) * t393 + t388 * t534 - t389 * t532 + t504 * t557 - t505 * t528) * t570 + (-mrSges(3,1) * t507 + mrSges(3,3) * t492 + Ifges(3,4) * t524 + Ifges(3,2) * t525 + Ifges(3,6) * t527 + pkin(2) * t548 + qJ(3) * t555 + t532 * t388 + t534 * t389 - t504 * t558 + t528 * t506) * t569 + t535 * t387 + pkin(1) * ((t391 * t542 + t392 * t538) * t535 + (-m(3) * t507 + t525 * mrSges(3,1) - t524 * mrSges(3,2) + (-t519 * t538 + t520 * t542) * t564 + t548) * t533) + (-t538 * t391 + t542 * t392) * t578; t387; -t548; t545; t579; t418;];
tauJ  = t1;
