% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-05-06 10:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:01:20
% EndTime: 2019-05-06 10:01:31
% DurationCPUTime: 11.17s
% Computational Cost: add. (149205->359), mult. (403014->471), div. (0->0), fcn. (321551->14), ass. (0->149)
t586 = -2 * qJD(3);
t552 = sin(pkin(6));
t557 = sin(qJ(2));
t561 = cos(qJ(2));
t577 = qJD(1) * qJD(2);
t541 = (qJDD(1) * t557 + t561 * t577) * t552;
t554 = cos(pkin(6));
t545 = qJDD(1) * t554 + qJDD(2);
t546 = qJD(1) * t554 + qJD(2);
t563 = qJD(1) ^ 2;
t558 = sin(qJ(1));
t562 = cos(qJ(1));
t568 = -g(1) * t562 - g(2) * t558;
t585 = t552 * pkin(8);
t539 = -pkin(1) * t563 + qJDD(1) * t585 + t568;
t573 = t558 * g(1) - g(2) * t562;
t538 = qJDD(1) * pkin(1) + t563 * t585 + t573;
t583 = t538 * t554;
t569 = -t539 * t557 + t561 * t583;
t582 = t552 ^ 2 * t563;
t479 = pkin(2) * t545 - qJ(3) * t541 + (pkin(2) * t557 * t582 + (qJ(3) * qJD(1) * t546 - g(3)) * t552) * t561 + t569;
t581 = t552 * t557;
t506 = -g(3) * t581 + t561 * t539 + t557 * t583;
t579 = qJD(1) * t552;
t575 = t557 * t579;
t535 = pkin(2) * t546 - qJ(3) * t575;
t542 = (qJDD(1) * t561 - t557 * t577) * t552;
t576 = t561 ^ 2 * t582;
t482 = -pkin(2) * t576 + qJ(3) * t542 - t535 * t546 + t506;
t551 = sin(pkin(11));
t584 = cos(pkin(11));
t532 = (t551 * t561 + t584 * t557) * t579;
t456 = t584 * t479 - t551 * t482 + t532 * t586;
t580 = t552 * t561;
t574 = t561 * t579;
t531 = t551 * t575 - t584 * t574;
t457 = t551 * t479 + t584 * t482 + t531 * t586;
t508 = mrSges(4,1) * t531 + mrSges(4,2) * t532;
t511 = t541 * t551 - t584 * t542;
t519 = mrSges(4,1) * t546 - mrSges(4,3) * t532;
t507 = pkin(3) * t531 - qJ(4) * t532;
t544 = t546 ^ 2;
t451 = -pkin(3) * t544 + qJ(4) * t545 - t507 * t531 + t457;
t523 = -g(3) * t554 - t538 * t552;
t490 = -pkin(2) * t542 - qJ(3) * t576 + t535 * t575 + qJDD(3) + t523;
t512 = t584 * t541 + t551 * t542;
t460 = (t531 * t546 - t512) * qJ(4) + (t532 * t546 + t511) * pkin(3) + t490;
t550 = sin(pkin(12));
t553 = cos(pkin(12));
t517 = t532 * t553 + t546 * t550;
t443 = -0.2e1 * qJD(4) * t517 - t451 * t550 + t553 * t460;
t500 = t512 * t553 + t545 * t550;
t516 = -t532 * t550 + t546 * t553;
t440 = (t516 * t531 - t500) * pkin(9) + (t516 * t517 + t511) * pkin(4) + t443;
t444 = 0.2e1 * qJD(4) * t516 + t553 * t451 + t550 * t460;
t497 = pkin(4) * t531 - pkin(9) * t517;
t499 = -t512 * t550 + t545 * t553;
t515 = t516 ^ 2;
t442 = -pkin(4) * t515 + pkin(9) * t499 - t497 * t531 + t444;
t556 = sin(qJ(5));
t560 = cos(qJ(5));
t437 = t556 * t440 + t560 * t442;
t491 = t516 * t560 - t517 * t556;
t492 = t516 * t556 + t517 * t560;
t473 = -pkin(5) * t491 - pkin(10) * t492;
t510 = qJDD(5) + t511;
t530 = qJD(5) + t531;
t529 = t530 ^ 2;
t435 = -pkin(5) * t529 + pkin(10) * t510 + t473 * t491 + t437;
t450 = -t545 * pkin(3) - t544 * qJ(4) + t532 * t507 + qJDD(4) - t456;
t445 = -t499 * pkin(4) - t515 * pkin(9) + t517 * t497 + t450;
t464 = -qJD(5) * t492 + t499 * t560 - t500 * t556;
t465 = qJD(5) * t491 + t499 * t556 + t500 * t560;
t438 = (-t491 * t530 - t465) * pkin(10) + (t492 * t530 - t464) * pkin(5) + t445;
t555 = sin(qJ(6));
t559 = cos(qJ(6));
t432 = -t435 * t555 + t438 * t559;
t475 = -t492 * t555 + t530 * t559;
t448 = qJD(6) * t475 + t465 * t559 + t510 * t555;
t476 = t492 * t559 + t530 * t555;
t461 = -mrSges(7,1) * t475 + mrSges(7,2) * t476;
t463 = qJDD(6) - t464;
t489 = qJD(6) - t491;
t466 = -mrSges(7,2) * t489 + mrSges(7,3) * t475;
t429 = m(7) * t432 + mrSges(7,1) * t463 - mrSges(7,3) * t448 - t461 * t476 + t466 * t489;
t433 = t435 * t559 + t438 * t555;
t447 = -qJD(6) * t476 - t465 * t555 + t510 * t559;
t467 = mrSges(7,1) * t489 - mrSges(7,3) * t476;
t430 = m(7) * t433 - mrSges(7,2) * t463 + mrSges(7,3) * t447 + t461 * t475 - t467 * t489;
t421 = -t429 * t555 + t559 * t430;
t472 = -mrSges(6,1) * t491 + mrSges(6,2) * t492;
t481 = mrSges(6,1) * t530 - mrSges(6,3) * t492;
t418 = m(6) * t437 - mrSges(6,2) * t510 + mrSges(6,3) * t464 + t472 * t491 - t481 * t530 + t421;
t436 = t440 * t560 - t442 * t556;
t434 = -pkin(5) * t510 - pkin(10) * t529 + t473 * t492 - t436;
t431 = -m(7) * t434 + t447 * mrSges(7,1) - mrSges(7,2) * t448 + t475 * t466 - t467 * t476;
t480 = -mrSges(6,2) * t530 + mrSges(6,3) * t491;
t425 = m(6) * t436 + mrSges(6,1) * t510 - mrSges(6,3) * t465 - t472 * t492 + t480 * t530 + t431;
t413 = t556 * t418 + t560 * t425;
t493 = -mrSges(5,1) * t516 + mrSges(5,2) * t517;
t495 = -mrSges(5,2) * t531 + mrSges(5,3) * t516;
t411 = m(5) * t443 + mrSges(5,1) * t511 - mrSges(5,3) * t500 - t493 * t517 + t495 * t531 + t413;
t496 = mrSges(5,1) * t531 - mrSges(5,3) * t517;
t570 = t560 * t418 - t425 * t556;
t412 = m(5) * t444 - mrSges(5,2) * t511 + mrSges(5,3) * t499 + t493 * t516 - t496 * t531 + t570;
t571 = -t411 * t550 + t553 * t412;
t403 = m(4) * t457 - mrSges(4,2) * t545 - mrSges(4,3) * t511 - t508 * t531 - t519 * t546 + t571;
t420 = t559 * t429 + t555 * t430;
t566 = m(6) * t445 - t464 * mrSges(6,1) + mrSges(6,2) * t465 - t491 * t480 + t481 * t492 + t420;
t419 = m(5) * t450 - t499 * mrSges(5,1) + mrSges(5,2) * t500 - t516 * t495 + t496 * t517 + t566;
t518 = -mrSges(4,2) * t546 - mrSges(4,3) * t531;
t415 = m(4) * t456 + mrSges(4,1) * t545 - mrSges(4,3) * t512 - t508 * t532 + t518 * t546 - t419;
t400 = t551 * t403 + t584 * t415;
t405 = t553 * t411 + t550 * t412;
t572 = t584 * t403 - t415 * t551;
t404 = m(4) * t490 + t511 * mrSges(4,1) + t512 * mrSges(4,2) + t531 * t518 + t532 * t519 + t405;
t453 = Ifges(7,4) * t476 + Ifges(7,2) * t475 + Ifges(7,6) * t489;
t454 = Ifges(7,1) * t476 + Ifges(7,4) * t475 + Ifges(7,5) * t489;
t565 = mrSges(7,1) * t432 - mrSges(7,2) * t433 + Ifges(7,5) * t448 + Ifges(7,6) * t447 + Ifges(7,3) * t463 + t453 * t476 - t454 * t475;
t452 = Ifges(7,5) * t476 + Ifges(7,6) * t475 + Ifges(7,3) * t489;
t422 = -mrSges(7,1) * t434 + mrSges(7,3) * t433 + Ifges(7,4) * t448 + Ifges(7,2) * t447 + Ifges(7,6) * t463 - t452 * t476 + t454 * t489;
t423 = mrSges(7,2) * t434 - mrSges(7,3) * t432 + Ifges(7,1) * t448 + Ifges(7,4) * t447 + Ifges(7,5) * t463 + t452 * t475 - t453 * t489;
t469 = Ifges(6,4) * t492 + Ifges(6,2) * t491 + Ifges(6,6) * t530;
t470 = Ifges(6,1) * t492 + Ifges(6,4) * t491 + Ifges(6,5) * t530;
t564 = mrSges(6,1) * t436 - mrSges(6,2) * t437 + Ifges(6,5) * t465 + Ifges(6,6) * t464 + Ifges(6,3) * t510 + pkin(5) * t431 + pkin(10) * t421 + t559 * t422 + t555 * t423 + t492 * t469 - t491 * t470;
t540 = (-mrSges(3,1) * t561 + mrSges(3,2) * t557) * t579;
t537 = -mrSges(3,2) * t546 + mrSges(3,3) * t574;
t536 = mrSges(3,1) * t546 - mrSges(3,3) * t575;
t522 = Ifges(3,5) * t546 + (Ifges(3,1) * t557 + Ifges(3,4) * t561) * t579;
t521 = Ifges(3,6) * t546 + (Ifges(3,4) * t557 + Ifges(3,2) * t561) * t579;
t520 = Ifges(3,3) * t546 + (Ifges(3,5) * t557 + Ifges(3,6) * t561) * t579;
t505 = -g(3) * t580 + t569;
t503 = Ifges(4,1) * t532 - Ifges(4,4) * t531 + Ifges(4,5) * t546;
t502 = Ifges(4,4) * t532 - Ifges(4,2) * t531 + Ifges(4,6) * t546;
t501 = Ifges(4,5) * t532 - Ifges(4,6) * t531 + Ifges(4,3) * t546;
t485 = Ifges(5,1) * t517 + Ifges(5,4) * t516 + Ifges(5,5) * t531;
t484 = Ifges(5,4) * t517 + Ifges(5,2) * t516 + Ifges(5,6) * t531;
t483 = Ifges(5,5) * t517 + Ifges(5,6) * t516 + Ifges(5,3) * t531;
t468 = Ifges(6,5) * t492 + Ifges(6,6) * t491 + Ifges(6,3) * t530;
t407 = -mrSges(6,1) * t445 + mrSges(6,3) * t437 + Ifges(6,4) * t465 + Ifges(6,2) * t464 + Ifges(6,6) * t510 - pkin(5) * t420 - t468 * t492 + t470 * t530 - t565;
t406 = mrSges(6,2) * t445 - mrSges(6,3) * t436 + Ifges(6,1) * t465 + Ifges(6,4) * t464 + Ifges(6,5) * t510 - pkin(10) * t420 - t422 * t555 + t423 * t559 + t468 * t491 - t469 * t530;
t399 = m(3) * t506 - mrSges(3,2) * t545 + mrSges(3,3) * t542 - t536 * t546 + t540 * t574 + t572;
t398 = m(3) * t505 + mrSges(3,1) * t545 - mrSges(3,3) * t541 + t537 * t546 - t540 * t575 + t400;
t397 = mrSges(5,2) * t450 - mrSges(5,3) * t443 + Ifges(5,1) * t500 + Ifges(5,4) * t499 + Ifges(5,5) * t511 - pkin(9) * t413 + t406 * t560 - t407 * t556 + t483 * t516 - t484 * t531;
t396 = -mrSges(5,1) * t450 + mrSges(5,3) * t444 + Ifges(5,4) * t500 + Ifges(5,2) * t499 + Ifges(5,6) * t511 - pkin(4) * t566 + pkin(9) * t570 + t556 * t406 + t560 * t407 - t517 * t483 + t531 * t485;
t395 = (-Ifges(5,3) - Ifges(4,2)) * t511 - pkin(4) * t413 - t564 - pkin(3) * t405 - mrSges(5,1) * t443 + mrSges(5,2) * t444 + mrSges(4,3) * t457 - mrSges(4,1) * t490 - Ifges(5,6) * t499 - Ifges(5,5) * t500 + Ifges(4,4) * t512 + t516 * t485 - t517 * t484 - t532 * t501 + Ifges(4,6) * t545 + t546 * t503;
t394 = mrSges(4,2) * t490 - mrSges(4,3) * t456 + Ifges(4,1) * t512 - Ifges(4,4) * t511 + Ifges(4,5) * t545 - qJ(4) * t405 - t396 * t550 + t397 * t553 - t501 * t531 - t502 * t546;
t393 = Ifges(3,5) * t541 + Ifges(3,6) * t542 + mrSges(3,1) * t505 - mrSges(3,2) * t506 + Ifges(4,5) * t512 - Ifges(4,6) * t511 + t532 * t502 + t531 * t503 + mrSges(4,1) * t456 - mrSges(4,2) * t457 + t550 * t397 + t553 * t396 - pkin(3) * t419 + qJ(4) * t571 + pkin(2) * t400 + (Ifges(3,3) + Ifges(4,3)) * t545 + (t521 * t557 - t522 * t561) * t579;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t573 - mrSges(2,2) * t568 + (mrSges(3,2) * t523 - mrSges(3,3) * t505 + Ifges(3,1) * t541 + Ifges(3,4) * t542 + Ifges(3,5) * t545 - qJ(3) * t400 + t394 * t584 - t395 * t551 + t520 * t574 - t521 * t546) * t581 + (-mrSges(3,1) * t523 + mrSges(3,3) * t506 + Ifges(3,4) * t541 + Ifges(3,2) * t542 + Ifges(3,6) * t545 - pkin(2) * t404 + qJ(3) * t572 + t551 * t394 + t395 * t584 - t520 * t575 + t546 * t522) * t580 + t554 * t393 + pkin(1) * ((t398 * t561 + t399 * t557) * t554 + (-m(3) * t523 + t542 * mrSges(3,1) - t541 * mrSges(3,2) + (-t536 * t557 + t537 * t561) * t579 - t404) * t552) + (-t398 * t557 + t399 * t561) * t585; t393; t404; t419; t564; t565;];
tauJ  = t1;
