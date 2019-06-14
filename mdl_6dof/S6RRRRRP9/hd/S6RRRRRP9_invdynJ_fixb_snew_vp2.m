% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 06:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRP9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:59:14
% EndTime: 2019-05-08 05:59:37
% DurationCPUTime: 9.67s
% Computational Cost: add. (135206->340), mult. (287561->429), div. (0->0), fcn. (230760->12), ass. (0->144)
t585 = Ifges(6,4) + Ifges(7,4);
t595 = Ifges(6,2) + Ifges(7,2);
t591 = Ifges(6,6) + Ifges(7,6);
t549 = cos(pkin(6));
t545 = qJD(1) * t549 + qJD(2);
t552 = sin(qJ(3));
t557 = cos(qJ(3));
t553 = sin(qJ(2));
t548 = sin(pkin(6));
t578 = qJD(1) * t548;
t573 = t553 * t578;
t526 = t545 * t557 - t552 * t573;
t558 = cos(qJ(2));
t576 = qJD(1) * qJD(2);
t538 = (qJDD(1) * t553 + t558 * t576) * t548;
t544 = qJDD(1) * t549 + qJDD(2);
t507 = qJD(3) * t526 + t538 * t557 + t544 * t552;
t527 = t545 * t552 + t557 * t573;
t577 = qJD(1) * t558;
t572 = t548 * t577;
t542 = qJD(3) - t572;
t551 = sin(qJ(4));
t556 = cos(qJ(4));
t514 = t527 * t556 + t542 * t551;
t539 = (-qJDD(1) * t558 + t553 * t576) * t548;
t531 = qJDD(3) + t539;
t470 = -qJD(4) * t514 - t507 * t551 + t531 * t556;
t513 = -t527 * t551 + t542 * t556;
t471 = qJD(4) * t513 + t507 * t556 + t531 * t551;
t550 = sin(qJ(5));
t555 = cos(qJ(5));
t489 = t513 * t555 - t514 * t550;
t445 = qJD(5) * t489 + t470 * t550 + t471 * t555;
t490 = t513 * t550 + t514 * t555;
t466 = -mrSges(7,1) * t489 + mrSges(7,2) * t490;
t537 = (-pkin(2) * t558 - pkin(9) * t553) * t578;
t543 = t545 ^ 2;
t560 = qJD(1) ^ 2;
t554 = sin(qJ(1));
t559 = cos(qJ(1));
t567 = -g(1) * t559 - g(2) * t554;
t587 = pkin(8) * t548;
t535 = -pkin(1) * t560 + qJDD(1) * t587 + t567;
t571 = t554 * g(1) - g(2) * t559;
t534 = qJDD(1) * pkin(1) + t560 * t587 + t571;
t584 = t534 * t549;
t579 = t558 * t535 + t553 * t584;
t486 = -t543 * pkin(2) + t544 * pkin(9) + (-g(3) * t553 + t537 * t577) * t548 + t579;
t586 = t549 * g(3);
t487 = t539 * pkin(2) - t538 * pkin(9) - t586 + (-t534 + (pkin(2) * t553 - pkin(9) * t558) * t545 * qJD(1)) * t548;
t456 = t557 * t486 + t552 * t487;
t511 = -pkin(3) * t526 - pkin(10) * t527;
t540 = t542 ^ 2;
t448 = -pkin(3) * t540 + pkin(10) * t531 + t511 * t526 + t456;
t582 = t548 * t558;
t508 = -g(3) * t582 - t553 * t535 + t558 * t584;
t485 = -t544 * pkin(2) - t543 * pkin(9) + t537 * t573 - t508;
t506 = -t527 * qJD(3) - t552 * t538 + t544 * t557;
t453 = (-t526 * t542 - t507) * pkin(10) + (t527 * t542 - t506) * pkin(3) + t485;
t433 = -t551 * t448 + t556 * t453;
t504 = qJDD(4) - t506;
t525 = qJD(4) - t526;
t430 = (t513 * t525 - t471) * pkin(11) + (t513 * t514 + t504) * pkin(4) + t433;
t434 = t556 * t448 + t551 * t453;
t495 = pkin(4) * t525 - pkin(11) * t514;
t512 = t513 ^ 2;
t432 = -pkin(4) * t512 + pkin(11) * t470 - t495 * t525 + t434;
t424 = t555 * t430 - t432 * t550;
t499 = qJDD(5) + t504;
t523 = qJD(5) + t525;
t419 = -0.2e1 * qJD(6) * t490 + (t489 * t523 - t445) * qJ(6) + (t489 * t490 + t499) * pkin(5) + t424;
t472 = -mrSges(7,2) * t523 + mrSges(7,3) * t489;
t575 = m(7) * t419 + t499 * mrSges(7,1) + t523 * t472;
t416 = -t445 * mrSges(7,3) - t466 * t490 + t575;
t425 = t550 * t430 + t555 * t432;
t444 = -qJD(5) * t490 + t470 * t555 - t471 * t550;
t474 = pkin(5) * t523 - qJ(6) * t490;
t488 = t489 ^ 2;
t421 = -pkin(5) * t488 + qJ(6) * t444 + 0.2e1 * qJD(6) * t489 - t474 * t523 + t425;
t592 = Ifges(6,5) + Ifges(7,5);
t593 = Ifges(6,1) + Ifges(7,1);
t580 = -t585 * t489 - t490 * t593 - t592 * t523;
t589 = t595 * t489 + t490 * t585 + t591 * t523;
t590 = Ifges(6,3) + Ifges(7,3);
t594 = mrSges(6,1) * t424 + mrSges(7,1) * t419 - mrSges(6,2) * t425 - mrSges(7,2) * t421 + pkin(5) * t416 + t444 * t591 + t445 * t592 + t580 * t489 + t490 * t589 + t499 * t590;
t467 = -mrSges(6,1) * t489 + mrSges(6,2) * t490;
t473 = -mrSges(6,2) * t523 + mrSges(6,3) * t489;
t408 = m(6) * t424 + mrSges(6,1) * t499 + t473 * t523 + (-t466 - t467) * t490 + (-mrSges(6,3) - mrSges(7,3)) * t445 + t575;
t475 = mrSges(7,1) * t523 - mrSges(7,3) * t490;
t476 = mrSges(6,1) * t523 - mrSges(6,3) * t490;
t574 = m(7) * t421 + t444 * mrSges(7,3) + t489 * t466;
t411 = m(6) * t425 + t444 * mrSges(6,3) + t467 * t489 + (-t475 - t476) * t523 + (-mrSges(6,2) - mrSges(7,2)) * t499 + t574;
t406 = t555 * t408 + t550 * t411;
t478 = Ifges(5,4) * t514 + Ifges(5,2) * t513 + Ifges(5,6) * t525;
t479 = Ifges(5,1) * t514 + Ifges(5,4) * t513 + Ifges(5,5) * t525;
t588 = mrSges(5,1) * t433 - mrSges(5,2) * t434 + Ifges(5,5) * t471 + Ifges(5,6) * t470 + Ifges(5,3) * t504 + pkin(4) * t406 + t514 * t478 - t513 * t479 + t594;
t583 = t548 * t553;
t491 = -mrSges(5,1) * t513 + mrSges(5,2) * t514;
t493 = -mrSges(5,2) * t525 + mrSges(5,3) * t513;
t403 = m(5) * t433 + mrSges(5,1) * t504 - mrSges(5,3) * t471 - t491 * t514 + t493 * t525 + t406;
t494 = mrSges(5,1) * t525 - mrSges(5,3) * t514;
t568 = -t408 * t550 + t555 * t411;
t404 = m(5) * t434 - mrSges(5,2) * t504 + mrSges(5,3) * t470 + t491 * t513 - t494 * t525 + t568;
t400 = -t403 * t551 + t556 * t404;
t510 = -mrSges(4,1) * t526 + mrSges(4,2) * t527;
t516 = mrSges(4,1) * t542 - mrSges(4,3) * t527;
t398 = m(4) * t456 - mrSges(4,2) * t531 + mrSges(4,3) * t506 + t510 * t526 - t516 * t542 + t400;
t455 = -t552 * t486 + t487 * t557;
t447 = -pkin(3) * t531 - pkin(10) * t540 + t527 * t511 - t455;
t435 = -pkin(4) * t470 - pkin(11) * t512 + t514 * t495 + t447;
t427 = -pkin(5) * t444 - qJ(6) * t488 + t474 * t490 + qJDD(6) + t435;
t422 = m(7) * t427 - t444 * mrSges(7,1) + t445 * mrSges(7,2) - t489 * t472 + t490 * t475;
t565 = m(6) * t435 - t444 * mrSges(6,1) + mrSges(6,2) * t445 - t489 * t473 + t476 * t490 + t422;
t414 = -m(5) * t447 + t470 * mrSges(5,1) - mrSges(5,2) * t471 + t513 * t493 - t494 * t514 - t565;
t515 = -mrSges(4,2) * t542 + mrSges(4,3) * t526;
t413 = m(4) * t455 + mrSges(4,1) * t531 - mrSges(4,3) * t507 - t510 * t527 + t515 * t542 + t414;
t394 = t552 * t398 + t557 * t413;
t581 = -t489 * t591 - t490 * t592 - t523 * t590;
t569 = t557 * t398 - t413 * t552;
t399 = t403 * t556 + t404 * t551;
t564 = -m(4) * t485 + t506 * mrSges(4,1) - mrSges(4,2) * t507 + t526 * t515 - t516 * t527 - t399;
t401 = -mrSges(6,1) * t435 + mrSges(6,3) * t425 - mrSges(7,1) * t427 + mrSges(7,3) * t421 - pkin(5) * t422 + qJ(6) * t574 + (-qJ(6) * t475 - t580) * t523 + (-mrSges(7,2) * qJ(6) + t591) * t499 + t581 * t490 + t585 * t445 + t595 * t444;
t405 = mrSges(6,2) * t435 + mrSges(7,2) * t427 - mrSges(6,3) * t424 - mrSges(7,3) * t419 - qJ(6) * t416 + t585 * t444 + t445 * t593 - t581 * t489 + t592 * t499 - t589 * t523;
t477 = Ifges(5,5) * t514 + Ifges(5,6) * t513 + Ifges(5,3) * t525;
t391 = -mrSges(5,1) * t447 + mrSges(5,3) * t434 + Ifges(5,4) * t471 + Ifges(5,2) * t470 + Ifges(5,6) * t504 - pkin(4) * t565 + pkin(11) * t568 + t555 * t401 + t550 * t405 - t514 * t477 + t525 * t479;
t392 = mrSges(5,2) * t447 - mrSges(5,3) * t433 + Ifges(5,1) * t471 + Ifges(5,4) * t470 + Ifges(5,5) * t504 - pkin(11) * t406 - t401 * t550 + t405 * t555 + t477 * t513 - t478 * t525;
t501 = Ifges(4,4) * t527 + Ifges(4,2) * t526 + Ifges(4,6) * t542;
t502 = Ifges(4,1) * t527 + Ifges(4,4) * t526 + Ifges(4,5) * t542;
t562 = mrSges(4,1) * t455 - mrSges(4,2) * t456 + Ifges(4,5) * t507 + Ifges(4,6) * t506 + Ifges(4,3) * t531 + pkin(3) * t414 + pkin(10) * t400 + t556 * t391 + t551 * t392 + t527 * t501 - t526 * t502;
t536 = (-mrSges(3,1) * t558 + mrSges(3,2) * t553) * t578;
t533 = -mrSges(3,2) * t545 + mrSges(3,3) * t572;
t532 = mrSges(3,1) * t545 - mrSges(3,3) * t573;
t520 = -t548 * t534 - t586;
t519 = Ifges(3,5) * t545 + (Ifges(3,1) * t553 + Ifges(3,4) * t558) * t578;
t518 = Ifges(3,6) * t545 + (Ifges(3,4) * t553 + Ifges(3,2) * t558) * t578;
t517 = Ifges(3,3) * t545 + (Ifges(3,5) * t553 + Ifges(3,6) * t558) * t578;
t509 = -g(3) * t583 + t579;
t500 = Ifges(4,5) * t527 + Ifges(4,6) * t526 + Ifges(4,3) * t542;
t395 = m(3) * t508 + mrSges(3,1) * t544 - mrSges(3,3) * t538 + t533 * t545 - t536 * t573 + t564;
t393 = m(3) * t509 - mrSges(3,2) * t544 - mrSges(3,3) * t539 - t532 * t545 + t536 * t572 + t569;
t390 = -mrSges(4,1) * t485 + mrSges(4,3) * t456 + Ifges(4,4) * t507 + Ifges(4,2) * t506 + Ifges(4,6) * t531 - pkin(3) * t399 - t527 * t500 + t542 * t502 - t588;
t389 = mrSges(4,2) * t485 - mrSges(4,3) * t455 + Ifges(4,1) * t507 + Ifges(4,4) * t506 + Ifges(4,5) * t531 - pkin(10) * t399 - t391 * t551 + t392 * t556 + t500 * t526 - t501 * t542;
t388 = Ifges(3,5) * t538 - Ifges(3,6) * t539 + Ifges(3,3) * t544 + mrSges(3,1) * t508 - mrSges(3,2) * t509 + t552 * t389 + t557 * t390 + pkin(2) * t564 + pkin(9) * t569 + (t518 * t553 - t519 * t558) * t578;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t571 - mrSges(2,2) * t567 + (mrSges(3,2) * t520 - mrSges(3,3) * t508 + Ifges(3,1) * t538 - Ifges(3,4) * t539 + Ifges(3,5) * t544 - pkin(9) * t394 + t389 * t557 - t390 * t552 + t517 * t572 - t518 * t545) * t583 + (-mrSges(3,1) * t520 + mrSges(3,3) * t509 + Ifges(3,4) * t538 - Ifges(3,2) * t539 + Ifges(3,6) * t544 - pkin(2) * t394 - t517 * t573 + t545 * t519 - t562) * t582 + t549 * t388 + pkin(1) * ((t393 * t553 + t395 * t558) * t549 + (-m(3) * t520 - t539 * mrSges(3,1) - t538 * mrSges(3,2) + (-t532 * t553 + t533 * t558) * t578 - t394) * t548) + (t393 * t558 - t395 * t553) * t587; t388; t562; t588; t594; t422;];
tauJ  = t1;
