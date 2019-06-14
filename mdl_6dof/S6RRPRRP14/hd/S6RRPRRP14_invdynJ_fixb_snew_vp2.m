% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRP14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 19:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRP14_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP14_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP14_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:18:01
% EndTime: 2019-05-06 19:18:11
% DurationCPUTime: 4.11s
% Computational Cost: add. (34925->322), mult. (77699->390), div. (0->0), fcn. (55962->10), ass. (0->142)
t617 = -2 * qJD(3);
t616 = Ifges(3,1) + Ifges(4,2);
t615 = Ifges(6,1) + Ifges(7,1);
t602 = Ifges(3,4) + Ifges(4,6);
t601 = Ifges(6,4) - Ifges(7,5);
t600 = Ifges(3,5) - Ifges(4,4);
t599 = -Ifges(6,5) - Ifges(7,4);
t614 = Ifges(3,2) + Ifges(4,3);
t613 = Ifges(6,2) + Ifges(7,3);
t598 = Ifges(3,6) - Ifges(4,5);
t597 = Ifges(6,6) - Ifges(7,6);
t612 = Ifges(3,3) + Ifges(4,1);
t611 = -Ifges(6,3) - Ifges(7,2);
t552 = cos(pkin(6));
t546 = qJD(1) * t552 + qJD(2);
t555 = sin(qJ(2));
t551 = sin(pkin(6));
t583 = qJD(1) * t551;
t576 = t555 * t583;
t610 = (pkin(2) * t546 + t617) * t576;
t560 = qJD(1) ^ 2;
t556 = sin(qJ(1));
t559 = cos(qJ(1));
t570 = -g(1) * t559 - g(2) * t556;
t580 = qJDD(1) * t551;
t528 = -pkin(1) * t560 + pkin(8) * t580 + t570;
t558 = cos(qJ(2));
t593 = t551 * t555;
t574 = t556 * g(1) - g(2) * t559;
t605 = t551 * pkin(8);
t527 = qJDD(1) * pkin(1) + t560 * t605 + t574;
t595 = t527 * t552;
t491 = -g(3) * t593 + t558 * t528 + t555 * t595;
t529 = (-t558 * pkin(2) - t555 * qJ(3)) * t583;
t544 = t546 ^ 2;
t545 = qJDD(1) * t552 + qJDD(2);
t582 = qJD(1) * t558;
t575 = t551 * t582;
t456 = pkin(2) * t544 - t545 * qJ(3) - t529 * t575 + t546 * t617 - t491;
t554 = sin(qJ(4));
t557 = cos(qJ(4));
t515 = -t546 * t554 - t557 * t575;
t534 = -qJD(2) * t576 + t558 * t580;
t489 = qJD(4) * t515 - t534 * t554 + t545 * t557;
t516 = t546 * t557 - t554 * t575;
t539 = qJD(4) + t576;
t553 = sin(qJ(5));
t607 = cos(qJ(5));
t494 = t516 * t553 - t607 * t539;
t533 = (qJD(2) * t582 + qJDD(1) * t555) * t551;
t522 = qJDD(4) + t533;
t454 = -t494 * qJD(5) + t607 * t489 + t553 * t522;
t495 = t607 * t516 + t553 * t539;
t470 = mrSges(7,1) * t494 - mrSges(7,3) * t495;
t532 = pkin(3) * t576 - pkin(9) * t546;
t594 = t551 ^ 2 * t560;
t578 = t558 ^ 2 * t594;
t606 = g(3) * t552;
t608 = -pkin(2) - pkin(9);
t448 = -pkin(3) * t578 - t606 - qJ(3) * t533 + t608 * t534 + (-t527 + (-qJ(3) * t546 * t558 - t532 * t555) * qJD(1)) * t551 + t610;
t592 = t551 * t558;
t584 = g(3) * t592 + t555 * t528;
t568 = -qJ(3) * t544 + t529 * t576 + qJDD(3) + t584;
t450 = pkin(3) * t533 + t608 * t545 + (-pkin(3) * t546 * t583 - pkin(9) * t555 * t594 - t595) * t558 + t568;
t443 = t557 * t448 + t554 * t450;
t493 = -pkin(4) * t515 - pkin(10) * t516;
t537 = t539 ^ 2;
t438 = -pkin(4) * t537 + pkin(10) * t522 + t493 * t515 + t443;
t447 = pkin(3) * t534 - pkin(9) * t578 + t546 * t532 - t456;
t488 = -qJD(4) * t516 - t534 * t557 - t545 * t554;
t440 = (-t515 * t539 - t489) * pkin(10) + (t516 * t539 - t488) * pkin(4) + t447;
t434 = -t553 * t438 + t607 * t440;
t469 = pkin(5) * t494 - qJ(6) * t495;
t486 = qJDD(5) - t488;
t513 = qJD(5) - t515;
t512 = t513 ^ 2;
t432 = -t486 * pkin(5) - t512 * qJ(6) + t495 * t469 + qJDD(6) - t434;
t474 = -mrSges(7,2) * t494 + mrSges(7,3) * t513;
t571 = -m(7) * t432 + t486 * mrSges(7,1) + t513 * t474;
t428 = mrSges(7,2) * t454 + t470 * t495 - t571;
t435 = t607 * t438 + t553 * t440;
t431 = -pkin(5) * t512 + qJ(6) * t486 + 0.2e1 * qJD(6) * t513 - t469 * t494 + t435;
t453 = qJD(5) * t495 + t489 * t553 - t607 * t522;
t477 = -mrSges(7,1) * t513 + mrSges(7,2) * t495;
t577 = m(7) * t431 + t486 * mrSges(7,3) + t513 * t477;
t589 = t601 * t494 - t615 * t495 + t599 * t513;
t590 = t613 * t494 - t601 * t495 - t597 * t513;
t609 = -t597 * t453 - t599 * t454 - t611 * t486 - t589 * t494 - t590 * t495 + mrSges(6,1) * t434 - mrSges(7,1) * t432 - mrSges(6,2) * t435 + mrSges(7,3) * t431 - pkin(5) * t428 + qJ(6) * (-mrSges(7,2) * t453 - t470 * t494 + t577);
t604 = mrSges(3,1) - mrSges(4,2);
t603 = -mrSges(6,3) - mrSges(7,2);
t476 = mrSges(6,1) * t513 - mrSges(6,3) * t495;
t588 = -mrSges(6,1) * t494 - mrSges(6,2) * t495 - t470;
t423 = m(6) * t435 - mrSges(6,2) * t486 + t603 * t453 - t476 * t513 + t588 * t494 + t577;
t475 = -mrSges(6,2) * t513 - mrSges(6,3) * t494;
t425 = m(6) * t434 + mrSges(6,1) * t486 + t603 * t454 + t475 * t513 + t588 * t495 + t571;
t419 = t553 * t423 + t607 * t425;
t591 = t597 * t494 + t599 * t495 + t611 * t513;
t587 = (t555 * t600 + t558 * t598) * t583 + t612 * t546;
t586 = (t555 * t602 + t558 * t614) * t583 + t598 * t546;
t585 = (t555 * t616 + t558 * t602) * t583 + t600 * t546;
t579 = t558 * t595;
t492 = -mrSges(5,1) * t515 + mrSges(5,2) * t516;
t497 = mrSges(5,1) * t539 - mrSges(5,3) * t516;
t572 = t607 * t423 - t425 * t553;
t413 = m(5) * t443 - mrSges(5,2) * t522 + mrSges(5,3) * t488 + t492 * t515 - t497 * t539 + t572;
t442 = -t554 * t448 + t450 * t557;
t496 = -mrSges(5,2) * t539 + mrSges(5,3) * t515;
t437 = -pkin(4) * t522 - pkin(10) * t537 + t516 * t493 - t442;
t433 = -0.2e1 * qJD(6) * t495 + (t494 * t513 - t454) * qJ(6) + (t495 * t513 + t453) * pkin(5) + t437;
t429 = m(7) * t433 + mrSges(7,1) * t453 - t454 * mrSges(7,3) + t474 * t494 - t495 * t477;
t561 = -m(6) * t437 - t453 * mrSges(6,1) - mrSges(6,2) * t454 - t494 * t475 - t476 * t495 - t429;
t420 = m(5) * t442 + mrSges(5,1) * t522 - mrSges(5,3) * t489 - t492 * t516 + t496 * t539 + t561;
t573 = t557 * t413 - t554 * t420;
t504 = -t527 * t551 - t606;
t410 = t413 * t554 + t420 * t557;
t457 = -pkin(2) * t534 + (-t546 * t575 - t533) * qJ(3) + t504 + t610;
t525 = -mrSges(4,1) * t575 - mrSges(4,3) * t546;
t569 = -m(4) * t457 + t533 * mrSges(4,3) - t525 * t575 - t573;
t466 = -pkin(2) * t545 + t568 - t579;
t567 = -m(4) * t466 - t533 * mrSges(4,1) - t410;
t565 = -m(5) * t447 + t488 * mrSges(5,1) - t489 * mrSges(5,2) + t515 * t496 - t516 * t497 - t419;
t415 = -mrSges(6,1) * t437 - mrSges(7,1) * t433 + mrSges(7,2) * t431 + mrSges(6,3) * t435 - pkin(5) * t429 - t613 * t453 + t601 * t454 + t597 * t486 + t591 * t495 - t589 * t513;
t417 = mrSges(6,2) * t437 + mrSges(7,2) * t432 - mrSges(6,3) * t434 - mrSges(7,3) * t433 - qJ(6) * t429 - t601 * t453 + t615 * t454 - t599 * t486 + t591 * t494 + t590 * t513;
t481 = Ifges(5,4) * t516 + Ifges(5,2) * t515 + Ifges(5,6) * t539;
t482 = Ifges(5,1) * t516 + Ifges(5,4) * t515 + Ifges(5,5) * t539;
t564 = mrSges(5,1) * t442 - mrSges(5,2) * t443 + Ifges(5,5) * t489 + Ifges(5,6) * t488 + Ifges(5,3) * t522 + pkin(4) * t561 + pkin(10) * t572 + t607 * t415 + t553 * t417 + t516 * t481 - t515 * t482;
t526 = mrSges(4,1) * t576 + mrSges(4,2) * t546;
t530 = (t558 * mrSges(4,2) - t555 * mrSges(4,3)) * t583;
t562 = -m(4) * t456 + t545 * mrSges(4,3) + t546 * t526 + t530 * t575 - t565;
t531 = (-t558 * mrSges(3,1) + t555 * mrSges(3,2)) * t583;
t524 = -mrSges(3,2) * t546 + mrSges(3,3) * t575;
t523 = mrSges(3,1) * t546 - mrSges(3,3) * t576;
t490 = t579 - t584;
t480 = Ifges(5,5) * t516 + Ifges(5,6) * t515 + Ifges(5,3) * t539;
t411 = (mrSges(3,3) + mrSges(4,1)) * t534 + t531 * t575 + t562 + m(3) * t491 - mrSges(3,2) * t545 - t523 * t546;
t409 = mrSges(4,2) * t545 + t525 * t546 + t530 * t576 - t567;
t408 = t534 * mrSges(4,2) - t526 * t576 - t569;
t407 = m(3) * t490 - mrSges(3,3) * t533 + (t524 - t525) * t546 + t604 * t545 + (-t530 - t531) * t576 + t567;
t406 = -mrSges(5,1) * t447 + mrSges(5,3) * t443 + Ifges(5,4) * t489 + Ifges(5,2) * t488 + Ifges(5,6) * t522 - pkin(4) * t419 - t516 * t480 + t539 * t482 - t609;
t405 = mrSges(5,2) * t447 - mrSges(5,3) * t442 + Ifges(5,1) * t489 + Ifges(5,4) * t488 + Ifges(5,5) * t522 - pkin(10) * t419 - t553 * t415 + t607 * t417 + t515 * t480 - t539 * t481;
t404 = mrSges(3,1) * t490 - mrSges(3,2) * t491 + mrSges(4,2) * t466 - mrSges(4,3) * t456 + t557 * t405 - t554 * t406 - pkin(9) * t410 - pkin(2) * t409 + qJ(3) * t562 + t612 * t545 + (qJ(3) * mrSges(4,1) + t598) * t534 + t600 * t533 + (t586 * t555 - t585 * t558) * t583;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t574 - mrSges(2,2) * t570 + (mrSges(4,1) * t466 + mrSges(3,2) * t504 - mrSges(3,3) * t490 - mrSges(4,3) * t457 + pkin(3) * t410 - qJ(3) * t408 + t616 * t533 + t602 * t534 + t600 * t545 - t586 * t546 + t587 * t575 + t564) * t593 + (-mrSges(3,1) * t504 - mrSges(4,1) * t456 + mrSges(4,2) * t457 + mrSges(3,3) * t491 - pkin(2) * t408 - pkin(3) * t565 - pkin(9) * t573 - t554 * t405 - t557 * t406 + t602 * t533 + t614 * t534 + t598 * t545 + t585 * t546 - t587 * t576) * t592 + t552 * t404 + pkin(1) * ((t558 * t407 + t555 * t411) * t552 + (-m(3) * t504 - t533 * mrSges(3,2) + t604 * t534 + (t524 * t558 + (-t523 + t526) * t555) * t583 + t569) * t551) + (-t407 * t555 + t411 * t558) * t605; t404; t409; t564; t609; t428;];
tauJ  = t1;
