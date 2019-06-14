% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 21:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPPRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:50:27
% EndTime: 2019-05-04 21:50:34
% DurationCPUTime: 4.39s
% Computational Cost: add. (57891->250), mult. (101778->312), div. (0->0), fcn. (67096->12), ass. (0->114)
t521 = sin(pkin(10));
t524 = cos(pkin(10));
t506 = t521 * g(1) - t524 * g(2);
t507 = -t524 * g(1) - t521 * g(2);
t517 = -g(3) + qJDD(1);
t528 = sin(qJ(2));
t525 = cos(pkin(6));
t531 = cos(qJ(2));
t554 = t525 * t531;
t522 = sin(pkin(6));
t556 = t522 * t531;
t473 = t506 * t554 - t528 * t507 + t517 * t556;
t471 = qJDD(2) * pkin(2) + t473;
t555 = t525 * t528;
t557 = t522 * t528;
t474 = t506 * t555 + t531 * t507 + t517 * t557;
t533 = qJD(2) ^ 2;
t472 = -t533 * pkin(2) + t474;
t520 = sin(pkin(11));
t523 = cos(pkin(11));
t467 = t520 * t471 + t523 * t472;
t562 = -qJDD(2) * qJ(4) - (2 * qJD(4) * qJD(2)) - t467;
t561 = -pkin(3) - pkin(8);
t560 = mrSges(4,1) - mrSges(5,2);
t559 = -Ifges(5,4) + Ifges(4,5);
t558 = Ifges(5,5) - Ifges(4,6);
t489 = -t522 * t506 + t525 * t517;
t488 = qJDD(3) + t489;
t527 = sin(qJ(5));
t553 = t527 * t488;
t466 = t523 * t471 - t520 * t472;
t538 = -t533 * qJ(4) + qJDD(4) - t466;
t463 = t561 * qJDD(2) + t538;
t530 = cos(qJ(5));
t459 = t527 * t463 + t530 * t488;
t502 = (mrSges(6,1) * t527 + mrSges(6,2) * t530) * qJD(2);
t549 = qJD(2) * qJD(5);
t545 = t530 * t549;
t504 = -t527 * qJDD(2) - t545;
t551 = qJD(2) * t530;
t509 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t551;
t503 = (pkin(5) * t527 - pkin(9) * t530) * qJD(2);
t532 = qJD(5) ^ 2;
t550 = t527 * qJD(2);
t456 = -t532 * pkin(5) + qJDD(5) * pkin(9) - t503 * t550 + t459;
t462 = t561 * t533 - t562;
t546 = t527 * t549;
t505 = t530 * qJDD(2) - t546;
t457 = (-t505 + t546) * pkin(9) + (-t504 + t545) * pkin(5) + t462;
t526 = sin(qJ(6));
t529 = cos(qJ(6));
t453 = -t526 * t456 + t529 * t457;
t500 = t529 * qJD(5) - t526 * t551;
t481 = t500 * qJD(6) + t526 * qJDD(5) + t529 * t505;
t501 = t526 * qJD(5) + t529 * t551;
t482 = -t500 * mrSges(7,1) + t501 * mrSges(7,2);
t512 = qJD(6) + t550;
t486 = -t512 * mrSges(7,2) + t500 * mrSges(7,3);
t498 = qJDD(6) - t504;
t451 = m(7) * t453 + t498 * mrSges(7,1) - t481 * mrSges(7,3) - t501 * t482 + t512 * t486;
t454 = t529 * t456 + t526 * t457;
t480 = -t501 * qJD(6) + t529 * qJDD(5) - t526 * t505;
t487 = t512 * mrSges(7,1) - t501 * mrSges(7,3);
t452 = m(7) * t454 - t498 * mrSges(7,2) + t480 * mrSges(7,3) + t500 * t482 - t512 * t487;
t541 = -t526 * t451 + t529 * t452;
t443 = m(6) * t459 - qJDD(5) * mrSges(6,2) + t504 * mrSges(6,3) - qJD(5) * t509 - t502 * t550 + t541;
t458 = t530 * t463 - t553;
t508 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t550;
t455 = -qJDD(5) * pkin(5) - t532 * pkin(9) + t553 + (qJD(2) * t503 - t463) * t530;
t535 = -m(7) * t455 + t480 * mrSges(7,1) - t481 * mrSges(7,2) + t500 * t486 - t501 * t487;
t447 = m(6) * t458 + qJDD(5) * mrSges(6,1) - t505 * mrSges(6,3) + qJD(5) * t508 - t502 * t551 + t535;
t438 = t527 * t443 + t530 * t447;
t465 = -qJDD(2) * pkin(3) + t538;
t537 = -m(5) * t465 + t533 * mrSges(5,3) - t438;
t434 = m(4) * t466 - t533 * mrSges(4,2) + t560 * qJDD(2) + t537;
t464 = t533 * pkin(3) + t562;
t444 = t529 * t451 + t526 * t452;
t536 = -m(6) * t462 + t504 * mrSges(6,1) - t505 * mrSges(6,2) - t508 * t550 - t509 * t551 - t444;
t534 = -m(5) * t464 + t533 * mrSges(5,2) + qJDD(2) * mrSges(5,3) - t536;
t441 = m(4) * t467 - t533 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t534;
t430 = t523 * t434 + t520 * t441;
t428 = m(3) * t473 + qJDD(2) * mrSges(3,1) - t533 * mrSges(3,2) + t430;
t543 = -t520 * t434 + t523 * t441;
t429 = m(3) * t474 - t533 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t543;
t542 = t530 * t443 - t527 * t447;
t437 = m(5) * t488 + t542;
t540 = m(4) * t488 + t437;
t436 = m(3) * t489 + t540;
t417 = t428 * t554 + t429 * t555 - t522 * t436;
t415 = m(2) * t506 + t417;
t422 = -t528 * t428 + t531 * t429;
t421 = m(2) * t507 + t422;
t552 = t524 * t415 + t521 * t421;
t416 = t428 * t556 + t429 * t557 + t525 * t436;
t544 = -t521 * t415 + t524 * t421;
t475 = Ifges(7,5) * t501 + Ifges(7,6) * t500 + Ifges(7,3) * t512;
t477 = Ifges(7,1) * t501 + Ifges(7,4) * t500 + Ifges(7,5) * t512;
t445 = -mrSges(7,1) * t455 + mrSges(7,3) * t454 + Ifges(7,4) * t481 + Ifges(7,2) * t480 + Ifges(7,6) * t498 - t501 * t475 + t512 * t477;
t476 = Ifges(7,4) * t501 + Ifges(7,2) * t500 + Ifges(7,6) * t512;
t446 = mrSges(7,2) * t455 - mrSges(7,3) * t453 + Ifges(7,1) * t481 + Ifges(7,4) * t480 + Ifges(7,5) * t498 + t500 * t475 - t512 * t476;
t492 = (Ifges(6,3) * qJD(5)) + (Ifges(6,5) * t530 - Ifges(6,6) * t527) * qJD(2);
t493 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t530 - Ifges(6,2) * t527) * qJD(2);
t431 = mrSges(6,2) * t462 - mrSges(6,3) * t458 + Ifges(6,1) * t505 + Ifges(6,4) * t504 + Ifges(6,5) * qJDD(5) - pkin(9) * t444 - qJD(5) * t493 - t526 * t445 + t529 * t446 - t492 * t550;
t494 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t530 - Ifges(6,4) * t527) * qJD(2);
t432 = -mrSges(6,1) * t462 - mrSges(7,1) * t453 + mrSges(7,2) * t454 + mrSges(6,3) * t459 + Ifges(6,4) * t505 - Ifges(7,5) * t481 + Ifges(6,2) * t504 + Ifges(6,6) * qJDD(5) - Ifges(7,6) * t480 - Ifges(7,3) * t498 - pkin(5) * t444 + qJD(5) * t494 - t501 * t476 + t500 * t477 - t492 * t551;
t413 = -mrSges(5,1) * t464 + mrSges(4,3) * t467 - pkin(3) * t437 - pkin(4) * t536 - pkin(8) * t542 - t558 * qJDD(2) - t527 * t431 - t530 * t432 - t560 * t488 + t559 * t533;
t418 = Ifges(6,5) * t505 + Ifges(6,6) * t504 + Ifges(6,3) * qJDD(5) + mrSges(6,1) * t458 - mrSges(6,2) * t459 + t526 * t446 + t529 * t445 + pkin(5) * t535 + pkin(9) * t541 + mrSges(5,1) * t465 + pkin(4) * t438 - qJ(4) * t437 - mrSges(4,3) * t466 + t558 * t533 + (mrSges(4,2) - mrSges(5,3)) * t488 + t559 * qJDD(2) + (t530 * t493 + t527 * t494) * qJD(2);
t410 = -mrSges(3,1) * t489 + mrSges(3,3) * t474 + t533 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t540 + qJ(3) * t543 + t523 * t413 + t520 * t418;
t411 = mrSges(3,2) * t489 - mrSges(3,3) * t473 + Ifges(3,5) * qJDD(2) - t533 * Ifges(3,6) - qJ(3) * t430 - t520 * t413 + t523 * t418;
t539 = pkin(7) * t422 + t410 * t531 + t411 * t528;
t412 = pkin(2) * t430 + mrSges(3,1) * t473 - mrSges(3,2) * t474 + pkin(3) * t537 + qJ(4) * t534 + t530 * t431 - t527 * t432 - pkin(8) * t438 + mrSges(4,1) * t466 - mrSges(4,2) * t467 + mrSges(5,2) * t465 - mrSges(5,3) * t464 + (-pkin(3) * mrSges(5,2) + Ifges(5,1) + Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t409 = mrSges(2,2) * t517 - mrSges(2,3) * t506 - t528 * t410 + t531 * t411 + (-t416 * t522 - t417 * t525) * pkin(7);
t408 = -mrSges(2,1) * t517 + mrSges(2,3) * t507 - pkin(1) * t416 - t522 * t412 + t539 * t525;
t1 = [-m(1) * g(1) + t544; -m(1) * g(2) + t552; -m(1) * g(3) + m(2) * t517 + t416; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t552 - t521 * t408 + t524 * t409; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t544 + t524 * t408 + t521 * t409; -mrSges(1,1) * g(2) + mrSges(2,1) * t506 + mrSges(1,2) * g(1) - mrSges(2,2) * t507 + pkin(1) * t417 + t525 * t412 + t539 * t522;];
tauB  = t1;
