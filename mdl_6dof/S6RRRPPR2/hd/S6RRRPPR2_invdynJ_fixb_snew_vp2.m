% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 04:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:20:29
% EndTime: 2019-05-07 04:20:36
% DurationCPUTime: 5.54s
% Computational Cost: add. (47414->332), mult. (108144->404), div. (0->0), fcn. (78474->10), ass. (0->136)
t584 = Ifges(5,4) + Ifges(6,6);
t596 = -Ifges(5,2) - Ifges(6,3);
t591 = Ifges(5,6) - Ifges(6,5);
t546 = sin(qJ(3));
t547 = sin(qJ(2));
t550 = cos(qJ(3));
t551 = cos(qJ(2));
t525 = (t546 * t551 + t547 * t550) * qJD(1);
t573 = qJD(1) * qJD(2);
t530 = qJDD(1) * t547 + t551 * t573;
t531 = qJDD(1) * t551 - t547 * t573;
t492 = -qJD(3) * t525 - t530 * t546 + t531 * t550;
t576 = qJD(1) * t547;
t534 = qJD(2) * pkin(2) - pkin(8) * t576;
t543 = t551 ^ 2;
t553 = qJD(1) ^ 2;
t548 = sin(qJ(1));
t552 = cos(qJ(1));
t571 = t548 * g(1) - t552 * g(2);
t564 = -qJDD(1) * pkin(1) - t571;
t494 = -t531 * pkin(2) + t534 * t576 + (-pkin(8) * t543 - pkin(7)) * t553 + t564;
t542 = qJD(2) + qJD(3);
t516 = pkin(3) * t542 - qJ(4) * t525;
t524 = (-t546 * t547 + t550 * t551) * qJD(1);
t520 = t524 ^ 2;
t443 = -t492 * pkin(3) - t520 * qJ(4) + t525 * t516 + qJDD(4) + t494;
t493 = qJD(3) * t524 + t530 * t550 + t531 * t546;
t544 = sin(pkin(10));
t583 = cos(pkin(10));
t464 = -t583 * t492 + t493 * t544;
t465 = t544 * t492 + t583 * t493;
t511 = t544 * t524 + t583 * t525;
t496 = mrSges(5,1) * t542 - mrSges(5,3) * t511;
t510 = -t583 * t524 + t525 * t544;
t582 = t510 * t542;
t587 = -2 * qJD(5);
t555 = (-t465 + t582) * qJ(5) + t443 + (t542 * pkin(4) + t587) * t511;
t427 = t464 * pkin(4) + t555;
t498 = mrSges(6,1) * t511 + mrSges(6,2) * t542;
t566 = -g(1) * t552 - g(2) * t548;
t527 = -pkin(1) * t553 + qJDD(1) * pkin(7) + t566;
t581 = t527 * t547;
t586 = pkin(2) * t553;
t484 = qJDD(2) * pkin(2) - pkin(8) * t530 - t581 + (pkin(8) * t573 + t547 * t586 - g(3)) * t551;
t514 = -g(3) * t547 + t551 * t527;
t485 = pkin(8) * t531 - qJD(2) * t534 - t543 * t586 + t514;
t453 = t550 * t484 - t485 * t546;
t541 = qJDD(2) + qJDD(3);
t432 = (t524 * t542 - t493) * qJ(4) + (t524 * t525 + t541) * pkin(3) + t453;
t454 = t546 * t484 + t550 * t485;
t435 = -pkin(3) * t520 + qJ(4) * t492 - t516 * t542 + t454;
t595 = -2 * qJD(4);
t429 = t583 * t432 - t544 * t435 + t511 * t595;
t477 = pkin(4) * t510 - qJ(5) * t511;
t540 = t542 ^ 2;
t425 = -t541 * pkin(4) - t540 * qJ(5) + t511 * t477 + qJDD(5) - t429;
t419 = (t510 * t511 - t541) * pkin(9) + (t465 + t582) * pkin(5) + t425;
t499 = pkin(5) * t511 - pkin(9) * t542;
t508 = t510 ^ 2;
t422 = (pkin(4) + pkin(9)) * t464 + t555 - t511 * t499 - t508 * pkin(5);
t545 = sin(qJ(6));
t549 = cos(qJ(6));
t417 = t419 * t549 - t422 * t545;
t490 = t510 * t549 - t542 * t545;
t441 = qJD(6) * t490 + t464 * t545 + t541 * t549;
t463 = qJDD(6) + t465;
t491 = t510 * t545 + t542 * t549;
t466 = -mrSges(7,1) * t490 + mrSges(7,2) * t491;
t504 = qJD(6) + t511;
t467 = -mrSges(7,2) * t504 + mrSges(7,3) * t490;
t414 = m(7) * t417 + mrSges(7,1) * t463 - mrSges(7,3) * t441 - t466 * t491 + t467 * t504;
t418 = t419 * t545 + t422 * t549;
t440 = -qJD(6) * t491 + t464 * t549 - t541 * t545;
t468 = mrSges(7,1) * t504 - mrSges(7,3) * t491;
t415 = m(7) * t418 - mrSges(7,2) * t463 + mrSges(7,3) * t440 + t466 * t490 - t468 * t504;
t567 = -t545 * t414 + t549 * t415;
t563 = -m(6) * t427 + t465 * mrSges(6,3) + t511 * t498 - t567;
t497 = mrSges(6,1) * t510 - mrSges(6,3) * t542;
t577 = -mrSges(5,2) * t542 - mrSges(5,3) * t510 - t497;
t585 = mrSges(6,2) - mrSges(5,1);
t401 = m(5) * t443 + t465 * mrSges(5,2) - t585 * t464 + t511 * t496 + t577 * t510 - t563;
t594 = Ifges(5,1) + Ifges(6,2);
t593 = -Ifges(6,1) - Ifges(5,3);
t592 = Ifges(5,5) - Ifges(6,4);
t590 = t596 * t510 + t511 * t584 + t591 * t542;
t515 = -mrSges(4,2) * t542 + mrSges(4,3) * t524;
t517 = mrSges(4,1) * t542 - mrSges(4,3) * t525;
t589 = m(4) * t494 - t492 * mrSges(4,1) + t493 * mrSges(4,2) - t524 * t515 + t525 * t517 + t401;
t478 = mrSges(5,1) * t510 + mrSges(5,2) * t511;
t405 = t414 * t549 + t415 * t545;
t479 = -mrSges(6,2) * t510 - mrSges(6,3) * t511;
t560 = -m(6) * t425 - t465 * mrSges(6,1) - t511 * t479 - t405;
t400 = m(5) * t429 - mrSges(5,3) * t465 - t478 * t511 - t585 * t541 + t577 * t542 + t560;
t502 = t510 * t595;
t580 = t544 * t432 + t583 * t435;
t430 = t502 + t580;
t562 = pkin(4) * t540 - qJ(5) * t541 - t580;
t423 = t542 * t587 + ((2 * qJD(4)) + t477) * t510 + t562;
t421 = -pkin(5) * t464 - pkin(9) * t508 - t477 * t510 + t502 + ((2 * qJD(5)) + t499) * t542 - t562;
t561 = -m(7) * t421 + mrSges(7,1) * t440 - t441 * mrSges(7,2) + t467 * t490 - t491 * t468;
t557 = -m(6) * t423 + t541 * mrSges(6,3) + t542 * t498 - t561;
t411 = m(5) * t430 - mrSges(5,2) * t541 - t496 * t542 + (-t478 - t479) * t510 + (-mrSges(5,3) - mrSges(6,1)) * t464 + t557;
t398 = t583 * t400 + t544 * t411;
t512 = -mrSges(4,1) * t524 + mrSges(4,2) * t525;
t395 = m(4) * t453 + mrSges(4,1) * t541 - mrSges(4,3) * t493 - t512 * t525 + t515 * t542 + t398;
t568 = -t400 * t544 + t583 * t411;
t396 = m(4) * t454 - mrSges(4,2) * t541 + mrSges(4,3) * t492 + t512 * t524 - t517 * t542 + t568;
t390 = t550 * t395 + t546 * t396;
t579 = t591 * t510 - t592 * t511 + t593 * t542;
t578 = -t584 * t510 + t594 * t511 + t592 * t542;
t575 = qJD(1) * t551;
t569 = -t395 * t546 + t550 * t396;
t445 = Ifges(7,4) * t491 + Ifges(7,2) * t490 + Ifges(7,6) * t504;
t446 = Ifges(7,1) * t491 + Ifges(7,4) * t490 + Ifges(7,5) * t504;
t558 = mrSges(7,1) * t417 - mrSges(7,2) * t418 + Ifges(7,5) * t441 + Ifges(7,6) * t440 + Ifges(7,3) * t463 + t491 * t445 - t490 * t446;
t403 = mrSges(6,2) * t541 + t497 * t542 - t560;
t444 = Ifges(7,5) * t491 + Ifges(7,6) * t490 + Ifges(7,3) * t504;
t407 = -mrSges(7,1) * t421 + mrSges(7,3) * t418 + Ifges(7,4) * t441 + Ifges(7,2) * t440 + Ifges(7,6) * t463 - t444 * t491 + t446 * t504;
t408 = mrSges(7,2) * t421 - mrSges(7,3) * t417 + Ifges(7,1) * t441 + Ifges(7,4) * t440 + Ifges(7,5) * t463 + t444 * t490 - t445 * t504;
t506 = Ifges(4,4) * t525 + Ifges(4,2) * t524 + Ifges(4,6) * t542;
t507 = Ifges(4,1) * t525 + Ifges(4,4) * t524 + Ifges(4,5) * t542;
t554 = t578 * t510 + qJ(5) * (-t479 * t510 + t557) + t549 * t408 - t545 * t407 - t524 * t507 + t525 * t506 + Ifges(4,5) * t493 + Ifges(4,6) * t492 + mrSges(4,1) * t453 - mrSges(4,2) * t454 + mrSges(5,1) * t429 - mrSges(5,2) * t430 - mrSges(6,3) * t423 + mrSges(6,2) * t425 - pkin(9) * t405 - pkin(4) * t403 + pkin(3) * t398 + t590 * t511 + t592 * t465 + (-qJ(5) * mrSges(6,1) - t591) * t464 + (Ifges(4,3) - t593) * t541;
t533 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t575;
t532 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t576;
t529 = (-mrSges(3,1) * t551 + mrSges(3,2) * t547) * qJD(1);
t526 = -pkin(7) * t553 + t564;
t523 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t547 + Ifges(3,4) * t551) * qJD(1);
t522 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t547 + Ifges(3,2) * t551) * qJD(1);
t513 = -g(3) * t551 - t581;
t505 = Ifges(4,5) * t525 + Ifges(4,6) * t524 + Ifges(4,3) * t542;
t404 = -t464 * mrSges(6,2) - t510 * t497 - t563;
t391 = mrSges(6,1) * t425 + mrSges(5,2) * t443 - mrSges(5,3) * t429 - mrSges(6,3) * t427 + pkin(5) * t405 - qJ(5) * t404 - t584 * t464 + t594 * t465 + t579 * t510 + t592 * t541 - t590 * t542 + t558;
t389 = -mrSges(5,1) * t443 - mrSges(6,1) * t423 + mrSges(6,2) * t427 + mrSges(5,3) * t430 - pkin(4) * t404 - pkin(5) * t561 - pkin(9) * t567 - t549 * t407 - t545 * t408 + t464 * t596 + t584 * t465 + t579 * t511 + t591 * t541 + t578 * t542;
t388 = mrSges(4,2) * t494 - mrSges(4,3) * t453 + Ifges(4,1) * t493 + Ifges(4,4) * t492 + Ifges(4,5) * t541 - qJ(4) * t398 - t544 * t389 + t583 * t391 + t524 * t505 - t542 * t506;
t387 = -mrSges(4,1) * t494 + mrSges(4,3) * t454 + Ifges(4,4) * t493 + Ifges(4,2) * t492 + Ifges(4,6) * t541 - pkin(3) * t401 + qJ(4) * t568 + t583 * t389 + t544 * t391 - t525 * t505 + t542 * t507;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t571 - mrSges(2,2) * t566 + t547 * (mrSges(3,2) * t526 - mrSges(3,3) * t513 + Ifges(3,1) * t530 + Ifges(3,4) * t531 + Ifges(3,5) * qJDD(2) - pkin(8) * t390 - qJD(2) * t522 - t546 * t387 + t550 * t388) + t551 * (-mrSges(3,1) * t526 + mrSges(3,3) * t514 + Ifges(3,4) * t530 + Ifges(3,2) * t531 + Ifges(3,6) * qJDD(2) - pkin(2) * t589 + pkin(8) * t569 + qJD(2) * t523 + t550 * t387 + t546 * t388) + pkin(1) * (-t530 * mrSges(3,2) + t531 * mrSges(3,1) - m(3) * t526 + (-t547 * t532 + t551 * t533) * qJD(1) - t589) + pkin(7) * (t551 * (m(3) * t514 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t531 - qJD(2) * t532 + t529 * t575 + t569) - t547 * (m(3) * t513 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t530 + qJD(2) * t533 - t529 * t576 + t390)); Ifges(3,3) * qJDD(2) + t554 + pkin(2) * t390 + (t547 * t522 - t551 * t523) * qJD(1) + Ifges(3,5) * t530 + Ifges(3,6) * t531 + mrSges(3,1) * t513 - mrSges(3,2) * t514; t554; t401; t403; t558;];
tauJ  = t1;
