% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 14:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRR11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR11_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR11_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 14:23:07
% EndTime: 2019-05-07 14:23:22
% DurationCPUTime: 7.31s
% Computational Cost: add. (81486->341), mult. (174559->427), div. (0->0), fcn. (136648->12), ass. (0->149)
t626 = Ifges(4,1) + Ifges(5,1);
t616 = Ifges(4,4) - Ifges(5,5);
t615 = Ifges(4,5) + Ifges(5,4);
t625 = Ifges(4,2) + Ifges(5,3);
t614 = Ifges(5,6) - Ifges(4,6);
t624 = Ifges(4,3) + Ifges(5,2);
t574 = cos(pkin(6));
t570 = qJD(1) * t574 + qJD(2);
t577 = sin(qJ(3));
t578 = sin(qJ(2));
t573 = sin(pkin(6));
t603 = qJD(1) * t573;
t600 = t578 * t603;
t621 = cos(qJ(3));
t541 = -t570 * t621 + t577 * t600;
t582 = cos(qJ(2));
t602 = qJD(1) * t582;
t556 = (qJD(2) * t602 + qJDD(1) * t578) * t573;
t569 = qJDD(1) * t574 + qJDD(2);
t514 = -t541 * qJD(3) + t556 * t621 + t577 * t569;
t542 = t577 * t570 + t600 * t621;
t521 = mrSges(5,1) * t541 - mrSges(5,3) * t542;
t555 = (-pkin(2) * t582 - pkin(9) * t578) * t603;
t568 = t570 ^ 2;
t584 = qJD(1) ^ 2;
t579 = sin(qJ(1));
t583 = cos(qJ(1));
t594 = -g(1) * t583 - g(2) * t579;
t601 = qJDD(1) * t573;
t553 = -pkin(1) * t584 + pkin(8) * t601 + t594;
t598 = t579 * g(1) - g(2) * t583;
t619 = pkin(8) * t573;
t552 = qJDD(1) * pkin(1) + t584 * t619 + t598;
t611 = t552 * t574;
t604 = t582 * t553 + t578 * t611;
t492 = -t568 * pkin(2) + t569 * pkin(9) + (-g(3) * t578 + t555 * t602) * t573 + t604;
t557 = -qJD(2) * t600 + t582 * t601;
t618 = t574 * g(3);
t493 = -t557 * pkin(2) - t556 * pkin(9) - t618 + (-t552 + (pkin(2) * t578 - pkin(9) * t582) * t570 * qJD(1)) * t573;
t464 = -t577 * t492 + t493 * t621;
t520 = pkin(3) * t541 - qJ(4) * t542;
t549 = qJDD(3) - t557;
t599 = t573 * t602;
t564 = qJD(3) - t599;
t563 = t564 ^ 2;
t459 = -t549 * pkin(3) - t563 * qJ(4) + t542 * t520 + qJDD(4) - t464;
t612 = t541 * t564;
t452 = (-t514 - t612) * pkin(10) + (t541 * t542 - t549) * pkin(4) + t459;
t465 = t621 * t492 + t577 * t493;
t622 = 2 * qJD(4);
t458 = -pkin(3) * t563 + t549 * qJ(4) - t541 * t520 + t564 * t622 + t465;
t513 = t542 * qJD(3) + t577 * t556 - t569 * t621;
t529 = -pkin(4) * t564 - pkin(10) * t542;
t540 = t541 ^ 2;
t455 = -pkin(4) * t540 + pkin(10) * t513 + t529 * t564 + t458;
t576 = sin(qJ(5));
t581 = cos(qJ(5));
t450 = t576 * t452 + t581 * t455;
t519 = t541 * t576 + t542 * t581;
t475 = -qJD(5) * t519 + t513 * t581 - t514 * t576;
t518 = t541 * t581 - t542 * t576;
t486 = -mrSges(6,1) * t518 + mrSges(6,2) * t519;
t560 = qJD(5) - t564;
t499 = mrSges(6,1) * t560 - mrSges(6,3) * t519;
t548 = qJDD(5) - t549;
t487 = -pkin(5) * t518 - pkin(11) * t519;
t559 = t560 ^ 2;
t446 = -pkin(5) * t559 + pkin(11) * t548 + t487 * t518 + t450;
t609 = t573 * t582;
t515 = -g(3) * t609 - t578 * t553 + t582 * t611;
t491 = -t569 * pkin(2) - t568 * pkin(9) + t555 * t600 - t515;
t591 = t513 * pkin(3) + t491 + (-t514 + t612) * qJ(4);
t620 = pkin(3) * t564;
t456 = -t513 * pkin(4) - t540 * pkin(10) - t591 + (t529 - t620 + t622) * t542;
t476 = qJD(5) * t518 + t513 * t576 + t514 * t581;
t447 = (-t518 * t560 - t476) * pkin(11) + (t519 * t560 - t475) * pkin(5) + t456;
t575 = sin(qJ(6));
t580 = cos(qJ(6));
t443 = -t446 * t575 + t447 * t580;
t496 = -t519 * t575 + t560 * t580;
t463 = qJD(6) * t496 + t476 * t580 + t548 * t575;
t474 = qJDD(6) - t475;
t497 = t519 * t580 + t560 * t575;
t478 = -mrSges(7,1) * t496 + mrSges(7,2) * t497;
t517 = qJD(6) - t518;
t479 = -mrSges(7,2) * t517 + mrSges(7,3) * t496;
t439 = m(7) * t443 + mrSges(7,1) * t474 - mrSges(7,3) * t463 - t478 * t497 + t479 * t517;
t444 = t446 * t580 + t447 * t575;
t462 = -qJD(6) * t497 - t476 * t575 + t548 * t580;
t480 = mrSges(7,1) * t517 - mrSges(7,3) * t497;
t440 = m(7) * t444 - mrSges(7,2) * t474 + mrSges(7,3) * t462 + t478 * t496 - t480 * t517;
t595 = -t439 * t575 + t580 * t440;
t427 = m(6) * t450 - mrSges(6,2) * t548 + mrSges(6,3) * t475 + t486 * t518 - t499 * t560 + t595;
t449 = t452 * t581 - t455 * t576;
t498 = -mrSges(6,2) * t560 + mrSges(6,3) * t518;
t445 = -pkin(5) * t548 - pkin(11) * t559 + t487 * t519 - t449;
t590 = -m(7) * t445 + t462 * mrSges(7,1) - mrSges(7,2) * t463 + t496 * t479 - t480 * t497;
t435 = m(6) * t449 + mrSges(6,1) * t548 - mrSges(6,3) * t476 - t486 * t519 + t498 * t560 + t590;
t423 = t576 * t427 + t581 * t435;
t528 = -mrSges(5,2) * t541 + mrSges(5,3) * t564;
t589 = -m(5) * t459 + t549 * mrSges(5,1) + t564 * t528 - t423;
t422 = t514 * mrSges(5,2) + t542 * t521 - t589;
t466 = Ifges(7,5) * t497 + Ifges(7,6) * t496 + Ifges(7,3) * t517;
t468 = Ifges(7,1) * t497 + Ifges(7,4) * t496 + Ifges(7,5) * t517;
t433 = -mrSges(7,1) * t445 + mrSges(7,3) * t444 + Ifges(7,4) * t463 + Ifges(7,2) * t462 + Ifges(7,6) * t474 - t466 * t497 + t468 * t517;
t467 = Ifges(7,4) * t497 + Ifges(7,2) * t496 + Ifges(7,6) * t517;
t434 = mrSges(7,2) * t445 - mrSges(7,3) * t443 + Ifges(7,1) * t463 + Ifges(7,4) * t462 + Ifges(7,5) * t474 + t466 * t496 - t467 * t517;
t482 = Ifges(6,4) * t519 + Ifges(6,2) * t518 + Ifges(6,6) * t560;
t483 = Ifges(6,1) * t519 + Ifges(6,4) * t518 + Ifges(6,5) * t560;
t588 = -mrSges(6,1) * t449 + mrSges(6,2) * t450 - Ifges(6,5) * t476 - Ifges(6,6) * t475 - Ifges(6,3) * t548 - pkin(5) * t590 - pkin(11) * t595 - t580 * t433 - t575 * t434 - t519 * t482 + t518 * t483;
t527 = -mrSges(5,1) * t564 + mrSges(5,2) * t542;
t596 = t581 * t427 - t576 * t435;
t593 = m(5) * t458 + t549 * mrSges(5,3) + t564 * t527 + t596;
t606 = -t616 * t541 + t626 * t542 + t615 * t564;
t608 = t625 * t541 - t616 * t542 + t614 * t564;
t623 = t513 * t614 + t514 * t615 + t541 * t606 - t542 * t608 + t624 * t549 + mrSges(4,1) * t464 - mrSges(5,1) * t459 - mrSges(4,2) * t465 + mrSges(5,3) * t458 - pkin(3) * t422 - pkin(4) * t423 + qJ(4) * (-t513 * mrSges(5,2) - t541 * t521 + t593) + t588;
t617 = -mrSges(4,3) - mrSges(5,2);
t610 = t573 * t578;
t526 = mrSges(4,1) * t564 - mrSges(4,3) * t542;
t605 = -mrSges(4,1) * t541 - mrSges(4,2) * t542 - t521;
t418 = m(4) * t465 - t549 * mrSges(4,2) + t513 * t617 - t564 * t526 + t541 * t605 + t593;
t525 = -mrSges(4,2) * t564 - mrSges(4,3) * t541;
t420 = m(4) * t464 + t549 * mrSges(4,1) + t514 * t617 + t564 * t525 + t542 * t605 + t589;
t413 = t577 * t418 + t621 * t420;
t429 = t580 * t439 + t575 * t440;
t607 = -t614 * t541 - t615 * t542 - t624 * t564;
t597 = t621 * t418 - t420 * t577;
t592 = -m(6) * t456 + t475 * mrSges(6,1) - t476 * mrSges(6,2) + t518 * t498 - t519 * t499 - t429;
t460 = (-(2 * qJD(4)) + t620) * t542 + t591;
t425 = m(5) * t460 + t513 * mrSges(5,1) - t514 * mrSges(5,3) - t542 * t527 + t541 * t528 + t592;
t587 = mrSges(7,1) * t443 - mrSges(7,2) * t444 + Ifges(7,5) * t463 + Ifges(7,6) * t462 + Ifges(7,3) * t474 + t467 * t497 - t468 * t496;
t586 = -m(4) * t491 - t513 * mrSges(4,1) - t514 * mrSges(4,2) - t541 * t525 - t542 * t526 - t425;
t554 = (-mrSges(3,1) * t582 + mrSges(3,2) * t578) * t603;
t551 = -mrSges(3,2) * t570 + mrSges(3,3) * t599;
t550 = mrSges(3,1) * t570 - mrSges(3,3) * t600;
t534 = -t573 * t552 - t618;
t533 = Ifges(3,5) * t570 + (Ifges(3,1) * t578 + Ifges(3,4) * t582) * t603;
t532 = Ifges(3,6) * t570 + (Ifges(3,4) * t578 + Ifges(3,2) * t582) * t603;
t531 = Ifges(3,3) * t570 + (Ifges(3,5) * t578 + Ifges(3,6) * t582) * t603;
t516 = -g(3) * t610 + t604;
t481 = Ifges(6,5) * t519 + Ifges(6,6) * t518 + Ifges(6,3) * t560;
t424 = m(3) * t515 + t569 * mrSges(3,1) - t556 * mrSges(3,3) + t570 * t551 - t554 * t600 + t586;
t415 = -mrSges(6,1) * t456 + mrSges(6,3) * t450 + Ifges(6,4) * t476 + Ifges(6,2) * t475 + Ifges(6,6) * t548 - pkin(5) * t429 - t481 * t519 + t483 * t560 - t587;
t414 = mrSges(6,2) * t456 - mrSges(6,3) * t449 + Ifges(6,1) * t476 + Ifges(6,4) * t475 + Ifges(6,5) * t548 - pkin(11) * t429 - t433 * t575 + t434 * t580 + t481 * t518 - t482 * t560;
t412 = m(3) * t516 - mrSges(3,2) * t569 + mrSges(3,3) * t557 - t550 * t570 + t554 * t599 + t597;
t411 = mrSges(4,2) * t491 + mrSges(5,2) * t459 - mrSges(4,3) * t464 - mrSges(5,3) * t460 - pkin(10) * t423 - qJ(4) * t425 + t581 * t414 - t576 * t415 - t616 * t513 + t626 * t514 + t607 * t541 + t615 * t549 + t608 * t564;
t410 = -mrSges(4,1) * t491 - mrSges(5,1) * t460 + mrSges(5,2) * t458 + mrSges(4,3) * t465 - pkin(3) * t425 - pkin(4) * t592 - pkin(10) * t596 - t576 * t414 - t581 * t415 - t625 * t513 + t616 * t514 + t607 * t542 - t614 * t549 + t606 * t564;
t409 = Ifges(3,5) * t556 + Ifges(3,6) * t557 + Ifges(3,3) * t569 + mrSges(3,1) * t515 - mrSges(3,2) * t516 + t577 * t411 + t621 * t410 + pkin(2) * t586 + pkin(9) * t597 + (t532 * t578 - t533 * t582) * t603;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t598 - mrSges(2,2) * t594 + (mrSges(3,2) * t534 - mrSges(3,3) * t515 + Ifges(3,1) * t556 + Ifges(3,4) * t557 + Ifges(3,5) * t569 - pkin(9) * t413 - t577 * t410 + t411 * t621 + t531 * t599 - t570 * t532) * t610 + (-mrSges(3,1) * t534 + mrSges(3,3) * t516 + Ifges(3,4) * t556 + Ifges(3,2) * t557 + Ifges(3,6) * t569 - pkin(2) * t413 - t531 * t600 + t570 * t533 - t623) * t609 + t574 * t409 + pkin(1) * ((t412 * t578 + t424 * t582) * t574 + (-m(3) * t534 + t557 * mrSges(3,1) - t556 * mrSges(3,2) + (-t550 * t578 + t551 * t582) * t603 - t413) * t573) + (t412 * t582 - t424 * t578) * t619; t409; t623; t422; -t588; t587;];
tauJ  = t1;
