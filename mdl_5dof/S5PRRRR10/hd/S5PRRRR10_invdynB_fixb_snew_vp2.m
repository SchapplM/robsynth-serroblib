% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR10_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR10_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:24:11
% EndTime: 2019-12-05 17:24:24
% DurationCPUTime: 11.07s
% Computational Cost: add. (175909->279), mult. (361309->377), div. (0->0), fcn. (282841->14), ass. (0->128)
t535 = sin(pkin(6));
t542 = sin(qJ(3));
t546 = cos(qJ(3));
t560 = qJD(2) * qJD(3);
t520 = (-qJDD(2) * t546 + t542 * t560) * t535;
t534 = sin(pkin(11));
t537 = cos(pkin(11));
t526 = t534 * g(1) - t537 * g(2);
t527 = -t537 * g(1) - t534 * g(2);
t533 = -g(3) + qJDD(1);
t543 = sin(qJ(2));
t539 = cos(pkin(5));
t547 = cos(qJ(2));
t563 = t539 * t547;
t536 = sin(pkin(5));
t566 = t536 * t547;
t498 = t526 * t563 - t543 * t527 + t533 * t566;
t548 = qJD(2) ^ 2;
t570 = pkin(8) * t535;
t494 = qJDD(2) * pkin(2) + t548 * t570 + t498;
t564 = t539 * t543;
t567 = t536 * t543;
t499 = t526 * t564 + t547 * t527 + t533 * t567;
t495 = -t548 * pkin(2) + qJDD(2) * t570 + t499;
t513 = -t536 * t526 + t539 * t533;
t538 = cos(pkin(6));
t470 = -t542 * t495 + (t494 * t538 + t513 * t535) * t546;
t532 = t538 * qJD(2) + qJD(3);
t561 = qJD(2) * t535;
t558 = t546 * t561;
t516 = -t532 * mrSges(4,2) + mrSges(4,3) * t558;
t517 = (-mrSges(4,1) * t546 + mrSges(4,2) * t542) * t561;
t519 = (qJDD(2) * t542 + t546 * t560) * t535;
t531 = t538 * qJDD(2) + qJDD(3);
t565 = t538 * t542;
t568 = t535 * t542;
t471 = t494 * t565 + t546 * t495 + t513 * t568;
t518 = (-pkin(3) * t546 - pkin(9) * t542) * t561;
t530 = t532 ^ 2;
t467 = -t530 * pkin(3) + t531 * pkin(9) + t518 * t558 + t471;
t509 = t538 * t513;
t469 = t520 * pkin(3) - t519 * pkin(9) + t509 + (-t494 + (pkin(3) * t542 - pkin(9) * t546) * t532 * qJD(2)) * t535;
t541 = sin(qJ(4));
t545 = cos(qJ(4));
t464 = t545 * t467 + t541 * t469;
t559 = t542 * t561;
t511 = t545 * t532 - t541 * t559;
t512 = t541 * t532 + t545 * t559;
t497 = -t511 * pkin(4) - t512 * pkin(10);
t514 = qJDD(4) + t520;
t525 = qJD(4) - t558;
t524 = t525 ^ 2;
t461 = -t524 * pkin(4) + t514 * pkin(10) + t511 * t497 + t464;
t466 = -t531 * pkin(3) - t530 * pkin(9) + t518 * t559 - t470;
t489 = -t512 * qJD(4) - t541 * t519 + t545 * t531;
t490 = t511 * qJD(4) + t545 * t519 + t541 * t531;
t462 = (-t511 * t525 - t490) * pkin(10) + (t512 * t525 - t489) * pkin(4) + t466;
t540 = sin(qJ(5));
t544 = cos(qJ(5));
t458 = -t540 * t461 + t544 * t462;
t500 = -t540 * t512 + t544 * t525;
t474 = t500 * qJD(5) + t544 * t490 + t540 * t514;
t501 = t544 * t512 + t540 * t525;
t479 = -t500 * mrSges(6,1) + t501 * mrSges(6,2);
t510 = qJD(5) - t511;
t481 = -t510 * mrSges(6,2) + t500 * mrSges(6,3);
t487 = qJDD(5) - t489;
t456 = m(6) * t458 + t487 * mrSges(6,1) - t474 * mrSges(6,3) - t501 * t479 + t510 * t481;
t459 = t544 * t461 + t540 * t462;
t473 = -t501 * qJD(5) - t540 * t490 + t544 * t514;
t482 = t510 * mrSges(6,1) - t501 * mrSges(6,3);
t457 = m(6) * t459 - t487 * mrSges(6,2) + t473 * mrSges(6,3) + t500 * t479 - t510 * t482;
t450 = t544 * t456 + t540 * t457;
t502 = -t525 * mrSges(5,2) + t511 * mrSges(5,3);
t503 = t525 * mrSges(5,1) - t512 * mrSges(5,3);
t549 = -m(5) * t466 + t489 * mrSges(5,1) - t490 * mrSges(5,2) + t511 * t502 - t512 * t503 - t450;
t446 = m(4) * t470 + t531 * mrSges(4,1) - t519 * mrSges(4,3) + t532 * t516 - t517 * t559 + t549;
t569 = t446 * t546;
t515 = t532 * mrSges(4,1) - mrSges(4,3) * t559;
t496 = -t511 * mrSges(5,1) + t512 * mrSges(5,2);
t555 = -t540 * t456 + t544 * t457;
t449 = m(5) * t464 - t514 * mrSges(5,2) + t489 * mrSges(5,3) + t511 * t496 - t525 * t503 + t555;
t463 = -t541 * t467 + t545 * t469;
t460 = -t514 * pkin(4) - t524 * pkin(10) + t512 * t497 - t463;
t550 = -m(6) * t460 + t473 * mrSges(6,1) - t474 * mrSges(6,2) + t500 * t481 - t501 * t482;
t454 = m(5) * t463 + t514 * mrSges(5,1) - t490 * mrSges(5,3) - t512 * t496 + t525 * t502 + t550;
t556 = t545 * t449 - t541 * t454;
t440 = m(4) * t471 - t531 * mrSges(4,2) - t520 * mrSges(4,3) - t532 * t515 + t517 * t558 + t556;
t443 = t541 * t449 + t545 * t454;
t480 = -t535 * t494 + t509;
t442 = m(4) * t480 + t520 * mrSges(4,1) + t519 * mrSges(4,2) + (t515 * t542 - t516 * t546) * t561 + t443;
t429 = t440 * t565 - t535 * t442 + t538 * t569;
t425 = m(3) * t498 + qJDD(2) * mrSges(3,1) - t548 * mrSges(3,2) + t429;
t428 = t440 * t568 + t538 * t442 + t535 * t569;
t427 = m(3) * t513 + t428;
t434 = t546 * t440 - t542 * t446;
t433 = m(3) * t499 - t548 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t434;
t415 = t425 * t563 - t536 * t427 + t433 * t564;
t413 = m(2) * t526 + t415;
t419 = -t543 * t425 + t547 * t433;
t418 = m(2) * t527 + t419;
t562 = t537 * t413 + t534 * t418;
t414 = t425 * t566 + t539 * t427 + t433 * t567;
t557 = -t534 * t413 + t537 * t418;
t475 = Ifges(6,5) * t501 + Ifges(6,6) * t500 + Ifges(6,3) * t510;
t477 = Ifges(6,1) * t501 + Ifges(6,4) * t500 + Ifges(6,5) * t510;
t451 = -mrSges(6,1) * t460 + mrSges(6,3) * t459 + Ifges(6,4) * t474 + Ifges(6,2) * t473 + Ifges(6,6) * t487 - t501 * t475 + t510 * t477;
t476 = Ifges(6,4) * t501 + Ifges(6,2) * t500 + Ifges(6,6) * t510;
t452 = mrSges(6,2) * t460 - mrSges(6,3) * t458 + Ifges(6,1) * t474 + Ifges(6,4) * t473 + Ifges(6,5) * t487 + t500 * t475 - t510 * t476;
t483 = Ifges(5,5) * t512 + Ifges(5,6) * t511 + Ifges(5,3) * t525;
t484 = Ifges(5,4) * t512 + Ifges(5,2) * t511 + Ifges(5,6) * t525;
t435 = mrSges(5,2) * t466 - mrSges(5,3) * t463 + Ifges(5,1) * t490 + Ifges(5,4) * t489 + Ifges(5,5) * t514 - pkin(10) * t450 - t540 * t451 + t544 * t452 + t511 * t483 - t525 * t484;
t485 = Ifges(5,1) * t512 + Ifges(5,4) * t511 + Ifges(5,5) * t525;
t436 = -mrSges(5,1) * t466 - mrSges(6,1) * t458 + mrSges(6,2) * t459 + mrSges(5,3) * t464 + Ifges(5,4) * t490 - Ifges(6,5) * t474 + Ifges(5,2) * t489 + Ifges(5,6) * t514 - Ifges(6,6) * t473 - Ifges(6,3) * t487 - pkin(4) * t450 - t501 * t476 + t500 * t477 - t512 * t483 + t525 * t485;
t506 = Ifges(4,6) * t532 + (Ifges(4,4) * t542 + Ifges(4,2) * t546) * t561;
t507 = Ifges(4,5) * t532 + (Ifges(4,1) * t542 + Ifges(4,4) * t546) * t561;
t420 = Ifges(4,5) * t519 - Ifges(4,6) * t520 + Ifges(4,3) * t531 + mrSges(4,1) * t470 - mrSges(4,2) * t471 + t541 * t435 + t545 * t436 + pkin(3) * t549 + pkin(9) * t556 + (t506 * t542 - t507 * t546) * t561;
t505 = Ifges(4,3) * t532 + (Ifges(4,5) * t542 + Ifges(4,6) * t546) * t561;
t421 = mrSges(4,2) * t480 - mrSges(4,3) * t470 + Ifges(4,1) * t519 - Ifges(4,4) * t520 + Ifges(4,5) * t531 - pkin(9) * t443 + t545 * t435 - t541 * t436 + t505 * t558 - t532 * t506;
t422 = Ifges(4,4) * t519 - Ifges(4,2) * t520 + Ifges(4,6) * t531 - t505 * t559 + t532 * t507 - mrSges(4,1) * t480 + mrSges(4,3) * t471 - Ifges(5,5) * t490 - Ifges(5,6) * t489 - Ifges(5,3) * t514 - t512 * t484 + t511 * t485 - mrSges(5,1) * t463 + mrSges(5,2) * t464 - t540 * t452 - t544 * t451 - pkin(4) * t550 - pkin(10) * t555 - pkin(3) * t443;
t551 = pkin(8) * t434 + t421 * t542 + t422 * t546;
t410 = -mrSges(3,1) * t513 + mrSges(3,3) * t499 + t548 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t428 - t535 * t420 + t551 * t538;
t411 = mrSges(3,2) * t513 - mrSges(3,3) * t498 + Ifges(3,5) * qJDD(2) - t548 * Ifges(3,6) + t546 * t421 - t542 * t422 + (-t428 * t535 - t429 * t538) * pkin(8);
t552 = pkin(7) * t419 + t410 * t547 + t411 * t543;
t409 = mrSges(3,1) * t498 - mrSges(3,2) * t499 + Ifges(3,3) * qJDD(2) + pkin(2) * t429 + t538 * t420 + t551 * t535;
t408 = mrSges(2,2) * t533 - mrSges(2,3) * t526 - t543 * t410 + t547 * t411 + (-t414 * t536 - t415 * t539) * pkin(7);
t407 = -mrSges(2,1) * t533 + mrSges(2,3) * t527 - pkin(1) * t414 - t536 * t409 + t552 * t539;
t1 = [-m(1) * g(1) + t557; -m(1) * g(2) + t562; -m(1) * g(3) + m(2) * t533 + t414; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t562 - t534 * t407 + t537 * t408; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t557 + t537 * t407 + t534 * t408; -mrSges(1,1) * g(2) + mrSges(2,1) * t526 + mrSges(1,2) * g(1) - mrSges(2,2) * t527 + pkin(1) * t415 + t539 * t409 + t552 * t536;];
tauB = t1;
