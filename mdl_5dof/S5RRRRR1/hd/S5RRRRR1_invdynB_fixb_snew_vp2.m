% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% Datum: 2019-12-05 18:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:50:13
% EndTime: 2019-12-05 18:50:21
% DurationCPUTime: 5.45s
% Computational Cost: add. (55552->294), mult. (121727->377), div. (0->0), fcn. (92430->10), ass. (0->114)
t535 = sin(qJ(1));
t540 = cos(qJ(1));
t521 = -t540 * g(1) - t535 * g(2);
t541 = qJD(1) ^ 2;
t515 = -t541 * pkin(1) + t521;
t534 = sin(qJ(2));
t539 = cos(qJ(2));
t502 = t539 * g(3) - t534 * t515;
t514 = (mrSges(3,1) * t539 - mrSges(3,2) * t534) * qJD(1);
t550 = qJD(1) * qJD(2);
t516 = -t534 * qJDD(1) - t539 * t550;
t551 = qJD(1) * t539;
t519 = -qJD(2) * mrSges(3,2) - mrSges(3,3) * t551;
t552 = qJD(1) * t534;
t497 = (t534 * t539 * t541 + qJDD(2)) * pkin(2) + t502;
t503 = t534 * g(3) + t539 * t515;
t498 = (-t539 ^ 2 * t541 - qJD(2) ^ 2) * pkin(2) + t503;
t533 = sin(qJ(3));
t538 = cos(qJ(3));
t475 = t538 * t497 - t533 * t498;
t507 = (t533 * t534 - t538 * t539) * qJD(1);
t549 = t534 * t550;
t517 = -t539 * qJDD(1) + t549;
t481 = t507 * qJD(3) + t538 * t516 + t533 * t517;
t508 = (-t533 * t539 - t534 * t538) * qJD(1);
t491 = -t507 * mrSges(4,1) + t508 * mrSges(4,2);
t529 = qJD(2) + qJD(3);
t499 = -t529 * mrSges(4,2) + t507 * mrSges(4,3);
t528 = qJDD(2) + qJDD(3);
t461 = (t507 * t508 + t528) * pkin(3) + t475;
t476 = t533 * t497 + t538 * t498;
t468 = (-t507 ^ 2 - t529 ^ 2) * pkin(3) + t476;
t532 = sin(qJ(4));
t537 = cos(qJ(4));
t451 = t532 * t461 + t537 * t468;
t480 = -t508 * qJD(3) - t533 * t516 + t538 * t517;
t490 = t532 * t507 + t537 * t508;
t458 = -t490 * qJD(4) + t537 * t480 - t532 * t481;
t489 = t537 * t507 - t532 * t508;
t473 = -t489 * mrSges(5,1) + t490 * mrSges(5,2);
t524 = qJD(4) + t529;
t483 = t524 * mrSges(5,1) - t490 * mrSges(5,3);
t523 = qJDD(4) + t528;
t459 = t489 * qJD(4) + t532 * t480 + t537 * t481;
t520 = t535 * g(1) - t540 * g(2);
t513 = qJDD(1) * pkin(1) + t520;
t495 = (-t517 - t549) * pkin(2) + t513;
t464 = t495 + (t508 * t529 - t480) * pkin(3);
t444 = (-t489 * t524 - t459) * pkin(6) + (t490 * t524 - t458) * pkin(4) + t464;
t474 = -t489 * pkin(4) - t490 * pkin(6);
t522 = t524 ^ 2;
t446 = -t522 * pkin(4) + t523 * pkin(6) + t489 * t474 + t451;
t531 = sin(qJ(5));
t536 = cos(qJ(5));
t442 = t536 * t444 - t531 * t446;
t477 = -t531 * t490 + t536 * t524;
t449 = t477 * qJD(5) + t536 * t459 + t531 * t523;
t456 = qJDD(5) - t458;
t478 = t536 * t490 + t531 * t524;
t462 = -t477 * mrSges(6,1) + t478 * mrSges(6,2);
t487 = qJD(5) - t489;
t465 = -t487 * mrSges(6,2) + t477 * mrSges(6,3);
t440 = m(6) * t442 + t456 * mrSges(6,1) - t449 * mrSges(6,3) - t478 * t462 + t487 * t465;
t443 = t531 * t444 + t536 * t446;
t448 = -t478 * qJD(5) - t531 * t459 + t536 * t523;
t466 = t487 * mrSges(6,1) - t478 * mrSges(6,3);
t441 = m(6) * t443 - t456 * mrSges(6,2) + t448 * mrSges(6,3) + t477 * t462 - t487 * t466;
t547 = -t531 * t440 + t536 * t441;
t431 = m(5) * t451 - t523 * mrSges(5,2) + t458 * mrSges(5,3) + t489 * t473 - t524 * t483 + t547;
t450 = t537 * t461 - t532 * t468;
t482 = -t524 * mrSges(5,2) + t489 * mrSges(5,3);
t445 = -t523 * pkin(4) - t522 * pkin(6) + t490 * t474 - t450;
t545 = -m(6) * t445 + t448 * mrSges(6,1) - t449 * mrSges(6,2) + t477 * t465 - t478 * t466;
t436 = m(5) * t450 + t523 * mrSges(5,1) - t459 * mrSges(5,3) - t490 * t473 + t524 * t482 + t545;
t553 = t532 * t431 + t537 * t436;
t426 = m(4) * t475 + t528 * mrSges(4,1) - t481 * mrSges(4,3) - t508 * t491 + t529 * t499 + t553;
t500 = t529 * mrSges(4,1) - t508 * mrSges(4,3);
t427 = m(4) * t476 - t528 * mrSges(4,2) + t480 * mrSges(4,3) + t537 * t431 - t532 * t436 + t507 * t491 - t529 * t500;
t554 = t538 * t426 + t533 * t427;
t420 = m(3) * t502 + qJDD(2) * mrSges(3,1) - t516 * mrSges(3,3) + qJD(2) * t519 + t514 * t552 + t554;
t518 = qJD(2) * mrSges(3,1) + mrSges(3,3) * t552;
t421 = m(3) * t503 - qJDD(2) * mrSges(3,2) + t517 * mrSges(3,3) - qJD(2) * t518 - t533 * t426 + t538 * t427 - t514 * t551;
t417 = m(2) * t521 - t541 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t534 * t420 + t539 * t421;
t432 = t536 * t440 + t531 * t441;
t544 = m(5) * t464 - t458 * mrSges(5,1) + t459 * mrSges(5,2) - t489 * t482 + t490 * t483 + t432;
t543 = m(4) * t495 - t480 * mrSges(4,1) + t481 * mrSges(4,2) - t507 * t499 + t508 * t500 + t544;
t542 = m(3) * t513 - t517 * mrSges(3,1) + t516 * mrSges(3,2) - t518 * t552 + t519 * t551 + t543;
t429 = m(2) * t520 + qJDD(1) * mrSges(2,1) - t541 * mrSges(2,2) + t542;
t555 = t535 * t417 + t540 * t429;
t548 = t540 * t417 - t429 * t535;
t546 = -t539 * t420 - t534 * t421;
t506 = Ifges(3,5) * qJD(2) + (-Ifges(3,1) * t534 - Ifges(3,4) * t539) * qJD(1);
t505 = Ifges(3,6) * qJD(2) + (-Ifges(3,4) * t534 - Ifges(3,2) * t539) * qJD(1);
t504 = Ifges(3,3) * qJD(2) + (-Ifges(3,5) * t534 - Ifges(3,6) * t539) * qJD(1);
t486 = Ifges(4,1) * t508 + Ifges(4,4) * t507 + Ifges(4,5) * t529;
t485 = Ifges(4,4) * t508 + Ifges(4,2) * t507 + Ifges(4,6) * t529;
t484 = Ifges(4,5) * t508 + Ifges(4,6) * t507 + Ifges(4,3) * t529;
t471 = Ifges(5,1) * t490 + Ifges(5,4) * t489 + Ifges(5,5) * t524;
t470 = Ifges(5,4) * t490 + Ifges(5,2) * t489 + Ifges(5,6) * t524;
t469 = Ifges(5,5) * t490 + Ifges(5,6) * t489 + Ifges(5,3) * t524;
t454 = Ifges(6,1) * t478 + Ifges(6,4) * t477 + Ifges(6,5) * t487;
t453 = Ifges(6,4) * t478 + Ifges(6,2) * t477 + Ifges(6,6) * t487;
t452 = Ifges(6,5) * t478 + Ifges(6,6) * t477 + Ifges(6,3) * t487;
t434 = mrSges(6,2) * t445 - mrSges(6,3) * t442 + Ifges(6,1) * t449 + Ifges(6,4) * t448 + Ifges(6,5) * t456 + t477 * t452 - t487 * t453;
t433 = -mrSges(6,1) * t445 + mrSges(6,3) * t443 + Ifges(6,4) * t449 + Ifges(6,2) * t448 + Ifges(6,6) * t456 - t478 * t452 + t487 * t454;
t423 = -mrSges(5,1) * t464 - mrSges(6,1) * t442 + mrSges(6,2) * t443 + mrSges(5,3) * t451 + Ifges(5,4) * t459 - Ifges(6,5) * t449 + Ifges(5,2) * t458 + Ifges(5,6) * t523 - Ifges(6,6) * t448 - Ifges(6,3) * t456 - pkin(4) * t432 - t478 * t453 + t477 * t454 - t490 * t469 + t524 * t471;
t422 = mrSges(5,2) * t464 - mrSges(5,3) * t450 + Ifges(5,1) * t459 + Ifges(5,4) * t458 + Ifges(5,5) * t523 - pkin(6) * t432 - t531 * t433 + t536 * t434 + t489 * t469 - t524 * t470;
t419 = mrSges(4,2) * t495 - mrSges(4,3) * t475 + Ifges(4,1) * t481 + Ifges(4,4) * t480 + Ifges(4,5) * t528 + t537 * t422 - t532 * t423 + t507 * t484 - t529 * t485;
t418 = -mrSges(4,1) * t495 + mrSges(4,3) * t476 + Ifges(4,4) * t481 + Ifges(4,2) * t480 + Ifges(4,6) * t528 - pkin(3) * t544 + t532 * t422 + t537 * t423 - t508 * t484 + t529 * t486;
t414 = mrSges(3,2) * t513 - mrSges(3,3) * t502 + Ifges(3,1) * t516 + Ifges(3,4) * t517 + Ifges(3,5) * qJDD(2) - qJD(2) * t505 - t533 * t418 + t538 * t419 - t504 * t551;
t413 = -mrSges(3,1) * t513 + mrSges(3,3) * t503 + Ifges(3,4) * t516 + Ifges(3,2) * t517 + Ifges(3,6) * qJDD(2) - pkin(2) * t543 + qJD(2) * t506 + t538 * t418 + t533 * t419 + t504 * t552;
t412 = mrSges(2,3) * t521 + Ifges(5,3) * t523 + Ifges(4,6) * t480 + Ifges(4,5) * t481 - t507 * t486 + t508 * t485 + Ifges(3,5) * t516 + Ifges(3,6) * t517 + (-t534 * t505 + t539 * t506) * qJD(1) + pkin(4) * t545 + t536 * t433 + t541 * Ifges(2,5) + mrSges(5,1) * t450 - mrSges(5,2) * t451 + pkin(6) * t547 + Ifges(2,6) * qJDD(1) + Ifges(3,3) * qJDD(2) + Ifges(5,6) * t458 + Ifges(5,5) * t459 + mrSges(2,1) * g(3) + Ifges(4,3) * t528 + t531 * t434 - pkin(1) * t546 + pkin(3) * t553 + pkin(2) * t554 + mrSges(4,1) * t475 - mrSges(4,2) * t476 - t489 * t471 + t490 * t470 + mrSges(3,1) * t502 - mrSges(3,2) * t503;
t411 = -mrSges(2,2) * g(3) - mrSges(2,3) * t520 + Ifges(2,5) * qJDD(1) - t541 * Ifges(2,6) - t534 * t413 + t539 * t414;
t1 = [-m(1) * g(1) + t548; -m(1) * g(2) + t555; (-m(1) - m(2)) * g(3) + t546; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t555 + t540 * t411 - t535 * t412; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t548 + t535 * t411 + t540 * t412; -mrSges(1,1) * g(2) + mrSges(2,1) * t520 + mrSges(1,2) * g(1) - mrSges(2,2) * t521 + Ifges(2,3) * qJDD(1) + pkin(1) * t542 - t539 * t413 - t534 * t414;];
tauB = t1;
