% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:22
% EndTime: 2019-12-05 19:00:28
% DurationCPUTime: 5.26s
% Computational Cost: add. (93727->273), mult. (120226->346), div. (0->0), fcn. (77730->10), ass. (0->110)
t521 = qJD(1) + qJD(2);
t517 = t521 ^ 2;
t548 = pkin(3) * t517;
t525 = sin(qJ(3));
t547 = t521 * t525;
t530 = cos(qJ(3));
t546 = t521 * t530;
t527 = sin(qJ(1));
t532 = cos(qJ(1));
t511 = t532 * g(2) + t527 * g(3);
t505 = qJDD(1) * pkin(1) + t511;
t510 = t527 * g(2) - t532 * g(3);
t533 = qJD(1) ^ 2;
t506 = -t533 * pkin(1) + t510;
t526 = sin(qJ(2));
t531 = cos(qJ(2));
t486 = t526 * t505 + t531 * t506;
t519 = qJDD(1) + qJDD(2);
t484 = -t517 * pkin(2) + t519 * pkin(7) + t486;
t545 = t525 * t484;
t544 = qJD(3) * t521;
t500 = t525 * t519 + t530 * t544;
t465 = qJDD(3) * pkin(3) - t500 * pkin(8) - t545 + (pkin(8) * t544 + t525 * t548 - g(1)) * t530;
t474 = -t525 * g(1) + t530 * t484;
t501 = t530 * t519 - t525 * t544;
t509 = qJD(3) * pkin(3) - pkin(8) * t547;
t522 = t530 ^ 2;
t466 = t501 * pkin(8) - qJD(3) * t509 - t522 * t548 + t474;
t524 = sin(qJ(4));
t529 = cos(qJ(4));
t450 = t529 * t465 - t524 * t466;
t494 = (-t524 * t525 + t529 * t530) * t521;
t470 = t494 * qJD(4) + t529 * t500 + t524 * t501;
t495 = (t524 * t530 + t525 * t529) * t521;
t518 = qJDD(3) + qJDD(4);
t520 = qJD(3) + qJD(4);
t446 = (t494 * t520 - t470) * pkin(9) + (t494 * t495 + t518) * pkin(4) + t450;
t451 = t524 * t465 + t529 * t466;
t469 = -t495 * qJD(4) - t524 * t500 + t529 * t501;
t489 = t520 * pkin(4) - t495 * pkin(9);
t490 = t494 ^ 2;
t447 = -t490 * pkin(4) + t469 * pkin(9) - t520 * t489 + t451;
t523 = sin(qJ(5));
t528 = cos(qJ(5));
t444 = t528 * t446 - t523 * t447;
t479 = t528 * t494 - t523 * t495;
t455 = t479 * qJD(5) + t523 * t469 + t528 * t470;
t480 = t523 * t494 + t528 * t495;
t461 = -t479 * mrSges(6,1) + t480 * mrSges(6,2);
t513 = qJD(5) + t520;
t471 = -t513 * mrSges(6,2) + t479 * mrSges(6,3);
t512 = qJDD(5) + t518;
t442 = m(6) * t444 + t512 * mrSges(6,1) - t455 * mrSges(6,3) - t480 * t461 + t513 * t471;
t445 = t523 * t446 + t528 * t447;
t454 = -t480 * qJD(5) + t528 * t469 - t523 * t470;
t472 = t513 * mrSges(6,1) - t480 * mrSges(6,3);
t443 = m(6) * t445 - t512 * mrSges(6,2) + t454 * mrSges(6,3) + t479 * t461 - t513 * t472;
t434 = t528 * t442 + t523 * t443;
t482 = -t494 * mrSges(5,1) + t495 * mrSges(5,2);
t487 = -t520 * mrSges(5,2) + t494 * mrSges(5,3);
t432 = m(5) * t450 + t518 * mrSges(5,1) - t470 * mrSges(5,3) - t495 * t482 + t520 * t487 + t434;
t488 = t520 * mrSges(5,1) - t495 * mrSges(5,3);
t539 = -t523 * t442 + t528 * t443;
t433 = m(5) * t451 - t518 * mrSges(5,2) + t469 * mrSges(5,3) + t494 * t482 - t520 * t488 + t539;
t428 = t529 * t432 + t524 * t433;
t473 = -t530 * g(1) - t545;
t499 = (-mrSges(4,1) * t530 + mrSges(4,2) * t525) * t521;
t508 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t546;
t426 = m(4) * t473 + qJDD(3) * mrSges(4,1) - t500 * mrSges(4,3) + qJD(3) * t508 - t499 * t547 + t428;
t507 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t547;
t540 = -t524 * t432 + t529 * t433;
t427 = m(4) * t474 - qJDD(3) * mrSges(4,2) + t501 * mrSges(4,3) - qJD(3) * t507 + t499 * t546 + t540;
t541 = -t525 * t426 + t530 * t427;
t419 = m(3) * t486 - t517 * mrSges(3,1) - t519 * mrSges(3,2) + t541;
t485 = t531 * t505 - t526 * t506;
t536 = -t519 * pkin(2) - t485;
t483 = -t517 * pkin(7) + t536;
t467 = -t501 * pkin(3) + t509 * t547 + (-pkin(8) * t522 - pkin(7)) * t517 + t536;
t449 = -t469 * pkin(4) - t490 * pkin(9) + t495 * t489 + t467;
t538 = m(6) * t449 - t454 * mrSges(6,1) + t455 * mrSges(6,2) - t479 * t471 + t480 * t472;
t535 = m(5) * t467 - t469 * mrSges(5,1) + t470 * mrSges(5,2) - t494 * t487 + t495 * t488 + t538;
t534 = -m(4) * t483 + t501 * mrSges(4,1) - t500 * mrSges(4,2) - t507 * t547 + t508 * t546 - t535;
t438 = m(3) * t485 + t519 * mrSges(3,1) - t517 * mrSges(3,2) + t534;
t416 = t526 * t419 + t531 * t438;
t420 = t530 * t426 + t525 * t427;
t542 = t531 * t419 - t526 * t438;
t414 = m(2) * t510 - t533 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t542;
t415 = m(2) * t511 + qJDD(1) * mrSges(2,1) - t533 * mrSges(2,2) + t416;
t543 = t532 * t414 - t527 * t415;
t537 = -t527 * t414 - t532 * t415;
t493 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t525 + Ifges(4,4) * t530) * t521;
t492 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t525 + Ifges(4,2) * t530) * t521;
t491 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t525 + Ifges(4,6) * t530) * t521;
t477 = Ifges(5,1) * t495 + Ifges(5,4) * t494 + Ifges(5,5) * t520;
t476 = Ifges(5,4) * t495 + Ifges(5,2) * t494 + Ifges(5,6) * t520;
t475 = Ifges(5,5) * t495 + Ifges(5,6) * t494 + Ifges(5,3) * t520;
t458 = Ifges(6,1) * t480 + Ifges(6,4) * t479 + Ifges(6,5) * t513;
t457 = Ifges(6,4) * t480 + Ifges(6,2) * t479 + Ifges(6,6) * t513;
t456 = Ifges(6,5) * t480 + Ifges(6,6) * t479 + Ifges(6,3) * t513;
t436 = mrSges(6,2) * t449 - mrSges(6,3) * t444 + Ifges(6,1) * t455 + Ifges(6,4) * t454 + Ifges(6,5) * t512 + t479 * t456 - t513 * t457;
t435 = -mrSges(6,1) * t449 + mrSges(6,3) * t445 + Ifges(6,4) * t455 + Ifges(6,2) * t454 + Ifges(6,6) * t512 - t480 * t456 + t513 * t458;
t422 = mrSges(5,2) * t467 - mrSges(5,3) * t450 + Ifges(5,1) * t470 + Ifges(5,4) * t469 + Ifges(5,5) * t518 - pkin(9) * t434 - t523 * t435 + t528 * t436 + t494 * t475 - t520 * t476;
t421 = -mrSges(5,1) * t467 + mrSges(5,3) * t451 + Ifges(5,4) * t470 + Ifges(5,2) * t469 + Ifges(5,6) * t518 - pkin(4) * t538 + pkin(9) * t539 + t528 * t435 + t523 * t436 - t495 * t475 + t520 * t477;
t412 = mrSges(4,2) * t483 - mrSges(4,3) * t473 + Ifges(4,1) * t500 + Ifges(4,4) * t501 + Ifges(4,5) * qJDD(3) - pkin(8) * t428 - qJD(3) * t492 - t524 * t421 + t529 * t422 + t491 * t546;
t411 = -mrSges(4,1) * t483 + mrSges(4,3) * t474 + Ifges(4,4) * t500 + Ifges(4,2) * t501 + Ifges(4,6) * qJDD(3) - pkin(3) * t535 + pkin(8) * t540 + qJD(3) * t493 + t529 * t421 + t524 * t422 - t491 * t547;
t410 = mrSges(3,1) * g(1) - Ifges(4,3) * qJDD(3) + (-t525 * t492 + t530 * t493) * t521 - Ifges(5,3) * t518 + Ifges(3,6) * t519 - Ifges(6,3) * t512 + t517 * Ifges(3,5) + t494 * t477 - t495 * t476 - Ifges(4,5) * t500 - Ifges(4,6) * t501 + mrSges(3,3) * t486 - Ifges(5,5) * t470 - mrSges(4,1) * t473 + mrSges(4,2) * t474 + t479 * t458 - t480 * t457 - Ifges(5,6) * t469 - Ifges(6,6) * t454 - Ifges(6,5) * t455 - mrSges(5,1) * t450 + mrSges(5,2) * t451 - mrSges(6,1) * t444 + mrSges(6,2) * t445 - pkin(4) * t434 - pkin(3) * t428 - pkin(2) * t420;
t409 = -mrSges(3,2) * g(1) - mrSges(3,3) * t485 + Ifges(3,5) * t519 - t517 * Ifges(3,6) - pkin(7) * t420 - t525 * t411 + t530 * t412;
t408 = -mrSges(2,2) * g(1) - mrSges(2,3) * t511 + Ifges(2,5) * qJDD(1) - t533 * Ifges(2,6) - pkin(6) * t416 + t531 * t409 - t526 * t410;
t407 = Ifges(2,6) * qJDD(1) + t533 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t510 + t526 * t409 + t531 * t410 - pkin(1) * (-m(3) * g(1) + t420) + pkin(6) * t542;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t420; -m(1) * g(2) + t537; -m(1) * g(3) + t543; mrSges(2,1) * t511 + mrSges(3,1) * t485 - mrSges(1,2) * g(3) - mrSges(2,2) * t510 - mrSges(3,2) * t486 + mrSges(1,3) * g(2) + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t519 + pkin(1) * t416 + pkin(2) * t534 + pkin(7) * t541 + t530 * t411 + t525 * t412; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t543 - t532 * t407 - t527 * t408; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t537 - t527 * t407 + t532 * t408;];
tauB = t1;
