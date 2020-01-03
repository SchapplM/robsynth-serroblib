% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:09:08
% EndTime: 2020-01-03 12:09:13
% DurationCPUTime: 5.00s
% Computational Cost: add. (86594->272), mult. (115951->347), div. (0->0), fcn. (73622->10), ass. (0->108)
t508 = qJD(1) + qJD(2);
t504 = t508 ^ 2;
t535 = pkin(3) * t504;
t513 = sin(qJ(3));
t534 = t508 * t513;
t517 = cos(qJ(3));
t533 = t508 * t517;
t515 = sin(qJ(1));
t519 = cos(qJ(1));
t502 = -t519 * g(2) - t515 * g(3);
t496 = qJDD(1) * pkin(1) + t502;
t501 = -t515 * g(2) + t519 * g(3);
t520 = qJD(1) ^ 2;
t497 = -t520 * pkin(1) + t501;
t514 = sin(qJ(2));
t518 = cos(qJ(2));
t475 = t514 * t496 + t518 * t497;
t506 = qJDD(1) + qJDD(2);
t470 = -t504 * pkin(2) + t506 * pkin(7) + t475;
t532 = t513 * t470;
t530 = qJD(3) * t508;
t491 = t513 * t506 + t517 * t530;
t454 = qJDD(3) * pkin(3) - t491 * qJ(4) - t532 + (qJ(4) * t530 + t513 * t535 - g(1)) * t517;
t460 = -t513 * g(1) + t517 * t470;
t492 = t517 * t506 - t513 * t530;
t498 = qJD(3) * pkin(3) - qJ(4) * t534;
t509 = t517 ^ 2;
t455 = t492 * qJ(4) - qJD(3) * t498 - t509 * t535 + t460;
t510 = sin(pkin(9));
t511 = cos(pkin(9));
t483 = (t510 * t517 + t511 * t513) * t508;
t437 = -0.2e1 * qJD(4) * t483 + t511 * t454 - t510 * t455;
t473 = t511 * t491 + t510 * t492;
t482 = (-t510 * t513 + t511 * t517) * t508;
t435 = (qJD(3) * t482 - t473) * pkin(8) + (t482 * t483 + qJDD(3)) * pkin(4) + t437;
t438 = 0.2e1 * qJD(4) * t482 + t510 * t454 + t511 * t455;
t472 = -t510 * t491 + t511 * t492;
t478 = qJD(3) * pkin(4) - t483 * pkin(8);
t481 = t482 ^ 2;
t436 = -t481 * pkin(4) + t472 * pkin(8) - qJD(3) * t478 + t438;
t512 = sin(qJ(5));
t516 = cos(qJ(5));
t433 = t516 * t435 - t512 * t436;
t464 = t516 * t482 - t512 * t483;
t444 = t464 * qJD(5) + t512 * t472 + t516 * t473;
t465 = t512 * t482 + t516 * t483;
t450 = -t464 * mrSges(6,1) + t465 * mrSges(6,2);
t507 = qJD(3) + qJD(5);
t457 = -t507 * mrSges(6,2) + t464 * mrSges(6,3);
t505 = qJDD(3) + qJDD(5);
t431 = m(6) * t433 + t505 * mrSges(6,1) - t444 * mrSges(6,3) - t465 * t450 + t507 * t457;
t434 = t512 * t435 + t516 * t436;
t443 = -t465 * qJD(5) + t516 * t472 - t512 * t473;
t458 = t507 * mrSges(6,1) - t465 * mrSges(6,3);
t432 = m(6) * t434 - t505 * mrSges(6,2) + t443 * mrSges(6,3) + t464 * t450 - t507 * t458;
t423 = t516 * t431 + t512 * t432;
t468 = -t482 * mrSges(5,1) + t483 * mrSges(5,2);
t476 = -qJD(3) * mrSges(5,2) + t482 * mrSges(5,3);
t421 = m(5) * t437 + qJDD(3) * mrSges(5,1) - t473 * mrSges(5,3) + qJD(3) * t476 - t483 * t468 + t423;
t477 = qJD(3) * mrSges(5,1) - t483 * mrSges(5,3);
t525 = -t512 * t431 + t516 * t432;
t422 = m(5) * t438 - qJDD(3) * mrSges(5,2) + t472 * mrSges(5,3) - qJD(3) * t477 + t482 * t468 + t525;
t417 = t511 * t421 + t510 * t422;
t459 = -t517 * g(1) - t532;
t490 = (-mrSges(4,1) * t517 + mrSges(4,2) * t513) * t508;
t500 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t533;
t415 = m(4) * t459 + qJDD(3) * mrSges(4,1) - t491 * mrSges(4,3) + qJD(3) * t500 - t490 * t534 + t417;
t499 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t534;
t526 = -t510 * t421 + t511 * t422;
t416 = m(4) * t460 - qJDD(3) * mrSges(4,2) + t492 * mrSges(4,3) - qJD(3) * t499 + t490 * t533 + t526;
t527 = -t513 * t415 + t517 * t416;
t408 = m(3) * t475 - t504 * mrSges(3,1) - t506 * mrSges(3,2) + t527;
t474 = t518 * t496 - t514 * t497;
t523 = -t506 * pkin(2) - t474;
t469 = -t504 * pkin(7) + t523;
t456 = -t492 * pkin(3) + qJDD(4) + t498 * t534 + (-qJ(4) * t509 - pkin(7)) * t504 + t523;
t440 = -t472 * pkin(4) - t481 * pkin(8) + t483 * t478 + t456;
t524 = m(6) * t440 - t443 * mrSges(6,1) + t444 * mrSges(6,2) - t464 * t457 + t465 * t458;
t522 = m(5) * t456 - t472 * mrSges(5,1) + t473 * mrSges(5,2) - t482 * t476 + t483 * t477 + t524;
t521 = -m(4) * t469 + t492 * mrSges(4,1) - t491 * mrSges(4,2) - t499 * t534 + t500 * t533 - t522;
t427 = m(3) * t474 + t506 * mrSges(3,1) - t504 * mrSges(3,2) + t521;
t528 = t518 * t408 - t514 * t427;
t403 = m(2) * t501 - t520 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t528;
t405 = t514 * t408 + t518 * t427;
t404 = m(2) * t502 + qJDD(1) * mrSges(2,1) - t520 * mrSges(2,2) + t405;
t531 = t515 * t403 + t519 * t404;
t409 = t517 * t415 + t513 * t416;
t529 = -t519 * t403 + t515 * t404;
t486 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t513 + Ifges(4,4) * t517) * t508;
t485 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t513 + Ifges(4,2) * t517) * t508;
t484 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t513 + Ifges(4,6) * t517) * t508;
t463 = Ifges(5,1) * t483 + Ifges(5,4) * t482 + Ifges(5,5) * qJD(3);
t462 = Ifges(5,4) * t483 + Ifges(5,2) * t482 + Ifges(5,6) * qJD(3);
t461 = Ifges(5,5) * t483 + Ifges(5,6) * t482 + Ifges(5,3) * qJD(3);
t447 = Ifges(6,1) * t465 + Ifges(6,4) * t464 + Ifges(6,5) * t507;
t446 = Ifges(6,4) * t465 + Ifges(6,2) * t464 + Ifges(6,6) * t507;
t445 = Ifges(6,5) * t465 + Ifges(6,6) * t464 + Ifges(6,3) * t507;
t425 = mrSges(6,2) * t440 - mrSges(6,3) * t433 + Ifges(6,1) * t444 + Ifges(6,4) * t443 + Ifges(6,5) * t505 + t464 * t445 - t507 * t446;
t424 = -mrSges(6,1) * t440 + mrSges(6,3) * t434 + Ifges(6,4) * t444 + Ifges(6,2) * t443 + Ifges(6,6) * t505 - t465 * t445 + t507 * t447;
t411 = mrSges(5,2) * t456 - mrSges(5,3) * t437 + Ifges(5,1) * t473 + Ifges(5,4) * t472 + Ifges(5,5) * qJDD(3) - pkin(8) * t423 - qJD(3) * t462 - t512 * t424 + t516 * t425 + t482 * t461;
t410 = -mrSges(5,1) * t456 + mrSges(5,3) * t438 + Ifges(5,4) * t473 + Ifges(5,2) * t472 + Ifges(5,6) * qJDD(3) - pkin(4) * t524 + pkin(8) * t525 + qJD(3) * t463 + t516 * t424 + t512 * t425 - t483 * t461;
t399 = mrSges(4,2) * t469 - mrSges(4,3) * t459 + Ifges(4,1) * t491 + Ifges(4,4) * t492 + Ifges(4,5) * qJDD(3) - qJ(4) * t417 - qJD(3) * t485 - t510 * t410 + t511 * t411 + t484 * t533;
t398 = -mrSges(4,1) * t469 + mrSges(4,3) * t460 + Ifges(4,4) * t491 + Ifges(4,2) * t492 + Ifges(4,6) * qJDD(3) - pkin(3) * t522 + qJ(4) * t526 + qJD(3) * t486 + t511 * t410 + t510 * t411 - t484 * t534;
t397 = mrSges(3,1) * g(1) + (-t513 * t485 + t517 * t486) * t508 + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + t504 * Ifges(3,5) - Ifges(6,3) * t505 + Ifges(3,6) * t506 - Ifges(4,5) * t491 - Ifges(4,6) * t492 + t482 * t463 - t483 * t462 - Ifges(5,6) * t472 - Ifges(5,5) * t473 + mrSges(3,3) * t475 + t464 * t447 - t465 * t446 - mrSges(4,1) * t459 + mrSges(4,2) * t460 - Ifges(6,5) * t444 - Ifges(6,6) * t443 - mrSges(5,1) * t437 + mrSges(5,2) * t438 + mrSges(6,2) * t434 - mrSges(6,1) * t433 - pkin(4) * t423 - pkin(3) * t417 - pkin(2) * t409;
t396 = -mrSges(3,2) * g(1) - mrSges(3,3) * t474 + Ifges(3,5) * t506 - t504 * Ifges(3,6) - pkin(7) * t409 - t513 * t398 + t517 * t399;
t395 = -mrSges(2,2) * g(1) - mrSges(2,3) * t502 + Ifges(2,5) * qJDD(1) - t520 * Ifges(2,6) - pkin(6) * t405 + t518 * t396 - t514 * t397;
t394 = Ifges(2,6) * qJDD(1) + t520 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t501 + t514 * t396 + t518 * t397 - pkin(1) * (-m(3) * g(1) + t409) + pkin(6) * t528;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t409; -m(1) * g(2) + t531; -m(1) * g(3) + t529; mrSges(2,1) * t502 + mrSges(3,1) * t474 - mrSges(1,2) * g(3) - mrSges(2,2) * t501 - mrSges(3,2) * t475 + mrSges(1,3) * g(2) + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t506 + pkin(1) * t405 + pkin(2) * t521 + pkin(7) * t527 + t517 * t398 + t513 * t399; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t529 + t519 * t394 + t515 * t395; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t531 + t515 * t394 - t519 * t395;];
tauB = t1;
