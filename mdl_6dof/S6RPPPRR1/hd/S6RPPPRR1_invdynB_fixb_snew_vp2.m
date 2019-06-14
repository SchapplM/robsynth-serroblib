% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPPRR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2019-05-05 13:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPPRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:32:03
% EndTime: 2019-05-05 13:32:04
% DurationCPUTime: 1.60s
% Computational Cost: add. (16458->245), mult. (28920->288), div. (0->0), fcn. (13970->8), ass. (0->95)
t505 = sin(qJ(1));
t508 = cos(qJ(1));
t480 = t505 * g(1) - g(2) * t508;
t472 = qJDD(1) * pkin(1) + t480;
t481 = -g(1) * t508 - g(2) * t505;
t510 = qJD(1) ^ 2;
t474 = -pkin(1) * t510 + t481;
t501 = sin(pkin(9));
t502 = cos(pkin(9));
t453 = t502 * t472 - t501 * t474;
t445 = -qJDD(1) * pkin(2) - t510 * qJ(3) + qJDD(3) - t453;
t442 = -qJDD(1) * qJ(4) - (2 * qJD(4) * qJD(1)) + t445;
t454 = t501 * t472 + t502 * t474;
t536 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t454;
t535 = mrSges(3,1) - mrSges(4,2);
t496 = -g(3) + qJDD(2);
t504 = sin(qJ(5));
t534 = t504 * t496;
t444 = pkin(2) * t510 - t536;
t443 = qJDD(4) + (-pkin(2) - qJ(4)) * t510 + t536;
t440 = -qJDD(1) * pkin(7) + t443;
t507 = cos(qJ(5));
t436 = t504 * t440 + t507 * t496;
t473 = (mrSges(6,1) * t504 + mrSges(6,2) * t507) * qJD(1);
t530 = qJD(1) * qJD(5);
t523 = t507 * t530;
t476 = -qJDD(1) * t504 - t523;
t531 = qJD(1) * t507;
t479 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t531;
t439 = -t510 * pkin(7) - t442;
t524 = t504 * t530;
t477 = qJDD(1) * t507 - t524;
t432 = (-t477 + t524) * pkin(8) + (-t476 + t523) * pkin(5) + t439;
t475 = (pkin(5) * t504 - pkin(8) * t507) * qJD(1);
t509 = qJD(5) ^ 2;
t532 = qJD(1) * t504;
t434 = -pkin(5) * t509 + qJDD(5) * pkin(8) - t475 * t532 + t436;
t503 = sin(qJ(6));
t506 = cos(qJ(6));
t430 = t432 * t506 - t434 * t503;
t470 = qJD(5) * t506 - t503 * t531;
t452 = qJD(6) * t470 + qJDD(5) * t503 + t477 * t506;
t471 = qJD(5) * t503 + t506 * t531;
t455 = -mrSges(7,1) * t470 + mrSges(7,2) * t471;
t482 = qJD(6) + t532;
t456 = -mrSges(7,2) * t482 + mrSges(7,3) * t470;
t469 = qJDD(6) - t476;
t428 = m(7) * t430 + mrSges(7,1) * t469 - mrSges(7,3) * t452 - t455 * t471 + t456 * t482;
t431 = t432 * t503 + t434 * t506;
t451 = -qJD(6) * t471 + qJDD(5) * t506 - t477 * t503;
t457 = mrSges(7,1) * t482 - mrSges(7,3) * t471;
t429 = m(7) * t431 - mrSges(7,2) * t469 + mrSges(7,3) * t451 + t455 * t470 - t457 * t482;
t519 = -t428 * t503 + t506 * t429;
t419 = m(6) * t436 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t476 - qJD(5) * t479 - t473 * t532 + t519;
t435 = t440 * t507 - t534;
t478 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t532;
t433 = -qJDD(5) * pkin(5) - t509 * pkin(8) + t534 + (qJD(1) * t475 - t440) * t507;
t512 = -m(7) * t433 + t451 * mrSges(7,1) - mrSges(7,2) * t452 + t470 * t456 - t457 * t471;
t424 = m(6) * t435 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t477 + qJD(5) * t478 - t473 * t531 + t512;
t413 = t504 * t419 + t507 * t424;
t518 = -m(5) * t443 - qJDD(1) * mrSges(5,2) - t413;
t513 = -m(4) * t444 + t510 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t518;
t411 = m(3) * t454 - qJDD(1) * mrSges(3,2) + (-mrSges(3,1) - mrSges(5,3)) * t510 + t513;
t420 = t506 * t428 + t503 * t429;
t514 = -m(6) * t439 + t476 * mrSges(6,1) - t477 * mrSges(6,2) - t478 * t532 - t479 * t531 - t420;
t416 = m(5) * t442 - t510 * mrSges(5,2) - qJDD(1) * mrSges(5,3) + t514;
t511 = -m(4) * t445 + t510 * mrSges(4,3) - t416;
t415 = m(3) * t453 - t510 * mrSges(3,2) + qJDD(1) * t535 + t511;
t406 = t501 * t411 + t502 * t415;
t404 = m(2) * t480 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t510 + t406;
t521 = t502 * t411 - t415 * t501;
t405 = m(2) * t481 - mrSges(2,1) * t510 - qJDD(1) * mrSges(2,2) + t521;
t533 = t508 * t404 + t505 * t405;
t527 = Ifges(3,5) - Ifges(4,4) + Ifges(5,5);
t526 = Ifges(3,6) - Ifges(4,5) - Ifges(5,4);
t522 = -t404 * t505 + t508 * t405;
t520 = t507 * t419 - t504 * t424;
t517 = m(5) * t496 + t520;
t412 = m(4) * t496 + t517;
t515 = m(3) * t496 + t412;
t464 = (Ifges(6,5) * qJD(5)) + (Ifges(6,1) * t507 - Ifges(6,4) * t504) * qJD(1);
t463 = (Ifges(6,6) * qJD(5)) + (Ifges(6,4) * t507 - Ifges(6,2) * t504) * qJD(1);
t462 = (Ifges(6,3) * qJD(5)) + (Ifges(6,5) * t507 - Ifges(6,6) * t504) * qJD(1);
t448 = Ifges(7,1) * t471 + Ifges(7,4) * t470 + Ifges(7,5) * t482;
t447 = Ifges(7,4) * t471 + Ifges(7,2) * t470 + Ifges(7,6) * t482;
t446 = Ifges(7,5) * t471 + Ifges(7,6) * t470 + Ifges(7,3) * t482;
t422 = mrSges(7,2) * t433 - mrSges(7,3) * t430 + Ifges(7,1) * t452 + Ifges(7,4) * t451 + Ifges(7,5) * t469 + t446 * t470 - t447 * t482;
t421 = -mrSges(7,1) * t433 + mrSges(7,3) * t431 + Ifges(7,4) * t452 + Ifges(7,2) * t451 + Ifges(7,6) * t469 - t446 * t471 + t448 * t482;
t408 = -mrSges(6,1) * t439 - mrSges(7,1) * t430 + mrSges(7,2) * t431 + mrSges(6,3) * t436 + Ifges(6,4) * t477 - Ifges(7,5) * t452 + Ifges(6,2) * t476 + Ifges(6,6) * qJDD(5) - Ifges(7,6) * t451 - Ifges(7,3) * t469 - pkin(5) * t420 + qJD(5) * t464 - t447 * t471 + t448 * t470 - t462 * t531;
t407 = mrSges(6,2) * t439 - mrSges(6,3) * t435 + Ifges(6,1) * t477 + Ifges(6,4) * t476 + Ifges(6,5) * qJDD(5) - pkin(8) * t420 - qJD(5) * t463 - t421 * t503 + t422 * t506 - t462 * t532;
t400 = Ifges(6,3) * qJDD(5) - pkin(3) * t518 + pkin(8) * t519 + t503 * t422 - qJ(4) * t517 + t506 * t421 + pkin(5) * t512 + Ifges(6,6) * t476 + Ifges(6,5) * t477 + mrSges(3,3) * t454 + mrSges(6,1) * t435 - mrSges(6,2) * t436 + mrSges(5,1) * t443 - mrSges(4,1) * t444 - pkin(2) * t412 + pkin(4) * t413 + (t463 * t507 + t464 * t504) * qJD(1) + (-mrSges(5,3) * pkin(3) + t527) * t510 + (-mrSges(5,3) - t535) * t496 + t526 * qJDD(1);
t399 = -qJ(3) * t412 - mrSges(3,3) * t453 + pkin(3) * t416 + mrSges(4,1) * t445 + t507 * t408 + pkin(4) * t514 + pkin(7) * t520 + t504 * t407 + mrSges(5,1) * t442 - t526 * t510 + (mrSges(3,2) - mrSges(5,2) - mrSges(4,3)) * t496 + t527 * qJDD(1);
t398 = -mrSges(2,2) * g(3) - mrSges(2,3) * t480 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t510 - qJ(2) * t406 + t399 * t502 - t400 * t501;
t397 = mrSges(2,1) * g(3) + mrSges(2,3) * t481 + t510 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t515 + qJ(2) * t521 + t501 * t399 + t502 * t400;
t1 = [-m(1) * g(1) + t522; -m(1) * g(2) + t533; (-m(1) - m(2)) * g(3) + t515; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t533 - t505 * t397 + t508 * t398; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t522 + t508 * t397 + t505 * t398; pkin(1) * t406 - mrSges(2,2) * t481 + mrSges(2,1) * t480 + pkin(2) * t511 + mrSges(3,1) * t453 - mrSges(3,2) * t454 + qJ(3) * (-mrSges(5,3) * t510 + t513) - qJ(4) * t416 + mrSges(4,2) * t445 - mrSges(4,3) * t444 + t507 * t407 - t504 * t408 - pkin(7) * t413 + mrSges(5,2) * t443 - mrSges(5,3) * t442 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(4,2) * pkin(2) + Ifges(4,1) + Ifges(5,1) + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
