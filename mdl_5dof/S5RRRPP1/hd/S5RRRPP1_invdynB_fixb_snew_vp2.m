% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:36
% EndTime: 2019-12-31 20:49:39
% DurationCPUTime: 2.47s
% Computational Cost: add. (32962->249), mult. (43897->304), div. (0->0), fcn. (25110->8), ass. (0->100)
t532 = Ifges(5,1) + Ifges(6,1);
t526 = Ifges(5,4) - Ifges(6,5);
t525 = Ifges(5,5) + Ifges(6,4);
t531 = Ifges(5,2) + Ifges(6,3);
t530 = -Ifges(6,2) - Ifges(5,3);
t524 = Ifges(5,6) - Ifges(6,6);
t529 = -2 * qJD(4);
t491 = qJD(1) + qJD(2);
t489 = t491 ^ 2;
t528 = pkin(3) * t489;
t527 = -mrSges(5,3) - mrSges(6,2);
t523 = cos(pkin(8));
t496 = sin(qJ(3));
t522 = t491 * t496;
t499 = cos(qJ(3));
t521 = t491 * t499;
t498 = sin(qJ(1));
t501 = cos(qJ(1));
t486 = t498 * g(1) - t501 * g(2);
t480 = qJDD(1) * pkin(1) + t486;
t487 = -t501 * g(1) - t498 * g(2);
t503 = qJD(1) ^ 2;
t481 = -t503 * pkin(1) + t487;
t497 = sin(qJ(2));
t500 = cos(qJ(2));
t457 = t497 * t480 + t500 * t481;
t490 = qJDD(1) + qJDD(2);
t450 = -t489 * pkin(2) + t490 * pkin(7) + t457;
t520 = t496 * t450;
t514 = qJD(3) * t491;
t475 = t496 * t490 + t499 * t514;
t431 = qJDD(3) * pkin(3) - t475 * qJ(4) - t520 + (qJ(4) * t514 + t496 * t528 - g(3)) * t499;
t435 = -t496 * g(3) + t499 * t450;
t476 = t499 * t490 - t496 * t514;
t482 = qJD(3) * pkin(3) - qJ(4) * t522;
t494 = t499 ^ 2;
t432 = t476 * qJ(4) - qJD(3) * t482 - t494 * t528 + t435;
t495 = sin(pkin(8));
t465 = t495 * t522 - t523 * t521;
t428 = t495 * t431 + t523 * t432 + t465 * t529;
t454 = t495 * t475 - t523 * t476;
t466 = (t495 * t499 + t523 * t496) * t491;
t461 = qJD(3) * mrSges(5,1) - t466 * mrSges(5,3);
t446 = t465 * pkin(4) - t466 * qJ(5);
t502 = qJD(3) ^ 2;
t423 = -t502 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t465 * t446 + t428;
t462 = -qJD(3) * mrSges(6,1) + t466 * mrSges(6,2);
t513 = m(6) * t423 + qJDD(3) * mrSges(6,3) + qJD(3) * t462;
t447 = t465 * mrSges(6,1) - t466 * mrSges(6,3);
t515 = -t465 * mrSges(5,1) - t466 * mrSges(5,2) - t447;
t419 = m(5) * t428 - qJDD(3) * mrSges(5,2) - qJD(3) * t461 + t527 * t454 + t515 * t465 + t513;
t506 = t523 * t431 - t495 * t432;
t427 = t466 * t529 + t506;
t455 = t523 * t475 + t495 * t476;
t460 = -qJD(3) * mrSges(5,2) - t465 * mrSges(5,3);
t424 = -qJDD(3) * pkin(4) - t502 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t446) * t466 - t506;
t463 = -t465 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t508 = -m(6) * t424 + qJDD(3) * mrSges(6,1) + qJD(3) * t463;
t420 = m(5) * t427 + qJDD(3) * mrSges(5,1) + qJD(3) * t460 + t527 * t455 + t515 * t466 + t508;
t413 = t495 * t419 + t523 * t420;
t434 = -t499 * g(3) - t520;
t474 = (-mrSges(4,1) * t499 + mrSges(4,2) * t496) * t491;
t484 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t521;
t409 = m(4) * t434 + qJDD(3) * mrSges(4,1) - t475 * mrSges(4,3) + qJD(3) * t484 - t474 * t522 + t413;
t483 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t522;
t509 = t523 * t419 - t495 * t420;
t410 = m(4) * t435 - qJDD(3) * mrSges(4,2) + t476 * mrSges(4,3) - qJD(3) * t483 + t474 * t521 + t509;
t510 = -t496 * t409 + t499 * t410;
t404 = m(3) * t457 - t489 * mrSges(3,1) - t490 * mrSges(3,2) + t510;
t456 = t500 * t480 - t497 * t481;
t507 = -t490 * pkin(2) - t456;
t449 = -t489 * pkin(7) + t507;
t433 = -t476 * pkin(3) + qJDD(4) + t482 * t522 + (-qJ(4) * t494 - pkin(7)) * t489 + t507;
t426 = -0.2e1 * qJD(5) * t466 + (qJD(3) * t465 - t455) * qJ(5) + (qJD(3) * t466 + t454) * pkin(4) + t433;
t421 = m(6) * t426 + t454 * mrSges(6,1) - t455 * mrSges(6,3) - t466 * t462 + t465 * t463;
t505 = m(5) * t433 + t454 * mrSges(5,1) + t455 * mrSges(5,2) + t465 * t460 + t466 * t461 + t421;
t504 = -m(4) * t449 + t476 * mrSges(4,1) - t475 * mrSges(4,2) - t483 * t522 + t484 * t521 - t505;
t415 = m(3) * t456 + t490 * mrSges(3,1) - t489 * mrSges(3,2) + t504;
t401 = t497 * t404 + t500 * t415;
t399 = m(2) * t486 + qJDD(1) * mrSges(2,1) - t503 * mrSges(2,2) + t401;
t511 = t500 * t404 - t497 * t415;
t400 = m(2) * t487 - t503 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t511;
t519 = t501 * t399 + t498 * t400;
t405 = t499 * t409 + t496 * t410;
t518 = -t524 * qJD(3) + t531 * t465 - t526 * t466;
t517 = t530 * qJD(3) + t524 * t465 - t525 * t466;
t516 = t525 * qJD(3) - t526 * t465 + t532 * t466;
t512 = -t498 * t399 + t501 * t400;
t470 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t496 + Ifges(4,4) * t499) * t491;
t469 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t496 + Ifges(4,2) * t499) * t491;
t468 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t496 + Ifges(4,6) * t499) * t491;
t412 = mrSges(5,2) * t433 + mrSges(6,2) * t424 - mrSges(5,3) * t427 - mrSges(6,3) * t426 - qJ(5) * t421 + t518 * qJD(3) + t525 * qJDD(3) - t526 * t454 + t532 * t455 + t517 * t465;
t411 = -mrSges(5,1) * t433 - mrSges(6,1) * t426 + mrSges(6,2) * t423 + mrSges(5,3) * t428 - pkin(4) * t421 + t516 * qJD(3) + t524 * qJDD(3) - t531 * t454 + t526 * t455 + t517 * t466;
t395 = mrSges(4,2) * t449 - mrSges(4,3) * t434 + Ifges(4,1) * t475 + Ifges(4,4) * t476 + Ifges(4,5) * qJDD(3) - qJ(4) * t413 - qJD(3) * t469 - t495 * t411 + t523 * t412 + t468 * t521;
t394 = -mrSges(4,1) * t449 + mrSges(4,3) * t435 + Ifges(4,4) * t475 + Ifges(4,2) * t476 + Ifges(4,6) * qJDD(3) - pkin(3) * t505 + qJ(4) * t509 + qJD(3) * t470 + t523 * t411 + t495 * t412 - t468 * t522;
t393 = mrSges(3,1) * g(3) - pkin(4) * t508 + t489 * Ifges(3,5) + Ifges(3,6) * t490 - qJ(5) * t513 - Ifges(4,5) * t475 - Ifges(4,6) * t476 + mrSges(3,3) * t457 + mrSges(4,2) * t435 - mrSges(5,1) * t427 + mrSges(5,2) * t428 - mrSges(4,1) * t434 - mrSges(6,3) * t423 + mrSges(6,1) * t424 - pkin(3) * t413 - pkin(2) * t405 + (-t496 * t469 + t499 * t470) * t491 + (pkin(4) * t447 + t518) * t466 + (qJ(5) * t447 - t516) * t465 + (pkin(4) * mrSges(6,2) - t525) * t455 + (qJ(5) * mrSges(6,2) + t524) * t454 + (-Ifges(4,3) + t530) * qJDD(3);
t392 = -mrSges(3,2) * g(3) - mrSges(3,3) * t456 + Ifges(3,5) * t490 - t489 * Ifges(3,6) - pkin(7) * t405 - t496 * t394 + t499 * t395;
t391 = -mrSges(2,2) * g(3) - mrSges(2,3) * t486 + Ifges(2,5) * qJDD(1) - t503 * Ifges(2,6) - pkin(6) * t401 + t500 * t392 - t497 * t393;
t390 = Ifges(2,6) * qJDD(1) + t503 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t487 + t497 * t392 + t500 * t393 - pkin(1) * (-m(3) * g(3) + t405) + pkin(6) * t511;
t1 = [-m(1) * g(1) + t512; -m(1) * g(2) + t519; (-m(1) - m(2) - m(3)) * g(3) + t405; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t519 - t498 * t390 + t501 * t391; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t512 + t501 * t390 + t498 * t391; -mrSges(1,1) * g(2) + mrSges(2,1) * t486 + mrSges(3,1) * t456 + mrSges(1,2) * g(1) - mrSges(2,2) * t487 - mrSges(3,2) * t457 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t490 + pkin(1) * t401 + pkin(2) * t504 + pkin(7) * t510 + t499 * t394 + t496 * t395;];
tauB = t1;
