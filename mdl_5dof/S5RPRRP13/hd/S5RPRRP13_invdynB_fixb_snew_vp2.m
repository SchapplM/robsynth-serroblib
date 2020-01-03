% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP13
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP13_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP13_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP13_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP13_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:41
% EndTime: 2019-12-31 18:58:43
% DurationCPUTime: 1.49s
% Computational Cost: add. (10688->240), mult. (20012->280), div. (0->0), fcn. (10784->6), ass. (0->95)
t515 = Ifges(5,1) + Ifges(6,1);
t505 = Ifges(5,4) - Ifges(6,5);
t514 = -Ifges(5,5) - Ifges(6,4);
t513 = Ifges(5,2) + Ifges(6,3);
t502 = Ifges(5,6) - Ifges(6,6);
t512 = -Ifges(5,3) - Ifges(6,2);
t475 = sin(qJ(1));
t477 = cos(qJ(1));
t462 = -t477 * g(1) - t475 * g(2);
t511 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t462;
t510 = -pkin(1) - pkin(6);
t509 = cos(qJ(4));
t508 = mrSges(2,1) - mrSges(3,2);
t507 = -mrSges(5,3) - mrSges(6,2);
t506 = -Ifges(3,4) + Ifges(2,5);
t503 = (Ifges(3,5) - Ifges(2,6));
t461 = t475 * g(1) - t477 * g(2);
t479 = qJD(1) ^ 2;
t484 = -t479 * qJ(2) + qJDD(2) - t461;
t440 = t510 * qJDD(1) + t484;
t474 = sin(qJ(3));
t476 = cos(qJ(3));
t432 = -t476 * g(3) + t474 * t440;
t455 = (mrSges(4,1) * t474 + mrSges(4,2) * t476) * qJD(1);
t494 = qJD(1) * qJD(3);
t490 = t476 * t494;
t457 = -t474 * qJDD(1) - t490;
t496 = qJD(1) * t476;
t460 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t496;
t439 = t510 * t479 - t511;
t491 = t474 * t494;
t458 = t476 * qJDD(1) - t491;
t410 = (-t458 + t491) * pkin(7) + (-t457 + t490) * pkin(3) + t439;
t456 = (pkin(3) * t474 - pkin(7) * t476) * qJD(1);
t478 = qJD(3) ^ 2;
t495 = t474 * qJD(1);
t413 = -t478 * pkin(3) + qJDD(3) * pkin(7) - t456 * t495 + t432;
t473 = sin(qJ(4));
t408 = t473 * t410 + t509 * t413;
t454 = t473 * qJD(3) + t509 * t496;
t424 = t454 * qJD(4) - t509 * qJDD(3) + t473 * t458;
t464 = qJD(4) + t495;
t434 = t464 * mrSges(5,1) - t454 * mrSges(5,3);
t452 = qJDD(4) - t457;
t453 = -t509 * qJD(3) + t473 * t496;
t428 = t453 * pkin(4) - t454 * qJ(5);
t463 = t464 ^ 2;
t404 = -t463 * pkin(4) + t452 * qJ(5) + 0.2e1 * qJD(5) * t464 - t453 * t428 + t408;
t435 = -t464 * mrSges(6,1) + t454 * mrSges(6,2);
t492 = m(6) * t404 + t452 * mrSges(6,3) + t464 * t435;
t429 = t453 * mrSges(6,1) - t454 * mrSges(6,3);
t497 = -t453 * mrSges(5,1) - t454 * mrSges(5,2) - t429;
t399 = m(5) * t408 - t452 * mrSges(5,2) + t507 * t424 - t464 * t434 + t497 * t453 + t492;
t407 = t509 * t410 - t473 * t413;
t425 = -t453 * qJD(4) + t473 * qJDD(3) + t509 * t458;
t433 = -t464 * mrSges(5,2) - t453 * mrSges(5,3);
t405 = -t452 * pkin(4) - t463 * qJ(5) + t454 * t428 + qJDD(5) - t407;
t436 = -t453 * mrSges(6,2) + t464 * mrSges(6,3);
t486 = -m(6) * t405 + t452 * mrSges(6,1) + t464 * t436;
t401 = m(5) * t407 + t452 * mrSges(5,1) + t507 * t425 + t464 * t433 + t497 * t454 + t486;
t487 = t509 * t399 - t473 * t401;
t394 = m(4) * t432 - qJDD(3) * mrSges(4,2) + t457 * mrSges(4,3) - qJD(3) * t460 - t455 * t495 + t487;
t431 = t474 * g(3) + t476 * t440;
t459 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t495;
t412 = -qJDD(3) * pkin(3) - t478 * pkin(7) + t456 * t496 - t431;
t406 = -0.2e1 * qJD(5) * t454 + (t453 * t464 - t425) * qJ(5) + (t454 * t464 + t424) * pkin(4) + t412;
t402 = m(6) * t406 + t424 * mrSges(6,1) - t425 * mrSges(6,3) - t454 * t435 + t453 * t436;
t480 = -m(5) * t412 - t424 * mrSges(5,1) - t425 * mrSges(5,2) - t453 * t433 - t454 * t434 - t402;
t396 = m(4) * t431 + qJDD(3) * mrSges(4,1) - t458 * mrSges(4,3) + qJD(3) * t459 - t455 * t496 + t480;
t387 = t474 * t394 + t476 * t396;
t442 = -qJDD(1) * pkin(1) + t484;
t483 = -m(3) * t442 + (t479 * mrSges(3,3)) - t387;
t385 = m(2) * t461 - (t479 * mrSges(2,2)) + t508 * qJDD(1) + t483;
t441 = t479 * pkin(1) + t511;
t395 = t473 * t399 + t509 * t401;
t482 = -m(4) * t439 + t457 * mrSges(4,1) - t458 * mrSges(4,2) - t459 * t495 - t460 * t496 - t395;
t481 = -m(3) * t441 + (t479 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t482;
t390 = m(2) * t462 - (t479 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t481;
t501 = t477 * t385 + t475 * t390;
t500 = t513 * t453 - t505 * t454 - t502 * t464;
t499 = t502 * t453 + t514 * t454 + t512 * t464;
t498 = -t505 * t453 + t515 * t454 - t514 * t464;
t489 = -t475 * t385 + t477 * t390;
t488 = t476 * t394 - t474 * t396;
t446 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t476 - Ifges(4,4) * t474) * qJD(1);
t445 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t476 - Ifges(4,2) * t474) * qJD(1);
t444 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t476 - Ifges(4,6) * t474) * qJD(1);
t392 = mrSges(5,2) * t412 + mrSges(6,2) * t405 - mrSges(5,3) * t407 - mrSges(6,3) * t406 - qJ(5) * t402 - t505 * t424 + t515 * t425 - t452 * t514 + t499 * t453 + t500 * t464;
t391 = -mrSges(5,1) * t412 - mrSges(6,1) * t406 + mrSges(6,2) * t404 + mrSges(5,3) * t408 - pkin(4) * t402 - t513 * t424 + t505 * t425 + t502 * t452 + t499 * t454 + t498 * t464;
t386 = -m(3) * g(3) + t488;
t383 = Ifges(4,4) * t458 + Ifges(4,2) * t457 + Ifges(4,6) * qJDD(3) - t444 * t496 + qJD(3) * t446 - mrSges(4,1) * t439 + mrSges(4,3) * t432 - mrSges(5,1) * t407 + mrSges(5,2) * t408 + mrSges(6,1) * t405 - mrSges(6,3) * t404 - pkin(4) * t486 - qJ(5) * t492 - pkin(3) * t395 + (pkin(4) * t429 + t500) * t454 + (qJ(5) * t429 - t498) * t453 + t512 * t452 + (pkin(4) * mrSges(6,2) + t514) * t425 + (qJ(5) * mrSges(6,2) + t502) * t424;
t382 = mrSges(4,2) * t439 - mrSges(4,3) * t431 + Ifges(4,1) * t458 + Ifges(4,4) * t457 + Ifges(4,5) * qJDD(3) - pkin(7) * t395 - qJD(3) * t445 - t473 * t391 + t509 * t392 - t444 * t495;
t381 = -qJ(2) * t386 - mrSges(2,3) * t461 + pkin(2) * t387 + mrSges(3,1) * t442 + pkin(7) * t487 + t473 * t392 + t509 * t391 + pkin(3) * t480 + Ifges(4,5) * t458 + Ifges(4,6) * t457 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t431 - mrSges(4,2) * t432 + (t503 * t479) + t506 * qJDD(1) + (t476 * t445 + t474 * t446) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t380 = -mrSges(3,1) * t441 + mrSges(2,3) * t462 - pkin(1) * t386 - pkin(2) * t482 - pkin(6) * t488 + t508 * g(3) - t503 * qJDD(1) - t474 * t382 - t476 * t383 + t506 * t479;
t1 = [-m(1) * g(1) + t489; -m(1) * g(2) + t501; (-m(1) - m(2) - m(3)) * g(3) + t488; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t501 - t475 * t380 + t477 * t381; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t489 + t477 * t380 + t475 * t381; pkin(1) * t483 + qJ(2) * t481 - t474 * t383 - pkin(6) * t387 + mrSges(2,1) * t461 - mrSges(2,2) * t462 + mrSges(3,2) * t442 - mrSges(3,3) * t441 + t476 * t382 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
