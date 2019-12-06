% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRPP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRPP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:57
% EndTime: 2019-12-05 16:09:01
% DurationCPUTime: 2.09s
% Computational Cost: add. (17720->236), mult. (37864->292), div. (0->0), fcn. (22938->8), ass. (0->95)
t505 = Ifges(5,1) + Ifges(6,1);
t500 = Ifges(5,4) - Ifges(6,5);
t499 = Ifges(5,5) + Ifges(6,4);
t504 = Ifges(5,2) + Ifges(6,3);
t503 = -Ifges(6,2) - Ifges(5,3);
t498 = Ifges(5,6) - Ifges(6,6);
t502 = -2 * qJD(4);
t501 = -mrSges(5,3) - mrSges(6,2);
t497 = cos(pkin(7));
t496 = cos(pkin(8));
t470 = sin(pkin(7));
t458 = -t497 * g(1) - t470 * g(2);
t468 = -g(3) + qJDD(1);
t472 = sin(qJ(2));
t474 = cos(qJ(2));
t439 = t474 * t458 + t472 * t468;
t476 = qJD(2) ^ 2;
t431 = -t476 * pkin(2) + qJDD(2) * pkin(6) + t439;
t457 = t470 * g(1) - t497 * g(2);
t471 = sin(qJ(3));
t473 = cos(qJ(3));
t410 = -t471 * t431 - t473 * t457;
t488 = qJD(2) * qJD(3);
t486 = t473 * t488;
t455 = t471 * qJDD(2) + t486;
t407 = (-t455 + t486) * qJ(4) + (t471 * t473 * t476 + qJDD(3)) * pkin(3) + t410;
t411 = t473 * t431 - t471 * t457;
t456 = t473 * qJDD(2) - t471 * t488;
t490 = qJD(2) * t471;
t459 = qJD(3) * pkin(3) - qJ(4) * t490;
t467 = t473 ^ 2;
t408 = -t467 * t476 * pkin(3) + t456 * qJ(4) - qJD(3) * t459 + t411;
t469 = sin(pkin(8));
t489 = qJD(2) * t473;
t441 = t469 * t490 - t496 * t489;
t404 = t469 * t407 + t496 * t408 + t441 * t502;
t427 = t469 * t455 - t496 * t456;
t442 = (t469 * t473 + t496 * t471) * qJD(2);
t435 = qJD(3) * mrSges(5,1) - t442 * mrSges(5,3);
t421 = t441 * pkin(4) - t442 * qJ(5);
t475 = qJD(3) ^ 2;
t399 = -t475 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t441 * t421 + t404;
t436 = -qJD(3) * mrSges(6,1) + t442 * mrSges(6,2);
t487 = m(6) * t399 + qJDD(3) * mrSges(6,3) + qJD(3) * t436;
t422 = t441 * mrSges(6,1) - t442 * mrSges(6,3);
t491 = -t441 * mrSges(5,1) - t442 * mrSges(5,2) - t422;
t395 = m(5) * t404 - qJDD(3) * mrSges(5,2) - qJD(3) * t435 + t501 * t427 + t491 * t441 + t487;
t479 = t496 * t407 - t469 * t408;
t403 = t442 * t502 + t479;
t428 = t496 * t455 + t469 * t456;
t434 = -qJD(3) * mrSges(5,2) - t441 * mrSges(5,3);
t400 = -qJDD(3) * pkin(4) - t475 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t421) * t442 - t479;
t437 = -t441 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t481 = -m(6) * t400 + qJDD(3) * mrSges(6,1) + qJD(3) * t437;
t396 = m(5) * t403 + qJDD(3) * mrSges(5,1) + qJD(3) * t434 + t501 * t428 + t491 * t442 + t481;
t389 = t469 * t395 + t496 * t396;
t454 = (-mrSges(4,1) * t473 + mrSges(4,2) * t471) * qJD(2);
t461 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t489;
t385 = m(4) * t410 + qJDD(3) * mrSges(4,1) - t455 * mrSges(4,3) + qJD(3) * t461 - t454 * t490 + t389;
t460 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t490;
t482 = t496 * t395 - t469 * t396;
t386 = m(4) * t411 - qJDD(3) * mrSges(4,2) + t456 * mrSges(4,3) - qJD(3) * t460 + t454 * t489 + t482;
t483 = -t471 * t385 + t473 * t386;
t380 = m(3) * t439 - t476 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t483;
t438 = -t472 * t458 + t474 * t468;
t480 = -qJDD(2) * pkin(2) - t438;
t430 = -t476 * pkin(6) + t480;
t409 = -t456 * pkin(3) + qJDD(4) + t459 * t490 + (-qJ(4) * t467 - pkin(6)) * t476 + t480;
t402 = -0.2e1 * qJD(5) * t442 + (qJD(3) * t441 - t428) * qJ(5) + (qJD(3) * t442 + t427) * pkin(4) + t409;
t397 = m(6) * t402 + t427 * mrSges(6,1) - t428 * mrSges(6,3) - t442 * t436 + t441 * t437;
t478 = m(5) * t409 + t427 * mrSges(5,1) + t428 * mrSges(5,2) + t441 * t434 + t442 * t435 + t397;
t477 = -m(4) * t430 + t456 * mrSges(4,1) - t455 * mrSges(4,2) - t460 * t490 + t461 * t489 - t478;
t394 = m(3) * t438 + qJDD(2) * mrSges(3,1) - t476 * mrSges(3,2) + t477;
t484 = t474 * t380 - t472 * t394;
t376 = m(2) * t458 + t484;
t383 = t473 * t385 + t471 * t386;
t382 = (m(2) + m(3)) * t457 - t383;
t495 = t470 * t376 + t497 * t382;
t377 = t472 * t380 + t474 * t394;
t494 = -t498 * qJD(3) + t504 * t441 - t500 * t442;
t493 = t503 * qJD(3) + t498 * t441 - t499 * t442;
t492 = t499 * qJD(3) - t500 * t441 + t505 * t442;
t485 = t497 * t376 - t470 * t382;
t445 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t471 + Ifges(4,4) * t473) * qJD(2);
t444 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t471 + Ifges(4,2) * t473) * qJD(2);
t443 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t471 + Ifges(4,6) * t473) * qJD(2);
t388 = mrSges(5,2) * t409 + mrSges(6,2) * t400 - mrSges(5,3) * t403 - mrSges(6,3) * t402 - qJ(5) * t397 + t494 * qJD(3) + t499 * qJDD(3) - t500 * t427 + t505 * t428 + t493 * t441;
t387 = -mrSges(5,1) * t409 - mrSges(6,1) * t402 + mrSges(6,2) * t399 + mrSges(5,3) * t404 - pkin(4) * t397 + t492 * qJD(3) + t498 * qJDD(3) - t504 * t427 + t500 * t428 + t493 * t442;
t373 = mrSges(4,2) * t430 - mrSges(4,3) * t410 + Ifges(4,1) * t455 + Ifges(4,4) * t456 + Ifges(4,5) * qJDD(3) - qJ(4) * t389 - qJD(3) * t444 - t469 * t387 + t496 * t388 + t443 * t489;
t372 = -mrSges(4,1) * t430 + mrSges(4,3) * t411 + Ifges(4,4) * t455 + Ifges(4,2) * t456 + Ifges(4,6) * qJDD(3) - pkin(3) * t478 + qJ(4) * t482 + qJD(3) * t445 + t496 * t387 + t469 * t388 - t443 * t490;
t371 = Ifges(3,6) * qJDD(2) + t476 * Ifges(3,5) - qJ(5) * t487 - pkin(4) * t481 - Ifges(4,5) * t455 - Ifges(4,6) * t456 + mrSges(3,1) * t457 + mrSges(3,3) * t439 - mrSges(4,1) * t410 + mrSges(4,2) * t411 + mrSges(6,1) * t400 - mrSges(5,1) * t403 + mrSges(5,2) * t404 - mrSges(6,3) * t399 - pkin(2) * t383 - pkin(3) * t389 + (pkin(4) * t422 + t494) * t442 + (qJ(5) * t422 - t492) * t441 + (pkin(4) * mrSges(6,2) - t499) * t428 + (qJ(5) * mrSges(6,2) + t498) * t427 + (-t471 * t444 + t473 * t445) * qJD(2) + (-Ifges(4,3) + t503) * qJDD(3);
t370 = -mrSges(3,2) * t457 - mrSges(3,3) * t438 + Ifges(3,5) * qJDD(2) - t476 * Ifges(3,6) - pkin(6) * t383 - t471 * t372 + t473 * t373;
t369 = -mrSges(2,1) * t468 - mrSges(3,1) * t438 + mrSges(3,2) * t439 + mrSges(2,3) * t458 - Ifges(3,3) * qJDD(2) - pkin(1) * t377 - pkin(2) * t477 - pkin(6) * t483 - t473 * t372 - t471 * t373;
t368 = mrSges(2,2) * t468 - mrSges(2,3) * t457 - pkin(5) * t377 + t474 * t370 - t472 * t371;
t1 = [-m(1) * g(1) + t485; -m(1) * g(2) + t495; -m(1) * g(3) + m(2) * t468 + t377; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t495 + t497 * t368 - t470 * t369; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t485 + t470 * t368 + t497 * t369; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t457 - mrSges(2,2) * t458 + t472 * t370 + t474 * t371 + pkin(1) * (m(3) * t457 - t383) + pkin(5) * t484;];
tauB = t1;
