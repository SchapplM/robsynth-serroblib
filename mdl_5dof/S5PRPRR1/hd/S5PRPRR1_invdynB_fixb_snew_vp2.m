% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRR1
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:54
% EndTime: 2019-12-05 15:42:58
% DurationCPUTime: 4.02s
% Computational Cost: add. (44310->237), mult. (100911->302), div. (0->0), fcn. (71884->10), ass. (0->106)
t494 = qJD(2) ^ 2;
t486 = cos(pkin(9));
t519 = pkin(3) * t486;
t484 = sin(pkin(9));
t518 = mrSges(4,2) * t484;
t481 = t486 ^ 2;
t517 = t481 * t494;
t485 = sin(pkin(8));
t487 = cos(pkin(8));
t469 = t485 * g(1) - t487 * g(2);
t470 = -t487 * g(1) - t485 * g(2);
t490 = sin(qJ(2));
t493 = cos(qJ(2));
t457 = t490 * t469 + t493 * t470;
t455 = -t494 * pkin(2) + qJDD(2) * qJ(3) + t457;
t483 = -g(3) + qJDD(1);
t512 = qJD(2) * qJD(3);
t515 = t486 * t483 - 0.2e1 * t484 * t512;
t436 = (-pkin(6) * qJDD(2) + t494 * t519 - t455) * t484 + t515;
t440 = t484 * t483 + (t455 + 0.2e1 * t512) * t486;
t511 = qJDD(2) * t486;
t437 = -pkin(3) * t517 + pkin(6) * t511 + t440;
t489 = sin(qJ(4));
t492 = cos(qJ(4));
t421 = t492 * t436 - t489 * t437;
t499 = t484 * t492 + t486 * t489;
t498 = -t484 * t489 + t486 * t492;
t462 = t498 * qJD(2);
t513 = t462 * qJD(4);
t454 = t499 * qJDD(2) + t513;
t463 = t499 * qJD(2);
t417 = (-t454 + t513) * pkin(7) + (t462 * t463 + qJDD(4)) * pkin(4) + t421;
t422 = t489 * t436 + t492 * t437;
t453 = -t463 * qJD(4) + t498 * qJDD(2);
t460 = qJD(4) * pkin(4) - t463 * pkin(7);
t461 = t462 ^ 2;
t418 = -t461 * pkin(4) + t453 * pkin(7) - qJD(4) * t460 + t422;
t488 = sin(qJ(5));
t491 = cos(qJ(5));
t415 = t491 * t417 - t488 * t418;
t446 = t491 * t462 - t488 * t463;
t426 = t446 * qJD(5) + t488 * t453 + t491 * t454;
t447 = t488 * t462 + t491 * t463;
t432 = -t446 * mrSges(6,1) + t447 * mrSges(6,2);
t482 = qJD(4) + qJD(5);
t441 = -t482 * mrSges(6,2) + t446 * mrSges(6,3);
t479 = qJDD(4) + qJDD(5);
t413 = m(6) * t415 + t479 * mrSges(6,1) - t426 * mrSges(6,3) - t447 * t432 + t482 * t441;
t416 = t488 * t417 + t491 * t418;
t425 = -t447 * qJD(5) + t491 * t453 - t488 * t454;
t442 = t482 * mrSges(6,1) - t447 * mrSges(6,3);
t414 = m(6) * t416 - t479 * mrSges(6,2) + t425 * mrSges(6,3) + t446 * t432 - t482 * t442;
t405 = t491 * t413 + t488 * t414;
t450 = -t462 * mrSges(5,1) + t463 * mrSges(5,2);
t458 = -qJD(4) * mrSges(5,2) + t462 * mrSges(5,3);
t403 = m(5) * t421 + qJDD(4) * mrSges(5,1) - t454 * mrSges(5,3) + qJD(4) * t458 - t463 * t450 + t405;
t459 = qJD(4) * mrSges(5,1) - t463 * mrSges(5,3);
t505 = -t488 * t413 + t491 * t414;
t404 = m(5) * t422 - qJDD(4) * mrSges(5,2) + t453 * mrSges(5,3) - qJD(4) * t459 + t462 * t450 + t505;
t399 = t492 * t403 + t489 * t404;
t439 = -t484 * t455 + t515;
t497 = mrSges(4,3) * qJDD(2) + t494 * (-mrSges(4,1) * t486 + t518);
t397 = m(4) * t439 - t497 * t484 + t399;
t506 = -t489 * t403 + t492 * t404;
t398 = m(4) * t440 + t497 * t486 + t506;
t507 = -t484 * t397 + t486 * t398;
t390 = m(3) * t457 - t494 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t507;
t456 = t493 * t469 - t490 * t470;
t501 = qJDD(3) - t456;
t452 = -qJDD(2) * pkin(2) - t494 * qJ(3) + t501;
t480 = t484 ^ 2;
t438 = (-pkin(2) - t519) * qJDD(2) + (-qJ(3) + (-t480 - t481) * pkin(6)) * t494 + t501;
t420 = -t453 * pkin(4) - t461 * pkin(7) + t463 * t460 + t438;
t500 = m(6) * t420 - t425 * mrSges(6,1) + t426 * mrSges(6,2) - t446 * t441 + t447 * t442;
t496 = m(5) * t438 - t453 * mrSges(5,1) + t454 * mrSges(5,2) - t462 * t458 + t463 * t459 + t500;
t495 = -m(4) * t452 + mrSges(4,1) * t511 - t496 + (t480 * t494 + t517) * mrSges(4,3);
t409 = t495 + (mrSges(3,1) - t518) * qJDD(2) - t494 * mrSges(3,2) + m(3) * t456;
t387 = t490 * t390 + t493 * t409;
t385 = m(2) * t469 + t387;
t508 = t493 * t390 - t490 * t409;
t386 = m(2) * t470 + t508;
t516 = t487 * t385 + t485 * t386;
t391 = t486 * t397 + t484 * t398;
t502 = Ifges(4,5) * t484 + Ifges(4,6) * t486;
t514 = t494 * t502;
t510 = m(3) * t483 + t391;
t509 = -t485 * t385 + t487 * t386;
t504 = Ifges(4,1) * t484 + Ifges(4,4) * t486;
t503 = Ifges(4,4) * t484 + Ifges(4,2) * t486;
t445 = Ifges(5,1) * t463 + Ifges(5,4) * t462 + Ifges(5,5) * qJD(4);
t444 = Ifges(5,4) * t463 + Ifges(5,2) * t462 + Ifges(5,6) * qJD(4);
t443 = Ifges(5,5) * t463 + Ifges(5,6) * t462 + Ifges(5,3) * qJD(4);
t429 = Ifges(6,1) * t447 + Ifges(6,4) * t446 + Ifges(6,5) * t482;
t428 = Ifges(6,4) * t447 + Ifges(6,2) * t446 + Ifges(6,6) * t482;
t427 = Ifges(6,5) * t447 + Ifges(6,6) * t446 + Ifges(6,3) * t482;
t407 = mrSges(6,2) * t420 - mrSges(6,3) * t415 + Ifges(6,1) * t426 + Ifges(6,4) * t425 + Ifges(6,5) * t479 + t446 * t427 - t482 * t428;
t406 = -mrSges(6,1) * t420 + mrSges(6,3) * t416 + Ifges(6,4) * t426 + Ifges(6,2) * t425 + Ifges(6,6) * t479 - t447 * t427 + t482 * t429;
t393 = mrSges(5,2) * t438 - mrSges(5,3) * t421 + Ifges(5,1) * t454 + Ifges(5,4) * t453 + Ifges(5,5) * qJDD(4) - pkin(7) * t405 - qJD(4) * t444 - t488 * t406 + t491 * t407 + t462 * t443;
t392 = -mrSges(5,1) * t438 + mrSges(5,3) * t422 + Ifges(5,4) * t454 + Ifges(5,2) * t453 + Ifges(5,6) * qJDD(4) - pkin(4) * t500 + pkin(7) * t505 + qJD(4) * t445 + t491 * t406 + t488 * t407 - t463 * t443;
t381 = mrSges(4,2) * t452 - mrSges(4,3) * t439 - pkin(6) * t399 + t504 * qJDD(2) - t489 * t392 + t492 * t393 + t486 * t514;
t380 = -mrSges(4,1) * t452 + mrSges(4,3) * t440 - pkin(3) * t496 + pkin(6) * t506 + t503 * qJDD(2) + t492 * t392 + t489 * t393 - t484 * t514;
t379 = -Ifges(5,3) * qJDD(4) - Ifges(6,3) * t479 - mrSges(3,1) * t483 + t462 * t445 - t463 * t444 - t447 * t428 - Ifges(5,6) * t453 - Ifges(5,5) * t454 + mrSges(3,3) * t457 - mrSges(4,1) * t439 + mrSges(4,2) * t440 + t446 * t429 - Ifges(6,5) * t426 - mrSges(5,1) * t421 + mrSges(5,2) * t422 - Ifges(6,6) * t425 - mrSges(6,1) * t415 + mrSges(6,2) * t416 - pkin(4) * t405 - pkin(3) * t399 - pkin(2) * t391 + (Ifges(3,6) - t502) * qJDD(2) + (-t484 * t503 + t486 * t504 + Ifges(3,5)) * t494;
t378 = mrSges(3,2) * t483 - mrSges(3,3) * t456 + Ifges(3,5) * qJDD(2) - t494 * Ifges(3,6) - qJ(3) * t391 - t484 * t380 + t486 * t381;
t377 = mrSges(2,2) * t483 - mrSges(2,3) * t469 - pkin(5) * t387 + t493 * t378 - t490 * t379;
t376 = -mrSges(2,1) * t483 + mrSges(2,3) * t470 - pkin(1) * t510 + pkin(5) * t508 + t490 * t378 + t493 * t379;
t1 = [-m(1) * g(1) + t509; -m(1) * g(2) + t516; -m(1) * g(3) + m(2) * t483 + t510; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t516 - t485 * t376 + t487 * t377; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t509 + t487 * t376 + t485 * t377; pkin(1) * t387 + mrSges(2,1) * t469 - mrSges(2,2) * t470 + t484 * t381 + t486 * t380 + pkin(2) * (-qJDD(2) * t518 + t495) + qJ(3) * t507 + mrSges(3,1) * t456 - mrSges(3,2) * t457 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(3,3) * qJDD(2);];
tauB = t1;
