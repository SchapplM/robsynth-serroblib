% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRR5
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:50
% EndTime: 2019-12-05 15:53:56
% DurationCPUTime: 3.97s
% Computational Cost: add. (41318->238), mult. (93831->302), div. (0->0), fcn. (66780->10), ass. (0->105)
t495 = qJD(2) ^ 2;
t488 = cos(pkin(9));
t520 = pkin(3) * t488;
t486 = sin(pkin(9));
t519 = mrSges(4,2) * t486;
t518 = cos(pkin(8));
t483 = t488 ^ 2;
t517 = t483 * t495;
t487 = sin(pkin(8));
t473 = -t518 * g(1) - t487 * g(2);
t485 = -g(3) + qJDD(1);
t491 = sin(qJ(2));
t494 = cos(qJ(2));
t462 = t494 * t473 + t491 * t485;
t457 = -t495 * pkin(2) + qJDD(2) * qJ(3) + t462;
t472 = t487 * g(1) - t518 * g(2);
t512 = qJD(2) * qJD(3);
t515 = -t488 * t472 - 0.2e1 * t486 * t512;
t438 = (-pkin(6) * qJDD(2) + t495 * t520 - t457) * t486 + t515;
t441 = -t486 * t472 + (t457 + 0.2e1 * t512) * t488;
t511 = qJDD(2) * t488;
t439 = -pkin(3) * t517 + pkin(6) * t511 + t441;
t490 = sin(qJ(4));
t493 = cos(qJ(4));
t423 = t493 * t438 - t490 * t439;
t500 = t486 * t493 + t488 * t490;
t499 = -t486 * t490 + t488 * t493;
t464 = t499 * qJD(2);
t513 = t464 * qJD(4);
t454 = t500 * qJDD(2) + t513;
t465 = t500 * qJD(2);
t419 = (-t454 + t513) * pkin(7) + (t464 * t465 + qJDD(4)) * pkin(4) + t423;
t424 = t490 * t438 + t493 * t439;
t453 = -t465 * qJD(4) + t499 * qJDD(2);
t460 = qJD(4) * pkin(4) - t465 * pkin(7);
t463 = t464 ^ 2;
t420 = -t463 * pkin(4) + t453 * pkin(7) - qJD(4) * t460 + t424;
t489 = sin(qJ(5));
t492 = cos(qJ(5));
t417 = t492 * t419 - t489 * t420;
t448 = t492 * t464 - t489 * t465;
t428 = t448 * qJD(5) + t489 * t453 + t492 * t454;
t449 = t489 * t464 + t492 * t465;
t434 = -t448 * mrSges(6,1) + t449 * mrSges(6,2);
t484 = qJD(4) + qJD(5);
t442 = -t484 * mrSges(6,2) + t448 * mrSges(6,3);
t481 = qJDD(4) + qJDD(5);
t415 = m(6) * t417 + t481 * mrSges(6,1) - t428 * mrSges(6,3) - t449 * t434 + t484 * t442;
t418 = t489 * t419 + t492 * t420;
t427 = -t449 * qJD(5) + t492 * t453 - t489 * t454;
t443 = t484 * mrSges(6,1) - t449 * mrSges(6,3);
t416 = m(6) * t418 - t481 * mrSges(6,2) + t427 * mrSges(6,3) + t448 * t434 - t484 * t443;
t407 = t492 * t415 + t489 * t416;
t451 = -t464 * mrSges(5,1) + t465 * mrSges(5,2);
t458 = -qJD(4) * mrSges(5,2) + t464 * mrSges(5,3);
t405 = m(5) * t423 + qJDD(4) * mrSges(5,1) - t454 * mrSges(5,3) + qJD(4) * t458 - t465 * t451 + t407;
t459 = qJD(4) * mrSges(5,1) - t465 * mrSges(5,3);
t506 = -t489 * t415 + t492 * t416;
t406 = m(5) * t424 - qJDD(4) * mrSges(5,2) + t453 * mrSges(5,3) - qJD(4) * t459 + t464 * t451 + t506;
t401 = t493 * t405 + t490 * t406;
t440 = -t486 * t457 + t515;
t498 = mrSges(4,3) * qJDD(2) + t495 * (-mrSges(4,1) * t488 + t519);
t399 = m(4) * t440 - t498 * t486 + t401;
t507 = -t490 * t405 + t493 * t406;
t400 = m(4) * t441 + t498 * t488 + t507;
t508 = -t486 * t399 + t488 * t400;
t392 = m(3) * t462 - t495 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t508;
t461 = -t491 * t473 + t494 * t485;
t502 = qJDD(3) - t461;
t456 = -qJDD(2) * pkin(2) - t495 * qJ(3) + t502;
t482 = t486 ^ 2;
t444 = (-pkin(2) - t520) * qJDD(2) + (-qJ(3) + (-t482 - t483) * pkin(6)) * t495 + t502;
t422 = -t453 * pkin(4) - t463 * pkin(7) + t465 * t460 + t444;
t501 = m(6) * t422 - t427 * mrSges(6,1) + t428 * mrSges(6,2) - t448 * t442 + t449 * t443;
t497 = m(5) * t444 - t453 * mrSges(5,1) + t454 * mrSges(5,2) - t464 * t458 + t465 * t459 + t501;
t496 = -m(4) * t456 + mrSges(4,1) * t511 - t497 + (t482 * t495 + t517) * mrSges(4,3);
t411 = t496 + (mrSges(3,1) - t519) * qJDD(2) - t495 * mrSges(3,2) + m(3) * t461;
t509 = t494 * t392 - t491 * t411;
t388 = m(2) * t473 + t509;
t395 = t488 * t399 + t486 * t400;
t394 = (m(2) + m(3)) * t472 - t395;
t516 = t487 * t388 + t518 * t394;
t389 = t491 * t392 + t494 * t411;
t503 = Ifges(4,5) * t486 + Ifges(4,6) * t488;
t514 = t495 * t503;
t510 = t518 * t388 - t487 * t394;
t505 = Ifges(4,1) * t486 + Ifges(4,4) * t488;
t504 = Ifges(4,4) * t486 + Ifges(4,2) * t488;
t447 = Ifges(5,1) * t465 + Ifges(5,4) * t464 + Ifges(5,5) * qJD(4);
t446 = Ifges(5,4) * t465 + Ifges(5,2) * t464 + Ifges(5,6) * qJD(4);
t445 = Ifges(5,5) * t465 + Ifges(5,6) * t464 + Ifges(5,3) * qJD(4);
t431 = Ifges(6,1) * t449 + Ifges(6,4) * t448 + Ifges(6,5) * t484;
t430 = Ifges(6,4) * t449 + Ifges(6,2) * t448 + Ifges(6,6) * t484;
t429 = Ifges(6,5) * t449 + Ifges(6,6) * t448 + Ifges(6,3) * t484;
t409 = mrSges(6,2) * t422 - mrSges(6,3) * t417 + Ifges(6,1) * t428 + Ifges(6,4) * t427 + Ifges(6,5) * t481 + t448 * t429 - t484 * t430;
t408 = -mrSges(6,1) * t422 + mrSges(6,3) * t418 + Ifges(6,4) * t428 + Ifges(6,2) * t427 + Ifges(6,6) * t481 - t449 * t429 + t484 * t431;
t397 = mrSges(5,2) * t444 - mrSges(5,3) * t423 + Ifges(5,1) * t454 + Ifges(5,4) * t453 + Ifges(5,5) * qJDD(4) - pkin(7) * t407 - qJD(4) * t446 - t489 * t408 + t492 * t409 + t464 * t445;
t396 = -mrSges(5,1) * t444 + mrSges(5,3) * t424 + Ifges(5,4) * t454 + Ifges(5,2) * t453 + Ifges(5,6) * qJDD(4) - pkin(4) * t501 + pkin(7) * t506 + qJD(4) * t447 + t492 * t408 + t489 * t409 - t465 * t445;
t385 = mrSges(4,2) * t456 - mrSges(4,3) * t440 - pkin(6) * t401 + t505 * qJDD(2) - t490 * t396 + t493 * t397 + t488 * t514;
t384 = -mrSges(4,1) * t456 + mrSges(4,3) * t441 - pkin(3) * t497 + pkin(6) * t507 + t504 * qJDD(2) + t493 * t396 + t490 * t397 - t486 * t514;
t383 = -Ifges(5,3) * qJDD(4) + mrSges(3,1) * t472 - Ifges(6,3) * t481 - t465 * t446 + mrSges(3,3) * t462 + t464 * t447 - Ifges(5,6) * t453 - Ifges(5,5) * t454 + t448 * t431 - t449 * t430 - mrSges(4,1) * t440 + mrSges(4,2) * t441 - mrSges(5,1) * t423 + mrSges(5,2) * t424 - Ifges(6,6) * t427 - Ifges(6,5) * t428 - mrSges(6,1) * t417 + mrSges(6,2) * t418 - pkin(4) * t407 - pkin(3) * t401 - pkin(2) * t395 + (Ifges(3,6) - t503) * qJDD(2) + (-t486 * t504 + t488 * t505 + Ifges(3,5)) * t495;
t382 = -mrSges(3,2) * t472 - mrSges(3,3) * t461 + Ifges(3,5) * qJDD(2) - t495 * Ifges(3,6) - qJ(3) * t395 - t486 * t384 + t488 * t385;
t381 = -mrSges(2,1) * t485 + mrSges(2,3) * t473 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t461 + mrSges(3,2) * t462 - t486 * t385 - t488 * t384 - pkin(2) * (-qJDD(2) * t519 + t496) - qJ(3) * t508 - pkin(1) * t389;
t380 = mrSges(2,2) * t485 - mrSges(2,3) * t472 - pkin(5) * t389 + t494 * t382 - t491 * t383;
t1 = [-m(1) * g(1) + t510; -m(1) * g(2) + t516; -m(1) * g(3) + m(2) * t485 + t389; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t516 + t518 * t380 - t487 * t381; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t510 + t487 * t380 + t518 * t381; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t472 - mrSges(2,2) * t473 + t491 * t382 + t494 * t383 + pkin(1) * (m(3) * t472 - t395) + pkin(5) * t509;];
tauB = t1;
