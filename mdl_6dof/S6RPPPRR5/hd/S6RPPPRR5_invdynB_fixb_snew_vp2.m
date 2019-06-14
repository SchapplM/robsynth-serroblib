% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2019-05-05 13:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPPRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:50:16
% EndTime: 2019-05-05 13:50:19
% DurationCPUTime: 1.73s
% Computational Cost: add. (18480->251), mult. (31061->289), div. (0->0), fcn. (12740->8), ass. (0->100)
t502 = -2 * qJD(1);
t468 = sin(qJ(1));
t471 = cos(qJ(1));
t449 = -t471 * g(1) - t468 * g(2);
t501 = -qJDD(1) * qJ(2) + (qJD(2) * t502) - t449;
t500 = -m(3) - m(4);
t499 = mrSges(2,1) - mrSges(3,2);
t498 = -pkin(1) - qJ(3);
t460 = -g(3) + qJDD(4);
t470 = cos(qJ(5));
t497 = t470 * t460;
t473 = qJD(1) ^ 2;
t496 = t473 * qJ(2);
t432 = (pkin(1) * t473) + t501;
t448 = t468 * g(1) - g(2) * t471;
t480 = qJDD(2) - t448;
t475 = (qJD(3) * t502) + qJDD(1) * t498 + t480;
t425 = (-pkin(3) - qJ(2)) * t473 + t475;
t429 = t473 * t498 + qJDD(3) - t501;
t426 = qJDD(1) * pkin(3) + t429;
t464 = sin(pkin(9));
t465 = cos(pkin(9));
t413 = t425 * t465 + t426 * t464;
t411 = -(pkin(4) * t473) + qJDD(1) * pkin(7) + t413;
t467 = sin(qJ(5));
t408 = t411 * t470 + t460 * t467;
t442 = (-mrSges(6,1) * t470 + mrSges(6,2) * t467) * qJD(1);
t491 = qJD(1) * qJD(5);
t486 = t467 * t491;
t445 = qJDD(1) * t470 - t486;
t493 = qJD(1) * t467;
t446 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t493;
t412 = -t425 * t464 + t465 * t426;
t410 = -qJDD(1) * pkin(4) - t473 * pkin(7) - t412;
t485 = t470 * t491;
t444 = qJDD(1) * t467 + t485;
t404 = (-t444 - t485) * pkin(8) + (-t445 + t486) * pkin(5) + t410;
t443 = (-pkin(5) * t470 - pkin(8) * t467) * qJD(1);
t472 = qJD(5) ^ 2;
t492 = t470 * qJD(1);
t406 = -pkin(5) * t472 + qJDD(5) * pkin(8) + t443 * t492 + t408;
t466 = sin(qJ(6));
t469 = cos(qJ(6));
t402 = t404 * t469 - t406 * t466;
t440 = qJD(5) * t469 - t466 * t493;
t420 = qJD(6) * t440 + qJDD(5) * t466 + t444 * t469;
t441 = qJD(5) * t466 + t469 * t493;
t424 = -mrSges(7,1) * t440 + mrSges(7,2) * t441;
t450 = qJD(6) - t492;
t430 = -mrSges(7,2) * t450 + mrSges(7,3) * t440;
t439 = qJDD(6) - t445;
t400 = m(7) * t402 + mrSges(7,1) * t439 - mrSges(7,3) * t420 - t424 * t441 + t430 * t450;
t403 = t404 * t466 + t406 * t469;
t419 = -qJD(6) * t441 + qJDD(5) * t469 - t444 * t466;
t431 = mrSges(7,1) * t450 - mrSges(7,3) * t441;
t401 = m(7) * t403 - mrSges(7,2) * t439 + mrSges(7,3) * t419 + t424 * t440 - t431 * t450;
t482 = -t400 * t466 + t401 * t469;
t393 = m(6) * t408 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t445 - qJD(5) * t446 + t442 * t492 + t482;
t407 = -t411 * t467 + t497;
t447 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t492;
t405 = -qJDD(5) * pkin(5) - t472 * pkin(8) - t497 + (qJD(1) * t443 + t411) * t467;
t477 = -m(7) * t405 + mrSges(7,1) * t419 - mrSges(7,2) * t420 + t430 * t440 - t431 * t441;
t398 = m(6) * t407 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t444 + qJD(5) * t447 - t442 * t493 + t477;
t483 = t393 * t470 - t398 * t467;
t385 = m(5) * t413 - (mrSges(5,1) * t473) - qJDD(1) * mrSges(5,2) + t483;
t394 = t400 * t469 + t401 * t466;
t474 = -m(6) * t410 + mrSges(6,1) * t445 - mrSges(6,2) * t444 - t446 * t493 + t447 * t492 - t394;
t390 = m(5) * t412 + qJDD(1) * mrSges(5,1) - mrSges(5,2) * t473 + t474;
t380 = t385 * t464 + t390 * t465;
t481 = -m(4) * t429 - qJDD(1) * mrSges(4,1) - t380;
t478 = -m(3) * t432 + (mrSges(3,2) * t473) + qJDD(1) * mrSges(3,3) - t481;
t377 = m(2) * t449 - qJDD(1) * mrSges(2,2) + ((-mrSges(2,1) - mrSges(4,3)) * t473) + t478;
t428 = t475 - t496;
t494 = t385 * t465 - t390 * t464;
t379 = m(4) * t428 - mrSges(4,1) * t473 - qJDD(1) * mrSges(4,3) + t494;
t433 = -qJDD(1) * pkin(1) + t480 - t496;
t476 = -m(3) * t433 + mrSges(3,3) * t473 - t379;
t378 = m(2) * t448 - t473 * mrSges(2,2) + qJDD(1) * t499 + t476;
t495 = t377 * t468 + t378 * t471;
t387 = t393 * t467 + t398 * t470;
t489 = -Ifges(3,4) + Ifges(2,5) - Ifges(4,6);
t488 = (Ifges(4,4) - Ifges(3,5) + Ifges(2,6));
t487 = -m(5) * t460 - t387;
t484 = t377 * t471 - t378 * t468;
t436 = (Ifges(6,5) * qJD(5)) + (Ifges(6,1) * t467 + Ifges(6,4) * t470) * qJD(1);
t435 = (Ifges(6,6) * qJD(5)) + (Ifges(6,4) * t467 + Ifges(6,2) * t470) * qJD(1);
t434 = (Ifges(6,3) * qJD(5)) + (Ifges(6,5) * t467 + Ifges(6,6) * t470) * qJD(1);
t416 = Ifges(7,1) * t441 + Ifges(7,4) * t440 + Ifges(7,5) * t450;
t415 = Ifges(7,4) * t441 + Ifges(7,2) * t440 + Ifges(7,6) * t450;
t414 = Ifges(7,5) * t441 + Ifges(7,6) * t440 + Ifges(7,3) * t450;
t396 = mrSges(7,2) * t405 - mrSges(7,3) * t402 + Ifges(7,1) * t420 + Ifges(7,4) * t419 + Ifges(7,5) * t439 + t414 * t440 - t415 * t450;
t395 = -mrSges(7,1) * t405 + mrSges(7,3) * t403 + Ifges(7,4) * t420 + Ifges(7,2) * t419 + Ifges(7,6) * t439 - t414 * t441 + t416 * t450;
t386 = g(3) * t500 - t487;
t382 = -mrSges(6,1) * t410 - mrSges(7,1) * t402 + mrSges(7,2) * t403 + mrSges(6,3) * t408 + Ifges(6,4) * t444 - Ifges(7,5) * t420 + Ifges(6,2) * t445 + Ifges(6,6) * qJDD(5) - Ifges(7,6) * t419 - Ifges(7,3) * t439 - pkin(5) * t394 + qJD(5) * t436 - t415 * t441 + t416 * t440 - t434 * t493;
t381 = mrSges(6,2) * t410 - mrSges(6,3) * t407 + Ifges(6,1) * t444 + Ifges(6,4) * t445 + Ifges(6,5) * qJDD(5) - pkin(8) * t394 - qJD(5) * t435 - t395 * t466 + t396 * t469 + t434 * t492;
t373 = Ifges(5,6) * qJDD(1) + (t473 * Ifges(5,5)) - mrSges(5,1) * t460 + mrSges(5,3) * t413 - Ifges(6,5) * t444 - Ifges(6,6) * t445 - Ifges(6,3) * qJDD(5) - mrSges(6,1) * t407 + mrSges(6,2) * t408 - t466 * t396 - t469 * t395 - pkin(5) * t477 - pkin(8) * t482 - pkin(4) * t387 + (-t435 * t467 + t436 * t470) * qJD(1);
t372 = mrSges(5,2) * t460 - mrSges(5,3) * t412 + Ifges(5,5) * qJDD(1) - Ifges(5,6) * t473 - pkin(7) * t387 + t381 * t470 - t382 * t467;
t371 = -qJ(2) * t386 - mrSges(2,3) * t448 + pkin(2) * t379 + mrSges(3,1) * t433 + t464 * t372 + t465 * t373 + pkin(3) * t487 + qJ(4) * t494 - mrSges(4,2) * t428 - (t488 * t473) + t489 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,1)) * g(3);
t370 = -pkin(1) * t386 + mrSges(2,3) * t449 - pkin(2) * t481 + qJ(3) * t487 - t465 * t372 + t464 * t373 + qJ(4) * t380 - mrSges(3,1) * t432 - mrSges(4,2) * t429 + (-mrSges(4,3) * pkin(2) + t489) * t473 + t488 * qJDD(1) + (m(4) * qJ(3) + mrSges(4,3) + t499) * g(3);
t1 = [-m(1) * g(1) + t484; -m(1) * g(2) + t495; (-m(1) - m(2) + t500) * g(3) - t487; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t495 - t370 * t468 + t371 * t471; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t484 + t471 * t370 + t468 * t371; pkin(1) * t476 - mrSges(2,2) * t449 + mrSges(2,1) * t448 + qJ(2) * (-(mrSges(4,3) * t473) + t478) - qJ(3) * t379 - mrSges(3,3) * t432 + mrSges(3,2) * t433 + pkin(3) * t380 + mrSges(4,1) * t429 - mrSges(4,3) * t428 + t467 * t381 + t470 * t382 + pkin(4) * t474 + pkin(7) * t483 + mrSges(5,1) * t412 - mrSges(5,2) * t413 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + Ifges(5,3)) * qJDD(1);];
tauB  = t1;
