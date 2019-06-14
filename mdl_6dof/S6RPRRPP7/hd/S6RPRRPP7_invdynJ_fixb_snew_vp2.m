% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPP7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:46:29
% EndTime: 2019-05-05 21:46:31
% DurationCPUTime: 1.32s
% Computational Cost: add. (6330->243), mult. (11967->279), div. (0->0), fcn. (6743->6), ass. (0->95)
t512 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t493 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t492 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t511 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t491 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t510 = Ifges(5,3) + Ifges(6,2) + Ifges(7,3);
t467 = sin(qJ(1));
t469 = cos(qJ(1));
t479 = -t469 * g(1) - t467 * g(2);
t509 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t479;
t468 = cos(qJ(3));
t496 = t468 * qJD(1);
t471 = qJD(1) ^ 2;
t503 = (-pkin(1) - pkin(7));
t435 = (t503 * t471) - t509;
t466 = sin(qJ(3));
t495 = qJD(1) * qJD(3);
t484 = t468 * t495;
t452 = -t466 * qJDD(1) - t484;
t485 = t466 * t495;
t453 = t468 * qJDD(1) - t485;
t388 = (-t453 + t485) * pkin(8) + (-t452 + t484) * pkin(3) + t435;
t483 = t467 * g(1) - t469 * g(2);
t475 = -t471 * qJ(2) + qJDD(2) - t483;
t436 = t503 * qJDD(1) + t475;
t425 = -t468 * g(3) + t466 * t436;
t451 = (t466 * pkin(3) - t468 * pkin(8)) * qJD(1);
t470 = qJD(3) ^ 2;
t497 = t466 * qJD(1);
t392 = -t470 * pkin(3) + qJDD(3) * pkin(8) - t451 * t497 + t425;
t465 = sin(qJ(4));
t502 = cos(qJ(4));
t385 = t502 * t388 - t465 * t392;
t448 = -t502 * qJD(3) + t465 * t496;
t449 = t465 * qJD(3) + t502 * t496;
t420 = t448 * pkin(4) - t449 * qJ(5);
t447 = qJDD(4) - t452;
t458 = qJD(4) + t497;
t457 = t458 ^ 2;
t383 = -t447 * pkin(4) - t457 * qJ(5) + t449 * t420 + qJDD(5) - t385;
t432 = -t448 * mrSges(6,2) + t458 * mrSges(6,3);
t508 = -m(6) * t383 + t447 * mrSges(6,1) + t458 * t432;
t417 = -t448 * qJD(4) + t465 * qJDD(3) + t502 * t453;
t424 = t466 * g(3) + t468 * t436;
t474 = qJDD(3) * pkin(3) + t470 * pkin(8) - t451 * t496 + t424;
t499 = t448 * t458;
t507 = (-t417 + t499) * qJ(5) - t474;
t426 = t458 * mrSges(7,2) + t448 * mrSges(7,3);
t505 = -0.2e1 * t449;
t376 = qJD(6) * t505 + (-t417 - t499) * qJ(6) + (t448 * t449 - t447) * pkin(5) + t383;
t422 = -t448 * mrSges(7,1) + t449 * mrSges(7,2);
t480 = -m(7) * t376 + t417 * mrSges(7,3) + t449 * t422;
t374 = -t447 * mrSges(7,1) - t458 * t426 - t480;
t421 = t448 * mrSges(6,1) - t449 * mrSges(6,3);
t372 = t417 * mrSges(6,2) + t449 * t421 + t374 - t508;
t386 = t465 * t388 + t502 * t392;
t504 = 2 * qJD(5);
t382 = -t457 * pkin(4) + t447 * qJ(5) - t448 * t420 + t458 * t504 + t386;
t416 = t449 * qJD(4) - t502 * qJDD(3) + t465 * t453;
t428 = -t458 * pkin(5) - t449 * qJ(6);
t446 = t448 ^ 2;
t378 = -t446 * pkin(5) + t416 * qJ(6) + 0.2e1 * qJD(6) * t448 + t458 * t428 + t382;
t429 = -t458 * mrSges(7,1) - t449 * mrSges(7,3);
t431 = -t458 * mrSges(6,1) + t449 * mrSges(6,2);
t489 = m(7) * t378 + t416 * mrSges(7,3) + t448 * t422;
t476 = m(6) * t382 + t447 * mrSges(6,3) + t458 * t431 + t489;
t486 = -t493 * t448 + t512 * t449 + t492 * t458;
t487 = t511 * t448 + t493 * t449 + t491 * t458;
t506 = -t491 * t416 + t492 * t417 + t510 * t447 + t486 * t448 + t487 * t449 + mrSges(5,1) * t385 - mrSges(6,1) * t383 - mrSges(7,1) * t376 - mrSges(5,2) * t386 + mrSges(7,2) * t378 + mrSges(6,3) * t382 - pkin(4) * t372 - pkin(5) * t374 + qJ(5) * (-t416 * mrSges(6,2) + t447 * mrSges(7,2) - t448 * t421 + t458 * t429 + t476);
t500 = -mrSges(5,3) - mrSges(6,2);
t427 = -t458 * mrSges(5,2) - t448 * mrSges(5,3);
t498 = -t448 * mrSges(5,1) - t449 * mrSges(5,2) - t421;
t369 = m(5) * t385 + (t426 + t427) * t458 + t498 * t449 + (mrSges(5,1) + mrSges(7,1)) * t447 + t500 * t417 + t480 + t508;
t430 = t458 * mrSges(5,1) - t449 * mrSges(5,3);
t370 = m(5) * t386 + (t429 - t430) * t458 + t498 * t448 + (-mrSges(5,2) + mrSges(7,2)) * t447 + t500 * t416 + t476;
t364 = t502 * t369 + t465 * t370;
t488 = t491 * t448 - t492 * t449 - t510 * t458;
t482 = -t465 * t369 + t502 * t370;
t380 = -t446 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t416 + (-pkin(4) * t458 + t428 + t504) * t449 - t507;
t375 = m(7) * t380 - t416 * mrSges(7,1) + t417 * mrSges(7,2) - t448 * t426 + t449 * t429;
t450 = (t466 * mrSges(4,1) + t468 * mrSges(4,2)) * qJD(1);
t454 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t497;
t455 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t496;
t384 = qJD(5) * t505 + (t449 * t458 + t416) * pkin(4) + t507;
t371 = m(6) * t384 + t416 * mrSges(6,1) - t417 * mrSges(6,3) - t449 * t431 + t448 * t432 - t375;
t472 = m(5) * t474 - t416 * mrSges(5,1) - t417 * mrSges(5,2) - t448 * t427 - t449 * t430 - t371;
t478 = t466 * (m(4) * t425 - qJDD(3) * mrSges(4,2) + t452 * mrSges(4,3) - qJD(3) * t455 - t450 * t497 + t482) + t468 * (m(4) * t424 + qJDD(3) * mrSges(4,1) - t453 * mrSges(4,3) + qJD(3) * t454 - t450 * t496 + t472);
t442 = (Ifges(4,5) * qJD(3)) + (t468 * Ifges(4,1) - t466 * Ifges(4,4)) * qJD(1);
t441 = (Ifges(4,6) * qJD(3)) + (t468 * Ifges(4,4) - t466 * Ifges(4,2)) * qJD(1);
t438 = -qJDD(1) * pkin(1) + t475;
t437 = t471 * pkin(1) + t509;
t362 = -mrSges(5,2) * t474 + mrSges(6,2) * t383 + mrSges(7,2) * t380 - mrSges(5,3) * t385 - mrSges(6,3) * t384 - mrSges(7,3) * t376 - qJ(5) * t371 - qJ(6) * t374 - t493 * t416 + t512 * t417 + t492 * t447 + t488 * t448 - t487 * t458;
t361 = m(3) * t438 + qJDD(1) * mrSges(3,2) - (t471 * mrSges(3,3)) + t478;
t360 = mrSges(5,1) * t474 + mrSges(5,3) * t386 - mrSges(6,1) * t384 + mrSges(6,2) * t382 + mrSges(7,1) * t380 - mrSges(7,3) * t378 + pkin(5) * t375 - qJ(6) * t489 - pkin(4) * t371 + (-qJ(6) * t429 + t486) * t458 + t488 * t449 + (-qJ(6) * mrSges(7,2) + t491) * t447 + t493 * t417 + t511 * t416;
t1 = [mrSges(2,1) * t483 - mrSges(2,2) * t479 + mrSges(3,2) * t438 - mrSges(3,3) * t437 + t468 * (mrSges(4,2) * t435 - mrSges(4,3) * t424 + Ifges(4,1) * t453 + Ifges(4,4) * t452 + Ifges(4,5) * qJDD(3) - pkin(8) * t364 - qJD(3) * t441 - t465 * t360 + t502 * t362) - pkin(7) * t478 - pkin(1) * t361 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t437 + m(4) * t435 - t452 * mrSges(4,1) + t471 * mrSges(3,2) + t453 * mrSges(4,2) + qJDD(1) * mrSges(3,3) + t455 * t496 + t364) * qJ(2) + (qJ(2) * t454 * qJD(1) + mrSges(4,1) * t435 - mrSges(4,3) * t425 - Ifges(4,4) * t453 - Ifges(4,2) * t452 - Ifges(4,6) * qJDD(3) + pkin(3) * t364 - qJD(3) * t442 + t506) * t466; t361; Ifges(4,5) * t453 + Ifges(4,6) * t452 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t424 - mrSges(4,2) * t425 + t465 * t362 + t502 * t360 + pkin(3) * t472 + pkin(8) * t482 + (t468 * t441 + t466 * t442) * qJD(1); t506; t372; t375;];
tauJ  = t1;
