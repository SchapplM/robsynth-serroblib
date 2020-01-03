% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR4_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR4_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR4_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:31
% EndTime: 2019-12-31 17:02:32
% DurationCPUTime: 1.15s
% Computational Cost: add. (15553->176), mult. (21577->224), div. (0->0), fcn. (12557->8), ass. (0->84)
t465 = qJD(1) + qJD(2);
t461 = t465 ^ 2;
t468 = cos(pkin(7));
t499 = pkin(3) * t468;
t467 = sin(pkin(7));
t498 = mrSges(4,2) * t467;
t464 = t468 ^ 2;
t497 = t461 * t464;
t462 = qJDD(1) + qJDD(2);
t496 = t462 * t468;
t485 = Ifges(4,5) * t467 + Ifges(4,6) * t468;
t495 = t461 * t485;
t471 = sin(qJ(1));
t474 = cos(qJ(1));
t454 = t471 * g(1) - t474 * g(2);
t449 = qJDD(1) * pkin(1) + t454;
t455 = -t474 * g(1) - t471 * g(2);
t475 = qJD(1) ^ 2;
t450 = -t475 * pkin(1) + t455;
t470 = sin(qJ(2));
t473 = cos(qJ(2));
t439 = t470 * t449 + t473 * t450;
t436 = -t461 * pkin(2) + t462 * qJ(3) + t439;
t493 = qJD(3) * t465;
t492 = -t468 * g(3) - 0.2e1 * t467 * t493;
t421 = (-pkin(6) * t462 + t461 * t499 - t436) * t467 + t492;
t425 = -t467 * g(3) + (t436 + 0.2e1 * t493) * t468;
t422 = -pkin(3) * t497 + pkin(6) * t496 + t425;
t469 = sin(qJ(4));
t472 = cos(qJ(4));
t419 = t472 * t421 - t469 * t422;
t481 = -t467 * t469 + t468 * t472;
t442 = t481 * t465;
t482 = t467 * t472 + t468 * t469;
t443 = t482 * t465;
t432 = -t442 * mrSges(5,1) + t443 * mrSges(5,2);
t435 = t442 * qJD(4) + t482 * t462;
t440 = -qJD(4) * mrSges(5,2) + t442 * mrSges(5,3);
t417 = m(5) * t419 + qJDD(4) * mrSges(5,1) - t435 * mrSges(5,3) + qJD(4) * t440 - t443 * t432;
t420 = t469 * t421 + t472 * t422;
t434 = -t443 * qJD(4) + t481 * t462;
t441 = qJD(4) * mrSges(5,1) - t443 * mrSges(5,3);
t418 = m(5) * t420 - qJDD(4) * mrSges(5,2) + t434 * mrSges(5,3) - qJD(4) * t441 + t442 * t432;
t407 = t472 * t417 + t469 * t418;
t424 = -t467 * t436 + t492;
t483 = mrSges(4,3) * t462 + (-mrSges(4,1) * t468 + t498) * t461;
t405 = m(4) * t424 - t483 * t467 + t407;
t488 = -t469 * t417 + t472 * t418;
t406 = m(4) * t425 + t483 * t468 + t488;
t489 = -t467 * t405 + t468 * t406;
t399 = m(3) * t439 - t461 * mrSges(3,1) - t462 * mrSges(3,2) + t489;
t438 = t473 * t449 - t470 * t450;
t484 = qJDD(3) - t438;
t433 = -t462 * pkin(2) - t461 * qJ(3) + t484;
t463 = t467 ^ 2;
t423 = (-pkin(2) - t499) * t462 + (-qJ(3) + (-t463 - t464) * pkin(6)) * t461 + t484;
t479 = m(5) * t423 - t434 * mrSges(5,1) + t435 * mrSges(5,2) - t442 * t440 + t443 * t441;
t477 = -m(4) * t433 + mrSges(4,1) * t496 - t479 + (t461 * t463 + t497) * mrSges(4,3);
t411 = m(3) * t438 - t461 * mrSges(3,2) + (mrSges(3,1) - t498) * t462 + t477;
t394 = t470 * t399 + t473 * t411;
t391 = m(2) * t454 + qJDD(1) * mrSges(2,1) - t475 * mrSges(2,2) + t394;
t490 = t473 * t399 - t470 * t411;
t392 = m(2) * t455 - t475 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t490;
t494 = t474 * t391 + t471 * t392;
t401 = t468 * t405 + t467 * t406;
t491 = -t471 * t391 + t474 * t392;
t487 = Ifges(4,1) * t467 + Ifges(4,4) * t468;
t486 = Ifges(4,4) * t467 + Ifges(4,2) * t468;
t426 = Ifges(5,5) * t443 + Ifges(5,6) * t442 + Ifges(5,3) * qJD(4);
t428 = Ifges(5,1) * t443 + Ifges(5,4) * t442 + Ifges(5,5) * qJD(4);
t408 = -mrSges(5,1) * t423 + mrSges(5,3) * t420 + Ifges(5,4) * t435 + Ifges(5,2) * t434 + Ifges(5,6) * qJDD(4) + qJD(4) * t428 - t443 * t426;
t427 = Ifges(5,4) * t443 + Ifges(5,2) * t442 + Ifges(5,6) * qJD(4);
t409 = mrSges(5,2) * t423 - mrSges(5,3) * t419 + Ifges(5,1) * t435 + Ifges(5,4) * t434 + Ifges(5,5) * qJDD(4) - qJD(4) * t427 + t442 * t426;
t387 = -mrSges(4,1) * t433 + mrSges(4,3) * t425 - pkin(3) * t479 + pkin(6) * t488 + t472 * t408 + t469 * t409 + t486 * t462 - t467 * t495;
t396 = mrSges(4,2) * t433 - mrSges(4,3) * t424 - pkin(6) * t407 - t469 * t408 + t472 * t409 + t487 * t462 + t468 * t495;
t413 = t462 * t498 - t477;
t480 = mrSges(3,1) * t438 - mrSges(3,2) * t439 + Ifges(3,3) * t462 - pkin(2) * t413 + qJ(3) * t489 + t468 * t387 + t467 * t396;
t478 = mrSges(2,1) * t454 - mrSges(2,2) * t455 + Ifges(2,3) * qJDD(1) + pkin(1) * t394 + t480;
t476 = mrSges(5,1) * t419 - mrSges(5,2) * t420 + Ifges(5,5) * t435 + Ifges(5,6) * t434 + Ifges(5,3) * qJDD(4) + t443 * t427 - t442 * t428;
t385 = mrSges(3,1) * g(3) - mrSges(4,1) * t424 + mrSges(4,2) * t425 + mrSges(3,3) * t439 - pkin(2) * t401 - pkin(3) * t407 + (Ifges(3,6) - t485) * t462 - t476 + (-t467 * t486 + t468 * t487 + Ifges(3,5)) * t461;
t384 = -mrSges(3,2) * g(3) - mrSges(3,3) * t438 + Ifges(3,5) * t462 - t461 * Ifges(3,6) - qJ(3) * t401 - t467 * t387 + t468 * t396;
t383 = -mrSges(2,2) * g(3) - mrSges(2,3) * t454 + Ifges(2,5) * qJDD(1) - t475 * Ifges(2,6) - pkin(5) * t394 + t473 * t384 - t470 * t385;
t382 = Ifges(2,6) * qJDD(1) + t475 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t455 + t470 * t384 + t473 * t385 - pkin(1) * (-m(3) * g(3) + t401) + pkin(5) * t490;
t1 = [-m(1) * g(1) + t491; -m(1) * g(2) + t494; (-m(1) - m(2) - m(3)) * g(3) + t401; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t494 - t471 * t382 + t474 * t383; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t491 + t474 * t382 + t471 * t383; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t478; t478; t480; t413; t476;];
tauJB = t1;
