% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPRP3
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP3_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP3_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP3_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:37
% EndTime: 2019-12-31 16:42:38
% DurationCPUTime: 0.80s
% Computational Cost: add. (5003->175), mult. (9472->210), div. (0->0), fcn. (4252->6), ass. (0->77)
t507 = Ifges(4,1) + Ifges(5,1);
t499 = Ifges(4,4) + Ifges(5,4);
t498 = Ifges(4,5) + Ifges(5,5);
t506 = Ifges(4,2) + Ifges(5,2);
t497 = Ifges(4,6) + Ifges(5,6);
t505 = Ifges(4,3) + Ifges(5,3);
t471 = sin(qJ(3));
t473 = cos(qJ(3));
t446 = (-mrSges(5,1) * t473 + mrSges(5,2) * t471) * qJD(1);
t488 = qJD(1) * qJD(3);
t484 = t473 * t488;
t449 = t471 * qJDD(1) + t484;
t472 = sin(qJ(1));
t474 = cos(qJ(1));
t458 = t472 * g(1) - t474 * g(2);
t445 = qJDD(1) * pkin(1) + t458;
t459 = -t474 * g(1) - t472 * g(2);
t475 = qJD(1) ^ 2;
t448 = -t475 * pkin(1) + t459;
t469 = sin(pkin(6));
t470 = cos(pkin(6));
t429 = t469 * t445 + t470 * t448;
t426 = -t475 * pkin(2) + qJDD(1) * pkin(5) + t429;
t468 = -g(3) + qJDD(2);
t461 = t473 * t468;
t487 = qJD(1) * qJD(4);
t501 = pkin(3) * t475;
t419 = qJDD(3) * pkin(3) + t461 + (-t449 + t484) * qJ(4) + (t473 * t501 - t426 - 0.2e1 * t487) * t471;
t489 = qJD(1) * t473;
t456 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t489;
t486 = m(5) * t419 + qJDD(3) * mrSges(5,1) + qJD(3) * t456;
t490 = qJD(1) * t471;
t415 = -t449 * mrSges(5,3) - t446 * t490 + t486;
t423 = t473 * t426 + t471 * t468;
t450 = t473 * qJDD(1) - t471 * t488;
t453 = qJD(3) * pkin(3) - qJ(4) * t490;
t467 = t473 ^ 2;
t420 = t450 * qJ(4) - qJD(3) * t453 - t467 * t501 + 0.2e1 * t473 * t487 + t423;
t422 = -t471 * t426 + t461;
t492 = t498 * qJD(3) + (t507 * t471 + t499 * t473) * qJD(1);
t493 = t497 * qJD(3) + (t499 * t471 + t506 * t473) * qJD(1);
t504 = mrSges(4,1) * t422 + mrSges(5,1) * t419 - mrSges(4,2) * t423 - mrSges(5,2) * t420 + pkin(3) * t415 + (t493 * t471 - t492 * t473) * qJD(1) + t505 * qJDD(3) + t498 * t449 + t497 * t450;
t500 = -mrSges(4,2) - mrSges(5,2);
t447 = (-mrSges(4,1) * t473 + mrSges(4,2) * t471) * qJD(1);
t457 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t489;
t412 = m(4) * t422 + qJDD(3) * mrSges(4,1) + qJD(3) * t457 + (-mrSges(4,3) - mrSges(5,3)) * t449 + (-t446 - t447) * t490 + t486;
t485 = m(5) * t420 + t450 * mrSges(5,3) + t446 * t489;
t454 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t490;
t491 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t490 - t454;
t413 = m(4) * t423 + t450 * mrSges(4,3) + t491 * qJD(3) + t500 * qJDD(3) + t447 * t489 + t485;
t481 = -t471 * t412 + t473 * t413;
t403 = m(3) * t429 - t475 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t481;
t428 = t470 * t445 - t469 * t448;
t479 = -qJDD(1) * pkin(2) - t428;
t425 = -t475 * pkin(5) + t479;
t421 = t453 * t490 - t450 * pkin(3) + qJDD(4) + (-qJ(4) * t467 - pkin(5)) * t475 + t479;
t480 = -m(5) * t421 + t450 * mrSges(5,1) + t456 * t489;
t476 = -m(4) * t425 + t450 * mrSges(4,1) + t500 * t449 + t457 * t489 + t491 * t490 + t480;
t408 = m(3) * t428 + qJDD(1) * mrSges(3,1) - t475 * mrSges(3,2) + t476;
t396 = t469 * t403 + t470 * t408;
t393 = m(2) * t458 + qJDD(1) * mrSges(2,1) - t475 * mrSges(2,2) + t396;
t482 = t470 * t403 - t469 * t408;
t394 = m(2) * t459 - t475 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t482;
t495 = t474 * t393 + t472 * t394;
t406 = t473 * t412 + t471 * t413;
t494 = t505 * qJD(3) + (t498 * t471 + t497 * t473) * qJD(1);
t404 = m(3) * t468 + t406;
t483 = -t472 * t393 + t474 * t394;
t416 = t449 * mrSges(5,2) + t454 * t490 - t480;
t398 = -mrSges(4,1) * t425 + mrSges(4,3) * t423 - mrSges(5,1) * t421 + mrSges(5,3) * t420 - pkin(3) * t416 + qJ(4) * t485 + t506 * t450 + t499 * t449 + (-qJ(4) * mrSges(5,2) + t497) * qJDD(3) + (-qJ(4) * t454 + t492) * qJD(3) - t494 * t490;
t400 = mrSges(4,2) * t425 + mrSges(5,2) * t421 - mrSges(4,3) * t422 - mrSges(5,3) * t419 - qJ(4) * t415 - t493 * qJD(3) + t498 * qJDD(3) + t507 * t449 + t499 * t450 + t494 * t489;
t477 = mrSges(2,1) * t458 + mrSges(3,1) * t428 - mrSges(2,2) * t459 - mrSges(3,2) * t429 + pkin(1) * t396 + pkin(2) * t476 + pkin(5) * t481 + t473 * t398 + t471 * t400 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t389 = -mrSges(3,1) * t468 + mrSges(3,3) * t429 + t475 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t406 - t504;
t388 = mrSges(3,2) * t468 - mrSges(3,3) * t428 + Ifges(3,5) * qJDD(1) - t475 * Ifges(3,6) - pkin(5) * t406 - t471 * t398 + t473 * t400;
t387 = -mrSges(2,2) * g(3) - mrSges(2,3) * t458 + Ifges(2,5) * qJDD(1) - t475 * Ifges(2,6) - qJ(2) * t396 + t470 * t388 - t469 * t389;
t386 = mrSges(2,1) * g(3) + mrSges(2,3) * t459 + t475 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t404 + qJ(2) * t482 + t469 * t388 + t470 * t389;
t1 = [-m(1) * g(1) + t483; -m(1) * g(2) + t495; (-m(1) - m(2)) * g(3) + t404; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t495 - t472 * t386 + t474 * t387; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t483 + t474 * t386 + t472 * t387; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t477; t477; t404; t504; t416;];
tauJB = t1;
