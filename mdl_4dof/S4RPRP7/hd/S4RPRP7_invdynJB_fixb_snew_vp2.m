% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPRP7
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPRP7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP7_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_invdynJB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP7_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP7_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:06
% EndTime: 2019-12-31 16:47:07
% DurationCPUTime: 0.75s
% Computational Cost: add. (2641->167), mult. (4843->197), div. (0->0), fcn. (1847->4), ass. (0->70)
t476 = Ifges(4,1) + Ifges(5,1);
t464 = Ifges(4,4) - Ifges(5,5);
t474 = Ifges(5,4) + Ifges(4,5);
t475 = Ifges(4,2) + Ifges(5,3);
t473 = Ifges(4,6) - Ifges(5,6);
t472 = Ifges(4,3) + Ifges(5,2);
t438 = sin(qJ(3));
t440 = cos(qJ(3));
t471 = -t473 * qJD(3) + (t475 * t438 - t464 * t440) * qJD(1);
t470 = t474 * qJD(3) + (-t464 * t438 + t476 * t440) * qJD(1);
t439 = sin(qJ(1));
t441 = cos(qJ(1));
t425 = -t441 * g(1) - t439 * g(2);
t469 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t425;
t468 = -pkin(1) - pkin(5);
t467 = t438 * g(3);
t466 = mrSges(2,1) - mrSges(3,2);
t465 = -mrSges(4,3) - mrSges(5,2);
t463 = Ifges(2,5) - Ifges(3,4);
t462 = Ifges(2,6) - Ifges(3,5);
t424 = t439 * g(1) - t441 * g(2);
t443 = qJD(1) ^ 2;
t449 = -t443 * qJ(2) + qJDD(2) - t424;
t390 = t468 * qJDD(1) + t449;
t386 = -t440 * g(3) + t438 * t390;
t458 = qJD(1) * qJD(3);
t415 = t438 * qJDD(1) + t440 * t458;
t459 = qJD(1) * t440;
t421 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t459;
t413 = (mrSges(5,1) * t438 - mrSges(5,3) * t440) * qJD(1);
t452 = qJD(1) * (-t413 - (mrSges(4,1) * t438 + mrSges(4,2) * t440) * qJD(1));
t412 = (pkin(3) * t438 - qJ(4) * t440) * qJD(1);
t442 = qJD(3) ^ 2;
t460 = qJD(1) * t438;
t382 = -t442 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) - t412 * t460 + t386;
t422 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t459;
t456 = m(5) * t382 + qJDD(3) * mrSges(5,3) + qJD(3) * t422;
t372 = m(4) * t386 - qJDD(3) * mrSges(4,2) - qJD(3) * t421 + t465 * t415 + t438 * t452 + t456;
t385 = t440 * t390 + t467;
t416 = t440 * qJDD(1) - t438 * t458;
t420 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t460;
t383 = -qJDD(3) * pkin(3) - t467 - t442 * qJ(4) + qJDD(4) + (qJD(1) * t412 - t390) * t440;
t423 = -mrSges(5,2) * t460 + qJD(3) * mrSges(5,3);
t451 = -m(5) * t383 + qJDD(3) * mrSges(5,1) + qJD(3) * t423;
t373 = m(4) * t385 + qJDD(3) * mrSges(4,1) + qJD(3) * t420 + t465 * t416 + t440 * t452 + t451;
t366 = t438 * t372 + t440 * t373;
t396 = -qJDD(1) * pkin(1) + t449;
t448 = -m(3) * t396 + t443 * mrSges(3,3) - t366;
t362 = m(2) * t424 - t443 * mrSges(2,2) + t466 * qJDD(1) + t448;
t394 = t443 * pkin(1) + t469;
t389 = t468 * t443 - t469;
t379 = t415 * pkin(3) - t416 * qJ(4) + (-0.2e1 * qJD(4) * t440 + (pkin(3) * t440 + qJ(4) * t438) * qJD(3)) * qJD(1) + t389;
t374 = m(5) * t379 + t415 * mrSges(5,1) - t416 * mrSges(5,3) - t422 * t459 + t423 * t460;
t447 = -m(4) * t389 - t415 * mrSges(4,1) - t416 * mrSges(4,2) - t420 * t460 - t421 * t459 - t374;
t445 = -m(3) * t394 + t443 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t447;
t369 = m(2) * t425 - t443 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t445;
t461 = t441 * t362 + t439 * t369;
t455 = -t439 * t362 + t441 * t369;
t454 = t440 * t372 - t438 * t373;
t453 = qJD(1) * (-t472 * qJD(3) + (t473 * t438 - t474 * t440) * qJD(1));
t359 = -mrSges(4,1) * t389 - mrSges(5,1) * t379 + mrSges(5,2) * t382 + mrSges(4,3) * t386 - pkin(3) * t374 + t470 * qJD(3) + t473 * qJDD(3) - t475 * t415 + t464 * t416 + t440 * t453;
t360 = mrSges(4,2) * t389 + mrSges(5,2) * t383 - mrSges(4,3) * t385 - mrSges(5,3) * t379 - qJ(4) * t374 + t471 * qJD(3) + t474 * qJDD(3) - t464 * t415 + t476 * t416 + t438 * t453;
t364 = qJDD(1) * mrSges(3,2) - t448;
t446 = mrSges(2,1) * t424 - mrSges(2,2) * t425 + mrSges(3,2) * t396 - mrSges(3,3) * t394 - pkin(1) * t364 - pkin(5) * t366 + qJ(2) * t445 - t438 * t359 + t440 * t360 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t377 = t416 * mrSges(5,2) + t413 * t459 - t451;
t444 = mrSges(4,1) * t385 - mrSges(5,1) * t383 - mrSges(4,2) * t386 + mrSges(5,3) * t382 - pkin(3) * t377 + qJ(4) * t456 + (-qJ(4) * t413 + t470) * t460 - t471 * t459 + t474 * t416 + (-qJ(4) * mrSges(5,2) - t473) * t415 + t472 * qJDD(3);
t365 = -m(3) * g(3) + t454;
t357 = t444 - t462 * t443 - mrSges(2,3) * t424 + mrSges(3,1) * t396 + pkin(2) * t366 - qJ(2) * t365 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t463 * qJDD(1);
t356 = -mrSges(3,1) * t394 + mrSges(2,3) * t425 - pkin(1) * t365 - pkin(2) * t447 - pkin(5) * t454 + t466 * g(3) + t462 * qJDD(1) - t440 * t359 - t438 * t360 + t463 * t443;
t1 = [-m(1) * g(1) + t455; -m(1) * g(2) + t461; (-m(1) - m(2) - m(3)) * g(3) + t454; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t461 - t439 * t356 + t441 * t357; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t455 + t441 * t356 + t439 * t357; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t446; t446; t364; t444; t377;];
tauJB = t1;
