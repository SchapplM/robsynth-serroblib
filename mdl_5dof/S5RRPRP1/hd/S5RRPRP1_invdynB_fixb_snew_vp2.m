% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:09
% EndTime: 2019-12-05 18:22:12
% DurationCPUTime: 1.76s
% Computational Cost: add. (20235->210), mult. (26070->250), div. (0->0), fcn. (12664->8), ass. (0->90)
t480 = Ifges(5,1) + Ifges(6,1);
t473 = Ifges(5,4) + Ifges(6,4);
t472 = Ifges(5,5) + Ifges(6,5);
t479 = Ifges(5,2) + Ifges(6,2);
t478 = Ifges(5,6) + Ifges(6,6);
t477 = Ifges(5,3) + Ifges(6,3);
t438 = qJD(1) + qJD(2);
t436 = t438 ^ 2;
t446 = sin(qJ(1));
t449 = cos(qJ(1));
t429 = t449 * g(2) + t446 * g(3);
t421 = qJDD(1) * pkin(1) + t429;
t428 = t446 * g(2) - t449 * g(3);
t450 = qJD(1) ^ 2;
t422 = -t450 * pkin(1) + t428;
t445 = sin(qJ(2));
t448 = cos(qJ(2));
t399 = t448 * t421 - t445 * t422;
t437 = qJDD(1) + qJDD(2);
t397 = t437 * pkin(2) + t399;
t400 = t445 * t421 + t448 * t422;
t398 = -t436 * pkin(2) + t400;
t442 = sin(pkin(8));
t443 = cos(pkin(8));
t392 = t443 * t397 - t442 * t398;
t452 = -t437 * pkin(3) - t392;
t390 = -t436 * pkin(7) + t452;
t444 = sin(qJ(4));
t447 = cos(qJ(4));
t464 = qJD(4) * t438;
t459 = t447 * t464;
t415 = t444 * t437 + t459;
t416 = t447 * t437 - t444 * t464;
t469 = t438 * t447;
t427 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t469;
t470 = t438 * t444;
t423 = qJD(4) * pkin(4) - qJ(5) * t470;
t440 = t447 ^ 2;
t386 = t423 * t470 - t416 * pkin(4) + qJDD(5) + (-qJ(5) * t440 - pkin(7)) * t436 + t452;
t426 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t469;
t454 = m(6) * t386 - t416 * mrSges(6,1) - t426 * t469;
t474 = -mrSges(5,2) - mrSges(6,2);
t476 = -m(5) * t390 + t416 * mrSges(5,1) + t474 * t415 + t427 * t469 - t454;
t475 = pkin(4) * t436;
t393 = t442 * t397 + t443 * t398;
t391 = -t436 * pkin(3) + t437 * pkin(7) + t393;
t441 = -g(1) + qJDD(3);
t388 = t447 * t391 + t444 * t441;
t414 = (-mrSges(5,1) * t447 + mrSges(5,2) * t444) * t438;
t463 = qJD(5) * t438;
t385 = t416 * qJ(5) - qJD(4) * t423 - t440 * t475 + 0.2e1 * t447 * t463 + t388;
t413 = (-mrSges(6,1) * t447 + mrSges(6,2) * t444) * t438;
t461 = m(6) * t385 + t416 * mrSges(6,3) + t413 * t469;
t424 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t470;
t465 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t470 - t424;
t380 = m(5) * t388 + t416 * mrSges(5,3) + t465 * qJD(4) + t474 * qJDD(4) + t414 * t469 + t461;
t378 = t447 * t380;
t431 = t447 * t441;
t387 = -t444 * t391 + t431;
t384 = qJDD(4) * pkin(4) + t431 + (-t415 + t459) * qJ(5) + (t447 * t475 - t391 - 0.2e1 * t463) * t444;
t462 = m(6) * t384 + qJDD(4) * mrSges(6,1) + qJD(4) * t426;
t379 = m(5) * t387 + qJDD(4) * mrSges(5,1) + qJD(4) * t427 + (-t413 - t414) * t470 + (-mrSges(5,3) - mrSges(6,3)) * t415 + t462;
t372 = m(4) * t393 - t436 * mrSges(4,1) - t437 * mrSges(4,2) - t444 * t379 + t378;
t458 = t438 * t465;
t375 = m(4) * t392 + t437 * mrSges(4,1) - t436 * mrSges(4,2) + t444 * t458 + t476;
t367 = t442 * t372 + t443 * t375;
t365 = m(3) * t399 + t437 * mrSges(3,1) - t436 * mrSges(3,2) + t367;
t455 = t443 * t372 - t442 * t375;
t366 = m(3) * t400 - t436 * mrSges(3,1) - t437 * mrSges(3,2) + t455;
t360 = t448 * t365 + t445 * t366;
t373 = t447 * t379 + t444 * t380;
t468 = (t472 * t444 + t478 * t447) * t438 + t477 * qJD(4);
t467 = (-t473 * t444 - t479 * t447) * t438 - t478 * qJD(4);
t466 = (t480 * t444 + t473 * t447) * t438 + t472 * qJD(4);
t460 = m(4) * t441 + t373;
t456 = -t445 * t365 + t448 * t366;
t358 = m(2) * t428 - t450 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t456;
t359 = m(2) * t429 + qJDD(1) * mrSges(2,1) - t450 * mrSges(2,2) + t360;
t457 = t449 * t358 - t446 * t359;
t453 = -t446 * t358 - t449 * t359;
t381 = -t415 * mrSges(6,3) - t413 * t470 + t462;
t369 = mrSges(5,2) * t390 + mrSges(6,2) * t386 - mrSges(5,3) * t387 - mrSges(6,3) * t384 - qJ(5) * t381 + t467 * qJD(4) + t472 * qJDD(4) + t480 * t415 + t473 * t416 + t468 * t469;
t368 = -mrSges(5,1) * t390 + mrSges(5,3) * t388 - mrSges(6,1) * t386 + mrSges(6,3) * t385 - pkin(4) * t454 + qJ(5) * t461 + (-pkin(4) * t424 - t468) * t470 + t479 * t416 + (-pkin(4) * mrSges(6,2) + t473) * t415 + (-qJ(5) * mrSges(6,2) + t478) * qJDD(4) + (-qJ(5) * t424 + t466) * qJD(4);
t361 = -mrSges(4,1) * t441 - mrSges(5,1) * t387 - mrSges(6,1) * t384 + mrSges(5,2) * t388 + mrSges(6,2) * t385 + mrSges(4,3) * t393 + t436 * Ifges(4,5) + Ifges(4,6) * t437 - pkin(3) * t373 - pkin(4) * t381 - t478 * t416 - t472 * t415 - t477 * qJDD(4) + (t467 * t444 + t466 * t447) * t438;
t356 = mrSges(4,2) * t441 - mrSges(4,3) * t392 + Ifges(4,5) * t437 - t436 * Ifges(4,6) - pkin(7) * t373 - t444 * t368 + t447 * t369;
t355 = -mrSges(3,2) * g(1) - mrSges(3,3) * t399 + Ifges(3,5) * t437 - t436 * Ifges(3,6) - qJ(3) * t367 + t443 * t356 - t442 * t361;
t354 = mrSges(3,1) * g(1) + mrSges(3,3) * t400 + t436 * Ifges(3,5) + Ifges(3,6) * t437 - pkin(2) * t460 + qJ(3) * t455 + t442 * t356 + t443 * t361;
t353 = -mrSges(2,2) * g(1) - mrSges(2,3) * t429 + Ifges(2,5) * qJDD(1) - t450 * Ifges(2,6) - pkin(6) * t360 - t445 * t354 + t448 * t355;
t352 = Ifges(2,6) * qJDD(1) + t450 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t428 + t445 * t355 + t448 * t354 - pkin(1) * (-m(3) * g(1) + t460) + pkin(6) * t456;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t460; -m(1) * g(2) + t453; -m(1) * g(3) + t457; pkin(1) * t360 + pkin(2) * t367 - mrSges(3,2) * t400 + mrSges(3,1) * t399 + t447 * t368 + pkin(3) * t476 + pkin(7) * t378 + mrSges(4,1) * t392 - mrSges(4,2) * t393 + mrSges(2,1) * t429 - mrSges(2,2) * t428 + Ifges(2,3) * qJDD(1) - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (pkin(3) * t458 - pkin(7) * t379 + t369) * t444 + (Ifges(4,3) + Ifges(3,3)) * t437; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t457 - t449 * t352 - t446 * t353; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t453 - t446 * t352 + t449 * t353;];
tauB = t1;
