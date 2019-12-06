% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:33
% EndTime: 2019-12-05 18:01:35
% DurationCPUTime: 1.70s
% Computational Cost: add. (19107->209), mult. (26070->250), div. (0->0), fcn. (12664->8), ass. (0->91)
t481 = Ifges(5,1) + Ifges(6,1);
t474 = Ifges(5,4) + Ifges(6,4);
t473 = Ifges(5,5) + Ifges(6,5);
t480 = Ifges(5,2) + Ifges(6,2);
t479 = Ifges(5,6) + Ifges(6,6);
t478 = Ifges(5,3) + Ifges(6,3);
t438 = qJD(1) + qJD(3);
t436 = t438 ^ 2;
t446 = sin(qJ(1));
t449 = cos(qJ(1));
t428 = t449 * g(2) + t446 * g(3);
t420 = qJDD(1) * pkin(1) + t428;
t427 = t446 * g(2) - t449 * g(3);
t450 = qJD(1) ^ 2;
t421 = -t450 * pkin(1) + t427;
t442 = sin(pkin(8));
t443 = cos(pkin(8));
t398 = t443 * t420 - t442 * t421;
t396 = qJDD(1) * pkin(2) + t398;
t399 = t442 * t420 + t443 * t421;
t397 = -t450 * pkin(2) + t399;
t445 = sin(qJ(3));
t448 = cos(qJ(3));
t391 = t448 * t396 - t445 * t397;
t437 = qJDD(1) + qJDD(3);
t452 = -t437 * pkin(3) - t391;
t389 = -t436 * pkin(7) + t452;
t444 = sin(qJ(4));
t447 = cos(qJ(4));
t465 = qJD(4) * t438;
t460 = t447 * t465;
t414 = t444 * t437 + t460;
t415 = t447 * t437 - t444 * t465;
t470 = t438 * t447;
t426 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t470;
t471 = t438 * t444;
t422 = qJD(4) * pkin(4) - qJ(5) * t471;
t440 = t447 ^ 2;
t385 = t422 * t471 - t415 * pkin(4) + qJDD(5) + (-qJ(5) * t440 - pkin(7)) * t436 + t452;
t425 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t470;
t454 = m(6) * t385 - t415 * mrSges(6,1) - t425 * t470;
t475 = -mrSges(5,2) - mrSges(6,2);
t477 = -m(5) * t389 + t415 * mrSges(5,1) + t475 * t414 + t426 * t470 - t454;
t476 = pkin(4) * t436;
t392 = t445 * t396 + t448 * t397;
t390 = -t436 * pkin(3) + t437 * pkin(7) + t392;
t441 = -g(1) + qJDD(2);
t387 = t447 * t390 + t444 * t441;
t413 = (-mrSges(5,1) * t447 + mrSges(5,2) * t444) * t438;
t464 = qJD(5) * t438;
t384 = t415 * qJ(5) - qJD(4) * t422 - t440 * t476 + 0.2e1 * t447 * t464 + t387;
t412 = (-mrSges(6,1) * t447 + mrSges(6,2) * t444) * t438;
t462 = m(6) * t384 + t415 * mrSges(6,3) + t412 * t470;
t423 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t471;
t466 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t471 - t423;
t379 = m(5) * t387 + t415 * mrSges(5,3) + t466 * qJD(4) + t475 * qJDD(4) + t413 * t470 + t462;
t377 = t447 * t379;
t430 = t447 * t441;
t386 = -t444 * t390 + t430;
t383 = qJDD(4) * pkin(4) + t430 + (-t414 + t460) * qJ(5) + (t447 * t476 - t390 - 0.2e1 * t464) * t444;
t463 = m(6) * t383 + qJDD(4) * mrSges(6,1) + qJD(4) * t425;
t378 = m(5) * t386 + qJDD(4) * mrSges(5,1) + qJD(4) * t426 + (-t412 - t413) * t471 + (-mrSges(5,3) - mrSges(6,3)) * t414 + t463;
t371 = m(4) * t392 - t436 * mrSges(4,1) - t437 * mrSges(4,2) - t444 * t378 + t377;
t459 = t438 * t466;
t374 = m(4) * t391 + t437 * mrSges(4,1) - t436 * mrSges(4,2) + t444 * t459 + t477;
t366 = t445 * t371 + t448 * t374;
t364 = m(3) * t398 + qJDD(1) * mrSges(3,1) - t450 * mrSges(3,2) + t366;
t456 = t448 * t371 - t445 * t374;
t365 = m(3) * t399 - t450 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t456;
t359 = t443 * t364 + t442 * t365;
t372 = t447 * t378 + t444 * t379;
t469 = (t473 * t444 + t479 * t447) * t438 + t478 * qJD(4);
t468 = (-t474 * t444 - t480 * t447) * t438 - t479 * qJD(4);
t467 = (t481 * t444 + t474 * t447) * t438 + t473 * qJD(4);
t461 = m(4) * t441 + t372;
t457 = -t442 * t364 + t443 * t365;
t357 = m(2) * t427 - t450 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t457;
t358 = m(2) * t428 + qJDD(1) * mrSges(2,1) - t450 * mrSges(2,2) + t359;
t458 = t449 * t357 - t446 * t358;
t455 = m(3) * t441 + t461;
t453 = -t446 * t357 - t449 * t358;
t380 = -t414 * mrSges(6,3) - t412 * t471 + t463;
t368 = mrSges(5,2) * t389 + mrSges(6,2) * t385 - mrSges(5,3) * t386 - mrSges(6,3) * t383 - qJ(5) * t380 + t468 * qJD(4) + t473 * qJDD(4) + t481 * t414 + t474 * t415 + t469 * t470;
t367 = -mrSges(5,1) * t389 + mrSges(5,3) * t387 - mrSges(6,1) * t385 + mrSges(6,3) * t384 - pkin(4) * t454 + qJ(5) * t462 + (-pkin(4) * t423 - t469) * t471 + t480 * t415 + (-pkin(4) * mrSges(6,2) + t474) * t414 + (-qJ(5) * mrSges(6,2) + t479) * qJDD(4) + (-qJ(5) * t423 + t467) * qJD(4);
t360 = -mrSges(4,1) * t441 - mrSges(5,1) * t386 - mrSges(6,1) * t383 + mrSges(5,2) * t387 + mrSges(6,2) * t384 + mrSges(4,3) * t392 + t436 * Ifges(4,5) + Ifges(4,6) * t437 - pkin(3) * t372 - pkin(4) * t380 - t479 * t415 - t473 * t414 - t478 * qJDD(4) + (t468 * t444 + t467 * t447) * t438;
t355 = mrSges(4,2) * t441 - mrSges(4,3) * t391 + Ifges(4,5) * t437 - t436 * Ifges(4,6) - pkin(7) * t372 - t444 * t367 + t447 * t368;
t354 = mrSges(3,2) * t441 - mrSges(3,3) * t398 + Ifges(3,5) * qJDD(1) - t450 * Ifges(3,6) - pkin(6) * t366 + t448 * t355 - t445 * t360;
t353 = -mrSges(3,1) * t441 + mrSges(3,3) * t399 + t450 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t461 + pkin(6) * t456 + t445 * t355 + t448 * t360;
t352 = -mrSges(2,2) * g(1) - mrSges(2,3) * t428 + Ifges(2,5) * qJDD(1) - t450 * Ifges(2,6) - qJ(2) * t359 - t442 * t353 + t443 * t354;
t351 = mrSges(2,1) * g(1) + mrSges(2,3) * t427 + t450 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t455 + qJ(2) * t457 + t443 * t353 + t442 * t354;
t1 = [(-m(1) - m(2)) * g(1) + t455; -m(1) * g(2) + t453; -m(1) * g(3) + t458; pkin(1) * t359 + pkin(2) * t366 + mrSges(3,1) * t398 - mrSges(3,2) * t399 + pkin(3) * t477 + pkin(7) * t377 + t447 * t367 + mrSges(4,1) * t391 - mrSges(4,2) * t392 + mrSges(2,1) * t428 - mrSges(2,2) * t427 + Ifges(4,3) * t437 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (pkin(3) * t459 - pkin(7) * t378 + t368) * t444 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t458 - t449 * t351 - t446 * t352; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t453 - t446 * t351 + t449 * t352;];
tauB = t1;
