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
% m [6x1]
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:19:44
% EndTime: 2022-01-20 10:19:46
% DurationCPUTime: 1.84s
% Computational Cost: add. (20235->210), mult. (26070->250), div. (0->0), fcn. (12664->8), ass. (0->90)
t474 = Ifges(5,1) + Ifges(6,1);
t467 = Ifges(5,4) + Ifges(6,4);
t466 = Ifges(5,5) + Ifges(6,5);
t473 = Ifges(5,2) + Ifges(6,2);
t472 = Ifges(5,6) + Ifges(6,6);
t471 = Ifges(5,3) + Ifges(6,3);
t432 = qJD(1) + qJD(2);
t430 = t432 ^ 2;
t440 = sin(qJ(1));
t443 = cos(qJ(1));
t424 = t440 * g(1) - t443 * g(2);
t417 = qJDD(1) * pkin(1) + t424;
t425 = -t443 * g(1) - t440 * g(2);
t444 = qJD(1) ^ 2;
t418 = -t444 * pkin(1) + t425;
t439 = sin(qJ(2));
t442 = cos(qJ(2));
t395 = t442 * t417 - t439 * t418;
t431 = qJDD(1) + qJDD(2);
t393 = t431 * pkin(2) + t395;
t396 = t439 * t417 + t442 * t418;
t394 = -t430 * pkin(2) + t396;
t436 = sin(pkin(8));
t437 = cos(pkin(8));
t388 = t437 * t393 - t436 * t394;
t446 = -t431 * pkin(3) - t388;
t386 = -t430 * pkin(7) + t446;
t438 = sin(qJ(4));
t441 = cos(qJ(4));
t457 = qJD(4) * t432;
t452 = t441 * t457;
t411 = t438 * t431 + t452;
t412 = t441 * t431 - t438 * t457;
t463 = t432 * t441;
t423 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t463;
t464 = t432 * t438;
t419 = qJD(4) * pkin(4) - qJ(5) * t464;
t434 = t441 ^ 2;
t382 = t419 * t464 - t412 * pkin(4) + qJDD(5) + (-qJ(5) * t434 - pkin(7)) * t430 + t446;
t422 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t463;
t447 = m(6) * t382 - t412 * mrSges(6,1) - t422 * t463;
t468 = -mrSges(5,2) - mrSges(6,2);
t470 = -m(5) * t386 + t412 * mrSges(5,1) + t468 * t411 + t423 * t463 - t447;
t469 = pkin(4) * t430;
t389 = t436 * t393 + t437 * t394;
t387 = -t430 * pkin(3) + t431 * pkin(7) + t389;
t435 = -g(3) + qJDD(3);
t384 = t441 * t387 + t438 * t435;
t410 = (-mrSges(5,1) * t441 + mrSges(5,2) * t438) * t432;
t456 = qJD(5) * t432;
t381 = t412 * qJ(5) - qJD(4) * t419 - t434 * t469 + 0.2e1 * t441 * t456 + t384;
t409 = (-mrSges(6,1) * t441 + mrSges(6,2) * t438) * t432;
t454 = m(6) * t381 + t412 * mrSges(6,3) + t409 * t463;
t420 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t464;
t458 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t464 - t420;
t376 = m(5) * t384 + t412 * mrSges(5,3) + t458 * qJD(4) + t468 * qJDD(4) + t410 * t463 + t454;
t374 = t441 * t376;
t427 = t441 * t435;
t383 = -t438 * t387 + t427;
t380 = qJDD(4) * pkin(4) + t427 + (-t411 + t452) * qJ(5) + (t441 * t469 - t387 - 0.2e1 * t456) * t438;
t455 = m(6) * t380 + qJDD(4) * mrSges(6,1) + qJD(4) * t422;
t375 = m(5) * t383 + qJDD(4) * mrSges(5,1) + qJD(4) * t423 + (-t409 - t410) * t464 + (-mrSges(5,3) - mrSges(6,3)) * t411 + t455;
t368 = m(4) * t389 - t430 * mrSges(4,1) - t431 * mrSges(4,2) - t438 * t375 + t374;
t451 = t432 * t458;
t371 = m(4) * t388 + t431 * mrSges(4,1) - t430 * mrSges(4,2) + t438 * t451 + t470;
t363 = t436 * t368 + t437 * t371;
t361 = m(3) * t395 + t431 * mrSges(3,1) - t430 * mrSges(3,2) + t363;
t448 = t437 * t368 - t436 * t371;
t362 = m(3) * t396 - t430 * mrSges(3,1) - t431 * mrSges(3,2) + t448;
t356 = t442 * t361 + t439 * t362;
t354 = m(2) * t424 + qJDD(1) * mrSges(2,1) - t444 * mrSges(2,2) + t356;
t449 = -t439 * t361 + t442 * t362;
t355 = m(2) * t425 - t444 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t449;
t462 = t443 * t354 + t440 * t355;
t369 = t441 * t375 + t438 * t376;
t461 = (t466 * t438 + t472 * t441) * t432 + t471 * qJD(4);
t460 = (-t467 * t438 - t473 * t441) * t432 - t472 * qJD(4);
t459 = (t474 * t438 + t467 * t441) * t432 + t466 * qJD(4);
t453 = m(4) * t435 + t369;
t450 = -t440 * t354 + t443 * t355;
t377 = -t411 * mrSges(6,3) - t409 * t464 + t455;
t365 = mrSges(5,2) * t386 + mrSges(6,2) * t382 - mrSges(5,3) * t383 - mrSges(6,3) * t380 - qJ(5) * t377 + t460 * qJD(4) + t466 * qJDD(4) + t474 * t411 + t467 * t412 + t461 * t463;
t364 = -mrSges(5,1) * t386 + mrSges(5,3) * t384 - mrSges(6,1) * t382 + mrSges(6,3) * t381 - pkin(4) * t447 + qJ(5) * t454 + (-pkin(4) * t420 - t461) * t464 + t473 * t412 + (-pkin(4) * mrSges(6,2) + t467) * t411 + (-qJ(5) * mrSges(6,2) + t472) * qJDD(4) + (-qJ(5) * t420 + t459) * qJD(4);
t357 = -mrSges(4,1) * t435 - mrSges(5,1) * t383 - mrSges(6,1) * t380 + mrSges(5,2) * t384 + mrSges(6,2) * t381 + mrSges(4,3) * t389 + t430 * Ifges(4,5) + Ifges(4,6) * t431 - pkin(3) * t369 - pkin(4) * t377 - t472 * t412 - t466 * t411 - t471 * qJDD(4) + (t460 * t438 + t459 * t441) * t432;
t350 = mrSges(4,2) * t435 - mrSges(4,3) * t388 + Ifges(4,5) * t431 - t430 * Ifges(4,6) - pkin(7) * t369 - t438 * t364 + t441 * t365;
t349 = -mrSges(3,2) * g(3) - mrSges(3,3) * t395 + Ifges(3,5) * t431 - t430 * Ifges(3,6) - qJ(3) * t363 + t437 * t350 - t436 * t357;
t348 = mrSges(3,1) * g(3) + mrSges(3,3) * t396 + t430 * Ifges(3,5) + Ifges(3,6) * t431 - pkin(2) * t453 + qJ(3) * t448 + t436 * t350 + t437 * t357;
t347 = -mrSges(2,2) * g(3) - mrSges(2,3) * t424 + Ifges(2,5) * qJDD(1) - t444 * Ifges(2,6) - pkin(6) * t356 - t439 * t348 + t442 * t349;
t346 = Ifges(2,6) * qJDD(1) + t444 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t425 + t439 * t349 + t442 * t348 - pkin(1) * (-m(3) * g(3) + t453) + pkin(6) * t449;
t1 = [-m(1) * g(1) + t450; -m(1) * g(2) + t462; (-m(1) - m(2) - m(3)) * g(3) + t453; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t462 - t440 * t346 + t443 * t347; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t450 + t443 * t346 + t440 * t347; pkin(1) * t356 + mrSges(2,1) * t424 - mrSges(2,2) * t425 + pkin(2) * t363 + mrSges(3,1) * t395 - mrSges(3,2) * t396 + t441 * t364 + pkin(3) * t470 + pkin(7) * t374 + mrSges(4,1) * t388 - mrSges(4,2) * t389 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + (pkin(3) * t451 - pkin(7) * t375 + t365) * t438 + (Ifges(4,3) + Ifges(3,3)) * t431;];
tauB = t1;
