% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-05-05 14:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:26:36
% EndTime: 2019-05-05 14:26:37
% DurationCPUTime: 1.18s
% Computational Cost: add. (2970->214), mult. (5661->248), div. (0->0), fcn. (2477->6), ass. (0->91)
t419 = sin(qJ(4));
t422 = cos(qJ(4));
t457 = Ifges(5,4) + Ifges(6,6);
t471 = t422 * (Ifges(5,1) + Ifges(6,2)) - t419 * t457;
t470 = t422 * t457 + t419 * (-Ifges(5,2) - Ifges(6,3));
t469 = 2 * qJD(1);
t456 = Ifges(5,5) - Ifges(6,4);
t455 = Ifges(5,6) - Ifges(6,5);
t466 = (t470 * qJD(1) + t455 * qJD(4)) * t422;
t420 = sin(qJ(1));
t423 = cos(qJ(1));
t440 = -g(1) * t423 - g(2) * t420;
t465 = qJDD(1) * qJ(2) + (qJD(2) * t469) + t440;
t425 = qJD(1) ^ 2;
t449 = t420 * g(1) - t423 * g(2);
t377 = -qJDD(1) * pkin(1) - t425 * qJ(2) + qJDD(2) - t449;
t435 = qJDD(1) * qJ(3) + (qJD(3) * t469) - t377;
t461 = pkin(7) * t425;
t369 = t435 - t461;
t446 = qJD(1) * qJD(4);
t443 = t422 * t446;
t394 = qJDD(1) * t419 + t443;
t406 = t419 * t446;
t395 = qJDD(1) * t422 - t406;
t447 = qJD(1) * t422;
t398 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t447;
t400 = mrSges(6,1) * t447 + qJD(4) * mrSges(6,2);
t463 = -2 * qJD(5);
t426 = pkin(4) * t443 + t447 * t463 + t435 + (-t395 + t406) * qJ(5);
t352 = pkin(4) * t394 + t426 - t461;
t401 = pkin(5) * t447 - qJD(4) * pkin(8);
t416 = t419 ^ 2;
t462 = pkin(4) + pkin(8);
t347 = t426 + t462 * t394 + (-pkin(5) * t416 - pkin(7)) * t425 - t401 * t447;
t391 = (pkin(4) * t419 - qJ(5) * t422) * qJD(1);
t424 = qJD(4) ^ 2;
t373 = qJDD(3) + (-pkin(1) - qJ(3)) * t425 + t465;
t370 = -qJDD(1) * pkin(7) + t373;
t454 = t370 * t422;
t434 = -qJ(5) * t424 + t391 * t447 + qJDD(5) - t454;
t460 = pkin(8) * t425;
t350 = pkin(5) * t395 - t462 * qJDD(4) + (pkin(5) * t446 + t422 * t460 - g(3)) * t419 + t434;
t418 = sin(qJ(6));
t421 = cos(qJ(6));
t345 = -t347 * t418 + t350 * t421;
t448 = qJD(1) * t419;
t389 = -qJD(4) * t418 + t421 * t448;
t364 = qJD(6) * t389 + qJDD(4) * t421 + t394 * t418;
t390 = qJD(4) * t421 + t418 * t448;
t368 = -mrSges(7,1) * t389 + mrSges(7,2) * t390;
t404 = qJD(6) + t447;
t374 = -mrSges(7,2) * t404 + mrSges(7,3) * t389;
t388 = qJDD(6) + t395;
t342 = m(7) * t345 + mrSges(7,1) * t388 - mrSges(7,3) * t364 - t368 * t390 + t374 * t404;
t346 = t347 * t421 + t350 * t418;
t363 = -qJD(6) * t390 - qJDD(4) * t418 + t394 * t421;
t375 = mrSges(7,1) * t404 - mrSges(7,3) * t390;
t343 = m(7) * t346 - mrSges(7,2) * t388 + mrSges(7,3) * t363 + t368 * t389 - t375 * t404;
t442 = -t342 * t418 + t421 * t343;
t438 = m(6) * t352 - t395 * mrSges(6,3) + t442;
t399 = mrSges(6,1) * t448 - qJD(4) * mrSges(6,3);
t450 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t448 - t399;
t458 = mrSges(5,1) - mrSges(6,2);
t464 = -m(4) * t435 - m(5) * t369 - mrSges(5,2) * t395 - t458 * t394 + ((-t398 + t400) * t422 - t450 * t419) * qJD(1) - t438;
t459 = g(3) * t419;
t365 = t454 + t459;
t366 = -g(3) * t422 + t419 * t370;
t428 = -pkin(4) * t424 + qJDD(4) * qJ(5) - t391 * t448 + t366;
t353 = qJD(4) * t463 - t428;
t349 = -t416 * t460 - pkin(5) * t394 + ((2 * qJD(5)) + t401) * qJD(4) + t428;
t433 = -m(7) * t349 + mrSges(7,1) * t363 - t364 * mrSges(7,2) + t374 * t389 - t390 * t375;
t427 = -m(6) * t353 + qJDD(4) * mrSges(6,3) + qJD(4) * t400 - t433;
t337 = t342 * t421 + t343 * t418;
t354 = -qJDD(4) * pkin(4) + t434 - t459;
t432 = -m(6) * t354 - t395 * mrSges(6,1) - t337;
t392 = (-mrSges(6,2) * t419 - mrSges(6,3) * t422) * qJD(1);
t441 = qJD(1) * (-t392 - (mrSges(5,1) * t419 + mrSges(5,2) * t422) * qJD(1));
t453 = t419 * (m(5) * t366 - qJDD(4) * mrSges(5,2) - qJD(4) * t398 + (-mrSges(5,3) - mrSges(6,1)) * t394 + t419 * t441 + t427) + t422 * (m(5) * t365 - mrSges(5,3) * t395 + t450 * qJD(4) + t458 * qJDD(4) + t422 * t441 + t432);
t451 = t471 * qJD(1) + t456 * qJD(4);
t436 = m(4) * t373 + qJDD(1) * mrSges(4,2) - (mrSges(4,3) * t425) + t453;
t357 = Ifges(7,4) * t390 + Ifges(7,2) * t389 + Ifges(7,6) * t404;
t358 = Ifges(7,1) * t390 + Ifges(7,4) * t389 + Ifges(7,5) * t404;
t430 = mrSges(7,1) * t345 - mrSges(7,2) * t346 + Ifges(7,5) * t364 + Ifges(7,6) * t363 + Ifges(7,3) * t388 + t390 * t357 - t389 * t358;
t376 = pkin(1) * t425 - t465;
t356 = Ifges(7,5) * t390 + Ifges(7,6) * t389 + Ifges(7,3) * t404;
t339 = mrSges(7,2) * t349 - mrSges(7,3) * t345 + Ifges(7,1) * t364 + Ifges(7,4) * t363 + Ifges(7,5) * t388 + t356 * t389 - t357 * t404;
t338 = -mrSges(7,1) * t349 + mrSges(7,3) * t346 + Ifges(7,4) * t364 + Ifges(7,2) * t363 + Ifges(7,6) * t388 - t356 * t390 + t358 * t404;
t336 = qJDD(4) * mrSges(6,2) + qJD(4) * t399 + t392 * t447 - t432;
t335 = -mrSges(6,2) * t394 + (-t399 * t419 - t400 * t422) * qJD(1) + t438;
t333 = m(3) * t377 + (-mrSges(4,2) - mrSges(3,3)) * t425 + (mrSges(3,2) - mrSges(4,3)) * qJDD(1) + t464;
t1 = [qJ(2) * (-m(3) * t376 + (mrSges(3,2) * t425) + t436) - pkin(1) * t333 + mrSges(2,1) * t449 - mrSges(2,2) * t440 + mrSges(3,2) * t377 - mrSges(3,3) * t376 - pkin(7) * t453 + mrSges(4,3) * t435 + t422 * (mrSges(6,1) * t354 + mrSges(5,2) * t369 - mrSges(5,3) * t365 - mrSges(6,3) * t352 + pkin(5) * t337 - qJ(5) * t335 + t430) - t419 * (-mrSges(5,1) * t369 - mrSges(6,1) * t353 + mrSges(6,2) * t352 + mrSges(5,3) * t366 - pkin(4) * t335 - pkin(5) * t433 - pkin(8) * t442 - t421 * t338 - t418 * t339) + mrSges(4,2) * t373 + t471 * t395 + (-t419 * t455 + t422 * t456) * qJDD(4) + (-t419 * t451 - t466) * qJD(4) - t470 * t394 + (mrSges(3,3) * qJ(2) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) + (mrSges(4,2) * t425 + qJDD(1) * mrSges(4,3) - t464) * qJ(3); t333; t436; mrSges(5,1) * t365 - mrSges(5,2) * t366 + mrSges(6,2) * t354 - mrSges(6,3) * t353 + t421 * t339 - t418 * t338 - pkin(8) * t337 - pkin(4) * t336 + qJ(5) * t427 + t456 * t395 + (-qJ(5) * mrSges(6,1) - t455) * t394 + (Ifges(5,3) + Ifges(6,1)) * qJDD(4) + (t466 + (-qJ(5) * t392 + t451) * t419) * qJD(1); t336; t430;];
tauJ  = t1;
