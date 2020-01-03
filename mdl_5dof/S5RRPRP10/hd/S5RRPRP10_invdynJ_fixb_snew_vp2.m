% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRP10
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRP10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP10_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP10_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:35
% EndTime: 2019-12-31 20:09:37
% DurationCPUTime: 1.73s
% Computational Cost: add. (4259->245), mult. (8686->277), div. (0->0), fcn. (4557->6), ass. (0->99)
t460 = Ifges(5,4) + Ifges(6,4);
t476 = Ifges(5,2) + Ifges(6,2);
t470 = Ifges(5,6) + Ifges(6,6);
t475 = -2 * qJD(3);
t474 = Ifges(3,1) + Ifges(4,2);
t473 = Ifges(5,1) + Ifges(6,1);
t461 = Ifges(3,4) + Ifges(4,6);
t459 = Ifges(3,5) - Ifges(4,4);
t472 = Ifges(5,5) + Ifges(6,5);
t471 = Ifges(3,2) + Ifges(4,3);
t458 = Ifges(3,6) - Ifges(4,5);
t469 = Ifges(3,3) + Ifges(4,1);
t468 = Ifges(5,3) + Ifges(6,3);
t421 = sin(qJ(4));
t424 = cos(qJ(4));
t425 = cos(qJ(2));
t447 = qJD(1) * t425;
t397 = -qJD(2) * t421 - t424 * t447;
t398 = qJD(2) * t424 - t421 * t447;
t422 = sin(qJ(2));
t448 = qJD(1) * t422;
t412 = qJD(4) + t448;
t467 = t476 * t397 + t460 * t398 + t470 * t412;
t428 = qJD(1) ^ 2;
t423 = sin(qJ(1));
t426 = cos(qJ(1));
t437 = -g(1) * t426 - g(2) * t423;
t390 = -pkin(1) * t428 + qJDD(1) * pkin(6) + t437;
t377 = -g(3) * t422 + t425 * t390;
t399 = (-t425 * pkin(2) - t422 * qJ(3)) * qJD(1);
t427 = qJD(2) ^ 2;
t348 = pkin(2) * t427 - qJDD(2) * qJ(3) + qJD(2) * t475 - t399 * t447 - t377;
t446 = qJD(1) * qJD(2);
t440 = t422 * t446;
t403 = qJDD(1) * t425 - t440;
t366 = -qJD(4) * t398 - qJDD(2) * t421 - t403 * t424;
t371 = -mrSges(6,2) * t412 + mrSges(6,3) * t397;
t466 = -t366 * mrSges(6,1) - t397 * t371;
t372 = -mrSges(5,2) * t412 + mrSges(5,3) * t397;
t465 = -(mrSges(5,1) + mrSges(6,1)) * t366 - (t371 + t372) * t397;
t464 = pkin(6) * t428;
t463 = mrSges(3,1) - mrSges(4,2);
t410 = pkin(3) * t448 - qJD(2) * pkin(7);
t420 = t425 ^ 2;
t441 = t425 * t446;
t402 = qJDD(1) * t422 + t441;
t439 = g(1) * t423 - t426 * g(2);
t435 = -qJDD(1) * pkin(1) - t439;
t430 = pkin(2) * t440 + t448 * t475 + (-t402 - t441) * qJ(3) + t435;
t339 = -t410 * t448 + (-pkin(3) * t420 - pkin(6)) * t428 + (-pkin(2) - pkin(7)) * t403 + t430;
t376 = -t425 * g(3) - t422 * t390;
t349 = -qJDD(2) * pkin(2) - qJ(3) * t427 + t399 * t448 + qJDD(3) - t376;
t344 = (-t422 * t425 * t428 - qJDD(2)) * pkin(7) + (t402 - t441) * pkin(3) + t349;
t335 = t424 * t339 + t421 * t344;
t455 = -t470 * t397 - t472 * t398 - t468 * t412;
t454 = -t460 * t397 - t473 * t398 - t472 * t412;
t452 = t469 * qJD(2) + (t422 * t459 + t425 * t458) * qJD(1);
t451 = t458 * qJD(2) + (t422 * t461 + t425 * t471) * qJD(1);
t450 = t459 * qJD(2) + (t422 * t474 + t425 * t461) * qJD(1);
t408 = -mrSges(4,1) * t447 - qJD(2) * mrSges(4,3);
t449 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t447 - t408;
t334 = -t339 * t421 + t424 * t344;
t367 = qJD(4) * t397 + qJDD(2) * t424 - t403 * t421;
t396 = qJDD(4) + t402;
t330 = -0.2e1 * qJD(5) * t398 + (t397 * t412 - t367) * qJ(5) + (t397 * t398 + t396) * pkin(4) + t334;
t444 = m(6) * t330 + t396 * mrSges(6,1) + t412 * t371;
t373 = pkin(4) * t412 - qJ(5) * t398;
t395 = t397 ^ 2;
t332 = -pkin(4) * t395 + qJ(5) * t366 + 0.2e1 * qJD(5) * t397 - t373 * t412 + t335;
t369 = -mrSges(6,1) * t397 + mrSges(6,2) * t398;
t443 = m(6) * t332 + t366 * mrSges(6,3) + t397 * t369;
t343 = -pkin(7) * t420 * t428 + pkin(3) * t403 + qJD(2) * t410 - t348;
t337 = -pkin(4) * t366 - qJ(5) * t395 + t373 * t398 + qJDD(5) + t343;
t374 = mrSges(6,1) * t412 - mrSges(6,3) * t398;
t442 = m(6) * t337 + t367 * mrSges(6,2) + t398 * t374;
t370 = -mrSges(5,1) * t397 + mrSges(5,2) * t398;
t322 = m(5) * t334 + mrSges(5,1) * t396 + t372 * t412 + (-t369 - t370) * t398 + (-mrSges(5,3) - mrSges(6,3)) * t367 + t444;
t375 = mrSges(5,1) * t412 - mrSges(5,3) * t398;
t324 = m(5) * t335 + mrSges(5,3) * t366 + t370 * t397 + (-t374 - t375) * t412 + (-mrSges(5,2) - mrSges(6,2)) * t396 + t443;
t438 = -t322 * t421 + t424 * t324;
t321 = t424 * t322 + t421 * t324;
t345 = -pkin(2) * t403 + t430 - t464;
t436 = m(4) * t345 + t438;
t434 = -m(5) * t343 - t367 * mrSges(5,2) - t398 * t375 - t442;
t432 = m(4) * t349 + t402 * mrSges(4,1) + t321;
t400 = (t425 * mrSges(4,2) - t422 * mrSges(4,3)) * qJD(1);
t409 = mrSges(4,1) * t448 + qJD(2) * mrSges(4,2);
t431 = -m(4) * t348 + qJDD(2) * mrSges(4,3) + qJD(2) * t409 + t400 * t447 - t434;
t326 = -mrSges(6,3) * t367 - t369 * t398 + t444;
t429 = mrSges(5,1) * t334 + mrSges(6,1) * t330 - mrSges(5,2) * t335 - mrSges(6,2) * t332 + pkin(4) * t326 + t470 * t366 + t472 * t367 + t468 * t396 + t454 * t397 + t467 * t398;
t406 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t448;
t401 = (-t425 * mrSges(3,1) + t422 * mrSges(3,2)) * qJD(1);
t389 = t435 - t464;
t327 = t442 + t466;
t320 = mrSges(5,2) * t343 + mrSges(6,2) * t337 - mrSges(5,3) * t334 - mrSges(6,3) * t330 - qJ(5) * t326 + t460 * t366 + t473 * t367 + t472 * t396 - t455 * t397 - t467 * t412;
t319 = qJDD(2) * mrSges(4,2) + qJD(2) * t408 + t400 * t448 + t432;
t318 = mrSges(4,2) * t403 - mrSges(4,3) * t402 + (t408 * t425 - t409 * t422) * qJD(1) + t436;
t317 = -mrSges(5,1) * t343 + mrSges(5,3) * t335 - mrSges(6,1) * t337 + mrSges(6,3) * t332 - pkin(4) * t327 + qJ(5) * t443 + (-qJ(5) * t374 - t454) * t412 + t455 * t398 + (-qJ(5) * mrSges(6,2) + t470) * t396 + t460 * t367 + t476 * t366;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t439 - mrSges(2,2) * t437 + t422 * (mrSges(4,1) * t349 + mrSges(3,2) * t389 - mrSges(3,3) * t376 - mrSges(4,3) * t345 + pkin(3) * t321 - qJ(3) * t318 - t451 * qJD(2) + t459 * qJDD(2) + t474 * t402 + t461 * t403 + t452 * t447 + t429) + t425 * (-mrSges(3,1) * t389 + mrSges(3,3) * t377 - mrSges(4,1) * t348 + mrSges(4,2) * t345 - t421 * t320 - t424 * t317 - pkin(3) * (t434 - t465) - pkin(7) * t438 - pkin(2) * t318 + t471 * t403 + t461 * t402 + t458 * qJDD(2) + t450 * qJD(2) - t452 * t448) + pkin(1) * (-m(3) * t389 + t463 * t403 + (-mrSges(3,2) + mrSges(4,3)) * t402 + (t449 * t425 + (-t406 + t409) * t422) * qJD(1) - t436) + pkin(6) * (t425 * (t431 + t401 * t447 - qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t403 + m(3) * t377 - qJD(2) * t406 + t465) + (-m(3) * t376 + t402 * mrSges(3,3) - t463 * qJDD(2) - t449 * qJD(2) + (t400 + t401) * t448 + t432) * t422); mrSges(3,1) * t376 - mrSges(3,2) * t377 + mrSges(4,2) * t349 - mrSges(4,3) * t348 + t424 * t320 - t421 * t317 - pkin(7) * t321 - pkin(2) * t319 + qJ(3) * (-t366 * mrSges(5,1) - t397 * t372 + t431 + t466) + (qJ(3) * mrSges(4,1) + t458) * t403 + t459 * t402 + t469 * qJDD(2) + (t451 * t422 - t450 * t425) * qJD(1); t319; t429; t327;];
tauJ = t1;
