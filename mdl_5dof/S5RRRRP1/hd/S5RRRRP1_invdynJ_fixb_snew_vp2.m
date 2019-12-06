% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:45:09
% EndTime: 2019-12-05 18:45:13
% DurationCPUTime: 2.40s
% Computational Cost: add. (15911->256), mult. (35165->313), div. (0->0), fcn. (24393->8), ass. (0->102)
t457 = Ifges(5,4) + Ifges(6,4);
t468 = Ifges(5,2) + Ifges(6,2);
t464 = Ifges(5,6) + Ifges(6,6);
t428 = sin(qJ(3));
t429 = sin(qJ(2));
t432 = cos(qJ(3));
t433 = cos(qJ(2));
t407 = (t433 * t428 + t429 * t432) * qJD(1);
t450 = qJD(1) * qJD(2);
t412 = t429 * qJDD(1) + t433 * t450;
t413 = t433 * qJDD(1) - t429 * t450;
t380 = -t407 * qJD(3) - t428 * t412 + t432 * t413;
t452 = qJD(1) * t429;
t416 = qJD(2) * pkin(2) - pkin(7) * t452;
t426 = t433 ^ 2;
t435 = qJD(1) ^ 2;
t430 = sin(qJ(1));
t434 = cos(qJ(1));
t445 = t430 * g(1) - t434 * g(2);
t440 = -qJDD(1) * pkin(1) - t445;
t382 = -t413 * pkin(2) + t416 * t452 + (-pkin(7) * t426 - pkin(6)) * t435 + t440;
t425 = qJD(2) + qJD(3);
t400 = t425 * pkin(3) - t407 * pkin(8);
t406 = (-t429 * t428 + t433 * t432) * qJD(1);
t402 = t406 ^ 2;
t345 = -t380 * pkin(3) - t402 * pkin(8) + t407 * t400 + t382;
t381 = t406 * qJD(3) + t432 * t412 + t428 * t413;
t427 = sin(qJ(4));
t431 = cos(qJ(4));
t394 = t427 * t406 + t431 * t407;
t353 = -t394 * qJD(4) + t431 * t380 - t427 * t381;
t393 = t431 * t406 - t427 * t407;
t354 = t393 * qJD(4) + t427 * t380 + t431 * t381;
t422 = qJD(4) + t425;
t383 = -t422 * mrSges(6,2) + t393 * mrSges(6,3);
t384 = -t422 * mrSges(5,2) + t393 * mrSges(5,3);
t387 = t422 * mrSges(5,1) - t394 * mrSges(5,3);
t385 = t422 * pkin(4) - t394 * qJ(5);
t392 = t393 ^ 2;
t338 = -t353 * pkin(4) - t392 * qJ(5) + t394 * t385 + qJDD(5) + t345;
t386 = t422 * mrSges(6,1) - t394 * mrSges(6,3);
t447 = m(6) * t338 + t354 * mrSges(6,2) + t394 * t386;
t467 = m(5) * t345 + t354 * mrSges(5,2) + t394 * t387 + t447 - (t384 + t383) * t393 - (mrSges(5,1) + mrSges(6,1)) * t353;
t466 = Ifges(5,1) + Ifges(6,1);
t465 = Ifges(5,5) + Ifges(6,5);
t463 = Ifges(5,3) + Ifges(6,3);
t462 = t468 * t393 + t457 * t394 + t464 * t422;
t398 = -t425 * mrSges(4,2) + t406 * mrSges(4,3);
t399 = t425 * mrSges(4,1) - t407 * mrSges(4,3);
t461 = m(4) * t382 - t380 * mrSges(4,1) + t381 * mrSges(4,2) - t406 * t398 + t407 * t399 + t467;
t459 = pkin(2) * t435;
t441 = -t434 * g(1) - t430 * g(2);
t409 = -t435 * pkin(1) + qJDD(1) * pkin(6) + t441;
t456 = t429 * t409;
t375 = qJDD(2) * pkin(2) - t412 * pkin(7) - t456 + (pkin(7) * t450 + t429 * t459 - g(3)) * t433;
t397 = -t429 * g(3) + t433 * t409;
t376 = t413 * pkin(7) - qJD(2) * t416 - t426 * t459 + t397;
t359 = t432 * t375 - t428 * t376;
t424 = qJDD(2) + qJDD(3);
t341 = (t406 * t425 - t381) * pkin(8) + (t406 * t407 + t424) * pkin(3) + t359;
t360 = t428 * t375 + t432 * t376;
t343 = -t402 * pkin(3) + t380 * pkin(8) - t425 * t400 + t360;
t335 = t431 * t341 - t427 * t343;
t369 = -t393 * mrSges(6,1) + t394 * mrSges(6,2);
t370 = -t393 * mrSges(5,1) + t394 * mrSges(5,2);
t421 = qJDD(4) + t424;
t331 = -0.2e1 * qJD(5) * t394 + (t393 * t422 - t354) * qJ(5) + (t393 * t394 + t421) * pkin(4) + t335;
t449 = m(6) * t331 + t421 * mrSges(6,1) + t422 * t383;
t322 = m(5) * t335 + t421 * mrSges(5,1) + t422 * t384 + (-t369 - t370) * t394 + (-mrSges(5,3) - mrSges(6,3)) * t354 + t449;
t336 = t427 * t341 + t431 * t343;
t333 = -t392 * pkin(4) + t353 * qJ(5) + 0.2e1 * qJD(5) * t393 - t422 * t385 + t336;
t448 = m(6) * t333 + t353 * mrSges(6,3) + t393 * t369;
t325 = m(5) * t336 + t353 * mrSges(5,3) + t393 * t370 + (-t386 - t387) * t422 + (-mrSges(5,2) - mrSges(6,2)) * t421 + t448;
t320 = t431 * t322 + t427 * t325;
t395 = -t406 * mrSges(4,1) + t407 * mrSges(4,2);
t316 = m(4) * t359 + t424 * mrSges(4,1) - t381 * mrSges(4,3) - t407 * t395 + t425 * t398 + t320;
t442 = -t427 * t322 + t431 * t325;
t317 = m(4) * t360 - t424 * mrSges(4,2) + t380 * mrSges(4,3) + t406 * t395 - t425 * t399 + t442;
t311 = t432 * t316 + t428 * t317;
t455 = -t464 * t393 - t465 * t394 - t463 * t422;
t454 = -t457 * t393 - t466 * t394 - t465 * t422;
t451 = qJD(1) * t433;
t443 = -t428 * t316 + t432 * t317;
t327 = -t354 * mrSges(6,3) - t394 * t369 + t449;
t437 = mrSges(5,1) * t335 + mrSges(6,1) * t331 - mrSges(5,2) * t336 - mrSges(6,2) * t333 + pkin(4) * t327 + t464 * t353 + t465 * t354 + t454 * t393 + t462 * t394 + t463 * t421;
t389 = Ifges(4,4) * t407 + Ifges(4,2) * t406 + Ifges(4,6) * t425;
t390 = Ifges(4,1) * t407 + Ifges(4,4) * t406 + Ifges(4,5) * t425;
t436 = mrSges(4,1) * t359 - mrSges(4,2) * t360 + Ifges(4,5) * t381 + Ifges(4,6) * t380 + Ifges(4,3) * t424 + pkin(3) * t320 + t407 * t389 - t406 * t390 + t437;
t415 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t451;
t414 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t452;
t411 = (-t433 * mrSges(3,1) + t429 * mrSges(3,2)) * qJD(1);
t408 = -t435 * pkin(6) + t440;
t405 = Ifges(3,5) * qJD(2) + (t429 * Ifges(3,1) + t433 * Ifges(3,4)) * qJD(1);
t404 = Ifges(3,6) * qJD(2) + (t429 * Ifges(3,4) + t433 * Ifges(3,2)) * qJD(1);
t396 = -t433 * g(3) - t456;
t388 = Ifges(4,5) * t407 + Ifges(4,6) * t406 + Ifges(4,3) * t425;
t328 = -t353 * mrSges(6,1) - t393 * t383 + t447;
t318 = mrSges(5,2) * t345 + mrSges(6,2) * t338 - mrSges(5,3) * t335 - mrSges(6,3) * t331 - qJ(5) * t327 + t457 * t353 + t466 * t354 - t455 * t393 + t465 * t421 - t462 * t422;
t312 = -mrSges(5,1) * t345 + mrSges(5,3) * t336 - mrSges(6,1) * t338 + mrSges(6,3) * t333 - pkin(4) * t328 + qJ(5) * t448 + (-qJ(5) * t386 - t454) * t422 + (-qJ(5) * mrSges(6,2) + t464) * t421 + t455 * t394 + t457 * t354 + t468 * t353;
t310 = mrSges(4,2) * t382 - mrSges(4,3) * t359 + Ifges(4,1) * t381 + Ifges(4,4) * t380 + Ifges(4,5) * t424 - pkin(8) * t320 - t427 * t312 + t431 * t318 + t406 * t388 - t425 * t389;
t309 = -mrSges(4,1) * t382 + mrSges(4,3) * t360 + Ifges(4,4) * t381 + Ifges(4,2) * t380 + Ifges(4,6) * t424 - pkin(3) * t467 + pkin(8) * t442 + t431 * t312 + t427 * t318 - t407 * t388 + t425 * t390;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t445 - mrSges(2,2) * t441 + t429 * (mrSges(3,2) * t408 - mrSges(3,3) * t396 + Ifges(3,1) * t412 + Ifges(3,4) * t413 + Ifges(3,5) * qJDD(2) - pkin(7) * t311 - qJD(2) * t404 - t428 * t309 + t432 * t310) + t433 * (-mrSges(3,1) * t408 + mrSges(3,3) * t397 + Ifges(3,4) * t412 + Ifges(3,2) * t413 + Ifges(3,6) * qJDD(2) - pkin(2) * t461 + pkin(7) * t443 + qJD(2) * t405 + t432 * t309 + t428 * t310) + pkin(1) * (-m(3) * t408 - t412 * mrSges(3,2) + t413 * mrSges(3,1) + (-t429 * t414 + t433 * t415) * qJD(1) - t461) + pkin(6) * (t433 * (m(3) * t397 - qJDD(2) * mrSges(3,2) + t413 * mrSges(3,3) - qJD(2) * t414 + t411 * t451 + t443) - t429 * (m(3) * t396 + qJDD(2) * mrSges(3,1) - t412 * mrSges(3,3) + qJD(2) * t415 - t411 * t452 + t311)); t436 + Ifges(3,3) * qJDD(2) + (t429 * t404 - t433 * t405) * qJD(1) + Ifges(3,5) * t412 + Ifges(3,6) * t413 + mrSges(3,1) * t396 - mrSges(3,2) * t397 + pkin(2) * t311; t436; t437; t328;];
tauJ = t1;
