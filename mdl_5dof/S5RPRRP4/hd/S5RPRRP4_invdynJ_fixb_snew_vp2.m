% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP4
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:33
% EndTime: 2020-01-03 11:49:37
% DurationCPUTime: 1.78s
% Computational Cost: add. (7217->220), mult. (17303->279), div. (0->0), fcn. (11278->8), ass. (0->98)
t456 = Ifges(5,4) + Ifges(6,4);
t465 = Ifges(5,2) + Ifges(6,2);
t460 = Ifges(5,6) + Ifges(6,6);
t427 = qJD(1) ^ 2;
t423 = sin(qJ(1));
t426 = cos(qJ(1));
t437 = -g(2) * t423 + t426 * g(3);
t464 = -pkin(1) * t427 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t437;
t421 = sin(qJ(4));
t422 = sin(qJ(3));
t424 = cos(qJ(4));
t425 = cos(qJ(3));
t419 = sin(pkin(8));
t450 = qJD(1) * t419;
t388 = (-t425 * t421 - t422 * t424) * t450;
t446 = qJD(1) * qJD(3);
t396 = (-qJDD(1) * t422 - t425 * t446) * t419;
t397 = (qJDD(1) * t425 - t422 * t446) * t419;
t359 = qJD(4) * t388 + t396 * t421 + t397 * t424;
t389 = (-t422 * t421 + t425 * t424) * t450;
t370 = -mrSges(6,1) * t388 + mrSges(6,2) * t389;
t420 = cos(pkin(8));
t381 = -g(1) * t419 + t464 * t420;
t434 = -pkin(2) * t420 - pkin(6) * t419;
t403 = t434 * qJD(1);
t448 = t420 * qJD(1);
t369 = t403 * t448 + t381;
t451 = -t426 * g(2) - t423 * g(3);
t431 = -qJ(2) * t427 + qJDD(2) - t451;
t382 = (-pkin(1) + t434) * qJDD(1) + t431;
t379 = t425 * t382;
t445 = qJDD(1) * t420;
t409 = qJDD(3) - t445;
t410 = qJD(3) - t448;
t417 = t419 ^ 2;
t455 = t417 * t427;
t344 = pkin(3) * t409 - pkin(7) * t397 + t379 + (-pkin(3) * t425 * t455 - pkin(7) * t410 * t450 - t369) * t422;
t348 = t425 * t369 + t422 * t382;
t439 = t425 * t450;
t395 = pkin(3) * t410 - pkin(7) * t439;
t444 = t422 ^ 2 * t455;
t345 = -pkin(3) * t444 + pkin(7) * t396 - t395 * t410 + t348;
t337 = t424 * t344 - t345 * t421;
t407 = qJDD(4) + t409;
t408 = qJD(4) + t410;
t333 = -0.2e1 * qJD(5) * t389 + (t388 * t408 - t359) * qJ(5) + (t388 * t389 + t407) * pkin(4) + t337;
t373 = -mrSges(6,2) * t408 + mrSges(6,3) * t388;
t443 = m(6) * t333 + t407 * mrSges(6,1) + t408 * t373;
t329 = -mrSges(6,3) * t359 - t370 * t389 + t443;
t338 = t421 * t344 + t424 * t345;
t358 = -qJD(4) * t389 + t396 * t424 - t397 * t421;
t375 = pkin(4) * t408 - qJ(5) * t389;
t387 = t388 ^ 2;
t335 = -pkin(4) * t387 + qJ(5) * t358 + 0.2e1 * qJD(5) * t388 - t375 * t408 + t338;
t461 = Ifges(5,5) + Ifges(6,5);
t462 = Ifges(5,1) + Ifges(6,1);
t453 = t456 * t388 + t462 * t389 + t461 * t408;
t458 = t465 * t388 + t456 * t389 + t460 * t408;
t459 = Ifges(5,3) + Ifges(6,3);
t463 = mrSges(5,1) * t337 + mrSges(6,1) * t333 - mrSges(5,2) * t338 - mrSges(6,2) * t335 + pkin(4) * t329 + t460 * t358 + t461 * t359 - t453 * t388 + t458 * t389 + t459 * t407;
t371 = -mrSges(5,1) * t388 + mrSges(5,2) * t389;
t374 = -mrSges(5,2) * t408 + mrSges(5,3) * t388;
t324 = m(5) * t337 + mrSges(5,1) * t407 + t374 * t408 + (-t370 - t371) * t389 + (-mrSges(5,3) - mrSges(6,3)) * t359 + t443;
t376 = mrSges(6,1) * t408 - mrSges(6,3) * t389;
t377 = mrSges(5,1) * t408 - mrSges(5,3) * t389;
t442 = m(6) * t335 + t358 * mrSges(6,3) + t388 * t370;
t327 = m(5) * t338 + mrSges(5,3) * t358 + t371 * t388 + (-t376 - t377) * t408 + (-mrSges(5,2) - mrSges(6,2)) * t407 + t442;
t322 = t424 * t324 + t421 * t327;
t347 = -t369 * t422 + t379;
t457 = -mrSges(4,1) * t347 + mrSges(4,2) * t348 - Ifges(4,5) * t397 - Ifges(4,6) * t396 - Ifges(4,3) * t409 - pkin(3) * t322 - t463;
t380 = -t420 * g(1) - t464 * t419;
t454 = -t460 * t388 - t461 * t389 - t459 * t408;
t449 = qJDD(1) * mrSges(3,3);
t368 = t403 * t450 - t380;
t346 = -pkin(3) * t396 - pkin(7) * t444 + t395 * t439 + t368;
t340 = -pkin(4) * t358 - qJ(5) * t387 + t375 * t389 + qJDD(5) + t346;
t441 = m(6) * t340 + t359 * mrSges(6,2) + t389 * t376;
t440 = t422 * t450;
t435 = -t324 * t421 + t424 * t327;
t433 = -mrSges(3,1) * t420 + mrSges(3,2) * t419;
t392 = -mrSges(4,2) * t410 - mrSges(4,3) * t440;
t394 = (t422 * mrSges(4,1) + t425 * mrSges(4,2)) * t450;
t319 = m(4) * t347 + mrSges(4,1) * t409 - mrSges(4,3) * t397 + t392 * t410 - t394 * t439 + t322;
t393 = mrSges(4,1) * t410 - mrSges(4,3) * t439;
t320 = m(4) * t348 - mrSges(4,2) * t409 + mrSges(4,3) * t396 - t393 * t410 - t394 * t440 + t435;
t317 = t319 * t425 + t320 * t422;
t384 = Ifges(4,6) * t410 + (t425 * Ifges(4,4) - t422 * Ifges(4,2)) * t450;
t385 = Ifges(4,5) * t410 + (t425 * Ifges(4,1) - t422 * Ifges(4,4)) * t450;
t432 = t425 * t384 + t422 * t385;
t429 = m(5) * t346 + t359 * mrSges(5,2) + t389 * t377 + (-t373 - t374) * t388 + (-mrSges(5,1) - mrSges(6,1)) * t358 + t441;
t402 = (Ifges(3,5) * t419 + Ifges(3,6) * t420) * qJD(1);
t401 = t433 * qJD(1);
t399 = -qJDD(1) * pkin(1) + t431;
t330 = -t358 * mrSges(6,1) - t388 * t373 + t441;
t321 = mrSges(5,2) * t346 + mrSges(6,2) * t340 - mrSges(5,3) * t337 - mrSges(6,3) * t333 - qJ(5) * t329 + t456 * t358 + t462 * t359 - t454 * t388 + t461 * t407 - t458 * t408;
t318 = -mrSges(5,1) * t346 + mrSges(5,3) * t338 - mrSges(6,1) * t340 + mrSges(6,3) * t335 - pkin(4) * t330 + qJ(5) * t442 + (-qJ(5) * t376 + t453) * t408 + (-qJ(5) * mrSges(6,2) + t460) * t407 + t454 * t389 + t456 * t359 + t465 * t358;
t316 = m(3) * t399 + t433 * qJDD(1) + (-t420 ^ 2 - t417) * t427 * mrSges(3,3) + t317;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t451 - mrSges(2,2) * t437 + t419 * (mrSges(3,2) * t399 - mrSges(3,3) * t380 + t425 * (mrSges(4,2) * t368 - mrSges(4,3) * t347 + Ifges(4,1) * t397 + Ifges(4,4) * t396 + Ifges(4,5) * t409 - pkin(7) * t322 - t318 * t421 + t321 * t424 - t384 * t410) - t422 * (-mrSges(4,1) * t368 + mrSges(4,3) * t348 + Ifges(4,4) * t397 + Ifges(4,2) * t396 + Ifges(4,6) * t409 - pkin(3) * t429 + pkin(7) * t435 + t424 * t318 + t421 * t321 + t410 * t385) - pkin(6) * t317 + (Ifges(3,1) * t419 + Ifges(3,4) * t420) * qJDD(1) + t402 * t448) + t420 * ((Ifges(3,4) * qJDD(1) + (-t402 - t432) * qJD(1)) * t419 + Ifges(3,2) * t445 - mrSges(3,1) * t399 + mrSges(3,3) * t381 - pkin(2) * t317 + t457) - pkin(1) * t316 + qJ(2) * ((m(3) * t381 - t422 * t319 + t425 * t320 + (qJD(1) * t401 + t449) * t420) * t420 + (t419 * t449 - m(3) * t380 + m(4) * t368 - t396 * mrSges(4,1) + t397 * mrSges(4,2) + (t392 * t422 + t393 * t425 + t401) * t450 + t429) * t419); t316; t432 * t450 - t457; t463; t330;];
tauJ = t1;
