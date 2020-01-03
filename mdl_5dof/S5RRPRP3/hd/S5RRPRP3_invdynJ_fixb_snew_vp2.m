% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRP3
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:04
% EndTime: 2019-12-31 19:51:05
% DurationCPUTime: 0.83s
% Computational Cost: add. (6613->167), mult. (9082->203), div. (0->0), fcn. (5517->8), ass. (0->81)
t401 = Ifges(5,1) + Ifges(6,1);
t394 = Ifges(5,4) - Ifges(6,5);
t393 = Ifges(5,5) + Ifges(6,4);
t400 = -Ifges(5,2) - Ifges(6,3);
t392 = Ifges(5,6) - Ifges(6,6);
t399 = Ifges(5,3) + Ifges(6,2);
t354 = qJD(1) + qJD(2);
t350 = t354 ^ 2;
t357 = sin(pkin(8));
t358 = cos(pkin(8));
t359 = sin(qJ(4));
t397 = cos(qJ(4));
t398 = t357 * t359 - t358 * t397;
t396 = pkin(3) * t358;
t395 = -mrSges(5,3) - mrSges(6,2);
t391 = mrSges(4,2) * t357;
t353 = t358 ^ 2;
t390 = t350 * t353;
t351 = qJDD(1) + qJDD(2);
t389 = t351 * t358;
t361 = sin(qJ(1));
t363 = cos(qJ(1));
t376 = t361 * g(1) - g(2) * t363;
t339 = qJDD(1) * pkin(1) + t376;
t371 = -g(1) * t363 - g(2) * t361;
t340 = -qJD(1) ^ 2 * pkin(1) + t371;
t360 = sin(qJ(2));
t362 = cos(qJ(2));
t325 = t360 * t339 + t362 * t340;
t322 = -pkin(2) * t350 + qJ(3) * t351 + t325;
t381 = qJD(3) * t354;
t375 = -t358 * g(3) - 0.2e1 * t357 * t381;
t298 = (-pkin(7) * t351 + t350 * t396 - t322) * t357 + t375;
t302 = -t357 * g(3) + (t322 + 0.2e1 * t381) * t358;
t299 = -pkin(3) * t390 + pkin(7) * t389 + t302;
t295 = t359 * t298 + t397 * t299;
t368 = t357 * t397 + t358 * t359;
t333 = t368 * t354;
t379 = t333 * qJD(4);
t320 = t351 * t398 + t379;
t329 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t333;
t332 = t398 * t354;
t316 = pkin(4) * t332 - qJ(5) * t333;
t364 = qJD(4) ^ 2;
t290 = -pkin(4) * t364 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t316 * t332 + t295;
t330 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t333;
t378 = m(6) * t290 + qJDD(4) * mrSges(6,3) + qJD(4) * t330;
t317 = mrSges(6,1) * t332 - mrSges(6,3) * t333;
t382 = -mrSges(5,1) * t332 - mrSges(5,2) * t333 - t317;
t285 = m(5) * t295 - qJDD(4) * mrSges(5,2) - qJD(4) * t329 + t320 * t395 + t332 * t382 + t378;
t294 = t298 * t397 - t299 * t359;
t380 = t332 * qJD(4);
t321 = t351 * t368 - t380;
t328 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t332;
t291 = -qJDD(4) * pkin(4) - qJ(5) * t364 + t316 * t333 + qJDD(5) - t294;
t331 = -mrSges(6,2) * t332 + qJD(4) * mrSges(6,3);
t372 = -m(6) * t291 + qJDD(4) * mrSges(6,1) + qJD(4) * t331;
t286 = m(5) * t294 + qJDD(4) * mrSges(5,1) + qJD(4) * t328 + t321 * t395 + t333 * t382 + t372;
t386 = t359 * t285 + t397 * t286;
t385 = -qJD(4) * t399 + t332 * t392 - t333 * t393;
t384 = qJD(4) * t392 + t332 * t400 + t333 * t394;
t383 = qJD(4) * t393 - t332 * t394 + t333 * t401;
t301 = -t357 * t322 + t375;
t369 = mrSges(4,3) * t351 + (-mrSges(4,1) * t358 + t391) * t350;
t373 = t397 * t285 - t359 * t286;
t374 = -t357 * (m(4) * t301 - t357 * t369 + t386) + t358 * (m(4) * t302 + t358 * t369 + t373);
t324 = t362 * t339 - t360 * t340;
t370 = qJDD(3) - t324;
t352 = t357 ^ 2;
t300 = (-pkin(2) - t396) * t351 + (-qJ(3) + (-t352 - t353) * pkin(7)) * t350 + t370;
t293 = -0.2e1 * qJD(5) * t333 + (-t321 + t380) * qJ(5) + (t320 + t379) * pkin(4) + t300;
t287 = m(6) * t293 + t320 * mrSges(6,1) - t321 * mrSges(6,3) - t333 * t330 + t332 * t331;
t276 = -mrSges(5,1) * t300 - mrSges(6,1) * t293 + mrSges(6,2) * t290 + mrSges(5,3) * t295 - pkin(4) * t287 + t383 * qJD(4) + t392 * qJDD(4) + t320 * t400 + t394 * t321 + t385 * t333;
t277 = mrSges(5,2) * t300 + mrSges(6,2) * t291 - mrSges(5,3) * t294 - mrSges(6,3) * t293 - qJ(5) * t287 - t384 * qJD(4) + t393 * qJDD(4) - t394 * t320 + t321 * t401 + t385 * t332;
t319 = -t351 * pkin(2) - t350 * qJ(3) + t370;
t366 = m(5) * t300 + t320 * mrSges(5,1) + t321 * mrSges(5,2) + t332 * t328 + t333 * t329 + t287;
t365 = -m(4) * t319 + mrSges(4,1) * t389 - t366 + (t350 * t352 + t390) * mrSges(4,3);
t281 = t351 * t391 - t365;
t367 = -mrSges(3,2) * t325 + t358 * (-mrSges(4,1) * t319 + mrSges(4,3) * t302 + t359 * t277 + t397 * t276 - pkin(3) * t366 + pkin(7) * t373 + (Ifges(4,4) * t357 + Ifges(4,2) * t358) * t351) + t357 * (mrSges(4,2) * t319 - mrSges(4,3) * t301 + t397 * t277 - t359 * t276 - pkin(7) * t386 + (Ifges(4,1) * t357 + Ifges(4,4) * t358) * t351) + qJ(3) * t374 - pkin(2) * t281 + mrSges(3,1) * t324 + Ifges(3,3) * t351;
t288 = t321 * mrSges(6,2) + t333 * t317 - t372;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t376 - mrSges(2,2) * t371 + pkin(1) * (t360 * (m(3) * t325 - mrSges(3,1) * t350 - mrSges(3,2) * t351 + t374) + t362 * (t365 + (mrSges(3,1) - t391) * t351 + m(3) * t324 - t350 * mrSges(3,2))) + t367; t367; t281; mrSges(5,1) * t294 - mrSges(5,2) * t295 - mrSges(6,1) * t291 + mrSges(6,3) * t290 - pkin(4) * t288 + qJ(5) * t378 + t384 * t333 + (-qJ(5) * t317 + t383) * t332 + t393 * t321 + (-mrSges(6,2) * qJ(5) - t392) * t320 + t399 * qJDD(4); t288;];
tauJ = t1;
