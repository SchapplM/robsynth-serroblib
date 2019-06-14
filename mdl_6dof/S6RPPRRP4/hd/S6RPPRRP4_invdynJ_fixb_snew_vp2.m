% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-05-05 14:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:54:08
% EndTime: 2019-05-05 14:54:10
% DurationCPUTime: 1.26s
% Computational Cost: add. (6804->213), mult. (12044->253), div. (0->0), fcn. (5859->8), ass. (0->93)
t414 = Ifges(6,1) + Ifges(7,1);
t402 = Ifges(6,4) - Ifges(7,5);
t410 = -Ifges(6,5) - Ifges(7,4);
t413 = Ifges(6,2) + Ifges(7,3);
t400 = Ifges(6,6) - Ifges(7,6);
t374 = cos(qJ(5));
t372 = sin(qJ(4));
t395 = t372 * qJD(1);
t404 = sin(qJ(5));
t352 = t374 * qJD(4) + t404 * t395;
t375 = cos(qJ(4));
t393 = qJD(1) * qJD(4);
t390 = t375 * t393;
t356 = -t372 * qJDD(1) - t390;
t327 = t352 * qJD(5) + t404 * qJDD(4) + t374 * t356;
t353 = -t404 * qJD(4) + t374 * t395;
t331 = -t352 * mrSges(7,1) + t353 * mrSges(7,3);
t378 = qJD(1) ^ 2;
t373 = sin(qJ(1));
t376 = cos(qJ(1));
t386 = -t376 * g(1) - t373 * g(2);
t383 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t386;
t405 = -pkin(1) - pkin(2);
t340 = t405 * t378 + t383;
t389 = t373 * g(1) - t376 * g(2);
t382 = -t378 * qJ(2) + qJDD(2) - t389;
t341 = t405 * qJDD(1) + t382;
t370 = sin(pkin(9));
t371 = cos(pkin(9));
t314 = -t370 * t340 + t371 * t341;
t312 = qJDD(1) * pkin(3) - t378 * pkin(7) - t314;
t391 = t372 * t393;
t357 = -t375 * qJDD(1) + t391;
t304 = (-t356 + t390) * pkin(8) + (-t357 - t391) * pkin(4) + t312;
t315 = t371 * t340 + t370 * t341;
t313 = -t378 * pkin(3) - qJDD(1) * pkin(7) + t315;
t368 = g(3) + qJDD(3);
t309 = t375 * t313 + t372 * t368;
t355 = (t375 * pkin(4) + t372 * pkin(8)) * qJD(1);
t377 = qJD(4) ^ 2;
t394 = t375 * qJD(1);
t307 = -t377 * pkin(4) + qJDD(4) * pkin(8) - t355 * t394 + t309;
t301 = t374 * t304 - t404 * t307;
t330 = -t352 * pkin(5) + t353 * qJ(6);
t351 = qJDD(5) - t357;
t361 = qJD(5) + t394;
t360 = t361 ^ 2;
t299 = -t351 * pkin(5) - t360 * qJ(6) - t353 * t330 + qJDD(6) - t301;
t339 = t352 * mrSges(7,2) + t361 * mrSges(7,3);
t387 = -m(7) * t299 + t351 * mrSges(7,1) + t361 * t339;
t295 = t327 * mrSges(7,2) - t353 * t331 - t387;
t302 = t404 * t304 + t374 * t307;
t406 = 2 * qJD(6);
t298 = -t360 * pkin(5) + t351 * qJ(6) + t352 * t330 + t361 * t406 + t302;
t326 = -t353 * qJD(5) - t374 * qJDD(4) + t404 * t356;
t338 = -t361 * mrSges(7,1) - t353 * mrSges(7,2);
t392 = m(7) * t298 + t351 * mrSges(7,3) + t361 * t338;
t397 = t402 * t352 - t414 * t353 - t410 * t361;
t407 = t413 * t352 - t402 * t353 + t400 * t361;
t409 = -Ifges(6,3) - Ifges(7,2);
t412 = -t410 * t327 - t407 * t353 - t400 * t326 - t409 * t351 + mrSges(6,1) * t301 - mrSges(7,1) * t299 - mrSges(6,2) * t302 + mrSges(7,3) * t298 - pkin(5) * t295 + qJ(6) * (-t326 * mrSges(7,2) + t352 * t331 + t392) - t397 * t352;
t403 = -mrSges(6,3) - mrSges(7,2);
t398 = t400 * t352 + t410 * t353 - t409 * t361;
t396 = -t352 * mrSges(6,1) - t353 * mrSges(6,2) + t331;
t354 = (t375 * mrSges(5,1) - t372 * mrSges(5,2)) * qJD(1);
t358 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t395;
t337 = t361 * mrSges(6,1) + t353 * mrSges(6,3);
t292 = m(6) * t302 - t351 * mrSges(6,2) + t403 * t326 - t361 * t337 + t396 * t352 + t392;
t336 = -t361 * mrSges(6,2) + t352 * mrSges(6,3);
t293 = m(6) * t301 + t351 * mrSges(6,1) + t403 * t327 + t361 * t336 + t396 * t353 + t387;
t385 = t374 * t292 - t404 * t293;
t286 = m(5) * t309 - qJDD(4) * mrSges(5,2) + t357 * mrSges(5,3) - qJD(4) * t358 - t354 * t394 + t385;
t308 = -t372 * t313 + t375 * t368;
t359 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t394;
t306 = -qJDD(4) * pkin(4) - t377 * pkin(8) - t355 * t395 - t308;
t300 = t353 * t406 + (-t352 * t361 - t327) * qJ(6) + (-t353 * t361 + t326) * pkin(5) + t306;
t296 = m(7) * t300 + t326 * mrSges(7,1) - t327 * mrSges(7,3) + t353 * t338 - t352 * t339;
t379 = -m(6) * t306 - t326 * mrSges(6,1) - t327 * mrSges(6,2) + t352 * t336 + t353 * t337 - t296;
t290 = m(5) * t308 + qJDD(4) * mrSges(5,1) - t356 * mrSges(5,3) + qJD(4) * t359 + t354 * t395 + t379;
t388 = t375 * t286 - t372 * t290;
t283 = m(4) * t315 - t378 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t388;
t289 = t404 * t292 + t374 * t293;
t380 = -m(5) * t312 + t357 * mrSges(5,1) - t356 * mrSges(5,2) + t358 * t395 - t359 * t394 - t289;
t284 = m(4) * t314 - qJDD(1) * mrSges(4,1) - t378 * mrSges(4,2) + t380;
t384 = t370 * t283 + t371 * t284;
t346 = Ifges(5,5) * qJD(4) + (-t372 * Ifges(5,1) - t375 * Ifges(5,4)) * qJD(1);
t345 = Ifges(5,6) * qJD(4) + (-t372 * Ifges(5,4) - t375 * Ifges(5,2)) * qJD(1);
t343 = -qJDD(1) * pkin(1) + t382;
t342 = -t378 * pkin(1) + t383;
t288 = mrSges(6,2) * t306 + mrSges(7,2) * t299 - mrSges(6,3) * t301 - mrSges(7,3) * t300 - qJ(6) * t296 - t402 * t326 + t414 * t327 - t410 * t351 + t398 * t352 - t407 * t361;
t287 = -mrSges(6,1) * t306 - mrSges(7,1) * t300 + mrSges(7,2) * t298 + mrSges(6,3) * t302 - pkin(5) * t296 - t413 * t326 + t402 * t327 + t400 * t351 + t398 * t353 + t397 * t361;
t282 = m(3) * t343 - qJDD(1) * mrSges(3,1) - t378 * mrSges(3,3) + t384;
t1 = [qJ(2) * (m(3) * t342 - t378 * mrSges(3,1) + t371 * t283 - t370 * t284) + mrSges(2,1) * t389 - pkin(1) * t282 - mrSges(2,2) * t386 - pkin(2) * t384 - mrSges(3,1) * t343 + mrSges(3,3) * t342 - t375 * (-mrSges(5,1) * t312 + mrSges(5,3) * t309 + Ifges(5,4) * t356 + Ifges(5,2) * t357 + Ifges(5,6) * qJDD(4) - pkin(4) * t289 + qJD(4) * t346 - t412) - pkin(3) * t380 - pkin(7) * t388 - mrSges(4,1) * t314 + mrSges(4,2) * t315 - t372 * (mrSges(5,2) * t312 - mrSges(5,3) * t308 + Ifges(5,1) * t356 + Ifges(5,4) * t357 + Ifges(5,5) * qJDD(4) - pkin(8) * t289 - qJD(4) * t345 - t404 * t287 + t374 * t288) + (qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1); t282; m(4) * t368 + t372 * t286 + t375 * t290; Ifges(5,5) * t356 + Ifges(5,6) * t357 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t308 - mrSges(5,2) * t309 + t404 * t288 + t374 * t287 + pkin(4) * t379 + pkin(8) * t385 + (-t372 * t345 + t375 * t346) * qJD(1); t412; t295;];
tauJ  = t1;
