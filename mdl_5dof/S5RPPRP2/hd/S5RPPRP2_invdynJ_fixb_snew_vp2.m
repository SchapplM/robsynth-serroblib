% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:12
% EndTime: 2019-12-31 17:49:13
% DurationCPUTime: 0.79s
% Computational Cost: add. (2881->172), mult. (6260->214), div. (0->0), fcn. (3789->8), ass. (0->84)
t395 = Ifges(5,1) + Ifges(6,1);
t389 = Ifges(5,4) - Ifges(6,5);
t388 = Ifges(5,5) + Ifges(6,4);
t394 = -Ifges(5,2) - Ifges(6,3);
t387 = Ifges(5,6) - Ifges(6,6);
t393 = Ifges(5,3) + Ifges(6,2);
t360 = qJD(1) ^ 2;
t392 = cos(qJ(4));
t354 = cos(pkin(8));
t391 = t354 * pkin(3);
t390 = -mrSges(5,3) - mrSges(6,2);
t352 = sin(pkin(8));
t386 = mrSges(4,2) * t352;
t348 = t354 ^ 2;
t385 = t348 * t360;
t357 = sin(qJ(1));
t358 = cos(qJ(1));
t370 = t357 * g(1) - g(2) * t358;
t335 = qJDD(1) * pkin(1) + t370;
t366 = -g(1) * t358 - g(2) * t357;
t336 = -pkin(1) * t360 + t366;
t353 = sin(pkin(7));
t355 = cos(pkin(7));
t321 = t353 * t335 + t355 * t336;
t308 = -pkin(2) * t360 + qJDD(1) * qJ(3) + t321;
t351 = -g(3) + qJDD(2);
t375 = qJD(1) * qJD(3);
t379 = t354 * t351 - 0.2e1 * t352 * t375;
t295 = (-pkin(6) * qJDD(1) + t360 * t391 - t308) * t352 + t379;
t299 = t352 * t351 + (t308 + 0.2e1 * t375) * t354;
t373 = t354 * qJDD(1);
t296 = -pkin(3) * t385 + pkin(6) * t373 + t299;
t356 = sin(qJ(4));
t292 = t356 * t295 + t392 * t296;
t371 = t354 * t392;
t374 = t352 * qJDD(1);
t363 = t392 * t352 + t354 * t356;
t329 = t363 * qJD(1);
t378 = qJD(4) * t329;
t318 = -qJDD(1) * t371 + t356 * t374 + t378;
t325 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t329;
t376 = t352 * qJD(1);
t328 = -qJD(1) * t371 + t356 * t376;
t312 = pkin(4) * t328 - qJ(5) * t329;
t359 = qJD(4) ^ 2;
t287 = -pkin(4) * t359 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t312 * t328 + t292;
t326 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t329;
t372 = m(6) * t287 + qJDD(4) * mrSges(6,3) + qJD(4) * t326;
t313 = mrSges(6,1) * t328 - mrSges(6,3) * t329;
t380 = -mrSges(5,1) * t328 - mrSges(5,2) * t329 - t313;
t282 = m(5) * t292 - qJDD(4) * mrSges(5,2) - qJD(4) * t325 + t390 * t318 + t380 * t328 + t372;
t291 = t392 * t295 - t356 * t296;
t377 = t328 * qJD(4);
t319 = t363 * qJDD(1) - t377;
t324 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t328;
t288 = -qJDD(4) * pkin(4) - t359 * qJ(5) + t329 * t312 + qJDD(5) - t291;
t327 = -mrSges(6,2) * t328 + qJD(4) * mrSges(6,3);
t367 = -m(6) * t288 + qJDD(4) * mrSges(6,1) + qJD(4) * t327;
t283 = m(5) * t291 + qJDD(4) * mrSges(5,1) + qJD(4) * t324 + t390 * t319 + t380 * t329 + t367;
t384 = t356 * t282 + t392 * t283;
t383 = -t393 * qJD(4) + t387 * t328 - t388 * t329;
t382 = t387 * qJD(4) + t394 * t328 + t389 * t329;
t381 = t388 * qJD(4) - t389 * t328 + t395 * t329;
t298 = -t308 * t352 + t379;
t364 = mrSges(4,3) * qJDD(1) + t360 * (-t354 * mrSges(4,1) + t386);
t276 = m(4) * t298 - t364 * t352 + t384;
t368 = t392 * t282 - t283 * t356;
t277 = m(4) * t299 + t364 * t354 + t368;
t369 = -t352 * t276 + t354 * t277;
t320 = t335 * t355 - t353 * t336;
t365 = qJDD(3) - t320;
t347 = t352 ^ 2;
t297 = (-pkin(2) - t391) * qJDD(1) + (-qJ(3) + (-t347 - t348) * pkin(6)) * t360 + t365;
t290 = -0.2e1 * qJD(5) * t329 + (-t319 + t377) * qJ(5) + (t318 + t378) * pkin(4) + t297;
t284 = m(6) * t290 + t318 * mrSges(6,1) - t319 * mrSges(6,3) - t329 * t326 + t328 * t327;
t362 = m(5) * t297 + t318 * mrSges(5,1) + mrSges(5,2) * t319 + t328 * t324 + t325 * t329 + t284;
t301 = -qJDD(1) * pkin(2) - qJ(3) * t360 + t365;
t361 = -m(4) * t301 + mrSges(4,1) * t373 - t362 + (t347 * t360 + t385) * mrSges(4,3);
t334 = (Ifges(4,5) * t352 + Ifges(4,6) * t354) * qJD(1);
t285 = mrSges(6,2) * t319 + t313 * t329 - t367;
t278 = mrSges(4,2) * t374 - t361;
t274 = mrSges(5,2) * t297 + mrSges(6,2) * t288 - mrSges(5,3) * t291 - mrSges(6,3) * t290 - qJ(5) * t284 - t382 * qJD(4) + t388 * qJDD(4) - t389 * t318 + t395 * t319 + t383 * t328;
t273 = -mrSges(5,1) * t297 - mrSges(6,1) * t290 + mrSges(6,2) * t287 + mrSges(5,3) * t292 - pkin(4) * t284 + t381 * qJD(4) + t387 * qJDD(4) + t394 * t318 + t389 * t319 + t383 * t329;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t370 - mrSges(2,2) * t366 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t320 - mrSges(3,2) * t321 + t352 * (t354 * qJD(1) * t334 + mrSges(4,2) * t301 - mrSges(4,3) * t298 + t392 * t274 - t356 * t273 - pkin(6) * t384 + (Ifges(4,1) * t352 + Ifges(4,4) * t354) * qJDD(1)) + t354 * (-t334 * t376 - mrSges(4,1) * t301 + mrSges(4,3) * t299 + t356 * t274 + t392 * t273 - pkin(3) * t362 + pkin(6) * t368 + (Ifges(4,4) * t352 + Ifges(4,2) * t354) * qJDD(1)) - pkin(2) * t278 + qJ(3) * t369 + pkin(1) * (t353 * (m(3) * t321 - mrSges(3,1) * t360 - qJDD(1) * mrSges(3,2) + t369) + t355 * (t361 + (mrSges(3,1) - t386) * qJDD(1) + m(3) * t320 - mrSges(3,2) * t360)); m(3) * t351 + t276 * t354 + t277 * t352; t278; mrSges(5,1) * t291 - mrSges(5,2) * t292 - mrSges(6,1) * t288 + mrSges(6,3) * t287 - pkin(4) * t285 + qJ(5) * t372 + t382 * t329 + (-qJ(5) * t313 + t381) * t328 + t388 * t319 + (-qJ(5) * mrSges(6,2) - t387) * t318 + t393 * qJDD(4); t285;];
tauJ = t1;
