% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRR12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR12_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR12_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:06
% EndTime: 2019-12-31 18:07:07
% DurationCPUTime: 0.94s
% Computational Cost: add. (5189->184), mult. (11507->235), div. (0->0), fcn. (7430->8), ass. (0->86)
t363 = qJD(1) ^ 2;
t358 = sin(qJ(1));
t361 = cos(qJ(1));
t379 = t358 * g(1) - t361 * g(2);
t369 = -t363 * qJ(2) + qJDD(2) - t379;
t389 = -pkin(1) - qJ(3);
t393 = -(2 * qJD(1) * qJD(3)) + t389 * qJDD(1) + t369;
t355 = cos(pkin(8));
t392 = t355 ^ 2;
t375 = -t361 * g(1) - t358 * g(2);
t391 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t375;
t354 = sin(pkin(8));
t357 = sin(qJ(4));
t360 = cos(qJ(4));
t373 = t354 * t360 + t355 * t357;
t339 = t373 * qJD(1);
t372 = -t354 * t357 + t355 * t360;
t340 = t372 * qJD(1);
t383 = t340 * qJD(4);
t325 = -t373 * qJDD(1) - t383;
t390 = pkin(3) * t363;
t388 = t355 * mrSges(4,2);
t322 = t354 * g(3) + t393 * t355;
t309 = (-pkin(6) * qJDD(1) - t354 * t390) * t355 + t322;
t323 = -t355 * g(3) + t393 * t354;
t351 = t354 ^ 2;
t381 = t354 * qJDD(1);
t310 = -pkin(6) * t381 - t351 * t390 + t323;
t298 = t357 * t309 + t360 * t310;
t319 = t339 * mrSges(5,1) + t340 * mrSges(5,2);
t334 = qJD(4) * mrSges(5,1) - t340 * mrSges(5,3);
t324 = t339 * pkin(4) - t340 * pkin(7);
t362 = qJD(4) ^ 2;
t295 = -t362 * pkin(4) + qJDD(4) * pkin(7) - t339 * t324 + t298;
t368 = qJDD(3) + t391;
t386 = -t351 - t392;
t314 = pkin(3) * t381 + (t386 * pkin(6) + t389) * t363 + t368;
t384 = t339 * qJD(4);
t326 = t372 * qJDD(1) - t384;
t296 = (-t326 + t384) * pkin(7) + (-t325 + t383) * pkin(4) + t314;
t356 = sin(qJ(5));
t359 = cos(qJ(5));
t292 = -t356 * t295 + t359 * t296;
t328 = t359 * qJD(4) - t356 * t340;
t305 = t328 * qJD(5) + t356 * qJDD(4) + t359 * t326;
t329 = t356 * qJD(4) + t359 * t340;
t307 = -t328 * mrSges(6,1) + t329 * mrSges(6,2);
t337 = qJD(5) + t339;
t311 = -t337 * mrSges(6,2) + t328 * mrSges(6,3);
t321 = qJDD(5) - t325;
t290 = m(6) * t292 + t321 * mrSges(6,1) - t305 * mrSges(6,3) - t329 * t307 + t337 * t311;
t293 = t359 * t295 + t356 * t296;
t304 = -t329 * qJD(5) + t359 * qJDD(4) - t356 * t326;
t312 = t337 * mrSges(6,1) - t329 * mrSges(6,3);
t291 = m(6) * t293 - t321 * mrSges(6,2) + t304 * mrSges(6,3) + t328 * t307 - t337 * t312;
t376 = -t356 * t290 + t359 * t291;
t281 = m(5) * t298 - qJDD(4) * mrSges(5,2) + t325 * mrSges(5,3) - qJD(4) * t334 - t339 * t319 + t376;
t297 = t360 * t309 - t357 * t310;
t333 = -qJD(4) * mrSges(5,2) - t339 * mrSges(5,3);
t294 = -qJDD(4) * pkin(4) - t362 * pkin(7) + t340 * t324 - t297;
t367 = -m(6) * t294 + t304 * mrSges(6,1) - t305 * mrSges(6,2) + t328 * t311 - t329 * t312;
t286 = m(5) * t297 + qJDD(4) * mrSges(5,1) - t326 * mrSges(5,3) + qJD(4) * t333 - t340 * t319 + t367;
t387 = t357 * t281 + t360 * t286;
t282 = t359 * t290 + t356 * t291;
t378 = t386 * mrSges(4,3);
t377 = t360 * t281 - t357 * t286;
t371 = -qJDD(1) * mrSges(4,3) - t363 * (t354 * mrSges(4,1) + t388);
t374 = t355 * (m(4) * t322 + t371 * t355 + t387) + t354 * (m(4) * t323 + t371 * t354 + t377);
t366 = m(5) * t314 - t325 * mrSges(5,1) + t326 * mrSges(5,2) + t339 * t333 + t340 * t334 + t282;
t332 = t389 * t363 + t368;
t365 = m(4) * t332 + mrSges(4,1) * t381 + qJDD(1) * t388 + t366;
t300 = Ifges(6,4) * t329 + Ifges(6,2) * t328 + Ifges(6,6) * t337;
t301 = Ifges(6,1) * t329 + Ifges(6,4) * t328 + Ifges(6,5) * t337;
t364 = mrSges(6,1) * t292 - mrSges(6,2) * t293 + Ifges(6,5) * t305 + Ifges(6,6) * t304 + Ifges(6,3) * t321 + t329 * t300 - t328 * t301;
t338 = -qJDD(1) * pkin(1) + t369;
t336 = t363 * pkin(1) - t391;
t317 = Ifges(5,1) * t340 - Ifges(5,4) * t339 + Ifges(5,5) * qJD(4);
t316 = Ifges(5,4) * t340 - Ifges(5,2) * t339 + Ifges(5,6) * qJD(4);
t315 = Ifges(5,5) * t340 - Ifges(5,6) * t339 + Ifges(5,3) * qJD(4);
t299 = Ifges(6,5) * t329 + Ifges(6,6) * t328 + Ifges(6,3) * t337;
t284 = mrSges(6,2) * t294 - mrSges(6,3) * t292 + Ifges(6,1) * t305 + Ifges(6,4) * t304 + Ifges(6,5) * t321 + t328 * t299 - t337 * t300;
t283 = -mrSges(6,1) * t294 + mrSges(6,3) * t293 + Ifges(6,4) * t305 + Ifges(6,2) * t304 + Ifges(6,6) * t321 - t329 * t299 + t337 * t301;
t276 = -mrSges(5,1) * t314 + mrSges(5,3) * t298 + Ifges(5,4) * t326 + Ifges(5,2) * t325 + Ifges(5,6) * qJDD(4) - pkin(4) * t282 + qJD(4) * t317 - t340 * t315 - t364;
t275 = mrSges(5,2) * t314 - mrSges(5,3) * t297 + Ifges(5,1) * t326 + Ifges(5,4) * t325 + Ifges(5,5) * qJDD(4) - pkin(7) * t282 - qJD(4) * t316 - t356 * t283 + t359 * t284 - t339 * t315;
t274 = m(3) * t338 + qJDD(1) * mrSges(3,2) - t363 * mrSges(3,3) + t374;
t1 = [mrSges(2,1) * t379 - mrSges(2,2) * t375 + mrSges(3,2) * t338 - mrSges(3,3) * t336 + t355 * (mrSges(4,2) * t332 - mrSges(4,3) * t322 - pkin(6) * t387 + t360 * t275 - t357 * t276) - t354 * (-mrSges(4,1) * t332 + mrSges(4,3) * t323 - pkin(3) * t366 + pkin(6) * t377 + t357 * t275 + t360 * t276) - qJ(3) * t374 - pkin(1) * t274 + qJ(2) * (-m(3) * t336 + (mrSges(3,2) + t378) * t363 + t365) + (Ifges(4,1) * t392 + qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (-0.2e1 * Ifges(4,4) * t355 + Ifges(4,2) * t354) * t354) * qJDD(1); t274; t363 * t378 + t365; mrSges(5,1) * t297 - mrSges(5,2) * t298 + Ifges(5,5) * t326 + Ifges(5,6) * t325 + Ifges(5,3) * qJDD(4) + pkin(4) * t367 + pkin(7) * t376 + t359 * t283 + t356 * t284 + t340 * t316 + t339 * t317; t364;];
tauJ = t1;
