% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRP1
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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:45
% EndTime: 2019-12-05 17:35:48
% DurationCPUTime: 1.30s
% Computational Cost: add. (2323->164), mult. (4820->209), div. (0->0), fcn. (2671->8), ass. (0->80)
t351 = sin(qJ(1));
t353 = cos(qJ(1));
t376 = t353 * g(2) + t351 * g(3);
t329 = qJDD(1) * pkin(1) + t376;
t354 = qJD(1) ^ 2;
t364 = t351 * g(2) - g(3) * t353;
t330 = -pkin(1) * t354 + t364;
t347 = sin(pkin(7));
t349 = cos(pkin(7));
t304 = t347 * t329 + t349 * t330;
t396 = -pkin(2) * t354 + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t304;
t350 = sin(qJ(4));
t352 = cos(qJ(4));
t386 = Ifges(5,4) + Ifges(6,4);
t395 = t352 * (Ifges(5,1) + Ifges(6,1)) - t350 * t386;
t394 = t352 * t386 - t350 * (Ifges(5,2) + Ifges(6,2));
t385 = Ifges(5,5) + Ifges(6,5);
t384 = Ifges(5,6) + Ifges(6,6);
t348 = cos(pkin(8));
t374 = qJD(1) * t348;
t333 = qJD(4) - t374;
t346 = sin(pkin(8));
t375 = qJD(1) * t346;
t363 = (t384 * t333 + t394 * t375) * t352;
t378 = t385 * t333 + t395 * t375;
t392 = t378 * t350 + t363;
t345 = -g(1) + qJDD(2);
t294 = t345 * t348 - t396 * t346;
t318 = (t350 * mrSges(6,1) + t352 * mrSges(6,2)) * t375;
t371 = qJD(1) * qJD(4);
t321 = (qJDD(1) * t352 - t350 * t371) * t346;
t366 = t352 * t375;
t295 = t346 * t345 + t396 * t348;
t359 = -t348 * pkin(3) - t346 * pkin(6);
t328 = t359 * qJD(1);
t293 = t328 * t374 + t295;
t303 = t329 * t349 - t347 * t330;
t357 = -qJ(3) * t354 + qJDD(3) - t303;
t298 = (-pkin(2) + t359) * qJDD(1) + t357;
t297 = t352 * t298;
t370 = qJDD(1) * t348;
t332 = qJDD(4) - t370;
t360 = -0.2e1 * qJD(5) * t375;
t381 = t346 ^ 2 * t354;
t285 = t352 * t360 + pkin(4) * t332 - qJ(5) * t321 + t297 + (-pkin(4) * t352 * t381 - qJ(5) * t333 * t375 - t293) * t350;
t367 = t350 * t375;
t313 = -mrSges(6,2) * t333 - mrSges(6,3) * t367;
t368 = m(6) * t285 + t332 * mrSges(6,1) + t333 * t313;
t282 = -mrSges(6,3) * t321 - t318 * t366 + t368;
t289 = t352 * t293 + t350 * t298;
t315 = pkin(4) * t333 - qJ(5) * t366;
t320 = (-qJDD(1) * t350 - t352 * t371) * t346;
t369 = t350 ^ 2 * t381;
t287 = -pkin(4) * t369 + qJ(5) * t320 - t315 * t333 + t350 * t360 + t289;
t288 = -t293 * t350 + t297;
t389 = mrSges(5,1) * t288 + mrSges(6,1) * t285 - mrSges(5,2) * t289 - mrSges(6,2) * t287 + pkin(4) * t282 + t384 * t320 + t385 * t321 + (Ifges(5,3) + Ifges(6,3)) * t332;
t292 = t328 * t375 - t294;
t290 = -pkin(4) * t320 - qJ(5) * t369 + t315 * t366 + qJDD(5) + t292;
t388 = m(6) * t290;
t387 = -mrSges(5,2) - mrSges(6,2);
t382 = t346 * mrSges(4,2);
t380 = m(6) * t287 + t320 * mrSges(6,3);
t316 = mrSges(6,1) * t333 - mrSges(6,3) * t366;
t377 = -mrSges(5,1) * t333 + mrSges(5,3) * t366 - t316;
t373 = qJDD(1) * mrSges(4,3);
t314 = -mrSges(5,2) * t333 - mrSges(5,3) * t367;
t358 = (-t318 - (t350 * mrSges(5,1) + t352 * mrSges(5,2)) * t375) * t375;
t279 = m(5) * t288 + mrSges(5,1) * t332 + t314 * t333 + (-mrSges(5,3) - mrSges(6,3)) * t321 + t352 * t358 + t368;
t280 = m(5) * t289 + mrSges(5,3) * t320 + t387 * t332 + t377 * t333 + t350 * t358 + t380;
t326 = (-t348 * mrSges(4,1) + t382) * qJD(1);
t276 = m(4) * t295 - t279 * t350 + t280 * t352 + (qJD(1) * t326 + t373) * t348;
t281 = m(4) * t294 - m(5) * t292 - t388 + t387 * t321 + (mrSges(5,1) + mrSges(6,1)) * t320 + (-t373 + (-t326 + t377 * t352 + (-t313 - t314) * t350) * qJD(1)) * t346;
t362 = t348 * t276 - t281 * t346;
t278 = t279 * t352 + t280 * t350;
t301 = -qJDD(1) * pkin(2) + t357;
t356 = -m(4) * t301 + mrSges(4,1) * t370 - t278 + (t348 ^ 2 * t354 + t381) * mrSges(4,3);
t327 = (Ifges(4,5) * t346 + Ifges(4,6) * t348) * qJD(1);
t283 = t388 - mrSges(6,1) * t320 + mrSges(6,2) * t321 + (t313 * t350 + t316 * t352) * t375;
t277 = qJDD(1) * t382 - t356;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t376 - mrSges(2,2) * t364 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t303 - mrSges(3,2) * t304 + t346 * (t327 * t374 + mrSges(4,2) * t301 - mrSges(4,3) * t294 + t352 * (mrSges(5,2) * t292 + mrSges(6,2) * t290 - mrSges(5,3) * t288 - mrSges(6,3) * t285 - qJ(5) * t282) - t350 * (-mrSges(5,1) * t292 + mrSges(5,3) * t289 - mrSges(6,1) * t290 + mrSges(6,3) * t287 - pkin(4) * t283 + qJ(5) * (-t318 * t367 + t380)) - pkin(6) * t278 + (-t363 - t350 * (-qJ(5) * t316 + t378)) * t333 + (t352 * t385 - t350 * (-qJ(5) * mrSges(6,2) + t384)) * t332 + t395 * t321 + t394 * t320 + (Ifges(4,1) * t346 + Ifges(4,4) * t348) * qJDD(1)) + t348 * (Ifges(4,2) * t370 - mrSges(4,1) * t301 + mrSges(4,3) * t295 - pkin(3) * t278 + (Ifges(4,4) * qJDD(1) + (-t327 - t392) * qJD(1)) * t346 - t389) - pkin(2) * t277 + qJ(3) * t362 + pkin(1) * (t347 * (m(3) * t304 - mrSges(3,1) * t354 - qJDD(1) * mrSges(3,2) + t362) + t349 * (m(3) * t303 - mrSges(3,2) * t354 + (mrSges(3,1) - t382) * qJDD(1) + t356)); m(3) * t345 + t276 * t346 + t281 * t348; t277; t375 * t392 + t389; t283;];
tauJ = t1;
