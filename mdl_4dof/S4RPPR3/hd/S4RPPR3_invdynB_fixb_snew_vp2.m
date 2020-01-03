% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR3_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR3_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR3_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:53
% EndTime: 2019-12-31 16:37:54
% DurationCPUTime: 1.01s
% Computational Cost: add. (8578->174), mult. (17882->222), div. (0->0), fcn. (10403->8), ass. (0->79)
t356 = qJD(1) ^ 2;
t350 = cos(pkin(7));
t378 = pkin(3) * t350;
t348 = sin(pkin(7));
t377 = mrSges(4,2) * t348;
t346 = t350 ^ 2;
t376 = t346 * t356;
t353 = sin(qJ(1));
t355 = cos(qJ(1));
t334 = t353 * g(1) - t355 * g(2);
t332 = qJDD(1) * pkin(1) + t334;
t335 = -t355 * g(1) - t353 * g(2);
t333 = -t356 * pkin(1) + t335;
t349 = sin(pkin(6));
t351 = cos(pkin(6));
t322 = t349 * t332 + t351 * t333;
t315 = -t356 * pkin(2) + qJDD(1) * qJ(3) + t322;
t347 = -g(3) + qJDD(2);
t372 = qJD(1) * qJD(3);
t374 = t350 * t347 - 0.2e1 * t348 * t372;
t305 = (-pkin(5) * qJDD(1) + t356 * t378 - t315) * t348 + t374;
t309 = t348 * t347 + (t315 + 0.2e1 * t372) * t350;
t371 = qJDD(1) * t350;
t306 = -pkin(3) * t376 + pkin(5) * t371 + t309;
t352 = sin(qJ(4));
t354 = cos(qJ(4));
t303 = t354 * t305 - t352 * t306;
t360 = -t348 * t352 + t350 * t354;
t325 = t360 * qJD(1);
t361 = t348 * t354 + t350 * t352;
t326 = t361 * qJD(1);
t317 = -t325 * mrSges(5,1) + t326 * mrSges(5,2);
t320 = t325 * qJD(4) + t361 * qJDD(1);
t323 = -qJD(4) * mrSges(5,2) + t325 * mrSges(5,3);
t301 = m(5) * t303 + qJDD(4) * mrSges(5,1) - t320 * mrSges(5,3) + qJD(4) * t323 - t326 * t317;
t304 = t352 * t305 + t354 * t306;
t319 = -t326 * qJD(4) + t360 * qJDD(1);
t324 = qJD(4) * mrSges(5,1) - t326 * mrSges(5,3);
t302 = m(5) * t304 - qJDD(4) * mrSges(5,2) + t319 * mrSges(5,3) - qJD(4) * t324 + t325 * t317;
t293 = t354 * t301 + t352 * t302;
t308 = -t348 * t315 + t374;
t359 = mrSges(4,3) * qJDD(1) + t356 * (-mrSges(4,1) * t350 + t377);
t291 = m(4) * t308 - t359 * t348 + t293;
t366 = -t352 * t301 + t354 * t302;
t292 = m(4) * t309 + t359 * t350 + t366;
t367 = -t348 * t291 + t350 * t292;
t286 = m(3) * t322 - t356 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t367;
t321 = t351 * t332 - t349 * t333;
t362 = qJDD(3) - t321;
t311 = -qJDD(1) * pkin(2) - t356 * qJ(3) + t362;
t345 = t348 ^ 2;
t307 = (-pkin(2) - t378) * qJDD(1) + (-qJ(3) + (-t345 - t346) * pkin(5)) * t356 + t362;
t358 = m(5) * t307 - t319 * mrSges(5,1) + t320 * mrSges(5,2) - t325 * t323 + t326 * t324;
t357 = -m(4) * t311 + mrSges(4,1) * t371 - t358 + (t345 * t356 + t376) * mrSges(4,3);
t297 = m(3) * t321 - t356 * mrSges(3,2) + (mrSges(3,1) - t377) * qJDD(1) + t357;
t282 = t349 * t286 + t351 * t297;
t280 = m(2) * t334 + qJDD(1) * mrSges(2,1) - t356 * mrSges(2,2) + t282;
t368 = t351 * t286 - t349 * t297;
t281 = m(2) * t335 - t356 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t368;
t375 = t355 * t280 + t353 * t281;
t287 = t350 * t291 + t348 * t292;
t363 = Ifges(4,5) * t348 + Ifges(4,6) * t350;
t373 = t356 * t363;
t370 = m(3) * t347 + t287;
t369 = -t353 * t280 + t355 * t281;
t365 = Ifges(4,1) * t348 + Ifges(4,4) * t350;
t364 = Ifges(4,4) * t348 + Ifges(4,2) * t350;
t314 = Ifges(5,1) * t326 + Ifges(5,4) * t325 + Ifges(5,5) * qJD(4);
t313 = Ifges(5,4) * t326 + Ifges(5,2) * t325 + Ifges(5,6) * qJD(4);
t312 = Ifges(5,5) * t326 + Ifges(5,6) * t325 + Ifges(5,3) * qJD(4);
t295 = mrSges(5,2) * t307 - mrSges(5,3) * t303 + Ifges(5,1) * t320 + Ifges(5,4) * t319 + Ifges(5,5) * qJDD(4) - qJD(4) * t313 + t325 * t312;
t294 = -mrSges(5,1) * t307 + mrSges(5,3) * t304 + Ifges(5,4) * t320 + Ifges(5,2) * t319 + Ifges(5,6) * qJDD(4) + qJD(4) * t314 - t326 * t312;
t283 = mrSges(4,2) * t311 - mrSges(4,3) * t308 - pkin(5) * t293 + t365 * qJDD(1) - t352 * t294 + t354 * t295 + t350 * t373;
t276 = -mrSges(4,1) * t311 + mrSges(4,3) * t309 - pkin(3) * t358 + pkin(5) * t366 + t364 * qJDD(1) + t354 * t294 + t352 * t295 - t348 * t373;
t275 = -mrSges(3,1) * t347 - mrSges(4,1) * t308 - mrSges(5,1) * t303 + mrSges(4,2) * t309 + mrSges(5,2) * t304 + mrSges(3,3) * t322 - Ifges(5,5) * t320 - Ifges(5,6) * t319 - Ifges(5,3) * qJDD(4) - pkin(2) * t287 - pkin(3) * t293 - t326 * t313 + t325 * t314 + (Ifges(3,6) - t363) * qJDD(1) + (-t348 * t364 + t350 * t365 + Ifges(3,5)) * t356;
t274 = mrSges(3,2) * t347 - mrSges(3,3) * t321 + Ifges(3,5) * qJDD(1) - t356 * Ifges(3,6) - qJ(3) * t287 - t348 * t276 + t350 * t283;
t273 = -mrSges(2,2) * g(3) - mrSges(2,3) * t334 + Ifges(2,5) * qJDD(1) - t356 * Ifges(2,6) - qJ(2) * t282 + t351 * t274 - t349 * t275;
t272 = mrSges(2,1) * g(3) + mrSges(2,3) * t335 + t356 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t370 + qJ(2) * t368 + t349 * t274 + t351 * t275;
t1 = [-m(1) * g(1) + t369; -m(1) * g(2) + t375; (-m(1) - m(2)) * g(3) + t370; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t375 - t353 * t272 + t355 * t273; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t369 + t355 * t272 + t353 * t273; pkin(1) * t282 + mrSges(2,1) * t334 - mrSges(2,2) * t335 + t348 * t283 + t350 * t276 + pkin(2) * t357 + qJ(3) * t367 + mrSges(3,1) * t321 - mrSges(3,2) * t322 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * t377 + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
