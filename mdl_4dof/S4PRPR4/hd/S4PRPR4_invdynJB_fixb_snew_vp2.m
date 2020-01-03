% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRPR4
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR4_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR4_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR4_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:58
% EndTime: 2019-12-31 16:21:59
% DurationCPUTime: 0.44s
% Computational Cost: add. (2805->135), mult. (4788->168), div. (0->0), fcn. (2384->6), ass. (0->63)
t354 = -pkin(2) - pkin(5);
t353 = mrSges(3,1) - mrSges(4,2);
t352 = -Ifges(4,4) + Ifges(3,5);
t351 = Ifges(4,5) - Ifges(3,6);
t328 = sin(pkin(6));
t329 = cos(pkin(6));
t314 = g(1) * t328 - g(2) * t329;
t315 = -g(1) * t329 - g(2) * t328;
t331 = sin(qJ(2));
t333 = cos(qJ(2));
t296 = t314 * t333 - t331 * t315;
t334 = qJD(2) ^ 2;
t340 = -qJ(3) * t334 + qJDD(3) - t296;
t291 = t354 * qJDD(2) + t340;
t325 = -g(3) + qJDD(1);
t330 = sin(qJ(4));
t332 = cos(qJ(4));
t287 = t291 * t332 - t325 * t330;
t311 = (mrSges(5,1) * t330 + mrSges(5,2) * t332) * qJD(2);
t347 = qJD(2) * qJD(4);
t313 = qJDD(2) * t332 - t330 * t347;
t349 = qJD(2) * t330;
t316 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t349;
t348 = qJD(2) * t332;
t284 = m(5) * t287 + qJDD(4) * mrSges(5,1) - t313 * mrSges(5,3) + qJD(4) * t316 - t311 * t348;
t288 = t291 * t330 + t325 * t332;
t312 = -qJDD(2) * t330 - t332 * t347;
t317 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t348;
t285 = m(5) * t288 - qJDD(4) * mrSges(5,2) + t312 * mrSges(5,3) - qJD(4) * t317 - t311 * t349;
t275 = t284 * t332 + t285 * t330;
t294 = -qJDD(2) * pkin(2) + t340;
t338 = -m(4) * t294 + t334 * mrSges(4,3) - t275;
t271 = m(3) * t296 - mrSges(3,2) * t334 + t353 * qJDD(2) + t338;
t297 = t331 * t314 + t333 * t315;
t339 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t297;
t292 = pkin(2) * t334 - t339;
t290 = t354 * t334 + t339;
t342 = -m(5) * t290 + mrSges(5,1) * t312 - t313 * mrSges(5,2) - t316 * t349 - t317 * t348;
t336 = -m(4) * t292 + t334 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t342;
t278 = m(3) * t297 - mrSges(3,1) * t334 - qJDD(2) * mrSges(3,2) + t336;
t269 = t333 * t271 + t331 * t278;
t267 = m(2) * t314 + t269;
t345 = -t331 * t271 + t333 * t278;
t268 = m(2) * t315 + t345;
t350 = t329 * t267 + t328 * t268;
t346 = -t267 * t328 + t329 * t268;
t344 = -t284 * t330 + t332 * t285;
t274 = m(4) * t325 + t344;
t343 = m(3) * t325 + t274;
t341 = m(2) * t325 + t343;
t301 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t332 - Ifges(5,2) * t330) * qJD(2);
t302 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t332 - Ifges(5,4) * t330) * qJD(2);
t337 = mrSges(5,1) * t287 - mrSges(5,2) * t288 + Ifges(5,5) * t313 + Ifges(5,6) * t312 + Ifges(5,3) * qJDD(4) + t301 * t348 + t302 * t349;
t273 = qJDD(2) * mrSges(4,2) - t338;
t300 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t332 - Ifges(5,6) * t330) * qJD(2);
t280 = -mrSges(5,1) * t290 + mrSges(5,3) * t288 + Ifges(5,4) * t313 + Ifges(5,2) * t312 + Ifges(5,6) * qJDD(4) + qJD(4) * t302 - t300 * t348;
t281 = mrSges(5,2) * t290 - mrSges(5,3) * t287 + Ifges(5,1) * t313 + Ifges(5,4) * t312 + Ifges(5,5) * qJDD(4) - qJD(4) * t301 - t300 * t349;
t335 = mrSges(3,1) * t296 - mrSges(3,2) * t297 + mrSges(4,2) * t294 - mrSges(4,3) * t292 - pkin(2) * t273 - pkin(5) * t275 + qJ(3) * t336 - t280 * t330 + t332 * t281 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2);
t263 = mrSges(4,1) * t294 - mrSges(3,3) * t296 + pkin(3) * t275 - qJ(3) * t274 + t351 * t334 + (mrSges(3,2) - mrSges(4,3)) * t325 + t352 * qJDD(2) + t337;
t262 = -mrSges(4,1) * t292 + mrSges(3,3) * t297 - pkin(2) * t274 - pkin(3) * t342 - pkin(5) * t344 - t351 * qJDD(2) - t332 * t280 - t330 * t281 - t353 * t325 + t352 * t334;
t261 = mrSges(2,2) * t325 - mrSges(2,3) * t314 - pkin(4) * t269 - t262 * t331 + t263 * t333;
t260 = -mrSges(2,1) * t325 + mrSges(2,3) * t315 - pkin(1) * t343 + pkin(4) * t345 + t333 * t262 + t331 * t263;
t1 = [-m(1) * g(1) + t346; -m(1) * g(2) + t350; -m(1) * g(3) + t341; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t350 - t328 * t260 + t329 * t261; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t346 + t329 * t260 + t328 * t261; -mrSges(1,1) * g(2) + mrSges(2,1) * t314 + mrSges(1,2) * g(1) - mrSges(2,2) * t315 + pkin(1) * t269 + t335; t341; t335; t273; t337;];
tauJB = t1;
