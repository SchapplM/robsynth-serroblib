% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRPP4
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRPP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP4_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP4_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:56
% EndTime: 2019-12-31 16:58:57
% DurationCPUTime: 0.72s
% Computational Cost: add. (775->163), mult. (1597->198), div. (0->0), fcn. (641->4), ass. (0->64)
t342 = Ifges(3,1) + Ifges(4,1) + Ifges(5,1);
t329 = Ifges(3,4) - Ifges(4,5) - Ifges(5,4);
t328 = Ifges(3,5) + Ifges(4,4) - Ifges(5,5);
t341 = Ifges(3,2) + Ifges(4,3) + Ifges(5,2);
t327 = Ifges(3,6) - Ifges(4,6) + Ifges(5,6);
t340 = Ifges(3,3) + Ifges(4,2) + Ifges(5,3);
t310 = cos(qJ(2));
t307 = t310 ^ 2;
t339 = 2 * qJD(3);
t313 = qJD(1) ^ 2;
t338 = pkin(5) * t313;
t337 = -mrSges(5,2) - mrSges(4,3);
t336 = qJ(3) * t310;
t309 = sin(qJ(1));
t311 = cos(qJ(1));
t319 = -g(1) * t311 - g(2) * t309;
t276 = -pkin(1) * t313 + qJDD(1) * pkin(5) + t319;
t308 = sin(qJ(2));
t258 = -t310 * g(3) - t308 * t276;
t332 = qJD(1) * t308;
t291 = -qJD(2) * mrSges(5,1) - mrSges(5,3) * t332;
t293 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t332;
t335 = -t291 - t293;
t331 = qJD(1) * t310;
t294 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t331;
t296 = mrSges(4,2) * t331 + qJD(2) * mrSges(4,3);
t334 = -t294 - t296;
t333 = t309 * g(1) - t311 * g(2);
t330 = qJD(1) * qJD(2);
t326 = -0.2e1 * qJD(1) * qJD(4);
t325 = -t340 * qJD(2) + (-t308 * t328 - t310 * t327) * qJD(1);
t324 = -t327 * qJD(2) + (-t308 * t329 - t341 * t310) * qJD(1);
t323 = t328 * qJD(2) + (t308 * t342 + t310 * t329) * qJD(1);
t322 = t310 * t330;
t285 = qJDD(1) * t310 - t308 * t330;
t290 = -qJD(2) * pkin(3) - qJ(4) * t332;
t284 = qJDD(1) * t308 + t322;
t318 = qJDD(1) * pkin(1) + t333;
t316 = -qJ(3) * t284 - t318;
t250 = qJDD(4) + (-qJ(4) * t307 + pkin(5)) * t313 + (pkin(2) + pkin(3)) * t285 + (qJD(2) * t336 + (-pkin(2) * qJD(2) + t290 + t339) * t308) * qJD(1) - t316;
t321 = m(5) * t250 + t285 * mrSges(5,1);
t259 = -g(3) * t308 + t310 * t276;
t280 = (-t310 * pkin(2) - t308 * qJ(3)) * qJD(1);
t312 = qJD(2) ^ 2;
t256 = -pkin(2) * t312 + qJDD(2) * qJ(3) + qJD(2) * t339 + t280 * t331 + t259;
t252 = -pkin(3) * t307 * t313 - qJ(4) * t285 + qJD(2) * t290 + t310 * t326 + t256;
t320 = m(5) * t252 + qJDD(2) * mrSges(5,2) - t285 * mrSges(5,3) + qJD(2) * t291;
t254 = -pkin(2) * t285 - t338 + (-0.2e1 * qJD(3) * t308 + (pkin(2) * t308 - t336) * qJD(2)) * qJD(1) + t316;
t317 = m(4) * t254 - t321;
t257 = -qJDD(2) * pkin(2) - qJ(3) * t312 + t280 * t332 + qJDD(3) - t258;
t253 = t308 * t326 + (-t284 + t322) * qJ(4) + (-t308 * t310 * t313 - qJDD(2)) * pkin(3) + t257;
t282 = (t310 * mrSges(5,1) + t308 * mrSges(5,2)) * qJD(1);
t249 = m(5) * t253 - qJDD(2) * mrSges(5,1) - t284 * mrSges(5,3) - qJD(2) * t294 - t282 * t332;
t281 = (-t310 * mrSges(4,1) - t308 * mrSges(4,3)) * qJD(1);
t315 = m(4) * t256 + qJDD(2) * mrSges(4,3) + qJD(2) * t293 + t281 * t331 + t320;
t314 = m(4) * t257 - qJDD(2) * mrSges(4,1) + t284 * mrSges(4,2) - qJD(2) * t296 + t249;
t295 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t331;
t292 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t332;
t283 = (-t310 * mrSges(3,1) + t308 * mrSges(3,2)) * qJD(1);
t275 = -t318 - t338;
t248 = mrSges(5,2) * t284 + (t291 * t308 + t294 * t310) * qJD(1) + t321;
t247 = t281 * t332 + t314;
t246 = -mrSges(4,1) * t285 + t337 * t284 + (t335 * t308 + t334 * t310) * qJD(1) + t317;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t333 - mrSges(2,2) * t319 + t308 * (mrSges(3,2) * t275 + mrSges(4,2) * t257 + mrSges(5,2) * t250 - mrSges(3,3) * t258 - mrSges(4,3) * t254 - mrSges(5,3) * t253 - qJ(3) * t246 - qJ(4) * t249 + t324 * qJD(2) + t328 * qJDD(2) + t342 * t284 + t329 * t285 - t325 * t331) + t310 * (-mrSges(3,1) * t275 + mrSges(3,3) * t259 - mrSges(4,1) * t254 + mrSges(4,2) * t256 + mrSges(5,1) * t250 - mrSges(5,3) * t252 + pkin(3) * t248 - qJ(4) * t320 - pkin(2) * t246 + t341 * t285 + t329 * t284 + t327 * qJDD(2) + t323 * qJD(2) + (qJ(4) * t282 * t310 + t325 * t308) * qJD(1)) + pkin(1) * (-m(3) * t275 + (mrSges(3,1) + mrSges(4,1)) * t285 + (-mrSges(3,2) - t337) * t284 + ((t295 - t334) * t310 + (-t292 - t335) * t308) * qJD(1) - t317) + pkin(5) * (t310 * (m(3) * t259 - qJDD(2) * mrSges(3,2) - qJD(2) * t292 + t315 + (mrSges(3,3) + mrSges(4,2)) * t285) - t308 * (m(3) * t258 + qJDD(2) * mrSges(3,1) - t284 * mrSges(3,3) + qJD(2) * t295 - t314) + ((-t282 + t283) * t307 + (t281 + t283) * t308 ^ 2) * qJD(1)); mrSges(3,1) * t258 - mrSges(3,2) * t259 - mrSges(4,1) * t257 + mrSges(4,3) * t256 - mrSges(5,1) * t253 + mrSges(5,2) * t252 - pkin(3) * t249 - pkin(2) * t247 + qJ(3) * t315 + (qJ(3) * mrSges(4,2) + t327) * t285 + t328 * t284 + t340 * qJDD(2) + (-t324 * t308 + (-qJ(3) * t282 - t323) * t310) * qJD(1); t247; t248;];
tauJ = t1;
