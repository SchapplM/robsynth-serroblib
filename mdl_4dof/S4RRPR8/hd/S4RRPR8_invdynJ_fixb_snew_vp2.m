% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRPR8
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRPR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR8_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR8_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:53
% EndTime: 2019-12-31 17:07:54
% DurationCPUTime: 0.87s
% Computational Cost: add. (1727->190), mult. (3588->236), div. (0->0), fcn. (1816->6), ass. (0->78)
t340 = Ifges(3,1) + Ifges(4,1);
t334 = Ifges(3,4) - Ifges(4,5);
t333 = Ifges(3,5) + Ifges(4,4);
t339 = Ifges(3,2) + Ifges(4,3);
t332 = Ifges(3,6) - Ifges(4,6);
t338 = Ifges(3,3) + Ifges(4,2);
t337 = 2 * qJD(3);
t313 = qJD(1) ^ 2;
t336 = pkin(5) * t313;
t335 = mrSges(3,3) + mrSges(4,2);
t310 = cos(qJ(2));
t331 = qJ(3) * t310;
t308 = sin(qJ(1));
t311 = cos(qJ(1));
t321 = -g(1) * t311 - g(2) * t308;
t282 = -pkin(1) * t313 + qJDD(1) * pkin(5) + t321;
t307 = sin(qJ(2));
t264 = -t310 * g(3) - t307 * t282;
t330 = t338 * qJD(2) + (t307 * t333 + t310 * t332) * qJD(1);
t329 = -t332 * qJD(2) + (-t307 * t334 - t339 * t310) * qJD(1);
t328 = t333 * qJD(2) + (t307 * t340 + t310 * t334) * qJD(1);
t327 = t308 * g(1) - t311 * g(2);
t326 = qJD(1) * t307;
t325 = qJD(1) * t310;
t324 = qJD(1) * qJD(2);
t323 = t310 * t324;
t265 = -g(3) * t307 + t310 * t282;
t283 = (-pkin(2) * t310 - qJ(3) * t307) * qJD(1);
t312 = qJD(2) ^ 2;
t254 = -pkin(2) * t312 + qJDD(2) * qJ(3) + qJD(2) * t337 + t283 * t325 + t265;
t287 = qJDD(1) * t310 - t307 * t324;
t294 = -qJD(2) * pkin(3) - pkin(6) * t326;
t305 = t310 ^ 2;
t250 = -pkin(3) * t305 * t313 - pkin(6) * t287 + qJD(2) * t294 + t254;
t256 = -qJDD(2) * pkin(2) - qJ(3) * t312 + t283 * t326 + qJDD(3) - t264;
t286 = qJDD(1) * t307 + t323;
t251 = (-t286 + t323) * pkin(6) + (-t307 * t310 * t313 - qJDD(2)) * pkin(3) + t256;
t306 = sin(qJ(4));
t309 = cos(qJ(4));
t247 = -t250 * t306 + t251 * t309;
t279 = (-t306 * t307 - t309 * t310) * qJD(1);
t258 = qJD(4) * t279 + t286 * t309 - t287 * t306;
t280 = (-t306 * t310 + t307 * t309) * qJD(1);
t263 = -mrSges(5,1) * t279 + mrSges(5,2) * t280;
t301 = -qJD(2) + qJD(4);
t266 = -mrSges(5,2) * t301 + mrSges(5,3) * t279;
t300 = -qJDD(2) + qJDD(4);
t245 = m(5) * t247 + mrSges(5,1) * t300 - mrSges(5,3) * t258 - t263 * t280 + t266 * t301;
t248 = t250 * t309 + t251 * t306;
t257 = -qJD(4) * t280 - t286 * t306 - t287 * t309;
t267 = mrSges(5,1) * t301 - mrSges(5,3) * t280;
t246 = m(5) * t248 - mrSges(5,2) * t300 + mrSges(5,3) * t257 + t263 * t279 - t267 * t301;
t322 = -t306 * t245 + t309 * t246;
t320 = qJDD(1) * pkin(1) + t327;
t240 = t309 * t245 + t306 * t246;
t319 = -qJ(3) * t286 - t320;
t284 = (-mrSges(4,1) * t310 - mrSges(4,3) * t307) * qJD(1);
t291 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t326;
t318 = m(4) * t254 + qJDD(2) * mrSges(4,3) + qJD(2) * t291 + t284 * t325 + t322;
t249 = (-pkin(6) * t305 + pkin(5)) * t313 + (pkin(2) + pkin(3)) * t287 + (qJD(2) * t331 + (-pkin(2) * qJD(2) + t294 + t337) * t307) * qJD(1) - t319;
t317 = m(5) * t249 - t257 * mrSges(5,1) + mrSges(5,2) * t258 - t279 * t266 + t267 * t280;
t293 = mrSges(4,2) * t325 + qJD(2) * mrSges(4,3);
t316 = m(4) * t256 - qJDD(2) * mrSges(4,1) - qJD(2) * t293 + t240;
t252 = -pkin(2) * t287 - t336 + (-0.2e1 * qJD(3) * t307 + (pkin(2) * t307 - t331) * qJD(2)) * qJD(1) + t319;
t315 = m(4) * t252 - t317;
t260 = Ifges(5,4) * t280 + Ifges(5,2) * t279 + Ifges(5,6) * t301;
t261 = Ifges(5,1) * t280 + Ifges(5,4) * t279 + Ifges(5,5) * t301;
t314 = mrSges(5,1) * t247 - mrSges(5,2) * t248 + Ifges(5,5) * t258 + Ifges(5,6) * t257 + Ifges(5,3) * t300 + t280 * t260 - t279 * t261;
t292 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t325;
t290 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t326;
t285 = (-mrSges(3,1) * t310 + mrSges(3,2) * t307) * qJD(1);
t281 = -t320 - t336;
t259 = Ifges(5,5) * t280 + Ifges(5,6) * t279 + Ifges(5,3) * t301;
t243 = -mrSges(4,1) * t287 - mrSges(4,3) * t286 + (-t291 * t307 - t293 * t310) * qJD(1) + t315;
t242 = mrSges(5,2) * t249 - mrSges(5,3) * t247 + Ifges(5,1) * t258 + Ifges(5,4) * t257 + Ifges(5,5) * t300 + t259 * t279 - t260 * t301;
t241 = -mrSges(5,1) * t249 + mrSges(5,3) * t248 + Ifges(5,4) * t258 + Ifges(5,2) * t257 + Ifges(5,6) * t300 - t259 * t280 + t261 * t301;
t239 = t286 * mrSges(4,2) + t284 * t326 + t316;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t327 - mrSges(2,2) * t321 + t307 * (mrSges(3,2) * t281 + mrSges(4,2) * t256 - mrSges(3,3) * t264 - mrSges(4,3) * t252 - pkin(6) * t240 - qJ(3) * t243 + t329 * qJD(2) + t333 * qJDD(2) - t306 * t241 + t309 * t242 + t340 * t286 + t334 * t287 + t330 * t325) + t310 * (-mrSges(3,1) * t281 - mrSges(4,1) * t252 + mrSges(4,2) * t254 + mrSges(3,3) * t265 - pkin(2) * t243 + pkin(3) * t317 - pkin(6) * t322 + t328 * qJD(2) + t332 * qJDD(2) - t309 * t241 - t306 * t242 + t334 * t286 + t339 * t287 - t330 * t326) + pkin(1) * (-m(3) * t281 + (mrSges(3,1) + mrSges(4,1)) * t287 + (-mrSges(3,2) + mrSges(4,3)) * t286 + ((t292 + t293) * t310 + (-t290 + t291) * t307) * qJD(1) - t315) + pkin(5) * (t310 * (m(3) * t265 - qJDD(2) * mrSges(3,2) - qJD(2) * t290 + t285 * t325 + t287 * t335 + t318) + (-m(3) * t264 - qJDD(2) * mrSges(3,1) - qJD(2) * t292 + t335 * t286 + (t284 + t285) * t326 + t316) * t307); -pkin(3) * t240 - pkin(2) * t239 - t314 + mrSges(4,3) * t254 - mrSges(4,1) * t256 + mrSges(3,1) * t264 - mrSges(3,2) * t265 + t338 * qJDD(2) + t333 * t286 + qJ(3) * t318 + (mrSges(4,2) * qJ(3) + t332) * t287 + (-t329 * t307 - t328 * t310) * qJD(1); t239; t314;];
tauJ = t1;
