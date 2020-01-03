% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRRP4
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP4_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP4_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:16
% EndTime: 2019-12-31 17:15:17
% DurationCPUTime: 0.87s
% Computational Cost: add. (2530->182), mult. (5362->223), div. (0->0), fcn. (3205->6), ass. (0->73)
t331 = Ifges(4,4) + Ifges(5,4);
t340 = Ifges(4,2) + Ifges(5,2);
t337 = Ifges(4,6) + Ifges(5,6);
t339 = Ifges(4,1) + Ifges(5,1);
t338 = Ifges(4,5) + Ifges(5,5);
t336 = Ifges(4,3) + Ifges(5,3);
t308 = sin(qJ(3));
t309 = sin(qJ(2));
t311 = cos(qJ(3));
t312 = cos(qJ(2));
t290 = (-t309 * t308 + t312 * t311) * qJD(1);
t291 = (t312 * t308 + t309 * t311) * qJD(1);
t306 = qJD(2) + qJD(3);
t335 = t340 * t290 + t331 * t291 + t337 * t306;
t324 = qJD(1) * qJD(2);
t296 = t309 * qJDD(1) + t312 * t324;
t297 = t312 * qJDD(1) - t309 * t324;
t265 = -t291 * qJD(3) - t308 * t296 + t311 * t297;
t266 = t290 * qJD(3) + t311 * t296 + t308 * t297;
t326 = qJD(1) * t309;
t300 = qJD(2) * pkin(2) - pkin(6) * t326;
t307 = t312 ^ 2;
t314 = qJD(1) ^ 2;
t310 = sin(qJ(1));
t313 = cos(qJ(1));
t320 = t310 * g(1) - t313 * g(2);
t317 = -qJDD(1) * pkin(1) - t320;
t267 = -t297 * pkin(2) + t300 * t326 + (-pkin(6) * t307 - pkin(5)) * t314 + t317;
t280 = -t306 * mrSges(5,2) + t290 * mrSges(5,3);
t281 = -t306 * mrSges(4,2) + t290 * mrSges(4,3);
t284 = t306 * mrSges(4,1) - t291 * mrSges(4,3);
t282 = t306 * pkin(3) - t291 * qJ(4);
t286 = t290 ^ 2;
t247 = -t265 * pkin(3) - t286 * qJ(4) + t291 * t282 + qJDD(4) + t267;
t283 = t306 * mrSges(5,1) - t291 * mrSges(5,3);
t321 = m(5) * t247 + t266 * mrSges(5,2) + t291 * t283;
t334 = -m(4) * t267 - t266 * mrSges(4,2) + (mrSges(4,1) + mrSges(5,1)) * t265 - t291 * t284 + (t280 + t281) * t290 - t321;
t333 = pkin(2) * t314;
t318 = -t313 * g(1) - t310 * g(2);
t293 = -t314 * pkin(1) + qJDD(1) * pkin(5) + t318;
t330 = t309 * t293;
t256 = qJDD(2) * pkin(2) - t296 * pkin(6) - t330 + (pkin(6) * t324 + t309 * t333 - g(3)) * t312;
t279 = -t309 * g(3) + t312 * t293;
t257 = t297 * pkin(6) - qJD(2) * t300 - t307 * t333 + t279;
t249 = t311 * t256 - t308 * t257;
t275 = -t290 * mrSges(5,1) + t291 * mrSges(5,2);
t276 = -t290 * mrSges(4,1) + t291 * mrSges(4,2);
t305 = qJDD(2) + qJDD(3);
t243 = -0.2e1 * qJD(4) * t291 + (t290 * t306 - t266) * qJ(4) + (t290 * t291 + t305) * pkin(3) + t249;
t323 = m(5) * t243 + t305 * mrSges(5,1) + t306 * t280;
t234 = m(4) * t249 + t305 * mrSges(4,1) + t306 * t281 + (-t275 - t276) * t291 + (-mrSges(4,3) - mrSges(5,3)) * t266 + t323;
t250 = t308 * t256 + t311 * t257;
t245 = -t286 * pkin(3) + t265 * qJ(4) + 0.2e1 * qJD(4) * t290 - t306 * t282 + t250;
t322 = m(5) * t245 + t265 * mrSges(5,3) + t290 * t275;
t237 = m(4) * t250 + t265 * mrSges(4,3) + t290 * t276 + (-t283 - t284) * t306 + (-mrSges(4,2) - mrSges(5,2)) * t305 + t322;
t232 = t311 * t234 + t308 * t237;
t329 = -t337 * t290 - t338 * t291 - t336 * t306;
t328 = -t331 * t290 - t339 * t291 - t338 * t306;
t325 = qJD(1) * t312;
t319 = -t308 * t234 + t311 * t237;
t239 = -t266 * mrSges(5,3) - t291 * t275 + t323;
t315 = mrSges(4,1) * t249 + mrSges(5,1) * t243 - mrSges(4,2) * t250 - mrSges(5,2) * t245 + pkin(3) * t239 + t337 * t265 + t338 * t266 + t290 * t328 + t335 * t291 + t336 * t305;
t299 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t325;
t298 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t326;
t295 = (-t312 * mrSges(3,1) + t309 * mrSges(3,2)) * qJD(1);
t292 = -t314 * pkin(5) + t317;
t289 = Ifges(3,5) * qJD(2) + (t309 * Ifges(3,1) + t312 * Ifges(3,4)) * qJD(1);
t288 = Ifges(3,6) * qJD(2) + (t309 * Ifges(3,4) + t312 * Ifges(3,2)) * qJD(1);
t278 = -t312 * g(3) - t330;
t240 = -t265 * mrSges(5,1) - t290 * t280 + t321;
t231 = mrSges(4,2) * t267 + mrSges(5,2) * t247 - mrSges(4,3) * t249 - mrSges(5,3) * t243 - qJ(4) * t239 + t331 * t265 + t339 * t266 - t329 * t290 + t338 * t305 - t335 * t306;
t230 = -mrSges(4,1) * t267 + mrSges(4,3) * t250 - mrSges(5,1) * t247 + mrSges(5,3) * t245 - pkin(3) * t240 + qJ(4) * t322 + (-qJ(4) * t283 - t328) * t306 + (-qJ(4) * mrSges(5,2) + t337) * t305 + t329 * t291 + t331 * t266 + t340 * t265;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t320 - mrSges(2,2) * t318 + t309 * (mrSges(3,2) * t292 - mrSges(3,3) * t278 + Ifges(3,1) * t296 + Ifges(3,4) * t297 + Ifges(3,5) * qJDD(2) - pkin(6) * t232 - qJD(2) * t288 - t308 * t230 + t311 * t231) + t312 * (-mrSges(3,1) * t292 + mrSges(3,3) * t279 + Ifges(3,4) * t296 + Ifges(3,2) * t297 + Ifges(3,6) * qJDD(2) + pkin(2) * t334 + pkin(6) * t319 + qJD(2) * t289 + t311 * t230 + t308 * t231) + pkin(1) * (-m(3) * t292 + t297 * mrSges(3,1) - t296 * mrSges(3,2) + (-t298 * t309 + t299 * t312) * qJD(1) + t334) + pkin(5) * (t312 * (m(3) * t279 - qJDD(2) * mrSges(3,2) + t297 * mrSges(3,3) - qJD(2) * t298 + t295 * t325 + t319) - t309 * (m(3) * t278 + qJDD(2) * mrSges(3,1) - t296 * mrSges(3,3) + qJD(2) * t299 - t295 * t326 + t232)); t315 + (t309 * t288 - t312 * t289) * qJD(1) + mrSges(3,1) * t278 - mrSges(3,2) * t279 + pkin(2) * t232 + Ifges(3,5) * t296 + Ifges(3,6) * t297 + Ifges(3,3) * qJDD(2); t315; t240;];
tauJ = t1;
