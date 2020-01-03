% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRPP5
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRPP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP5_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP5_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:14
% EndTime: 2019-12-31 17:00:15
% DurationCPUTime: 0.74s
% Computational Cost: add. (773->166), mult. (1585->191), div. (0->0), fcn. (630->4), ass. (0->69)
t334 = Ifges(3,1) + Ifges(4,2) + Ifges(5,3);
t317 = Ifges(3,4) + Ifges(4,6) - Ifges(5,6);
t316 = Ifges(3,5) - Ifges(4,4) + Ifges(5,5);
t333 = Ifges(3,2) + Ifges(4,3) + Ifges(5,2);
t315 = Ifges(3,6) - Ifges(4,5) - Ifges(5,4);
t332 = Ifges(3,3) + Ifges(4,1) + Ifges(5,1);
t295 = sin(qJ(2));
t318 = qJD(1) * qJD(2);
t310 = t295 * t318;
t320 = qJD(1) * t295;
t330 = -2 * qJD(3);
t331 = pkin(2) * t310 + t320 * t330;
t329 = -2 * qJD(4);
t328 = mrSges(3,1) - mrSges(4,2);
t327 = mrSges(5,1) + mrSges(3,3);
t326 = -mrSges(5,2) - mrSges(4,3);
t297 = cos(qJ(2));
t311 = t297 * t318;
t274 = qJDD(1) * t295 + t311;
t325 = t274 * mrSges(5,1);
t275 = qJDD(1) * t297 - t310;
t281 = pkin(3) * t320 - qJD(2) * qJ(4);
t294 = t297 ^ 2;
t300 = qJD(1) ^ 2;
t296 = sin(qJ(1));
t298 = cos(qJ(1));
t309 = g(1) * t296 - t298 * g(2);
t304 = -qJDD(1) * pkin(1) - t309;
t243 = -qJ(3) * t274 + (-pkin(3) * t294 - pkin(5)) * t300 + (-pkin(2) - qJ(4)) * t275 + (-t281 * t295 + (-qJ(3) * qJD(2) + t329) * t297) * qJD(1) + t304 + t331;
t324 = m(5) * t243 - t275 * mrSges(5,3);
t305 = -g(1) * t298 - g(2) * t296;
t267 = -pkin(1) * t300 + qJDD(1) * pkin(5) + t305;
t250 = -t297 * g(3) - t295 * t267;
t271 = (mrSges(4,2) * t297 - mrSges(4,3) * t295) * qJD(1);
t273 = (-mrSges(5,2) * t295 - mrSges(5,3) * t297) * qJD(1);
t323 = t271 + t273;
t282 = mrSges(5,1) * t320 - qJD(2) * mrSges(5,3);
t285 = mrSges(4,1) * t320 + qJD(2) * mrSges(4,2);
t322 = -t282 - t285;
t319 = qJD(1) * t297;
t283 = -mrSges(4,1) * t319 - qJD(2) * mrSges(4,3);
t284 = mrSges(5,1) * t319 + qJD(2) * mrSges(5,2);
t321 = t283 - t284;
t314 = t332 * qJD(2) + (t295 * t316 + t297 * t315) * qJD(1);
t313 = t315 * qJD(2) + (t295 * t317 + t297 * t333) * qJD(1);
t312 = t316 * qJD(2) + (t295 * t334 + t297 * t317) * qJD(1);
t251 = -g(3) * t295 + t297 * t267;
t270 = (-pkin(2) * t297 - qJ(3) * t295) * qJD(1);
t299 = qJD(2) ^ 2;
t301 = -pkin(2) * t299 + qJDD(2) * qJ(3) + t270 * t319 + t251;
t246 = -qJ(4) * t294 * t300 + pkin(3) * t275 + qJDD(4) + ((2 * qJD(3)) + t281) * qJD(2) + t301;
t308 = m(5) * t246 + qJDD(2) * mrSges(5,2) + qJD(2) * t282 + t273 * t319;
t266 = -pkin(5) * t300 + t304;
t247 = -pkin(2) * t275 + (-t274 - t311) * qJ(3) + t266 + t331;
t307 = m(4) * t247 + t324;
t249 = -qJDD(2) * pkin(2) - qJ(3) * t299 + t270 * t320 + qJDD(3) - t250;
t245 = qJD(2) * t329 + (-t295 * t297 * t300 - qJDD(2)) * qJ(4) + (t274 - t311) * pkin(3) + t249;
t306 = m(5) * t245 - qJDD(2) * mrSges(5,3) - qJD(2) * t284;
t303 = m(4) * t249 + t274 * mrSges(4,1) + t306;
t248 = qJD(2) * t330 - t301;
t302 = -m(4) * t248 + qJDD(2) * mrSges(4,3) + qJD(2) * t285 + t271 * t319 + t308;
t280 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t319;
t279 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t320;
t272 = (-mrSges(3,1) * t297 + mrSges(3,2) * t295) * qJD(1);
t241 = mrSges(5,1) * t275 + t308;
t240 = t273 * t320 + t306 + t325;
t239 = qJDD(2) * mrSges(4,2) + qJD(2) * t283 + t323 * t320 + t303 + t325;
t238 = mrSges(4,2) * t275 + t326 * t274 + (t322 * t295 + t321 * t297) * qJD(1) + t307;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t309 - mrSges(2,2) * t305 + t295 * (mrSges(4,1) * t249 + mrSges(5,1) * t245 + mrSges(3,2) * t266 - mrSges(5,2) * t243 - mrSges(3,3) * t250 - mrSges(4,3) * t247 + pkin(3) * t240 - qJ(3) * t238 - t313 * qJD(2) + t316 * qJDD(2) + t334 * t274 + t317 * t275 + t314 * t319) + t297 * (-mrSges(3,1) * t266 + mrSges(3,3) * t251 - mrSges(4,1) * t248 + mrSges(4,2) * t247 + mrSges(5,1) * t246 - mrSges(5,3) * t243 + pkin(3) * t241 - qJ(4) * t324 - pkin(2) * t238 + t333 * t275 + (mrSges(5,2) * qJ(4) + t317) * t274 + t315 * qJDD(2) + t312 * qJD(2) + (qJ(4) * t284 * t297 + (qJ(4) * t282 - t314) * t295) * qJD(1)) + pkin(1) * (-m(3) * t266 + t328 * t275 + (-mrSges(3,2) - t326) * t274 + ((t280 - t321) * t297 + (-t279 - t322) * t295) * qJD(1) - t307) + pkin(5) * (t297 * (t272 * t319 + m(3) * t251 - qJDD(2) * mrSges(3,2) - qJD(2) * t279 + (mrSges(4,1) + t327) * t275 + t302) + (-m(3) * t250 + t327 * t274 - t328 * qJDD(2) + (-t280 + t283) * qJD(2) + (t272 + t323) * t320 + t303) * t295); mrSges(3,1) * t250 - mrSges(3,2) * t251 + mrSges(4,2) * t249 - mrSges(4,3) * t248 + mrSges(5,2) * t246 - mrSges(5,3) * t245 - qJ(4) * t240 - pkin(2) * t239 + qJ(3) * t302 + t316 * t274 + t332 * qJDD(2) + (qJ(3) * (mrSges(4,1) + mrSges(5,1)) + t315) * t275 + (t313 * t295 - t312 * t297) * qJD(1); t239; t241;];
tauJ = t1;
