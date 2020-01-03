% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:07
% EndTime: 2019-12-31 17:52:08
% DurationCPUTime: 0.63s
% Computational Cost: add. (1474->148), mult. (2586->176), div. (0->0), fcn. (936->6), ass. (0->68)
t318 = Ifges(5,1) + Ifges(6,1);
t312 = Ifges(5,4) + Ifges(6,4);
t311 = Ifges(5,5) + Ifges(6,5);
t317 = Ifges(5,2) + Ifges(6,2);
t310 = Ifges(5,6) + Ifges(6,6);
t316 = (Ifges(5,3) + Ifges(6,3));
t315 = -pkin(1) - pkin(2);
t290 = qJD(1) ^ 2;
t314 = pkin(4) * t290;
t313 = (-mrSges(5,2) - mrSges(6,2));
t287 = sin(qJ(1));
t289 = cos(qJ(1));
t296 = -t289 * g(1) - t287 * g(2);
t293 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t296;
t249 = t290 * t315 + t293;
t300 = t287 * g(1) - t289 * g(2);
t292 = -t290 * qJ(2) + qJDD(2) - t300;
t250 = qJDD(1) * t315 + t292;
t284 = sin(pkin(7));
t285 = cos(pkin(7));
t245 = t285 * t249 + t284 * t250;
t243 = -t290 * pkin(3) - qJDD(1) * pkin(6) + t245;
t282 = g(3) + qJDD(3);
t286 = sin(qJ(4));
t288 = cos(qJ(4));
t240 = t288 * t243 + t286 * t282;
t303 = qJD(1) * qJD(4);
t268 = -t288 * qJDD(1) + t286 * t303;
t304 = t286 * qJD(1);
t270 = (qJD(4) * pkin(4)) + qJ(5) * t304;
t281 = t288 ^ 2;
t302 = qJD(1) * qJD(5);
t237 = t268 * qJ(5) - qJD(4) * t270 - t281 * t314 - 0.2e1 * t288 * t302 + t240;
t309 = m(6) * t237 + t268 * mrSges(6,3);
t308 = (t316 * qJD(4)) + (-t286 * t311 - t310 * t288) * qJD(1);
t307 = -t310 * qJD(4) + (t286 * t312 + t288 * t317) * qJD(1);
t306 = t311 * qJD(4) + (-t318 * t286 - t288 * t312) * qJD(1);
t305 = qJD(1) * t288;
t301 = t288 * t303;
t276 = t288 * t282;
t239 = -t286 * t243 + t276;
t266 = (t288 * mrSges(5,1) - t286 * mrSges(5,2)) * qJD(1);
t267 = -t286 * qJDD(1) - t301;
t274 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t305;
t236 = (qJDD(4) * pkin(4)) + t276 + (-t267 - t301) * qJ(5) + (t288 * t314 - t243 + (2 * t302)) * t286;
t265 = (t288 * mrSges(6,1) - t286 * mrSges(6,2)) * qJD(1);
t273 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t305;
t298 = m(6) * t236 + qJDD(4) * mrSges(6,1) + qJD(4) * t273 + t265 * t304;
t230 = t266 * t304 + m(5) * t239 + qJDD(4) * mrSges(5,1) + qJD(4) * t274 + (-mrSges(5,3) - mrSges(6,3)) * t267 + t298;
t271 = qJD(4) * mrSges(6,1) + mrSges(6,3) * t304;
t272 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t304;
t231 = m(5) * t240 + t268 * mrSges(5,3) + (t313 * qJDD(4)) + (-t271 - t272) * qJD(4) + (-t265 - t266) * t305 + t309;
t299 = -t286 * t230 + t288 * t231;
t244 = -t284 * t249 + t285 * t250;
t295 = qJDD(1) * pkin(3) - t244;
t238 = -t270 * t304 - t268 * pkin(4) + qJDD(5) + (-qJ(5) * t281 - pkin(6)) * t290 + t295;
t297 = -m(6) * t238 + t268 * mrSges(6,1) + t271 * t304;
t227 = m(4) * t245 - t290 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t299;
t242 = -t290 * pkin(6) + t295;
t291 = -m(5) * t242 + t272 * t304 + t268 * mrSges(5,1) + t313 * t267 + (-t273 - t274) * t305 + t297;
t228 = m(4) * t244 - qJDD(1) * mrSges(4,1) - t290 * mrSges(4,2) + t291;
t294 = t284 * t227 + t285 * t228;
t252 = -qJDD(1) * pkin(1) + t292;
t251 = -t290 * pkin(1) + t293;
t233 = t267 * mrSges(6,2) + t273 * t305 - t297;
t232 = -t267 * mrSges(6,3) + t298;
t226 = m(3) * t252 - qJDD(1) * mrSges(3,1) - t290 * mrSges(3,3) + t294;
t1 = [qJ(2) * (m(3) * t251 - t290 * mrSges(3,1) + t285 * t227 - t284 * t228) - pkin(1) * t226 + mrSges(2,1) * t300 - mrSges(2,2) * t296 - pkin(2) * t294 - mrSges(3,1) * t252 + mrSges(3,3) * t251 - t286 * (mrSges(5,2) * t242 + mrSges(6,2) * t238 - mrSges(5,3) * t239 - mrSges(6,3) * t236 - qJ(5) * t232 + t307 * qJD(4) + t311 * qJDD(4) + t318 * t267 + t312 * t268 - t308 * t305) - t288 * (-mrSges(5,1) * t242 + mrSges(5,3) * t240 - mrSges(6,1) * t238 + mrSges(6,3) * t237 - pkin(4) * t233 + qJ(5) * t309 + t317 * t268 + t312 * t267 + (-qJ(5) * mrSges(6,2) + t310) * qJDD(4) + (-qJ(5) * t271 + t306) * qJD(4) + (-qJ(5) * t265 * t288 + t286 * t308) * qJD(1)) - pkin(3) * t291 - pkin(6) * t299 - mrSges(4,1) * t244 + mrSges(4,2) * t245 + (mrSges(3,3) * qJ(2) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1); t226; m(4) * t282 + t288 * t230 + t286 * t231; mrSges(5,1) * t239 + mrSges(6,1) * t236 - mrSges(5,2) * t240 - mrSges(6,2) * t237 + pkin(4) * t232 + t310 * t268 + t311 * t267 + (t316 * qJDD(4)) + (t286 * t307 + t306 * t288) * qJD(1); t233;];
tauJ = t1;
