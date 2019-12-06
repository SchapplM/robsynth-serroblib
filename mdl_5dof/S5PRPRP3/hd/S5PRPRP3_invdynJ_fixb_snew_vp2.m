% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRP3
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:49
% EndTime: 2019-12-05 15:32:50
% DurationCPUTime: 0.61s
% Computational Cost: add. (1204->139), mult. (2096->169), div. (0->0), fcn. (1100->8), ass. (0->66)
t297 = cos(qJ(4));
t315 = Ifges(5,4) + Ifges(6,4);
t322 = t297 * t315;
t321 = Ifges(5,1) + Ifges(6,1);
t314 = Ifges(5,5) + Ifges(6,5);
t320 = Ifges(5,2) + Ifges(6,2);
t319 = Ifges(5,6) + Ifges(6,6);
t295 = sin(qJ(4));
t318 = -t314 * qJD(4) + (-t321 * t295 - t322) * qJD(2);
t299 = qJD(2) ^ 2;
t317 = pkin(4) * t299;
t316 = -mrSges(5,2) - mrSges(6,2);
t292 = sin(pkin(7));
t294 = cos(pkin(7));
t280 = -g(1) * t294 - g(2) * t292;
t290 = -g(3) + qJDD(1);
t296 = sin(qJ(2));
t298 = cos(qJ(2));
t257 = -t280 * t296 + t298 * t290;
t255 = qJDD(2) * pkin(2) + t257;
t258 = t298 * t280 + t296 * t290;
t256 = -pkin(2) * t299 + t258;
t291 = sin(pkin(8));
t293 = cos(pkin(8));
t251 = t291 * t255 + t293 * t256;
t249 = -pkin(3) * t299 + qJDD(2) * pkin(6) + t251;
t279 = -g(1) * t292 + g(2) * t294 + qJDD(3);
t246 = t297 * t249 + t295 * t279;
t275 = (-mrSges(5,1) * t297 + mrSges(5,2) * t295) * qJD(2);
t308 = qJD(2) * qJD(4);
t277 = qJDD(2) * t297 - t295 * t308;
t310 = qJD(2) * t295;
t281 = qJD(4) * pkin(4) - qJ(5) * t310;
t289 = t297 ^ 2;
t307 = qJD(2) * qJD(5);
t243 = qJ(5) * t277 - qJD(4) * t281 - t289 * t317 + 0.2e1 * t297 * t307 + t246;
t274 = (-mrSges(6,1) * t297 + mrSges(6,2) * t295) * qJD(2);
t309 = qJD(2) * t297;
t305 = m(6) * t243 + t277 * mrSges(6,3) + t274 * t309;
t282 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t310;
t311 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t310 - t282;
t237 = m(5) * t246 + t277 * mrSges(5,3) + t311 * qJD(4) + t316 * qJDD(4) + t275 * t309 + t305;
t235 = t297 * t237;
t267 = t297 * t279;
t245 = -t295 * t249 + t267;
t304 = t297 * t308;
t276 = qJDD(2) * t295 + t304;
t285 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t309;
t242 = qJDD(4) * pkin(4) + t267 + (-t276 + t304) * qJ(5) + (t297 * t317 - t249 - 0.2e1 * t307) * t295;
t284 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t309;
t306 = m(6) * t242 + qJDD(4) * mrSges(6,1) + qJD(4) * t284;
t236 = m(5) * t245 + qJDD(4) * mrSges(5,1) + qJD(4) * t285 + (-mrSges(5,3) - mrSges(6,3)) * t276 + (-t274 - t275) * t310 + t306;
t232 = m(4) * t251 - mrSges(4,1) * t299 - qJDD(2) * mrSges(4,2) - t236 * t295 + t235;
t250 = t293 * t255 - t291 * t256;
t301 = -qJDD(2) * pkin(3) - t250;
t248 = -pkin(6) * t299 + t301;
t244 = t281 * t310 - t277 * pkin(4) + qJDD(5) + (-qJ(5) * t289 - pkin(6)) * t299 + t301;
t302 = m(6) * t244 - t277 * mrSges(6,1) - t284 * t309;
t300 = -m(5) * t248 + t277 * mrSges(5,1) + t285 * t309 - t302;
t303 = qJD(2) * t311;
t234 = m(4) * t250 + qJDD(2) * mrSges(4,1) - t299 * mrSges(4,2) + t316 * t276 + t295 * t303 + t300;
t313 = t291 * t232 + t293 * t234;
t312 = t319 * qJD(4) + (t315 * t295 + t297 * t320) * qJD(2);
t239 = t276 * mrSges(6,2) + t282 * t310 + t302;
t238 = -t276 * mrSges(6,3) - t274 * t310 + t306;
t1 = [m(2) * t290 + t296 * (m(3) * t258 - mrSges(3,1) * t299 - qJDD(2) * mrSges(3,2) + t232 * t293 - t234 * t291) + t298 * (m(3) * t257 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t299 + t313); mrSges(3,1) * t257 - mrSges(3,2) * t258 + mrSges(4,1) * t250 - mrSges(4,2) * t251 + t297 * (-mrSges(5,1) * t248 - mrSges(6,1) * t244 + mrSges(5,3) * t246 + mrSges(6,3) * t243 - pkin(4) * t239 + qJ(5) * t305 + t320 * t277 + (-qJ(5) * mrSges(6,2) + t319) * qJDD(4) + (-qJ(5) * t282 - t318) * qJD(4)) + pkin(3) * t300 + pkin(6) * t235 + pkin(2) * t313 + (pkin(3) * t316 + t322) * t276 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (mrSges(5,2) * t248 + mrSges(6,2) * t244 - mrSges(5,3) * t245 - mrSges(6,3) * t242 + pkin(3) * t303 - pkin(6) * t236 - qJ(5) * t238 - t312 * qJD(4) + t314 * qJDD(4) + t321 * t276 + t315 * t277) * t295; m(4) * t279 + t236 * t297 + t237 * t295; mrSges(5,1) * t245 + mrSges(6,1) * t242 - mrSges(5,2) * t246 - mrSges(6,2) * t243 + pkin(4) * t238 + t319 * t277 + t314 * t276 + (Ifges(5,3) + Ifges(6,3)) * qJDD(4) + (t312 * t295 + t318 * t297) * qJD(2); t239;];
tauJ = t1;
