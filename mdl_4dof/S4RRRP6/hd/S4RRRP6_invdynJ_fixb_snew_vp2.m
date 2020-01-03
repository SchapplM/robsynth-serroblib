% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRRP6
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP6_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:15
% EndTime: 2019-12-31 17:18:17
% DurationCPUTime: 0.81s
% Computational Cost: add. (2442->178), mult. (4805->217), div. (0->0), fcn. (2766->6), ass. (0->73)
t317 = Ifges(4,1) + Ifges(5,1);
t309 = Ifges(4,4) + Ifges(5,4);
t308 = Ifges(4,5) + Ifges(5,5);
t316 = Ifges(4,2) + Ifges(5,2);
t307 = Ifges(4,6) + Ifges(5,6);
t282 = sin(qJ(3));
t285 = cos(qJ(3));
t283 = sin(qJ(2));
t302 = qJD(1) * t283;
t271 = qJD(2) * t285 - t282 * t302;
t286 = cos(qJ(2));
t300 = qJD(1) * qJD(2);
t296 = t286 * t300;
t275 = qJDD(1) * t283 + t296;
t249 = qJD(3) * t271 + qJDD(2) * t282 + t275 * t285;
t272 = qJD(2) * t282 + t285 * t302;
t251 = -mrSges(5,1) * t271 + mrSges(5,2) * t272;
t289 = qJD(1) ^ 2;
t284 = sin(qJ(1));
t287 = cos(qJ(1));
t295 = g(1) * t284 - t287 * g(2);
t266 = -qJDD(1) * pkin(1) - pkin(5) * t289 - t295;
t297 = t283 * t300;
t276 = qJDD(1) * t286 - t297;
t232 = (-t275 - t296) * pkin(6) + (-t276 + t297) * pkin(2) + t266;
t292 = -g(1) * t287 - g(2) * t284;
t267 = -pkin(1) * t289 + qJDD(1) * pkin(5) + t292;
t259 = -g(3) * t283 + t286 * t267;
t274 = (-pkin(2) * t286 - pkin(6) * t283) * qJD(1);
t288 = qJD(2) ^ 2;
t301 = qJD(1) * t286;
t235 = -pkin(2) * t288 + qJDD(2) * pkin(6) + t274 * t301 + t259;
t228 = t285 * t232 - t235 * t282;
t270 = qJDD(3) - t276;
t279 = qJD(3) - t301;
t224 = -0.2e1 * qJD(4) * t272 + (t271 * t279 - t249) * qJ(4) + (t271 * t272 + t270) * pkin(3) + t228;
t253 = -mrSges(5,2) * t279 + mrSges(5,3) * t271;
t299 = m(5) * t224 + t270 * mrSges(5,1) + t279 * t253;
t221 = -mrSges(5,3) * t249 - t251 * t272 + t299;
t229 = t282 * t232 + t285 * t235;
t248 = -qJD(3) * t272 + qJDD(2) * t285 - t275 * t282;
t255 = pkin(3) * t279 - qJ(4) * t272;
t269 = t271 ^ 2;
t226 = -pkin(3) * t269 + qJ(4) * t248 + 0.2e1 * qJD(4) * t271 - t255 * t279 + t229;
t305 = -t316 * t271 - t309 * t272 - t307 * t279;
t311 = t309 * t271 + t317 * t272 + t308 * t279;
t313 = Ifges(4,3) + Ifges(5,3);
t315 = mrSges(4,1) * t228 + mrSges(5,1) * t224 - mrSges(4,2) * t229 - mrSges(5,2) * t226 + pkin(3) * t221 + t307 * t248 + t308 * t249 + t313 * t270 - t311 * t271 - t305 * t272;
t310 = -mrSges(4,2) - mrSges(5,2);
t306 = -t307 * t271 - t308 * t272 - t313 * t279;
t256 = mrSges(5,1) * t279 - mrSges(5,3) * t272;
t303 = -mrSges(4,1) * t279 + mrSges(4,3) * t272 - t256;
t258 = -t286 * g(3) - t283 * t267;
t298 = m(5) * t226 + t248 * mrSges(5,3) + t271 * t251;
t252 = -mrSges(4,1) * t271 + mrSges(4,2) * t272;
t254 = -mrSges(4,2) * t279 + mrSges(4,3) * t271;
t218 = m(4) * t228 + mrSges(4,1) * t270 + t254 * t279 + (-t251 - t252) * t272 + (-mrSges(4,3) - mrSges(5,3)) * t249 + t299;
t220 = m(4) * t229 + mrSges(4,3) * t248 + t252 * t271 + t270 * t310 + t303 * t279 + t298;
t294 = -t218 * t282 + t285 * t220;
t234 = -qJDD(2) * pkin(2) - pkin(6) * t288 + t274 * t302 - t258;
t227 = -pkin(3) * t248 - qJ(4) * t269 + t255 * t272 + qJDD(4) + t234;
t293 = -m(5) * t227 + t248 * mrSges(5,1) + t271 * t253;
t217 = t218 * t285 + t220 * t282;
t290 = -m(4) * t234 + t248 * mrSges(4,1) + t249 * t310 + t271 * t254 + t303 * t272 + t293;
t278 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t301;
t277 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t302;
t273 = (-mrSges(3,1) * t286 + mrSges(3,2) * t283) * qJD(1);
t265 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t283 + Ifges(3,4) * t286) * qJD(1);
t264 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t283 + Ifges(3,2) * t286) * qJD(1);
t222 = mrSges(5,2) * t249 + t256 * t272 - t293;
t216 = mrSges(4,2) * t234 + mrSges(5,2) * t227 - mrSges(4,3) * t228 - mrSges(5,3) * t224 - qJ(4) * t221 + t309 * t248 + t317 * t249 + t308 * t270 - t306 * t271 + t305 * t279;
t215 = -mrSges(4,1) * t234 + mrSges(4,3) * t229 - mrSges(5,1) * t227 + mrSges(5,3) * t226 - pkin(3) * t222 + qJ(4) * t298 + (-qJ(4) * t256 + t311) * t279 + t306 * t272 + (-mrSges(5,2) * qJ(4) + t307) * t270 + t309 * t249 + t316 * t248;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t295 - mrSges(2,2) * t292 + t283 * (mrSges(3,2) * t266 - mrSges(3,3) * t258 + Ifges(3,1) * t275 + Ifges(3,4) * t276 + Ifges(3,5) * qJDD(2) - pkin(6) * t217 - qJD(2) * t264 - t282 * t215 + t285 * t216) + t286 * (-mrSges(3,1) * t266 + mrSges(3,3) * t259 + Ifges(3,4) * t275 + Ifges(3,2) * t276 + Ifges(3,6) * qJDD(2) - pkin(2) * t217 + qJD(2) * t265 - t315) + pkin(1) * (-m(3) * t266 + mrSges(3,1) * t276 - mrSges(3,2) * t275 + (-t277 * t283 + t278 * t286) * qJD(1) - t217) + pkin(5) * (t286 * (m(3) * t259 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t276 - qJD(2) * t277 + t273 * t301 + t294) - t283 * (m(3) * t258 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t275 + qJD(2) * t278 - t273 * t302 + t290)); Ifges(3,5) * t275 + Ifges(3,6) * t276 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t258 - mrSges(3,2) * t259 + t282 * t216 + t285 * t215 + pkin(2) * t290 + pkin(6) * t294 + (t283 * t264 - t286 * t265) * qJD(1); t315; t222;];
tauJ = t1;
