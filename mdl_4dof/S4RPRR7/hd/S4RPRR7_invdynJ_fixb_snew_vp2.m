% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPRR7
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR7_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR7_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:48
% EndTime: 2019-12-31 16:53:49
% DurationCPUTime: 0.66s
% Computational Cost: add. (3644->171), mult. (8564->224), div. (0->0), fcn. (5736->8), ass. (0->79)
t296 = qJD(1) ^ 2;
t287 = sin(pkin(7));
t288 = cos(pkin(7));
t290 = sin(qJ(3));
t293 = cos(qJ(3));
t301 = t287 * t290 - t288 * t293;
t277 = t301 * qJD(1);
t302 = t287 * t293 + t288 * t290;
t278 = t302 * qJD(1);
t311 = t278 * qJD(3);
t266 = -qJDD(1) * t301 - t311;
t317 = pkin(2) * t296;
t316 = pkin(5) * qJDD(1);
t291 = sin(qJ(1));
t294 = cos(qJ(1));
t305 = -t294 * g(1) - t291 * g(2);
t279 = -t296 * pkin(1) + qJDD(1) * qJ(2) + t305;
t310 = qJD(1) * qJD(2);
t308 = -t288 * g(3) - 0.2e1 * t287 * t310;
t257 = (t288 * t317 - t279 - t316) * t287 + t308;
t269 = -t287 * g(3) + (t279 + 0.2e1 * t310) * t288;
t286 = t288 ^ 2;
t258 = -t286 * t317 + t288 * t316 + t269;
t244 = t290 * t257 + t293 * t258;
t262 = t277 * mrSges(4,1) + t278 * mrSges(4,2);
t273 = qJD(3) * mrSges(4,1) - t278 * mrSges(4,3);
t264 = t277 * pkin(3) - t278 * pkin(6);
t295 = qJD(3) ^ 2;
t241 = -t295 * pkin(3) + qJDD(3) * pkin(6) - t277 * t264 + t244;
t309 = t291 * g(1) - t294 * g(2);
t304 = qJDD(2) - t309;
t314 = -t287 ^ 2 - t286;
t265 = (-pkin(2) * t288 - pkin(1)) * qJDD(1) + (pkin(5) * t314 - qJ(2)) * t296 + t304;
t312 = t277 * qJD(3);
t267 = qJDD(1) * t302 - t312;
t242 = (-t267 + t312) * pkin(6) + (-t266 + t311) * pkin(3) + t265;
t289 = sin(qJ(4));
t292 = cos(qJ(4));
t238 = -t289 * t241 + t292 * t242;
t270 = t292 * qJD(3) - t289 * t278;
t251 = t270 * qJD(4) + t289 * qJDD(3) + t292 * t267;
t271 = t289 * qJD(3) + t292 * t278;
t252 = -t270 * mrSges(5,1) + t271 * mrSges(5,2);
t275 = qJD(4) + t277;
t254 = -t275 * mrSges(5,2) + t270 * mrSges(5,3);
t263 = qJDD(4) - t266;
t236 = m(5) * t238 + t263 * mrSges(5,1) - t251 * mrSges(5,3) - t271 * t252 + t275 * t254;
t239 = t292 * t241 + t289 * t242;
t250 = -t271 * qJD(4) + t292 * qJDD(3) - t289 * t267;
t255 = t275 * mrSges(5,1) - t271 * mrSges(5,3);
t237 = m(5) * t239 - t263 * mrSges(5,2) + t250 * mrSges(5,3) + t270 * t252 - t275 * t255;
t306 = -t289 * t236 + t292 * t237;
t227 = m(4) * t244 - qJDD(3) * mrSges(4,2) + t266 * mrSges(4,3) - qJD(3) * t273 - t277 * t262 + t306;
t243 = t293 * t257 - t290 * t258;
t272 = -qJD(3) * mrSges(4,2) - t277 * mrSges(4,3);
t240 = -qJDD(3) * pkin(3) - t295 * pkin(6) + t278 * t264 - t243;
t299 = -m(5) * t240 + t250 * mrSges(5,1) - t251 * mrSges(5,2) + t270 * t254 - t271 * t255;
t232 = m(4) * t243 + qJDD(3) * mrSges(4,1) - t267 * mrSges(4,3) + qJD(3) * t272 - t278 * t262 + t299;
t315 = t290 * t227 + t293 * t232;
t228 = t292 * t236 + t289 * t237;
t307 = t293 * t227 - t290 * t232;
t303 = -t288 * mrSges(3,1) + t287 * mrSges(3,2);
t300 = mrSges(3,3) * qJDD(1) + t296 * t303;
t298 = m(4) * t265 - t266 * mrSges(4,1) + t267 * mrSges(4,2) + t277 * t272 + t278 * t273 + t228;
t246 = Ifges(5,4) * t271 + Ifges(5,2) * t270 + Ifges(5,6) * t275;
t247 = Ifges(5,1) * t271 + Ifges(5,4) * t270 + Ifges(5,5) * t275;
t297 = mrSges(5,1) * t238 - mrSges(5,2) * t239 + Ifges(5,5) * t251 + Ifges(5,6) * t250 + Ifges(5,3) * t263 + t271 * t246 - t270 * t247;
t276 = -qJDD(1) * pkin(1) - t296 * qJ(2) + t304;
t268 = -t287 * t279 + t308;
t261 = Ifges(4,1) * t278 - Ifges(4,4) * t277 + Ifges(4,5) * qJD(3);
t260 = Ifges(4,4) * t278 - Ifges(4,2) * t277 + Ifges(4,6) * qJD(3);
t259 = Ifges(4,5) * t278 - Ifges(4,6) * t277 + Ifges(4,3) * qJD(3);
t245 = Ifges(5,5) * t271 + Ifges(5,6) * t270 + Ifges(5,3) * t275;
t230 = mrSges(5,2) * t240 - mrSges(5,3) * t238 + Ifges(5,1) * t251 + Ifges(5,4) * t250 + Ifges(5,5) * t263 + t270 * t245 - t275 * t246;
t229 = -mrSges(5,1) * t240 + mrSges(5,3) * t239 + Ifges(5,4) * t251 + Ifges(5,2) * t250 + Ifges(5,6) * t263 - t271 * t245 + t275 * t247;
t224 = mrSges(3,3) * t296 * t314 + m(3) * t276 + qJDD(1) * t303 + t298;
t223 = -mrSges(4,1) * t265 + mrSges(4,3) * t244 + Ifges(4,4) * t267 + Ifges(4,2) * t266 + Ifges(4,6) * qJDD(3) - pkin(3) * t228 + qJD(3) * t261 - t278 * t259 - t297;
t222 = mrSges(4,2) * t265 - mrSges(4,3) * t243 + Ifges(4,1) * t267 + Ifges(4,4) * t266 + Ifges(4,5) * qJDD(3) - pkin(6) * t228 - qJD(3) * t260 - t289 * t229 + t292 * t230 - t277 * t259;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t309 - mrSges(2,2) * t305 + t287 * (mrSges(3,2) * t276 - mrSges(3,3) * t268 + t293 * t222 - t290 * t223 - pkin(5) * t315 + (Ifges(3,1) * t287 + Ifges(3,4) * t288) * qJDD(1)) + t288 * (-mrSges(3,1) * t276 + mrSges(3,3) * t269 + t290 * t222 + t293 * t223 - pkin(2) * t298 + pkin(5) * t307 + (Ifges(3,4) * t287 + Ifges(3,2) * t288) * qJDD(1)) - pkin(1) * t224 + qJ(2) * ((m(3) * t269 + t288 * t300 + t307) * t288 + (-m(3) * t268 + t287 * t300 - t315) * t287); t224; mrSges(4,1) * t243 - mrSges(4,2) * t244 + Ifges(4,5) * t267 + Ifges(4,6) * t266 + Ifges(4,3) * qJDD(3) + pkin(3) * t299 + pkin(6) * t306 + t292 * t229 + t289 * t230 + t278 * t260 + t277 * t261; t297;];
tauJ = t1;
