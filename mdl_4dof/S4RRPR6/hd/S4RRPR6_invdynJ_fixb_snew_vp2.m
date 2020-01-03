% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRPR6
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR6_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR6_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:31
% EndTime: 2019-12-31 17:04:33
% DurationCPUTime: 0.81s
% Computational Cost: add. (5379->201), mult. (12462->266), div. (0->0), fcn. (8068->8), ass. (0->80)
t311 = qJD(1) ^ 2;
t323 = pkin(2) * t311;
t307 = sin(qJ(1));
t310 = cos(qJ(1));
t315 = -g(1) * t310 - g(2) * t307;
t290 = -pkin(1) * t311 + qJDD(1) * pkin(5) + t315;
t306 = sin(qJ(2));
t322 = t290 * t306;
t309 = cos(qJ(2));
t319 = qJD(1) * qJD(2);
t293 = qJDD(1) * t306 + t309 * t319;
t262 = qJDD(2) * pkin(2) - qJ(3) * t293 - t322 + (qJ(3) * t319 + t306 * t323 - g(3)) * t309;
t276 = -g(3) * t306 + t309 * t290;
t294 = qJDD(1) * t309 - t306 * t319;
t321 = qJD(1) * t306;
t295 = qJD(2) * pkin(2) - qJ(3) * t321;
t302 = t309 ^ 2;
t263 = qJ(3) * t294 - qJD(2) * t295 - t302 * t323 + t276;
t303 = sin(pkin(7));
t304 = cos(pkin(7));
t285 = (t303 * t309 + t304 * t306) * qJD(1);
t246 = -0.2e1 * qJD(3) * t285 + t304 * t262 - t263 * t303;
t274 = t293 * t304 + t294 * t303;
t284 = (-t303 * t306 + t304 * t309) * qJD(1);
t242 = (qJD(2) * t284 - t274) * pkin(6) + (t284 * t285 + qJDD(2)) * pkin(3) + t246;
t247 = 0.2e1 * qJD(3) * t284 + t303 * t262 + t304 * t263;
t273 = -t293 * t303 + t294 * t304;
t279 = qJD(2) * pkin(3) - pkin(6) * t285;
t283 = t284 ^ 2;
t243 = -pkin(3) * t283 + pkin(6) * t273 - qJD(2) * t279 + t247;
t305 = sin(qJ(4));
t308 = cos(qJ(4));
t240 = t242 * t308 - t243 * t305;
t270 = t284 * t308 - t285 * t305;
t253 = qJD(4) * t270 + t273 * t305 + t274 * t308;
t271 = t284 * t305 + t285 * t308;
t258 = -mrSges(5,1) * t270 + mrSges(5,2) * t271;
t301 = qJD(2) + qJD(4);
t265 = -t301 * mrSges(5,2) + t270 * mrSges(5,3);
t300 = qJDD(2) + qJDD(4);
t236 = m(5) * t240 + mrSges(5,1) * t300 - mrSges(5,3) * t253 - t258 * t271 + t265 * t301;
t241 = t242 * t305 + t243 * t308;
t252 = -qJD(4) * t271 + t273 * t308 - t274 * t305;
t266 = mrSges(5,1) * t301 - mrSges(5,3) * t271;
t237 = m(5) * t241 - mrSges(5,2) * t300 + mrSges(5,3) * t252 + t258 * t270 - t266 * t301;
t230 = t308 * t236 + t305 * t237;
t272 = -mrSges(4,1) * t284 + mrSges(4,2) * t285;
t277 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t284;
t228 = m(4) * t246 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t274 + qJD(2) * t277 - t272 * t285 + t230;
t278 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t285;
t316 = -t236 * t305 + t308 * t237;
t229 = m(4) * t247 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t273 - qJD(2) * t278 + t272 * t284 + t316;
t224 = t304 * t228 + t303 * t229;
t320 = qJD(1) * t309;
t318 = g(1) * t307 - t310 * g(2);
t317 = -t228 * t303 + t304 * t229;
t314 = -qJDD(1) * pkin(1) - t318;
t264 = -pkin(2) * t294 + qJDD(3) + t295 * t321 + (-qJ(3) * t302 - pkin(5)) * t311 + t314;
t245 = -pkin(3) * t273 - pkin(6) * t283 + t279 * t285 + t264;
t313 = -m(5) * t245 + mrSges(5,1) * t252 - t253 * mrSges(5,2) + t265 * t270 - t271 * t266;
t255 = Ifges(5,4) * t271 + Ifges(5,2) * t270 + Ifges(5,6) * t301;
t256 = Ifges(5,1) * t271 + Ifges(5,4) * t270 + Ifges(5,5) * t301;
t312 = mrSges(5,1) * t240 - mrSges(5,2) * t241 + Ifges(5,5) * t253 + Ifges(5,6) * t252 + Ifges(5,3) * t300 + t271 * t255 - t270 * t256;
t238 = m(4) * t264 - mrSges(4,1) * t273 + mrSges(4,2) * t274 - t277 * t284 + t278 * t285 - t313;
t297 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t320;
t296 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t321;
t292 = (-mrSges(3,1) * t309 + mrSges(3,2) * t306) * qJD(1);
t289 = -pkin(5) * t311 + t314;
t288 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t306 + Ifges(3,4) * t309) * qJD(1);
t287 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t306 + Ifges(3,2) * t309) * qJD(1);
t275 = -g(3) * t309 - t322;
t269 = Ifges(4,1) * t285 + Ifges(4,4) * t284 + Ifges(4,5) * qJD(2);
t268 = Ifges(4,4) * t285 + Ifges(4,2) * t284 + Ifges(4,6) * qJD(2);
t267 = Ifges(4,5) * t285 + Ifges(4,6) * t284 + Ifges(4,3) * qJD(2);
t254 = Ifges(5,5) * t271 + Ifges(5,6) * t270 + Ifges(5,3) * t301;
t232 = mrSges(5,2) * t245 - mrSges(5,3) * t240 + Ifges(5,1) * t253 + Ifges(5,4) * t252 + Ifges(5,5) * t300 + t254 * t270 - t255 * t301;
t231 = -mrSges(5,1) * t245 + mrSges(5,3) * t241 + Ifges(5,4) * t253 + Ifges(5,2) * t252 + Ifges(5,6) * t300 - t254 * t271 + t256 * t301;
t223 = mrSges(4,2) * t264 - mrSges(4,3) * t246 + Ifges(4,1) * t274 + Ifges(4,4) * t273 + Ifges(4,5) * qJDD(2) - pkin(6) * t230 - qJD(2) * t268 - t231 * t305 + t232 * t308 + t267 * t284;
t222 = -mrSges(4,1) * t264 + mrSges(4,3) * t247 + Ifges(4,4) * t274 + Ifges(4,2) * t273 + Ifges(4,6) * qJDD(2) + pkin(3) * t313 + pkin(6) * t316 + qJD(2) * t269 + t308 * t231 + t305 * t232 - t285 * t267;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t318 - mrSges(2,2) * t315 + t306 * (mrSges(3,2) * t289 - mrSges(3,3) * t275 + Ifges(3,1) * t293 + Ifges(3,4) * t294 + Ifges(3,5) * qJDD(2) - qJ(3) * t224 - qJD(2) * t287 - t303 * t222 + t304 * t223) + t309 * (-mrSges(3,1) * t289 + mrSges(3,3) * t276 + Ifges(3,4) * t293 + Ifges(3,2) * t294 + Ifges(3,6) * qJDD(2) - pkin(2) * t238 + qJ(3) * t317 + qJD(2) * t288 + t304 * t222 + t303 * t223) + pkin(1) * (-t238 - m(3) * t289 + mrSges(3,1) * t294 - mrSges(3,2) * t293 + (-t296 * t306 + t297 * t309) * qJD(1)) + pkin(5) * (t309 * (m(3) * t276 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t294 - qJD(2) * t296 + t292 * t320 + t317) - t306 * (m(3) * t275 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t293 + qJD(2) * t297 - t292 * t321 + t224)); pkin(2) * t224 + t312 + Ifges(3,5) * t293 + Ifges(3,6) * t294 - t284 * t269 + t285 * t268 + Ifges(4,6) * t273 + Ifges(4,5) * t274 + mrSges(3,1) * t275 - mrSges(3,2) * t276 + pkin(3) * t230 + (t306 * t287 - t309 * t288) * qJD(1) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + mrSges(4,1) * t246 - mrSges(4,2) * t247; t238; t312;];
tauJ = t1;
