% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRPP3
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRPP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP3_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:27
% EndTime: 2019-12-31 16:57:28
% DurationCPUTime: 0.70s
% Computational Cost: add. (2015->180), mult. (4575->223), div. (0->0), fcn. (2637->6), ass. (0->73)
t323 = Ifges(4,1) + Ifges(5,1);
t316 = Ifges(4,4) - Ifges(5,5);
t315 = Ifges(4,5) + Ifges(5,4);
t322 = -Ifges(4,2) - Ifges(5,3);
t321 = -Ifges(5,2) - Ifges(4,3);
t314 = Ifges(4,6) - Ifges(5,6);
t289 = sin(qJ(2));
t291 = cos(qJ(2));
t304 = qJD(1) * qJD(2);
t279 = t291 * qJDD(1) - t289 * t304;
t306 = qJD(1) * t289;
t280 = qJD(2) * pkin(2) - qJ(3) * t306;
t287 = t291 ^ 2;
t294 = qJD(1) ^ 2;
t290 = sin(qJ(1));
t292 = cos(qJ(1));
t301 = t290 * g(1) - t292 * g(2);
t297 = -qJDD(1) * pkin(1) - t301;
t244 = -t279 * pkin(2) + qJDD(3) + t280 * t306 + (-qJ(3) * t287 - pkin(5)) * t294 + t297;
t278 = t289 * qJDD(1) + t291 * t304;
t288 = sin(pkin(6));
t313 = cos(pkin(6));
t256 = t288 * t278 - t313 * t279;
t257 = t313 * t278 + t288 * t279;
t305 = qJD(1) * t291;
t268 = t288 * t306 - t313 * t305;
t262 = -qJD(2) * mrSges(4,2) - t268 * mrSges(4,3);
t269 = (t288 * t291 + t313 * t289) * qJD(1);
t263 = qJD(2) * mrSges(4,1) - t269 * mrSges(4,3);
t264 = -qJD(2) * mrSges(5,1) + t269 * mrSges(5,2);
t237 = -0.2e1 * qJD(4) * t269 + (qJD(2) * t268 - t257) * qJ(4) + (qJD(2) * t269 + t256) * pkin(3) + t244;
t265 = -t268 * mrSges(5,2) + qJD(2) * mrSges(5,3);
t302 = m(5) * t237 + t256 * mrSges(5,1) + t268 * t265;
t230 = m(4) * t244 + t256 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t257 + t268 * t262 + (t263 - t264) * t269 + t302;
t320 = -2 * qJD(3);
t319 = pkin(2) * t294;
t317 = -mrSges(4,3) - mrSges(5,2);
t298 = -t292 * g(1) - t290 * g(2);
t275 = -t294 * pkin(1) + qJDD(1) * pkin(5) + t298;
t312 = t289 * t275;
t242 = qJDD(2) * pkin(2) - t278 * qJ(3) - t312 + (qJ(3) * t304 + t289 * t319 - g(3)) * t291;
t261 = -t289 * g(3) + t291 * t275;
t243 = t279 * qJ(3) - qJD(2) * t280 - t287 * t319 + t261;
t239 = t288 * t242 + t313 * t243 + t268 * t320;
t252 = t268 * pkin(3) - t269 * qJ(4);
t293 = qJD(2) ^ 2;
t234 = -t293 * pkin(3) + qJDD(2) * qJ(4) + 0.2e1 * qJD(4) * qJD(2) - t268 * t252 + t239;
t303 = m(5) * t234 + qJDD(2) * mrSges(5,3) + qJD(2) * t264;
t253 = t268 * mrSges(5,1) - t269 * mrSges(5,3);
t308 = -t268 * mrSges(4,1) - t269 * mrSges(4,2) - t253;
t228 = m(4) * t239 - qJDD(2) * mrSges(4,2) - qJD(2) * t263 + t317 * t256 + t308 * t268 + t303;
t296 = t313 * t242 - t288 * t243;
t238 = t269 * t320 + t296;
t235 = -qJDD(2) * pkin(3) - t293 * qJ(4) + qJDD(4) + ((2 * qJD(3)) + t252) * t269 - t296;
t299 = -m(5) * t235 + qJDD(2) * mrSges(5,1) + qJD(2) * t265;
t229 = m(4) * t238 + qJDD(2) * mrSges(4,1) + qJD(2) * t262 + t317 * t257 + t308 * t269 + t299;
t224 = t288 * t228 + t313 * t229;
t311 = t321 * qJD(2) + t314 * t268 - t315 * t269;
t310 = t314 * qJD(2) + t322 * t268 + t316 * t269;
t309 = t315 * qJD(2) - t316 * t268 + t323 * t269;
t300 = t313 * t228 - t288 * t229;
t282 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t305;
t281 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t306;
t277 = (-t291 * mrSges(3,1) + t289 * mrSges(3,2)) * qJD(1);
t274 = -t294 * pkin(5) + t297;
t272 = Ifges(3,5) * qJD(2) + (t289 * Ifges(3,1) + t291 * Ifges(3,4)) * qJD(1);
t271 = Ifges(3,6) * qJD(2) + (t289 * Ifges(3,4) + t291 * Ifges(3,2)) * qJD(1);
t260 = -t291 * g(3) - t312;
t232 = t257 * mrSges(5,2) + t269 * t253 - t299;
t231 = -t257 * mrSges(5,3) - t269 * t264 + t302;
t223 = mrSges(4,2) * t244 + mrSges(5,2) * t235 - mrSges(4,3) * t238 - mrSges(5,3) * t237 - qJ(4) * t231 - t310 * qJD(2) + t315 * qJDD(2) - t316 * t256 + t323 * t257 + t311 * t268;
t222 = -mrSges(4,1) * t244 - mrSges(5,1) * t237 + mrSges(5,2) * t234 + mrSges(4,3) * t239 - pkin(3) * t231 + t309 * qJD(2) + t314 * qJDD(2) + t322 * t256 + t316 * t257 + t311 * t269;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t301 - mrSges(2,2) * t298 + t289 * (mrSges(3,2) * t274 - mrSges(3,3) * t260 + Ifges(3,1) * t278 + Ifges(3,4) * t279 + Ifges(3,5) * qJDD(2) - qJ(3) * t224 - qJD(2) * t271 - t288 * t222 + t313 * t223) + t291 * (-mrSges(3,1) * t274 + mrSges(3,3) * t261 + Ifges(3,4) * t278 + Ifges(3,2) * t279 + Ifges(3,6) * qJDD(2) - pkin(2) * t230 + qJ(3) * t300 + qJD(2) * t272 + t313 * t222 + t288 * t223) + pkin(1) * (-m(3) * t274 + t279 * mrSges(3,1) - t278 * mrSges(3,2) + (-t281 * t289 + t282 * t291) * qJD(1) - t230) + pkin(5) * (t291 * (m(3) * t261 - qJDD(2) * mrSges(3,2) + t279 * mrSges(3,3) - qJD(2) * t281 + t277 * t305 + t300) - t289 * (m(3) * t260 + qJDD(2) * mrSges(3,1) - t278 * mrSges(3,3) + qJD(2) * t282 - t277 * t306 + t224)); Ifges(3,5) * t278 + Ifges(3,6) * t279 + mrSges(3,1) * t260 - mrSges(3,2) * t261 + mrSges(4,1) * t238 - mrSges(4,2) * t239 - mrSges(5,1) * t235 + mrSges(5,3) * t234 - pkin(3) * t232 + qJ(4) * t303 + pkin(2) * t224 + t310 * t269 + (-qJ(4) * t253 + t309) * t268 + t315 * t257 + (-qJ(4) * mrSges(5,2) - t314) * t256 + (t289 * t271 - t291 * t272) * qJD(1) + (Ifges(3,3) - t321) * qJDD(2); t230; t232;];
tauJ = t1;
