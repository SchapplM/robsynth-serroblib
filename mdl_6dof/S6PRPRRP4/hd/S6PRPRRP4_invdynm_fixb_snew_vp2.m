% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 23:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:47:26
% EndTime: 2019-05-04 23:47:46
% DurationCPUTime: 12.57s
% Computational Cost: add. (219164->316), mult. (481223->390), div. (0->0), fcn. (356395->12), ass. (0->136)
t307 = qJD(2) ^ 2;
t296 = sin(pkin(10));
t299 = cos(pkin(10));
t281 = t296 * g(1) - t299 * g(2);
t282 = -t299 * g(1) - t296 * g(2);
t294 = -g(3) + qJDD(1);
t305 = cos(qJ(2));
t300 = cos(pkin(6));
t303 = sin(qJ(2));
t341 = t300 * t303;
t297 = sin(pkin(6));
t342 = t297 * t303;
t248 = t281 * t341 + t305 * t282 + t294 * t342;
t243 = -t307 * pkin(2) + qJDD(2) * qJ(3) + t248;
t295 = sin(pkin(11));
t268 = -t297 * t281 + t300 * t294;
t298 = cos(pkin(11));
t333 = qJD(2) * qJD(3);
t337 = t298 * t268 - 0.2e1 * t295 * t333;
t349 = pkin(3) * t298;
t208 = (-pkin(8) * qJDD(2) + t307 * t349 - t243) * t295 + t337;
t212 = t295 * t268 + (t243 + 0.2e1 * t333) * t298;
t332 = qJDD(2) * t298;
t292 = t298 ^ 2;
t343 = t292 * t307;
t209 = -pkin(3) * t343 + pkin(8) * t332 + t212;
t302 = sin(qJ(4));
t304 = cos(qJ(4));
t201 = t302 * t208 + t304 * t209;
t319 = t295 * t302 - t298 * t304;
t271 = t319 * qJD(2);
t320 = t295 * t304 + t298 * t302;
t272 = t320 * qJD(2);
t254 = t271 * mrSges(5,1) + t272 * mrSges(5,2);
t334 = t272 * qJD(4);
t260 = -t319 * qJDD(2) - t334;
t267 = qJD(4) * mrSges(5,1) - t272 * mrSges(5,3);
t259 = t271 * pkin(4) - t272 * pkin(9);
t306 = qJD(4) ^ 2;
t198 = -t306 * pkin(4) + qJDD(4) * pkin(9) - t271 * t259 + t201;
t291 = t295 ^ 2;
t247 = -t303 * t282 + (t281 * t300 + t294 * t297) * t305;
t314 = qJDD(3) - t247;
t229 = (-pkin(2) - t349) * qJDD(2) + (-qJ(3) + (-t291 - t292) * pkin(8)) * t307 + t314;
t335 = t271 * qJD(4);
t261 = t320 * qJDD(2) - t335;
t203 = (-t261 + t335) * pkin(9) + (-t260 + t334) * pkin(4) + t229;
t301 = sin(qJ(5));
t350 = cos(qJ(5));
t195 = t350 * t198 + t301 * t203;
t263 = t301 * qJD(4) + t350 * t272;
t227 = t263 * qJD(5) - t350 * qJDD(4) + t301 * t261;
t270 = qJD(5) + t271;
t241 = t270 * mrSges(6,1) - t263 * mrSges(6,3);
t258 = qJDD(5) - t260;
t262 = -t350 * qJD(4) + t301 * t272;
t233 = t262 * pkin(5) - t263 * qJ(6);
t269 = t270 ^ 2;
t190 = -t269 * pkin(5) + t258 * qJ(6) + 0.2e1 * qJD(6) * t270 - t262 * t233 + t195;
t242 = -t270 * mrSges(7,1) + t263 * mrSges(7,2);
t331 = m(7) * t190 + t258 * mrSges(7,3) + t270 * t242;
t234 = t262 * mrSges(7,1) - t263 * mrSges(7,3);
t338 = -t262 * mrSges(6,1) - t263 * mrSges(6,2) - t234;
t347 = -mrSges(6,3) - mrSges(7,2);
t180 = m(6) * t195 - t258 * mrSges(6,2) + t347 * t227 - t270 * t241 + t338 * t262 + t331;
t194 = -t301 * t198 + t350 * t203;
t228 = -t262 * qJD(5) + t301 * qJDD(4) + t350 * t261;
t240 = -t270 * mrSges(6,2) - t262 * mrSges(6,3);
t192 = -t258 * pkin(5) - t269 * qJ(6) + t263 * t233 + qJDD(6) - t194;
t239 = -t262 * mrSges(7,2) + t270 * mrSges(7,3);
t327 = -m(7) * t192 + t258 * mrSges(7,1) + t270 * t239;
t182 = m(6) * t194 + t258 * mrSges(6,1) + t347 * t228 + t270 * t240 + t338 * t263 + t327;
t328 = t350 * t180 - t301 * t182;
t168 = m(5) * t201 - qJDD(4) * mrSges(5,2) + t260 * mrSges(5,3) - qJD(4) * t267 - t271 * t254 + t328;
t200 = t304 * t208 - t302 * t209;
t266 = -qJD(4) * mrSges(5,2) - t271 * mrSges(5,3);
t197 = -qJDD(4) * pkin(4) - t306 * pkin(9) + t272 * t259 - t200;
t193 = -0.2e1 * qJD(6) * t263 + (t262 * t270 - t228) * qJ(6) + (t263 * t270 + t227) * pkin(5) + t197;
t187 = m(7) * t193 + t227 * mrSges(7,1) - t228 * mrSges(7,3) + t262 * t239 - t263 * t242;
t310 = -m(6) * t197 - t227 * mrSges(6,1) - t228 * mrSges(6,2) - t262 * t240 - t263 * t241 - t187;
t177 = m(5) * t200 + qJDD(4) * mrSges(5,1) - t261 * mrSges(5,3) + qJD(4) * t266 - t272 * t254 + t310;
t163 = t302 * t168 + t304 * t177;
t211 = -t295 * t243 + t337;
t346 = mrSges(4,2) * t295;
t318 = mrSges(4,3) * qJDD(2) + t307 * (-mrSges(4,1) * t298 + t346);
t161 = m(4) * t211 - t318 * t295 + t163;
t329 = t304 * t168 - t302 * t177;
t162 = m(4) * t212 + t318 * t298 + t329;
t155 = t298 * t161 + t295 * t162;
t323 = Ifges(4,5) * t295 + Ifges(4,6) * t298;
t217 = Ifges(7,1) * t263 + Ifges(7,4) * t270 + Ifges(7,5) * t262;
t218 = Ifges(6,1) * t263 - Ifges(6,4) * t262 + Ifges(6,5) * t270;
t326 = -mrSges(7,1) * t193 + mrSges(7,2) * t190;
t215 = Ifges(7,4) * t263 + Ifges(7,2) * t270 + Ifges(7,6) * t262;
t340 = -Ifges(6,5) * t263 + Ifges(6,6) * t262 - Ifges(6,3) * t270 - t215;
t170 = -mrSges(6,1) * t197 + mrSges(6,3) * t195 - pkin(5) * t187 + (t217 + t218) * t270 + t340 * t263 + (Ifges(6,6) - Ifges(7,6)) * t258 + (Ifges(6,4) - Ifges(7,5)) * t228 + (-Ifges(6,2) - Ifges(7,3)) * t227 + t326;
t216 = Ifges(6,4) * t263 - Ifges(6,2) * t262 + Ifges(6,6) * t270;
t213 = Ifges(7,5) * t263 + Ifges(7,6) * t270 + Ifges(7,3) * t262;
t316 = mrSges(7,2) * t192 - mrSges(7,3) * t193 + Ifges(7,1) * t228 + Ifges(7,4) * t258 + Ifges(7,5) * t227 + t270 * t213;
t172 = mrSges(6,2) * t197 - mrSges(6,3) * t194 + Ifges(6,1) * t228 - Ifges(6,4) * t227 + Ifges(6,5) * t258 - qJ(6) * t187 - t270 * t216 + t340 * t262 + t316;
t245 = Ifges(5,4) * t272 - Ifges(5,2) * t271 + Ifges(5,6) * qJD(4);
t246 = Ifges(5,1) * t272 - Ifges(5,4) * t271 + Ifges(5,5) * qJD(4);
t312 = -mrSges(5,1) * t200 + mrSges(5,2) * t201 - Ifges(5,5) * t261 - Ifges(5,6) * t260 - Ifges(5,3) * qJDD(4) - pkin(4) * t310 - pkin(9) * t328 - t350 * t170 - t301 * t172 - t272 * t245 - t271 * t246;
t324 = Ifges(4,4) * t295 + Ifges(4,2) * t298;
t325 = Ifges(4,1) * t295 + Ifges(4,4) * t298;
t351 = -mrSges(4,1) * t211 + mrSges(4,2) * t212 - pkin(3) * t163 - (t295 * t324 - t298 * t325) * t307 + t312;
t141 = (Ifges(3,6) - t323) * qJDD(2) + t307 * Ifges(3,5) - mrSges(3,1) * t268 + mrSges(3,3) * t248 - pkin(2) * t155 + t351;
t330 = -t295 * t161 + t298 * t162;
t153 = m(3) * t248 - t307 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t330;
t238 = -qJDD(2) * pkin(2) - t307 * qJ(3) + t314;
t174 = t301 * t180 + t350 * t182;
t313 = m(5) * t229 - t260 * mrSges(5,1) + t261 * mrSges(5,2) + t271 * t266 + t272 * t267 + t174;
t311 = -m(4) * t238 + mrSges(4,1) * t332 - t313 + (t291 * t307 + t343) * mrSges(4,3);
t165 = -t307 * mrSges(3,2) + m(3) * t247 + t311 + (mrSges(3,1) - t346) * qJDD(2);
t150 = t305 * t153 - t303 * t165;
t353 = pkin(7) * t150 + t141 * t305;
t317 = mrSges(7,1) * t192 - mrSges(7,3) * t190 - Ifges(7,4) * t228 - Ifges(7,2) * t258 - Ifges(7,6) * t227 - t262 * t217;
t352 = -(-t216 + t213) * t263 + mrSges(6,1) * t194 - mrSges(6,2) * t195 + Ifges(6,5) * t228 - Ifges(6,6) * t227 + Ifges(6,3) * t258 + pkin(5) * (-t228 * mrSges(7,2) - t263 * t234 + t327) + qJ(6) * (-t227 * mrSges(7,2) - t262 * t234 + t331) + t262 * t218 - t317;
t344 = t165 * t305;
t336 = t307 * t323;
t154 = m(3) * t268 + t155;
t146 = t153 * t341 - t297 * t154 + t300 * t344;
t244 = Ifges(5,5) * t272 - Ifges(5,6) * t271 + Ifges(5,3) * qJD(4);
t156 = mrSges(5,2) * t229 - mrSges(5,3) * t200 + Ifges(5,1) * t261 + Ifges(5,4) * t260 + Ifges(5,5) * qJDD(4) - pkin(9) * t174 - qJD(4) * t245 - t301 * t170 + t350 * t172 - t271 * t244;
t157 = -mrSges(5,1) * t229 + mrSges(5,3) * t201 + Ifges(5,4) * t261 + Ifges(5,2) * t260 + Ifges(5,6) * qJDD(4) - pkin(4) * t174 + qJD(4) * t246 - t272 * t244 - t352;
t142 = -mrSges(4,1) * t238 + mrSges(4,3) * t212 - pkin(3) * t313 + pkin(8) * t329 + t324 * qJDD(2) + t302 * t156 + t304 * t157 - t295 * t336;
t147 = mrSges(4,2) * t238 - mrSges(4,3) * t211 - pkin(8) * t163 + t325 * qJDD(2) + t304 * t156 - t302 * t157 + t298 * t336;
t137 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t247 - mrSges(3,2) * t248 + t295 * t147 + t298 * t142 + pkin(2) * (-qJDD(2) * t346 + t311) + qJ(3) * t330;
t139 = mrSges(3,2) * t268 - mrSges(3,3) * t247 + Ifges(3,5) * qJDD(2) - t307 * Ifges(3,6) - qJ(3) * t155 - t295 * t142 + t298 * t147;
t315 = mrSges(2,1) * t281 - mrSges(2,2) * t282 + pkin(1) * t146 + t300 * t137 + t139 * t342 + t353 * t297;
t148 = m(2) * t282 + t150;
t145 = t300 * t154 + (t153 * t303 + t344) * t297;
t143 = m(2) * t281 + t146;
t135 = mrSges(2,2) * t294 - mrSges(2,3) * t281 + t305 * t139 - t303 * t141 + (-t145 * t297 - t146 * t300) * pkin(7);
t134 = -mrSges(2,1) * t294 + mrSges(2,3) * t282 - pkin(1) * t145 - t297 * t137 + (t139 * t303 + t353) * t300;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t299 * t135 - t296 * t134 - qJ(1) * (t299 * t143 + t296 * t148), t135, t139, t147, t156, t172, -t262 * t215 + t316; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t296 * t135 + t299 * t134 + qJ(1) * (-t296 * t143 + t299 * t148), t134, t141, t142, t157, t170, -t263 * t213 - t317; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t315, t315, t137, t323 * qJDD(2) - t351, -t312, t352, Ifges(7,5) * t228 + Ifges(7,6) * t258 + Ifges(7,3) * t227 + t263 * t215 - t270 * t217 - t326;];
m_new  = t1;
