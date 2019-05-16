% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-05-06 01:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRP8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:46:21
% EndTime: 2019-05-06 01:46:35
% DurationCPUTime: 6.21s
% Computational Cost: add. (112632->339), mult. (220128->402), div. (0->0), fcn. (143802->8), ass. (0->127)
t298 = sin(qJ(1));
t301 = cos(qJ(1));
t277 = -t301 * g(1) - t298 * g(2);
t317 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t277;
t296 = sin(qJ(4));
t297 = sin(qJ(3));
t299 = cos(qJ(4));
t300 = cos(qJ(3));
t261 = (t296 * t300 + t297 * t299) * qJD(1);
t327 = qJD(1) * qJD(3);
t270 = -t297 * qJDD(1) - t300 * t327;
t328 = qJD(1) * t300;
t275 = (qJD(3) * pkin(3)) - pkin(8) * t328;
t292 = t297 ^ 2;
t302 = qJD(1) ^ 2;
t338 = -pkin(1) - pkin(7);
t217 = -t270 * pkin(3) + t275 * t328 + (-pkin(8) * t292 + t338) * t302 + t317;
t262 = (-t296 * t297 + t299 * t300) * qJD(1);
t324 = t297 * t327;
t271 = t300 * qJDD(1) - t324;
t227 = -t262 * qJD(4) + t299 * t270 - t296 * t271;
t228 = -t261 * qJD(4) + t296 * t270 + t299 * t271;
t286 = qJD(3) + qJD(4);
t180 = (t261 * t286 - t228) * pkin(9) + (t262 * t286 - t227) * pkin(4) + t217;
t276 = t298 * g(1) - t301 * g(2);
t316 = -t302 * qJ(2) + qJDD(2) - t276;
t250 = t338 * qJDD(1) + t316;
t240 = t297 * g(3) + t300 * t250;
t210 = (-t271 - t324) * pkin(8) + (-t297 * t300 * t302 + qJDD(3)) * pkin(3) + t240;
t241 = -t300 * g(3) + t297 * t250;
t211 = -t292 * t302 * pkin(3) + t270 * pkin(8) - qJD(3) * t275 + t241;
t186 = t296 * t210 + t299 * t211;
t238 = t261 * pkin(4) - t262 * pkin(9);
t284 = t286 ^ 2;
t285 = qJDD(3) + qJDD(4);
t183 = -t284 * pkin(4) + t285 * pkin(9) - t261 * t238 + t186;
t295 = sin(qJ(5));
t337 = cos(qJ(5));
t177 = t337 * t180 - t295 * t183;
t178 = t295 * t180 + t337 * t183;
t243 = t337 * t262 + t295 * t286;
t195 = t243 * qJD(5) + t295 * t228 - t337 * t285;
t242 = t295 * t262 - t337 * t286;
t196 = -t242 * qJD(5) + t337 * t228 + t295 * t285;
t257 = qJD(5) + t261;
t197 = Ifges(7,5) * t243 + Ifges(7,6) * t257 + Ifges(7,3) * t242;
t200 = Ifges(6,4) * t243 - Ifges(6,2) * t242 + Ifges(6,6) * t257;
t202 = Ifges(6,1) * t243 - Ifges(6,4) * t242 + Ifges(6,5) * t257;
t213 = t242 * mrSges(7,1) - t243 * mrSges(7,3);
t226 = qJDD(5) - t227;
t212 = t242 * pkin(5) - t243 * qJ(6);
t255 = t257 ^ 2;
t173 = -t255 * pkin(5) + t226 * qJ(6) + 0.2e1 * qJD(6) * t257 - t242 * t212 + t178;
t175 = -t226 * pkin(5) - t255 * qJ(6) + t243 * t212 + qJDD(6) - t177;
t201 = Ifges(7,1) * t243 + Ifges(7,4) * t257 + Ifges(7,5) * t242;
t315 = mrSges(7,1) * t175 - mrSges(7,3) * t173 - Ifges(7,4) * t196 - Ifges(7,2) * t226 - Ifges(7,6) * t195 - t242 * t201;
t229 = -t242 * mrSges(7,2) + t257 * mrSges(7,3);
t321 = -m(7) * t175 + t226 * mrSges(7,1) + t257 * t229;
t232 = -t257 * mrSges(7,1) + t243 * mrSges(7,2);
t325 = m(7) * t173 + t226 * mrSges(7,3) + t257 * t232;
t339 = -(-t200 + t197) * t243 + mrSges(6,1) * t177 - mrSges(6,2) * t178 + Ifges(6,5) * t196 - Ifges(6,6) * t195 + Ifges(6,3) * t226 + pkin(5) * (-t196 * mrSges(7,2) - t243 * t213 + t321) + qJ(6) * (-t195 * mrSges(7,2) - t242 * t213 + t325) + t242 * t202 - t315;
t336 = mrSges(2,1) - mrSges(3,2);
t335 = -mrSges(6,3) - mrSges(7,2);
t334 = Ifges(2,5) - Ifges(3,4);
t333 = (-Ifges(2,6) + Ifges(3,5));
t237 = t261 * mrSges(5,1) + t262 * mrSges(5,2);
t248 = t286 * mrSges(5,1) - t262 * mrSges(5,3);
t231 = t257 * mrSges(6,1) - t243 * mrSges(6,3);
t330 = -t242 * mrSges(6,1) - t243 * mrSges(6,2) - t213;
t163 = m(6) * t178 - t226 * mrSges(6,2) + t335 * t195 - t257 * t231 + t330 * t242 + t325;
t230 = -t257 * mrSges(6,2) - t242 * mrSges(6,3);
t165 = m(6) * t177 + t226 * mrSges(6,1) + t335 * t196 + t257 * t230 + t330 * t243 + t321;
t322 = t337 * t163 - t295 * t165;
t151 = m(5) * t186 - t285 * mrSges(5,2) + t227 * mrSges(5,3) - t261 * t237 - t286 * t248 + t322;
t185 = t299 * t210 - t296 * t211;
t247 = -t286 * mrSges(5,2) - t261 * mrSges(5,3);
t182 = -t285 * pkin(4) - t284 * pkin(9) + t262 * t238 - t185;
t176 = -0.2e1 * qJD(6) * t243 + (t242 * t257 - t196) * qJ(6) + (t243 * t257 + t195) * pkin(5) + t182;
t170 = m(7) * t176 + t195 * mrSges(7,1) - t196 * mrSges(7,3) + t242 * t229 - t243 * t232;
t307 = -m(6) * t182 - t195 * mrSges(6,1) - t196 * mrSges(6,2) - t242 * t230 - t243 * t231 - t170;
t160 = m(5) * t185 + t285 * mrSges(5,1) - t228 * mrSges(5,3) - t262 * t237 + t286 * t247 + t307;
t145 = t296 * t151 + t299 * t160;
t157 = t295 * t163 + t337 * t165;
t199 = Ifges(7,4) * t243 + Ifges(7,2) * t257 + Ifges(7,6) * t242;
t332 = -Ifges(6,5) * t243 + Ifges(6,6) * t242 - Ifges(6,3) * t257 - t199;
t329 = qJD(1) * t297;
t269 = (mrSges(4,1) * t297 + mrSges(4,2) * t300) * qJD(1);
t273 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t329;
t142 = m(4) * t240 + qJDD(3) * mrSges(4,1) - t271 * mrSges(4,3) + qJD(3) * t273 - t269 * t328 + t145;
t274 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t328;
t323 = t299 * t151 - t296 * t160;
t143 = m(4) * t241 - qJDD(3) * mrSges(4,2) + t270 * mrSges(4,3) - qJD(3) * t274 - t269 * t329 + t323;
t139 = -t297 * t142 + t300 * t143;
t320 = -mrSges(7,1) * t176 + mrSges(7,2) * t173;
t138 = t300 * t142 + t297 * t143;
t314 = mrSges(7,2) * t175 - mrSges(7,3) * t176 + Ifges(7,1) * t196 + Ifges(7,4) * t226 + Ifges(7,5) * t195 + t257 * t197;
t256 = -qJDD(1) * pkin(1) + t316;
t313 = -m(3) * t256 + (t302 * mrSges(3,3)) - t138;
t312 = m(5) * t217 - t227 * mrSges(5,1) + t228 * mrSges(5,2) + t261 * t247 + t262 * t248 + t157;
t153 = -mrSges(6,1) * t182 + mrSges(6,3) * t178 - pkin(5) * t170 + (t201 + t202) * t257 + t332 * t243 + (Ifges(6,6) - Ifges(7,6)) * t226 + (Ifges(6,4) - Ifges(7,5)) * t196 + (-Ifges(6,2) - Ifges(7,3)) * t195 + t320;
t155 = mrSges(6,2) * t182 - mrSges(6,3) * t177 + Ifges(6,1) * t196 - Ifges(6,4) * t195 + Ifges(6,5) * t226 - qJ(6) * t170 - t257 * t200 + t332 * t242 + t314;
t233 = Ifges(5,5) * t262 - Ifges(5,6) * t261 + Ifges(5,3) * t286;
t234 = Ifges(5,4) * t262 - Ifges(5,2) * t261 + Ifges(5,6) * t286;
t135 = mrSges(5,2) * t217 - mrSges(5,3) * t185 + Ifges(5,1) * t228 + Ifges(5,4) * t227 + Ifges(5,5) * t285 - pkin(9) * t157 - t295 * t153 + t337 * t155 - t261 * t233 - t286 * t234;
t235 = Ifges(5,1) * t262 - Ifges(5,4) * t261 + Ifges(5,5) * t286;
t140 = -mrSges(5,1) * t217 + mrSges(5,3) * t186 + Ifges(5,4) * t228 + Ifges(5,2) * t227 + Ifges(5,6) * t285 - pkin(4) * t157 - t262 * t233 + t286 * t235 - t339;
t249 = t338 * t302 + t317;
t258 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t300 - Ifges(4,6) * t297) * qJD(1);
t260 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t300 - Ifges(4,4) * t297) * qJD(1);
t131 = -mrSges(4,1) * t249 + mrSges(4,3) * t241 + Ifges(4,4) * t271 + Ifges(4,2) * t270 + Ifges(4,6) * qJDD(3) - pkin(3) * t312 + pkin(8) * t323 + qJD(3) * t260 + t296 * t135 + t299 * t140 - t258 * t328;
t259 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t300 - Ifges(4,2) * t297) * qJD(1);
t133 = mrSges(4,2) * t249 - mrSges(4,3) * t240 + Ifges(4,1) * t271 + Ifges(4,4) * t270 + Ifges(4,5) * qJDD(3) - pkin(8) * t145 - qJD(3) * t259 + t299 * t135 - t296 * t140 - t258 * t329;
t253 = t302 * pkin(1) - t317;
t311 = mrSges(3,2) * t256 - mrSges(3,3) * t253 + Ifges(3,1) * qJDD(1) - pkin(7) * t138 - t297 * t131 + t300 * t133;
t148 = -m(4) * t249 + t270 * mrSges(4,1) - t271 * mrSges(4,2) - t273 * t329 - t274 * t328 - t312;
t310 = -mrSges(3,1) * t253 - pkin(2) * t148 - pkin(7) * t139 - t300 * t131 - t297 * t133;
t309 = -mrSges(5,1) * t185 + mrSges(5,2) * t186 - Ifges(5,5) * t228 - Ifges(5,6) * t227 - Ifges(5,3) * t285 - pkin(4) * t307 - pkin(9) * t322 - t337 * t153 - t295 * t155 - t262 * t234 - t261 * t235;
t306 = -m(3) * t253 + t302 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t148;
t308 = -mrSges(2,2) * t277 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t313) + qJ(2) * t306 + mrSges(2,1) * t276 + Ifges(2,3) * qJDD(1) + t311;
t305 = -mrSges(4,1) * t240 + mrSges(4,2) * t241 - Ifges(4,5) * t271 - Ifges(4,6) * t270 - Ifges(4,3) * qJDD(3) - pkin(3) * t145 - t259 * t328 - t260 * t329 + t309;
t303 = -mrSges(3,1) * t256 - pkin(2) * t138 + t305;
t146 = m(2) * t277 - t302 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t306;
t137 = -m(3) * g(3) + t139;
t134 = m(2) * t276 - t302 * mrSges(2,2) + t336 * qJDD(1) + t313;
t130 = -t303 - qJ(2) * t137 + (t333 * t302) + t334 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t276;
t129 = mrSges(2,3) * t277 - pkin(1) * t137 + t336 * g(3) - t333 * qJDD(1) + t334 * t302 + t310;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t301 * t130 - t298 * t129 - pkin(6) * (t301 * t134 + t298 * t146), t130, t311, t133, t135, t155, -t242 * t199 + t314; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t298 * t130 + t301 * t129 + pkin(6) * (-t298 * t134 + t301 * t146), t129, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t302 * Ifges(3,5)) + t303, t131, t140, t153, -t243 * t197 - t315; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t308, t308, mrSges(3,2) * g(3) + t302 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t310, -t305, -t309, t339, Ifges(7,5) * t196 + Ifges(7,6) * t226 + Ifges(7,3) * t195 + t243 * t199 - t257 * t201 - t320;];
m_new  = t1;
