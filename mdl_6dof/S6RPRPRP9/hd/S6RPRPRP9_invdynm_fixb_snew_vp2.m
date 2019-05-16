% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-05-05 18:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRP9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:04:53
% EndTime: 2019-05-05 18:05:07
% DurationCPUTime: 6.32s
% Computational Cost: add. (104342->337), mult. (216732->399), div. (0->0), fcn. (136556->8), ass. (0->126)
t295 = sin(qJ(1));
t297 = cos(qJ(1));
t279 = -t297 * g(1) - t295 * g(2);
t334 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t279;
t333 = (-pkin(1) - pkin(7));
t332 = cos(qJ(5));
t331 = mrSges(2,1) - mrSges(3,2);
t330 = -mrSges(6,3) - mrSges(7,2);
t329 = -Ifges(3,4) + Ifges(2,5);
t328 = (Ifges(3,5) - Ifges(2,6));
t299 = qJD(1) ^ 2;
t238 = (t333 * t299) - t334;
t294 = sin(qJ(3));
t296 = cos(qJ(3));
t323 = qJD(1) * qJD(3);
t319 = t296 * t323;
t273 = t294 * qJDD(1) + t319;
t320 = t294 * t323;
t274 = t296 * qJDD(1) - t320;
t213 = (-t274 + t320) * qJ(4) + (t273 + t319) * pkin(3) + t238;
t278 = t295 * g(1) - t297 * g(2);
t313 = -t299 * qJ(2) + qJDD(2) - t278;
t244 = t333 * qJDD(1) + t313;
t233 = -t296 * g(3) + t294 * t244;
t271 = (pkin(3) * t294 - qJ(4) * t296) * qJD(1);
t298 = qJD(3) ^ 2;
t324 = t294 * qJD(1);
t218 = -t298 * pkin(3) + qJDD(3) * qJ(4) - t271 * t324 + t233;
t291 = sin(pkin(9));
t292 = cos(pkin(9));
t325 = qJD(1) * t296;
t266 = t291 * qJD(3) + t292 * t325;
t182 = -0.2e1 * qJD(4) * t266 + t292 * t213 - t291 * t218;
t242 = t291 * qJDD(3) + t292 * t274;
t265 = t292 * qJD(3) - t291 * t325;
t179 = (t265 * t324 - t242) * pkin(8) + (t265 * t266 + t273) * pkin(4) + t182;
t183 = 0.2e1 * qJD(4) * t265 + t291 * t213 + t292 * t218;
t241 = t292 * qJDD(3) - t291 * t274;
t243 = pkin(4) * t324 - t266 * pkin(8);
t264 = t265 ^ 2;
t181 = -t264 * pkin(4) + t241 * pkin(8) - t243 * t324 + t183;
t293 = sin(qJ(5));
t175 = t293 * t179 + t332 * t181;
t229 = t293 * t265 + t332 * t266;
t197 = t229 * qJD(5) - t332 * t241 + t293 * t242;
t281 = qJD(5) + t324;
t221 = t281 * mrSges(6,1) - t229 * mrSges(6,3);
t228 = -t332 * t265 + t293 * t266;
t270 = qJDD(5) + t273;
t208 = t228 * pkin(5) - t229 * qJ(6);
t280 = t281 ^ 2;
t170 = -t280 * pkin(5) + t270 * qJ(6) + 0.2e1 * qJD(6) * t281 - t228 * t208 + t175;
t222 = -t281 * mrSges(7,1) + t229 * mrSges(7,2);
t321 = m(7) * t170 + t270 * mrSges(7,3) + t281 * t222;
t209 = t228 * mrSges(7,1) - t229 * mrSges(7,3);
t326 = -t228 * mrSges(6,1) - t229 * mrSges(6,2) - t209;
t157 = m(6) * t175 - t270 * mrSges(6,2) + t330 * t197 - t281 * t221 + t326 * t228 + t321;
t174 = t332 * t179 - t293 * t181;
t198 = -t228 * qJD(5) + t293 * t241 + t332 * t242;
t219 = -t281 * mrSges(6,2) - t228 * mrSges(6,3);
t172 = -t270 * pkin(5) - t280 * qJ(6) + t229 * t208 + qJDD(6) - t174;
t220 = -t228 * mrSges(7,2) + t281 * mrSges(7,3);
t316 = -m(7) * t172 + t270 * mrSges(7,1) + t281 * t220;
t159 = m(6) * t174 + t270 * mrSges(6,1) + t330 * t198 + t281 * t219 + t326 * t229 + t316;
t154 = t293 * t157 + t332 * t159;
t230 = -t265 * mrSges(5,1) + t266 * mrSges(5,2);
t239 = -mrSges(5,2) * t324 + t265 * mrSges(5,3);
t150 = m(5) * t182 + t273 * mrSges(5,1) - t242 * mrSges(5,3) - t266 * t230 + t239 * t324 + t154;
t240 = mrSges(5,1) * t324 - t266 * mrSges(5,3);
t317 = t332 * t157 - t293 * t159;
t151 = m(5) * t183 - t273 * mrSges(5,2) + t241 * mrSges(5,3) + t265 * t230 - t240 * t324 + t317;
t145 = t292 * t150 + t291 * t151;
t201 = Ifges(7,4) * t229 + Ifges(7,2) * t281 + Ifges(7,6) * t228;
t327 = -Ifges(6,5) * t229 + Ifges(6,6) * t228 - Ifges(6,3) * t281 - t201;
t272 = (mrSges(4,1) * t294 + mrSges(4,2) * t296) * qJD(1);
t277 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t325;
t318 = -t291 * t150 + t292 * t151;
t143 = m(4) * t233 - qJDD(3) * mrSges(4,2) - t273 * mrSges(4,3) - qJD(3) * t277 - t272 * t324 + t318;
t232 = t294 * g(3) + t296 * t244;
t276 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t324;
t215 = -qJDD(3) * pkin(3) - t298 * qJ(4) + t271 * t325 + qJDD(4) - t232;
t184 = -t241 * pkin(4) - t264 * pkin(8) + t266 * t243 + t215;
t177 = -0.2e1 * qJD(6) * t229 + (t228 * t281 - t198) * qJ(6) + (t229 * t281 + t197) * pkin(5) + t184;
t167 = m(7) * t177 + t197 * mrSges(7,1) - t198 * mrSges(7,3) + t228 * t220 - t229 * t222;
t305 = m(6) * t184 + t197 * mrSges(6,1) + t198 * mrSges(6,2) + t228 * t219 + t229 * t221 + t167;
t301 = -m(5) * t215 + t241 * mrSges(5,1) - t242 * mrSges(5,2) + t265 * t239 - t266 * t240 - t305;
t160 = m(4) * t232 + qJDD(3) * mrSges(4,1) - t274 * mrSges(4,3) + qJD(3) * t276 - t272 * t325 + t301;
t138 = t296 * t143 - t294 * t160;
t315 = -mrSges(7,1) * t177 + mrSges(7,2) * t170;
t137 = t294 * t143 + t296 * t160;
t199 = Ifges(7,5) * t229 + Ifges(7,6) * t281 + Ifges(7,3) * t228;
t312 = mrSges(7,2) * t172 - mrSges(7,3) * t177 + Ifges(7,1) * t198 + Ifges(7,4) * t270 + Ifges(7,5) * t197 + t281 * t199;
t249 = -qJDD(1) * pkin(1) + t313;
t311 = -m(3) * t249 + (t299 * mrSges(3,3)) - t137;
t141 = -m(4) * t238 - t273 * mrSges(4,1) - t274 * mrSges(4,2) - t276 * t324 - t277 * t325 - t145;
t203 = Ifges(7,1) * t229 + Ifges(7,4) * t281 + Ifges(7,5) * t228;
t310 = mrSges(7,1) * t172 - mrSges(7,3) * t170 - Ifges(7,4) * t198 - Ifges(7,2) * t270 - Ifges(7,6) * t197 + t229 * t199 - t228 * t203;
t204 = Ifges(6,1) * t229 - Ifges(6,4) * t228 + Ifges(6,5) * t281;
t152 = -mrSges(6,1) * t184 + mrSges(6,3) * t175 - pkin(5) * t167 + (t203 + t204) * t281 + (Ifges(6,6) - Ifges(7,6)) * t270 + t327 * t229 + (Ifges(6,4) - Ifges(7,5)) * t198 + (-Ifges(6,2) - Ifges(7,3)) * t197 + t315;
t202 = Ifges(6,4) * t229 - Ifges(6,2) * t228 + Ifges(6,6) * t281;
t153 = mrSges(6,2) * t184 - mrSges(6,3) * t174 + Ifges(6,1) * t198 - Ifges(6,4) * t197 + Ifges(6,5) * t270 - qJ(6) * t167 - t281 * t202 + t327 * t228 + t312;
t223 = Ifges(5,5) * t266 + Ifges(5,6) * t265 + Ifges(5,3) * t324;
t225 = Ifges(5,1) * t266 + Ifges(5,4) * t265 + Ifges(5,5) * t324;
t131 = -mrSges(5,1) * t215 + mrSges(5,3) * t183 + Ifges(5,4) * t242 + Ifges(5,2) * t241 + Ifges(5,6) * t273 - pkin(4) * t305 + pkin(8) * t317 + t332 * t152 + t293 * t153 - t266 * t223 + t225 * t324;
t224 = Ifges(5,4) * t266 + Ifges(5,2) * t265 + Ifges(5,6) * t324;
t133 = mrSges(5,2) * t215 - mrSges(5,3) * t182 + Ifges(5,1) * t242 + Ifges(5,4) * t241 + Ifges(5,5) * t273 - pkin(8) * t154 - t293 * t152 + t332 * t153 + t265 * t223 - t224 * t324;
t254 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t296 - Ifges(4,6) * t294) * qJD(1);
t255 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t296 - Ifges(4,2) * t294) * qJD(1);
t128 = mrSges(4,2) * t238 - mrSges(4,3) * t232 + Ifges(4,1) * t274 - Ifges(4,4) * t273 + Ifges(4,5) * qJDD(3) - qJ(4) * t145 - qJD(3) * t255 - t131 * t291 + t133 * t292 - t254 * t324;
t256 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t296 - Ifges(4,4) * t294) * qJD(1);
t302 = mrSges(6,2) * t175 - t228 * t204 - qJ(6) * (-t197 * mrSges(7,2) - t228 * t209 + t321) - pkin(5) * (-t198 * mrSges(7,2) - t229 * t209 + t316) - mrSges(6,1) * t174 - t229 * t202 + Ifges(6,6) * t197 - Ifges(6,5) * t198 - Ifges(6,3) * t270 + t310;
t300 = mrSges(5,1) * t182 - mrSges(5,2) * t183 + Ifges(5,5) * t242 + Ifges(5,6) * t241 + pkin(4) * t154 + t266 * t224 - t265 * t225 - t302;
t129 = Ifges(4,6) * qJDD(3) + (-Ifges(5,3) - Ifges(4,2)) * t273 + Ifges(4,4) * t274 + qJD(3) * t256 + mrSges(4,3) * t233 - mrSges(4,1) * t238 - pkin(3) * t145 - t254 * t325 - t300;
t247 = t299 * pkin(1) + t334;
t309 = mrSges(3,2) * t249 - mrSges(3,3) * t247 + Ifges(3,1) * qJDD(1) - pkin(7) * t137 + t296 * t128 - t129 * t294;
t308 = -mrSges(3,1) * t247 - pkin(2) * t141 - pkin(7) * t138 - t294 * t128 - t296 * t129;
t307 = mrSges(4,1) * t232 - mrSges(4,2) * t233 + Ifges(4,5) * t274 - Ifges(4,6) * t273 + Ifges(4,3) * qJDD(3) + pkin(3) * t301 + qJ(4) * t318 + t292 * t131 + t291 * t133 + t255 * t325 + t256 * t324;
t306 = -m(3) * t247 + t299 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t141;
t304 = -mrSges(2,2) * t279 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t311) + qJ(2) * t306 + mrSges(2,1) * t278 + Ifges(2,3) * qJDD(1) + t309;
t303 = mrSges(3,1) * t249 + pkin(2) * t137 + t307;
t139 = m(2) * t279 - t299 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t306;
t136 = -m(3) * g(3) + t138;
t134 = m(2) * t278 - t299 * mrSges(2,2) + t331 * qJDD(1) + t311;
t126 = -mrSges(2,3) * t278 - qJ(2) * t136 + (t328 * t299) + t329 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t303;
t125 = mrSges(2,3) * t279 - pkin(1) * t136 + t331 * g(3) - t328 * qJDD(1) + t329 * t299 + t308;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t297 * t126 - t295 * t125 - pkin(6) * (t134 * t297 + t139 * t295), t126, t309, t128, t133, t153, -t201 * t228 + t312; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t295 * t126 + t297 * t125 + pkin(6) * (-t134 * t295 + t139 * t297), t125, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t299 * Ifges(3,5)) - t303, t129, t131, t152, -t310; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t304, t304, mrSges(3,2) * g(3) + t299 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t308, t307, Ifges(5,3) * t273 + t300, -t302, Ifges(7,5) * t198 + Ifges(7,6) * t270 + Ifges(7,3) * t197 + t229 * t201 - t281 * t203 - t315;];
m_new  = t1;
