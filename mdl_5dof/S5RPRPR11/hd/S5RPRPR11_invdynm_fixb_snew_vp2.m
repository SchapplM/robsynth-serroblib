% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR11_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR11_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR11_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:15
% EndTime: 2019-12-31 18:27:23
% DurationCPUTime: 4.68s
% Computational Cost: add. (47525->290), mult. (115174->350), div. (0->0), fcn. (78532->8), ass. (0->120)
t273 = sin(pkin(8));
t263 = t273 ^ 2;
t274 = cos(pkin(8));
t264 = t274 ^ 2;
t281 = qJD(1) ^ 2;
t277 = sin(qJ(1));
t279 = cos(qJ(1));
t249 = t277 * g(1) - t279 * g(2);
t303 = -qJDD(2) + t249;
t318 = pkin(2) * t274;
t226 = -(qJ(2) + (t263 + t264) * pkin(6)) * t281 - (pkin(1) + t318) * qJDD(1) - t303;
t276 = sin(qJ(3));
t319 = cos(qJ(3));
t293 = t319 * t273 + t274 * t276;
t304 = t274 * t319;
t311 = qJD(1) * t273;
t241 = -qJD(1) * t304 + t276 * t311;
t310 = t241 * qJD(3);
t228 = t293 * qJDD(1) - t310;
t326 = t226 + (-t228 + t310) * qJ(4);
t324 = qJD(1) * t274;
t250 = -t279 * g(1) - t277 * g(2);
t243 = -t281 * pkin(1) + qJDD(1) * qJ(2) + t250;
t308 = qJD(1) * qJD(2);
t302 = -t274 * g(3) - 0.2e1 * t273 * t308;
t200 = (-pkin(6) * qJDD(1) + t281 * t318 - t243) * t273 + t302;
t230 = -t273 * g(3) + (t243 + 0.2e1 * t308) * t274;
t306 = qJDD(1) * t274;
t315 = t264 * t281;
t201 = -pkin(2) * t315 + pkin(6) * t306 + t230;
t184 = t276 * t200 + t319 * t201;
t307 = qJDD(1) * t273;
t242 = t293 * qJD(1);
t309 = t242 * qJD(3);
t227 = -qJDD(1) * t304 + t276 * t307 + t309;
t234 = qJD(3) * mrSges(4,1) - t242 * mrSges(4,3);
t183 = t319 * t200 - t276 * t201;
t214 = t241 * pkin(3) - t242 * qJ(4);
t280 = qJD(3) ^ 2;
t173 = -qJDD(3) * pkin(3) - t280 * qJ(4) + t242 * t214 + qJDD(4) - t183;
t165 = (-t228 - t310) * pkin(7) + (t241 * t242 - qJDD(3)) * pkin(4) + t173;
t320 = 2 * qJD(4);
t171 = -t280 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t320 - t241 * t214 + t184;
t237 = -qJD(3) * pkin(4) - t242 * pkin(7);
t240 = t241 ^ 2;
t166 = -t240 * pkin(4) + t227 * pkin(7) + qJD(3) * t237 + t171;
t275 = sin(qJ(5));
t278 = cos(qJ(5));
t161 = t278 * t165 - t275 * t166;
t209 = t278 * t241 - t275 * t242;
t181 = t209 * qJD(5) + t275 * t227 + t278 * t228;
t210 = t275 * t241 + t278 * t242;
t190 = -t209 * mrSges(6,1) + t210 * mrSges(6,2);
t265 = -qJD(3) + qJD(5);
t196 = -t265 * mrSges(6,2) + t209 * mrSges(6,3);
t262 = -qJDD(3) + qJDD(5);
t156 = m(6) * t161 + t262 * mrSges(6,1) - t181 * mrSges(6,3) - t210 * t190 + t265 * t196;
t162 = t275 * t165 + t278 * t166;
t180 = -t210 * qJD(5) + t278 * t227 - t275 * t228;
t197 = t265 * mrSges(6,1) - t210 * mrSges(6,3);
t157 = m(6) * t162 - t262 * mrSges(6,2) + t180 * mrSges(6,3) + t209 * t190 - t265 * t197;
t148 = -t275 * t156 + t278 * t157;
t235 = -qJD(3) * mrSges(5,1) + t242 * mrSges(5,2);
t292 = m(5) * t171 + qJDD(3) * mrSges(5,3) + qJD(3) * t235 + t148;
t215 = t241 * mrSges(5,1) - t242 * mrSges(5,3);
t313 = -t241 * mrSges(4,1) - t242 * mrSges(4,2) - t215;
t317 = -mrSges(4,3) - mrSges(5,2);
t141 = m(4) * t184 - qJDD(3) * mrSges(4,2) - qJD(3) * t234 + t317 * t227 + t313 * t241 + t292;
t233 = -qJD(3) * mrSges(4,2) - t241 * mrSges(4,3);
t147 = t278 * t156 + t275 * t157;
t236 = -t241 * mrSges(5,2) + qJD(3) * mrSges(5,3);
t289 = -m(5) * t173 + qJDD(3) * mrSges(5,1) + qJD(3) * t236 - t147;
t142 = m(4) * t183 + qJDD(3) * mrSges(4,1) + qJD(3) * t233 + t317 * t228 + t313 * t242 + t289;
t137 = t276 * t141 + t319 * t142;
t229 = -t273 * t243 + t302;
t206 = Ifges(4,4) * t242 - Ifges(4,2) * t241 + Ifges(4,6) * qJD(3);
t208 = Ifges(4,1) * t242 - Ifges(4,4) * t241 + Ifges(4,5) * qJD(3);
t203 = Ifges(5,5) * t242 + Ifges(5,6) * qJD(3) + Ifges(5,3) * t241;
t207 = Ifges(5,1) * t242 + Ifges(5,4) * qJD(3) + Ifges(5,5) * t241;
t186 = Ifges(6,4) * t210 + Ifges(6,2) * t209 + Ifges(6,6) * t265;
t187 = Ifges(6,1) * t210 + Ifges(6,4) * t209 + Ifges(6,5) * t265;
t291 = mrSges(6,1) * t161 - mrSges(6,2) * t162 + Ifges(6,5) * t181 + Ifges(6,6) * t180 + Ifges(6,3) * t262 + t210 * t186 - t209 * t187;
t285 = mrSges(5,1) * t173 - mrSges(5,3) * t171 - Ifges(5,4) * t228 - Ifges(5,2) * qJDD(3) - Ifges(5,6) * t227 + pkin(4) * t147 + t242 * t203 - t241 * t207 + t291;
t283 = mrSges(4,2) * t184 - t241 * t208 - qJ(4) * (-t227 * mrSges(5,2) - t241 * t215 + t292) - pkin(3) * (-t228 * mrSges(5,2) - t242 * t215 + t289) - mrSges(4,1) * t183 - t242 * t206 + Ifges(4,6) * t227 - Ifges(4,5) * t228 - Ifges(4,3) * qJDD(3) + t285;
t297 = Ifges(3,4) * t273 + Ifges(3,2) * t274;
t298 = Ifges(3,1) * t273 + Ifges(3,4) * t274;
t322 = -mrSges(3,1) * t229 + mrSges(3,2) * t230 - pkin(2) * t137 - (t297 * t311 - t298 * t324) * qJD(1) + t283;
t316 = mrSges(3,2) * t273;
t205 = Ifges(5,4) * t242 + Ifges(5,2) * qJD(3) + Ifges(5,6) * t241;
t314 = -Ifges(4,5) * t242 + Ifges(4,6) * t241 - Ifges(4,3) * qJD(3) - t205;
t294 = mrSges(3,3) * qJDD(1) + t281 * (-mrSges(3,1) * t274 + t316);
t135 = m(3) * t229 - t294 * t273 + t137;
t300 = t319 * t141 - t276 * t142;
t136 = m(3) * t230 + t294 * t274 + t300;
t301 = -t273 * t135 + t274 * t136;
t296 = Ifges(3,5) * t273 + Ifges(3,6) * t274;
t164 = -t240 * pkin(7) + (-pkin(3) - pkin(4)) * t227 + (-pkin(3) * qJD(3) + t237 + t320) * t242 - t326;
t158 = -m(6) * t164 + t180 * mrSges(6,1) - t181 * mrSges(6,2) + t209 * t196 - t210 * t197;
t168 = -0.2e1 * qJD(4) * t242 + (t227 + t309) * pkin(3) + t326;
t154 = m(5) * t168 + t227 * mrSges(5,1) - t228 * mrSges(5,3) - t242 * t235 + t241 * t236 + t158;
t185 = Ifges(6,5) * t210 + Ifges(6,6) * t209 + Ifges(6,3) * t265;
t151 = -mrSges(6,1) * t164 + mrSges(6,3) * t162 + Ifges(6,4) * t181 + Ifges(6,2) * t180 + Ifges(6,6) * t262 - t210 * t185 + t265 * t187;
t152 = mrSges(6,2) * t164 - mrSges(6,3) * t161 + Ifges(6,1) * t181 + Ifges(6,4) * t180 + Ifges(6,5) * t262 + t209 * t185 - t265 * t186;
t287 = -mrSges(5,1) * t168 + mrSges(5,2) * t171 - pkin(4) * t158 - pkin(7) * t148 - t278 * t151 - t275 * t152;
t132 = -mrSges(4,1) * t226 + mrSges(4,3) * t184 - pkin(3) * t154 + t314 * t242 + (Ifges(4,4) - Ifges(5,5)) * t228 + (-Ifges(4,2) - Ifges(5,3)) * t227 + (Ifges(4,6) - Ifges(5,6)) * qJDD(3) + (t208 + t207) * qJD(3) + t287;
t288 = mrSges(5,2) * t173 - mrSges(5,3) * t168 + Ifges(5,1) * t228 + Ifges(5,4) * qJDD(3) + Ifges(5,5) * t227 - pkin(7) * t147 + qJD(3) * t203 - t275 * t151 + t278 * t152;
t133 = mrSges(4,2) * t226 - mrSges(4,3) * t183 + Ifges(4,1) * t228 - Ifges(4,4) * t227 + Ifges(4,5) * qJDD(3) - qJ(4) * t154 - qJD(3) * t206 + t314 * t241 + t288;
t239 = -qJDD(1) * pkin(1) - t281 * qJ(2) - t303;
t245 = t296 * qJD(1);
t286 = m(4) * t226 + t227 * mrSges(4,1) + t228 * mrSges(4,2) + t241 * t233 + t242 * t234 + t154;
t126 = -mrSges(3,1) * t239 + mrSges(3,3) * t230 - pkin(2) * t286 + pkin(6) * t300 + t297 * qJDD(1) + t319 * t132 + t276 * t133 - t245 * t311;
t128 = mrSges(3,2) * t239 - mrSges(3,3) * t229 - pkin(6) * t137 + t298 * qJDD(1) - t276 * t132 + t319 * t133 + t245 * t324;
t284 = -m(3) * t239 + mrSges(3,1) * t306 - t286 + (t263 * t281 + t315) * mrSges(3,3);
t290 = -mrSges(2,2) * t250 + qJ(2) * t301 + t274 * t126 + t273 * t128 + pkin(1) * (-mrSges(3,2) * t307 + t284) + mrSges(2,1) * t249 + Ifges(2,3) * qJDD(1);
t149 = (mrSges(2,1) - t316) * qJDD(1) - t281 * mrSges(2,2) + m(2) * t249 + t284;
t131 = t274 * t135 + t273 * t136;
t129 = m(2) * t250 - t281 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t301;
t124 = (Ifges(2,6) - t296) * qJDD(1) + mrSges(2,1) * g(3) + t281 * Ifges(2,5) + mrSges(2,3) * t250 - pkin(1) * t131 + t322;
t123 = -mrSges(2,2) * g(3) - mrSges(2,3) * t249 + Ifges(2,5) * qJDD(1) - t281 * Ifges(2,6) - qJ(2) * t131 - t273 * t126 + t274 * t128;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t279 * t123 - t277 * t124 - pkin(5) * (t277 * t129 + t279 * t149), t123, t128, t133, -t241 * t205 + t288, t152; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t277 * t123 + t279 * t124 + pkin(5) * (t279 * t129 - t277 * t149), t124, t126, t132, -t285, t151; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t290, t290, qJDD(1) * t296 - t322, -t283, Ifges(5,5) * t228 + Ifges(5,6) * qJDD(3) + Ifges(5,3) * t227 - qJD(3) * t207 + t242 * t205 - t287, t291;];
m_new = t1;
