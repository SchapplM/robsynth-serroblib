% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:41:44
% EndTime: 2020-01-03 11:42:10
% DurationCPUTime: 8.84s
% Computational Cost: add. (113962->288), mult. (293736->378), div. (0->0), fcn. (200705->10), ass. (0->125)
t268 = sin(qJ(1));
t271 = cos(qJ(1));
t245 = -t268 * g(2) + t271 * g(3);
t272 = qJD(1) ^ 2;
t306 = -t272 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t245;
t263 = sin(pkin(8));
t265 = cos(pkin(8));
t212 = -t265 * g(1) - t306 * t263;
t301 = qJD(1) * t263;
t299 = t265 * qJD(1);
t213 = -t263 * g(1) + t306 * t265;
t288 = -pkin(2) * t265 - pkin(6) * t263;
t240 = t288 * qJD(1);
t200 = t240 * t299 + t213;
t246 = -t271 * g(2) - t268 * g(3);
t281 = -t272 * qJ(2) + qJDD(2) - t246;
t214 = (-pkin(1) + t288) * qJDD(1) + t281;
t270 = cos(qJ(3));
t211 = t270 * t214;
t267 = sin(qJ(3));
t297 = qJD(1) * qJD(3);
t232 = (qJDD(1) * t270 - t267 * t297) * t263;
t296 = t265 * qJDD(1);
t248 = qJDD(3) - t296;
t249 = qJD(3) - t299;
t303 = t263 ^ 2 * t272;
t178 = t248 * pkin(3) - t232 * qJ(4) + t211 + (-pkin(3) * t270 * t303 - qJ(4) * t249 * t301 - t200) * t267;
t188 = t270 * t200 + t267 * t214;
t293 = t270 * t301;
t228 = t249 * pkin(3) - qJ(4) * t293;
t231 = (-qJDD(1) * t267 - t270 * t297) * t263;
t295 = t267 ^ 2 * t303;
t179 = -pkin(3) * t295 + t231 * qJ(4) - t249 * t228 + t188;
t262 = sin(pkin(9));
t264 = cos(pkin(9));
t223 = (-t262 * t267 + t264 * t270) * t301;
t164 = -0.2e1 * qJD(4) * t223 + t264 * t178 - t262 * t179;
t206 = t262 * t231 + t264 * t232;
t222 = (-t262 * t270 - t264 * t267) * t301;
t161 = (t222 * t249 - t206) * pkin(7) + (t222 * t223 + t248) * pkin(4) + t164;
t165 = 0.2e1 * qJD(4) * t222 + t262 * t178 + t264 * t179;
t205 = t264 * t231 - t262 * t232;
t209 = t249 * pkin(4) - t223 * pkin(7);
t221 = t222 ^ 2;
t162 = -t221 * pkin(4) + t205 * pkin(7) - t249 * t209 + t165;
t266 = sin(qJ(5));
t269 = cos(qJ(5));
t160 = t266 * t161 + t269 * t162;
t199 = t240 * t301 - t212;
t186 = -t231 * pkin(3) - qJ(4) * t295 + t228 * t293 + qJDD(4) + t199;
t167 = -t205 * pkin(4) - t221 * pkin(7) + t223 * t209 + t186;
t198 = t266 * t222 + t269 * t223;
t173 = -t198 * qJD(5) + t269 * t205 - t266 * t206;
t197 = t269 * t222 - t266 * t223;
t174 = t197 * qJD(5) + t266 * t205 + t269 * t206;
t247 = qJD(5) + t249;
t180 = Ifges(6,5) * t198 + Ifges(6,6) * t197 + Ifges(6,3) * t247;
t182 = Ifges(6,1) * t198 + Ifges(6,4) * t197 + Ifges(6,5) * t247;
t244 = qJDD(5) + t248;
t149 = -mrSges(6,1) * t167 + mrSges(6,3) * t160 + Ifges(6,4) * t174 + Ifges(6,2) * t173 + Ifges(6,6) * t244 - t198 * t180 + t247 * t182;
t159 = t269 * t161 - t266 * t162;
t181 = Ifges(6,4) * t198 + Ifges(6,2) * t197 + Ifges(6,6) * t247;
t150 = mrSges(6,2) * t167 - mrSges(6,3) * t159 + Ifges(6,1) * t174 + Ifges(6,4) * t173 + Ifges(6,5) * t244 + t197 * t180 - t247 * t181;
t192 = Ifges(5,5) * t223 + Ifges(5,6) * t222 + Ifges(5,3) * t249;
t194 = Ifges(5,1) * t223 + Ifges(5,4) * t222 + Ifges(5,5) * t249;
t190 = -t247 * mrSges(6,2) + t197 * mrSges(6,3);
t191 = t247 * mrSges(6,1) - t198 * mrSges(6,3);
t285 = m(6) * t167 - t173 * mrSges(6,1) + t174 * mrSges(6,2) - t197 * t190 + t198 * t191;
t185 = -t197 * mrSges(6,1) + t198 * mrSges(6,2);
t155 = m(6) * t159 + mrSges(6,1) * t244 - mrSges(6,3) * t174 - t185 * t198 + t190 * t247;
t156 = m(6) * t160 - mrSges(6,2) * t244 + mrSges(6,3) * t173 + t185 * t197 - t191 * t247;
t289 = -t155 * t266 + t156 * t269;
t136 = -mrSges(5,1) * t186 + mrSges(5,3) * t165 + Ifges(5,4) * t206 + Ifges(5,2) * t205 + Ifges(5,6) * t248 - pkin(4) * t285 + pkin(7) * t289 + t269 * t149 + t266 * t150 - t223 * t192 + t249 * t194;
t148 = t155 * t269 + t156 * t266;
t193 = Ifges(5,4) * t223 + Ifges(5,2) * t222 + Ifges(5,6) * t249;
t137 = mrSges(5,2) * t186 - mrSges(5,3) * t164 + Ifges(5,1) * t206 + Ifges(5,4) * t205 + Ifges(5,5) * t248 - pkin(7) * t148 - t149 * t266 + t150 * t269 + t192 * t222 - t193 * t249;
t215 = Ifges(4,3) * t249 + (Ifges(4,5) * t270 - Ifges(4,6) * t267) * t301;
t217 = Ifges(4,5) * t249 + (Ifges(4,1) * t270 - Ifges(4,4) * t267) * t301;
t207 = -mrSges(5,2) * t249 + mrSges(5,3) * t222;
t208 = mrSges(5,1) * t249 - mrSges(5,3) * t223;
t277 = m(5) * t186 - mrSges(5,1) * t205 + t206 * mrSges(5,2) - t207 * t222 + t223 * t208 + t285;
t201 = -mrSges(5,1) * t222 + mrSges(5,2) * t223;
t145 = m(5) * t164 + mrSges(5,1) * t248 - mrSges(5,3) * t206 - t201 * t223 + t207 * t249 + t148;
t146 = m(5) * t165 - mrSges(5,2) * t248 + mrSges(5,3) * t205 + t201 * t222 - t208 * t249 + t289;
t290 = -t145 * t262 + t146 * t264;
t125 = -mrSges(4,1) * t199 + mrSges(4,3) * t188 + Ifges(4,4) * t232 + Ifges(4,2) * t231 + Ifges(4,6) * t248 - pkin(3) * t277 + qJ(4) * t290 + t264 * t136 + t262 * t137 - t215 * t293 + t249 * t217;
t141 = t145 * t264 + t146 * t262;
t187 = -t267 * t200 + t211;
t216 = Ifges(4,6) * t249 + (Ifges(4,4) * t270 - Ifges(4,2) * t267) * t301;
t294 = t267 * t301;
t126 = mrSges(4,2) * t199 - mrSges(4,3) * t187 + Ifges(4,1) * t232 + Ifges(4,4) * t231 + Ifges(4,5) * t248 - qJ(4) * t141 - t136 * t262 + t137 * t264 - t215 * t294 - t216 * t249;
t227 = -mrSges(4,2) * t249 - mrSges(4,3) * t294;
t230 = (mrSges(4,1) * t267 + mrSges(4,2) * t270) * t301;
t139 = m(4) * t187 + mrSges(4,1) * t248 - mrSges(4,3) * t232 + t227 * t249 - t230 * t293 + t141;
t229 = mrSges(4,1) * t249 - mrSges(4,3) * t293;
t140 = m(4) * t188 - mrSges(4,2) * t248 + mrSges(4,3) * t231 - t229 * t249 - t230 * t294 + t290;
t135 = -t267 * t139 + t140 * t270;
t275 = -m(4) * t199 + t231 * mrSges(4,1) - t232 * mrSges(4,2) - t277;
t283 = -t227 * t267 - t229 * t270;
t287 = Ifges(3,1) * t263 + Ifges(3,4) * t265;
t305 = -((Ifges(3,4) * t263 + Ifges(3,2) * t265) * t301 - t287 * t299) * qJD(1) - mrSges(3,1) * t212 + mrSges(3,2) * t213 - pkin(2) * (t283 * t301 + t275) - pkin(6) * t135 - t270 * t125 - t267 * t126;
t304 = mrSges(3,2) * t263;
t300 = qJDD(1) * mrSges(3,3);
t236 = (-mrSges(3,1) * t265 + t304) * qJD(1);
t132 = m(3) * t213 + (qJD(1) * t236 + t300) * t265 + t135;
t154 = t275 + (-t300 + (-t236 + t283) * qJD(1)) * t263 + m(3) * t212;
t291 = t132 * t265 - t154 * t263;
t286 = Ifges(3,5) * t263 + Ifges(3,6) * t265;
t134 = t270 * t139 + t267 * t140;
t284 = t216 * t270 + t217 * t267;
t234 = -qJDD(1) * pkin(1) + t281;
t237 = t286 * qJD(1);
t122 = mrSges(3,2) * t234 - mrSges(3,3) * t212 - pkin(6) * t134 + qJDD(1) * t287 - t267 * t125 + t270 * t126 + t237 * t299;
t279 = -mrSges(6,1) * t159 + mrSges(6,2) * t160 - Ifges(6,5) * t174 - Ifges(6,6) * t173 - Ifges(6,3) * t244 - t181 * t198 + t197 * t182;
t274 = -mrSges(5,1) * t164 + mrSges(5,2) * t165 - Ifges(5,5) * t206 - Ifges(5,6) * t205 - Ifges(5,3) * t248 - pkin(4) * t148 - t193 * t223 + t222 * t194 + t279;
t273 = mrSges(4,1) * t187 - mrSges(4,2) * t188 + Ifges(4,5) * t232 + Ifges(4,6) * t231 + Ifges(4,3) * t248 + pkin(3) * t141 - t274;
t124 = Ifges(3,2) * t296 - t273 - pkin(2) * t134 + (Ifges(3,4) * qJDD(1) + (-t237 - t284) * qJD(1)) * t263 - mrSges(3,1) * t234 + mrSges(3,3) * t213;
t278 = -m(3) * t234 + mrSges(3,1) * t296 - t134 + (t265 ^ 2 * t272 + t303) * mrSges(3,3);
t280 = -mrSges(2,2) * t245 + qJ(2) * t291 + t263 * t122 + t265 * t124 + pkin(1) * (-qJDD(1) * t304 + t278) + mrSges(2,1) * t246 + Ifges(2,3) * qJDD(1);
t130 = m(2) * t246 - t272 * mrSges(2,2) + (mrSges(2,1) - t304) * qJDD(1) + t278;
t129 = t132 * t263 + t154 * t265;
t127 = m(2) * t245 - mrSges(2,1) * t272 - qJDD(1) * mrSges(2,2) + t291;
t120 = mrSges(2,1) * g(1) + mrSges(2,3) * t245 + t272 * Ifges(2,5) - pkin(1) * t129 + (Ifges(2,6) - t286) * qJDD(1) + t305;
t119 = -mrSges(2,2) * g(1) - mrSges(2,3) * t246 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t272 - qJ(2) * t129 + t122 * t265 - t124 * t263;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t280, t119, t122, t126, t137, t150; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t268 * t119 + t271 * t120 - pkin(5) * (-t127 * t271 + t130 * t268), t120, t124, t125, t136, t149; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t271 * t119 + t268 * t120 + pkin(5) * (t127 * t268 + t130 * t271), t280, qJDD(1) * t286 - t305, t284 * t301 + t273, -t274, -t279;];
m_new = t1;
