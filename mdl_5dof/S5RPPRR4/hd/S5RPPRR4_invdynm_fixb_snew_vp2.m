% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:19
% EndTime: 2019-12-05 17:44:35
% DurationCPUTime: 8.15s
% Computational Cost: add. (104853->279), mult. (287853->380), div. (0->0), fcn. (198151->10), ass. (0->133)
t267 = sin(qJ(1));
t270 = cos(qJ(1));
t245 = t267 * g(2) - t270 * g(3);
t271 = qJD(1) ^ 2;
t317 = -t271 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t245;
t246 = t270 * g(2) + t267 * g(3);
t282 = -t271 * qJ(2) + qJDD(2) - t246;
t262 = sin(pkin(8));
t264 = cos(pkin(8));
t292 = -pkin(2) * t264 - qJ(3) * t262;
t308 = qJD(1) * t262;
t316 = (-pkin(1) + t292) * qJDD(1) + t282 - 0.2e1 * qJD(3) * t308;
t215 = -t264 * g(1) - t317 * t262;
t307 = t264 * qJD(1);
t216 = -t262 * g(1) + t317 * t264;
t234 = t292 * qJD(1);
t203 = t234 * t307 + t216;
t259 = t262 ^ 2;
t261 = sin(pkin(9));
t263 = cos(pkin(9));
t311 = t262 * t263;
t287 = -pkin(3) * t264 - pkin(6) * t311;
t310 = t316 * t263;
t181 = t287 * qJDD(1) + (-t203 + (-pkin(3) * t259 * t263 + pkin(6) * t262 * t264) * t271) * t261 + t310;
t190 = t263 * t203 + t316 * t261;
t233 = t287 * qJD(1);
t305 = qJDD(1) * t262;
t301 = t261 * t305;
t313 = t259 * t271;
t303 = t261 ^ 2 * t313;
t182 = -pkin(3) * t303 - pkin(6) * t301 + t233 * t307 + t190;
t266 = sin(qJ(4));
t269 = cos(qJ(4));
t167 = t269 * t181 - t266 * t182;
t284 = (-t261 * t269 - t263 * t266) * t262;
t223 = qJD(1) * t284;
t283 = (-t261 * t266 + t263 * t269) * t262;
t209 = t223 * qJD(4) + qJDD(1) * t283;
t224 = qJD(1) * t283;
t304 = t264 * qJDD(1);
t249 = qJDD(4) - t304;
t250 = qJD(4) - t307;
t164 = (t223 * t250 - t209) * pkin(7) + (t223 * t224 + t249) * pkin(4) + t167;
t168 = t266 * t181 + t269 * t182;
t208 = -t224 * qJD(4) + qJDD(1) * t284;
t214 = t250 * pkin(4) - t224 * pkin(7);
t222 = t223 ^ 2;
t165 = -t222 * pkin(4) + t208 * pkin(7) - t250 * t214 + t168;
t265 = sin(qJ(5));
t268 = cos(qJ(5));
t163 = t265 * t164 + t268 * t165;
t202 = t234 * t308 + qJDD(3) - t215;
t191 = t263 * t233 * t308 + pkin(3) * t301 - pkin(6) * t303 + t202;
t170 = -t208 * pkin(4) - t222 * pkin(7) + t224 * t214 + t191;
t201 = t265 * t223 + t268 * t224;
t176 = -t201 * qJD(5) + t268 * t208 - t265 * t209;
t200 = t268 * t223 - t265 * t224;
t177 = t200 * qJD(5) + t265 * t208 + t268 * t209;
t248 = qJD(5) + t250;
t183 = Ifges(6,5) * t201 + Ifges(6,6) * t200 + Ifges(6,3) * t248;
t185 = Ifges(6,1) * t201 + Ifges(6,4) * t200 + Ifges(6,5) * t248;
t244 = qJDD(5) + t249;
t152 = -mrSges(6,1) * t170 + mrSges(6,3) * t163 + Ifges(6,4) * t177 + Ifges(6,2) * t176 + Ifges(6,6) * t244 - t201 * t183 + t248 * t185;
t162 = t268 * t164 - t265 * t165;
t184 = Ifges(6,4) * t201 + Ifges(6,2) * t200 + Ifges(6,6) * t248;
t153 = mrSges(6,2) * t170 - mrSges(6,3) * t162 + Ifges(6,1) * t177 + Ifges(6,4) * t176 + Ifges(6,5) * t244 + t200 * t183 - t248 * t184;
t195 = Ifges(5,5) * t224 + Ifges(5,6) * t223 + Ifges(5,3) * t250;
t197 = Ifges(5,1) * t224 + Ifges(5,4) * t223 + Ifges(5,5) * t250;
t193 = -t248 * mrSges(6,2) + t200 * mrSges(6,3);
t194 = t248 * mrSges(6,1) - t201 * mrSges(6,3);
t291 = m(6) * t170 - t176 * mrSges(6,1) + t177 * mrSges(6,2) - t200 * t193 + t201 * t194;
t188 = -t200 * mrSges(6,1) + t201 * mrSges(6,2);
t158 = m(6) * t162 + t244 * mrSges(6,1) - t177 * mrSges(6,3) - t201 * t188 + t248 * t193;
t159 = m(6) * t163 - t244 * mrSges(6,2) + t176 * mrSges(6,3) + t200 * t188 - t248 * t194;
t298 = -t265 * t158 + t268 * t159;
t139 = -mrSges(5,1) * t191 + mrSges(5,3) * t168 + Ifges(5,4) * t209 + Ifges(5,2) * t208 + Ifges(5,6) * t249 - pkin(4) * t291 + pkin(7) * t298 + t268 * t152 + t265 * t153 - t224 * t195 + t250 * t197;
t151 = t268 * t158 + t265 * t159;
t196 = Ifges(5,4) * t224 + Ifges(5,2) * t223 + Ifges(5,6) * t250;
t140 = mrSges(5,2) * t191 - mrSges(5,3) * t167 + Ifges(5,1) * t209 + Ifges(5,4) * t208 + Ifges(5,5) * t249 - pkin(7) * t151 - t265 * t152 + t268 * t153 + t223 * t195 - t250 * t196;
t293 = Ifges(4,5) * t263 - Ifges(4,6) * t261;
t218 = (-Ifges(4,3) * t264 + t293 * t262) * qJD(1);
t280 = -Ifges(4,5) * t264 + (Ifges(4,1) * t263 - Ifges(4,4) * t261) * t262;
t220 = t280 * qJD(1);
t210 = -t250 * mrSges(5,2) + t223 * mrSges(5,3);
t211 = t250 * mrSges(5,1) - t224 * mrSges(5,3);
t276 = m(5) * t191 - t208 * mrSges(5,1) + t209 * mrSges(5,2) - t223 * t210 + t224 * t211 + t291;
t279 = -Ifges(4,6) * t264 + (Ifges(4,4) * t263 - Ifges(4,2) * t261) * t262;
t204 = -t223 * mrSges(5,1) + t224 * mrSges(5,2);
t148 = m(5) * t167 + t249 * mrSges(5,1) - t209 * mrSges(5,3) - t224 * t204 + t250 * t210 + t151;
t149 = m(5) * t168 - t249 * mrSges(5,2) + t208 * mrSges(5,3) + t223 * t204 - t250 * t211 + t298;
t299 = -t266 * t148 + t269 * t149;
t128 = -mrSges(4,1) * t202 + mrSges(4,3) * t190 + t266 * t140 + t269 * t139 - pkin(3) * t276 + pkin(6) * t299 + (-t218 * t311 - t264 * t220) * qJD(1) + t279 * qJDD(1);
t144 = t269 * t148 + t266 * t149;
t189 = -t261 * t203 + t310;
t219 = t279 * qJD(1);
t312 = t261 * t262;
t129 = mrSges(4,2) * t202 - mrSges(4,3) * t189 - pkin(6) * t144 - t266 * t139 + t269 * t140 + (-t218 * t312 + t219 * t264) * qJD(1) + t280 * qJDD(1);
t296 = mrSges(4,1) * t261 + mrSges(4,2) * t263;
t227 = t296 * t308;
t285 = mrSges(4,2) * t264 - mrSges(4,3) * t312;
t230 = t285 * qJD(1);
t286 = -mrSges(4,1) * t264 - mrSges(4,3) * t311;
t142 = m(4) * t189 + t286 * qJDD(1) + (-t227 * t311 - t230 * t264) * qJD(1) + t144;
t231 = t286 * qJD(1);
t143 = m(4) * t190 + t285 * qJDD(1) + (-t227 * t312 + t231 * t264) * qJD(1) + t299;
t138 = -t261 * t142 + t263 * t143;
t274 = -m(4) * t202 - t276;
t289 = -t230 * t261 - t231 * t263;
t295 = Ifges(3,1) * t262 + Ifges(3,4) * t264;
t315 = -((Ifges(3,4) * t262 + Ifges(3,2) * t264) * t308 - t295 * t307) * qJD(1) - mrSges(3,1) * t215 + mrSges(3,2) * t216 - pkin(2) * ((t289 * qJD(1) - t296 * qJDD(1)) * t262 + t274) - qJ(3) * t138 - t263 * t128 - t261 * t129;
t314 = mrSges(3,2) * t262;
t235 = (-mrSges(3,1) * t264 + t314) * qJD(1);
t135 = m(3) * t216 + (qJDD(1) * mrSges(3,3) + qJD(1) * t235) * t264 + t138;
t154 = t274 + ((-mrSges(3,3) - t296) * qJDD(1) + (-t235 + t289) * qJD(1)) * t262 + m(3) * t215;
t300 = t264 * t135 - t262 * t154;
t294 = Ifges(3,5) * t262 + Ifges(3,6) * t264;
t137 = t263 * t142 + t261 * t143;
t290 = t219 * t263 + t220 * t261;
t229 = -qJDD(1) * pkin(1) + t282;
t236 = t294 * qJD(1);
t125 = mrSges(3,2) * t229 - mrSges(3,3) * t215 - qJ(3) * t137 + t295 * qJDD(1) - t261 * t128 + t263 * t129 + t236 * t307;
t278 = -mrSges(6,1) * t162 + mrSges(6,2) * t163 - Ifges(6,5) * t177 - Ifges(6,6) * t176 - Ifges(6,3) * t244 - t201 * t184 + t200 * t185;
t273 = -mrSges(5,1) * t167 + mrSges(5,2) * t168 - Ifges(5,5) * t209 - Ifges(5,6) * t208 - Ifges(5,3) * t249 - pkin(4) * t151 - t224 * t196 + t223 * t197 + t278;
t272 = mrSges(4,1) * t189 - mrSges(4,2) * t190 + pkin(3) * t144 - t273;
t127 = ((Ifges(3,4) - t293) * qJDD(1) + (-t236 - t290) * qJD(1)) * t262 - t272 - mrSges(3,1) * t229 + mrSges(3,3) * t216 - pkin(2) * t137 + (Ifges(3,2) + Ifges(4,3)) * t304;
t277 = -m(3) * t229 + mrSges(3,1) * t304 - t137 + (t264 ^ 2 * t271 + t313) * mrSges(3,3);
t281 = -mrSges(2,2) * t245 + qJ(2) * t300 + t262 * t125 + t264 * t127 + pkin(1) * (-mrSges(3,2) * t305 + t277) + mrSges(2,1) * t246 + Ifges(2,3) * qJDD(1);
t133 = m(2) * t246 - t271 * mrSges(2,2) + (mrSges(2,1) - t314) * qJDD(1) + t277;
t132 = t262 * t135 + t264 * t154;
t130 = m(2) * t245 - t271 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t300;
t123 = mrSges(2,1) * g(1) + mrSges(2,3) * t245 + t271 * Ifges(2,5) - pkin(1) * t132 + (Ifges(2,6) - t294) * qJDD(1) + t315;
t122 = -mrSges(2,2) * g(1) - mrSges(2,3) * t246 + Ifges(2,5) * qJDD(1) - t271 * Ifges(2,6) - qJ(2) * t132 + t264 * t125 - t262 * t127;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t281, t122, t125, t129, t140, t153; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t267 * t122 - t270 * t123 - pkin(5) * (t270 * t130 - t267 * t133), t123, t127, t128, t139, t152; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t270 * t122 - t267 * t123 + pkin(5) * (-t267 * t130 - t270 * t133), t281, t294 * qJDD(1) - t315, t272 + (t290 * qJD(1) + t293 * qJDD(1)) * t262 - Ifges(4,3) * t304, -t273, -t278;];
m_new = t1;
