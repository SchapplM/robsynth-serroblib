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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:30:42
% EndTime: 2020-01-03 11:31:01
% DurationCPUTime: 8.24s
% Computational Cost: add. (104853->279), mult. (287853->380), div. (0->0), fcn. (198151->10), ass. (0->133)
t271 = sin(qJ(1));
t274 = cos(qJ(1));
t247 = -t271 * g(2) + t274 * g(3);
t275 = qJD(1) ^ 2;
t321 = -t275 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t247;
t248 = -t274 * g(2) - t271 * g(3);
t288 = -t275 * qJ(2) + qJDD(2) - t248;
t266 = sin(pkin(8));
t268 = cos(pkin(8));
t296 = -pkin(2) * t268 - qJ(3) * t266;
t312 = qJD(1) * t266;
t320 = (-pkin(1) + t296) * qJDD(1) + t288 - 0.2e1 * qJD(3) * t312;
t217 = -t268 * g(1) - t321 * t266;
t311 = t268 * qJD(1);
t218 = -t266 * g(1) + t321 * t268;
t236 = t296 * qJD(1);
t205 = t236 * t311 + t218;
t263 = t266 ^ 2;
t265 = sin(pkin(9));
t267 = cos(pkin(9));
t315 = t266 * t267;
t291 = -pkin(3) * t268 - pkin(6) * t315;
t314 = t320 * t267;
t183 = t291 * qJDD(1) + (-t205 + (-pkin(3) * t263 * t267 + pkin(6) * t266 * t268) * t275) * t265 + t314;
t192 = t267 * t205 + t320 * t265;
t235 = t291 * qJD(1);
t309 = qJDD(1) * t266;
t305 = t265 * t309;
t317 = t263 * t275;
t307 = t265 ^ 2 * t317;
t184 = -pkin(3) * t307 - pkin(6) * t305 + t235 * t311 + t192;
t270 = sin(qJ(4));
t273 = cos(qJ(4));
t169 = t273 * t183 - t270 * t184;
t287 = (-t265 * t273 - t267 * t270) * t266;
t225 = qJD(1) * t287;
t286 = (-t265 * t270 + t267 * t273) * t266;
t211 = t225 * qJD(4) + qJDD(1) * t286;
t226 = qJD(1) * t286;
t308 = t268 * qJDD(1);
t251 = qJDD(4) - t308;
t252 = qJD(4) - t311;
t166 = (t225 * t252 - t211) * pkin(7) + (t225 * t226 + t251) * pkin(4) + t169;
t170 = t270 * t183 + t273 * t184;
t210 = -t226 * qJD(4) + qJDD(1) * t287;
t216 = t252 * pkin(4) - t226 * pkin(7);
t224 = t225 ^ 2;
t167 = -t224 * pkin(4) + t210 * pkin(7) - t252 * t216 + t170;
t269 = sin(qJ(5));
t272 = cos(qJ(5));
t165 = t269 * t166 + t272 * t167;
t204 = t236 * t312 + qJDD(3) - t217;
t193 = t267 * t235 * t312 + pkin(3) * t305 - pkin(6) * t307 + t204;
t172 = -t210 * pkin(4) - t224 * pkin(7) + t226 * t216 + t193;
t203 = t269 * t225 + t272 * t226;
t178 = -t203 * qJD(5) + t272 * t210 - t269 * t211;
t202 = t272 * t225 - t269 * t226;
t179 = t202 * qJD(5) + t269 * t210 + t272 * t211;
t250 = qJD(5) + t252;
t185 = Ifges(6,5) * t203 + Ifges(6,6) * t202 + Ifges(6,3) * t250;
t187 = Ifges(6,1) * t203 + Ifges(6,4) * t202 + Ifges(6,5) * t250;
t246 = qJDD(5) + t251;
t154 = -mrSges(6,1) * t172 + mrSges(6,3) * t165 + Ifges(6,4) * t179 + Ifges(6,2) * t178 + Ifges(6,6) * t246 - t203 * t185 + t250 * t187;
t164 = t272 * t166 - t269 * t167;
t186 = Ifges(6,4) * t203 + Ifges(6,2) * t202 + Ifges(6,6) * t250;
t155 = mrSges(6,2) * t172 - mrSges(6,3) * t164 + Ifges(6,1) * t179 + Ifges(6,4) * t178 + Ifges(6,5) * t246 + t202 * t185 - t250 * t186;
t197 = Ifges(5,5) * t226 + Ifges(5,6) * t225 + Ifges(5,3) * t252;
t199 = Ifges(5,1) * t226 + Ifges(5,4) * t225 + Ifges(5,5) * t252;
t195 = -t250 * mrSges(6,2) + t202 * mrSges(6,3);
t196 = t250 * mrSges(6,1) - t203 * mrSges(6,3);
t295 = m(6) * t172 - t178 * mrSges(6,1) + t179 * mrSges(6,2) - t202 * t195 + t203 * t196;
t190 = -t202 * mrSges(6,1) + t203 * mrSges(6,2);
t160 = m(6) * t164 + t246 * mrSges(6,1) - t179 * mrSges(6,3) - t203 * t190 + t250 * t195;
t161 = m(6) * t165 - t246 * mrSges(6,2) + t178 * mrSges(6,3) + t202 * t190 - t250 * t196;
t302 = -t269 * t160 + t272 * t161;
t141 = -mrSges(5,1) * t193 + mrSges(5,3) * t170 + Ifges(5,4) * t211 + Ifges(5,2) * t210 + Ifges(5,6) * t251 - pkin(4) * t295 + pkin(7) * t302 + t272 * t154 + t269 * t155 - t226 * t197 + t252 * t199;
t153 = t272 * t160 + t269 * t161;
t198 = Ifges(5,4) * t226 + Ifges(5,2) * t225 + Ifges(5,6) * t252;
t142 = mrSges(5,2) * t193 - mrSges(5,3) * t169 + Ifges(5,1) * t211 + Ifges(5,4) * t210 + Ifges(5,5) * t251 - pkin(7) * t153 - t269 * t154 + t272 * t155 + t225 * t197 - t252 * t198;
t297 = Ifges(4,5) * t267 - Ifges(4,6) * t265;
t220 = (-Ifges(4,3) * t268 + t297 * t266) * qJD(1);
t284 = -Ifges(4,5) * t268 + (Ifges(4,1) * t267 - Ifges(4,4) * t265) * t266;
t222 = t284 * qJD(1);
t212 = -t252 * mrSges(5,2) + t225 * mrSges(5,3);
t213 = t252 * mrSges(5,1) - t226 * mrSges(5,3);
t280 = m(5) * t193 - t210 * mrSges(5,1) + t211 * mrSges(5,2) - t225 * t212 + t226 * t213 + t295;
t283 = -Ifges(4,6) * t268 + (Ifges(4,4) * t267 - Ifges(4,2) * t265) * t266;
t206 = -t225 * mrSges(5,1) + t226 * mrSges(5,2);
t150 = m(5) * t169 + t251 * mrSges(5,1) - t211 * mrSges(5,3) - t226 * t206 + t252 * t212 + t153;
t151 = m(5) * t170 - t251 * mrSges(5,2) + t210 * mrSges(5,3) + t225 * t206 - t252 * t213 + t302;
t303 = -t270 * t150 + t273 * t151;
t130 = -mrSges(4,1) * t204 + mrSges(4,3) * t192 + t270 * t142 + t273 * t141 - pkin(3) * t280 + pkin(6) * t303 + (-t220 * t315 - t268 * t222) * qJD(1) + t283 * qJDD(1);
t146 = t273 * t150 + t270 * t151;
t191 = -t265 * t205 + t314;
t221 = t283 * qJD(1);
t316 = t265 * t266;
t131 = mrSges(4,2) * t204 - mrSges(4,3) * t191 - pkin(6) * t146 - t270 * t141 + t273 * t142 + (-t220 * t316 + t221 * t268) * qJD(1) + t284 * qJDD(1);
t300 = mrSges(4,1) * t265 + mrSges(4,2) * t267;
t229 = t300 * t312;
t289 = mrSges(4,2) * t268 - mrSges(4,3) * t316;
t232 = t289 * qJD(1);
t290 = -mrSges(4,1) * t268 - mrSges(4,3) * t315;
t144 = m(4) * t191 + t290 * qJDD(1) + (-t229 * t315 - t232 * t268) * qJD(1) + t146;
t233 = t290 * qJD(1);
t145 = m(4) * t192 + t289 * qJDD(1) + (-t229 * t316 + t233 * t268) * qJD(1) + t303;
t140 = -t265 * t144 + t267 * t145;
t278 = -m(4) * t204 - t280;
t293 = -t232 * t265 - t233 * t267;
t299 = Ifges(3,1) * t266 + Ifges(3,4) * t268;
t319 = -((Ifges(3,4) * t266 + Ifges(3,2) * t268) * t312 - t299 * t311) * qJD(1) - mrSges(3,1) * t217 + mrSges(3,2) * t218 - pkin(2) * ((t293 * qJD(1) - t300 * qJDD(1)) * t266 + t278) - qJ(3) * t140 - t267 * t130 - t265 * t131;
t318 = mrSges(3,2) * t266;
t237 = (-mrSges(3,1) * t268 + t318) * qJD(1);
t137 = m(3) * t218 + (qJDD(1) * mrSges(3,3) + qJD(1) * t237) * t268 + t140;
t156 = m(3) * t217 + ((-mrSges(3,3) - t300) * qJDD(1) + (-t237 + t293) * qJD(1)) * t266 + t278;
t304 = t268 * t137 - t266 * t156;
t298 = Ifges(3,5) * t266 + Ifges(3,6) * t268;
t139 = t267 * t144 + t265 * t145;
t294 = t221 * t267 + t222 * t265;
t231 = -qJDD(1) * pkin(1) + t288;
t238 = t298 * qJD(1);
t127 = mrSges(3,2) * t231 - mrSges(3,3) * t217 - qJ(3) * t139 + t299 * qJDD(1) - t265 * t130 + t267 * t131 + t238 * t311;
t282 = -mrSges(6,1) * t164 + mrSges(6,2) * t165 - Ifges(6,5) * t179 - Ifges(6,6) * t178 - Ifges(6,3) * t246 - t203 * t186 + t202 * t187;
t277 = -mrSges(5,1) * t169 + mrSges(5,2) * t170 - Ifges(5,5) * t211 - Ifges(5,6) * t210 - Ifges(5,3) * t251 - pkin(4) * t153 - t226 * t198 + t225 * t199 + t282;
t276 = mrSges(4,1) * t191 - mrSges(4,2) * t192 + pkin(3) * t146 - t277;
t129 = ((Ifges(3,4) - t297) * qJDD(1) + (-t238 - t294) * qJD(1)) * t266 - mrSges(3,1) * t231 + mrSges(3,3) * t218 - pkin(2) * t139 - t276 + (Ifges(3,2) + Ifges(4,3)) * t308;
t281 = -m(3) * t231 + mrSges(3,1) * t308 - t139 + (t268 ^ 2 * t275 + t317) * mrSges(3,3);
t285 = -mrSges(2,2) * t247 + qJ(2) * t304 + t266 * t127 + t268 * t129 + pkin(1) * (-mrSges(3,2) * t309 + t281) + mrSges(2,1) * t248 + Ifges(2,3) * qJDD(1);
t135 = m(2) * t248 - t275 * mrSges(2,2) + (mrSges(2,1) - t318) * qJDD(1) + t281;
t134 = t266 * t137 + t268 * t156;
t132 = m(2) * t247 - t275 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t304;
t125 = mrSges(2,1) * g(1) + mrSges(2,3) * t247 + t275 * Ifges(2,5) - pkin(1) * t134 + (Ifges(2,6) - t298) * qJDD(1) + t319;
t124 = -mrSges(2,2) * g(1) - mrSges(2,3) * t248 + Ifges(2,5) * qJDD(1) - t275 * Ifges(2,6) - qJ(2) * t134 + t268 * t127 - t266 * t129;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t285, t124, t127, t131, t142, t155; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t271 * t124 + t274 * t125 - pkin(5) * (-t274 * t132 + t271 * t135), t125, t129, t130, t141, t154; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t274 * t124 + t271 * t125 + pkin(5) * (t271 * t132 + t274 * t135), t285, t298 * qJDD(1) - t319, (t294 * qJD(1) + t297 * qJDD(1)) * t266 + t276 - Ifges(4,3) * t308, -t277, -t282;];
m_new = t1;
