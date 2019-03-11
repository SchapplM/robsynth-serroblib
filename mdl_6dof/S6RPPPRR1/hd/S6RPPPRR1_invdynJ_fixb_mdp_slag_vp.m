% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPPPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:19
% EndTime: 2019-03-09 01:30:21
% DurationCPUTime: 1.82s
% Computational Cost: add. (972->280), mult. (1611->353), div. (0->0), fcn. (934->10), ass. (0->123)
t225 = sin(pkin(9));
t201 = pkin(1) * t225 + qJ(3);
t231 = cos(qJ(5));
t273 = qJDD(1) * t231;
t315 = qJD(5) * qJD(6) + t273;
t308 = pkin(1) * qJDD(1);
t216 = qJ(1) + pkin(9);
t209 = sin(t216);
t210 = cos(t216);
t314 = -g(1) * t209 + g(2) * t210;
t219 = qJD(3) * qJD(1);
t293 = qJDD(1) * qJ(3) + t225 * t308;
t267 = -t219 - t293;
t256 = qJDD(4) - t267;
t177 = -qJDD(1) * pkin(7) + t256;
t228 = sin(qJ(5));
t193 = t201 * qJD(1);
t190 = qJD(4) + t193;
t181 = -qJD(1) * pkin(7) + t190;
t172 = qJD(2) * t231 + t181 * t228;
t287 = qJD(5) * t172;
t162 = -qJDD(5) * pkin(5) + qJDD(2) * t228 - t177 * t231 + t287;
t171 = -qJD(2) * t228 + t181 * t231;
t168 = -qJD(5) * pkin(5) - t171;
t226 = cos(pkin(9));
t204 = -pkin(1) * t226 - pkin(2);
t197 = qJ(4) - t204;
t254 = pkin(5) * t228 - pkin(8) * t231;
t180 = t197 + t254;
t278 = qJD(1) * qJD(5);
t264 = t231 * t278;
t274 = qJDD(1) * t228;
t184 = qJDD(6) + t264 + t274;
t196 = -pkin(7) + t201;
t198 = qJD(1) * t228 + qJD(6);
t170 = qJD(1) * t180 - qJD(3);
t288 = qJD(5) * t171;
t259 = -qJDD(5) * pkin(8) - qJD(6) * t170 - qJDD(2) * t231 - t177 * t228 - t288;
t284 = qJD(5) * t231;
t313 = -t198 * (qJD(6) * t180 + t196 * t284) + t162 * t231 + (-qJD(3) * t198 - t168 * qJD(5) - t196 * t184 + t259) * t228;
t290 = qJD(1) * t197;
t182 = -qJD(3) + t290;
t312 = qJD(5) * (qJD(3) + t182 + t290) + qJDD(5) * t196;
t253 = g(1) * t210 + g(2) * t209;
t255 = pkin(5) * t231 + pkin(8) * t228;
t311 = t198 * (pkin(8) * qJD(6) + t255 * qJD(1)) + t231 * t253 - g(3) * t228 + t162;
t227 = sin(qJ(6));
t230 = cos(qJ(6));
t279 = t230 * qJD(5);
t265 = t228 * t279;
t281 = qJD(6) * t231;
t243 = -t227 * t281 - t265;
t269 = t227 * qJDD(5) + t315 * t230;
t166 = qJD(1) * t243 + t269;
t307 = t166 * t227;
t280 = t227 * qJD(5);
t289 = qJD(1) * t231;
t188 = t230 * t289 + t280;
t266 = t228 * t280;
t294 = -qJD(1) * t266 - t230 * qJDD(5);
t167 = qJD(6) * t188 + t227 * t273 + t294;
t306 = t167 * t228;
t305 = t180 * t184;
t186 = t227 * t289 - t279;
t304 = t186 * t198;
t303 = t186 * t231;
t302 = t188 * t198;
t301 = t188 * t231;
t300 = t227 * t228;
t299 = t228 * t230;
t298 = t230 * t198;
t297 = t230 * t231;
t296 = qJDD(2) - g(3);
t295 = t166 * t228 + t188 * t284;
t222 = t231 ^ 2;
t292 = t228 ^ 2 - t222;
t233 = qJD(5) ^ 2;
t234 = qJD(1) ^ 2;
t291 = -t233 - t234;
t286 = qJD(5) * t186;
t285 = qJD(5) * t228;
t169 = qJD(5) * pkin(8) + t172;
t283 = qJD(6) * t169;
t282 = qJD(6) * t198;
t277 = qJD(4) * qJD(1);
t275 = qJDD(1) * t197;
t271 = t198 * t300;
t270 = t228 * t298;
t232 = cos(qJ(1));
t268 = t232 * pkin(1) + t210 * pkin(2) + t209 * qJ(3);
t263 = -qJDD(1) * pkin(2) - t226 * t308 + qJDD(3);
t229 = sin(qJ(1));
t262 = -pkin(1) * t229 + t210 * qJ(3);
t185 = qJD(5) * t255 + qJD(4);
t217 = qJDD(1) * qJ(4);
t257 = t217 - t263;
t164 = qJD(1) * t185 + qJDD(1) * t254 + t257;
t258 = -t164 + t283;
t252 = g(1) * t229 - g(2) * t232;
t251 = -qJD(5) * t181 - t296;
t249 = g(3) * t231 + t259;
t248 = t184 * t227 + t230 * t282;
t247 = -t184 * t230 + t227 * t282;
t245 = -t253 + t293;
t244 = t263 + t314;
t242 = -t217 + t244;
t241 = qJDD(1) * t201 + 0.2e1 * t219 + t245;
t240 = (-t248 + t286) * MDP(23);
t239 = -pkin(8) * t184 + (t168 + t171) * t198;
t238 = qJD(1) * t182 + qJD(2) * qJD(5) - t177 + t253;
t178 = t257 + t277;
t237 = -t196 * t233 + t178 + t275 + t277 - t314;
t192 = qJDD(5) * t231 - t228 * t233;
t191 = -qJDD(5) * t228 - t231 * t233;
t183 = t198 * t266;
t176 = -t209 * t227 + t210 * t299;
t175 = -t209 * t230 - t210 * t300;
t174 = -t209 * t299 - t210 * t227;
t173 = t209 * t300 - t210 * t230;
t163 = t230 * t164;
t160 = t169 * t230 + t170 * t227;
t159 = -t169 * t227 + t170 * t230;
t1 = [qJDD(1) * MDP(1) + t252 * MDP(2) + (g(1) * t232 + g(2) * t229) * MDP(3) + (t252 + (t225 ^ 2 + t226 ^ 2) * t308) * pkin(1) * MDP(4) + (qJDD(1) * t204 + t244) * MDP(5) + t241 * MDP(6) + (-t267 * t201 + t193 * qJD(3) + t263 * t204 - g(1) * (-pkin(2) * t209 + t262) - g(2) * t268) * MDP(7) + (qJDD(4) + t241) * MDP(8) + (-t242 + t275 + 0.2e1 * t277) * MDP(9) + (t178 * t197 + t182 * qJD(4) + t256 * t201 + t190 * qJD(3) - g(1) * ((-pkin(2) - qJ(4)) * t209 + t262) - g(2) * (qJ(4) * t210 + t268)) * MDP(10) + (qJDD(1) * t222 - 0.2e1 * t228 * t264) * MDP(11) + 0.2e1 * (-t228 * t273 + t292 * t278) * MDP(12) + t192 * MDP(13) + t191 * MDP(14) + (t237 * t228 + t231 * t312) * MDP(16) + (-t228 * t312 + t237 * t231) * MDP(17) + (t166 * t297 + t188 * t243) * MDP(18) + ((t186 * t230 + t188 * t227) * t285 + (-t307 - t167 * t230 + (t186 * t227 - t188 * t230) * qJD(6)) * t231) * MDP(19) + (t184 * t297 + t198 * t243 + t295) * MDP(20) + (-t306 + t183 + (-t248 - t286) * t231) * MDP(21) + (t184 * t228 + t198 * t284) * MDP(22) + (-g(1) * t174 - g(2) * t176 + (t196 * t286 + t163) * t228 + (-qJD(3) * t186 + t159 * qJD(5) - t196 * t167) * t231 + (t305 + t185 * t198 + (t168 * t231 + (-t196 * t198 - t169) * t228) * qJD(6)) * t230 + t313 * t227) * MDP(23) + (t196 * t188 * t285 - g(1) * t173 - g(2) * t175 + (-qJD(3) * t188 - t160 * qJD(5) - t196 * t166) * t231 + (-(-qJD(6) * t196 * t228 + t185) * t198 - t305 + t258 * t228 - t168 * t281) * t227 + t313 * t230) * MDP(24); t191 * MDP(16) - t192 * MDP(17) + (t183 + t306) * MDP(23) + (t198 * t265 + t295) * MDP(24) + (MDP(24) * t247 + t240) * t231 + (MDP(4) + MDP(7) + MDP(10)) * t296; t244 * MDP(7) + t242 * MDP(10) + t247 * MDP(23) + t248 * MDP(24) + (-MDP(6) - MDP(8)) * t234 + (-t228 * MDP(16) - MDP(17) * t231 + MDP(5) - MDP(9)) * qJDD(1) + (-t193 * MDP(7) + (-qJD(4) - t190) * MDP(10) + (t271 + t303) * MDP(23) + (t270 + t301) * MDP(24) + 0.2e1 * (-MDP(16) * t231 + MDP(17) * t228) * qJD(5)) * qJD(1); qJDD(1) * MDP(8) - t234 * MDP(9) + (qJDD(4) + t219 + t245) * MDP(10) + (-t182 * MDP(10) + (-MDP(23) * t230 + MDP(24) * t227) * t198) * qJD(1) + (qJDD(5) * MDP(16) + t291 * MDP(17) + (-t198 * t280 - t167) * MDP(23) + (-t198 * t279 - t166) * MDP(24)) * t231 + (t291 * MDP(16) - qJDD(5) * MDP(17) + t240 + (qJD(5) * t188 + t247) * MDP(24)) * t228; MDP(13) * t273 - MDP(14) * t274 + qJDD(5) * MDP(15) + (t228 * t251 - t231 * t238 + t287) * MDP(16) + (t228 * t238 + t231 * t251 + t288) * MDP(17) + (t188 * t298 + t307) * MDP(18) + ((t166 - t304) * t230 + (-t167 - t302) * t227) * MDP(19) + ((t270 - t301) * qJD(1) + t248) * MDP(20) + ((-t271 + t303) * qJD(1) - t247) * MDP(21) - t198 * MDP(22) * t289 + (-pkin(5) * t167 - t159 * t289 - t172 * t186 + t239 * t227 - t230 * t311) * MDP(23) + (-pkin(5) * t166 + t160 * t289 - t172 * t188 + t227 * t311 + t239 * t230) * MDP(24) + (t231 * t228 * MDP(11) - t292 * MDP(12)) * t234; t188 * t186 * MDP(18) + (-t186 ^ 2 + t188 ^ 2) * MDP(19) + (t269 + t304) * MDP(20) + (-t294 + t302) * MDP(21) + t184 * MDP(22) + (-g(1) * t175 + g(2) * t173 + t160 * t198 - t168 * t188 + t163) * MDP(23) + (g(1) * t176 - g(2) * t174 + t159 * t198 + t168 * t186) * MDP(24) + (-MDP(23) * t283 + t249 * MDP(24) + (-MDP(20) * t285 - MDP(21) * t281) * qJD(1)) * t230 + (-qJD(1) * MDP(20) * t281 - t315 * MDP(21) + t249 * MDP(23) + t258 * MDP(24)) * t227;];
tau  = t1;
