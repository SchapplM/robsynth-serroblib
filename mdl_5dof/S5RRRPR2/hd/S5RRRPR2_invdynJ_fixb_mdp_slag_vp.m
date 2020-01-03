% Calculate vector of inverse dynamics joint torques for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:40
% EndTime: 2020-01-03 12:07:43
% DurationCPUTime: 1.17s
% Computational Cost: add. (1290->199), mult. (2175->265), div. (0->0), fcn. (1241->16), ass. (0->124)
t238 = qJD(1) + qJD(2);
t250 = cos(qJ(2));
t310 = pkin(1) * qJD(1);
t206 = pkin(2) * t238 + t250 * t310;
t249 = cos(qJ(3));
t293 = qJD(2) * t250;
t282 = qJD(1) * t293;
t246 = sin(qJ(2));
t286 = qJDD(1) * t246;
t263 = t282 + t286;
t259 = t263 * pkin(1);
t318 = (qJD(3) * t206 + t259) * t249;
t314 = pkin(1) * t250;
t229 = qJDD(1) * t314;
t237 = qJDD(1) + qJDD(2);
t285 = t246 * t310;
t192 = pkin(2) * t237 - qJD(2) * t285 + t229;
t189 = t249 * t192;
t232 = qJDD(3) + t237;
t245 = sin(qJ(3));
t289 = qJD(3) * t249;
t283 = t246 * t289;
t291 = qJD(3) * t245;
t284 = t206 * t291;
t165 = -t284 + pkin(3) * t232 + t189 + (-t245 * t286 + (-t245 * t293 - t283) * qJD(1)) * pkin(1);
t290 = qJD(3) * t246;
t281 = qJD(1) * t290;
t277 = pkin(1) * t281;
t209 = t245 * t277;
t307 = t192 * t245;
t171 = -t209 + t307 + t318;
t242 = sin(pkin(9));
t243 = cos(pkin(9));
t160 = t165 * t243 - t171 * t242;
t241 = qJ(1) + qJ(2);
t236 = qJ(3) + t241;
t224 = pkin(9) + t236;
t212 = sin(t224);
t213 = cos(t224);
t316 = -pkin(4) * t232 + g(2) * t213 + g(3) * t212 - t160;
t227 = pkin(2) * t249 + pkin(3);
t306 = t242 * t245;
t265 = -pkin(2) * t306 + t227 * t243;
t193 = -pkin(4) - t265;
t304 = t243 * t245;
t295 = pkin(2) * t304 + t242 * t227;
t194 = pkin(8) + t295;
t233 = qJD(3) + t238;
t252 = qJD(5) ^ 2;
t302 = t246 * t249;
t270 = -t245 * t250 - t302;
t198 = t270 * t310;
t303 = t245 * t246;
t269 = t249 * t250 - t303;
t199 = t269 * t310;
t309 = pkin(2) * qJD(3);
t299 = -t198 * t243 + t199 * t242 - (t242 * t249 + t304) * t309;
t315 = t193 * t232 + t194 * t252 - t299 * t233;
t313 = pkin(2) * t232;
t188 = t206 * t245 + t249 * t285;
t308 = t188 * t242;
t305 = t243 * t188;
t301 = qJDD(4) - g(1);
t161 = t242 * t165 + t243 * t171;
t228 = pkin(2) + t314;
t211 = t249 * t228;
t197 = -pkin(1) * t303 + pkin(3) + t211;
t200 = pkin(1) * t302 + t228 * t245;
t300 = t242 * t197 + t243 * t200;
t298 = -t198 * t242 - t199 * t243 + (t243 * t249 - t306) * t309;
t225 = sin(t236);
t234 = sin(t241);
t297 = pkin(2) * t234 + pkin(3) * t225;
t226 = cos(t236);
t235 = cos(t241);
t296 = pkin(2) * t235 + pkin(3) * t226;
t244 = sin(qJ(5));
t239 = t244 ^ 2;
t248 = cos(qJ(5));
t294 = -t248 ^ 2 + t239;
t288 = qJD(5) * t233;
t287 = qJD(5) * t248;
t280 = g(2) * t234 - g(3) * t235;
t279 = qJD(1) * (-qJD(2) + t238);
t278 = qJD(2) * (-qJD(1) - t238);
t187 = t249 * t206 - t245 * t285;
t184 = pkin(3) * t233 + t187;
t172 = t184 * t243 - t308;
t169 = -pkin(4) * t233 - t172;
t276 = t169 * t287 + t316 * t244;
t275 = g(2) * t225 - g(3) * t226 + t209;
t274 = -g(2) * t226 - g(3) * t225;
t207 = qJDD(5) * t244 + t248 * t252;
t208 = qJDD(5) * t248 - t244 * t252;
t273 = 0.2e1 * (t232 * t244 * t248 - t294 * t288) * MDP(12) + (0.2e1 * t233 * t244 * t287 + t232 * t239) * MDP(11) + t207 * MDP(13) + t208 * MDP(14) + t232 * MDP(7);
t271 = t197 * t243 - t200 * t242;
t267 = t189 + t274;
t266 = -g(2) * t235 - g(3) * t234 + t229;
t264 = t237 * MDP(4) + t273;
t262 = t275 - t307;
t180 = t228 * t289 + (t269 * qJD(2) - t245 * t290) * pkin(1);
t181 = -t228 * t291 + (t270 * qJD(2) - t283) * pkin(1);
t162 = t180 * t242 - t181 * t243;
t176 = -pkin(4) - t271;
t177 = pkin(8) + t300;
t261 = t162 * t233 + t176 * t232 + t177 * t252;
t174 = t187 * t242 + t305;
t219 = pkin(3) * t242 + pkin(8);
t220 = -pkin(3) * t243 - pkin(4);
t260 = -t174 * t233 + t219 * t252 + t220 * t232;
t258 = -pkin(8) * t232 + g(2) * t212 - g(3) * t213 - t169 * t233 - t161;
t163 = t180 * t243 + t181 * t242;
t257 = -qJDD(5) * t177 + (t176 * t233 - t163) * qJD(5);
t175 = t187 * t243 - t308;
t256 = qJD(5) * t175 - qJDD(5) * t219 + t220 * t288;
t255 = t267 - t284;
t254 = -qJDD(5) * t194 + (t193 * t233 - t298) * qJD(5);
t253 = (-pkin(2) * t233 - t206) * qJD(3) - t259;
t251 = cos(qJ(1));
t247 = sin(qJ(1));
t231 = t233 ^ 2;
t173 = t242 * t184 + t305;
t166 = t169 * qJD(5) * t244;
t1 = [qJDD(1) * MDP(1) + (-g(2) * t251 - g(3) * t247) * MDP(2) + (g(2) * t247 - g(3) * t251) * MDP(3) + ((t237 * t250 + t246 * t278) * pkin(1) + t266) * MDP(5) + (((-qJDD(1) - t237) * t246 + t250 * t278) * pkin(1) + t280) * MDP(6) + (t181 * t233 + t211 * t232 + (-t249 * t281 + (-t282 + (-qJDD(1) - t232) * t246) * t245) * pkin(1) + t255) * MDP(8) + (-t180 * t233 - t200 * t232 + t262 - t318) * MDP(9) + (t161 * t300 + t173 * t163 + t160 * t271 - t172 * t162 - g(2) * (pkin(1) * t251 + t296) - g(3) * (pkin(1) * t247 + t297)) * MDP(10) + (t166 + t257 * t244 + (-t261 - t316) * t248) * MDP(16) + (t261 * t244 + t257 * t248 + t276) * MDP(17) + t264; (t246 * pkin(1) * t279 + t266) * MDP(5) + ((t250 * t279 - t286) * pkin(1) + t280) * MDP(6) + (-t198 * t233 + (-t277 + t313) * t249 + t253 * t245 + t267) * MDP(8) + (t199 * t233 + (-t192 - t313) * t245 + t253 * t249 + t275) * MDP(9) + (-g(2) * t296 - g(3) * t297 + t160 * t265 + t161 * t295 + t299 * t172 + t298 * t173) * MDP(10) + (t166 + t254 * t244 + (-t316 - t315) * t248) * MDP(16) + (t315 * t244 + t254 * t248 + t276) * MDP(17) + t264; (t188 * t233 + t255) * MDP(8) + (t187 * t233 - t206 * t289 + t262) * MDP(9) + t166 * MDP(16) + t276 * MDP(17) + (-t263 * MDP(8) * t245 + (-MDP(8) * t281 - t263 * MDP(9)) * t249) * pkin(1) + (t256 * MDP(16) + t260 * MDP(17)) * t244 + ((-t260 - t316) * MDP(16) + t256 * MDP(17)) * t248 + t273 + (t172 * t174 - t173 * t175 + (t160 * t243 + t161 * t242 + t274) * pkin(3)) * MDP(10); t301 * MDP(10) + t208 * MDP(16) - t207 * MDP(17); qJDD(5) * MDP(15) + t294 * MDP(12) * t231 + (t232 * MDP(14) + t301 * MDP(16) + t258 * MDP(17)) * t248 + (-t231 * t248 * MDP(11) + t232 * MDP(13) + t258 * MDP(16) - t301 * MDP(17)) * t244;];
tau = t1;
