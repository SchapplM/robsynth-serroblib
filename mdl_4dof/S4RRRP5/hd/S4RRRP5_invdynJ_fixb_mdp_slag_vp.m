% Calculate vector of inverse dynamics joint torques for
% S4RRRP5
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:09
% EndTime: 2019-12-31 17:17:11
% DurationCPUTime: 1.50s
% Computational Cost: add. (1209->237), mult. (2655->301), div. (0->0), fcn. (1674->8), ass. (0->117)
t236 = qJD(2) + qJD(3);
t304 = pkin(6) + pkin(5);
t239 = qJ(2) + qJ(3);
t233 = sin(t239);
t234 = cos(t239);
t290 = t234 * pkin(3) + t233 * qJ(4);
t235 = qJDD(2) + qJDD(3);
t230 = t235 * qJ(4);
t231 = t236 * qJD(4);
t308 = t230 + t231;
t240 = sin(qJ(3));
t241 = sin(qJ(2));
t243 = cos(qJ(2));
t303 = cos(qJ(3));
t197 = t240 * t243 + t303 * t241;
t193 = t197 * qJD(1);
t232 = t235 * pkin(3);
t307 = qJDD(4) - t232;
t281 = qJD(2) * t304;
t199 = t241 * t281;
t201 = t243 * t281;
t203 = t304 * t241;
t204 = t304 * t243;
t260 = -t303 * t203 - t240 * t204;
t161 = t260 * qJD(3) - t303 * t199 - t240 * t201;
t185 = -t240 * t203 + t303 * t204;
t242 = sin(qJ(1));
t244 = cos(qJ(1));
t266 = g(1) * t242 - g(2) * t244;
t306 = t161 * t236 + t185 * t235 + t266 * t233;
t305 = t193 ^ 2;
t302 = pkin(2) * t243;
t301 = qJD(2) * pkin(2);
t198 = qJD(1) * t203;
t195 = -t198 + t301;
t200 = qJD(1) * t204;
t293 = t240 * t200;
t172 = t303 * t195 - t293;
t300 = t172 * t236;
t280 = t303 * t200;
t173 = t240 * t195 + t280;
t299 = t173 * t236;
t279 = t303 * t243;
t270 = qJD(1) * t279;
t288 = qJD(1) * t241;
t278 = t240 * t288;
t191 = -t270 + t278;
t298 = t191 * t193;
t297 = t233 * t242;
t296 = t233 * t244;
t295 = t234 * t242;
t294 = t234 * t244;
t292 = t240 * t241;
t177 = -t303 * t198 - t293;
t277 = qJD(3) * t303;
t291 = pkin(2) * t277 + qJD(4) - t177;
t237 = t241 ^ 2;
t289 = -t243 ^ 2 + t237;
t287 = qJD(3) * t240;
t286 = qJD(4) - t172;
t285 = qJD(1) * qJD(2);
t284 = qJDD(1) * t241;
t283 = qJDD(1) * t243;
t282 = t241 * t301;
t229 = pkin(1) + t302;
t276 = t241 * t285;
t275 = t243 * t285;
t274 = qJDD(1) * t303;
t182 = qJDD(2) * pkin(2) + t304 * (-t275 - t284);
t183 = t304 * (-t276 + t283);
t273 = t240 * t182 + t303 * t183 + t195 * t277 - t200 * t287;
t272 = -t303 * t182 + t240 * t183 + t195 * t287 + t200 * t277;
t271 = -t236 * t270 - t240 * t283 - t241 * t274;
t176 = -t240 * t198 + t280;
t269 = pkin(2) * t287 - t176;
t268 = -pkin(2) * t241 - pkin(3) * t233;
t267 = g(1) * t244 + g(2) * t242;
t202 = t229 * qJD(1);
t265 = t236 * t292;
t170 = pkin(3) * t193 + qJ(4) * t191;
t264 = t240 * t284 - t243 * t274;
t262 = t229 + t290;
t261 = -0.2e1 * pkin(1) * t285 - pkin(5) * qJDD(2);
t188 = pkin(2) * t276 - t229 * qJDD(1);
t257 = -g(1) * t294 - g(2) * t295 - g(3) * t233 + t273;
t160 = -t271 + (t191 - t278) * t236;
t256 = MDP(11) * t298 + t160 * MDP(13) - t264 * MDP(14) + (-t191 ^ 2 + t305) * MDP(12) + t235 * MDP(15);
t255 = g(1) * t296 + g(2) * t297 - g(3) * t234 - t272;
t246 = qJD(2) ^ 2;
t254 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t246 + t266;
t247 = qJD(1) ^ 2;
t253 = pkin(1) * t247 - pkin(5) * qJDD(1) + t267;
t162 = t185 * qJD(3) - t240 * t199 + t303 * t201;
t252 = g(1) * t295 - g(2) * t294 - t162 * t236 + t235 * t260;
t166 = pkin(3) * t191 - qJ(4) * t193 - t202;
t251 = -t166 * t191 + t257;
t250 = -t202 * t191 - t257;
t249 = t202 * t193 + t255;
t248 = t166 * t193 - t255 + t307;
t179 = t236 * t197;
t228 = -t303 * pkin(2) - pkin(3);
t223 = pkin(2) * t240 + qJ(4);
t206 = qJ(4) * t294;
t205 = qJ(4) * t295;
t196 = -t279 + t292;
t178 = -qJD(2) * t279 - t243 * t277 + t265;
t171 = pkin(3) * t196 - qJ(4) * t197 - t229;
t169 = t236 * qJ(4) + t173;
t168 = pkin(2) * t288 + t170;
t167 = -t236 * pkin(3) + t286;
t164 = t179 * qJD(1) + t264;
t163 = qJD(1) * t265 + t271;
t157 = pkin(3) * t179 + qJ(4) * t178 - qJD(4) * t197 + t282;
t156 = t272 + t307;
t155 = t273 + t308;
t154 = pkin(3) * t164 + qJ(4) * t163 - qJD(4) * t193 + t188;
t1 = [qJDD(1) * MDP(1) + t266 * MDP(2) + t267 * MDP(3) + (qJDD(1) * t237 + 0.2e1 * t241 * t275) * MDP(4) + 0.2e1 * (t241 * t283 - t289 * t285) * MDP(5) + (qJDD(2) * t241 + t243 * t246) * MDP(6) + (qJDD(2) * t243 - t241 * t246) * MDP(7) + (t261 * t241 + t254 * t243) * MDP(9) + (-t254 * t241 + t261 * t243) * MDP(10) + (-t163 * t197 - t178 * t193) * MDP(11) + (t163 * t196 - t164 * t197 + t178 * t191 - t179 * t193) * MDP(12) + (-t178 * t236 + t197 * t235) * MDP(13) + (-t179 * t236 - t196 * t235) * MDP(14) + (-t164 * t229 - t179 * t202 + t188 * t196 + t191 * t282 + t252) * MDP(16) + (t163 * t229 + t178 * t202 + t188 * t197 + t193 * t282 - t306) * MDP(17) + (t154 * t196 + t157 * t191 + t164 * t171 + t166 * t179 + t252) * MDP(18) + (-t155 * t196 + t156 * t197 - t161 * t191 + t162 * t193 + t163 * t260 - t164 * t185 - t167 * t178 - t169 * t179 - t267) * MDP(19) + (-t154 * t197 - t157 * t193 + t163 * t171 + t166 * t178 + t306) * MDP(20) + (t154 * t171 + t155 * t185 - t156 * t260 + t166 * t157 + t169 * t161 + t167 * t162 + (-g(1) * t304 - g(2) * t262) * t244 + (g(1) * t262 - g(2) * t304) * t242) * MDP(21); MDP(6) * t284 + MDP(7) * t283 + qJDD(2) * MDP(8) + (-g(3) * t243 + t253 * t241) * MDP(9) + (g(3) * t241 + t253 * t243) * MDP(10) + (t176 * t236 + (-t191 * t288 + t303 * t235 - t236 * t287) * pkin(2) + t249) * MDP(16) + (t177 * t236 + (-t193 * t288 - t235 * t240 - t236 * t277) * pkin(2) + t250) * MDP(17) + (-t168 * t191 - t228 * t235 - t269 * t236 - t248) * MDP(18) + (-t163 * t228 - t164 * t223 + (t169 + t269) * t193 + (t167 - t291) * t191) * MDP(19) + (t168 * t193 + t223 * t235 + t291 * t236 + t251 + t308) * MDP(20) + (t155 * t223 + t156 * t228 - t166 * t168 - g(1) * (t268 * t244 + t206) - g(2) * (t268 * t242 + t205) - g(3) * (t290 + t302) + t291 * t169 + t269 * t167) * MDP(21) + t256 + (-t241 * t243 * MDP(4) + t289 * MDP(5)) * t247; (t249 + t299) * MDP(16) + (t250 + t300) * MDP(17) + (-t170 * t191 + t232 - t248 + t299) * MDP(18) + (pkin(3) * t163 - qJ(4) * t164 + (t169 - t173) * t193 + (t167 - t286) * t191) * MDP(19) + (t170 * t193 + 0.2e1 * t230 + 0.2e1 * t231 + t251 - t300) * MDP(20) + (t155 * qJ(4) - t156 * pkin(3) - t166 * t170 - t167 * t173 - g(1) * (-pkin(3) * t296 + t206) - g(2) * (-pkin(3) * t297 + t205) - g(3) * t290 + t286 * t169) * MDP(21) + t256; (-t235 + t298) * MDP(18) + t160 * MDP(19) + (-t236 ^ 2 - t305) * MDP(20) + (-t169 * t236 + t248) * MDP(21);];
tau = t1;
