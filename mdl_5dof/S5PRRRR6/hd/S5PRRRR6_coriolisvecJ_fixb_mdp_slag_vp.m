% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRRRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:16
% EndTime: 2019-12-05 17:10:20
% DurationCPUTime: 1.22s
% Computational Cost: add. (967->165), mult. (1821->241), div. (0->0), fcn. (1280->8), ass. (0->103)
t230 = qJD(2) + qJD(3);
t233 = sin(qJ(5));
t234 = sin(qJ(4));
t237 = cos(qJ(5));
t238 = cos(qJ(4));
t203 = t233 * t238 + t234 * t237;
t229 = qJD(4) + qJD(5);
t304 = t203 * t229;
t168 = t304 * t230;
t308 = pkin(7) + pkin(8);
t236 = sin(qJ(2));
t276 = qJD(1) * qJD(2);
t281 = qJD(1) * t236;
t307 = qJD(3) * t281 + t236 * t276;
t301 = qJD(5) - t229;
t306 = MDP(8) * t234;
t288 = t237 * t238;
t292 = t233 * t234;
t201 = -t288 + t292;
t248 = t201 * t229;
t235 = sin(qJ(3));
t239 = cos(qJ(3));
t240 = cos(qJ(2));
t202 = t235 * t236 - t239 * t240;
t305 = t202 * t230;
t204 = t235 * t240 + t236 * t239;
t279 = qJD(3) * t235;
t260 = pkin(2) * t279 - t204 * qJD(1);
t252 = t260 * t230;
t303 = (t234 ^ 2 - t238 ^ 2) * MDP(9);
t280 = qJD(1) * t240;
t221 = qJD(2) * pkin(2) + t280;
t266 = t240 * t276;
t172 = t221 * t279 + t235 * t266 + t307 * t239;
t302 = -t252 - t172;
t283 = t307 * t235;
t171 = t239 * (qJD(3) * t221 + t266) - t283;
t299 = pkin(2) * t239;
t298 = pkin(3) * t230;
t224 = pkin(2) * t235 + pkin(7);
t297 = -pkin(8) - t224;
t182 = t230 * t204;
t296 = t182 * t230;
t191 = t221 * t235 + t239 * t281;
t295 = t191 * t230;
t294 = t204 * t229;
t241 = qJD(4) ^ 2;
t293 = t224 * t241;
t263 = t308 * t230 + t191;
t175 = t263 * t238;
t289 = t237 * t175;
t278 = qJD(4) * t234;
t272 = pkin(4) * t278;
t166 = t230 * t272 + t172;
t222 = t235 * t281;
t190 = t221 * t239 - t222;
t226 = -pkin(4) * t238 - pkin(3);
t178 = t226 * t230 - t190;
t287 = t166 * t203 - t178 * t248;
t286 = t166 * t201 + t178 * t304;
t187 = -t190 - t298;
t277 = qJD(4) * t238;
t285 = t172 * t234 + t187 * t277;
t284 = t272 + t260;
t275 = pkin(4) * t230 * t234;
t273 = qJD(3) * t299;
t271 = t230 * t292;
t270 = t230 * t288;
t269 = qJD(4) * t308;
t268 = t230 * t277;
t174 = t263 * t234;
t173 = qJD(4) * pkin(4) - t174;
t264 = -pkin(4) * t229 - t173;
t262 = qJD(4) * t297;
t259 = -t191 + t272;
t258 = pkin(7) * t241 - t295;
t257 = qJD(4) * (t190 - t298);
t256 = qJD(4) * t263;
t254 = t204 * t241 + t296;
t253 = 0.2e1 * qJD(4) * t305;
t157 = t171 * t238 - t234 * t256;
t158 = -t171 * t234 - t238 * t256;
t194 = t203 * t230;
t251 = -t233 * t157 + t237 * t158 - t178 * t194;
t167 = qJD(5) * t270 - t229 * t271 + t237 * t268;
t192 = -t270 + t271;
t247 = t194 * t192 * MDP(15) + (t192 * t229 + t167) * MDP(17) + (t194 * t229 - t168) * MDP(18) + (-t192 ^ 2 + t194 ^ 2) * MDP(16);
t196 = t239 * t280 - t222;
t246 = qJD(4) * ((-pkin(3) - t299) * t230 + t196 - t273);
t245 = t178 * t192 + (t301 * t175 - t158) * t233;
t244 = (-t167 * t201 - t168 * t203 + t192 * t248 - t194 * t304) * MDP(16) + (t167 * t203 - t194 * t248) * MDP(15) - 0.2e1 * t230 * qJD(4) * t303 + 0.2e1 * t268 * t306 + (t238 * MDP(10) - t234 * MDP(11)) * t241 + (-MDP(17) * t248 - MDP(18) * t304) * t229;
t227 = t238 * pkin(8);
t215 = pkin(7) * t238 + t227;
t214 = t308 * t234;
t210 = t226 - t299;
t206 = t238 * t269;
t205 = t234 * t269;
t198 = t224 * t238 + t227;
t197 = t297 * t234;
t186 = -t234 * t273 + t238 * t262;
t185 = t234 * t262 + t238 * t273;
t183 = t187 * t278;
t1 = [-MDP(6) * t296 + t305 * t230 * MDP(7) + (t234 * t253 - t238 * t254) * MDP(13) + (t234 * t254 + t238 * t253) * MDP(14) + (t202 * t168 + t182 * t192 + t294 * t248 + t304 * t305) * MDP(20) + (t202 * t167 + t182 * t194 - t248 * t305 + t294 * t304) * MDP(21) + (-t236 * MDP(3) - t240 * MDP(4)) * qJD(2) ^ 2; t302 * MDP(6) + (t196 * t230 + (-t266 + (-pkin(2) * t230 - t221) * qJD(3)) * t239 + t283) * MDP(7) + (t183 + t234 * t246 + (-t293 + t302) * t238) * MDP(13) + ((t293 + t252) * t234 + t238 * t246 + t285) * MDP(14) + ((-t185 * t233 + t237 * t186 + (-t197 * t233 - t198 * t237) * qJD(5)) * t229 + t210 * t168 + t196 * t304 + t284 * t192 + t286) * MDP(20) + (-(t185 * t237 + t186 * t233 + (t197 * t237 - t198 * t233) * qJD(5)) * t229 + t210 * t167 - t196 * t248 + t284 * t194 + t287) * MDP(21) + t244; (-t172 + t295) * MDP(6) + (t190 * t230 - t171) * MDP(7) + (t183 + t234 * t257 + (-t172 - t258) * t238) * MDP(13) + (t258 * t234 + t238 * t257 + t285) * MDP(14) + ((t205 * t233 - t206 * t237 + (t214 * t233 - t215 * t237) * qJD(5)) * t229 + t226 * t168 + t259 * t192 + t190 * t304 + t286) * MDP(20) + (-(-t205 * t237 - t206 * t233 + (-t214 * t237 - t215 * t233) * qJD(5)) * t229 + t226 * t167 + t259 * t194 - t190 * t248 + t287) * MDP(21) + t244; (-(t174 * t233 - t289) * t229 - t192 * t275 + (t264 * t233 - t289) * qJD(5) + t251) * MDP(20) + (-t194 * t275 + (t264 * qJD(5) - t174 * t229 - t157) * t237 + t245) * MDP(21) + t247 + (t234 * MDP(13) + t238 * MDP(14)) * (-t187 * t230 - t171) + (-t238 * t306 + t303) * t230 ^ 2; (t251 + t301 * (-t173 * t233 - t289)) * MDP(20) + ((-t301 * t173 - t157) * t237 + t245) * MDP(21) + t247;];
tauc = t1;
