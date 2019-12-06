% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:06:58
% EndTime: 2019-12-05 18:07:04
% DurationCPUTime: 1.69s
% Computational Cost: add. (1604->210), mult. (4205->298), div. (0->0), fcn. (2876->6), ass. (0->122)
t256 = sin(pkin(8));
t295 = qJD(3) + qJD(4);
t341 = t256 * t295;
t259 = sin(qJ(3));
t261 = cos(qJ(3));
t340 = (t259 * MDP(10) + t261 * MDP(11)) * t256;
t252 = t256 ^ 2;
t330 = 0.2e1 * t252;
t258 = sin(qJ(4));
t260 = cos(qJ(4));
t231 = t258 * t261 + t259 * t260;
t222 = t231 * t256;
t339 = (t259 ^ 2 - t261 ^ 2) * MDP(9);
t318 = t260 * t261;
t319 = t258 * t259;
t230 = t318 - t319;
t257 = cos(pkin(8));
t305 = qJD(1) * t257;
t315 = (t295 - t305) * t230;
t210 = t295 * t231;
t266 = qJD(1) * t231;
t338 = t257 * t266 - t210;
t306 = qJD(1) * t256;
t291 = t259 * t306;
t278 = t258 * t291;
t290 = t261 * t306;
t279 = t260 * t290;
t218 = -t278 + t279;
t212 = t218 * qJ(5);
t234 = -pkin(2) * t257 - pkin(6) * t256 - pkin(1);
t225 = qJD(1) * t234 + qJD(2);
t221 = t261 * t225;
t320 = t257 * t259;
t329 = pkin(7) * t256;
t265 = -qJ(2) * t320 - t261 * t329;
t206 = qJD(1) * t265 + t221;
t245 = -qJD(3) + t305;
t192 = -pkin(3) * t245 + t206;
t326 = qJ(2) * t261;
t293 = t257 * t326;
t207 = -pkin(7) * t291 + qJD(1) * t293 + t225 * t259;
t199 = t258 * t207;
t285 = t192 * t260 - t199;
t337 = t212 - t285;
t253 = t257 ^ 2;
t336 = t330 + t253;
t335 = -MDP(13) * t259 - MDP(14) * t261;
t287 = -t234 + t329;
t332 = t259 * t287 - t293;
t331 = t218 ^ 2;
t328 = MDP(7) * qJ(2);
t327 = qJ(2) * t259;
t215 = t256 * t266;
t325 = qJ(5) * t215;
t226 = pkin(3) * t291 + qJ(2) * t306;
t204 = pkin(4) * t215 + qJD(5) + t226;
t324 = t204 * t218;
t323 = t226 * t218;
t262 = qJD(1) ^ 2;
t322 = t252 * t262;
t321 = t256 * t259;
t201 = t260 * t207;
t240 = -qJD(4) + t245;
t173 = -pkin(4) * t240 - t337;
t317 = t173 + t337;
t316 = t206 * t260 - t199;
t297 = qJD(1) * qJD(2);
t288 = t257 * t297;
t301 = qJD(3) * t261;
t313 = t225 * t301 + t261 * t288;
t303 = qJD(2) * t261;
t312 = t234 * t301 + t257 * t303;
t311 = t295 * t278;
t243 = t256 * pkin(3) * t301;
t224 = qJD(1) * t243 + t256 * t297;
t228 = qJD(2) * t256 + t243;
t229 = pkin(3) * t321 + qJ(2) * t256;
t310 = t252 + t253;
t304 = qJD(2) * t259;
t302 = qJD(3) * t259;
t300 = qJD(4) * t258;
t299 = qJD(4) * t260;
t298 = qJD(3) + t245;
t292 = qJ(2) * t302;
t289 = t257 * t304;
t264 = t265 * qJD(3);
t187 = qJD(1) * t264 + t313;
t188 = -t225 * t302 + (-t289 + (pkin(7) * t321 - t293) * qJD(3)) * qJD(1);
t286 = -t258 * t187 + t188 * t260;
t202 = t264 + t312;
t203 = qJD(3) * t332 - t289;
t284 = -t202 * t258 + t203 * t260;
t283 = -t206 * t258 - t201;
t282 = qJD(1) * t298;
t281 = t261 * t252 * t259 * MDP(8);
t280 = -t187 * t260 - t188 * t258 - t192 * t299 + t207 * t300;
t272 = t318 * t341;
t191 = qJD(1) * t272 - t311;
t277 = pkin(4) * t191 + t224;
t275 = -t192 * t258 - t201;
t208 = -t287 * t261 + (-pkin(3) - t327) * t257;
t274 = -t208 * t258 + t260 * t332;
t269 = t226 * t215 + t280;
t213 = t215 ^ 2;
t268 = t218 * t215 * MDP(15) + (-t210 * t306 - t215 * t240) * MDP(17) + (-t218 * t240 - t279 * t295 + t311) * MDP(18) + (-t213 + t331) * MDP(16);
t267 = t202 * t260 + t203 * t258 + t208 * t299 + t300 * t332;
t263 = qJD(4) * t275 + t286;
t197 = t295 * t222;
t249 = pkin(3) * t260 + pkin(4);
t223 = t230 * t256;
t198 = -t319 * t341 + t272;
t190 = qJD(1) * t197;
t179 = -qJ(5) * t222 - t274;
t178 = -pkin(4) * t257 - qJ(5) * t223 + t208 * t260 + t258 * t332;
t177 = -t212 + t316;
t176 = t283 + t325;
t175 = -t275 - t325;
t172 = qJ(5) * t197 + qJD(4) * t274 - qJD(5) * t223 + t284;
t171 = -qJ(5) * t198 - qJD(5) * t222 + t267;
t170 = qJ(5) * t190 - qJD(5) * t218 + t263;
t169 = -qJ(5) * t191 - qJD(5) * t215 - t280;
t1 = [(t245 * t289 + t336 * qJD(1) * (qJ(2) * t301 + t304)) * MDP(13) + ((-t257 * t292 + t312) * t245 + t313 * t257 + (-t292 * t336 + t303 * t330) * qJD(1)) * MDP(14) + (-t190 * t223 - t197 * t218) * MDP(15) + (t190 * t222 - t191 * t223 + t197 * t215 - t198 * t218) * MDP(16) + (t190 * t257 + t197 * t240) * MDP(17) + (t191 * t257 + t198 * t240) * MDP(18) + (-t284 * t240 - t286 * t257 + t228 * t215 + t229 * t191 + t224 * t222 + t226 * t198 + (-t240 * t274 - t257 * t275) * qJD(4)) * MDP(20) + (-t229 * t190 - t226 * t197 + t228 * t218 + t224 * t223 + t240 * t267 - t257 * t280) * MDP(21) + (-t169 * t222 - t170 * t223 - t171 * t215 - t172 * t218 + t173 * t197 - t175 * t198 + t178 * t190 - t179 * t191) * MDP(22) + (t169 * t179 + t175 * t171 + t170 * t178 + t173 * t172 + t277 * (pkin(4) * t222 + t229) + t204 * (pkin(4) * t198 + t228)) * MDP(23) + 0.2e1 * (MDP(6) + t328) * t310 * t297 + ((-(-t234 * t259 - t293) * t245 + t225 * t320) * MDP(13) + (t330 * t339 - 0.2e1 * t281) * qJD(1) + (t245 + t305) * t340) * qJD(3); (t190 * t230 - t191 * t231 - t215 * t315 - t218 * t338) * MDP(22) + (t169 * t231 + t170 * t230 + t173 * t338 + t175 * t315) * MDP(23) + (-MDP(20) * t215 - MDP(21) * t218 - MDP(23) * t204) * t306 + (-MDP(20) * t338 + MDP(21) * t315) * t240 + (-t253 * MDP(6) + (-MDP(6) + t335) * t252 - t310 * t328) * t262 + t335 * t245 ^ 2; t262 * t281 - t322 * t339 + ((-t225 * t298 - t288) * t259 + (-t257 * t282 - t322) * t326) * MDP(13) + (-t221 * t245 + (t298 * t305 + t322) * t327 - t313) * MDP(14) + (t283 * t240 - pkin(3) * t215 * t290 - t323 + (-t201 + (pkin(3) * t240 - t192) * t258) * qJD(4) + t286) * MDP(20) + (-t316 * t240 + (-t218 * t290 + t240 * t299) * pkin(3) + t269) * MDP(21) + (t190 * t249 + (t175 + t176) * t218 + (-t173 + t177) * t215 + (-t191 * t258 + (-t215 * t260 + t218 * t258) * qJD(4)) * pkin(3)) * MDP(22) + (-pkin(4) * t324 + t170 * t249 - t173 * t176 - t175 * t177 + (-t204 * t290 + t169 * t258 + (-t173 * t258 + t175 * t260) * qJD(4)) * pkin(3)) * MDP(23) + t268 - t282 * t340; (t240 * t275 + t263 - t323) * MDP(20) + (-t240 * t285 + t269) * MDP(21) + (pkin(4) * t190 - t215 * t317) * MDP(22) + (t317 * t175 + (t170 - t324) * pkin(4)) * MDP(23) + t268; (-t213 - t331) * MDP(22) + (t173 * t218 + t175 * t215 + t277) * MDP(23);];
tauc = t1;
