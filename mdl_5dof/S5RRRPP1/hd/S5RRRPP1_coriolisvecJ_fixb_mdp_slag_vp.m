% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:50:01
% EndTime: 2019-12-31 20:50:04
% DurationCPUTime: 1.28s
% Computational Cost: add. (1842->220), mult. (3192->284), div. (0->0), fcn. (1904->6), ass. (0->121)
t340 = qJ(4) + pkin(7);
t272 = cos(qJ(3));
t262 = t272 * qJD(4);
t270 = sin(qJ(3));
t292 = qJD(3) * t340;
t234 = -t270 * t292 + t262;
t269 = cos(pkin(8));
t322 = t269 * t272;
t268 = sin(pkin(8));
t325 = t268 * t270;
t241 = -t322 + t325;
t278 = -qJD(4) * t270 - t272 * t292;
t273 = cos(qJ(2));
t332 = pkin(1) * qJD(1);
t303 = t273 * t332;
t315 = t269 * t234 + t241 * t303 + t268 * t278;
t242 = t268 * t272 + t269 * t270;
t265 = qJD(1) + qJD(2);
t231 = t242 * t265;
t339 = t231 * MDP(16);
t338 = MDP(7) * t270;
t337 = (t270 ^ 2 - t272 ^ 2) * MDP(8);
t335 = MDP(12) * t270 + MDP(13) * t272;
t226 = t231 ^ 2;
t334 = pkin(1) * t273;
t331 = pkin(1) * qJD(2);
t300 = qJD(1) * t331;
t288 = t273 * t300;
t279 = qJD(4) * t265 + t288;
t271 = sin(qJ(2));
t301 = t271 * t332;
t290 = t340 * t265 + t301;
t284 = qJD(3) * t290;
t203 = -t270 * t284 + t279 * t272;
t275 = -t279 * t270 - t272 * t284;
t175 = t203 * t268 - t269 * t275;
t258 = pkin(1) * t271 + pkin(7);
t263 = t272 * qJ(4);
t240 = t258 * t272 + t263;
t319 = -qJ(4) - t258;
t291 = t319 * t270;
t207 = t240 * t268 - t269 * t291;
t330 = t175 * t207;
t250 = pkin(7) * t272 + t263;
t293 = t340 * t270;
t215 = t250 * t268 + t269 * t293;
t329 = t175 * t215;
t328 = t175 * t242;
t219 = t290 * t272;
t327 = t219 * t268;
t326 = t265 * t270;
t213 = t269 * t219;
t274 = qJD(3) ^ 2;
t321 = t270 * t274;
t320 = t272 * t274;
t236 = t242 * qJD(3);
t222 = t265 * t236;
t310 = qJD(3) * t270;
t294 = t268 * t310;
t247 = t265 * t294;
t309 = qJD(3) * t272;
t295 = t269 * t309;
t223 = t265 * t295 - t247;
t253 = t271 * t300;
t297 = t265 * t310;
t235 = pkin(3) * t297 + t253;
t282 = pkin(4) * t222 - qJ(5) * t223 + t235;
t178 = -qJD(5) * t231 + t282;
t298 = -pkin(3) * t272 - pkin(2);
t228 = t298 * t265 + qJD(4) - t303;
t299 = t265 * t322;
t229 = t265 * t325 - t299;
t187 = pkin(4) * t229 - qJ(5) * t231 + t228;
t318 = t178 * t241 + t187 * t236;
t237 = -t294 + t295;
t317 = -t178 * t242 - t187 * t237;
t176 = t269 * t203 + t268 * t275;
t316 = t234 * t268 - t242 * t303 - t269 * t278;
t218 = t290 * t270;
t217 = qJD(3) * pkin(3) - t218;
t194 = t268 * t217 + t213;
t246 = -pkin(2) * t265 - t303;
t314 = t246 * t309 + t270 * t253;
t308 = qJD(3) * t273;
t307 = t272 * MDP(12);
t306 = -qJD(2) + t265;
t197 = -t218 * t269 - t327;
t305 = qJD(5) - t197;
t304 = qJD(3) * qJD(5);
t302 = t273 * t331;
t260 = pkin(3) * t310;
t296 = t265 * t309;
t289 = qJD(3) * t319;
t174 = t304 + t176;
t193 = t217 * t269 - t327;
t190 = -qJD(3) * pkin(4) + qJD(5) - t193;
t191 = qJD(3) * qJ(5) + t194;
t287 = -t174 * t241 + t190 * t237 - t191 * t236 + t328;
t286 = -t176 * t241 - t193 * t237 - t194 * t236 + t328;
t195 = pkin(4) * t236 - qJ(5) * t237 - qJD(5) * t242 + t260;
t285 = -t195 + t301;
t283 = -0.2e1 * t265 * qJD(3) * t337 - MDP(10) * t321 + MDP(9) * t320 + 0.2e1 * t296 * t338;
t211 = t270 * t289 + t272 * t302 + t262;
t276 = (-qJD(4) - t302) * t270 + t272 * t289;
t185 = t211 * t268 - t269 * t276;
t186 = t269 * t211 + t268 * t276;
t208 = t269 * t240 + t268 * t291;
t280 = t185 * t231 - t186 * t229 + t207 * t223 - t208 * t222;
t209 = pkin(4) * t241 - qJ(5) * t242 + t298;
t216 = t269 * t250 - t268 * t293;
t277 = t215 * t223 - t216 * t222 - t315 * t229 + t316 * t231;
t261 = t271 * t331;
t259 = -pkin(2) - t334;
t256 = -pkin(3) * t269 - pkin(4);
t254 = pkin(3) * t268 + qJ(5);
t238 = t246 * t310;
t206 = t209 - t334;
t198 = pkin(3) * t326 + pkin(4) * t231 + qJ(5) * t229;
t196 = -t218 * t268 + t213;
t188 = t195 + t261;
t1 = [-t253 * MDP(5) + (-t258 * t320 + t259 * t297 + t238) * MDP(12) + (t258 * t321 + t259 * t296 + t314) * MDP(13) + (t280 + t286) * MDP(14) + (t176 * t208 + t194 * t186 + t330 - t193 * t185 + t235 * (t298 - t334) + t228 * (t261 + t260)) * MDP(15) + (-qJD(3) * t185 + t188 * t229 + t206 * t222 + t318) * MDP(16) + (t280 + t287) * MDP(17) + (qJD(3) * t186 - t188 * t231 - t206 * t223 + t317) * MDP(18) + (t174 * t208 + t178 * t206 + t185 * t190 + t186 * t191 + t187 * t188 + t330) * MDP(19) + (((-qJD(1) - t265) * MDP(6) - t335 * qJD(3)) * t273 + (-qJD(1) * t307 + (t270 * MDP(13) - MDP(5) - t307) * t265) * t271) * t331 + t283; (t265 * t301 - t253) * MDP(5) + t306 * MDP(6) * t303 + (-pkin(2) * t297 - pkin(7) * t320 + t238 + (t306 * t272 * t271 + t270 * t308) * t332) * MDP(12) + (-pkin(2) * t296 + pkin(7) * t321 + (-t271 * t326 + t272 * t308) * t332 + t314) * MDP(13) + (t277 + t286) * MDP(14) + (t176 * t216 + t329 + t235 * t298 + (-t301 + t260) * t228 + t315 * t194 - t316 * t193) * MDP(15) + (-t316 * qJD(3) + t209 * t222 - t285 * t229 + t318) * MDP(16) + (t277 + t287) * MDP(17) + (t315 * qJD(3) - t209 * t223 + t285 * t231 + t317) * MDP(18) + (t174 * t216 + t178 * t209 - t285 * t187 + t316 * t190 + t315 * t191 + t329) * MDP(19) + t283; (t193 * t196 - t194 * t197) * MDP(15) + (qJD(3) * t196 - t175) * MDP(16) + (-t222 * t254 + t223 * t256) * MDP(17) + (-qJD(3) * t197 + t176 + 0.2e1 * t304) * MDP(18) + (t174 * t254 + t175 * t256 - t187 * t198 - t190 * t196 + t305 * t191) * MDP(19) + (-t272 * t338 + t337) * t265 ^ 2 + ((t194 - t196) * MDP(14) - t187 * MDP(16) + (t191 - t196) * MDP(17) + t198 * MDP(18)) * t231 + ((-t193 + t197) * MDP(14) - t198 * MDP(16) + (t190 - t305) * MDP(17) - t187 * MDP(18)) * t229 + ((-t222 * t268 - t223 * t269) * MDP(14) + (-t175 * t269 + t176 * t268 - t228 * t326) * MDP(15)) * pkin(3) + t335 * (-t246 * t265 - t288); (t193 * t231 + t194 * t229 + t235) * MDP(15) + t247 * MDP(18) + (t191 * t229 + (-qJD(5) - t190) * t231 + t282) * MDP(19) + (0.2e1 * t339 + (t229 - t299) * MDP(18)) * qJD(3) + (MDP(14) + MDP(17)) * (-t229 ^ 2 - t226); t229 * t339 + (-t247 + (t229 + t299) * qJD(3)) * MDP(17) + (-t226 - t274) * MDP(18) + (-qJD(3) * t191 + t187 * t231 + t175) * MDP(19);];
tauc = t1;
