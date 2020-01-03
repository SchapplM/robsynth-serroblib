% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:16
% EndTime: 2019-12-31 18:22:21
% DurationCPUTime: 1.92s
% Computational Cost: add. (1287->264), mult. (3130->395), div. (0->0), fcn. (2067->8), ass. (0->123)
t271 = sin(qJ(3));
t331 = MDP(5) * t271;
t264 = t271 ^ 2;
t273 = cos(qJ(3));
t330 = (-t273 ^ 2 + t264) * MDP(6);
t258 = sin(pkin(8)) * pkin(1) + pkin(6);
t254 = t258 * qJD(1);
t329 = qJD(2) * t273 - t271 * t254;
t328 = pkin(7) + qJ(4);
t266 = sin(pkin(9));
t268 = cos(pkin(9));
t270 = sin(qJ(5));
t272 = cos(qJ(5));
t245 = t266 * t272 + t268 * t270;
t326 = t266 * t270;
t244 = -t272 * t268 + t326;
t280 = t244 * t273;
t307 = qJD(5) * t271;
t188 = -qJD(3) * t280 - t245 * t307;
t313 = qJD(1) * t273;
t256 = -qJD(5) + t313;
t327 = t188 * t256;
t325 = t266 * t271;
t324 = t266 * t273;
t323 = t268 * t271;
t322 = t268 * t273;
t274 = qJD(3) ^ 2;
t321 = t271 * t274;
t320 = t273 * t274;
t281 = t245 * t273;
t277 = qJD(3) * t281;
t306 = qJD(5) * t272;
t189 = t306 * t323 - t307 * t326 + t277;
t223 = t245 * t271;
t304 = qJD(1) * qJD(3);
t300 = t271 * t304;
t319 = t189 * t256 - t223 * t300;
t210 = (qJD(4) + t329) * qJD(3);
t292 = pkin(3) * t271 - qJ(4) * t273;
t231 = qJD(3) * t292 - qJD(4) * t271;
t218 = t231 * qJD(1);
t179 = t268 * t210 + t266 * t218;
t263 = t271 * qJD(2);
t228 = t273 * t254 + t263;
t214 = qJD(3) * qJ(4) + t228;
t259 = -cos(pkin(8)) * pkin(1) - pkin(2);
t237 = -pkin(3) * t273 - qJ(4) * t271 + t259;
t217 = t237 * qJD(1);
t182 = t268 * t214 + t266 * t217;
t249 = t292 * qJD(1);
t193 = t266 * t249 + t268 * t329;
t318 = qJD(1) * t280 - t244 * qJD(5);
t317 = -qJD(1) * t281 + t245 * qJD(5);
t309 = qJD(3) * t273;
t219 = qJD(3) * t263 + t254 * t309;
t310 = qJD(3) * t271;
t301 = t258 * t310;
t195 = t268 * t231 + t266 * t301;
t248 = t258 * t322;
t201 = t266 * t237 + t248;
t305 = t268 * qJD(3);
t314 = qJD(1) * t271;
t240 = t266 * t314 - t305;
t299 = t273 * t304;
t294 = t268 * t299;
t316 = -t240 * t306 + t272 * t294;
t255 = qJD(1) * t259;
t311 = qJD(3) * t266;
t174 = -pkin(7) * t240 + t182;
t308 = qJD(5) * t174;
t303 = t258 * t324;
t302 = t266 * t313;
t298 = MDP(20) * t310;
t297 = pkin(4) * t266 + t258;
t296 = -qJD(3) * pkin(3) + qJD(4);
t181 = -t214 * t266 + t268 * t217;
t178 = -t210 * t266 + t268 * t218;
t192 = t268 * t249 - t266 * t329;
t242 = t268 * t314 + t311;
t172 = -pkin(4) * t313 - pkin(7) * t242 + t181;
t293 = t266 * t299;
t175 = -pkin(7) * t293 + t179;
t295 = -qJD(5) * t172 - t175;
t169 = t172 * t272 - t174 * t270;
t170 = t172 * t270 + t174 * t272;
t291 = -t178 * t266 + t179 * t268;
t290 = -t181 * t266 + t182 * t268;
t226 = t268 * t237;
t187 = -pkin(7) * t323 + t226 + (-t258 * t266 - pkin(4)) * t273;
t191 = -pkin(7) * t325 + t201;
t289 = t187 * t272 - t191 * t270;
t288 = t187 * t270 + t191 * t272;
t287 = t240 * t270 - t242 * t272;
t285 = 0.2e1 * qJD(3) * t255;
t284 = pkin(4) * t271 - pkin(7) * t322;
t253 = t328 * t268;
t283 = qJD(1) * t284 + qJD(4) * t266 + qJD(5) * t253 + t192;
t252 = t328 * t266;
t282 = pkin(7) * t302 + qJD(4) * t268 - qJD(5) * t252 - t193;
t279 = t284 * qJD(3);
t212 = -t329 + t296;
t278 = -qJD(5) * t242 - t293;
t276 = -qJ(4) * t310 + (-t212 + t296) * t273;
t177 = qJD(1) * t277 - qJD(5) * t287;
t260 = -pkin(4) * t268 - pkin(3);
t232 = t272 * t240;
t230 = t297 * t271;
t224 = t244 * t271;
t222 = t297 * t309;
t220 = t266 * t231;
t204 = pkin(4) * t302 + t228;
t202 = pkin(4) * t293 + t219;
t200 = t226 - t303;
t197 = t242 * t270 + t232;
t196 = -t268 * t301 + t220;
t194 = pkin(4) * t240 + t212;
t190 = t287 * t310;
t186 = t220 + (-pkin(7) * t324 - t258 * t323) * qJD(3);
t184 = t279 + t195;
t176 = t270 * t278 + t316;
t173 = qJD(1) * t279 + t178;
t171 = t272 * t173;
t1 = [0.2e1 * t299 * t331 - 0.2e1 * t304 * t330 + MDP(7) * t320 - MDP(8) * t321 + (-t258 * t320 + t271 * t285) * MDP(10) + (t258 * t321 + t273 * t285) * MDP(11) + (t219 * t325 + (-qJD(1) * t195 - t178) * t273 + ((t212 * t266 + t240 * t258) * t273 + (t181 + (t200 + t303) * qJD(1)) * t271) * qJD(3)) * MDP(12) + (t219 * t323 + (qJD(1) * t196 + t179) * t273 + ((t212 * t268 + t242 * t258) * t273 + (-t182 + (-t201 + t248) * qJD(1)) * t271) * qJD(3)) * MDP(13) + (-t195 * t242 - t196 * t240 + (-t178 * t268 - t179 * t266) * t271 + (-t181 * t268 - t182 * t266 + (-t200 * t268 - t201 * t266) * qJD(1)) * t309) * MDP(14) + (t178 * t200 + t179 * t201 + t181 * t195 + t182 * t196 + (t212 * t309 + t219 * t271) * t258) * MDP(15) + (-t176 * t224 - t188 * t287) * MDP(16) + (-t176 * t223 + t177 * t224 - t188 * t197 + t189 * t287) * MDP(17) + (-t176 * t273 - t224 * t300 - t190 - t327) * MDP(18) + (t177 * t273 - t197 * t310 + t319) * MDP(19) + (-t256 - t313) * t298 + (-(t184 * t272 - t186 * t270) * t256 - (-t175 * t270 + t171) * t273 + t222 * t197 + t230 * t177 + t202 * t223 + t194 * t189 + (t170 * t273 + t256 * t288) * qJD(5) + (qJD(1) * t289 + t169) * t310) * MDP(21) + ((t184 * t270 + t186 * t272) * t256 + (t173 * t270 + t175 * t272) * t273 - t222 * t287 + t230 * t176 - t202 * t224 + t194 * t188 + (t169 * t273 + t256 * t289) * qJD(5) + (-qJD(1) * t288 - t170) * t310) * MDP(22); t319 * MDP(21) + (-t190 + t327) * MDP(22) + (-t274 * MDP(11) - t219 * MDP(15) - t177 * MDP(21) - t176 * MDP(22)) * t273 + (-t274 * MDP(10) + MDP(15) * t291) * t271 + ((-t266 * MDP(12) - t268 * MDP(13)) * t264 * qJD(1) + (MDP(22) * qJD(1) * t224 + t240 * MDP(12) + t242 * MDP(13) + t212 * MDP(15) + t197 * MDP(21)) * t271 + ((-t240 * t268 + t242 * t266) * MDP(14) + t290 * MDP(15)) * t273) * qJD(3); (qJD(3) * t228 - t255 * t314 - t219) * MDP(10) - t255 * t313 * MDP(11) + (-t219 * t268 - t228 * t240 + (-t181 * t271 + t192 * t273 + t266 * t276) * qJD(1)) * MDP(12) + (t219 * t266 - t228 * t242 + (t182 * t271 - t193 * t273 + t268 * t276) * qJD(1)) * MDP(13) + (t192 * t242 + t193 * t240 + (-qJD(4) * t240 + t181 * t313 + t179) * t268 + (qJD(4) * t242 + t182 * t313 - t178) * t266) * MDP(14) + (-pkin(3) * t219 + t291 * qJ(4) + t290 * qJD(4) - t181 * t192 - t182 * t193 - t212 * t228) * MDP(15) + (t176 * t245 - t287 * t318) * MDP(16) + (-t176 * t244 - t177 * t245 - t318 * t197 + t287 * t317) * MDP(17) + (-t318 * t256 + (qJD(3) * t245 + t287) * t314) * MDP(18) + (t317 * t256 + (-qJD(3) * t244 + t197) * t314) * MDP(19) + t256 * MDP(20) * t314 + (t260 * t177 - t204 * t197 + t202 * t244 + (t270 * t282 + t272 * t283) * t256 + t317 * t194 + ((-t252 * t272 - t253 * t270) * qJD(3) - t169) * t314) * MDP(21) + (t260 * t176 + t204 * t287 + t202 * t245 + (-t270 * t283 + t272 * t282) * t256 + t318 * t194 + (-(-t252 * t270 + t253 * t272) * qJD(3) + t170) * t314) * MDP(22) + (-t273 * t331 + t330) * qJD(1) ^ 2; (-t240 ^ 2 - t242 ^ 2) * MDP(14) + (t181 * t242 + t182 * t240 + t219) * MDP(15) + (t256 * t287 + t177) * MDP(21) + (t232 * t256 + (-t293 + (-qJD(5) + t256) * t242) * t270 + t316) * MDP(22) + ((-t242 + t311) * MDP(12) + (t240 + t305) * MDP(13)) * t313; -t197 ^ 2 * MDP(17) + (-t197 * t256 + t316) * MDP(18) + qJD(1) * t298 + (-t170 * t256 + t171) * MDP(21) + (-t169 * t256 + t194 * t197) * MDP(22) - (MDP(16) * t197 - MDP(17) * t287 - MDP(19) * t256 - t194 * MDP(21)) * t287 + (t278 * MDP(19) - MDP(21) * t308 + t295 * MDP(22)) * t272 + (t278 * MDP(18) + (qJD(5) * t240 - t294) * MDP(19) + t295 * MDP(21) + (-t173 + t308) * MDP(22)) * t270;];
tauc = t1;
