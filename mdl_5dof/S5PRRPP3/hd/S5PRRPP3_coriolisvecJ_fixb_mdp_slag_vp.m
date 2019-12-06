% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:39
% EndTime: 2019-12-05 16:13:46
% DurationCPUTime: 2.03s
% Computational Cost: add. (1149->258), mult. (2821->351), div. (0->0), fcn. (1766->6), ass. (0->122)
t248 = sin(qJ(3));
t250 = cos(qJ(3));
t261 = pkin(3) * t248 - qJ(4) * t250;
t213 = t261 * qJD(3) - qJD(4) * t248;
t246 = sin(pkin(8));
t249 = sin(qJ(2));
t247 = cos(pkin(8));
t251 = cos(qJ(2));
t312 = t247 * t251;
t257 = t246 * t249 + t250 * t312;
t327 = -t257 * qJD(1) + t246 * t213;
t303 = qJD(1) * t251;
t304 = qJD(1) * t249;
t313 = t246 * t250;
t326 = t303 * t313 + (t213 - t304) * t247;
t289 = MDP(12) + MDP(16);
t325 = MDP(14) + MDP(17);
t245 = t250 ^ 2;
t324 = (t248 ^ 2 - t245) * MDP(6);
t294 = qJD(5) * t250;
t298 = qJD(3) * t248;
t319 = pkin(6) * t247;
t323 = -t294 + (qJ(5) - t319) * t298 + t327;
t284 = pkin(6) * t246 + pkin(4);
t309 = -t284 * t298 - t326;
t287 = pkin(6) * t298;
t322 = t246 * t287 + t326;
t307 = -t247 * t287 + t327;
t288 = MDP(13) - MDP(18);
t299 = qJD(3) * t246;
t302 = qJD(2) * t248;
t225 = t247 * t302 + t299;
t221 = t225 ^ 2;
t320 = -0.2e1 * t251;
t318 = qJD(2) * pkin(2);
t317 = qJD(3) * pkin(3);
t228 = t261 * qJD(2);
t315 = t228 * t247;
t232 = -pkin(3) * t250 - qJ(4) * t248 - pkin(2);
t314 = t232 * t247;
t238 = qJD(2) * pkin(6) + t304;
t229 = t248 * t238;
t311 = t249 * t250;
t230 = t250 * t238;
t195 = (t213 + t304) * qJD(2);
t300 = qJD(2) * t251;
t277 = qJD(1) * t300;
t198 = t250 * t277 + (qJD(4) - t229) * qJD(3);
t171 = t246 * t195 + t247 * t198;
t206 = t232 * qJD(2) - t303;
t218 = qJD(3) * qJ(4) + t230;
t177 = t246 * t206 + t247 * t218;
t237 = t248 * t277;
t297 = qJD(3) * t250;
t200 = t238 * t297 + t237;
t202 = t246 * t232 + t250 * t319;
t252 = qJD(3) ^ 2;
t253 = qJD(2) ^ 2;
t305 = t252 + t253;
t301 = qJD(2) * t250;
t296 = qJD(4) * t225;
t295 = qJD(5) * t225;
t274 = qJD(4) - t317;
t214 = t229 + t274;
t291 = t247 * qJD(3);
t223 = t246 * t302 - t291;
t175 = pkin(4) * t223 - qJ(5) * t225 + t214;
t293 = t175 * MDP(19);
t292 = t214 * MDP(15);
t290 = qJD(2) * qJD(3);
t286 = t247 * t311;
t283 = qJ(4) * t299;
t282 = qJ(4) * t291;
t281 = t223 * t303;
t280 = t225 * t303;
t279 = t246 * t301;
t278 = t249 * t298;
t276 = t248 * t290;
t275 = t250 * t290;
t170 = t195 * t247 - t246 * t198;
t239 = -t303 - t318;
t272 = -t239 - t303;
t176 = t206 * t247 - t218 * t246;
t173 = pkin(4) * t301 + qJD(5) - t176;
t197 = t284 * t250 - t314;
t271 = qJD(2) * t197 + t173;
t174 = -qJ(5) * t301 + t177;
t196 = -qJ(5) * t250 + t202;
t270 = qJD(2) * t196 + t174;
t260 = pkin(4) * t246 - qJ(5) * t247;
t258 = pkin(6) + t260;
t207 = t258 * t248;
t269 = qJD(2) * t207 + t175;
t201 = -pkin(6) * t313 + t314;
t268 = qJD(2) * t201 + t176;
t267 = -qJD(2) * t202 - t177;
t231 = -pkin(4) * t247 - qJ(5) * t246 - pkin(3);
t266 = qJD(3) * t231 - t175;
t265 = t247 * t275;
t236 = t246 * t275;
t264 = pkin(6) * t302 + t214;
t256 = -t272 - t318;
t255 = pkin(4) * t236 - qJ(5) * t265 + t200;
t254 = t288 * t225 + t292 + t293;
t235 = qJD(4) * t279;
t217 = t246 * t228;
t216 = -t246 * t251 + t286;
t215 = t246 * t311 + t312;
t212 = t223 * t301;
t211 = t247 * qJD(4) * t223;
t193 = t257 * qJD(2) - t247 * t278;
t192 = -qJD(2) * t249 * t247 - t246 * t278 + t251 * t279;
t187 = t260 * t301 + t230;
t185 = -t247 * t229 + t217;
t184 = t246 * t229 + t315;
t183 = -qJD(5) * t247 * t248 + t258 * t297;
t181 = -t315 + (-pkin(4) * qJD(2) - t238 * t246) * t248;
t180 = t217 + (qJ(5) * qJD(2) - t238 * t247) * t248;
t172 = t255 - t295;
t169 = -pkin(4) * t276 - t170;
t168 = (qJ(5) * t298 - t294) * qJD(2) + t171;
t1 = [-t253 * t251 * MDP(4) + (-t170 * t215 + t171 * t216 - t176 * t192 + t177 * t193) * MDP(15) + (t168 * t216 + t169 * t215 + t173 * t192 + t174 * t193) * MDP(19) + (-t253 * MDP(3) + (t305 * MDP(11) + t200 * MDP(15) + t172 * MDP(19)) * t248 + (-t305 * MDP(10) + t254 * qJD(3)) * t250) * t249 + ((t254 * t251 + (MDP(10) * t320 - t289 * t215 + t288 * (-t216 + t286)) * qJD(3)) * t248 + (t288 * t193 + t289 * t192 + (MDP(11) * t320 + t325 * (t215 * t247 - t216 * t246)) * qJD(3)) * t250) * qJD(2) + t289 * (t248 * t249 * t236 + (t248 * t300 + t249 * t297) * t223) + t325 * (t192 * t225 - t193 * t223); (t170 * t201 + t171 * t202 + t322 * t176 + t307 * t177) * MDP(15) + (t168 * t196 + t169 * t197 + t172 * t207 + t173 * t309 + t174 * t323 + t175 * t183) * MDP(19) - 0.2e1 * t290 * t324 + (-MDP(14) * t322 + t309 * MDP(17) - t183 * MDP(18)) * t225 + (-t307 * MDP(14) + t183 * MDP(16) - MDP(17) * t323) * t223 + ((t200 * t246 - t281) * MDP(12) + (t200 * t247 - t280) * MDP(13) + (-t170 * t247 - t171 * t246) * MDP(14) + (pkin(6) * t200 - t214 * t303) * MDP(15) + (t172 * t246 - t281) * MDP(16) + (-t168 * t246 + t169 * t247) * MDP(17) + (-t172 * t247 + t280) * MDP(18) - t293 * t303 + (pkin(6) * MDP(11) - MDP(8)) * t252 + (t256 * MDP(10) + t268 * MDP(12) + t267 * MDP(13) - t271 * MDP(16) + t270 * MDP(18)) * qJD(3)) * t248 + (-t170 * MDP(12) + t171 * MDP(13) + t169 * MDP(16) - t168 * MDP(18) + (-pkin(6) * MDP(10) + MDP(7)) * t252 + (-MDP(12) * t322 + t307 * MDP(13) + t309 * MDP(16) - MDP(18) * t323) * qJD(2) + (0.2e1 * MDP(5) * t302 + t256 * MDP(11) + (t223 * MDP(12) + t225 * MDP(13) + t292) * pkin(6) + (t264 * MDP(13) - t268 * MDP(14) + t271 * MDP(17) - t269 * MDP(18)) * t247 + (t264 * MDP(12) + t267 * MDP(14) + t269 * MDP(16) - t270 * MDP(17)) * t246) * qJD(3)) * t250; -t237 * MDP(10) + (-t223 * t230 + t235) * MDP(12) + (t185 * t223 - t211) * MDP(14) + (-pkin(3) * t200 - t176 * t184 - t177 * t185 - t214 * t230) * MDP(15) + (-t187 * t223 + t235) * MDP(16) + (t180 * t223 - t211) * MDP(17) + (t172 * t231 - t173 * t181 - t174 * t180 - t175 * t187) * MDP(19) + (-t200 * MDP(12) + t171 * MDP(14) + (qJ(4) * t171 + qJD(4) * t177) * MDP(15) - t172 * MDP(16) + t168 * MDP(17) + (qJ(4) * t168 + qJD(4) * t174) * MDP(19)) * t247 + (-MDP(13) * t230 + MDP(14) * t184 - MDP(17) * t181 + MDP(18) * t187) * t225 + (-t248 * t250 * MDP(5) + t324) * t253 + (t200 * MDP(13) + (-t170 + t296) * MDP(14) + (-qJ(4) * t170 - qJD(4) * t176) * MDP(15) - qJD(5) * t223 * MDP(16) + (t169 + t296) * MDP(17) + (-t172 + t295) * MDP(18) + (qJ(4) * t169 + qJD(4) * t173 - qJD(5) * t175) * MDP(19)) * t246 + ((-t239 * MDP(10) + (-t176 - t283) * MDP(12) + (t177 - t282) * MDP(13) + (t173 - t283) * MDP(16) + (-t174 + t282) * MDP(18)) * t248 + (t272 * MDP(11) + t184 * MDP(12) - t185 * MDP(13) - t181 * MDP(16) + t180 * MDP(18) + ((-t214 - t317) * MDP(12) + t177 * MDP(14) + t266 * MDP(16) + t174 * MDP(17)) * t246 + ((-t214 + t274) * MDP(13) + t176 * MDP(14) - t173 * MDP(17) + (-qJD(4) - t266) * MDP(18)) * t247) * t250) * qJD(2); (t176 * t225 + t177 * t223 + t200) * MDP(15) + (t174 * t223 + (-qJD(5) - t173) * t225 + t255) * MDP(19) + t288 * (t212 + t265) + t289 * (-t225 * t301 + t236) + t325 * (-t223 ^ 2 - t221); (t223 * t225 - t276) * MDP(16) + (-t212 + t265) * MDP(17) + (-t245 * t253 - t221) * MDP(18) + (t175 * t225 + (-pkin(4) * t298 + t174 * t250) * qJD(2) - t170) * MDP(19);];
tauc = t1;
