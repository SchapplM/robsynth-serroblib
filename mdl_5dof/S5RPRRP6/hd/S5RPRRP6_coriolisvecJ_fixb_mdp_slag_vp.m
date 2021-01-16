% Calculate Coriolis joint torque vector for
% S5RPRRP6
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
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 18:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 18:08:33
% EndTime: 2021-01-15 18:08:39
% DurationCPUTime: 1.96s
% Computational Cost: add. (1514->282), mult. (3509->389), div. (0->0), fcn. (2061->6), ass. (0->124)
t248 = sin(qJ(3));
t333 = MDP(5) * t248;
t243 = t248 ^ 2;
t250 = cos(qJ(3));
t332 = (-t250 ^ 2 + t243) * MDP(6);
t236 = sin(pkin(8)) * pkin(1) + pkin(6);
t230 = t236 * qJD(1);
t241 = t248 * qJD(2);
t203 = t250 * t230 + t241;
t196 = qJD(3) * pkin(7) + t203;
t237 = -cos(pkin(8)) * pkin(1) - pkin(2);
t217 = -pkin(3) * t250 - pkin(7) * t248 + t237;
t199 = t217 * qJD(1);
t247 = sin(qJ(4));
t249 = cos(qJ(4));
t178 = -t196 * t247 + t249 * t199;
t289 = qJD(3) * t247;
t292 = qJD(1) * t248;
t225 = t249 * t292 + t289;
t175 = -qJ(5) * t225 + t178;
t291 = qJD(1) * t250;
t235 = -qJD(4) + t291;
t174 = -pkin(4) * t235 + t175;
t331 = t174 - t175;
t286 = qJD(4) * t247;
t275 = t248 * t286;
t280 = qJD(1) * qJD(3);
t272 = t250 * t280;
t279 = qJD(3) * qJD(4);
t294 = (t272 + t279) * t249;
t191 = qJD(1) * t275 - t294;
t278 = t247 * t292;
t281 = t249 * qJD(3);
t223 = t278 - t281;
t330 = t223 * t235 - t191;
t285 = qJD(4) * t249;
t274 = t248 * t285;
t287 = qJD(3) * t250;
t253 = t247 * t287 + t274;
t192 = t253 * qJD(1) + t247 * t279;
t329 = -t225 * t235 + t192;
t328 = qJD(2) * t250 - t248 * t230;
t327 = MDP(17) + MDP(19);
t260 = qJD(1) * t243 - t235 * t250;
t326 = -t235 * t275 - t260 * t281;
t325 = pkin(4) * t223;
t324 = pkin(4) * t248;
t323 = -qJ(5) - pkin(7);
t322 = qJ(5) * t248;
t307 = t247 * t199;
t179 = t196 * t249 + t307;
t176 = -qJ(5) * t223 + t179;
t321 = t176 * t235;
t198 = qJD(3) * t241 + t230 * t287;
t180 = pkin(4) * t192 + t198;
t320 = t180 * t247;
t319 = t180 * t249;
t318 = t192 * t250;
t195 = -qJD(3) * pkin(3) - t328;
t317 = t195 * t247;
t316 = t195 * t249;
t315 = t198 * t247;
t314 = t198 * t249;
t312 = t223 * t248;
t310 = t225 * t247;
t309 = t235 * t249;
t308 = t236 * t247;
t306 = t247 * t250;
t305 = t248 * t249;
t251 = qJD(3) ^ 2;
t304 = t248 * t251;
t303 = t249 * t250;
t302 = t250 * t251;
t257 = -qJ(5) * t303 + t324;
t263 = pkin(3) * t248 - pkin(7) * t250;
t228 = t263 * qJD(1);
t267 = t249 * t228 - t247 * t328;
t271 = qJD(4) * t323;
t301 = t257 * qJD(1) + qJD(5) * t247 - t249 * t271 + t267;
t273 = t247 * t291;
t284 = qJD(5) * t249;
t298 = t247 * t228 + t249 * t328;
t300 = -qJ(5) * t273 - t247 * t271 - t284 + t298;
t277 = t250 * t281;
t299 = -t192 * t305 - t223 * t277;
t229 = t263 * qJD(3);
t297 = t217 * t285 + t247 * t229;
t288 = qJD(3) * t248;
t296 = t249 * t229 + t288 * t308;
t227 = t236 * t303;
t295 = t247 * t217 + t227;
t231 = qJD(1) * t237;
t283 = t191 * MDP(21);
t282 = t195 * qJD(4);
t276 = t235 * t286;
t270 = -qJD(5) - t325;
t269 = t191 * t250 + t225 * t288;
t197 = t328 * qJD(3);
t216 = qJD(1) * t229;
t268 = t247 * t197 - t249 * t216;
t266 = t196 * t286 - t249 * t197 - t199 * t285 - t247 * t216;
t264 = t235 * t274;
t262 = -t174 * t249 - t176 * t247;
t261 = t174 * t247 - t176 * t249;
t258 = 0.2e1 * qJD(3) * t231;
t256 = qJ(5) * t191 - t268;
t255 = -qJ(5) * t192 - t266;
t254 = t260 * t247;
t240 = -pkin(4) * t249 - pkin(3);
t233 = t323 * t249;
t232 = t323 * t247;
t220 = t223 ^ 2;
t211 = (pkin(4) * t247 + t236) * t248;
t207 = t249 * t217;
t187 = pkin(4) * t273 + t203;
t186 = t253 * pkin(4) + t236 * t287;
t184 = t195 - t270;
t183 = -t247 * t322 + t295;
t182 = -qJ(5) * t305 + t207 + (-pkin(4) - t308) * t250;
t171 = (-qJ(5) * qJD(4) - qJD(3) * t236) * t305 + (-qJD(5) * t248 + (-qJ(5) * qJD(3) - qJD(4) * t236) * t250) * t247 + t297;
t170 = -t248 * t284 + t257 * qJD(3) + (-t227 + (-t217 + t322) * t247) * qJD(4) + t296;
t169 = -qJD(5) * t223 + t255;
t168 = -t179 * qJD(4) - qJD(5) * t225 + t280 * t324 + t256;
t1 = [0.2e1 * t272 * t333 - 0.2e1 * t280 * t332 + MDP(7) * t302 - MDP(8) * t304 + (-t236 * t302 + t248 * t258) * MDP(10) + (t236 * t304 + t250 * t258) * MDP(11) + (-t191 * t305 + (-t275 + t277) * t225) * MDP(12) + (-t225 * t274 + (-t225 * t287 + (qJD(4) * t223 + t191) * t248) * t247 + t299) * MDP(13) + (t269 - t326) * MDP(14) + (t264 + t318 + (-t254 - t312) * qJD(3)) * MDP(15) + (-t235 - t291) * MDP(16) * t288 + (-(-t217 * t286 + t296) * t235 + ((t223 * t236 + t317) * qJD(3) + (t307 + (t235 * t236 + t196) * t249) * qJD(4) + t268) * t250 + (t249 * t282 + t236 * t192 + t315 + ((-t236 * t306 + t207) * qJD(1) + t178) * qJD(3)) * t248) * MDP(17) + (t297 * t235 + (-t236 * t276 + (t225 * t236 + t316) * qJD(3) - t266) * t250 + (-t247 * t282 - t236 * t191 + t314 + (-t295 * qJD(1) - t236 * t309 - t179) * qJD(3)) * t248) * MDP(18) + (-t170 * t235 + t186 * t223 + t192 * t211 + (t184 * t289 - t168) * t250 + (t184 * t285 + t320 + (qJD(1) * t182 + t174) * qJD(3)) * t248) * MDP(19) + (t171 * t235 + t186 * t225 - t191 * t211 + (t184 * t281 + t169) * t250 + (-t184 * t286 + t319 + (-qJD(1) * t183 - t176) * qJD(3)) * t248) * MDP(20) + (-t170 * t225 - t171 * t223 + t182 * t191 - t183 * t192 + t262 * t287 + (t261 * qJD(4) - t168 * t249 - t169 * t247) * t248) * MDP(21) + (t168 * t182 + t169 * t183 + t170 * t174 + t171 * t176 + t180 * t211 + t184 * t186) * MDP(22); t299 * MDP(21) + t327 * (t264 - t318 + (-t254 + t312) * qJD(3)) + (MDP(18) + MDP(20)) * (t269 + t326) + (-t251 * MDP(11) - t180 * MDP(22) + (MDP(21) * t310 - t261 * MDP(22)) * qJD(3)) * t250 + (-t251 * MDP(10) - t247 * t283 + (qJD(3) * t184 - t168 * t247 + t169 * t249) * MDP(22) + ((t223 * t247 + t225 * t249) * MDP(21) + t262 * MDP(22)) * qJD(4)) * t248; (qJD(3) * t203 - t231 * t292 - t198) * MDP(10) - t231 * t291 * MDP(11) + (-t191 * t247 - t225 * t309) * MDP(12) + (-t329 * t247 + t330 * t249) * MDP(13) + (-t235 * t285 + (t235 * t303 + (-t225 + t289) * t248) * qJD(1)) * MDP(14) + (t276 + (-t235 * t306 + (t223 + t281) * t248) * qJD(1)) * MDP(15) + t235 * MDP(16) * t292 + (-pkin(3) * t192 - t314 + t267 * t235 - t203 * t223 + (pkin(7) * t309 + t317) * qJD(4) + (-t178 * t248 + (-pkin(7) * t288 - t195 * t250) * t247) * qJD(1)) * MDP(17) + (pkin(3) * t191 + t315 - t298 * t235 - t203 * t225 + (-pkin(7) * t235 * t247 + t316) * qJD(4) + (-t195 * t303 + (-pkin(7) * t281 + t179) * t248) * qJD(1)) * MDP(18) + (-t319 - t187 * t223 + t192 * t240 + t301 * t235 + (t184 + t325) * t286 + (-t184 * t306 + (qJD(3) * t232 - t174) * t248) * qJD(1)) * MDP(19) + (t320 - t187 * t225 - t191 * t240 - t300 * t235 + (pkin(4) * t310 + t184 * t249) * qJD(4) + (-t184 * t303 + (qJD(3) * t233 + t176) * t248) * qJD(1)) * MDP(20) + (t191 * t232 + t192 * t233 + t301 * t225 + t300 * t223 + (t235 * t174 + t169) * t249 + (-t168 + t321) * t247) * MDP(21) + (t168 * t232 - t169 * t233 + t180 * t240 + (pkin(4) * t286 - t187) * t184 - t300 * t176 - t301 * t174) * MDP(22) + (-t250 * t333 + t332) * qJD(1) ^ 2; -t220 * MDP(13) + t294 * MDP(14) + (-t179 * t235 - t268) * MDP(17) + (-t178 * t235 + t266) * MDP(18) + (t256 - t321) * MDP(19) + (-t175 * t235 - t255) * MDP(20) + pkin(4) * t283 + (pkin(4) * t168 + t331 * t176) * MDP(22) + (-t235 * MDP(14) + t195 * MDP(18) + (qJD(5) + t184) * MDP(20) - t331 * MDP(21)) * t223 + (-MDP(15) * t306 + (0.2e1 * MDP(19) * pkin(4) + MDP(16)) * t248) * t280 + (-MDP(14) * t278 - t225 * MDP(15) - t327 * t179) * qJD(4) + (t223 * MDP(12) - t235 * MDP(15) - t195 * MDP(17) + (-t184 + t270) * MDP(19) - pkin(4) * t184 * MDP(22) + (-MDP(20) * pkin(4) + MDP(13)) * t225) * t225; t329 * MDP(19) + t330 * MDP(20) + (-t225 ^ 2 - t220) * MDP(21) + (t174 * t225 + t176 * t223 + t180) * MDP(22);];
tauc = t1;
