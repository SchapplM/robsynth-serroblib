% Calculate joint inertia matrix for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP8_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:56:00
% EndTime: 2019-03-10 01:56:06
% DurationCPUTime: 2.09s
% Computational Cost: add. (2944->328), mult. (6370->447), div. (0->0), fcn. (7131->10), ass. (0->140)
t257 = sin(qJ(4));
t241 = pkin(3) * t257 + pkin(11);
t256 = sin(qJ(5));
t251 = t256 ^ 2;
t260 = cos(qJ(5));
t252 = t260 ^ 2;
t314 = t251 + t252;
t316 = t314 * t241;
t351 = t256 * MDP(27) + t260 * MDP(28);
t258 = sin(qJ(3));
t337 = cos(qJ(4));
t338 = cos(qJ(3));
t226 = t257 * t338 + t258 * t337;
t350 = 0.2e1 * t226;
t349 = pkin(10) + pkin(9);
t255 = cos(pkin(6));
t254 = sin(pkin(6));
t261 = cos(qJ(2));
t322 = t254 * t261;
t296 = pkin(8) * t322;
t259 = sin(qJ(2));
t336 = pkin(1) * t259;
t213 = t296 + (pkin(9) + t336) * t255;
t214 = (-pkin(2) * t261 - pkin(9) * t259 - pkin(1)) * t254;
t191 = -t213 * t258 + t338 * t214;
t323 = t254 * t259;
t218 = t255 * t258 + t323 * t338;
t297 = pkin(3) * t322;
t182 = -pkin(10) * t218 + t191 - t297;
t192 = t213 * t338 + t258 * t214;
t234 = t258 * t323;
t283 = -t255 * t338 + t234;
t186 = -pkin(10) * t283 + t192;
t291 = t337 * t186;
t174 = t257 * t182 + t291;
t172 = -pkin(11) * t322 + t174;
t197 = t218 * t257 + t283 * t337;
t198 = t218 * t337 - t257 * t283;
t237 = pkin(8) * t323;
t243 = -pkin(3) * t338 - pkin(2);
t332 = t261 * pkin(1);
t201 = t234 * pkin(3) + t237 + (t243 - t332) * t255;
t177 = t197 * pkin(4) - t198 * pkin(11) + t201;
t166 = t260 * t172 + t256 * t177;
t330 = qJ(6) * t197;
t163 = t166 + t330;
t287 = t172 * t256 - t260 * t177;
t335 = pkin(5) * t197;
t164 = t287 - t335;
t348 = t163 * t260 + t164 * t256;
t347 = t198 * MDP(20) - t197 * MDP(21);
t230 = t349 * t258;
t231 = t349 * t338;
t208 = t337 * t230 + t231 * t257;
t209 = -t257 * t230 + t231 * t337;
t346 = -t208 * MDP(23) - t209 * MDP(24);
t345 = MDP(31) * t256;
t344 = -(t258 * MDP(16) + MDP(17) * t338) * pkin(9) + t258 * MDP(13) + MDP(14) * t338;
t343 = 0.2e1 * MDP(23);
t342 = 0.2e1 * MDP(24);
t341 = 2 * MDP(32);
t340 = 2 * MDP(33);
t339 = 2 * MDP(34);
t225 = t257 * t258 - t337 * t338;
t334 = pkin(5) * t225;
t333 = pkin(11) * t225;
t295 = t337 * pkin(3);
t242 = -t295 - pkin(4);
t331 = pkin(4) - t242;
t329 = qJ(6) * t225;
t318 = -t337 * t182 + t257 * t186;
t171 = pkin(4) * t322 + t318;
t189 = t198 * t256 + t260 * t322;
t190 = t198 * t260 - t256 * t322;
t167 = pkin(5) * t189 - qJ(6) * t190 + t171;
t328 = t167 * t256;
t327 = t167 * t260;
t326 = t171 * t260;
t325 = t197 * t241;
t324 = t225 * t241;
t321 = t256 * t260;
t320 = t259 * MDP(6);
t203 = t225 * pkin(4) - t226 * pkin(11) + t243;
t185 = t256 * t203 + t260 * t209;
t282 = -t260 * pkin(5) - t256 * qJ(6);
t229 = -pkin(4) + t282;
t223 = -t295 + t229;
t317 = -t223 - t229;
t315 = t314 * pkin(11);
t313 = MDP(35) * t223;
t312 = MDP(35) * t229;
t311 = MDP(35) * t256;
t179 = t185 + t329;
t310 = t179 * MDP(35);
t286 = -t260 * t203 + t209 * t256;
t180 = t286 - t334;
t309 = t180 * MDP(35);
t281 = pkin(5) * t256 - qJ(6) * t260;
t187 = t226 * t281 + t208;
t308 = t187 * MDP(34);
t307 = t189 * MDP(26);
t306 = t189 * MDP(28);
t305 = t190 * MDP(25);
t304 = t190 * MDP(27);
t303 = t197 * MDP(29);
t302 = t198 * MDP(19);
t301 = t218 * MDP(11);
t300 = t225 * MDP(29);
t299 = 0.2e1 * t338;
t298 = -MDP(31) + MDP(34);
t294 = t256 * t325;
t293 = t260 * t325;
t292 = -t281 * MDP(33) + t351;
t290 = MDP(26) * t321;
t289 = t251 * MDP(25) + MDP(22) + 0.2e1 * t290;
t288 = -MDP(35) * pkin(5) - MDP(32);
t285 = MDP(35) * t314;
t284 = -pkin(4) * t226 - t333;
t280 = -t226 * t229 + t333;
t278 = -t189 * t260 + t190 * t256;
t277 = t223 * t226 - t324;
t276 = t226 * t242 - t324;
t275 = t260 * MDP(27) - t256 * MDP(28);
t274 = -MDP(30) * t286 - t185 * MDP(31);
t273 = MDP(30) * t260 - t345;
t272 = -t208 * MDP(30) - t187 * MDP(32);
t271 = -0.2e1 * MDP(32) * t260 - 0.2e1 * MDP(34) * t256;
t268 = (MDP(23) * t337 - t257 * MDP(24)) * pkin(3);
t267 = -t218 * MDP(13) + MDP(14) * t283;
t266 = -MDP(22) * t322 + (-t189 * t256 + t190 * t260) * MDP(26) + t256 * t305 + t347 + t351 * t197;
t265 = t208 * t345 + (t179 * t260 + t180 * t256) * MDP(33) + t346 + ((-t251 + t252) * MDP(26) + MDP(25) * t321 + MDP(20)) * t226 + (-MDP(21) + t351) * t225;
t264 = -MDP(30) * t287 - t166 * MDP(31) + t303 + t304 - t306;
t263 = (MDP(35) * qJ(6) + t298) * t260 + (-MDP(30) + t288) * t256;
t250 = t254 ^ 2;
t245 = t256 * MDP(33);
t220 = t255 * t336 + t296;
t219 = t255 * t332 - t237;
t212 = t237 + (-pkin(2) - t332) * t255;
t170 = t171 * t256;
t1 = [MDP(1) + (t163 ^ 2 + t164 ^ 2 + t167 ^ 2) * MDP(35) + t198 ^ 2 * MDP(18) + t255 ^ 2 * MDP(8) + (-0.2e1 * MDP(12) * t283 + t301) * t218 + (t305 - 0.2e1 * t307) * t190 + ((MDP(4) * t259 + 0.2e1 * MDP(5) * t261) * t259 + (MDP(15) + MDP(22)) * t261 ^ 2) * t250 + (-0.2e1 * t302 + t303 + 0.2e1 * t304 - 0.2e1 * t306) * t197 + 0.2e1 * (t255 * t320 + (t255 * MDP(7) + t267 - t347) * t261) * t254 + 0.2e1 * (t192 * t322 + t212 * t218) * MDP(17) + (t197 * t201 + t318 * t322) * t343 + (t174 * t322 + t198 * t201) * t342 + 0.2e1 * (t219 * t255 + t250 * t332) * MDP(9) + 0.2e1 * (-t220 * t255 - t250 * t336) * MDP(10) + (t163 * t197 - t167 * t190) * t339 + 0.2e1 * (t171 * t189 - t197 * t287) * MDP(30) + 0.2e1 * (-t166 * t197 + t171 * t190) * MDP(31) + (-t164 * t197 + t167 * t189) * t341 + (-t163 * t189 + t164 * t190) * t340 + 0.2e1 * (-t191 * t322 + t212 * t283) * MDP(16); (-t179 * t189 + t180 * t190) * MDP(33) + (-t180 * t197 + t187 * t189) * MDP(32) + (t179 * t197 - t187 * t190) * MDP(34) + (t189 * t208 - t197 * t286) * MDP(30) + (-t185 * t197 + t190 * t208) * MDP(31) + t258 * t301 + (t163 * t179 + t164 * t180 + t167 * t187) * MDP(35) + t219 * MDP(9) - t220 * MDP(10) + t255 * MDP(8) + (t218 * t338 - t258 * t283) * MDP(12) + (-pkin(2) * t218 + t212 * t258) * MDP(17) + (-pkin(2) * t283 - t212 * t338) * MDP(16) + (t197 * MDP(23) + t198 * MDP(24)) * t243 + (t320 + (MDP(7) - t344 - t346) * t261) * t254 + (MDP(21) * t322 + t201 * MDP(23) - t164 * MDP(32) + t163 * MDP(34) + t264 - t302) * t225 + (-MDP(20) * t322 + t198 * MDP(18) - t197 * MDP(19) + t201 * MDP(24) + (-t190 * MDP(26) - t197 * MDP(28) + t171 * MDP(30) + t167 * MDP(32) - t163 * MDP(33)) * t256 + (t197 * MDP(27) + t171 * MDP(31) + t164 * MDP(33) - t167 * MDP(34) + t305 - t307) * t260) * t226; MDP(8) + pkin(2) * MDP(16) * t299 + (t179 ^ 2 + t180 ^ 2 + t187 ^ 2) * MDP(35) + (MDP(11) * t258 + MDP(12) * t299 - 0.2e1 * pkin(2) * MDP(17)) * t258 + (t243 * t343 + t300 + (-MDP(19) + t275) * t350) * t225 + 0.2e1 * (-t180 * MDP(32) + t179 * MDP(34) + t274) * t225 + ((-t179 * t256 + t180 * t260) * MDP(33) + (t256 * MDP(30) + t260 * MDP(31)) * t208 + (t256 * MDP(32) - t260 * MDP(34)) * t187) * t350 + (t243 * t342 + (MDP(25) * t252 + MDP(18) - 0.2e1 * t290) * t226) * t226; (t167 * t223 + t241 * t348) * MDP(35) + t266 - MDP(15) * t322 + t191 * MDP(16) - t192 * MDP(17) + (t189 * t242 - t294 - t326) * MDP(30) + (t190 * t242 + t170 - t293) * MDP(31) + (t189 * t223 - t294 - t327) * MDP(32) + (-t190 * t223 + t293 - t328) * MDP(34) + (-t295 * t322 - t318) * MDP(23) + (-t291 + (-t182 + t297) * t257) * MDP(24) + (t241 * t278 + t348) * MDP(33) - t267; (MDP(31) * t276 - MDP(34) * t277 + t241 * t310 + t272) * t260 + t265 + (MDP(30) * t276 + MDP(32) * t277 + t241 * t309 - t308) * t256 + t187 * t313 + t344; MDP(15) + t316 * t340 + t241 ^ 2 * t285 + (t271 + t313) * t223 + t289 - 0.2e1 * t242 * t273 + 0.2e1 * t268; -t318 * MDP(23) - t174 * MDP(24) + (-pkin(4) * t189 - t326) * MDP(30) + (-pkin(4) * t190 + t170) * MDP(31) + (t189 * t229 - t327) * MDP(32) + t348 * MDP(33) + (-t190 * t229 - t328) * MDP(34) + t167 * t312 + (t278 * MDP(33) + t348 * MDP(35) + (t298 * t260 + (-MDP(30) - MDP(32)) * t256) * t197) * pkin(11) + t266; t187 * t312 + (MDP(31) * t284 + MDP(34) * t280 + pkin(11) * t310 + t272) * t260 + (MDP(30) * t284 - MDP(32) * t280 + pkin(11) * t309 - t308) * t256 + t265; (t315 + t316) * MDP(33) + (pkin(11) * t316 + t223 * t229) * MDP(35) + t268 + (MDP(30) * t331 + MDP(32) * t317) * t260 + (-MDP(31) * t331 + MDP(34) * t317) * t256 + t289; t315 * t340 + pkin(11) ^ 2 * t285 + (t271 + t312) * t229 + 0.2e1 * t273 * pkin(4) + t289; (-t287 + 0.2e1 * t335) * MDP(32) + (-pkin(5) * t190 - qJ(6) * t189) * MDP(33) + (t166 + 0.2e1 * t330) * MDP(34) + (-pkin(5) * t164 + qJ(6) * t163) * MDP(35) + t264; t300 + (-t286 + 0.2e1 * t334) * MDP(32) + (t185 + 0.2e1 * t329) * MDP(34) + (-pkin(5) * t180 + qJ(6) * t179) * MDP(35) + (MDP(33) * t282 + t275) * t226 + t274; t241 * t263 + t292; pkin(11) * t263 + t292; MDP(29) + pkin(5) * t341 + qJ(6) * t339 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(35); -t197 * MDP(32) + t190 * MDP(33) + t164 * MDP(35); MDP(33) * t226 * t260 - t225 * MDP(32) + t309; t241 * t311 + t245; pkin(11) * t311 + t245; t288; MDP(35);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
