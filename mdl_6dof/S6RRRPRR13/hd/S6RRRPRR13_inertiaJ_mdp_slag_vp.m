% Calculate joint inertia matrix for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR13_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR13_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR13_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:05:16
% EndTime: 2019-03-09 20:05:24
% DurationCPUTime: 2.75s
% Computational Cost: add. (4745->373), mult. (12275->549), div. (0->0), fcn. (14190->14), ass. (0->165)
t251 = cos(pkin(7));
t248 = sin(pkin(7));
t258 = cos(qJ(3));
t316 = t248 * t258;
t284 = pkin(10) * t316;
t255 = sin(qJ(3));
t325 = pkin(2) * t255;
t217 = t284 + (qJ(4) + t325) * t251;
t218 = (-pkin(3) * t258 - qJ(4) * t255 - pkin(2)) * t248;
t247 = sin(pkin(13));
t250 = cos(pkin(13));
t197 = -t217 * t247 + t250 * t218;
t317 = t248 * t255;
t223 = t247 * t251 + t250 * t317;
t188 = -pkin(4) * t316 - pkin(11) * t223 + t197;
t198 = t250 * t217 + t247 * t218;
t234 = t247 * t317;
t273 = t250 * t251 - t234;
t191 = t273 * pkin(11) + t198;
t254 = sin(qJ(5));
t326 = cos(qJ(5));
t176 = t254 * t188 + t326 * t191;
t174 = -pkin(12) * t316 + t176;
t199 = t223 * t254 - t326 * t273;
t200 = t326 * t223 + t254 * t273;
t236 = pkin(10) * t317;
t240 = -pkin(4) * t250 - pkin(3);
t324 = pkin(2) * t258;
t201 = t234 * pkin(4) + t236 + (t240 - t324) * t251;
t181 = t199 * pkin(5) - t200 * pkin(12) + t201;
t253 = sin(qJ(6));
t257 = cos(qJ(6));
t166 = -t174 * t253 + t181 * t257;
t167 = t174 * t257 + t181 * t253;
t298 = t199 * MDP(33);
t196 = t200 * t257 - t253 * t316;
t299 = t196 * MDP(31);
t195 = t200 * t253 + t257 * t316;
t301 = t195 * MDP(32);
t346 = t166 * MDP(34) - t167 * MDP(35) + t298 + t299 - t301;
t252 = cos(pkin(6));
t327 = cos(qJ(2));
t283 = pkin(1) * t327;
t237 = t252 * t283;
t249 = sin(pkin(6));
t256 = sin(qJ(2));
t315 = t249 * t256;
t209 = pkin(2) * t252 + t237 + (-pkin(10) * t251 - pkin(9)) * t315;
t216 = (-pkin(10) * t248 * t256 - t327 * pkin(2) - pkin(1)) * t249;
t192 = -t209 * t248 + t251 * t216;
t280 = t251 * t327;
t272 = t249 * t280;
t207 = -t252 * t316 + t255 * t315 - t258 * t272;
t208 = t252 * t317 + (t255 * t280 + t256 * t258) * t249;
t179 = pkin(3) * t207 - qJ(4) * t208 + t192;
t281 = t249 * t327;
t227 = t252 * t256 * pkin(1) + pkin(9) * t281;
t206 = (t248 * t252 + t272) * pkin(10) + t227;
t269 = t209 * t251 + t216 * t248;
t187 = t206 * t258 + t269 * t255;
t221 = t248 * t281 - t252 * t251;
t182 = -qJ(4) * t221 + t187;
t168 = t250 * t179 - t182 * t247;
t194 = t208 * t250 - t221 * t247;
t164 = pkin(4) * t207 - pkin(11) * t194 + t168;
t169 = t247 * t179 + t250 * t182;
t274 = t208 * t247 + t221 * t250;
t165 = -t274 * pkin(11) + t169;
t162 = t254 * t164 + t326 * t165;
t160 = t207 * pkin(12) + t162;
t186 = -t255 * t206 + t269 * t258;
t183 = t221 * pkin(3) - t186;
t170 = t274 * pkin(4) + t183;
t184 = t194 * t254 + t326 * t274;
t185 = t326 * t194 - t254 * t274;
t163 = t184 * pkin(5) - t185 * pkin(12) + t170;
t157 = -t160 * t253 + t163 * t257;
t158 = t160 * t257 + t163 * t253;
t303 = t184 * MDP(33);
t172 = t185 * t257 + t207 * t253;
t304 = t172 * MDP(31);
t171 = t185 * t253 - t207 * t257;
t306 = t171 * MDP(32);
t345 = t157 * MDP(34) - t158 * MDP(35) + t303 + t304 - t306;
t323 = pkin(11) + qJ(4);
t233 = t323 * t250;
t275 = t323 * t247;
t211 = t326 * t233 - t254 * t275;
t344 = t211 * MDP(28);
t230 = t247 * t254 - t326 * t250;
t231 = t326 * t247 + t254 * t250;
t203 = pkin(5) * t230 - pkin(12) * t231 + t240;
t189 = t203 * t257 - t211 * t253;
t190 = t203 * t253 + t211 * t257;
t343 = MDP(34) * t189 - MDP(35) * t190;
t210 = t233 * t254 + t326 * t275;
t340 = -t231 * MDP(24) + t210 * MDP(27) + MDP(14) + t344;
t319 = t231 * t257;
t320 = t231 * t253;
t339 = -t231 * MDP(23) + t240 * MDP(27) + t319 * MDP(31) - t320 * MDP(32) + t343;
t338 = 2 * MDP(16);
t337 = 2 * MDP(17);
t336 = 2 * MDP(18);
t335 = 2 * MDP(19);
t334 = 2 * MDP(20);
t333 = -0.2e1 * MDP(23);
t332 = 0.2e1 * MDP(27);
t331 = 0.2e1 * MDP(28);
t330 = -2 * MDP(30);
t329 = 0.2e1 * MDP(34);
t328 = 0.2e1 * MDP(35);
t322 = pkin(3) * MDP(21);
t321 = qJ(4) * t207;
t243 = t249 ^ 2;
t318 = t243 * t256;
t314 = t252 * MDP(8);
t312 = MDP(18) * t250;
t311 = MDP(19) * t247;
t310 = MDP(28) * t231;
t309 = MDP(29) * t257;
t308 = MDP(30) * t231;
t307 = MDP(33) * t230;
t305 = t172 * MDP(29);
t302 = t185 * MDP(24);
t300 = t196 * MDP(29);
t297 = t200 * MDP(22);
t296 = t200 * MDP(24);
t295 = t207 * MDP(25);
t294 = t207 * MDP(26);
t293 = t208 * MDP(11);
t292 = t208 * MDP(12);
t291 = t208 * MDP(13);
t290 = t221 * MDP(13);
t289 = t221 * MDP(14);
t287 = t251 * MDP(14);
t286 = t251 * MDP(15);
t285 = t258 * MDP(26);
t282 = qJ(4) * t316;
t279 = MDP(25) * t316;
t278 = t231 * t309;
t277 = MDP(30) * t253 * t257;
t276 = MDP(13) * t317;
t271 = -t168 * t247 + t169 * t250;
t270 = -t197 * t247 + t198 * t250;
t268 = MDP(31) * t257 - MDP(32) * t253;
t267 = MDP(34) * t257 - MDP(35) * t253;
t266 = -MDP(34) * t253 - MDP(35) * t257;
t161 = t326 * t164 - t254 * t165;
t175 = t326 * t188 - t254 * t191;
t265 = MDP(27) + t267;
t264 = (t256 * MDP(6) + t327 * MDP(7)) * t249;
t263 = t175 * MDP(27) - t176 * MDP(28) + t296;
t262 = t161 * MDP(27) - t162 * MDP(28) + t294 + t302;
t261 = t253 * MDP(31) + t257 * MDP(32) + t266 * pkin(12);
t260 = -MDP(25) + t261;
t246 = t257 ^ 2;
t245 = t253 ^ 2;
t242 = t248 ^ 2;
t226 = t251 * t325 + t284;
t225 = -pkin(9) * t315 + t237;
t224 = t251 * t324 - t236;
t220 = t236 + (-pkin(3) - t324) * t251;
t173 = pkin(5) * t316 - t175;
t159 = -t207 * pkin(5) - t161;
t1 = [t221 ^ 2 * MDP(15) + t185 ^ 2 * MDP(22) + (t168 ^ 2 + t169 ^ 2 + t183 ^ 2) * MDP(21) + MDP(1) + (MDP(4) * t256 + 0.2e1 * MDP(5) * t327) * t318 + (-0.2e1 * t290 + t293) * t208 + (t171 * t330 + t305) * t172 + (0.2e1 * t264 + t314) * t252 + (0.2e1 * t289 - 0.2e1 * t292 + t294 + 0.2e1 * t302) * t207 + (t185 * t333 - 0.2e1 * t295 + t303 + 0.2e1 * t304 - 0.2e1 * t306) * t184 + (t157 * t184 + t159 * t171) * t329 + (-t158 * t184 + t159 * t172) * t328 + (t161 * t207 + t170 * t184) * t332 + (-t162 * t207 + t170 * t185) * t331 + (-t169 * t207 + t183 * t194) * t335 + (t187 * t221 + t192 * t208) * t337 + (-t186 * t221 + t192 * t207) * t338 + (-t168 * t194 - t169 * t274) * t334 + (t168 * t207 + t183 * t274) * t336 + 0.2e1 * (-pkin(1) * t318 - t227 * t252) * MDP(10) + 0.2e1 * (t225 * t252 + t243 * t283) * MDP(9); t225 * MDP(9) + t314 - t227 * MDP(10) + (-t183 * t273 + t220 * t274) * MDP(18) + (-t168 * t223 + t169 * t273 - t197 * t194 - t198 * t274) * MDP(20) + (t186 * t251 - t221 * t224) * MDP(16) + (-t187 * t251 + t221 * t226) * MDP(17) + (t172 * t199 + t184 * t196) * MDP(31) + (t157 * t199 + t159 * t195 + t166 * t184 + t171 * t173) * MDP(34) + (-t171 * t199 - t184 * t195) * MDP(32) + (-t158 * t199 + t159 * t196 - t167 * t184 + t172 * t173) * MDP(35) + (-t184 * t200 - t185 * t199) * MDP(23) + (-t171 * t196 - t172 * t195) * MDP(30) + (t168 * t197 + t169 * t198 + t183 * t220) * MDP(21) + (t170 * t200 + t185 * t201) * MDP(28) + (t170 * t199 + t184 * t201) * MDP(27) + (t183 * t223 + t194 * t220) * MDP(19) + t251 * t291 + t185 * t297 + t184 * t298 + t172 * t300 + ((-t207 * MDP(16) - t208 * MDP(17)) * pkin(2) + (-t207 * MDP(12) + t192 * MDP(17) - t290 + t293) * t255 + (-t192 * MDP(16) - t168 * MDP(18) + t169 * MDP(19) + t184 * MDP(25) - t262 - t289 + t292) * t258) * t248 - t221 * t286 + (t197 * MDP(18) - t198 * MDP(19) - t199 * MDP(25) + t263 - t287) * t207 + t264; t242 * t255 ^ 2 * MDP(11) + t200 ^ 2 * MDP(22) + (t197 ^ 2 + t198 ^ 2 + t220 ^ 2) * MDP(21) + MDP(8) + (0.2e1 * t276 + t286) * t251 + (t195 * t330 + t300) * t196 + (0.2e1 * (t287 - t296) * t248 + (0.2e1 * t255 * MDP(12) + t285) * t242) * t258 + (t200 * t333 + 0.2e1 * t279 + t298 + 0.2e1 * t299 - 0.2e1 * t301) * t199 + (t166 * t199 + t173 * t195) * t329 + (-t167 * t199 + t173 * t196) * t328 + (-t197 * t223 + t198 * t273) * t334 + (-t226 * t251 - t242 * t325) * t337 + (-t175 * t316 + t199 * t201) * t332 + (t176 * t316 + t200 * t201) * t331 + (t198 * t316 + t220 * t223) * t335 + (t224 * t251 + t242 * t324) * t338 + (-t197 * t316 - t220 * t273) * t336; t291 - t221 * MDP(15) + t186 * MDP(16) - t187 * MDP(17) + (-pkin(3) * t274 - t183 * t250 - t247 * t321) * MDP(18) + (-pkin(3) * t194 + t183 * t247 - t250 * t321) * MDP(19) + ((t247 * t194 - t250 * t274) * qJ(4) + t271) * MDP(20) + (-pkin(3) * t183 + t271 * qJ(4)) * MDP(21) + t185 * t231 * MDP(22) + (t170 * t231 + t185 * t240) * MDP(28) + t172 * t278 + (-t171 * t257 - t172 * t253) * t308 + (t159 * t320 + t171 * t210) * MDP(34) + (t159 * t319 + t172 * t210) * MDP(35) + t339 * t184 + (-t185 * MDP(23) + t170 * MDP(27) - t295 + t345) * t230 - t340 * t207; t276 + t286 + t224 * MDP(16) - t226 * MDP(17) + (pkin(3) * t273 - t220 * t250 + t247 * t282) * MDP(18) + (-pkin(3) * t223 + t220 * t247 + t250 * t282) * MDP(19) + ((t247 * t223 + t250 * t273) * qJ(4) + t270) * MDP(20) + (-pkin(3) * t220 + t270 * qJ(4)) * MDP(21) + t231 * t297 + (t200 * t240 + t201 * t231) * MDP(28) + t196 * t278 + (-t195 * t257 - t196 * t253) * t308 + (t173 * t320 + t195 * t210) * MDP(34) + (t173 * t319 + t196 * t210) * MDP(35) + t340 * t316 + t339 * t199 + (-t200 * MDP(23) + t201 * MDP(27) + t279 + t346) * t230; 0.2e1 * t240 * t310 + MDP(15) + (-0.2e1 * t311 + 0.2e1 * t312 + t322) * pkin(3) + (MDP(29) * t246 + MDP(22) - 0.2e1 * t277) * t231 ^ 2 + (t240 * t332 + t307 + 0.2e1 * (-MDP(23) + t268) * t231) * t230 + (t189 * t230 + t210 * t320) * t329 + (-t190 * t230 + t210 * t319) * t328 + (MDP(21) * qJ(4) + t334) * (t247 ^ 2 + t250 ^ 2) * qJ(4); t274 * MDP(18) + t194 * MDP(19) + t183 * MDP(21) + t185 * MDP(28) + t265 * t184; -t273 * MDP(18) + t223 * MDP(19) + t220 * MDP(21) + t200 * MDP(28) + t265 * t199; t265 * t230 + t310 + t311 - t312 - t322; MDP(21); t253 * t305 + (-t171 * t253 + t172 * t257) * MDP(30) + (-pkin(5) * t171 - t159 * t257) * MDP(34) + (-pkin(5) * t172 + t159 * t253) * MDP(35) + t260 * t184 + t262; -t248 * t285 + t253 * t300 + (-t195 * t253 + t196 * t257) * MDP(30) + (-pkin(5) * t195 - t173 * t257) * MDP(34) + (-pkin(5) * t196 + t173 * t253) * MDP(35) + t260 * t199 + t263; -t344 - t265 * t210 + t260 * t230 + (MDP(24) + t253 * t309 + (-t245 + t246) * MDP(30) + t266 * pkin(5)) * t231; 0; MDP(29) * t245 + 0.2e1 * pkin(5) * t267 + MDP(26) + 0.2e1 * t277; t345; t346; t268 * t231 + t307 + t343; t267; t261; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
