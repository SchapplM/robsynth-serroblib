% Calculate joint inertia matrix for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR9_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:10:23
% EndTime: 2019-03-09 19:10:29
% DurationCPUTime: 1.83s
% Computational Cost: add. (4212->332), mult. (11309->495), div. (0->0), fcn. (13043->14), ass. (0->155)
t246 = sin(qJ(6));
t250 = cos(qJ(6));
t264 = t246 * MDP(32) + t250 * MDP(33);
t259 = t246 * MDP(29) + t250 * MDP(30) - t264 * pkin(12);
t328 = 2 * MDP(16);
t327 = 2 * MDP(17);
t326 = 2 * MDP(18);
t325 = 2 * MDP(22);
t324 = -2 * MDP(23);
t323 = 2 * MDP(25);
t322 = 2 * MDP(26);
t321 = -2 * MDP(28);
t320 = 0.2e1 * MDP(32);
t319 = 0.2e1 * MDP(33);
t253 = cos(qJ(2));
t318 = pkin(1) * t253;
t317 = pkin(10) + qJ(4);
t245 = cos(pkin(6));
t249 = sin(qJ(2));
t242 = sin(pkin(6));
t304 = t242 * t253;
t223 = pkin(1) * t245 * t249 + pkin(9) * t304;
t241 = sin(pkin(7));
t244 = cos(pkin(7));
t302 = t244 * t253;
t275 = t242 * t302;
t198 = (t241 * t245 + t275) * pkin(10) + t223;
t230 = t245 * t318;
t305 = t242 * t249;
t204 = pkin(2) * t245 + t230 + (-pkin(10) * t244 - pkin(9)) * t305;
t213 = (-pkin(10) * t241 * t249 - pkin(2) * t253 - pkin(1)) * t242;
t248 = sin(qJ(3));
t252 = cos(qJ(3));
t303 = t244 * t252;
t306 = t241 * t252;
t180 = -t198 * t248 + t204 * t303 + t213 * t306;
t307 = t241 * t248;
t200 = t245 * t307 + (t248 * t302 + t249 * t252) * t242;
t218 = t241 * t304 - t245 * t244;
t171 = -pkin(3) * t218 - qJ(4) * t200 + t180;
t181 = t198 * t252 + (t204 * t244 + t213 * t241) * t248;
t199 = -t245 * t306 + t248 * t305 - t252 * t275;
t176 = -qJ(4) * t199 + t181;
t240 = sin(pkin(13));
t243 = cos(pkin(13));
t164 = t240 * t171 + t243 * t176;
t162 = -pkin(11) * t218 + t164;
t190 = -t204 * t241 + t244 * t213;
t184 = pkin(3) * t199 + t190;
t186 = t243 * t199 + t200 * t240;
t187 = -t199 * t240 + t200 * t243;
t167 = pkin(4) * t186 - pkin(11) * t187 + t184;
t247 = sin(qJ(5));
t251 = cos(qJ(5));
t158 = -t162 * t247 + t167 * t251;
t156 = -pkin(5) * t186 - t158;
t316 = t156 * t246;
t315 = t156 * t250;
t229 = pkin(2) * t303;
t208 = pkin(3) * t244 - t317 * t307 + t229;
t276 = pkin(2) * t244 * t248;
t212 = t317 * t306 + t276;
t192 = t240 * t208 + t243 * t212;
t189 = pkin(11) * t244 + t192;
t216 = t240 * t307 - t243 * t306;
t217 = (t240 * t252 + t243 * t248) * t241;
t225 = (-pkin(3) * t252 - pkin(2)) * t241;
t195 = pkin(4) * t216 - pkin(11) * t217 + t225;
t178 = -t189 * t247 + t195 * t251;
t172 = -pkin(5) * t216 - t178;
t314 = t172 * t246;
t313 = t172 * t250;
t233 = pkin(3) * t240 + pkin(11);
t312 = t233 * t246;
t311 = t233 * t250;
t310 = t233 * t251;
t235 = t241 ^ 2;
t309 = t235 * t248;
t236 = t242 ^ 2;
t308 = t236 * t249;
t301 = t245 * MDP(8);
t300 = t246 * t247;
t299 = t247 * t250;
t183 = t187 * t251 - t218 * t247;
t168 = t183 * t246 - t250 * t186;
t298 = t168 * MDP(30);
t169 = t183 * t250 + t186 * t246;
t297 = t169 * MDP(27);
t296 = t169 * MDP(29);
t182 = t187 * t247 + t218 * t251;
t295 = t182 * MDP(31);
t294 = t183 * MDP(20);
t293 = t183 * MDP(21);
t206 = t217 * t251 + t244 * t247;
t193 = t206 * t246 - t250 * t216;
t292 = t193 * MDP(30);
t194 = t206 * t250 + t216 * t246;
t291 = t194 * MDP(27);
t290 = t194 * MDP(29);
t289 = t199 * MDP(14);
t288 = t200 * MDP(11);
t205 = t217 * t247 - t251 * t244;
t287 = t205 * MDP(31);
t286 = t206 * MDP(20);
t285 = t206 * MDP(21);
t284 = t216 * MDP(24);
t283 = t218 * MDP(15);
t234 = -pkin(3) * t243 - pkin(4);
t282 = t234 * MDP(25);
t281 = t234 * MDP(26);
t280 = t244 * MDP(15);
t279 = t247 * MDP(26);
t278 = t250 * MDP(27);
t277 = t251 * MDP(31);
t274 = MDP(28) * t246 * t250;
t163 = t171 * t243 - t240 * t176;
t191 = t208 * t243 - t240 * t212;
t273 = -t233 * MDP(25) + MDP(22);
t272 = -t233 * MDP(26) + MDP(23);
t159 = t162 * t251 + t167 * t247;
t179 = t189 * t251 + t195 * t247;
t271 = t248 * MDP(13) + t252 * MDP(14);
t220 = -pkin(10) * t307 + t229;
t222 = pkin(10) * t306 + t276;
t270 = -t220 * MDP(16) + t222 * MDP(17);
t269 = t251 * MDP(25) - t279;
t268 = t250 * MDP(29) - t246 * MDP(30);
t224 = -pkin(5) * t251 - pkin(12) * t247 + t234;
t202 = t224 * t250 - t246 * t310;
t203 = t224 * t246 + t250 * t310;
t266 = t202 * MDP(32) - t203 * MDP(33);
t265 = MDP(32) * t250 - MDP(33) * t246;
t161 = pkin(4) * t218 - t163;
t188 = -pkin(4) * t244 - t191;
t263 = -MDP(21) + t268;
t262 = (t249 * MDP(6) + t253 * MDP(7)) * t242;
t261 = t271 * t241;
t260 = t266 + t282;
t258 = t200 * MDP(13) + t180 * MDP(16) - t181 * MDP(17) - t283 - t289;
t157 = pkin(12) * t186 + t159;
t160 = pkin(5) * t182 - pkin(12) * t183 + t161;
t154 = -t157 * t246 + t160 * t250;
t155 = t157 * t250 + t160 * t246;
t257 = t154 * MDP(32) - t155 * MDP(33) + t295 + t296 - t298;
t173 = pkin(12) * t216 + t179;
t177 = pkin(5) * t205 - pkin(12) * t206 + t188;
t165 = -t173 * t246 + t177 * t250;
t166 = t173 * t250 + t177 * t246;
t256 = t165 * MDP(32) - t166 * MDP(33) + t287 + t290 - t292;
t255 = -MDP(23) + t259;
t239 = t250 ^ 2;
t238 = t247 ^ 2;
t237 = t246 ^ 2;
t221 = -pkin(9) * t305 + t230;
t1 = [(t163 ^ 2 + t164 ^ 2 + t184 ^ 2) * MDP(19) + t186 ^ 2 * MDP(24) + MDP(1) + (MDP(4) * t249 + 0.2e1 * MDP(5) * t253) * t308 + (t283 + 0.2e1 * t289) * t218 + (t186 * t325 + t294) * t183 + (t168 * t321 + t297) * t169 + (0.2e1 * t262 + t301) * t245 + (-0.2e1 * t199 * MDP(12) - 0.2e1 * t218 * MDP(13) + t288) * t200 + (t186 * t324 - 0.2e1 * t293 + t295 + 0.2e1 * t296 - 0.2e1 * t298) * t182 + 0.2e1 * (t221 * t245 + t236 * t318) * MDP(9) + 0.2e1 * (-pkin(1) * t308 - t223 * t245) * MDP(10) + (t181 * t218 + t190 * t200) * t327 + (-t180 * t218 + t190 * t199) * t328 + (t158 * t186 + t161 * t182) * t323 + (-t159 * t186 + t161 * t183) * t322 + (-t163 * t187 - t164 * t186) * t326 + (t154 * t182 + t156 * t168) * t320 + (-t155 * t182 + t156 * t169) * t319; t183 * t286 + t186 * t284 + t182 * t287 + t169 * t291 + (-t168 * t194 - t169 * t193) * MDP(28) + (t154 * t205 + t156 * t193 + t165 * t182 + t168 * t172) * MDP(32) + (-t168 * t205 - t182 * t193) * MDP(30) + (t169 * t205 + t182 * t194) * MDP(29) + (-t155 * t205 + t156 * t194 - t166 * t182 + t169 * t172) * MDP(33) + (-t182 * t206 - t183 * t205) * MDP(21) + (t158 * t216 + t161 * t205 + t178 * t186 + t182 * t188) * MDP(25) + (-t182 * t216 - t186 * t205) * MDP(23) + (t183 * t216 + t186 * t206) * MDP(22) + (-t159 * t216 + t161 * t206 - t179 * t186 + t183 * t188) * MDP(26) + (-t163 * t217 - t164 * t216 - t186 * t192 - t187 * t191) * MDP(18) + t221 * MDP(9) - t223 * MDP(10) + (t163 * t191 + t164 * t192 + t184 * t225) * MDP(19) + t301 + t262 + t270 * t218 + t258 * t244 + ((-pkin(2) * t200 + t190 * t248) * MDP(17) + (-pkin(2) * t199 - t190 * t252) * MDP(16) + (-t199 * t248 + t200 * t252) * MDP(12) + t248 * t288 - t271 * t218) * t241; (t191 ^ 2 + t192 ^ 2 + t225 ^ 2) * MDP(19) + t216 ^ 2 * MDP(24) + MDP(8) + (MDP(11) * t248 + 0.2e1 * MDP(12) * t252) * t309 + (t216 * t325 + t286) * t206 + (t193 * t321 + t291) * t194 + (0.2e1 * t261 + t280) * t244 + (t216 * t324 - 0.2e1 * t285 + t287 + 0.2e1 * t290 - 0.2e1 * t292) * t205 + (pkin(2) * t235 * t252 + t220 * t244) * t328 + (-pkin(2) * t309 - t222 * t244) * t327 + (-t191 * t217 - t192 * t216) * t326 + (t178 * t216 + t188 * t205) * t323 + (-t179 * t216 + t188 * t206) * t322 + (t165 * t205 + t172 * t193) * t320 + (-t166 * t205 + t172 * t194) * t319; t183 * t281 + t260 * t182 + ((-t186 * t240 - t187 * t243) * MDP(18) + (t163 * t243 + t164 * t240) * MDP(19)) * pkin(3) + (-t161 * MDP(25) + t272 * t186 - t257 + t293) * t251 + (t294 + t161 * MDP(26) + t169 * t278 + (-t168 * t250 - t169 * t246) * MDP(28) + (t168 * t233 + t316) * MDP(32) + (t169 * t233 + t315) * MDP(33) + t273 * t186 + t263 * t182) * t247 + t258; t206 * t281 + t280 + t261 + t260 * t205 + ((-t216 * t240 - t217 * t243) * MDP(18) + (t191 * t243 + t192 * t240) * MDP(19)) * pkin(3) + (-t188 * MDP(25) + t272 * t216 - t256 + t285) * t251 + (t286 + t188 * MDP(26) + t194 * t278 + (-t193 * t250 - t194 * t246) * MDP(28) + (t193 * t233 + t314) * MDP(32) + (t194 * t233 + t313) * MDP(33) + t273 * t216 + t263 * t205) * t247 - t270; MDP(15) + (t240 ^ 2 + t243 ^ 2) * MDP(19) * pkin(3) ^ 2 + (t277 - 0.2e1 * t282) * t251 + (MDP(27) * t239 + MDP(20) - 0.2e1 * t274) * t238 + (-t202 * t251 + t238 * t312) * t320 + (t203 * t251 + t238 * t311) * t319 + 0.2e1 * (-t251 * t263 + t281) * t247; t184 * MDP(19) + (-t168 * t251 - t182 * t300) * MDP(32) + (-t169 * t251 - t182 * t299) * MDP(33) + t269 * t186; t225 * MDP(19) + (-t193 * t251 - t205 * t300) * MDP(32) + (-t194 * t251 - t205 * t299) * MDP(33) + t269 * t216; 0; MDP(19); t183 * MDP(22) + t186 * MDP(24) + t158 * MDP(25) - t159 * MDP(26) + t246 * t297 + (-t168 * t246 + t169 * t250) * MDP(28) + (-pkin(5) * t168 - t315) * MDP(32) + (-pkin(5) * t169 + t316) * MDP(33) + t255 * t182; t206 * MDP(22) + t284 + t178 * MDP(25) - t179 * MDP(26) + t246 * t291 + (-t193 * t246 + t194 * t250) * MDP(28) + (-pkin(5) * t193 - t313) * MDP(32) + (-pkin(5) * t194 + t314) * MDP(33) + t255 * t205; (-t259 + t272) * t251 + (t246 * t278 + (-t237 + t239) * MDP(28) + (-pkin(5) * t246 - t311) * MDP(32) + (-pkin(5) * t250 + t312) * MDP(33) + t273) * t247; -t279 + (MDP(25) + t265) * t251; MDP(27) * t237 + 0.2e1 * pkin(5) * t265 + MDP(24) + 0.2e1 * t274; t257; t256; t268 * t247 + t266 - t277; -t264 * t247; t259; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
