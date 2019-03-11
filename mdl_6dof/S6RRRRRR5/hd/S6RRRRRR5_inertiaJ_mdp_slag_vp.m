% Calculate joint inertia matrix for
% S6RRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR5_inertiaJ_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:04:53
% EndTime: 2019-03-10 04:04:58
% DurationCPUTime: 1.76s
% Computational Cost: add. (2932->279), mult. (6473->383), div. (0->0), fcn. (7578->12), ass. (0->141)
t239 = sin(qJ(6));
t244 = cos(qJ(6));
t293 = t239 * MDP(34) + t244 * MDP(35);
t291 = MDP(37) * t244;
t253 = -0.2e1 * MDP(38) * t239 + 0.2e1 * t291;
t241 = sin(qJ(4));
t246 = cos(qJ(4));
t275 = t246 * MDP(23);
t324 = pkin(3) * (-MDP(24) * t241 + t275);
t245 = cos(qJ(5));
t240 = sin(qJ(5));
t292 = MDP(31) * t240;
t254 = (MDP(30) * t245 - t292) * pkin(4);
t238 = cos(pkin(6));
t242 = sin(qJ(3));
t247 = cos(qJ(3));
t237 = sin(pkin(6));
t243 = sin(qJ(2));
t296 = t237 * t243;
t208 = -t238 * t247 + t242 * t296;
t209 = t238 * t242 + t247 * t296;
t184 = t208 * t246 + t209 * t241;
t185 = -t208 * t241 + t209 * t246;
t169 = t245 * t184 + t185 * t240;
t170 = -t184 * t240 + t185 * t245;
t323 = t170 * MDP(27) - t169 * MDP(28);
t322 = t185 * MDP(20) - t184 * MDP(21);
t234 = t239 ^ 2;
t272 = t239 * t244 * MDP(33);
t321 = t234 * MDP(32) + 0.2e1 * t272;
t248 = cos(qJ(2));
t295 = t237 * t248;
t273 = pkin(8) * t295;
t306 = pkin(1) * t243;
t202 = t273 + (pkin(9) + t306) * t238;
t203 = (-pkin(2) * t248 - pkin(9) * t243 - pkin(1)) * t237;
t176 = -t242 * t202 + t247 * t203;
t174 = -pkin(3) * t295 - t209 * pkin(10) + t176;
t177 = t202 * t247 + t203 * t242;
t175 = -pkin(10) * t208 + t177;
t154 = t246 * t174 - t241 * t175;
t155 = t241 * t174 + t246 * t175;
t150 = -pkin(4) * t295 - t185 * pkin(11) + t154;
t151 = -pkin(11) * t184 + t155;
t145 = t245 * t150 - t240 * t151;
t146 = t240 * t150 + t245 * t151;
t157 = t239 * t170 + t244 * t295;
t158 = t244 * t170 - t239 * t295;
t285 = t158 * MDP(32);
t319 = (-t157 * t239 + t158 * t244) * MDP(33) + t239 * t285 + t323 + t293 * t169 + t145 * MDP(30) - t146 * MDP(31);
t320 = t154 * MDP(23) - t155 * MDP(24) + t319 + t322;
t318 = -(MDP(16) * t242 + MDP(17) * t247) * pkin(9) + t242 * MDP(13) + t247 * MDP(14);
t307 = pkin(9) + pkin(10);
t218 = t307 * t242;
t219 = t307 * t247;
t197 = -t218 * t241 + t219 * t246;
t215 = t241 * t242 - t246 * t247;
t183 = -pkin(11) * t215 + t197;
t196 = -t246 * t218 - t219 * t241;
t216 = t241 * t247 + t242 * t246;
t258 = -pkin(11) * t216 + t196;
t167 = t183 * t240 - t245 * t258;
t168 = t245 * t183 + t240 * t258;
t192 = t245 * t215 + t216 * t240;
t193 = -t215 * t240 + t216 * t245;
t317 = t193 * MDP(27) - t192 * MDP(28) - t167 * MDP(30) - t168 * MDP(31);
t316 = t216 * MDP(20) - t215 * MDP(21) + t196 * MDP(23) - t197 * MDP(24);
t315 = 0.2e1 * MDP(16);
t314 = -2 * MDP(19);
t313 = 0.2e1 * MDP(23);
t312 = 0.2e1 * MDP(24);
t311 = 0.2e1 * MDP(30);
t310 = 0.2e1 * MDP(31);
t309 = 0.2e1 * MDP(37);
t308 = 0.2e1 * MDP(38);
t305 = pkin(1) * t248;
t304 = pkin(3) * t241;
t303 = pkin(4) * t245;
t302 = pkin(5) * t239;
t143 = pkin(5) * t295 - t145;
t301 = t143 * t244;
t300 = t167 * t244;
t299 = t169 * t239;
t298 = t169 * t244;
t297 = t193 * t239;
t294 = t238 * t243;
t286 = t157 * MDP(35);
t284 = t158 * MDP(34);
t283 = t169 * MDP(36);
t282 = t170 * MDP(26);
t281 = t192 * MDP(36);
t227 = pkin(3) * t246 + pkin(4);
t220 = t245 * t227;
t206 = -t240 * t304 + t220;
t280 = t206 * MDP(30);
t207 = t227 * t240 + t245 * t304;
t279 = t207 * MDP(31);
t278 = t216 * MDP(18);
t277 = t242 * MDP(11);
t276 = t244 * MDP(32);
t274 = MDP(22) + MDP(29);
t228 = -pkin(3) * t247 - pkin(2);
t271 = MDP(29) + t321;
t270 = MDP(22) + t271;
t269 = -pkin(5) * t193 - pkin(12) * t192;
t204 = -pkin(5) - t206;
t205 = pkin(12) + t207;
t268 = -t192 * t205 + t193 * t204;
t225 = pkin(4) * t240 + pkin(12);
t226 = -pkin(5) - t303;
t267 = -t192 * t225 + t193 * t226;
t200 = pkin(4) * t215 + t228;
t266 = t209 * MDP(13) - t208 * MDP(14);
t261 = t244 * MDP(34) - t239 * MDP(35);
t259 = t239 * MDP(37) + t244 * MDP(38);
t223 = pkin(8) * t296;
t201 = t223 + (-pkin(2) - t305) * t238;
t256 = -MDP(26) + t261;
t255 = -MDP(29) + t279 - t280;
t235 = t244 ^ 2;
t252 = (-t234 + t235) * t193 * MDP(33) + t276 * t297 + t317 + t293 * t192;
t188 = t208 * pkin(3) + t201;
t171 = t184 * pkin(4) + t188;
t251 = t252 + t316;
t144 = -pkin(12) * t295 + t146;
t147 = t169 * pkin(5) - t170 * pkin(12) + t171;
t140 = -t144 * t239 + t147 * t244;
t141 = t144 * t244 + t147 * t239;
t250 = t140 * MDP(37) - t141 * MDP(38) + t283 + t284 - t286;
t233 = t237 ^ 2;
t232 = pkin(5) * t244;
t221 = t226 * t239;
t212 = pkin(1) * t294 + t273;
t211 = t238 * t305 - t223;
t199 = t204 * t239;
t164 = pkin(5) * t192 - pkin(12) * t193 + t200;
t161 = t167 * t239;
t153 = t164 * t239 + t168 * t244;
t152 = t164 * t244 - t168 * t239;
t142 = t143 * t239;
t1 = [t170 ^ 2 * MDP(25) + t238 ^ 2 * MDP(8) + MDP(1) + (t209 * MDP(11) - 0.2e1 * t208 * MDP(12)) * t209 + (t185 * MDP(18) + t184 * t314) * t185 + (-0.2e1 * t157 * MDP(33) + t285) * t158 + (-0.2e1 * t282 + t283 + 0.2e1 * t284 - 0.2e1 * t286) * t169 + ((MDP(4) * t243 + 0.2e1 * MDP(5) * t248) * t243 + (MDP(15) + t274) * t248 ^ 2) * t233 + 0.2e1 * (MDP(6) * t294 + (MDP(7) * t238 - t266 - t322 - t323) * t248) * t237 + (-t145 * t295 + t171 * t169) * t311 + (t146 * t295 + t171 * t170) * t310 + 0.2e1 * (t211 * t238 + t233 * t305) * MDP(9) + (-t176 * t295 + t201 * t208) * t315 + 0.2e1 * (t177 * t295 + t201 * t209) * MDP(17) + (-t154 * t295 + t188 * t184) * t313 + (t155 * t295 + t188 * t185) * t312 + 0.2e1 * (-t212 * t238 - t233 * t306) * MDP(10) + (t140 * t169 + t143 * t157) * t309 + (-t141 * t169 + t143 * t158) * t308; (-t208 * t242 + t209 * t247) * MDP(12) + (-pkin(2) * t208 - t201 * t247) * MDP(16) + (-pkin(2) * t209 + t201 * t242) * MDP(17) + (t228 * t185 + t188 * t216) * MDP(24) + (t228 * t184 + t188 * t215) * MDP(23) + (-t153 * t169 + t158 * t167) * MDP(38) + t209 * t277 + t238 * MDP(8) + (t152 * t169 + t157 * t167) * MDP(37) + (-t184 * t216 - t185 * t215) * MDP(19) + t211 * MDP(9) - t212 * MDP(10) + t185 * t278 + (t169 * MDP(30) + t170 * MDP(31)) * t200 + (t171 * MDP(30) + t250 - t282) * t192 + (t171 * MDP(31) + (-t157 * t244 - t158 * t239) * MDP(33) + t170 * MDP(25) + t158 * t276 + t259 * t143 + t256 * t169) * t193 + (MDP(6) * t243 + (MDP(7) - t316 - t317 - t318) * t248) * t237; pkin(2) * t247 * t315 + t228 * t215 * t313 + MDP(8) + t167 * t297 * t309 + (0.2e1 * t247 * MDP(12) - 0.2e1 * pkin(2) * MDP(17) + t277) * t242 + (t215 * t314 + t228 * t312 + t278) * t216 + (t152 * t309 - t153 * t308 + t200 * t311 + t281) * t192 + (0.2e1 * t256 * t192 + t200 * t310 + t300 * t308 + (MDP(32) * t235 + MDP(25) - 0.2e1 * t272) * t193) * t193; (-MDP(15) - MDP(22) + t255 - t324) * t295 + (t157 * t204 - t205 * t299 - t301) * MDP(37) + (t158 * t204 - t205 * t298 + t142) * MDP(38) + t176 * MDP(16) - t177 * MDP(17) + t266 + t320; (t239 * t268 - t300) * MDP(37) + (t244 * t268 + t161) * MDP(38) + t251 + t318; -t204 * t253 + MDP(15) + t270 - 0.2e1 * t279 + 0.2e1 * t280 + 0.2e1 * t324; (-t274 - t254) * t295 + (t157 * t226 - t225 * t299 - t301) * MDP(37) + (t158 * t226 - t225 * t298 + t142) * MDP(38) + t320; (t239 * t267 - t300) * MDP(37) + (t244 * t267 + t161) * MDP(38) + t251; (t220 + t303) * MDP(30) + (t221 + t199) * MDP(38) + (-t204 - t226) * t291 + (-pkin(4) - t227) * t292 + (t275 + (-MDP(30) * t240 - MDP(31) * t245 - MDP(24)) * t241) * pkin(3) + t270; -t226 * t253 + 0.2e1 * t254 + t270; -MDP(29) * t295 + (-pkin(5) * t157 - pkin(12) * t299 - t301) * MDP(37) + (-pkin(5) * t158 - pkin(12) * t298 + t142) * MDP(38) + t319; (t239 * t269 - t300) * MDP(37) + (t244 * t269 + t161) * MDP(38) + t252; (-t204 * t244 + t232) * MDP(37) + (t199 - t302) * MDP(38) - t255 + t321; (-t226 * t244 + t232) * MDP(37) + (t221 - t302) * MDP(38) + t254 + t271; pkin(5) * t253 + t271; t250; t152 * MDP(37) - t153 * MDP(38) + t193 * t261 + t281; -t205 * t259 + t293; -t225 * t259 + t293; -pkin(12) * t259 + t293; MDP(36);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
