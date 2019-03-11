% Calculate joint inertia matrix for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR14_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR14_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR14_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:20:39
% EndTime: 2019-03-09 20:20:43
% DurationCPUTime: 1.31s
% Computational Cost: add. (1507->268), mult. (3305->367), div. (0->0), fcn. (3557->10), ass. (0->131)
t295 = pkin(4) + pkin(9);
t211 = sin(pkin(6));
t220 = cos(qJ(2));
t266 = t211 * t220;
t212 = cos(pkin(6));
t215 = sin(qJ(3));
t219 = cos(qJ(3));
t216 = sin(qJ(2));
t267 = t211 * t216;
t181 = -t212 * t219 + t215 * t267;
t214 = sin(qJ(5));
t218 = cos(qJ(5));
t164 = t181 * t218 + t214 * t266;
t165 = t181 * t214 - t218 * t266;
t213 = sin(qJ(6));
t217 = cos(qJ(6));
t149 = -t217 * t164 + t165 * t213;
t150 = t164 * t213 + t165 * t217;
t294 = t150 * MDP(31) - t149 * MDP(32);
t263 = t218 * t219;
t264 = t214 * t219;
t177 = t213 * t264 - t217 * t263;
t190 = t213 * t218 + t214 * t217;
t178 = t190 * t219;
t293 = -t178 * MDP(31) + t177 * MDP(32);
t191 = -t213 * t214 + t217 * t218;
t281 = pkin(3) + pkin(10);
t275 = -pkin(11) - t281;
t193 = t275 * t214;
t194 = t275 * t218;
t237 = t191 * MDP(31) - t190 * MDP(32) + (-t193 * t213 + t194 * t217) * MDP(34) - (t193 * t217 + t194 * t213) * MDP(35);
t292 = (-MDP(27) * t281 + MDP(24)) * t218 + (MDP(28) * t281 - MDP(25)) * t214 + t237;
t246 = pkin(8) * t266;
t280 = pkin(1) * t216;
t174 = t246 + (pkin(9) + t280) * t212;
t175 = (-pkin(2) * t220 - pkin(9) * t216 - pkin(1)) * t211;
t154 = -t215 * t174 + t175 * t219;
t200 = pkin(3) * t266;
t153 = -t154 + t200;
t182 = t212 * t215 + t219 * t267;
t144 = pkin(4) * t182 + pkin(10) * t266 + t153;
t199 = pkin(8) * t267;
t279 = pkin(1) * t220;
t173 = t199 + (-pkin(2) - t279) * t212;
t227 = -qJ(4) * t182 + t173;
t145 = t281 * t181 + t227;
t138 = t218 * t144 - t145 * t214;
t139 = t144 * t214 + t145 * t218;
t291 = t138 * MDP(27) - t139 * MDP(28);
t241 = -qJ(4) * t215 - pkin(2);
t195 = -pkin(3) * t219 + t241;
t290 = -pkin(2) * MDP(17) - t195 * MDP(20) + t293;
t289 = 2 * MDP(18);
t288 = 2 * MDP(19);
t287 = 0.2e1 * MDP(20);
t286 = 0.2e1 * MDP(27);
t285 = 0.2e1 * MDP(28);
t284 = -2 * MDP(30);
t283 = 0.2e1 * MDP(34);
t282 = 0.2e1 * MDP(35);
t278 = pkin(5) * t182;
t277 = pkin(5) * t215;
t276 = pkin(5) * t217;
t274 = MDP(16) * pkin(2);
t273 = MDP(21) * pkin(3);
t271 = pkin(9) * MDP(21);
t137 = pkin(11) * t164 + t139;
t270 = t137 * t217;
t189 = -t281 * t219 + t241;
t240 = pkin(11) * t219 - t189;
t196 = t295 * t215;
t268 = t196 * t214;
t158 = -t240 * t218 + t268;
t269 = t158 * t217;
t265 = t212 * MDP(8);
t155 = t219 * t174 + t215 * t175;
t262 = t191 * MDP(34) - t190 * MDP(35);
t197 = t295 * t219;
t261 = MDP(15) * t220;
t260 = MDP(21) * t195;
t259 = MDP(21) * pkin(9) ^ 2;
t258 = MDP(27) * t214;
t257 = MDP(28) * t218;
t256 = MDP(29) * t178;
t255 = MDP(29) * t191;
t252 = t153 * MDP(21);
t251 = t165 * MDP(22);
t249 = t214 * MDP(22);
t248 = MDP(17) - MDP(20);
t247 = MDP(26) + MDP(33);
t245 = qJ(4) * t266;
t244 = t182 * MDP(33) + t294;
t243 = t215 * MDP(33) + t293;
t242 = t218 * t214 * MDP(23);
t239 = MDP(11) + t247;
t238 = MDP(19) - t273;
t136 = -pkin(11) * t165 + t138 + t278;
t133 = t217 * t136 - t137 * t213;
t192 = t218 * t196;
t157 = t240 * t214 + t192 + t277;
t141 = t217 * t157 - t158 * t213;
t236 = t165 * MDP(24) + t164 * MDP(25);
t235 = -MDP(24) * t214 - MDP(25) * t218;
t160 = -t189 * t214 + t192;
t161 = t189 * t218 + t268;
t234 = t160 * MDP(27) - t161 * MDP(28);
t233 = MDP(27) * t218 - MDP(28) * t214;
t134 = t136 * t213 + t270;
t232 = t133 * MDP(34) - t134 * MDP(35);
t142 = t157 * t213 + t269;
t231 = MDP(34) * t141 - MDP(35) * t142;
t152 = t245 - t155;
t230 = MDP(12) + t235;
t229 = MDP(18) + t233;
t228 = (MDP(34) * t217 - MDP(35) * t213) * pkin(5);
t146 = -pkin(4) * t181 - t152;
t225 = t229 + t262;
t224 = -t181 * MDP(12) - MDP(13) * t266 + t236 + t294;
t223 = -pkin(3) * MDP(18) + MDP(13) + t292;
t210 = t219 ^ 2;
t209 = t218 ^ 2;
t208 = t215 ^ 2;
t207 = t214 ^ 2;
t206 = t211 ^ 2;
t202 = pkin(5) * t214 + qJ(4);
t184 = t212 * t280 + t246;
t183 = t212 * t279 - t199;
t180 = pkin(5) * t263 + t197;
t151 = pkin(3) * t181 + t227;
t140 = -pkin(5) * t164 + t146;
t1 = [MDP(1) + (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) * MDP(21) + t206 * t216 ^ 2 * MDP(4) + (0.2e1 * MDP(6) * t267 + t265) * t212 + (0.2e1 * t164 * MDP(23) + t251) * t165 + (MDP(29) * t150 + t149 * t284) * t150 + t239 * t182 ^ 2 + (0.2e1 * MDP(5) * t216 + t261) * t206 * t220 + (t133 * t182 + t140 * t149) * t283 + (t138 * t182 - t146 * t164) * t286 + (t152 * t181 + t153 * t182) * t289 + (-t134 * t182 + t140 * t150) * t282 + (-t139 * t182 + t146 * t165) * t285 + 0.2e1 * (-t184 * t212 - t206 * t280) * MDP(10) + 0.2e1 * (-t154 * t266 + t173 * t181) * MDP(16) + (-t151 * t181 - t153 * t266) * t288 + 0.2e1 * (t155 * t266 + t173 * t182) * MDP(17) + (-t151 * t182 + t152 * t266) * t287 + 0.2e1 * (t183 * t212 + t206 * t279) * MDP(9) + 0.2e1 * (t181 * MDP(14) + MDP(7) * t212) * t266 + 0.2e1 * t224 * t182; (-t140 * t177 + t180 * t149) * MDP(34) + (-t140 * t178 + t180 * t150) * MDP(35) + t151 * t260 + (t178 * t149 + t150 * t177) * MDP(30) + t183 * MDP(9) - t184 * MDP(10) + t265 - t150 * t256 + (MDP(6) * t216 + MDP(7) * t220) * t211 + (-t164 * MDP(27) + t165 * MDP(28)) * t197 + (-MDP(19) * t195 - t274) * t181 + (t231 + t234 + t290) * t182 + (-t165 * t249 - MDP(14) * t266 - t152 * MDP(18) + (-t164 * t214 - t165 * t218) * MDP(23) - t173 * MDP(16) + t151 * MDP(19) + t233 * t146 + t230 * t182 + (-MDP(18) * t181 - MDP(21) * t152 + t248 * t266) * pkin(9)) * t219 + (t173 * MDP(17) + t153 * MDP(18) - t151 * MDP(20) + t239 * t182 + (t182 * MDP(18) + t252 + (MDP(16) - MDP(19)) * t266) * pkin(9) + t224 + t232 + t291) * t215; 0.2e1 * t219 * t274 + MDP(8) + (t219 * t288 + t260) * t195 + (MDP(22) * t207 + 0.2e1 * t242 + t259) * t210 - (0.2e1 * MDP(30) * t177 - t256) * t178 + (t239 + t259) * t208 + 0.2e1 * (t230 * t219 + t290) * t215 + (t160 * t215 + t197 * t263) * t286 + (-t161 * t215 - t197 * t264) * t285 + (t141 * t215 - t177 * t180) * t283 + (-t142 * t215 - t178 * t180) * t282 + (t208 + t210) * pkin(9) * t289; -t211 * t261 + t154 * MDP(16) - t155 * MDP(17) + (-t154 + 0.2e1 * t200) * MDP(19) + (-0.2e1 * t245 + t155) * MDP(20) + (-pkin(3) * t153 - qJ(4) * t152) * MDP(21) + t218 * t251 + (t164 * t218 - t165 * t214) * MDP(23) + (-qJ(4) * t164 + t146 * t214) * MDP(27) + (qJ(4) * t165 + t146 * t218) * MDP(28) + t150 * t255 + (-t149 * t191 - t150 * t190) * MDP(30) + (t140 * t190 + t149 * t202) * MDP(34) + (t140 * t191 + t150 * t202) * MDP(35) + (-qJ(4) * MDP(18) - MDP(14)) * t181 + t223 * t182; -t178 * t255 + (t177 * t191 + t178 * t190) * MDP(30) + (-t177 * t202 + t180 * t190) * MDP(34) + (-t178 * t202 + t180 * t191) * MDP(35) + (t257 + t258) * t197 + (MDP(14) - t218 * t249 + (t207 - t209) * MDP(23) - t248 * pkin(9) + (t229 + t271) * qJ(4)) * t219 + ((-MDP(16) + t238) * pkin(9) + t223) * t215; -0.2e1 * t242 + t202 * t190 * t283 + t209 * MDP(22) + MDP(15) + (-(2 * MDP(19)) + t273) * pkin(3) + (t190 * t284 + t202 * t282 + t255) * t191 + (MDP(21) * qJ(4) + 0.2e1 * t257 + 0.2e1 * t258 + t287) * qJ(4); -MDP(19) * t266 + t225 * t182 + t252; (t225 + t271) * t215; t238; MDP(21); t182 * MDP(26) + (t182 * t276 + t133) * MDP(34) + (-t270 + (-t136 - t278) * t213) * MDP(35) + t236 + t244 + t291; t215 * MDP(26) + (t215 * t276 + t141) * MDP(34) + (-t269 + (-t157 - t277) * t213) * MDP(35) + t235 * t219 + t234 + t243; t292; t233 + t262; 0.2e1 * t228 + t247; t232 + t244; t231 + t243; t237; t262; MDP(33) + t228; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
