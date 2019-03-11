% Calculate joint inertia matrix for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR11_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:46:09
% EndTime: 2019-03-09 07:46:16
% DurationCPUTime: 2.00s
% Computational Cost: add. (3198->275), mult. (8503->412), div. (0->0), fcn. (9941->14), ass. (0->128)
t233 = sin(qJ(5));
t237 = cos(qJ(5));
t232 = sin(qJ(6));
t236 = cos(qJ(6));
t206 = t232 * t233 - t236 * t237;
t207 = t232 * t237 + t233 * t236;
t293 = pkin(11) + pkin(12);
t213 = t293 * t233;
t214 = t293 * t237;
t255 = t207 * MDP(31) - t206 * MDP(32) + (-t213 * t236 - t214 * t232) * MDP(34) - (-t213 * t232 + t214 * t236) * MDP(35);
t306 = -(t233 * MDP(27) + t237 * MDP(28)) * pkin(11) + t233 * MDP(24) + t237 * MDP(25) + t255;
t234 = sin(qJ(4));
t238 = cos(qJ(4));
t210 = -pkin(4) * t238 - pkin(11) * t234 - pkin(3);
t205 = t237 * t210;
t276 = t234 * t237;
t289 = pkin(10) * t233;
t178 = -pkin(12) * t276 + t205 + (-pkin(5) - t289) * t238;
t287 = pkin(10) * t238;
t260 = t237 * t287;
t184 = t260 + (-pkin(12) * t234 + t210) * t233;
t165 = t236 * t178 - t184 * t232;
t166 = t178 * t232 + t184 * t236;
t196 = t207 * t234;
t194 = t196 * MDP(32);
t197 = t206 * t234;
t195 = t197 * MDP(31);
t305 = MDP(34) * t165 - MDP(35) * t166 - t194 - t195;
t304 = -MDP(18) + t306;
t191 = -t233 * t287 + t205;
t192 = t210 * t233 + t260;
t303 = t191 * MDP(27) - t192 * MDP(28) + t305;
t242 = (MDP(34) * t236 - MDP(35) * t232) * pkin(5);
t226 = sin(pkin(13));
t228 = sin(pkin(6));
t231 = cos(pkin(6));
t235 = sin(qJ(3));
t239 = cos(qJ(3));
t229 = cos(pkin(13));
t230 = cos(pkin(7));
t277 = t229 * t230;
t227 = sin(pkin(7));
t280 = t227 * t235;
t180 = t231 * t280 + (t226 * t239 + t235 * t277) * t228;
t278 = t228 * t229;
t198 = -t227 * t278 + t230 * t231;
t173 = t180 * t238 + t198 * t234;
t259 = t228 * t277;
t279 = t227 * t239;
t281 = t226 * t228;
t179 = -t231 * t279 + t235 * t281 - t239 * t259;
t161 = t173 * t233 - t179 * t237;
t162 = t173 * t237 + t179 * t233;
t152 = t236 * t161 + t162 * t232;
t153 = -t161 * t232 + t162 * t236;
t301 = t153 * MDP(31) - t152 * MDP(32);
t291 = pkin(1) * t231;
t217 = t229 * t291;
t181 = pkin(2) * t231 + t217 + (-pkin(9) * t230 - qJ(2)) * t281;
t190 = (-pkin(9) * t226 * t227 - pkin(2) * t229 - pkin(1)) * t228;
t169 = -t181 * t227 + t230 * t190;
t154 = pkin(3) * t179 - pkin(10) * t180 + t169;
t200 = qJ(2) * t278 + t226 * t291;
t177 = (t227 * t231 + t259) * pkin(9) + t200;
t254 = t181 * t230 + t190 * t227;
t159 = t177 * t239 + t235 * t254;
t157 = t198 * pkin(10) + t159;
t147 = t154 * t234 + t157 * t238;
t145 = pkin(11) * t179 + t147;
t158 = -t235 * t177 + t239 * t254;
t156 = -t198 * pkin(3) - t158;
t172 = t180 * t234 - t198 * t238;
t149 = t172 * pkin(4) - t173 * pkin(11) + t156;
t141 = -t145 * t233 + t237 * t149;
t142 = t145 * t237 + t149 * t233;
t299 = -t141 * MDP(27) + t142 * MDP(28);
t298 = 0.2e1 * MDP(27);
t297 = 0.2e1 * MDP(28);
t296 = -2 * MDP(30);
t295 = 0.2e1 * MDP(34);
t294 = 0.2e1 * MDP(35);
t221 = t228 ^ 2;
t292 = pkin(1) * t221;
t290 = pkin(5) * t172;
t288 = pkin(10) * t237;
t286 = pkin(3) * MDP(20);
t285 = pkin(10) * MDP(21);
t140 = -pkin(12) * t161 + t142;
t284 = t140 * t236;
t146 = t154 * t238 - t234 * t157;
t144 = -pkin(4) * t179 - t146;
t283 = t144 * t233;
t282 = t144 * t237;
t202 = t230 * t234 + t238 * t280;
t182 = -t233 * t202 - t237 * t279;
t183 = t237 * t202 - t233 * t279;
t167 = t182 * t236 - t183 * t232;
t168 = t182 * t232 + t183 * t236;
t275 = t167 * MDP(34) - t168 * MDP(35);
t273 = MDP(21) * t234;
t272 = MDP(29) * t197;
t271 = MDP(29) * t207;
t269 = MDP(34) * t206;
t265 = t162 * MDP(22);
t264 = t173 * MDP(17);
t263 = t179 * MDP(19);
t262 = t237 * MDP(22);
t261 = MDP(26) + MDP(33);
t258 = t172 * MDP(33) + t301;
t257 = MDP(23) * t233 * t237;
t256 = -pkin(10) * MDP(20) + MDP(17);
t139 = -pkin(12) * t162 + t141 + t290;
t136 = t236 * t139 - t140 * t232;
t253 = t162 * MDP(24) - t161 * MDP(25);
t252 = t237 * MDP(24) - t233 * MDP(25);
t249 = MDP(27) * t237 - MDP(28) * t233;
t137 = t139 * t232 + t284;
t247 = -t136 * MDP(34) + t137 * MDP(35);
t244 = -MDP(16) + t252;
t241 = t173 * MDP(16) - t253 - t301;
t224 = t237 ^ 2;
t222 = t233 ^ 2;
t220 = -pkin(5) * t237 - pkin(4);
t209 = (pkin(5) * t233 + pkin(10)) * t234;
t201 = -t238 * t230 + t234 * t280;
t199 = -qJ(2) * t281 + t217;
t143 = pkin(5) * t161 + t144;
t1 = [MDP(1) + t173 ^ 2 * MDP(15) + t198 ^ 2 * MDP(12) + (pkin(1) ^ 2 * t221 + t199 ^ 2 + t200 ^ 2) * MDP(7) + (0.2e1 * t198 * MDP(10) + MDP(8) * t180) * t180 + t261 * t172 ^ 2 + (-0.2e1 * t161 * MDP(23) + t265) * t162 + (MDP(29) * t153 + t152 * t296) * t153 + (-0.2e1 * t198 * MDP(11) - 0.2e1 * t180 * MDP(9) + t263 + 0.2e1 * t264) * t179 + 0.2e1 * (-t179 * MDP(18) - t241) * t172 + 0.2e1 * (t199 * t231 + t229 * t292) * MDP(4) + 0.2e1 * (-t200 * t231 - t226 * t292) * MDP(5) + 0.2e1 * (t158 * t198 + t169 * t179) * MDP(13) + 0.2e1 * (-t159 * t198 + t169 * t180) * MDP(14) + 0.2e1 * (t146 * t179 + t156 * t172) * MDP(20) + 0.2e1 * (-t147 * t179 + t156 * t173) * MDP(21) + (t141 * t172 + t144 * t161) * t298 + (t136 * t172 + t143 * t152) * t295 + (-t142 * t172 + t144 * t162) * t297 + (-t137 * t172 + t143 * t153) * t294 + 0.2e1 * (-t199 * t226 + t200 * t229) * MDP(6) * t228; (t230 * t179 + t198 * t279) * MDP(13) + (t180 * t230 - t198 * t280) * MDP(14) + (-t172 * t279 - t201 * t179) * MDP(20) + (-t173 * t279 - t202 * t179) * MDP(21) + (t161 * t201 + t172 * t182) * MDP(27) + (t162 * t201 - t172 * t183) * MDP(28) + (t152 * t201 + t167 * t172) * MDP(34) + (t153 * t201 - t168 * t172) * MDP(35) + (-MDP(4) * t229 + MDP(5) * t226 - MDP(7) * pkin(1)) * t228; MDP(7); -t153 * t272 + t158 * MDP(13) - t159 * MDP(14) - t179 * MDP(11) + t180 * MDP(10) + (t197 * t152 - t153 * t196) * MDP(30) + t198 * MDP(12) + (t143 * t196 + t209 * t152) * MDP(34) + (-t143 * t197 + t209 * t153) * MDP(35) - pkin(3) * t173 * MDP(21) + (-t286 + t303) * t172 + (-t156 * MDP(20) + (MDP(18) - t285) * t179 - t261 * t172 + t241 + t247 + t299) * t238 + (t162 * t262 + t173 * MDP(15) + (-t161 * t237 - t162 * t233) * MDP(23) + t156 * MDP(21) + (pkin(10) * t161 + t283) * MDP(27) + (pkin(10) * t162 + t282) * MDP(28) + t256 * t179 + t244 * t172) * t234; (t201 * t233 * t234 - t182 * t238) * MDP(27) + (t183 * t238 + t201 * t276) * MDP(28) + (-t167 * t238 + t196 * t201) * MDP(34) + (t168 * t238 - t197 * t201) * MDP(35) + (-t235 * MDP(14) + (MDP(20) * t238 + MDP(13) - t273) * t239) * t227; t196 * t209 * t295 - 0.2e1 * pkin(3) * t273 + MDP(12) - (t196 * t296 + t209 * t294 - t272) * t197 + (MDP(22) * t224 + t288 * t297 + t289 * t298 + MDP(15) - 0.2e1 * t257) * t234 ^ 2 + (-t165 * t295 + t166 * t294 - t191 * t298 + t192 * t297 - 0.2e1 * t244 * t234 + t261 * t238 + 0.2e1 * t194 + 0.2e1 * t195 + 0.2e1 * t286) * t238; t264 + t263 + t146 * MDP(20) - t147 * MDP(21) + t233 * t265 + (-t161 * t233 + t162 * t237) * MDP(23) + (-pkin(4) * t161 - t282) * MDP(27) + (-pkin(4) * t162 + t283) * MDP(28) + t153 * t271 + (-t152 * t207 - t153 * t206) * MDP(30) + (t143 * t206 + t152 * t220) * MDP(34) + (t143 * t207 + t153 * t220) * MDP(35) + t304 * t172; -MDP(21) * t202 + (MDP(35) * t207 - MDP(20) - t249 + t269) * t201; -t197 * t271 + (-t196 * t207 + t197 * t206) * MDP(30) + (t196 * t220 + t206 * t209) * MDP(34) + (-t197 * t220 + t207 * t209) * MDP(35) + (-t285 - t304) * t238 + (t233 * t262 + (-t222 + t224) * MDP(23) + (-pkin(4) * t233 - t288) * MDP(27) + (-pkin(4) * t237 + t289) * MDP(28) + t256) * t234; 0.2e1 * t257 + 0.2e1 * t220 * t269 + MDP(22) * t222 + MDP(19) + 0.2e1 * t249 * pkin(4) + (t206 * t296 + t220 * t294 + t271) * t207; t172 * MDP(26) + (t236 * t290 + t136) * MDP(34) + (-t284 + (-t139 - t290) * t232) * MDP(35) + t253 + t258 - t299; MDP(27) * t182 - MDP(28) * t183 + t275; (-t261 - t242) * t238 + t252 * t234 + t303; t306; 0.2e1 * t242 + t261; -t247 + t258; t275; -t238 * MDP(33) + t305; t255; MDP(33) + t242; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
