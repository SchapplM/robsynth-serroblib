% Calculate joint inertia matrix for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR11_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:23:38
% EndTime: 2019-03-09 23:23:42
% DurationCPUTime: 1.77s
% Computational Cost: add. (2833->292), mult. (6386->428), div. (0->0), fcn. (7234->12), ass. (0->135)
t250 = sin(qJ(3));
t253 = cos(qJ(3));
t228 = -pkin(3) * t253 - pkin(10) * t250 - pkin(2);
t252 = cos(qJ(4));
t223 = t252 * t228;
t290 = t250 * t252;
t249 = sin(qJ(4));
t303 = pkin(9) * t249;
t194 = -qJ(5) * t290 + t223 + (-pkin(4) - t303) * t253;
t301 = pkin(9) * t253;
t272 = t252 * t301;
t199 = t272 + (-qJ(5) * t250 + t228) * t249;
t244 = sin(pkin(12));
t246 = cos(pkin(12));
t172 = t194 * t246 - t199 * t244;
t291 = t249 * t250;
t211 = -t244 * t291 + t246 * t290;
t166 = -pkin(5) * t253 - pkin(11) * t211 + t172;
t173 = t194 * t244 + t199 * t246;
t222 = t244 * t252 + t246 * t249;
t210 = t222 * t250;
t169 = -pkin(11) * t210 + t173;
t248 = sin(qJ(6));
t307 = cos(qJ(6));
t149 = t166 * t307 - t169 * t248;
t150 = t166 * t248 + t169 * t307;
t184 = t210 * t307 + t211 * t248;
t182 = t184 * MDP(30);
t185 = -t210 * t248 + t211 * t307;
t183 = t185 * MDP(29);
t322 = MDP(32) * t149 - MDP(33) * t150 - t182 + t183;
t300 = -qJ(5) - pkin(10);
t229 = t300 * t249;
t230 = t300 * t252;
t200 = t229 * t246 + t230 * t244;
t186 = -pkin(11) * t222 + t200;
t201 = t229 * t244 - t230 * t246;
t221 = -t244 * t249 + t246 * t252;
t187 = pkin(11) * t221 + t201;
t192 = -t221 * t307 + t222 * t248;
t193 = t221 * t248 + t222 * t307;
t268 = t193 * MDP(29) - t192 * MDP(30) + (t186 * t307 - t187 * t248) * MDP(32) - (t186 * t248 + t187 * t307) * MDP(33);
t321 = t249 * MDP(20) + t252 * MDP(21) + t268 - (MDP(23) * t249 + MDP(24) * t252) * pkin(10);
t320 = -MDP(14) + t321;
t204 = -t249 * t301 + t223;
t205 = t228 * t249 + t272;
t319 = t204 * MDP(23) - t205 * MDP(24) + t322;
t317 = -0.2e1 * t250;
t245 = sin(pkin(6));
t254 = cos(qJ(2));
t293 = t245 * t254;
t247 = cos(pkin(6));
t251 = sin(qJ(2));
t294 = t245 * t251;
t218 = t247 * t250 + t253 * t294;
t197 = t218 * t249 + t252 * t293;
t198 = t218 * t252 - t249 * t293;
t174 = -t197 * t246 - t198 * t244;
t175 = -t197 * t244 + t198 * t246;
t159 = -t174 * t307 + t175 * t248;
t160 = t174 * t248 + t175 * t307;
t316 = t160 * MDP(29) - t159 * MDP(30);
t232 = pkin(8) * t294;
t305 = pkin(1) * t254;
t207 = t232 + (-pkin(2) - t305) * t247;
t217 = -t247 * t253 + t250 * t294;
t177 = pkin(3) * t217 - pkin(10) * t218 + t207;
t273 = pkin(8) * t293;
t306 = pkin(1) * t251;
t208 = t273 + (pkin(9) + t306) * t247;
t209 = (-pkin(2) * t254 - pkin(9) * t251 - pkin(1)) * t245;
t181 = t208 * t253 + t209 * t250;
t179 = -pkin(10) * t293 + t181;
t161 = t177 * t252 - t179 * t249;
t162 = t177 * t249 + t179 * t252;
t314 = -t161 * MDP(23) + t162 * MDP(24);
t313 = 0.2e1 * MDP(23);
t312 = 0.2e1 * MDP(24);
t311 = 2 * MDP(25);
t310 = -2 * MDP(28);
t309 = 0.2e1 * MDP(32);
t308 = 0.2e1 * MDP(33);
t304 = pkin(4) * t244;
t302 = pkin(9) * t252;
t299 = pkin(2) * MDP(16);
t298 = pkin(2) * MDP(17);
t297 = pkin(9) * MDP(17);
t180 = -t208 * t250 + t209 * t253;
t178 = pkin(3) * t293 - t180;
t296 = t178 * t249;
t295 = t178 * t252;
t292 = t247 * MDP(8);
t289 = t251 * MDP(6);
t153 = pkin(4) * t217 - qJ(5) * t198 + t161;
t156 = -qJ(5) * t197 + t162;
t148 = t153 * t244 + t156 * t246;
t227 = pkin(4) * t291 + pkin(9) * t250;
t287 = MDP(15) * t254;
t236 = pkin(4) * t246 + pkin(5);
t213 = t236 * t307 - t248 * t304;
t286 = MDP(32) * t213;
t283 = t160 * MDP(27);
t280 = t192 * MDP(32);
t279 = t193 * MDP(27);
t278 = t198 * MDP(18);
t214 = t236 * t248 + t304 * t307;
t277 = t214 * MDP(33);
t276 = t218 * MDP(13);
t275 = t252 * MDP(18);
t274 = MDP(22) + MDP(31);
t271 = t217 * MDP(31) + t316;
t237 = -pkin(4) * t252 - pkin(3);
t270 = t249 * t252 * MDP(19);
t269 = pkin(9) * MDP(16) - MDP(13);
t147 = t153 * t246 - t156 * t244;
t145 = pkin(5) * t217 - pkin(11) * t175 + t147;
t146 = pkin(11) * t174 + t148;
t142 = t145 * t307 - t146 * t248;
t267 = t198 * MDP(20) - t197 * MDP(21);
t266 = MDP(20) * t252 - MDP(21) * t249;
t143 = t145 * t248 + t146 * t307;
t262 = -t142 * MDP(32) + t143 * MDP(33);
t260 = -MDP(12) + t266;
t259 = MDP(31) - t277 + t286;
t170 = pkin(4) * t197 + t178;
t257 = t218 * MDP(12) - t267 - t316;
t242 = t252 ^ 2;
t240 = t249 ^ 2;
t239 = t245 ^ 2;
t220 = t247 * t306 + t273;
t219 = t247 * t305 - t232;
t206 = -pkin(5) * t221 + t237;
t195 = pkin(5) * t210 + t227;
t154 = -pkin(5) * t174 + t170;
t1 = [(t147 ^ 2 + t148 ^ 2 + t170 ^ 2) * MDP(26) + t218 ^ 2 * MDP(11) + t239 * t251 ^ 2 * MDP(4) + MDP(1) + (0.2e1 * t245 * t289 + t292) * t247 + t274 * t217 ^ 2 + (-0.2e1 * MDP(19) * t197 + t278) * t198 + (t159 * t310 + t283) * t160 + (0.2e1 * MDP(5) * t251 + t287) * t239 * t254 + 0.2e1 * (t219 * t247 + t239 * t305) * MDP(9) + 0.2e1 * (-t180 * t293 + t207 * t217) * MDP(16) + 0.2e1 * (t181 * t293 + t207 * t218) * MDP(17) + 0.2e1 * (-t220 * t247 - t239 * t306) * MDP(10) + (t161 * t217 + t178 * t197) * t313 + (t142 * t217 + t154 * t159) * t309 + (-t162 * t217 + t178 * t198) * t312 + (-t143 * t217 + t154 * t160) * t308 + (-t147 * t175 + t148 * t174) * t311 + 0.2e1 * (MDP(7) * t247 - t276) * t293 + 0.2e1 * (MDP(14) * t293 - t257) * t217; (-t159 * t185 - t160 * t184) * MDP(28) + (-t147 * t211 - t148 * t210 - t172 * t175 + t173 * t174) * MDP(25) + t219 * MDP(9) - t220 * MDP(10) + (t147 * t172 + t148 * t173 + t170 * t227) * MDP(26) + t292 + (t154 * t184 + t159 * t195) * MDP(32) + (t154 * t185 + t160 * t195) * MDP(33) - t218 * t298 + t185 * t283 + (MDP(7) * t254 + t289) * t245 + (-t299 + t319) * t217 + (-t207 * MDP(16) + (-MDP(14) + t297) * t293 - t274 * t217 + t257 + t262 + t314) * t253 + (t198 * t275 + (-t197 * t252 - t198 * t249) * MDP(19) + (pkin(9) * t197 + t296) * MDP(23) + (pkin(9) * t198 + t295) * MDP(24) + t207 * MDP(17) + t218 * MDP(11) + t269 * t293 + t260 * t217) * t250; MDP(8) + t298 * t317 + (t172 ^ 2 + t173 ^ 2 + t227 ^ 2) * MDP(26) + (-t172 * t211 - t173 * t210) * t311 + t184 * t195 * t309 + (MDP(27) * t185 + t184 * t310 + t195 * t308) * t185 + (MDP(18) * t242 + t302 * t312 + t303 * t313 + MDP(11) - 0.2e1 * t270) * t250 ^ 2 + (-t149 * t309 + t150 * t308 - t204 * t313 + t205 * t312 + t253 * t274 + t260 * t317 + 0.2e1 * t182 - 0.2e1 * t183 + 0.2e1 * t299) * t253; t276 - t245 * t287 + t180 * MDP(16) - t181 * MDP(17) + t249 * t278 + (-t197 * t249 + t198 * t252) * MDP(19) + (-pkin(3) * t197 - t295) * MDP(23) + (-pkin(3) * t198 + t296) * MDP(24) + (-t147 * t222 + t148 * t221 + t174 * t201 - t175 * t200) * MDP(25) + (t147 * t200 + t148 * t201 + t170 * t237) * MDP(26) + t160 * t279 + (-t159 * t193 - t160 * t192) * MDP(28) + (t154 * t192 + t159 * t206) * MDP(32) + (t154 * t193 + t160 * t206) * MDP(33) + t320 * t217; (-t172 * t222 + t173 * t221 - t200 * t211 - t201 * t210) * MDP(25) + (t172 * t200 + t173 * t201 + t227 * t237) * MDP(26) + t185 * t279 + (-t184 * t193 - t185 * t192) * MDP(28) + (t184 * t206 + t192 * t195) * MDP(32) + (t185 * t206 + t193 * t195) * MDP(33) + (-t297 - t320) * t253 + (t249 * t275 + (-t240 + t242) * MDP(19) + (-pkin(3) * t249 - t302) * MDP(23) + (-pkin(3) * t252 + t303) * MDP(24) - t269) * t250; MDP(15) + t240 * MDP(18) + 0.2e1 * t270 + (-t200 * t222 + t201 * t221) * t311 + (t200 ^ 2 + t201 ^ 2 + t237 ^ 2) * MDP(26) + 0.2e1 * t206 * t280 + 0.2e1 * (MDP(23) * t252 - MDP(24) * t249) * pkin(3) + (t192 * t310 + t206 * t308 + t279) * t193; t217 * MDP(22) + (t213 * t217 + t142) * MDP(32) + (-t214 * t217 - t143) * MDP(33) + ((t174 * t244 - t175 * t246) * MDP(25) + (t147 * t246 + t148 * t244) * MDP(26)) * pkin(4) + t267 + t271 - t314; (-MDP(22) - t259) * t253 + t266 * t250 + ((-t210 * t244 - t211 * t246) * MDP(25) + (t172 * t246 + t173 * t244) * MDP(26)) * pkin(4) + t319; ((t221 * t244 - t222 * t246) * MDP(25) + (t200 * t246 + t201 * t244) * MDP(26)) * pkin(4) + t321; (t244 ^ 2 + t246 ^ 2) * MDP(26) * pkin(4) ^ 2 + 0.2e1 * t286 - 0.2e1 * t277 + t274; MDP(26) * t170 + MDP(32) * t159 + MDP(33) * t160; MDP(26) * t227 + MDP(32) * t184 + MDP(33) * t185; MDP(26) * t237 + MDP(33) * t193 + t280; 0; MDP(26); -t262 + t271; -t253 * MDP(31) + t322; t268; t259; 0; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
