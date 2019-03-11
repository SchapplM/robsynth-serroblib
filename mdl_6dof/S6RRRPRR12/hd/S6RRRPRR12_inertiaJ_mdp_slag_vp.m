% Calculate joint inertia matrix for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR12_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:46:34
% EndTime: 2019-03-09 19:46:40
% DurationCPUTime: 1.82s
% Computational Cost: add. (2762->301), mult. (6304->428), div. (0->0), fcn. (7197->12), ass. (0->133)
t245 = sin(qJ(3));
t249 = cos(qJ(3));
t222 = -pkin(3) * t249 - qJ(4) * t245 - pkin(2);
t241 = cos(pkin(12));
t217 = t241 * t222;
t239 = sin(pkin(12));
t301 = pkin(9) * t239;
t191 = -pkin(10) * t241 * t245 + t217 + (-pkin(4) - t301) * t249;
t300 = pkin(9) * t249;
t201 = t239 * t222 + t241 * t300;
t294 = t239 * t245;
t196 = -pkin(10) * t294 + t201;
t244 = sin(qJ(5));
t248 = cos(qJ(5));
t168 = t248 * t191 - t196 * t244;
t169 = t191 * t244 + t196 * t248;
t219 = t239 * t244 - t248 * t241;
t208 = t219 * t245;
t162 = -pkin(5) * t249 + pkin(11) * t208 + t168;
t220 = t239 * t248 + t241 * t244;
t207 = t220 * t245;
t165 = -pkin(11) * t207 + t169;
t243 = sin(qJ(6));
t247 = cos(qJ(6));
t147 = t247 * t162 - t165 * t243;
t148 = t162 * t243 + t165 * t247;
t260 = t147 * MDP(34) - t148 * MDP(35);
t322 = t168 * MDP(27) - t169 * MDP(28) + t260;
t258 = (MDP(34) * t247 - MDP(35) * t243) * pkin(5);
t321 = pkin(9) * MDP(21);
t240 = sin(pkin(6));
t250 = cos(qJ(2));
t292 = t240 * t250;
t242 = cos(pkin(6));
t246 = sin(qJ(2));
t293 = t240 * t246;
t213 = t242 * t245 + t249 * t293;
t194 = t213 * t239 + t241 * t292;
t195 = t213 * t241 - t239 * t292;
t170 = t248 * t194 + t195 * t244;
t171 = -t194 * t244 + t195 * t248;
t155 = t247 * t170 + t171 * t243;
t156 = -t170 * t243 + t171 * t247;
t320 = t156 * MDP(31) - t155 * MDP(32);
t299 = pkin(10) + qJ(4);
t223 = t299 * t239;
t224 = t299 * t241;
t197 = -t248 * t223 - t224 * t244;
t198 = -t223 * t244 + t224 * t248;
t184 = -pkin(11) * t220 + t197;
t185 = -pkin(11) * t219 + t198;
t189 = t247 * t219 + t220 * t243;
t190 = -t219 * t243 + t220 * t247;
t267 = t190 * MDP(31) - t189 * MDP(32) + (t184 * t247 - t185 * t243) * MDP(34) - (t184 * t243 + t185 * t247) * MDP(35);
t319 = t220 * MDP(24) - t219 * MDP(25) + t197 * MDP(27) - t198 * MDP(28) + t267;
t182 = t247 * t207 - t208 * t243;
t183 = -t207 * t243 - t208 * t247;
t289 = t183 * MDP(31) - t182 * MDP(32);
t318 = -t208 * MDP(24) - t207 * MDP(25) + t289;
t227 = pkin(8) * t293;
t303 = pkin(1) * t250;
t204 = t227 + (-pkin(2) - t303) * t242;
t212 = -t242 * t249 + t245 * t293;
t175 = pkin(3) * t212 - qJ(4) * t213 + t204;
t269 = pkin(8) * t292;
t304 = pkin(1) * t246;
t205 = t269 + (pkin(9) + t304) * t242;
t206 = (-pkin(2) * t250 - pkin(9) * t246 - pkin(1)) * t240;
t179 = t205 * t249 + t206 * t245;
t176 = -qJ(4) * t292 + t179;
t157 = t241 * t175 - t176 * t239;
t150 = pkin(4) * t212 - pkin(10) * t195 + t157;
t158 = t239 * t175 + t241 * t176;
t152 = -pkin(10) * t194 + t158;
t145 = t248 * t150 - t152 * t244;
t146 = t150 * t244 + t152 * t248;
t316 = t145 * MDP(27) - t146 * MDP(28);
t315 = MDP(21) * qJ(4);
t264 = t239 * MDP(18) + t241 * MDP(19);
t314 = -qJ(4) * t264 - MDP(14) + t319;
t313 = 0.2e1 * MDP(18);
t312 = 0.2e1 * MDP(19);
t311 = 2 * MDP(20);
t310 = -2 * MDP(23);
t309 = 0.2e1 * MDP(27);
t308 = 0.2e1 * MDP(28);
t307 = -2 * MDP(30);
t306 = 0.2e1 * MDP(34);
t305 = 0.2e1 * MDP(35);
t302 = pkin(5) * t212;
t298 = pkin(2) * MDP(17);
t297 = pkin(3) * MDP(21);
t296 = pkin(9) * MDP(17);
t144 = -pkin(11) * t170 + t146;
t295 = t144 * t247;
t291 = t242 * MDP(8);
t290 = t246 * MDP(6);
t221 = pkin(4) * t294 + t245 * pkin(9);
t287 = MDP(15) * t250;
t286 = MDP(18) * t241;
t285 = MDP(19) * t239;
t284 = MDP(27) * t219;
t283 = MDP(34) * t189;
t278 = t156 * MDP(29);
t275 = t171 * MDP(22);
t178 = -t245 * t205 + t206 * t249;
t177 = pkin(3) * t292 - t178;
t274 = t177 * MDP(21);
t273 = t183 * MDP(29);
t272 = t208 * MDP(22);
t271 = t213 * MDP(13);
t270 = MDP(26) + MDP(33);
t268 = t212 * MDP(33) + t320;
t232 = -pkin(4) * t241 - pkin(3);
t143 = -pkin(11) * t171 + t145 + t302;
t140 = t247 * t143 - t144 * t243;
t266 = -t157 * t239 + t158 * t241;
t263 = t171 * MDP(24) - t170 * MDP(25);
t141 = t143 * t243 + t295;
t261 = -t140 * MDP(34) + t141 * MDP(35);
t166 = pkin(4) * t194 + t177;
t257 = t285 - t286 - t297;
t256 = t194 * MDP(18) + t195 * MDP(19) + t274;
t254 = pkin(2) * MDP(16) - t318;
t253 = t213 * MDP(12) - t263 - t320;
t235 = t240 ^ 2;
t215 = t242 * t304 + t269;
t214 = t242 * t303 - t227;
t203 = pkin(5) * t219 + t232;
t200 = -t239 * t300 + t217;
t193 = pkin(5) * t207 + t221;
t151 = pkin(5) * t170 + t166;
t1 = [t235 * t246 ^ 2 * MDP(4) + MDP(1) + (t157 ^ 2 + t158 ^ 2 + t177 ^ 2) * MDP(21) + t213 ^ 2 * MDP(11) + (0.2e1 * t240 * t290 + t291) * t242 + t270 * t212 ^ 2 + (t170 * t310 + t275) * t171 + (t155 * t307 + t278) * t156 + (0.2e1 * MDP(5) * t246 + t287) * t235 * t250 + 0.2e1 * (-t178 * t292 + t204 * t212) * MDP(16) + 0.2e1 * (t179 * t292 + t204 * t213) * MDP(17) + 0.2e1 * (t214 * t242 + t235 * t303) * MDP(9) + 0.2e1 * (-t215 * t242 - t235 * t304) * MDP(10) + (t157 * t212 + t177 * t194) * t313 + (t145 * t212 + t166 * t170) * t309 + (t140 * t212 + t151 * t155) * t306 + (-t158 * t212 + t177 * t195) * t312 + (-t146 * t212 + t166 * t171) * t308 + (-t141 * t212 + t151 * t156) * t305 + (-t157 * t195 - t158 * t194) * t311 + 0.2e1 * (MDP(7) * t242 - t271) * t292 + 0.2e1 * (MDP(14) * t292 - t253) * t212; -t171 * t272 + t156 * t273 - t213 * t298 + (t157 * t200 + t158 * t201) * MDP(21) + (-t155 * t183 - t156 * t182) * MDP(30) + (t170 * t208 - t171 * t207) * MDP(23) + t214 * MDP(9) - t215 * MDP(10) + t291 + (t151 * t182 + t155 * t193) * MDP(34) + (t151 * t183 + t156 * t193) * MDP(35) + (t166 * t207 + t170 * t221) * MDP(27) + (-t166 * t208 + t171 * t221) * MDP(28) + (-t194 * t201 - t195 * t200) * MDP(20) + (t250 * MDP(7) + t290) * t240 + (t200 * MDP(18) - t201 * MDP(19) - t254 + t322) * t212 + (-t204 * MDP(16) - t157 * MDP(18) + t158 * MDP(19) + (-MDP(14) + t296) * t292 - t270 * t212 + t253 + t261 - t316) * t249 + (t213 * MDP(11) - MDP(13) * t292 + (-t157 * t241 - t158 * t239) * MDP(20) - t212 * MDP(12) + t204 * MDP(17) + t264 * t177 + (MDP(16) * t292 + t256) * pkin(9)) * t245; MDP(8) + (t200 ^ 2 + t201 ^ 2) * MDP(21) + t207 * t221 * t309 + t182 * t193 * t306 - (t207 * t310 + t221 * t308 - t272) * t208 + (t182 * t307 + t193 * t305 + t273) * t183 + (-t147 * t306 + t148 * t305 - t168 * t309 + t169 * t308 - t200 * t313 + t201 * t312 + t270 * t249 + 0.2e1 * t254) * t249 + (-0.2e1 * t298 + (-t200 * t241 - t201 * t239) * t311 + 0.2e1 * t249 * MDP(12) + (t301 * t313 + MDP(11) + (t241 * t312 + t321) * pkin(9)) * t245) * t245; t271 - t240 * t287 + t178 * MDP(16) - t179 * MDP(17) + (-pkin(3) * t194 - t177 * t241) * MDP(18) + (-pkin(3) * t195 + t177 * t239) * MDP(19) + t266 * MDP(20) - pkin(3) * t274 + t220 * t275 + (-t170 * t220 - t171 * t219) * MDP(23) + (t166 * t219 + t170 * t232) * MDP(27) + (t166 * t220 + t171 * t232) * MDP(28) + t190 * t278 + (-t155 * t190 - t156 * t189) * MDP(30) + (t151 * t189 + t155 * t203) * MDP(34) + (t151 * t190 + t156 * t203) * MDP(35) + ((-t194 * t241 + t195 * t239) * MDP(20) + t266 * MDP(21)) * qJ(4) + t314 * t212; -t220 * t272 + (-t207 * t220 + t208 * t219) * MDP(23) + (t207 * t232 + t219 * t221) * MDP(27) + (-t208 * t232 + t220 * t221) * MDP(28) + t190 * t273 + (-t182 * t190 - t183 * t189) * MDP(30) + (t182 * t203 + t189 * t193) * MDP(34) + (t183 * t203 + t190 * t193) * MDP(35) + (MDP(13) - t264 * pkin(3) + (-MDP(16) + t257) * pkin(9)) * t245 + (-t296 - t314) * t249 + (MDP(20) + t315) * (-t200 * t239 + t201 * t241); 0.2e1 * t232 * t284 + 0.2e1 * t203 * t283 + MDP(15) + (-0.2e1 * t285 + 0.2e1 * t286 + t297) * pkin(3) + (MDP(22) * t220 + t219 * t310 + t232 * t308) * t220 + (MDP(29) * t190 + t189 * t307 + t203 * t305) * t190 + (t311 + t315) * (t239 ^ 2 + t241 ^ 2) * qJ(4); t170 * MDP(27) + t171 * MDP(28) + t155 * MDP(34) + t156 * MDP(35) + t256; t207 * MDP(27) - t208 * MDP(28) + t182 * MDP(34) + t183 * MDP(35) + (t264 + t321) * t245; MDP(28) * t220 + MDP(35) * t190 + t257 + t283 + t284; MDP(21); t212 * MDP(26) + (t247 * t302 + t140) * MDP(34) + (-t295 + (-t143 - t302) * t243) * MDP(35) + t263 + t268 + t316; (-t270 - t258) * t249 + t318 + t322; t319; 0; 0.2e1 * t258 + t270; -t261 + t268; -t249 * MDP(33) + t260 + t289; t267; 0; MDP(33) + t258; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
