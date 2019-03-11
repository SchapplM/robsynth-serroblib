% Calculate joint inertia matrix for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR13_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR13_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR13_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:53:23
% EndTime: 2019-03-09 14:53:28
% DurationCPUTime: 1.55s
% Computational Cost: add. (1358->262), mult. (3003->361), div. (0->0), fcn. (3225->10), ass. (0->119)
t205 = sin(pkin(6));
t214 = cos(qJ(2));
t258 = t205 * t214;
t210 = sin(qJ(2));
t259 = t205 * t210;
t191 = pkin(8) * t259;
t206 = cos(pkin(6));
t269 = pkin(1) * t214;
t235 = -pkin(2) - t269;
t151 = pkin(3) * t259 + t191 + (-pkin(9) + t235) * t206;
t215 = -pkin(2) - pkin(9);
t233 = -qJ(3) * t210 - pkin(1);
t159 = (t214 * t215 + t233) * t205;
t209 = sin(qJ(4));
t213 = cos(qJ(4));
t142 = t151 * t213 - t209 * t159;
t143 = t151 * t209 + t159 * t213;
t178 = t206 * t213 - t209 * t258;
t284 = t178 * MDP(17) + t142 * MDP(20) - t143 * MDP(21);
t283 = 0.2e1 * t206;
t208 = sin(qJ(5));
t212 = cos(qJ(5));
t154 = t178 * t208 - t212 * t259;
t155 = t178 * t212 + t208 * t259;
t207 = sin(qJ(6));
t211 = cos(qJ(6));
t145 = t211 * t154 + t155 * t207;
t146 = -t154 * t207 + t155 * t211;
t282 = t146 * MDP(31) - t145 * MDP(32);
t186 = t207 * t212 + t208 * t211;
t173 = t186 * t213;
t185 = t207 * t208 - t211 * t212;
t175 = t185 * t213;
t281 = -t175 * MDP(31) - t173 * MDP(32);
t270 = pkin(10) + pkin(11);
t188 = t270 * t208;
t189 = t270 * t212;
t231 = t186 * MDP(31) - t185 * MDP(32) + (-t188 * t211 - t189 * t207) * MDP(34) - (-t188 * t207 + t189 * t211) * MDP(35);
t278 = MDP(27) * t208 + MDP(28) * t212;
t280 = t208 * MDP(24) + t212 * MDP(25) - pkin(10) * t278 + t231;
t137 = pkin(10) * t259 + t143;
t180 = t206 * t210 * pkin(1) + pkin(8) * t258;
t197 = t206 * qJ(3);
t164 = -t197 - t180;
t158 = pkin(3) * t258 - t164;
t177 = t206 * t209 + t213 * t258;
t147 = pkin(4) * t177 - pkin(10) * t178 + t158;
t133 = -t137 * t208 + t212 * t147;
t134 = t137 * t212 + t147 * t208;
t279 = t133 * MDP(27) - t134 * MDP(28);
t277 = qJ(3) * MDP(20) + t281;
t275 = 0.2e1 * MDP(27);
t274 = 0.2e1 * MDP(28);
t273 = -2 * MDP(30);
t272 = 0.2e1 * MDP(34);
t271 = 0.2e1 * MDP(35);
t268 = pkin(5) * t177;
t267 = pkin(11) * t213;
t266 = t209 * pkin(5);
t265 = pkin(2) * MDP(14);
t132 = -pkin(11) * t154 + t134;
t264 = t132 * t211;
t136 = -pkin(4) * t259 - t142;
t263 = t136 * t208;
t262 = t136 * t212;
t187 = pkin(4) * t209 - pkin(10) * t213 + qJ(3);
t255 = t212 * t215;
t238 = t209 * t255;
t153 = t238 + (t187 - t267) * t208;
t261 = t153 * t211;
t260 = t177 * t209;
t257 = t208 * t212;
t256 = t208 * t215;
t172 = t186 * t209;
t174 = t185 * t209;
t254 = -t172 * MDP(34) + t174 * MDP(35);
t252 = MDP(21) * t178;
t249 = MDP(29) * t175;
t248 = MDP(29) * t186;
t247 = MDP(34) * t185;
t243 = t155 * MDP(22);
t241 = t209 * MDP(21);
t240 = t215 * MDP(21);
t239 = MDP(26) + MDP(33);
t237 = t177 * MDP(33) + t282;
t236 = t209 * MDP(33) + t281;
t234 = MDP(23) * t257;
t131 = -pkin(11) * t155 + t133 + t268;
t128 = t211 * t131 - t132 * t207;
t183 = t212 * t187;
t149 = -t212 * t267 + t183 + (pkin(5) - t256) * t209;
t138 = t211 * t149 - t153 * t207;
t232 = t215 * MDP(20) + MDP(17);
t230 = t210 * MDP(6) + t214 * MDP(7);
t229 = -t180 * MDP(10) + (t206 * t269 - t191) * MDP(9);
t227 = t155 * MDP(24) - t154 * MDP(25);
t226 = MDP(24) * t212 - MDP(25) * t208;
t162 = -t209 * t256 + t183;
t163 = t208 * t187 + t238;
t225 = t162 * MDP(27) - t163 * MDP(28);
t224 = MDP(27) * t212 - MDP(28) * t208;
t129 = t131 * t207 + t264;
t222 = t128 * MDP(34) - t129 * MDP(35);
t139 = t149 * t207 + t261;
t221 = MDP(34) * t138 - MDP(35) * t139;
t220 = -MDP(16) + t226;
t219 = (MDP(34) * t211 - MDP(35) * t207) * pkin(5);
t217 = -t178 * MDP(16) + t227 + t282;
t216 = -MDP(18) + t280;
t204 = t213 ^ 2;
t203 = t212 ^ 2;
t201 = t209 ^ 2;
t200 = t208 ^ 2;
t196 = -pkin(5) * t212 - pkin(4);
t184 = (pkin(5) * t208 - t215) * t213;
t170 = t206 * t235 + t191;
t165 = (-pkin(2) * t214 + t233) * t205;
t135 = pkin(5) * t154 + t136;
t1 = [t206 ^ 2 * MDP(8) + t178 ^ 2 * MDP(15) + (t164 ^ 2 + t165 ^ 2 + t170 ^ 2) * MDP(14) + MDP(1) + t239 * t177 ^ 2 + (-0.2e1 * t154 * MDP(23) + t243) * t155 + (MDP(29) * t146 + t145 * t273) * t146 + (-t134 * t177 + t136 * t155) * t274 + (-t129 * t177 + t135 * t146) * t271 + (t133 * t177 + t136 * t154) * t275 + (t128 * t177 + t135 * t145) * t272 + (MDP(12) * t170 - MDP(13) * t164 + t229) * t283 + 0.2e1 * (MDP(20) * t177 + t252) * t158 + 0.2e1 * (-MDP(18) * t259 + t217) * t177 + (0.2e1 * (-MDP(11) * t164 + MDP(12) * t165) * t214 + t230 * t283 + ((MDP(19) + MDP(4)) * t210 ^ 2 + 0.2e1 * (-MDP(10) * t210 + MDP(9) * t214) * pkin(1)) * t205 + 0.2e1 * (MDP(11) * t170 - MDP(13) * t165 + MDP(5) * t258 + t284) * t210) * t205; t191 * MDP(12) + (-t135 * t175 + t146 * t184) * MDP(35) + (t135 * t173 + t145 * t184) * MDP(34) + (0.2e1 * t197 + t180) * MDP(13) + (t145 * t175 - t146 * t173) * MDP(30) + (-pkin(2) * t170 - qJ(3) * t164) * MDP(14) - t146 * t249 + qJ(3) * t252 + ((-0.2e1 * pkin(2) - t269) * MDP(12) + MDP(8)) * t206 + ((-pkin(2) * t210 + qJ(3) * t214) * MDP(11) + t230) * t205 + (t221 + t225 + t277) * t177 + (t158 * MDP(20) + (-MDP(18) - t240) * t259 + t239 * t177 + t217 + t222 + t279) * t209 + ((-t154 * t212 - t155 * t208) * MDP(23) + (-t155 * t215 + t262) * MDP(28) + (-t154 * t215 + t263) * MDP(27) + t158 * MDP(21) + t178 * MDP(15) + t212 * t243 + t232 * t259 + t220 * t177) * t213 + t229; MDP(8) + t239 * t201 - (t173 * t273 - t249) * t175 + (-0.2e1 * MDP(12) + t265) * pkin(2) + (MDP(14) * qJ(3) + 0.2e1 * MDP(21) * t213 + 0.2e1 * MDP(13)) * qJ(3) + (MDP(22) * t203 + MDP(15) - 0.2e1 * t234) * t204 + 0.2e1 * (t213 * t220 + t277) * t209 + (t162 * t209 - t204 * t256) * t275 + (-t163 * t209 - t204 * t255) * t274 + (t138 * t209 + t173 * t184) * t272 + (-t139 * t209 - t175 * t184) * t271; t206 * MDP(12) + t170 * MDP(14) + (-t154 * t213 - t208 * t260) * MDP(27) + (-t155 * t213 - t212 * t260) * MDP(28) + (-t145 * t213 - t172 * t177) * MDP(34) + (-t146 * t213 + t174 * t177) * MDP(35) + (MDP(20) * t213 + MDP(11) - t241) * t259; MDP(12) - t265 + (-t172 * t209 - t173 * t213) * MDP(34) + (t174 * t209 + t175 * t213) * MDP(35) + t278 * (-t201 - t204); MDP(14); MDP(19) * t259 + t208 * t243 + (-t154 * t208 + t155 * t212) * MDP(23) + (-pkin(4) * t154 - t262) * MDP(27) + (-pkin(4) * t155 + t263) * MDP(28) + t146 * t248 + (-t145 * t186 - t146 * t185) * MDP(30) + (t135 * t185 + t145 * t196) * MDP(34) + (t135 * t186 + t146 * t196) * MDP(35) + t216 * t177 + t284; -t175 * t248 + (-t173 * t186 + t175 * t185) * MDP(30) + (t173 * t196 + t184 * t185) * MDP(34) + (-t175 * t196 + t184 * t186) * MDP(35) + (t216 - t240) * t209 + (MDP(22) * t257 + (-t200 + t203) * MDP(23) + (-pkin(4) * t208 + t255) * MDP(27) + (-pkin(4) * t212 - t256) * MDP(28) + t232) * t213; -t241 + (-MDP(35) * t186 + MDP(20) + t224 - t247) * t213; 0.2e1 * t234 + 0.2e1 * t196 * t247 + MDP(22) * t200 + MDP(19) + 0.2e1 * t224 * pkin(4) + (t185 * t273 + t196 * t271 + t248) * t186; t177 * MDP(26) + (t211 * t268 + t128) * MDP(34) + (-t264 + (-t131 - t268) * t207) * MDP(35) + t227 + t237 + t279; t209 * MDP(26) + (t211 * t266 + t138) * MDP(34) + (-t261 + (-t149 - t266) * t207) * MDP(35) + t226 * t213 + t225 + t236; -t209 * t278 + t254; t280; 0.2e1 * t219 + t239; t222 + t237; t221 + t236; t254; t231; MDP(33) + t219; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
