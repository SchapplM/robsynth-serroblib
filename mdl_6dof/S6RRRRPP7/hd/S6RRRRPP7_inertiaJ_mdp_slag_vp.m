% Calculate joint inertia matrix for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRRPP7_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:26:35
% EndTime: 2019-03-09 21:26:39
% DurationCPUTime: 1.51s
% Computational Cost: add. (2810->309), mult. (6115->442), div. (0->0), fcn. (6657->10), ass. (0->116)
t204 = sin(pkin(6));
t270 = 0.2e1 * t204;
t208 = sin(qJ(3));
t269 = -0.2e1 * t208;
t209 = sin(qJ(2));
t249 = t204 * t209;
t189 = pkin(8) * t249;
t206 = cos(pkin(6));
t212 = cos(qJ(2));
t259 = pkin(1) * t212;
t171 = t189 + (-pkin(2) - t259) * t206;
t211 = cos(qJ(3));
t177 = -t206 * t211 + t208 * t249;
t178 = t206 * t208 + t211 * t249;
t147 = t177 * pkin(3) - t178 * pkin(10) + t171;
t248 = t204 * t212;
t229 = pkin(8) * t248;
t260 = pkin(1) * t209;
t172 = t229 + (pkin(9) + t260) * t206;
t173 = (-pkin(2) * t212 - pkin(9) * t209 - pkin(1)) * t204;
t152 = t172 * t211 + t173 * t208;
t149 = -pkin(10) * t248 + t152;
t207 = sin(qJ(4));
t210 = cos(qJ(4));
t136 = t210 * t147 - t149 * t207;
t162 = t178 * t210 - t207 * t248;
t133 = t177 * pkin(4) - t162 * qJ(5) + t136;
t137 = t147 * t207 + t149 * t210;
t161 = t178 * t207 + t210 * t248;
t135 = -t161 * qJ(5) + t137;
t203 = sin(pkin(11));
t205 = cos(pkin(11));
t129 = t203 * t133 + t205 * t135;
t126 = t177 * qJ(6) + t129;
t238 = t162 * MDP(20);
t240 = t161 * MDP(21);
t268 = t136 * MDP(23) - t137 * MDP(24) + t126 * MDP(29) + t238 - t240;
t254 = -qJ(5) - pkin(10);
t187 = t254 * t210;
t226 = t254 * t207;
t165 = -t187 * t203 - t205 * t226;
t167 = -t205 * t187 + t203 * t226;
t214 = t207 * MDP(20) + t210 * MDP(21) - t165 * MDP(27) + t167 * MDP(29) - pkin(10) * (t207 * MDP(23) + t210 * MDP(24));
t195 = pkin(4) * t205 + pkin(5);
t267 = -(pkin(5) + t195) * MDP(27) - MDP(22);
t266 = 0.2e1 * MDP(23);
t265 = 0.2e1 * MDP(24);
t264 = 2 * MDP(25);
t263 = 0.2e1 * MDP(27);
t262 = 2 * MDP(28);
t261 = 0.2e1 * MDP(29);
t258 = pkin(9) * t207;
t257 = pkin(9) * t210;
t256 = pkin(9) * t211;
t253 = MDP(17) * pkin(2);
t252 = pkin(2) * MDP(16);
t151 = -t208 * t172 + t173 * t211;
t148 = pkin(3) * t248 - t151;
t251 = t148 * t207;
t250 = t148 * t210;
t247 = t206 * MDP(8);
t246 = t207 * t208;
t245 = t208 * t210;
t244 = t209 * MDP(6);
t186 = -pkin(3) * t211 - pkin(10) * t208 - pkin(2);
t183 = t210 * t186;
t159 = -qJ(5) * t245 + t183 + (-pkin(4) - t258) * t211;
t228 = t210 * t256;
t163 = t228 + (-qJ(5) * t208 + t186) * t207;
t143 = t203 * t159 + t205 * t163;
t185 = pkin(4) * t246 + t208 * pkin(9);
t243 = MDP(15) * t212;
t181 = t203 * t207 - t205 * t210;
t182 = t203 * t210 + t205 * t207;
t197 = -pkin(4) * t210 - pkin(3);
t155 = t181 * pkin(5) - t182 * qJ(6) + t197;
t241 = t155 * MDP(30);
t239 = t162 * MDP(18);
t237 = t177 * MDP(22);
t236 = t178 * MDP(12);
t235 = t178 * MDP(13);
t234 = t181 * MDP(27);
t233 = t182 * MDP(29);
t192 = pkin(4) * t203 + qJ(6);
t232 = t192 * MDP(29);
t231 = t210 * MDP(18);
t230 = t165 ^ 2 + t167 ^ 2;
t227 = t207 * t210 * MDP(19);
t225 = pkin(9) * MDP(16) - MDP(13);
t224 = pkin(9) * MDP(17) - MDP(14);
t128 = t205 * t133 - t135 * t203;
t144 = t205 * t161 + t162 * t203;
t145 = -t161 * t203 + t162 * t205;
t223 = -t167 * t144 + t145 * t165;
t174 = t182 * t208;
t175 = -t203 * t246 + t205 * t245;
t222 = t165 * t175 - t167 * t174;
t142 = t205 * t159 - t163 * t203;
t220 = t210 * MDP(20) - t207 * MDP(21);
t169 = -t207 * t256 + t183;
t170 = t186 * t207 + t228;
t219 = t169 * MDP(23) - t170 * MDP(24);
t217 = -MDP(12) + t220;
t139 = t161 * pkin(4) + t148;
t202 = t210 ^ 2;
t201 = t208 ^ 2;
t200 = t207 ^ 2;
t199 = t204 ^ 2;
t180 = t206 * t260 + t229;
t179 = t206 * t259 - t189;
t150 = t174 * pkin(5) - t175 * qJ(6) + t185;
t141 = pkin(5) * t211 - t142;
t140 = -qJ(6) * t211 + t143;
t130 = t144 * pkin(5) - t145 * qJ(6) + t139;
t127 = -t177 * pkin(5) - t128;
t1 = [t199 * t209 ^ 2 * MDP(4) + MDP(1) + (t126 ^ 2 + t127 ^ 2 + t130 ^ 2) * MDP(30) + (t128 ^ 2 + t129 ^ 2 + t139 ^ 2) * MDP(26) + t178 ^ 2 * MDP(11) + (t244 * t270 + t247) * t206 + (-0.2e1 * t161 * MDP(19) + t239) * t162 + ((MDP(7) * t206 - t235) * t270 + (0.2e1 * MDP(5) * t209 + t243) * t199) * t212 + (0.2e1 * MDP(14) * t248 - 0.2e1 * t236 + t237 + 0.2e1 * t238 - 0.2e1 * t240) * t177 + (t136 * t177 + t148 * t161) * t266 + (t126 * t177 - t130 * t145) * t261 + (-t137 * t177 + t148 * t162) * t265 + (-t127 * t177 + t130 * t144) * t263 + (-t126 * t144 + t127 * t145) * t262 + (-t128 * t145 - t129 * t144) * t264 + 0.2e1 * (-t180 * t206 - t199 * t260) * MDP(10) + 0.2e1 * (t179 * t206 + t199 * t259) * MDP(9) + 0.2e1 * (t152 * t248 + t171 * t178) * MDP(17) + 0.2e1 * (-t151 * t248 + t171 * t177) * MDP(16); t247 + t179 * MDP(9) - t180 * MDP(10) - t178 * t253 + (-t128 * t175 - t129 * t174 - t142 * t145 - t143 * t144) * MDP(25) + (t128 * t142 + t129 * t143 + t139 * t185) * MDP(26) + (t130 * t174 + t150 * t144) * MDP(27) + (-t126 * t174 + t127 * t175 - t140 * t144 + t141 * t145) * MDP(28) + (-t130 * t175 - t150 * t145) * MDP(29) + (t126 * t140 + t127 * t141 + t130 * t150) * MDP(30) + (t212 * MDP(7) + t244) * t204 + (-t141 * MDP(27) + t140 * MDP(29) + t219 - t252) * t177 + (-t171 * MDP(16) + t127 * MDP(27) + t224 * t248 + t236 - t237 - t268) * t211 + (t171 * MDP(17) + (-t161 * t210 - t162 * t207) * MDP(19) + (pkin(9) * t161 + t251) * MDP(23) + (pkin(9) * t162 + t250) * MDP(24) + t162 * t231 + t178 * MDP(11) + t225 * t248 + t217 * t177) * t208; MDP(8) + t253 * t269 + (t142 ^ 2 + t143 ^ 2 + t185 ^ 2) * MDP(26) + (t140 ^ 2 + t141 ^ 2 + t150 ^ 2) * MDP(30) + (MDP(18) * t202 + MDP(11) - 0.2e1 * t227) * t201 + (t211 * MDP(22) + t217 * t269 + 0.2e1 * t252) * t211 + (-t169 * t211 + t201 * t258) * t266 + (t170 * t211 + t201 * t257) * t265 + (-t142 * t175 - t143 * t174) * t264 + (t141 * t211 + t150 * t174) * t263 + (-t140 * t174 + t141 * t175) * t262 + (-t140 * t211 - t150 * t175) * t261; t235 - t204 * t243 + t151 * MDP(16) - t152 * MDP(17) + t207 * t239 + (-t161 * t207 + t162 * t210) * MDP(19) + (-pkin(3) * t161 - t250) * MDP(23) + (-pkin(3) * t162 + t251) * MDP(24) + (-t128 * t182 - t129 * t181 + t223) * MDP(25) + (-t128 * t165 + t129 * t167 + t139 * t197) * MDP(26) + (t130 * t181 + t144 * t155) * MDP(27) + (-t126 * t181 + t127 * t182 + t223) * MDP(28) + (-t130 * t182 - t145 * t155) * MDP(29) + (t126 * t167 + t127 * t165 + t130 * t155) * MDP(30) + (-MDP(14) + t214) * t177; (-t142 * t182 - t143 * t181 + t222) * MDP(25) + (-t142 * t165 + t143 * t167 + t185 * t197) * MDP(26) + (t150 * t181 + t155 * t174) * MDP(27) + (-t140 * t181 + t141 * t182 + t222) * MDP(28) + (-t150 * t182 - t155 * t175) * MDP(29) + (t140 * t167 + t141 * t165 + t150 * t155) * MDP(30) + (-t214 - t224) * t211 + (t207 * t231 + (-t200 + t202) * MDP(19) + (-pkin(3) * t207 - t257) * MDP(23) + (-pkin(3) * t210 + t258) * MDP(24) - t225) * t208; MDP(15) + t200 * MDP(18) + 0.2e1 * t227 + (t197 ^ 2 + t230) * MDP(26) + t230 * MDP(30) + (-0.2e1 * t233 + 0.2e1 * t234 + t241) * t155 + 0.2e1 * (MDP(23) * t210 - MDP(24) * t207) * pkin(3) + (t264 + t262) * (t165 * t182 - t167 * t181); t128 * MDP(27) + (-t144 * t192 - t145 * t195) * MDP(28) + (t126 * t192 - t127 * t195) * MDP(30) + (t232 - t267) * t177 + ((-t144 * t203 - t145 * t205) * MDP(25) + (t128 * t205 + t129 * t203) * MDP(26)) * pkin(4) + t268; t142 * MDP(27) + (-t174 * t192 - t175 * t195) * MDP(28) + t143 * MDP(29) + (t140 * t192 - t141 * t195) * MDP(30) + t220 * t208 + ((-qJ(6) - t192) * MDP(29) + t267) * t211 + ((-t174 * t203 - t175 * t205) * MDP(25) + (t142 * t205 + t143 * t203) * MDP(26)) * pkin(4) + t219; (-t181 * t192 - t182 * t195) * MDP(28) + (-t165 * t195 + t167 * t192) * MDP(30) + ((-t181 * t203 - t182 * t205) * MDP(25) + (-t165 * t205 + t167 * t203) * MDP(26)) * pkin(4) + t214; MDP(22) + (t192 ^ 2 + t195 ^ 2) * MDP(30) + (t203 ^ 2 + t205 ^ 2) * MDP(26) * pkin(4) ^ 2 + t195 * t263 + 0.2e1 * t232; MDP(26) * t139 + t144 * MDP(27) - t145 * MDP(29) + t130 * MDP(30); MDP(26) * t185 + t174 * MDP(27) - t175 * MDP(29) + t150 * MDP(30); MDP(26) * t197 - t233 + t234 + t241; 0; MDP(26) + MDP(30); -t177 * MDP(27) + t145 * MDP(28) + t127 * MDP(30); MDP(27) * t211 + t175 * MDP(28) + t141 * MDP(30); MDP(28) * t182 + MDP(30) * t165; -MDP(30) * t195 - MDP(27); 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
