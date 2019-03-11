% Calculate joint inertia matrix for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP10_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:37:01
% EndTime: 2019-03-09 17:37:06
% DurationCPUTime: 1.69s
% Computational Cost: add. (2730->327), mult. (6033->470), div. (0->0), fcn. (6616->10), ass. (0->114)
t269 = (MDP(21) * qJ(4));
t209 = cos(qJ(3));
t268 = t209 * MDP(12);
t204 = cos(pkin(11));
t249 = pkin(10) + qJ(4);
t186 = t249 * t204;
t206 = sin(qJ(5));
t202 = sin(pkin(11));
t221 = t249 * t202;
t256 = cos(qJ(5));
t164 = t186 * t206 + t256 * t221;
t165 = t256 * t186 - t206 * t221;
t222 = t256 * t204;
t182 = t202 * t206 - t222;
t183 = t256 * t202 + t206 * t204;
t215 = t183 * MDP(24) - t182 * MDP(25);
t216 = t202 * MDP(18) + t204 * MDP(19);
t224 = MDP(28) - MDP(31);
t225 = MDP(27) + MDP(29);
t267 = -t216 * qJ(4) - t225 * t164 - t224 * t165 - MDP(14) + t215;
t266 = 2 * MDP(16);
t265 = 0.2e1 * MDP(18);
t264 = 0.2e1 * MDP(19);
t263 = 2 * MDP(20);
t262 = -2 * MDP(23);
t261 = 2 * MDP(27);
t260 = 2 * MDP(28);
t259 = 2 * MDP(29);
t258 = 2 * MDP(30);
t257 = 2 * MDP(31);
t208 = sin(qJ(2));
t255 = pkin(1) * t208;
t210 = cos(qJ(2));
t254 = pkin(1) * t210;
t205 = cos(pkin(6));
t207 = sin(qJ(3));
t203 = sin(pkin(6));
t242 = t203 * t208;
t175 = -t205 * t209 + t207 * t242;
t253 = pkin(5) * t175;
t252 = pkin(5) * t209;
t251 = pkin(9) * t202;
t250 = pkin(9) * t209;
t248 = pkin(2) * MDP(17);
t247 = pkin(3) * MDP(21);
t246 = pkin(9) * MDP(17);
t245 = qJ(6) * t175;
t244 = qJ(6) * t209;
t243 = t202 * t207;
t241 = t203 * t210;
t240 = t205 * MDP(8);
t191 = pkin(8) * t242;
t169 = t191 + (-pkin(2) - t254) * t205;
t176 = t205 * t207 + t209 * t242;
t149 = pkin(3) * t175 - qJ(4) * t176 + t169;
t223 = pkin(8) * t241;
t170 = t223 + (pkin(9) + t255) * t205;
t171 = (-pkin(2) * t210 - pkin(9) * t208 - pkin(1)) * t203;
t154 = t170 * t209 + t171 * t207;
t150 = -qJ(4) * t241 + t154;
t137 = t204 * t149 - t150 * t202;
t162 = t176 * t204 - t202 * t241;
t134 = pkin(4) * t175 - pkin(10) * t162 + t137;
t138 = t202 * t149 + t204 * t150;
t161 = t176 * t202 + t204 * t241;
t136 = -pkin(10) * t161 + t138;
t130 = t206 * t134 + t256 * t136;
t185 = -pkin(3) * t209 - qJ(4) * t207 - pkin(2);
t180 = t204 * t185;
t158 = -pkin(10) * t204 * t207 + t180 + (-pkin(4) - t251) * t209;
t168 = t202 * t185 + t204 * t250;
t163 = -pkin(10) * t243 + t168;
t143 = t206 * t158 + t256 * t163;
t184 = pkin(4) * t243 + t207 * pkin(9);
t238 = MDP(15) * t210;
t237 = MDP(18) * t204;
t236 = MDP(19) * t202;
t144 = t256 * t161 + t162 * t206;
t235 = t144 * MDP(25);
t145 = -t206 * t161 + t256 * t162;
t234 = t145 * MDP(22);
t233 = t145 * MDP(24);
t153 = -t207 * t170 + t171 * t209;
t151 = pkin(3) * t241 - t153;
t232 = t151 * MDP(21);
t173 = t183 * t207;
t231 = t173 * MDP(25);
t174 = -t206 * t243 + t207 * t222;
t230 = t174 * MDP(22);
t229 = t174 * MDP(24);
t228 = t175 * MDP(26);
t227 = t176 * MDP(13);
t226 = t209 * MDP(26);
t196 = -pkin(4) * t204 - pkin(3);
t220 = -MDP(32) * pkin(5) - MDP(29);
t219 = -t256 * t134 + t136 * t206;
t142 = t256 * t158 - t163 * t206;
t218 = -t137 * t202 + t138 * t204;
t139 = pkin(4) * t161 + t151;
t213 = t236 - t237 - t247;
t212 = t161 * MDP(18) + t162 * MDP(19) + t232;
t201 = t207 ^ 2;
t199 = t203 ^ 2;
t178 = t205 * t255 + t223;
t177 = t205 * t254 - t191;
t167 = -t202 * t250 + t180;
t157 = pkin(5) * t182 - qJ(6) * t183 + t196;
t152 = pkin(5) * t173 - qJ(6) * t174 + t184;
t141 = -t142 + t252;
t140 = t143 - t244;
t131 = pkin(5) * t144 - qJ(6) * t145 + t139;
t128 = t219 - t253;
t127 = t130 + t245;
t1 = [(t127 ^ 2 + t128 ^ 2 + t131 ^ 2) * MDP(32) + (t137 ^ 2 + t138 ^ 2 + t151 ^ 2) * MDP(21) + t176 ^ 2 * MDP(11) + t199 * t208 ^ 2 * MDP(4) + MDP(1) + (0.2e1 * MDP(6) * t242 + t240) * t205 + (t144 * t262 + t234) * t145 + (0.2e1 * (MDP(7) * t205 - t227) * t203 + (0.2e1 * MDP(5) * t208 + t238) * t199) * t210 + (-0.2e1 * t176 * MDP(12) + 0.2e1 * MDP(14) * t241 + t228 + 0.2e1 * t233 - 0.2e1 * t235) * t175 + 0.2e1 * (t177 * t205 + t199 * t254) * MDP(9) + 0.2e1 * (t154 * t241 + t169 * t176) * MDP(17) + (-t153 * t241 + t169 * t175) * t266 + 0.2e1 * (-t178 * t205 - t199 * t255) * MDP(10) + (-t128 * t175 + t131 * t144) * t259 + (-t130 * t175 + t139 * t145) * t260 + (t137 * t175 + t151 * t161) * t265 + (t139 * t144 - t175 * t219) * t261 + (t127 * t175 - t131 * t145) * t257 + (-t138 * t175 + t151 * t162) * t264 + (-t137 * t162 - t138 * t161) * t263 + (-t127 * t144 + t128 * t145) * t258; (t127 * t140 + t128 * t141 + t131 * t152) * MDP(32) + (-t127 * t173 + t128 * t174 - t140 * t144 + t141 * t145) * MDP(30) + (-t144 * t174 - t145 * t173) * MDP(23) + t177 * MDP(9) - t178 * MDP(10) + t240 + (t137 * t167 + t138 * t168) * MDP(21) + (-t161 * t168 - t162 * t167) * MDP(20) + (t144 * t209 - t173 * t175) * MDP(25) + (-t145 * t209 + t174 * t175) * MDP(24) + (-t127 * t209 - t131 * t174 + t140 * t175 - t152 * t145) * MDP(31) + (t128 * t209 + t131 * t173 - t141 * t175 + t144 * t152) * MDP(29) + (t139 * t173 + t142 * t175 + t144 * t184 + t209 * t219) * MDP(27) + (t130 * t209 + t139 * t174 - t143 * t175 + t145 * t184) * MDP(28) + (-t137 * t209 + t167 * t175) * MDP(18) + (t138 * t209 - t168 * t175) * MDP(19) + (-pkin(2) * t175 - t169 * t209) * MDP(16) - t175 * t226 + t145 * t230 + (-t248 + t268) * t176 + (MDP(6) * t208 + (MDP(7) + (-MDP(14) + t246) * t209) * t210) * t203 + (-MDP(13) * t241 + (-t137 * t204 - t138 * t202) * MDP(20) - t175 * MDP(12) + t169 * MDP(17) + t176 * MDP(11) + t216 * t151 + (MDP(16) * t241 + t212) * pkin(9)) * t207; MDP(8) + t201 * MDP(11) + (pkin(9) ^ 2 * t201 + t167 ^ 2 + t168 ^ 2) * MDP(21) + (t140 ^ 2 + t141 ^ 2 + t152 ^ 2) * MDP(32) + (t173 * t262 + t230) * t174 + (pkin(2) * t266 + t226 - 0.2e1 * t229 + 0.2e1 * t231) * t209 + (-t167 * t209 + t201 * t251) * t265 + (pkin(9) * t201 * t204 + t168 * t209) * t264 + (-t142 * t209 + t173 * t184) * t261 + (t143 * t209 + t174 * t184) * t260 + (t141 * t209 + t152 * t173) * t259 + (-t140 * t173 + t141 * t174) * t258 + (-t140 * t209 - t152 * t174) * t257 + (-0.2e1 * t248 + 0.2e1 * t268 + (-t167 * t204 - t168 * t202) * t263) * t207; t227 - t203 * t238 + t153 * MDP(16) - t154 * MDP(17) + (-pkin(3) * t161 - t151 * t204) * MDP(18) + (-pkin(3) * t162 + t151 * t202) * MDP(19) + t218 * MDP(20) - pkin(3) * t232 + t183 * t234 + (-t144 * t183 - t145 * t182) * MDP(23) + (t139 * t182 + t144 * t196) * MDP(27) + (t139 * t183 + t145 * t196) * MDP(28) + (t131 * t182 + t144 * t157) * MDP(29) + (-t127 * t182 + t128 * t183 - t144 * t165 + t145 * t164) * MDP(30) + (-t131 * t183 - t145 * t157) * MDP(31) + (t127 * t165 + t128 * t164 + t131 * t157) * MDP(32) + ((-t161 * t204 + t162 * t202) * MDP(20) + t218 * MDP(21)) * qJ(4) + t267 * t175; t183 * t230 + (-t173 * t183 - t174 * t182) * MDP(23) + (t173 * t196 + t182 * t184) * MDP(27) + (t174 * t196 + t183 * t184) * MDP(28) + (t152 * t182 + t157 * t173) * MDP(29) + (-t140 * t182 + t141 * t183 + t164 * t174 - t165 * t173) * MDP(30) + (-t152 * t183 - t157 * t174) * MDP(31) + (t140 * t165 + t141 * t164 + t152 * t157) * MDP(32) + (-t246 - t267) * t209 + (MDP(13) - t216 * pkin(3) + (-MDP(16) + t213) * pkin(9)) * t207 + (MDP(20) + t269) * (-t167 * t202 + t168 * t204); MDP(15) + (t157 ^ 2 + t164 ^ 2 + t165 ^ 2) * MDP(32) + (-0.2e1 * t236 + 0.2e1 * t237 + t247) * pkin(3) + 0.2e1 * (MDP(27) * t196 + MDP(29) * t157 - MDP(30) * t165) * t182 + (MDP(22) * t183 - 0.2e1 * MDP(31) * t157 + t164 * t258 + t182 * t262 + t196 * t260) * t183 + (t263 + t269) * (t202 ^ 2 + t204 ^ 2) * qJ(4); t131 * MDP(32) + t225 * t144 + t224 * t145 + t212; t152 * MDP(32) + t224 * t174 + t225 * t173 + (pkin(9) * MDP(21) + t216) * t207; MDP(32) * t157 + t225 * t182 + t224 * t183 + t213; MDP(21) + MDP(32); t233 - t235 + t228 - t219 * MDP(27) - t130 * MDP(28) + (-t219 + 0.2e1 * t253) * MDP(29) + (-pkin(5) * t145 - qJ(6) * t144) * MDP(30) + (t130 + 0.2e1 * t245) * MDP(31) + (-pkin(5) * t128 + qJ(6) * t127) * MDP(32); t229 - t231 - t226 + t142 * MDP(27) - t143 * MDP(28) + (t142 - 0.2e1 * t252) * MDP(29) + (-pkin(5) * t174 - qJ(6) * t173) * MDP(30) + (t143 - 0.2e1 * t244) * MDP(31) + (-pkin(5) * t141 + qJ(6) * t140) * MDP(32); (-pkin(5) * t183 - qJ(6) * t182) * MDP(30) + (MDP(32) * qJ(6) - t224) * t165 + (-MDP(27) + t220) * t164 + t215; 0; MDP(26) + pkin(5) * t259 + qJ(6) * t257 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(32); -t175 * MDP(29) + t145 * MDP(30) + t128 * MDP(32); t209 * MDP(29) + t174 * MDP(30) + t141 * MDP(32); MDP(30) * t183 + MDP(32) * t164; 0; t220; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
