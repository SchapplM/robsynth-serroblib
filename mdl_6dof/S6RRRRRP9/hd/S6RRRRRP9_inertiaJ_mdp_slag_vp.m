% Calculate joint inertia matrix for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP9_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:11:43
% EndTime: 2019-03-10 02:11:47
% DurationCPUTime: 1.63s
% Computational Cost: add. (2180->288), mult. (4876->409), div. (0->0), fcn. (5373->10), ass. (0->125)
t205 = sin(qJ(4));
t208 = cos(qJ(4));
t263 = pkin(10) + pkin(11);
t186 = t263 * t205;
t187 = t263 * t208;
t204 = sin(qJ(5));
t262 = cos(qJ(5));
t159 = -t262 * t186 - t187 * t204;
t160 = -t204 * t186 + t187 * t262;
t226 = t262 * t208;
t181 = t204 * t205 - t226;
t182 = t204 * t208 + t205 * t262;
t222 = t182 * MDP(27) - t181 * MDP(28) + t159 * MDP(30) - t160 * MDP(31);
t277 = t205 * MDP(20) + t208 * MDP(21) + t222 - (t205 * MDP(23) + t208 * MDP(24)) * pkin(10);
t276 = -MDP(14) + t277;
t206 = sin(qJ(3));
t274 = -0.2e1 * t206;
t202 = sin(pkin(6));
t210 = cos(qJ(2));
t247 = t202 * t210;
t203 = cos(pkin(6));
t209 = cos(qJ(3));
t207 = sin(qJ(2));
t248 = t202 * t207;
t175 = t203 * t206 + t209 * t248;
t154 = t175 * t205 + t208 * t247;
t155 = t175 * t208 - t205 * t247;
t140 = t154 * t262 + t155 * t204;
t141 = -t204 * t154 + t155 * t262;
t273 = t141 * MDP(27) - t140 * MDP(28);
t171 = t182 * t206;
t165 = t171 * MDP(28);
t245 = t205 * t206;
t172 = -t204 * t245 + t206 * t226;
t166 = t172 * MDP(27);
t272 = t166 - t165;
t271 = t208 * MDP(20) - t205 * MDP(21);
t189 = pkin(8) * t248;
t260 = pkin(1) * t210;
t167 = t189 + (-pkin(2) - t260) * t203;
t174 = -t203 * t209 + t206 * t248;
t143 = t174 * pkin(3) - t175 * pkin(10) + t167;
t232 = pkin(8) * t247;
t261 = pkin(1) * t207;
t168 = t232 + (pkin(9) + t261) * t203;
t169 = (-pkin(2) * t210 - pkin(9) * t207 - pkin(1)) * t202;
t147 = t209 * t168 + t206 * t169;
t145 = -pkin(10) * t247 + t147;
t131 = t208 * t143 - t145 * t205;
t132 = t143 * t205 + t145 * t208;
t270 = -t131 * MDP(23) + t132 * MDP(24);
t269 = 0.2e1 * MDP(23);
t268 = 0.2e1 * MDP(24);
t267 = -2 * MDP(26);
t266 = 0.2e1 * MDP(30);
t265 = 0.2e1 * MDP(31);
t264 = 2 * MDP(32);
t259 = pkin(4) * t174;
t258 = pkin(4) * t204;
t257 = pkin(9) * t205;
t256 = pkin(9) * t208;
t255 = pkin(9) * t209;
t254 = pkin(11) * t206;
t253 = pkin(2) * MDP(16);
t252 = pkin(2) * MDP(17);
t251 = pkin(9) * MDP(17);
t146 = -t206 * t168 + t169 * t209;
t144 = pkin(3) * t247 - t146;
t250 = t144 * t205;
t249 = t144 * t208;
t246 = t203 * MDP(8);
t244 = t207 * MDP(6);
t184 = pkin(4) * t245 + t206 * pkin(9);
t243 = MDP(15) * t210;
t240 = t141 * MDP(25);
t239 = t155 * MDP(18);
t238 = t172 * MDP(25);
t237 = t175 * MDP(13);
t235 = t208 * MDP(18);
t233 = MDP(22) + MDP(29);
t231 = t208 * t255;
t230 = t262 * pkin(4);
t229 = t174 * MDP(29) + t273;
t195 = -pkin(4) * t208 - pkin(3);
t130 = -pkin(11) * t154 + t132;
t228 = t262 * t130;
t185 = -pkin(3) * t209 - pkin(10) * t206 - pkin(2);
t156 = t231 + (t185 - t254) * t205;
t227 = t262 * t156;
t225 = t205 * t208 * MDP(19);
t224 = MDP(30) * t262;
t223 = pkin(9) * MDP(16) - MDP(13);
t128 = -pkin(11) * t155 + t131 + t259;
t125 = t262 * t128 - t204 * t130;
t180 = t208 * t185;
t152 = -t208 * t254 + t180 + (-pkin(4) - t257) * t209;
t138 = t262 * t152 - t204 * t156;
t221 = -t209 * MDP(29) + t272;
t220 = t155 * MDP(20) - t154 * MDP(21);
t162 = -t205 * t255 + t180;
t163 = t185 * t205 + t231;
t218 = t162 * MDP(23) - t163 * MDP(24);
t126 = t204 * t128 + t228;
t216 = t125 * MDP(30) - t126 * MDP(31);
t139 = t204 * t152 + t227;
t215 = t138 * MDP(30) - t139 * MDP(31);
t214 = -MDP(12) + t271;
t135 = pkin(4) * t154 + t144;
t212 = t175 * MDP(12) - t220 - t273;
t200 = t208 ^ 2;
t198 = t205 ^ 2;
t197 = t202 ^ 2;
t194 = t230 + pkin(5);
t177 = t203 * t261 + t232;
t176 = t203 * t260 - t189;
t164 = pkin(5) * t181 + t195;
t153 = pkin(5) * t171 + t184;
t149 = -t181 * qJ(6) + t160;
t148 = -qJ(6) * t182 + t159;
t134 = -t171 * qJ(6) + t139;
t133 = -pkin(5) * t209 - qJ(6) * t172 + t138;
t129 = pkin(5) * t140 + t135;
t124 = -t140 * qJ(6) + t126;
t123 = pkin(5) * t174 - qJ(6) * t141 + t125;
t1 = [t175 ^ 2 * MDP(11) + (t123 ^ 2 + t124 ^ 2 + t129 ^ 2) * MDP(33) + t197 * t207 ^ 2 * MDP(4) + MDP(1) + (0.2e1 * t202 * t244 + t246) * t203 + t233 * t174 ^ 2 + (-0.2e1 * t154 * MDP(19) + t239) * t155 + (t140 * t267 + t240) * t141 + (0.2e1 * MDP(5) * t207 + t243) * t197 * t210 + 0.2e1 * (t147 * t247 + t167 * t175) * MDP(17) + 0.2e1 * (-t146 * t247 + t167 * t174) * MDP(16) + 0.2e1 * (t176 * t203 + t197 * t260) * MDP(9) + 0.2e1 * (-t177 * t203 - t197 * t261) * MDP(10) + (-t132 * t174 + t144 * t155) * t268 + (-t126 * t174 + t135 * t141) * t265 + (t131 * t174 + t144 * t154) * t269 + (t125 * t174 + t135 * t140) * t266 + (-t123 * t141 - t124 * t140) * t264 + 0.2e1 * (MDP(7) * t203 - t237) * t247 + 0.2e1 * (MDP(14) * t247 - t212) * t174; (t135 * t171 + t184 * t140) * MDP(30) + (t135 * t172 + t184 * t141) * MDP(31) - t175 * t252 + t246 + t176 * MDP(9) - t177 * MDP(10) + (-t172 * t140 - t141 * t171) * MDP(26) + (-t123 * t172 - t124 * t171 - t133 * t141 - t134 * t140) * MDP(32) + (t123 * t133 + t124 * t134 + t129 * t153) * MDP(33) + t141 * t238 + (t210 * MDP(7) + t244) * t202 + (t215 + t218 - t253 + t272) * t174 + (-t167 * MDP(16) + (-MDP(14) + t251) * t247 - t233 * t174 + t212 - t216 + t270) * t209 + ((pkin(9) * t154 + t250) * MDP(23) + (pkin(9) * t155 + t249) * MDP(24) + t167 * MDP(17) + (-t154 * t208 - t155 * t205) * MDP(19) + t155 * t235 + t175 * MDP(11) + t223 * t247 + t214 * t174) * t206; MDP(8) + t252 * t274 + (t133 ^ 2 + t134 ^ 2 + t153 ^ 2) * MDP(33) + (-t134 * t264 + t184 * t266) * t171 + (-t133 * t264 + t171 * t267 + t184 * t265 + t238) * t172 + (MDP(18) * t200 + t256 * t268 + t257 * t269 + MDP(11) - 0.2e1 * t225) * t206 ^ 2 + (-t138 * t266 + t139 * t265 - t162 * t269 + t163 * t268 + t209 * t233 + t214 * t274 + 0.2e1 * t165 - 0.2e1 * t166 + 0.2e1 * t253) * t209; t237 - t202 * t243 + t146 * MDP(16) - t147 * MDP(17) + t205 * t239 + (-t154 * t205 + t155 * t208) * MDP(19) + (-pkin(3) * t154 - t249) * MDP(23) + (-pkin(3) * t155 + t250) * MDP(24) + t182 * t240 + (-t140 * t182 - t141 * t181) * MDP(26) + (t135 * t181 + t140 * t195) * MDP(30) + (t135 * t182 + t141 * t195) * MDP(31) + (-t123 * t182 - t124 * t181 - t140 * t149 - t141 * t148) * MDP(32) + (t123 * t148 + t124 * t149 + t129 * t164) * MDP(33) + t276 * t174; t182 * t238 + (-t171 * t182 - t172 * t181) * MDP(26) + (t171 * t195 + t181 * t184) * MDP(30) + (t172 * t195 + t182 * t184) * MDP(31) + (-t133 * t182 - t134 * t181 - t148 * t172 - t149 * t171) * MDP(32) + (t133 * t148 + t134 * t149 + t153 * t164) * MDP(33) + (-t251 - t276) * t209 + (t205 * t235 + (-t198 + t200) * MDP(19) + (-pkin(3) * t205 - t256) * MDP(23) + (-pkin(3) * t208 + t257) * MDP(24) - t223) * t206; MDP(15) + t198 * MDP(18) + 0.2e1 * t225 - t149 * t181 * t264 + (t148 ^ 2 + t149 ^ 2 + t164 ^ 2) * MDP(33) + (MDP(25) * t182 - t148 * t264 + t181 * t267) * t182 + 0.2e1 * (MDP(30) * t181 + MDP(31) * t182) * t195 + 0.2e1 * (MDP(23) * t208 - MDP(24) * t205) * pkin(3); t174 * MDP(22) + (t174 * t230 + t125) * MDP(30) + (-t228 + (-t128 - t259) * t204) * MDP(31) + (-t140 * t258 - t141 * t194) * MDP(32) + (t123 * t194 + t124 * t258) * MDP(33) + t220 + t229 - t270; -t209 * MDP(22) + (-t209 * t230 + t138) * MDP(30) + (-t227 + (pkin(4) * t209 - t152) * t204) * MDP(31) + (-t171 * t258 - t172 * t194) * MDP(32) + (t133 * t194 + t134 * t258) * MDP(33) + t218 + t221 + t271 * t206; (-t181 * t258 - t182 * t194) * MDP(32) + (t148 * t194 + t149 * t258) * MDP(33) + t277; t194 ^ 2 * MDP(33) + (0.2e1 * t224 + (MDP(33) * t258 - 0.2e1 * MDP(31)) * t204) * pkin(4) + t233; (-MDP(32) * t141 + MDP(33) * t123) * pkin(5) + t216 + t229; (-t172 * MDP(32) + t133 * MDP(33)) * pkin(5) + t215 + t221; (-MDP(32) * t182 + MDP(33) * t148) * pkin(5) + t222; t194 * pkin(5) * MDP(33) + MDP(29) + (-MDP(31) * t204 + t224) * pkin(4); MDP(33) * pkin(5) ^ 2 + MDP(29); t129 * MDP(33); t153 * MDP(33); t164 * MDP(33); 0; 0; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
