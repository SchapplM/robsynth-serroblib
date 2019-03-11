% Calculate joint inertia matrix for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR4_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:00:04
% EndTime: 2019-03-09 01:00:08
% DurationCPUTime: 1.17s
% Computational Cost: add. (1357->235), mult. (3253->352), div. (0->0), fcn. (3784->14), ass. (0->114)
t205 = sin(qJ(6));
t210 = cos(qJ(6));
t247 = t205 * MDP(28) + t210 * MDP(29);
t221 = t210 * MDP(31) - t205 * MDP(32);
t203 = cos(pkin(7));
t207 = sin(qJ(4));
t211 = cos(qJ(4));
t201 = sin(pkin(7));
t208 = sin(qJ(3));
t252 = t201 * t208;
t175 = -t203 * t211 + t207 * t252;
t176 = t203 * t207 + t211 * t252;
t206 = sin(qJ(5));
t262 = cos(qJ(5));
t154 = t262 * t175 + t176 * t206;
t155 = -t206 * t175 + t262 * t176;
t269 = t155 * MDP(21) - t154 * MDP(22);
t268 = -(t207 * MDP(17) + t211 * MDP(18)) * pkin(10) + t207 * MDP(14) + t211 * MDP(15);
t263 = pkin(10) + pkin(11);
t187 = t263 * t211;
t230 = t262 * t207;
t164 = t187 * t206 + t263 * t230;
t249 = t206 * t207;
t165 = t262 * t187 - t263 * t249;
t182 = -t262 * t211 + t249;
t183 = t206 * t211 + t230;
t267 = t183 * MDP(21) - t182 * MDP(22) - t164 * MDP(24) - t165 * MDP(25);
t266 = 0.2e1 * MDP(24);
t265 = 0.2e1 * MDP(31);
t264 = 0.2e1 * MDP(32);
t261 = pkin(2) * t208;
t212 = cos(qJ(3));
t260 = pkin(2) * t212;
t251 = t201 * t212;
t233 = pkin(9) * t251;
t169 = t233 + (pkin(10) + t261) * t203;
t170 = (-pkin(3) * t212 - pkin(10) * t208 - pkin(2)) * t201;
t148 = -t169 * t207 + t211 * t170;
t234 = pkin(4) * t251;
t138 = -pkin(11) * t176 + t148 - t234;
t149 = t169 * t211 + t170 * t207;
t141 = -pkin(11) * t175 + t149;
t128 = t262 * t138 - t206 * t141;
t126 = pkin(5) * t251 - t128;
t258 = t126 * t210;
t257 = t154 * t205;
t256 = t154 * t210;
t255 = t164 * t210;
t254 = t183 * t205;
t253 = t183 * t210;
t213 = cos(qJ(2));
t250 = t203 * t213;
t248 = t208 * MDP(7);
t246 = MDP(30) * t182;
t146 = t155 * t205 + t210 * t251;
t245 = t146 * MDP(29);
t147 = t155 * t210 - t205 * t251;
t244 = t147 * MDP(26);
t243 = t147 * MDP(28);
t242 = t154 * MDP(30);
t241 = t155 * MDP(20);
t240 = t176 * MDP(12);
t239 = t183 * MDP(25);
t237 = t210 * MDP(26);
t235 = t211 * MDP(17);
t232 = t262 * pkin(4);
t193 = -pkin(4) * t211 - pkin(3);
t231 = t262 * t141;
t229 = t205 * t210 * MDP(27);
t198 = t205 ^ 2;
t228 = t198 * MDP(26) + MDP(23) + 0.2e1 * t229;
t227 = -pkin(5) * t183 - pkin(12) * t182;
t191 = pkin(4) * t206 + pkin(12);
t192 = -t232 - pkin(5);
t226 = -t182 * t191 + t183 * t192;
t225 = t176 * MDP(14) - t175 * MDP(15);
t222 = MDP(28) * t210 - MDP(29) * t205;
t220 = t205 * MDP(31) + t210 * MDP(32);
t189 = pkin(9) * t252;
t168 = t189 + (-pkin(3) - t260) * t203;
t129 = t206 * t138 + t231;
t219 = -MDP(20) + t222;
t202 = sin(pkin(6));
t204 = cos(pkin(6));
t209 = sin(qJ(2));
t160 = t204 * t252 + (t208 * t250 + t209 * t212) * t202;
t174 = -t201 * t202 * t213 + t204 * t203;
t144 = -t160 * t207 + t174 * t211;
t145 = t160 * t211 + t174 * t207;
t134 = -t262 * t144 + t145 * t206;
t135 = t206 * t144 + t262 * t145;
t218 = -t135 * MDP(25) + (-MDP(24) - t221) * t134;
t217 = (t262 * MDP(24) - t206 * MDP(25)) * pkin(4);
t199 = t210 ^ 2;
t216 = (-t198 + t199) * t183 * MDP(27) + t237 * t254 + t267 + t247 * t182;
t156 = pkin(4) * t175 + t168;
t215 = -MDP(23) * t251 + (-t146 * t205 + t147 * t210) * MDP(27) + t205 * t244 + t269 + t247 * t154;
t127 = -pkin(12) * t251 + t129;
t133 = pkin(5) * t154 - pkin(12) * t155 + t156;
t121 = -t127 * t205 + t133 * t210;
t122 = t127 * t210 + t133 * t205;
t214 = t121 * MDP(31) - t122 * MDP(32) + t242 + t243 - t245;
t197 = t201 ^ 2;
t178 = t203 * t261 + t233;
t177 = t203 * t260 - t189;
t161 = t164 * t205;
t159 = -t204 * t251 + (t208 * t209 - t212 * t250) * t202;
t158 = pkin(5) * t182 - pkin(12) * t183 + t193;
t140 = t158 * t205 + t165 * t210;
t139 = t158 * t210 - t165 * t205;
t125 = t126 * t205;
t124 = t135 * t210 + t159 * t205;
t123 = -t135 * t205 + t159 * t210;
t1 = [MDP(1); (-t159 * t203 - t174 * t251) * MDP(10) + (-t160 * t203 + t174 * t252) * MDP(11) + (-t144 * t251 + t159 * t175) * MDP(17) + (t145 * t251 + t159 * t176) * MDP(18) + (t134 * t251 + t154 * t159) * MDP(24) + (t135 * t251 + t155 * t159) * MDP(25) + (t123 * t154 + t134 * t146) * MDP(31) + (-t124 * t154 + t134 * t147) * MDP(32) + (t213 * MDP(3) - t209 * MDP(4)) * t202; t155 ^ 2 * MDP(19) + t203 ^ 2 * MDP(9) + MDP(2) + (-0.2e1 * t175 * MDP(13) + t240) * t176 + (-0.2e1 * t146 * MDP(27) + t244) * t147 + ((MDP(5) * t208 + 0.2e1 * MDP(6) * t212) * t208 + (MDP(16) + MDP(23)) * t212 ^ 2) * t197 + (-0.2e1 * t241 + t242 + 0.2e1 * t243 - 0.2e1 * t245) * t154 + 0.2e1 * (t203 * t248 + (MDP(8) * t203 - t225 - t269) * t212) * t201 + 0.2e1 * (t149 * t251 + t168 * t176) * MDP(18) + 0.2e1 * (t129 * t251 + t155 * t156) * MDP(25) + (-t128 * t251 + t154 * t156) * t266 + 0.2e1 * (t177 * t203 + t197 * t260) * MDP(10) + 0.2e1 * (-t148 * t251 + t168 * t175) * MDP(17) + 0.2e1 * (-t178 * t203 - t197 * t261) * MDP(11) + (-t122 * t154 + t126 * t147) * t264 + (t121 * t154 + t126 * t146) * t265; -t160 * MDP(11) + (t123 * t182 + t134 * t254) * MDP(31) + (-t124 * t182 + t134 * t253) * MDP(32) + (t207 * MDP(18) + t182 * MDP(24) - MDP(10) - t235 + t239) * t159; t203 * MDP(9) + t177 * MDP(10) - t178 * MDP(11) + (-t175 * t207 + t176 * t211) * MDP(13) + (-pkin(3) * t175 - t168 * t211) * MDP(17) + (-pkin(3) * t176 + t168 * t207) * MDP(18) + (t139 * t154 + t146 * t164) * MDP(31) + (-t140 * t154 + t147 * t164) * MDP(32) + t207 * t240 + (t154 * MDP(24) + t155 * MDP(25)) * t193 + (t156 * MDP(24) + t214 - t241) * t182 + (t156 * MDP(25) + (-t146 * t210 - t147 * t205) * MDP(27) + t155 * MDP(19) + t147 * t237 + t220 * t126 + t219 * t154) * t183 + (t248 + (MDP(8) - t267 - t268) * t212) * t201; 0.2e1 * pkin(3) * t235 + 0.2e1 * t193 * t239 + MDP(9) + (0.2e1 * t219 * t183 + t193 * t266 + t246) * t182 + (t139 * t182 + t164 * t254) * t265 + (-t140 * t182 + t164 * t253) * t264 + (MDP(12) * t207 + 0.2e1 * t211 * MDP(13) - 0.2e1 * pkin(3) * MDP(18)) * t207 + (t199 * MDP(26) + MDP(19) - 0.2e1 * t229) * t183 ^ 2; MDP(17) * t144 - MDP(18) * t145 + t218; -MDP(16) * t251 + t148 * MDP(17) - t149 * MDP(18) + (-t232 * t251 + t128) * MDP(24) + (-t231 + (-t138 + t234) * t206) * MDP(25) + (t146 * t192 - t191 * t257 - t258) * MDP(31) + (t147 * t192 - t191 * t256 + t125) * MDP(32) + t215 + t225; (t226 * t205 - t255) * MDP(31) + (t226 * t210 + t161) * MDP(32) + t216 + t268; -0.2e1 * t192 * t221 + MDP(16) + 0.2e1 * t217 + t228; t218; t128 * MDP(24) - t129 * MDP(25) + (-pkin(5) * t146 - pkin(12) * t257 - t258) * MDP(31) + (-pkin(5) * t147 - pkin(12) * t256 + t125) * MDP(32) + t215; (t227 * t205 - t255) * MDP(31) + (t227 * t210 + t161) * MDP(32) + t216; t217 + t228 + t221 * (pkin(5) - t192); 0.2e1 * pkin(5) * t221 + t228; t123 * MDP(31) - t124 * MDP(32); t214; t139 * MDP(31) - t140 * MDP(32) + t222 * t183 + t246; -t220 * t191 + t247; -t220 * pkin(12) + t247; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
