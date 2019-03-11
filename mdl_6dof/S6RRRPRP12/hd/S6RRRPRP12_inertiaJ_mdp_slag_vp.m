% Calculate joint inertia matrix for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP12_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:59:30
% EndTime: 2019-03-09 17:59:34
% DurationCPUTime: 1.56s
% Computational Cost: add. (1569->300), mult. (3298->413), div. (0->0), fcn. (3389->8), ass. (0->116)
t252 = pkin(4) + pkin(9);
t183 = cos(qJ(3));
t251 = 0.2e1 * t183;
t179 = sin(qJ(5));
t172 = t179 ^ 2;
t182 = cos(qJ(5));
t174 = t182 ^ 2;
t165 = t172 + t174;
t250 = t165 * MDP(32);
t178 = cos(pkin(6));
t180 = sin(qJ(3));
t177 = sin(pkin(6));
t181 = sin(qJ(2));
t230 = t177 * t181;
t194 = -t178 * t183 + t180 * t230;
t249 = t194 * MDP(14);
t248 = -MDP(17) + MDP(20);
t152 = t178 * t180 + t183 * t230;
t247 = -0.2e1 * t152;
t246 = 2 * MDP(16);
t245 = 2 * MDP(18);
t244 = 2 * MDP(19);
t243 = -2 * MDP(20);
t242 = 2 * MDP(20);
t241 = 2 * MDP(27);
t240 = 2 * MDP(28);
t239 = 2 * MDP(29);
t238 = 2 * MDP(31);
t185 = -pkin(3) - pkin(10);
t237 = pkin(1) * t181;
t236 = t152 * pkin(5);
t235 = t180 * pkin(5);
t184 = cos(qJ(2));
t234 = t184 * pkin(1);
t233 = (pkin(3) * MDP(21));
t229 = t177 * t184;
t140 = -t194 * t179 + t182 * t229;
t232 = t140 * t182;
t231 = t152 * qJ(6);
t228 = t178 * MDP(8);
t139 = t179 * t229 + t194 * t182;
t227 = t179 * t139;
t226 = t179 * t183;
t225 = t180 * qJ(6);
t224 = t180 * t152;
t223 = t180 * t185;
t222 = t182 * t183;
t208 = pkin(8) * t229;
t147 = t208 + (pkin(9) + t237) * t178;
t148 = (-pkin(2) * t184 - pkin(9) * t181 - pkin(1)) * t177;
t133 = -t180 * t147 + t183 * t148;
t167 = pkin(3) * t229;
t132 = -t133 + t167;
t125 = t152 * pkin(4) + pkin(10) * t229 + t132;
t166 = pkin(8) * t230;
t146 = t166 + (-pkin(2) - t234) * t178;
t130 = t194 * pkin(3) - t152 * qJ(4) + t146;
t128 = t194 * pkin(10) + t130;
t121 = t179 * t125 + t182 * t128;
t134 = t183 * t147 + t180 * t148;
t203 = -t180 * qJ(4) - pkin(2);
t155 = t185 * t183 + t203;
t162 = t252 * t180;
t138 = t182 * t155 + t179 * t162;
t163 = t252 * t183;
t221 = MDP(21) * qJ(4);
t193 = t179 * pkin(5) - t182 * qJ(6);
t159 = qJ(4) + t193;
t220 = MDP(32) * t159;
t219 = t132 * MDP(21);
t218 = t139 * MDP(23);
t217 = t139 * MDP(25);
t216 = t140 * MDP(22);
t215 = t179 * MDP(23);
t214 = t180 * MDP(26);
t213 = t184 * MDP(15);
t212 = t185 * MDP(32);
t211 = pkin(9) ^ 2 * MDP(21);
t210 = MDP(11) + MDP(26);
t209 = -MDP(28) + MDP(31);
t207 = t180 * t229;
t206 = MDP(6) * t230;
t205 = t183 * t229;
t204 = t182 * t215;
t202 = -pkin(3) * MDP(18) + MDP(13);
t201 = MDP(19) - t233;
t200 = -MDP(32) * pkin(5) - MDP(29);
t199 = -t182 * t125 + t179 * t128;
t198 = t179 * t155 - t182 * t162;
t197 = pkin(9) * t207;
t196 = pkin(9) * t205;
t195 = (MDP(27) + MDP(29)) * t182;
t161 = pkin(5) * t182 + t179 * qJ(6);
t131 = qJ(4) * t229 - t134;
t192 = -t131 * t183 + t132 * t180;
t191 = t182 * MDP(24) - t179 * MDP(25);
t190 = -t179 * MDP(24) - t182 * MDP(25);
t189 = -MDP(27) * t198 - t138 * MDP(28);
t188 = t209 * t179 + t195;
t187 = MDP(18) + t188;
t129 = -t194 * pkin(4) - t131;
t175 = t183 ^ 2;
t173 = t180 ^ 2;
t171 = t177 ^ 2;
t160 = -t183 * pkin(3) + t203;
t158 = t165 * t185;
t154 = t178 * t237 + t208;
t153 = t178 * t234 - t166;
t142 = t161 * t183 + t163;
t136 = t198 - t235;
t135 = t138 + t225;
t126 = t135 * t179 - t136 * t182;
t122 = -t139 * pkin(5) + t140 * qJ(6) + t129;
t119 = t199 - t236;
t118 = t121 + t231;
t1 = [(t118 ^ 2 + t119 ^ 2 + t122 ^ 2) * MDP(32) + (t130 ^ 2 + t131 ^ 2 + t132 ^ 2) * MDP(21) + t171 * t181 ^ 2 * MDP(4) + MDP(1) + (0.2e1 * t206 + t228) * t178 + 0.2e1 * (-t194 * MDP(12) + t217) * t152 + t210 * t152 ^ 2 + (MDP(24) * t247 + t216 - 0.2e1 * t218) * t140 + ((MDP(13) * t247 + 0.2e1 * t178 * MDP(7) + 0.2e1 * t249) * t177 + (0.2e1 * t181 * MDP(5) + t213) * t171) * t184 + (-t133 * t229 + t146 * t194) * t246 + (-t130 * t194 - t132 * t229) * t244 + (-t130 * t152 + t131 * t229) * t242 + 0.2e1 * (t134 * t229 + t146 * t152) * MDP(17) + 0.2e1 * (t153 * t178 + t171 * t234) * MDP(9) + (t131 * t194 + t132 * t152) * t245 + 0.2e1 * (-t154 * t178 - t171 * t237) * MDP(10) + (-t119 * t152 - t122 * t139) * t239 + (-t121 * t152 - t129 * t140) * t240 + (-t129 * t139 - t152 * t199) * t241 + (t118 * t152 + t122 * t140) * t238 + 0.2e1 * (t118 * t139 - t119 * t140) * MDP(30); t228 + t206 + (t129 * t222 - t163 * t139 - t152 * t198 - t180 * t199) * MDP(27) + (t118 * t135 + t119 * t136 + t122 * t142) * MDP(32) + (t192 * pkin(9) + t130 * t160) * MDP(21) + (-t227 + t232) * t183 * MDP(23) + (t152 * t183 - t180 * t194) * MDP(12) + (-t130 * t180 - t160 * t152 - t196) * MDP(20) + (-pkin(2) * t152 + t146 * t180 + t196) * MDP(17) + (t130 * t183 - t160 * t194 - t197) * MDP(19) + (-pkin(2) * t194 - t146 * t183 + t197) * MDP(16) + t153 * MDP(9) - t154 * MDP(10) + (t135 * t139 - t136 * t140 + (-t118 * t182 - t119 * t179) * t183) * MDP(30) + (-t140 * t180 - t152 * t226) * MDP(24) + (-t121 * t180 - t129 * t226 - t138 * t152 - t163 * t140) * MDP(28) + (t118 * t180 + t122 * t226 + t135 * t152 + t142 * t140) * MDP(31) - MDP(14) * t205 - MDP(13) * t207 + (t139 * t180 - t152 * t222) * MDP(25) + (-t119 * t180 + t122 * t222 - t136 * t152 - t142 * t139) * MDP(29) + ((-t183 * t194 + t224) * pkin(9) + t192) * MDP(18) + MDP(11) * t224 + t216 * t226 + MDP(7) * t229 + t152 * t214; MDP(8) + pkin(2) * t183 * t246 + (t135 ^ 2 + t136 ^ 2 + t142 ^ 2) * MDP(32) + (MDP(21) * t160 + t183 * t244) * t160 + (t172 * MDP(22) + 0.2e1 * t204 + t211) * t175 + (t210 + t211) * t173 + (-0.2e1 * pkin(2) * MDP(17) + t160 * t243 + (MDP(12) + t190) * t251) * t180 + 0.2e1 * (-t136 * MDP(29) + t135 * MDP(31) + t189) * t180 + ((-t135 * t182 - t136 * t179) * MDP(30) + (t182 * MDP(27) - t179 * MDP(28)) * t163 + (t182 * MDP(29) + t179 * MDP(31)) * t142) * t251 + (t173 + t175) * pkin(9) * t245; -t249 - t177 * t213 + t133 * MDP(16) + (-t133 + 0.2e1 * t167) * MDP(19) - pkin(3) * t219 + (-t139 * MDP(29) + t122 * MDP(32)) * t159 + (-t194 * MDP(18) - t131 * MDP(21) - t139 * MDP(27) + t229 * t243) * qJ(4) + (t129 * MDP(27) + t122 * MDP(29) + (t139 * t185 - t118) * MDP(30) + t118 * t212) * t179 + (-qJ(4) * MDP(28) + t159 * MDP(31) + t215) * t140 + (-t216 + t218 + t129 * MDP(28) + (t140 * t185 + t119) * MDP(30) - t122 * MDP(31) - t119 * t212) * t182 + (t188 * t185 + t191 + t202) * t152 + t248 * t134; t142 * t220 - t126 * MDP(30) + t202 * t180 + t195 * t223 + (MDP(14) + qJ(4) * MDP(18) + (t172 - t174) * MDP(23)) * t183 + (-t136 * t212 + t180 * MDP(24) + t163 * MDP(28) - t142 * MDP(31) + (qJ(4) * MDP(27) + t159 * MDP(29)) * t183) * t182 + (-MDP(22) * t222 - t180 * MDP(25) + t163 * MDP(27) + (-qJ(4) * t183 - t223) * MDP(28) + t142 * MDP(29) + (t159 * t183 + t223) * MDP(31) + t135 * t212) * t179 + ((t221 + t248) * t183 + (-MDP(16) + t201) * t180) * pkin(9); -0.2e1 * t204 - 0.2e1 * t158 * MDP(30) + t174 * MDP(22) + MDP(15) + t185 ^ 2 * t250 + (-0.2e1 * t182 * MDP(31) + t179 * t239 + t220) * t159 + (-2 * MDP(19) + t233) * pkin(3) + (t179 * t241 + t182 * t240 + t221 + t242) * qJ(4); -MDP(19) * t229 + t219 + (t227 + t232) * MDP(30) + (t118 * t179 - t119 * t182) * MDP(32) + t187 * t152; t126 * MDP(32) + (pkin(9) * MDP(21) + t187) * t180; -t165 * MDP(30) + t158 * MDP(32) + t201; MDP(21) + t250; -t140 * MDP(24) + t217 + t152 * MDP(26) - t199 * MDP(27) - t121 * MDP(28) + (-t199 + 0.2e1 * t236) * MDP(29) + (pkin(5) * t140 + qJ(6) * t139) * MDP(30) + (t121 + 0.2e1 * t231) * MDP(31) + (-pkin(5) * t119 + qJ(6) * t118) * MDP(32); t214 + (-t198 + 0.2e1 * t235) * MDP(29) + (t138 + 0.2e1 * t225) * MDP(31) + (-pkin(5) * t136 + qJ(6) * t135) * MDP(32) + (t193 * MDP(30) + t190) * t183 + t189; -t161 * MDP(30) + ((MDP(27) - t200) * t182 + (MDP(32) * qJ(6) + t209) * t179) * t185 + t191; t161 * MDP(32) + t188; MDP(26) + pkin(5) * t239 + qJ(6) * t238 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(32); -t152 * MDP(29) - t140 * MDP(30) + t119 * MDP(32); -t180 * MDP(29) - MDP(30) * t226 + t136 * MDP(32); (MDP(30) - t212) * t182; -t182 * MDP(32); t200; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
