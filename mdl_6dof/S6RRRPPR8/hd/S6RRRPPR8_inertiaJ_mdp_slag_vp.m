% Calculate joint inertia matrix for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPPR8_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:09:33
% EndTime: 2019-03-09 16:09:36
% DurationCPUTime: 1.08s
% Computational Cost: add. (987->239), mult. (2107->318), div. (0->0), fcn. (2101->8), ass. (0->98)
t167 = cos(pkin(6));
t170 = sin(qJ(3));
t173 = cos(qJ(3));
t166 = sin(pkin(6));
t171 = sin(qJ(2));
t212 = t166 * t171;
t139 = t167 * t170 + t173 * t212;
t227 = 0.2e1 * t139;
t226 = pkin(9) - qJ(5);
t174 = cos(qJ(2));
t211 = t166 * t174;
t196 = MDP(19) - MDP(24);
t225 = -t196 * qJ(4) - MDP(14);
t224 = 2 * MDP(18);
t223 = 2 * MDP(23);
t222 = 2 * MDP(24);
t221 = 2 * MDP(31);
t220 = 2 * MDP(32);
t219 = -pkin(4) - pkin(10);
t218 = pkin(1) * t171;
t217 = pkin(1) * t174;
t168 = qJ(4) + pkin(5);
t216 = MDP(19) * pkin(9);
t215 = MDP(21) * pkin(9);
t214 = pkin(2) * MDP(17);
t146 = t226 * t173;
t213 = t146 * t173;
t210 = t167 * MDP(8);
t209 = t171 * MDP(6);
t194 = pkin(8) * t211;
t132 = t194 + (pkin(9) + t218) * t167;
t133 = (-pkin(2) * t174 - pkin(9) * t171 - pkin(1)) * t166;
t208 = t170 * t132 - t173 * t133;
t121 = t173 * t132 + t170 * t133;
t140 = -pkin(8) * t212 + t167 * t217;
t207 = MDP(15) * t174;
t206 = MDP(21) * pkin(9) ^ 2;
t138 = -t167 * t173 + t170 * t212;
t169 = sin(qJ(6));
t172 = cos(qJ(6));
t125 = t138 * t172 + t169 * t211;
t205 = MDP(26) * t125;
t204 = MDP(26) * t172;
t124 = t138 * t169 - t172 * t211;
t203 = t124 * MDP(29);
t202 = t138 * MDP(12);
t201 = t138 * MDP(23);
t143 = -t173 * pkin(3) - t170 * qJ(4) - pkin(2);
t142 = t173 * pkin(4) - t143;
t200 = t142 * MDP(22);
t199 = t143 * MDP(20);
t198 = MDP(11) + MDP(30);
t197 = MDP(17) - MDP(20);
t195 = 0.2e1 * t173;
t193 = qJ(4) * t211;
t152 = pkin(3) * t211;
t119 = t152 + t208;
t192 = -0.2e1 * t193 + t121;
t131 = -t167 * pkin(2) - t140;
t191 = t169 * t172 * MDP(27);
t190 = -pkin(3) * MDP(21) - MDP(18);
t175 = -pkin(3) - pkin(4);
t189 = MDP(25) * t175 + MDP(23);
t188 = pkin(4) * t211 + t119;
t117 = t138 * pkin(3) - t139 * qJ(4) + t131;
t187 = -MDP(28) * t172 + MDP(29) * t169;
t126 = pkin(5) * t170 + pkin(10) * t173 + t142;
t145 = t226 * t170;
t122 = t126 * t172 - t145 * t169;
t123 = t126 * t169 + t145 * t172;
t186 = t122 * MDP(31) - t123 * MDP(32);
t185 = MDP(31) * t172 - MDP(32) * t169;
t184 = -MDP(31) * t169 - MDP(32) * t172;
t118 = -t193 + t121;
t183 = MDP(12) + t187;
t182 = MDP(22) + t185;
t114 = -qJ(5) * t139 + t188;
t181 = t184 + t196;
t111 = t139 * pkin(5) + t219 * t138 - t117;
t112 = pkin(10) * t211 + t114;
t109 = t111 * t172 - t112 * t169;
t110 = t111 * t169 + t112 * t172;
t180 = t125 * MDP(28) + t109 * MDP(31) - t110 * MDP(32) - t203;
t179 = -t169 * MDP(28) - t172 * MDP(29) + (-pkin(3) + t219) * t184;
t178 = -pkin(3) * MDP(19) - t175 * MDP(24) + MDP(13) + t179;
t176 = qJ(4) ^ 2;
t165 = t173 ^ 2;
t164 = t172 ^ 2;
t163 = t170 ^ 2;
t162 = t169 ^ 2;
t160 = t166 ^ 2;
t144 = t170 * pkin(9) * t211;
t141 = t167 * t218 + t194;
t134 = t138 * qJ(5);
t116 = -t118 - t134;
t115 = -pkin(4) * t138 - t117;
t113 = -t168 * t211 + t121 + t134;
t1 = [(t114 ^ 2 + t115 ^ 2 + t116 ^ 2) * MDP(25) + (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) * MDP(21) + MDP(1) + t160 * t171 ^ 2 * MDP(4) + (0.2e1 * t166 * t209 + t210) * t167 + t198 * t139 ^ 2 + (-0.2e1 * MDP(27) * t124 + MDP(28) * t227 + t205) * t125 + (0.2e1 * MDP(5) * t171 + t207) * t160 * t174 + (t109 * t139 + t113 * t124) * t221 + 0.2e1 * (-t118 * t138 + t119 * t139) * MDP(19) + (-t110 * t139 + t113 * t125) * t220 + (-t114 * t139 - t116 * t138) * t222 + 0.2e1 * (-t141 * t167 - t160 * t218) * MDP(10) + (-t114 * t211 + t115 * t138) * t223 + (t117 * t138 + t119 * t211) * t224 + 0.2e1 * (t131 * t138 + t208 * t211) * MDP(16) + 0.2e1 * (t115 * t139 + t116 * t211) * MDP(22) + 0.2e1 * (t121 * t211 + t131 * t139) * MDP(17) + 0.2e1 * (-t117 * t139 - t118 * t211) * MDP(20) + 0.2e1 * (t140 * t167 + t160 * t217) * MDP(9) + (-t202 - t203) * t227 + 0.2e1 * (-t139 * MDP(13) + t138 * MDP(14) + MDP(7) * t167) * t211; t140 * MDP(9) - t141 * MDP(10) + (t114 * t145 + t115 * t142) * MDP(25) + t210 + (-pkin(2) * t138 + t144) * MDP(16) + (t138 * t143 + t144) * MDP(18) + t117 * t143 * MDP(21) + t142 * t201 + (t209 + (-t145 * MDP(23) + MDP(7)) * t174) * t166 + (-MDP(22) * t211 + MDP(24) * t138 - MDP(25) * t116 + MDP(31) * t124 + MDP(32) * t125) * t146 + (-MDP(24) * t145 + t186 - t199 + t200 - t214) * t139 + (-MDP(13) * t211 - t202 + t131 * MDP(17) - t117 * MDP(20) + t115 * MDP(22) - t114 * MDP(24) + (MDP(19) + t215) * t119 + (t198 + t216) * t139 + t180) * t170 + (-t131 * MDP(16) - t117 * MDP(18) + t118 * MDP(19) + t116 * MDP(24) + (t124 * t172 + t125 * t169) * MDP(27) - t115 * MDP(23) - t125 * t204 - MDP(14) * t211 + t184 * t113 + t183 * t139 + (-t138 * MDP(19) + t118 * MDP(21) + t197 * t211) * pkin(9)) * t173; MDP(8) + t143 ^ 2 * MDP(21) + (t142 ^ 2 + t145 ^ 2 + t146 ^ 2) * MDP(25) + (MDP(26) * t164 - 0.2e1 * t191 + t206) * t165 + (t198 + t206) * t163 + (MDP(16) * pkin(2) - MDP(18) * t143 - MDP(23) * t142) * t195 + (t183 * t195 - 0.2e1 * t199 + 0.2e1 * t200 - 0.2e1 * t214) * t170 + (-t145 * t170 - t213) * t222 + (t122 * t170 - t169 * t213) * t221 + (-t123 * t170 - t172 * t213) * t220 + 0.2e1 * (t163 + t165) * t216; -t166 * t207 - t208 * MDP(16) - t121 * MDP(17) + (-0.2e1 * t152 - t208) * MDP(18) + t192 * MDP(20) + (-pkin(3) * t119 + qJ(4) * t118) * MDP(21) + (t134 + t192) * MDP(22) + (-t175 * t211 + t188) * MDP(23) + (-qJ(4) * t116 + t114 * t175) * MDP(25) - t169 * t205 + (t124 * t169 - t125 * t172) * MDP(27) + (t113 * t172 + t124 * t168) * MDP(31) + (-t113 * t169 + t125 * t168) * MDP(32) + t225 * t138 + (-qJ(5) * MDP(23) + t178) * t139; t189 * t145 + (MDP(25) * qJ(4) + t182) * t146 + t178 * t170 + (t169 * t204 + (-t162 + t164) * MDP(27) + t184 * t168 - t225) * t173 + ((MDP(21) * qJ(4) - t197) * t173 + (-MDP(16) + t190) * t170) * pkin(9); MDP(15) + pkin(3) * t224 + (pkin(3) ^ 2 + t176) * MDP(21) + t175 * t223 + (t175 ^ 2 + t176) * MDP(25) + t162 * MDP(26) + 0.2e1 * t191 + 0.2e1 * t185 * t168 + 0.2e1 * (MDP(20) + MDP(22)) * qJ(4); t119 * MDP(21) + t114 * MDP(25) + (MDP(18) - MDP(23)) * t211 + t181 * t139; t145 * MDP(25) + (t181 + t215) * t170; t189 + t190; MDP(21) + MDP(25); t115 * MDP(25) + t139 * t182 + t201; -t173 * MDP(23) + t142 * MDP(25) + t170 * t182; 0; 0; MDP(25); t139 * MDP(30) + t180; t170 * MDP(30) + t173 * t187 + t186; t179; t184; t185; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
