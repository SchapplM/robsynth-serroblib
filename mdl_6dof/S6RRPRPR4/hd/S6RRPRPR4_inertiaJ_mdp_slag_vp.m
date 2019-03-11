% Calculate joint inertia matrix for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRPR4_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:26:20
% EndTime: 2019-03-09 10:26:22
% DurationCPUTime: 0.77s
% Computational Cost: add. (1932->205), mult. (4553->325), div. (0->0), fcn. (5160->12), ass. (0->97)
t170 = sin(pkin(12));
t173 = cos(pkin(12));
t177 = sin(qJ(4));
t180 = cos(qJ(4));
t155 = t170 * t177 - t173 * t180;
t176 = sin(qJ(6));
t179 = cos(qJ(6));
t190 = t179 * MDP(27) - t176 * MDP(28);
t219 = t190 * t155;
t218 = 2 * MDP(19);
t217 = 2 * MDP(20);
t216 = 0.2e1 * MDP(27);
t215 = 0.2e1 * MDP(28);
t181 = cos(qJ(2));
t214 = pkin(1) * t181;
t213 = pkin(8) + qJ(3);
t157 = t170 * t180 + t173 * t177;
t212 = t157 * t176;
t211 = t157 * t179;
t172 = sin(pkin(6));
t167 = t172 ^ 2;
t178 = sin(qJ(2));
t210 = t167 * t178;
t209 = t172 * t178;
t208 = t172 * t181;
t175 = cos(pkin(6));
t207 = t175 * MDP(8);
t171 = sin(pkin(11));
t164 = t171 * pkin(2) + pkin(9);
t206 = qJ(5) + t164;
t161 = t175 * t214;
t143 = t175 * pkin(2) - t213 * t209 + t161;
t194 = t175 * t178 * pkin(1);
t146 = t213 * t208 + t194;
t174 = cos(pkin(11));
t132 = t171 * t143 + t174 * t146;
t130 = t175 * pkin(9) + t132;
t148 = t171 * t209 - t174 * t208;
t149 = (t171 * t181 + t174 * t178) * t172;
t158 = (-pkin(2) * t181 - pkin(1)) * t172;
t134 = t148 * pkin(3) - t149 * pkin(9) + t158;
t119 = -t177 * t130 + t180 * t134;
t142 = t149 * t180 + t175 * t177;
t116 = t148 * pkin(4) - t142 * qJ(5) + t119;
t120 = t180 * t130 + t177 * t134;
t141 = t149 * t177 - t175 * t180;
t118 = -t141 * qJ(5) + t120;
t113 = t170 * t116 + t173 * t118;
t128 = -t170 * t141 + t173 * t142;
t122 = t176 * t128 - t148 * t179;
t205 = t122 * MDP(23);
t204 = t122 * MDP(25);
t123 = t179 * t128 + t148 * t176;
t203 = t123 * MDP(22);
t127 = t173 * t141 + t170 * t142;
t202 = t127 * MDP(24);
t201 = t127 * MDP(26);
t200 = t141 * MDP(16);
t199 = t142 * MDP(13);
t198 = t148 * MDP(17);
t197 = t155 * MDP(26);
t196 = t176 * MDP(22);
t195 = t180 * MDP(18);
t166 = -t174 * pkin(2) - pkin(3);
t193 = t176 * t179 * MDP(23);
t192 = t206 * t177;
t131 = t174 * t143 - t171 * t146;
t112 = t173 * t116 - t170 * t118;
t159 = -t180 * pkin(4) + t166;
t191 = -t177 * MDP(19) + t195;
t189 = MDP(27) * t176 + MDP(28) * t179;
t129 = -t175 * pkin(3) - t131;
t188 = (t178 * MDP(6) + t181 * MDP(7)) * t172;
t187 = (MDP(24) * t179 - MDP(25) * t176) * t157;
t124 = t141 * pkin(4) + t129;
t111 = t148 * pkin(10) + t113;
t114 = t127 * pkin(5) - t128 * pkin(10) + t124;
t108 = -t176 * t111 + t179 * t114;
t109 = t179 * t111 + t176 * t114;
t186 = t108 * MDP(27) - t109 * MDP(28) + t201 - t204;
t185 = t176 * MDP(24) + t179 * MDP(25) - t189 * (t170 * pkin(4) + pkin(10));
t184 = t177 * MDP(15) + t180 * MDP(16) + (-t177 * MDP(18) - t180 * MDP(19)) * t164;
t169 = t179 ^ 2;
t168 = t176 ^ 2;
t165 = -t173 * pkin(4) - pkin(5);
t154 = t157 ^ 2;
t153 = t206 * t180;
t152 = pkin(8) * t208 + t194;
t151 = -pkin(8) * t209 + t161;
t138 = t173 * t153 - t170 * t192;
t136 = t170 * t153 + t173 * t192;
t135 = t155 * pkin(5) - t157 * pkin(10) + t159;
t126 = t176 * t135 + t179 * t138;
t125 = t179 * t135 - t176 * t138;
t121 = t123 * t155;
t110 = -t148 * pkin(5) - t112;
t1 = [(t112 ^ 2 + t113 ^ 2 + t124 ^ 2) * MDP(21) + (t131 ^ 2 + t132 ^ 2 + t158 ^ 2) * MDP(12) + MDP(1) + (MDP(4) * t178 + 0.2e1 * MDP(5) * t181) * t210 + (t198 - 0.2e1 * t200) * t148 + (t201 - 0.2e1 * t204) * t127 + (0.2e1 * t188 + t207) * t175 + (-0.2e1 * t141 * MDP(14) + 0.2e1 * t148 * MDP(15) + t199) * t142 + (0.2e1 * t202 + t203 - 0.2e1 * t205) * t123 + 0.2e1 * (t151 * t175 + t167 * t214) * MDP(9) + 0.2e1 * (-pkin(1) * t210 - t152 * t175) * MDP(10) + 0.2e1 * (-t131 * t149 - t132 * t148) * MDP(11) + 0.2e1 * (t119 * t148 + t129 * t141) * MDP(18) + (-t120 * t148 + t129 * t142) * t218 + (-t109 * t127 + t110 * t123) * t215 + (-t112 * t128 - t113 * t127) * t217 + (t108 * t127 + t110 * t122) * t216; t207 + t151 * MDP(9) - t152 * MDP(10) + t177 * t199 + (-t177 * t141 + t142 * t180) * MDP(14) + (-t129 * t180 + t166 * t141) * MDP(18) + (t129 * t177 + t166 * t142) * MDP(19) + (-t138 * t127 + t136 * t128) * MDP(20) + (-t112 * t136 + t113 * t138 + t124 * t159) * MDP(21) + t121 * MDP(24) + (t136 * t122 + t125 * t127) * MDP(27) + (t136 * t123 - t126 * t127) * MDP(28) + t188 + (-t113 * MDP(20) + t186) * t155 + t184 * t148 + (-t112 * MDP(20) + (-t123 * MDP(23) - t127 * MDP(25) + t110 * MDP(27)) * t176 + (t110 * MDP(28) + t202 + t203 - t205) * t179) * t157 + ((-t148 * t171 - t149 * t174) * MDP(11) + (t131 * t174 + t132 * t171) * MDP(12)) * pkin(2); MDP(8) - 0.2e1 * t166 * t195 + (t136 ^ 2 + t138 ^ 2 + t159 ^ 2) * MDP(21) + (t171 ^ 2 + t174 ^ 2) * MDP(12) * pkin(2) ^ 2 + (t169 * MDP(22) - 0.2e1 * t193) * t154 + (MDP(13) * t177 + 0.2e1 * t180 * MDP(14) + t166 * t218) * t177 + (0.2e1 * t187 + t197) * t155 + (t136 * t157 - t138 * t155) * t217 + (t125 * t155 + t136 * t212) * t216 + (-t126 * t155 + t136 * t211) * t215; t158 * MDP(12) + (-t157 * t127 + t155 * t128) * MDP(20) + (-t112 * t155 + t113 * t157) * MDP(21) + (t155 * t122 - t127 * t212) * MDP(27) + (-t127 * t211 + t121) * MDP(28) + t191 * t148; (t136 * t155 + t138 * t157) * MDP(21); MDP(12) + (t155 ^ 2 + t154) * MDP(21); t142 * MDP(15) - t200 + t198 + t119 * MDP(18) - t120 * MDP(19) + t123 * t196 + (-t176 * t122 + t123 * t179) * MDP(23) + (-t110 * t179 + t165 * t122) * MDP(27) + (t110 * t176 + t165 * t123) * MDP(28) + t185 * t127 + ((-t127 * t170 - t128 * t173) * MDP(20) + (t112 * t173 + t113 * t170) * MDP(21)) * pkin(4); -t190 * t136 + t185 * t155 + (t179 * t196 + (-t168 + t169) * MDP(23) + t189 * t165) * t157 + ((-t155 * t170 - t157 * t173) * MDP(20) + (-t136 * t173 + t138 * t170) * MDP(21)) * pkin(4) + t184; -t219 + (-t155 * t173 + t157 * t170) * MDP(21) * pkin(4) + t191; 0.2e1 * t193 + t168 * MDP(22) + MDP(17) + (t170 ^ 2 + t173 ^ 2) * MDP(21) * pkin(4) ^ 2 - 0.2e1 * t190 * t165; t124 * MDP(21) + t190 * t127; t159 * MDP(21) + t219; 0; 0; MDP(21); t123 * MDP(24) + t186; t125 * MDP(27) - t126 * MDP(28) + t187 + t197; -t189 * t157; t185; t190; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
