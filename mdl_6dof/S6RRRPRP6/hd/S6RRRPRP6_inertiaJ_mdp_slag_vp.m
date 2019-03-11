% Calculate joint inertia matrix for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRRPRP6_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:01:11
% EndTime: 2019-03-09 17:01:15
% DurationCPUTime: 0.99s
% Computational Cost: add. (1944->232), mult. (4258->351), div. (0->0), fcn. (4746->10), ass. (0->100)
t171 = sin(pkin(11));
t173 = cos(pkin(11));
t176 = sin(qJ(3));
t179 = cos(qJ(3));
t156 = t171 * t179 + t173 * t176;
t224 = 0.2e1 * t156;
t172 = sin(pkin(6));
t223 = 0.2e1 * t172;
t222 = -(MDP(16) * t176 + MDP(17) * t179) * pkin(9) + t176 * MDP(13) + t179 * MDP(14);
t221 = 0.2e1 * MDP(16);
t220 = 2 * MDP(27);
t177 = sin(qJ(2));
t219 = pkin(1) * t177;
t180 = cos(qJ(2));
t218 = pkin(1) * t180;
t217 = -qJ(4) - pkin(9);
t216 = MDP(28) * pkin(5);
t215 = qJ(6) * t156;
t174 = cos(pkin(6));
t212 = t172 * t177;
t148 = -t174 * t179 + t176 * t212;
t149 = t174 * t176 + t179 * t212;
t136 = -t148 * t171 + t149 * t173;
t175 = sin(qJ(5));
t178 = cos(qJ(5));
t211 = t172 * t180;
t131 = t136 * t178 - t175 * t211;
t214 = t131 * t178;
t159 = t217 * t179;
t193 = t217 * t176;
t142 = -t173 * t159 + t171 * t193;
t213 = t142 * t178;
t210 = t174 * MDP(8);
t209 = t175 * t178;
t208 = t177 * MDP(6);
t165 = pkin(3) * t171 + pkin(10);
t207 = qJ(6) + t165;
t195 = pkin(8) * t211;
t145 = t195 + (pkin(9) + t219) * t174;
t146 = (-pkin(2) * t180 - pkin(9) * t177 - pkin(1)) * t172;
t132 = -t145 * t176 + t179 * t146;
t123 = -pkin(3) * t211 - qJ(4) * t149 + t132;
t133 = t145 * t179 + t146 * t176;
t128 = -qJ(4) * t148 + t133;
t117 = t171 * t123 + t173 * t128;
t169 = t175 ^ 2;
t170 = t178 ^ 2;
t206 = t169 + t170;
t205 = MDP(15) * t180;
t204 = MDP(23) * t175;
t155 = t171 * t176 - t173 * t179;
t203 = MDP(24) * t155;
t130 = t136 * t175 + t178 * t211;
t202 = t130 * MDP(21);
t201 = t130 * MDP(23);
t200 = t131 * MDP(20);
t135 = t173 * t148 + t149 * t171;
t199 = t135 * MDP(22);
t198 = t135 * MDP(24);
t197 = t149 * MDP(11);
t196 = t175 * MDP(26);
t167 = -pkin(3) * t179 - pkin(2);
t166 = -pkin(3) * t173 - pkin(4);
t194 = MDP(21) * t209;
t192 = -(MDP(27) * pkin(5)) + MDP(22);
t115 = -pkin(10) * t211 + t117;
t161 = pkin(8) * t212;
t144 = t161 + (-pkin(2) - t218) * t174;
t138 = pkin(3) * t148 + t144;
t119 = pkin(4) * t135 - pkin(10) * t136 + t138;
t111 = -t115 * t175 + t178 * t119;
t116 = t123 * t173 - t171 * t128;
t139 = pkin(4) * t155 - pkin(10) * t156 + t167;
t124 = t178 * t139 - t142 * t175;
t140 = -t159 * t171 - t173 * t193;
t114 = pkin(4) * t211 - t116;
t112 = t115 * t178 + t119 * t175;
t120 = pkin(5) * t155 - t178 * t215 + t124;
t121 = t213 + (t139 - t215) * t175;
t191 = t120 * t178 + t121 * t175;
t152 = t207 * t175;
t153 = t207 * t178;
t190 = -t152 * t178 + t153 * t175;
t189 = t149 * MDP(13) - t148 * MDP(14);
t125 = t139 * t175 + t213;
t186 = t124 * MDP(25) - t125 * MDP(26);
t185 = MDP(25) * t178 - t196;
t184 = t175 * MDP(25) + t178 * MDP(26);
t183 = t175 * MDP(22) + t178 * MDP(23) - t184 * t165;
t182 = t131 * MDP(22) + t111 * MDP(25) - t112 * MDP(26) + t198 - t201;
t168 = t172 ^ 2;
t158 = -pkin(5) * t178 + t166;
t151 = t174 * t219 + t195;
t150 = t174 * t218 - t161;
t134 = pkin(5) * t156 * t175 + t140;
t129 = t175 * t130;
t113 = pkin(5) * t130 + t114;
t110 = -qJ(6) * t130 + t112;
t109 = pkin(5) * t135 - qJ(6) * t131 + t111;
t1 = [(t109 ^ 2 + t110 ^ 2 + t113 ^ 2) * MDP(28) + (t116 ^ 2 + t117 ^ 2 + t138 ^ 2) * MDP(19) + t168 * t177 ^ 2 * MDP(4) + MDP(1) + (t208 * t223 + t210) * t174 + (-0.2e1 * t148 * MDP(12) + t197) * t149 + (t198 - 0.2e1 * t201) * t135 + (0.2e1 * t199 + t200 - 0.2e1 * t202) * t131 + ((0.2e1 * MDP(5) * t177 + t205) * t168 + (MDP(7) * t174 - t189) * t223) * t180 + 0.2e1 * (t150 * t174 + t168 * t218) * MDP(9) + (-t132 * t211 + t144 * t148) * t221 + 0.2e1 * (t133 * t211 + t144 * t149) * MDP(17) + 0.2e1 * (-t151 * t174 - t168 * t219) * MDP(10) + 0.2e1 * (t111 * t135 + t114 * t130) * MDP(25) + 0.2e1 * (-t112 * t135 + t114 * t131) * MDP(26) + 0.2e1 * (-t116 * t136 - t117 * t135) * MDP(18) + (-t109 * t131 - t110 * t130) * t220; t210 + t150 * MDP(9) - t151 * MDP(10) + t176 * t197 + (-t148 * t176 + t149 * t179) * MDP(12) + (-pkin(2) * t148 - t144 * t179) * MDP(16) + (-pkin(2) * t149 + t144 * t176) * MDP(17) + (-t135 * t142 + t136 * t140) * MDP(18) + (-t116 * t140 + t117 * t142 + t138 * t167) * MDP(19) + (t124 * t135 + t130 * t140) * MDP(25) + (-t125 * t135 + t131 * t140) * MDP(26) + (-t120 * t131 - t121 * t130) * MDP(27) + (t109 * t120 + t110 * t121 + t113 * t134) * MDP(28) + (-t117 * MDP(18) + t182) * t155 + (t208 + (MDP(7) - t222) * t180) * t172 + (-t116 * MDP(18) + (-t131 * MDP(21) - t135 * MDP(23) + t114 * MDP(25) - t110 * MDP(27)) * t175 + (t114 * MDP(26) - t109 * MDP(27) + t199 + t200 - t202) * t178) * t156; MDP(8) + pkin(2) * t179 * t221 + (t140 ^ 2 + t142 ^ 2 + t167 ^ 2) * MDP(19) + (t120 ^ 2 + t121 ^ 2 + t134 ^ 2) * MDP(28) + (t170 * MDP(20) - 0.2e1 * t194) * t156 ^ 2 + (MDP(11) * t176 + 0.2e1 * t179 * MDP(12) - 0.2e1 * pkin(2) * MDP(17)) * t176 + (t203 + (MDP(22) * t178 - t204) * t224) * t155 + 0.2e1 * (-t142 * MDP(18) + t186) * t155 + (-t191 * MDP(27) + (MDP(18) + t184) * t140) * t224; -t172 * t205 + t132 * MDP(16) - t133 * MDP(17) + t175 * t200 + (-t129 + t214) * MDP(21) + (-t114 * t178 + t130 * t166) * MDP(25) + (t114 * t175 + t131 * t166) * MDP(26) + (-t109 * t175 + t110 * t178 - t130 * t153 + t131 * t152) * MDP(27) + (-t109 * t152 + t110 * t153 + t113 * t158) * MDP(28) + t183 * t135 + ((-t135 * t171 - t136 * t173) * MDP(18) + (t116 * t173 + t117 * t171) * MDP(19)) * pkin(3) + t189; (-t120 * t175 + t121 * t178) * MDP(27) + (-t120 * t152 + t121 * t153 + t134 * t158) * MDP(28) - t185 * t140 + t183 * t155 + (MDP(20) * t209 + (-t169 + t170) * MDP(21) - t190 * MDP(27) + t184 * t166) * t156 + ((-t155 * t171 - t156 * t173) * MDP(18) + (-t140 * t173 + t142 * t171) * MDP(19)) * pkin(3) + t222; MDP(15) + t169 * MDP(20) + 0.2e1 * t194 + (t152 * t175 + t153 * t178) * t220 + (t152 ^ 2 + t153 ^ 2 + t158 ^ 2) * MDP(28) + (t171 ^ 2 + t173 ^ 2) * MDP(19) * pkin(3) ^ 2 - 0.2e1 * t185 * t166; t138 * MDP(19) + (-t129 - t214) * MDP(27) + (t109 * t178 + t110 * t175) * MDP(28) + t185 * t135; -t206 * MDP(27) * t156 + t167 * MDP(19) + t191 * MDP(28) + t185 * t155; t190 * MDP(28); t206 * MDP(28) + MDP(19); (-t131 * MDP(27) + t109 * MDP(28)) * pkin(5) + t182; t120 * t216 + t203 + (t192 * t178 - t204) * t156 + t186; -t152 * t216 + (-MDP(26) * t165 + MDP(23)) * t178 + (-MDP(25) * t165 + t192) * t175; -t196 + (MDP(25) + t216) * t178; MDP(28) * (pkin(5) ^ 2) + MDP(24); t113 * MDP(28); t134 * MDP(28); t158 * MDP(28); 0; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
