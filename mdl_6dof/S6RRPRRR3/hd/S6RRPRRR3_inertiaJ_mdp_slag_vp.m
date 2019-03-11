% Calculate joint inertia matrix for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR3_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:25:21
% EndTime: 2019-03-09 13:25:24
% DurationCPUTime: 0.80s
% Computational Cost: add. (1475->178), mult. (2788->259), div. (0->0), fcn. (3253->10), ass. (0->99)
t175 = sin(pkin(11));
t168 = pkin(2) * t175 + pkin(8);
t179 = sin(qJ(4));
t183 = cos(qJ(4));
t213 = pkin(9) + t168;
t154 = t213 * t179;
t155 = t213 * t183;
t178 = sin(qJ(5));
t182 = cos(qJ(5));
t130 = -t182 * t154 - t155 * t178;
t131 = -t154 * t178 + t155 * t182;
t161 = t178 * t179 - t182 * t183;
t162 = t178 * t183 + t179 * t182;
t119 = -pkin(10) * t162 + t130;
t120 = -pkin(10) * t161 + t131;
t177 = sin(qJ(6));
t181 = cos(qJ(6));
t137 = t181 * t161 + t162 * t177;
t138 = -t161 * t177 + t162 * t181;
t193 = t138 * MDP(29) - t137 * MDP(30) + (t119 * t181 - t120 * t177) * MDP(32) - (t119 * t177 + t120 * t181) * MDP(33);
t186 = t162 * MDP(22) - t161 * MDP(23) + t130 * MDP(25) - t131 * MDP(26) + t193;
t190 = t179 * MDP(18) + t183 * MDP(19);
t229 = MDP(15) * t179 + MDP(16) * t183 - t168 * t190 + t186;
t228 = -qJ(3) - pkin(7);
t176 = cos(pkin(11));
t180 = sin(qJ(2));
t219 = cos(qJ(2));
t159 = t175 * t219 + t176 * t180;
t125 = t162 * t159;
t126 = t161 * t159;
t114 = t181 * t125 - t126 * t177;
t115 = -t125 * t177 - t126 * t181;
t227 = t115 * MDP(29) - t114 * MDP(30);
t226 = -t126 * MDP(22) - t125 * MDP(23);
t191 = MDP(18) * t183 - MDP(19) * t179;
t150 = t161 * MDP(25);
t133 = t137 * MDP(32);
t205 = -t138 * MDP(33) - t133;
t192 = -t162 * MDP(26) - t150 + t205;
t225 = t191 + t192;
t223 = -2 * MDP(21);
t222 = 0.2e1 * MDP(26);
t221 = -2 * MDP(28);
t220 = 0.2e1 * MDP(33);
t158 = t175 * t180 - t176 * t219;
t218 = pkin(4) * t158;
t217 = pkin(4) * t178;
t216 = pkin(4) * t182;
t215 = pkin(5) * t158;
t214 = pkin(5) * t181;
t212 = MDP(11) * pkin(2);
t172 = -pkin(2) * t219 - pkin(1);
t132 = t158 * pkin(3) - t159 * pkin(8) + t172;
t164 = t228 * t180;
t165 = t228 * t219;
t141 = t164 * t175 - t165 * t176;
t117 = t183 * t132 - t141 * t179;
t207 = t159 * t183;
t110 = -pkin(9) * t207 + t117 + t218;
t209 = t141 * t183;
t113 = t209 + (-pkin(9) * t159 + t132) * t179;
t210 = t113 * t182;
t104 = t110 * t178 + t210;
t102 = -pkin(10) * t125 + t104;
t211 = t102 * t181;
t208 = t159 * t179;
t206 = t179 * t183;
t171 = pkin(5) + t216;
t166 = t181 * t171;
t145 = -t177 * t217 + t166;
t204 = MDP(32) * t145;
t146 = t171 * t177 + t181 * t217;
t203 = MDP(33) * t146;
t202 = MDP(33) * t177;
t201 = t115 * MDP(27);
t200 = t126 * MDP(20);
t199 = t182 * MDP(25);
t198 = 0.2e1 * t219;
t197 = MDP(24) + MDP(31);
t196 = t158 * MDP(31) + t227;
t169 = -pkin(2) * t176 - pkin(3);
t195 = MDP(14) * t206;
t194 = MDP(17) + t197;
t103 = t182 * t110 - t113 * t178;
t101 = pkin(10) * t126 + t103 + t215;
t98 = t181 * t101 - t102 * t177;
t139 = -t176 * t164 - t165 * t175;
t121 = pkin(4) * t208 + t139;
t99 = t101 * t177 + t211;
t163 = -pkin(4) * t183 + t169;
t189 = t158 * MDP(24) + t196 + t226;
t188 = (MDP(32) * t181 - t202) * pkin(5);
t187 = (MDP(15) * t183 - MDP(16) * t179) * t159;
t174 = t183 ^ 2;
t173 = t179 ^ 2;
t142 = pkin(5) * t161 + t163;
t118 = t132 * t179 + t209;
t116 = pkin(5) * t125 + t121;
t1 = [(t139 ^ 2 + t141 ^ 2 + t172 ^ 2) * MDP(12) + pkin(1) * MDP(9) * t198 + MDP(1) + (t174 * MDP(13) - 0.2e1 * t195) * t159 ^ 2 - (t125 * t223 - t200) * t126 + (t114 * t221 + t201) * t115 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t180 + MDP(5) * t198) * t180 + t194 * t158 ^ 2 + 0.2e1 * (t187 + t226 + t227) * t158 + 0.2e1 * (t114 * t116 + t158 * t98) * MDP(32) + 0.2e1 * (t103 * t158 + t121 * t125) * MDP(25) + (t115 * t116 - t158 * t99) * t220 + (-t104 * t158 - t121 * t126) * t222 + 0.2e1 * (t139 * t159 - t141 * t158) * MDP(11) + 0.2e1 * (t117 * t158 + t139 * t208) * MDP(18) + 0.2e1 * (-t118 * t158 + t139 * t207) * MDP(19); t180 * MDP(6) + t219 * MDP(7) + (-t125 * t162 + t126 * t161) * MDP(21) + (t121 * t161 + t125 * t163) * MDP(25) + (t121 * t162 - t126 * t163) * MDP(26) + (-t114 * t138 - t115 * t137) * MDP(28) + (t114 * t142 + t116 * t137) * MDP(32) + (t115 * t142 + t116 * t138) * MDP(33) - t162 * t200 + t138 * t201 - t191 * t139 + (-MDP(10) * t219 - t180 * MDP(9)) * pkin(7) + (-t139 * t176 + t141 * t175) * MDP(12) * pkin(2) + (-t176 * t212 + (-t173 + t174) * MDP(14) + MDP(13) * t206 + t190 * t169) * t159 + (-t175 * t212 + t229) * t158; 0.2e1 * t195 + 0.2e1 * t163 * t150 + 0.2e1 * t142 * t133 + t173 * MDP(13) + MDP(8) + (t175 ^ 2 + t176 ^ 2) * MDP(12) * pkin(2) ^ 2 - 0.2e1 * t191 * t169 + (MDP(20) * t162 + t161 * t223 + t163 * t222) * t162 + (MDP(27) * t138 + t137 * t221 + t142 * t220) * t138; MDP(12) * t172 + t225 * t158; 0; MDP(12); t158 * MDP(17) + t117 * MDP(18) - t118 * MDP(19) + (t158 * t216 + t103) * MDP(25) + (-t210 + (-t110 - t218) * t178) * MDP(26) + (t145 * t158 + t98) * MDP(32) + (-t146 * t158 - t99) * MDP(33) + t187 + t189; t229; t225; 0.2e1 * (-MDP(26) * t178 + t199) * pkin(4) + 0.2e1 * t204 - 0.2e1 * t203 + t194; t103 * MDP(25) - t104 * MDP(26) + (t158 * t214 + t98) * MDP(32) + (-t211 + (-t101 - t215) * t177) * MDP(33) + t189; t186; t192; (t166 + t214) * MDP(32) + (-pkin(5) - t171) * t202 + (t199 + (-MDP(32) * t177 - MDP(33) * t181 - MDP(26)) * t178) * pkin(4) + t197; 0.2e1 * t188 + t197; t98 * MDP(32) - t99 * MDP(33) + t196; t193; t205; MDP(31) - t203 + t204; MDP(31) + t188; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
