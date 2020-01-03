% Calculate joint inertia matrix for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR12_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:40:37
% EndTime: 2019-12-31 21:40:40
% DurationCPUTime: 0.88s
% Computational Cost: add. (1117->216), mult. (2574->326), div. (0->0), fcn. (2753->10), ass. (0->97)
t155 = sin(pkin(5));
t211 = 0.2e1 * t155;
t210 = (MDP(21) * qJ(4));
t162 = cos(qJ(3));
t209 = t162 * MDP(12);
t154 = sin(pkin(10));
t156 = cos(pkin(10));
t158 = sin(qJ(5));
t161 = cos(qJ(5));
t137 = t154 * t158 - t161 * t156;
t138 = t154 * t161 + t156 * t158;
t195 = pkin(9) + qJ(4);
t141 = t195 * t154;
t142 = t195 * t156;
t166 = MDP(24) * t138 - MDP(25) * t137 + MDP(27) * (-t141 * t161 - t142 * t158) - MDP(28) * (-t141 * t158 + t142 * t161);
t169 = MDP(18) * t154 + MDP(19) * t156;
t208 = -t169 * qJ(4) - MDP(14) + t166;
t207 = 2 * MDP(16);
t206 = 0.2e1 * MDP(18);
t205 = 0.2e1 * MDP(19);
t204 = 2 * MDP(20);
t203 = -2 * MDP(23);
t202 = 0.2e1 * MDP(27);
t201 = 0.2e1 * MDP(28);
t160 = sin(qJ(2));
t200 = pkin(1) * t160;
t163 = cos(qJ(2));
t199 = pkin(1) * t163;
t198 = pkin(8) * t154;
t197 = pkin(8) * t162;
t159 = sin(qJ(3));
t196 = pkin(9) * t159;
t194 = MDP(17) * pkin(8);
t193 = pkin(2) * MDP(17);
t192 = pkin(3) * MDP(21);
t191 = t155 * t160;
t190 = t155 * t163;
t157 = cos(pkin(5));
t189 = t157 * MDP(8);
t188 = t160 * MDP(6);
t144 = pkin(7) * t191;
t126 = t144 + (-pkin(2) - t199) * t157;
t131 = -t157 * t162 + t159 * t191;
t132 = t157 * t159 + t162 * t191;
t111 = pkin(3) * t131 - qJ(4) * t132 + t126;
t172 = pkin(7) * t190;
t127 = t172 + (pkin(8) + t200) * t157;
t128 = (-pkin(2) * t163 - pkin(8) * t160 - pkin(1)) * t155;
t115 = t127 * t162 + t128 * t159;
t112 = -qJ(4) * t190 + t115;
t102 = t154 * t111 + t156 * t112;
t140 = -pkin(3) * t162 - qJ(4) * t159 - pkin(2);
t125 = t154 * t140 + t156 * t197;
t186 = MDP(15) * t163;
t185 = MDP(26) * t162;
t118 = t132 * t154 + t156 * t190;
t119 = t132 * t156 - t154 * t190;
t106 = t161 * t118 + t119 * t158;
t184 = t106 * MDP(25);
t107 = -t118 * t158 + t119 * t161;
t183 = t107 * MDP(24);
t114 = -t159 * t127 + t128 * t162;
t113 = pkin(3) * t190 - t114;
t182 = t113 * MDP(21);
t129 = t138 * t159;
t181 = t129 * MDP(25);
t130 = t137 * t159;
t180 = t130 * MDP(22);
t179 = t130 * MDP(24);
t178 = t131 * MDP(26);
t177 = t132 * MDP(13);
t176 = t137 * MDP(27);
t175 = t138 * MDP(22);
t174 = t154 * MDP(19);
t173 = t156 * MDP(18);
t101 = t156 * t111 - t112 * t154;
t171 = -t101 * t154 + t102 * t156;
t168 = -t173 + t174 - t192;
t167 = t118 * MDP(18) + t119 * MDP(19) + t182;
t153 = t159 ^ 2;
t151 = t155 ^ 2;
t149 = -pkin(4) * t156 - pkin(3);
t139 = (pkin(4) * t154 + pkin(8)) * t159;
t136 = t156 * t140;
t134 = t157 * t200 + t172;
t133 = t157 * t199 - t144;
t124 = -t154 * t197 + t136;
t120 = -t154 * t196 + t125;
t116 = -t156 * t196 + t136 + (-pkin(4) - t198) * t162;
t105 = t116 * t158 + t120 * t161;
t104 = t116 * t161 - t120 * t158;
t103 = pkin(4) * t118 + t113;
t100 = -pkin(9) * t118 + t102;
t99 = pkin(4) * t131 - pkin(9) * t119 + t101;
t98 = t100 * t161 + t158 * t99;
t97 = -t100 * t158 + t161 * t99;
t1 = [(t101 ^ 2 + t102 ^ 2 + t113 ^ 2) * MDP(21) + t132 ^ 2 * MDP(11) + MDP(1) + t151 * t160 ^ 2 * MDP(4) + (t188 * t211 + t189) * t157 + (t107 * MDP(22) + t106 * t203) * t107 + ((MDP(7) * t157 - t177) * t211 + (0.2e1 * MDP(5) * t160 + t186) * t151) * t163 + (-0.2e1 * t132 * MDP(12) + 0.2e1 * MDP(14) * t190 + t178 + 0.2e1 * t183 - 0.2e1 * t184) * t131 + (-t101 * t119 - t102 * t118) * t204 + (t103 * t106 + t131 * t97) * t202 + (t101 * t131 + t113 * t118) * t206 + (t103 * t107 - t131 * t98) * t201 + (-t102 * t131 + t113 * t119) * t205 + 0.2e1 * (-t134 * t157 - t151 * t200) * MDP(10) + (-t114 * t190 + t126 * t131) * t207 + 0.2e1 * (t115 * t190 + t126 * t132) * MDP(17) + 0.2e1 * (t133 * t157 + t151 * t199) * MDP(9); t189 + t133 * MDP(9) - t134 * MDP(10) + (-pkin(2) * t131 - t126 * t162) * MDP(16) + (-t101 * t162 + t124 * t131) * MDP(18) + (t102 * t162 - t125 * t131) * MDP(19) + (-t118 * t125 - t119 * t124) * MDP(20) + (t101 * t124 + t102 * t125) * MDP(21) - t107 * t180 + (t106 * t130 - t107 * t129) * MDP(23) + (-t107 * t162 - t130 * t131) * MDP(24) + (t106 * t162 - t129 * t131) * MDP(25) - t162 * t178 + (t103 * t129 + t104 * t131 + t106 * t139 - t162 * t97) * MDP(27) + (-t103 * t130 - t105 * t131 + t107 * t139 + t162 * t98) * MDP(28) + (-t193 + t209) * t132 + (t188 + (MDP(7) + (-MDP(14) + t194) * t162) * t163) * t155 + (t132 * MDP(11) - t131 * MDP(12) - MDP(13) * t190 + t126 * MDP(17) + (-t101 * t156 - t102 * t154) * MDP(20) + t169 * t113 + (MDP(16) * t190 + t167) * pkin(8)) * t159; MDP(8) + t153 * MDP(11) + (pkin(8) ^ 2 * t153 + t124 ^ 2 + t125 ^ 2) * MDP(21) - (t129 * t203 - t180) * t130 + (pkin(2) * t207 + 0.2e1 * t179 + 0.2e1 * t181 + t185) * t162 + (-t124 * t162 + t153 * t198) * t206 + (pkin(8) * t153 * t156 + t125 * t162) * t205 + (-t104 * t162 + t129 * t139) * t202 + (t105 * t162 - t130 * t139) * t201 + (-0.2e1 * t193 + 0.2e1 * t209 + (-t124 * t156 - t125 * t154) * t204) * t159; t177 - t155 * t186 + t114 * MDP(16) - t115 * MDP(17) + (-pkin(3) * t118 - t113 * t156) * MDP(18) + (-pkin(3) * t119 + t113 * t154) * MDP(19) + t171 * MDP(20) - pkin(3) * t182 + t107 * t175 + (-t106 * t138 - t107 * t137) * MDP(23) + (t103 * t137 + t106 * t149) * MDP(27) + (t103 * t138 + t107 * t149) * MDP(28) + ((-t118 * t156 + t119 * t154) * MDP(20) + t171 * MDP(21)) * qJ(4) + t208 * t131; -t130 * t175 + (-t129 * t138 + t130 * t137) * MDP(23) + (t129 * t149 + t137 * t139) * MDP(27) + (-t130 * t149 + t138 * t139) * MDP(28) + (-t194 - t208) * t162 + (MDP(13) - t169 * pkin(3) + (-MDP(16) + t168) * pkin(8)) * t159 + (MDP(20) + t210) * (-t124 * t154 + t125 * t156); 0.2e1 * t149 * t176 + MDP(15) + (0.2e1 * t173 - 0.2e1 * t174 + t192) * pkin(3) + (t137 * t203 + t149 * t201 + t175) * t138 + (t204 + t210) * (t154 ^ 2 + t156 ^ 2) * qJ(4); t106 * MDP(27) + t107 * MDP(28) + t167; t129 * MDP(27) - t130 * MDP(28) + (pkin(8) * MDP(21) + t169) * t159; t138 * MDP(28) + t168 + t176; MDP(21); t97 * MDP(27) - t98 * MDP(28) + t178 + t183 - t184; t104 * MDP(27) - t105 * MDP(28) - t179 - t181 - t185; t166; 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
