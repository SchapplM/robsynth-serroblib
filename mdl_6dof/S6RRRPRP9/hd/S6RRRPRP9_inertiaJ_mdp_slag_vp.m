% Calculate joint inertia matrix for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP9_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:27:18
% EndTime: 2019-03-09 17:27:22
% DurationCPUTime: 1.26s
% Computational Cost: add. (997->237), mult. (1724->318), div. (0->0), fcn. (1592->6), ass. (0->82)
t144 = sin(qJ(5));
t146 = sin(qJ(2));
t148 = cos(qJ(3));
t184 = t146 * t148;
t145 = sin(qJ(3));
t147 = cos(qJ(5));
t186 = t145 * t147;
t114 = t144 * t184 - t146 * t186;
t119 = t144 * t145 + t147 * t148;
t115 = t119 * t146;
t149 = cos(qJ(2));
t139 = t149 * pkin(3);
t127 = -pkin(2) * t149 - pkin(8) * t146 - pkin(1);
t185 = t145 * t149;
t181 = pkin(7) * t185 - t148 * t127;
t110 = t139 + t181;
t102 = pkin(4) * t149 - pkin(9) * t184 + t110;
t188 = pkin(7) * t148;
t112 = t145 * t127 + t149 * t188;
t109 = -qJ(4) * t149 + t112;
t105 = pkin(9) * t145 * t146 + t109;
t96 = t147 * t102 - t144 * t105;
t97 = t144 * t102 + t147 * t105;
t201 = t115 * MDP(24) - t114 * MDP(25) + t96 * MDP(27) - t97 * MDP(28);
t200 = 2 * MDP(19);
t199 = -t181 * MDP(16) - t112 * MDP(17) - t201;
t198 = 0.2e1 * t149;
t196 = -t148 * pkin(3) - t145 * qJ(4);
t195 = -2 * MDP(23);
t194 = 2 * MDP(29);
t193 = 2 * MDP(30);
t192 = -2 * MDP(31);
t150 = -pkin(3) - pkin(4);
t191 = pkin(8) - pkin(9);
t190 = pkin(3) * t145;
t189 = pkin(5) * t149;
t187 = qJ(6) * t149;
t140 = t145 ^ 2;
t142 = t148 ^ 2;
t180 = t140 + t142;
t179 = qJ(4) * MDP(20);
t177 = t115 * MDP(22);
t124 = qJ(4) * t144 - t147 * t150;
t175 = t124 * MDP(27);
t125 = qJ(4) * t147 + t144 * t150;
t174 = t125 * MDP(28);
t173 = t148 * MDP(12);
t172 = MDP(15) + MDP(26);
t171 = MDP(17) - MDP(20);
t170 = MDP(27) + MDP(29);
t169 = -MDP(28) + MDP(31);
t126 = -pkin(2) + t196;
t168 = t191 * t145;
t167 = -pkin(3) * MDP(21) - MDP(18);
t166 = -MDP(32) * pkin(5) - MDP(29);
t116 = t148 * pkin(4) - t126;
t123 = pkin(5) + t124;
t165 = t123 * MDP(32) + MDP(29);
t164 = MDP(27) - t166;
t163 = MDP(27) + t165;
t162 = MDP(32) * qJ(6) + t169;
t122 = -qJ(6) + t125;
t161 = t122 * MDP(32) - t169;
t159 = t109 * t148 + t110 * t145;
t158 = t148 * MDP(13) - t145 * MDP(14);
t156 = t145 * MDP(18) - t148 * MDP(20);
t120 = -t144 * t148 + t186;
t154 = t120 * MDP(24) - t119 * MDP(25);
t153 = -MDP(26) - t174 - t175;
t130 = qJ(4) * t184;
t106 = t130 + (t150 * t145 - pkin(7)) * t146;
t152 = -t145 * MDP(13) + t154;
t132 = pkin(8) * t185;
t128 = t191 * t148;
t113 = -t130 + (pkin(7) + t190) * t146;
t108 = t147 * t128 + t144 * t168;
t107 = t128 * t144 - t147 * t168;
t99 = t119 * pkin(5) - t120 * qJ(6) + t116;
t98 = t114 * pkin(5) - t115 * qJ(6) + t106;
t95 = -t96 - t189;
t94 = t97 + t187;
t1 = [MDP(1) + (t109 ^ 2 + t110 ^ 2 + t113 ^ 2) * MDP(21) + (t94 ^ 2 + t95 ^ 2 + t98 ^ 2) * MDP(32) + t172 * t149 ^ 2 + (t114 * t195 + t177) * t115 + (-t114 * t94 + t115 * t95) * t193 + 0.2e1 * (t114 * MDP(29) - t115 * MDP(31)) * t98 + 0.2e1 * (t114 * MDP(27) + t115 * MDP(28)) * t106 + (t110 * MDP(18) - t109 * MDP(20) - t95 * MDP(29) + t94 * MDP(31) + pkin(1) * MDP(9) - t199) * t198 + (-0.2e1 * pkin(1) * MDP(10) + (MDP(5) - t158) * t198 + (-t109 * t145 + t110 * t148) * t200 + 0.2e1 * t156 * t113 + (t142 * MDP(11) - 0.2e1 * t145 * t173 + MDP(4) + 0.2e1 * (t145 * MDP(16) + t148 * MDP(17)) * pkin(7)) * t146) * t146; t120 * t177 + t132 * MDP(16) + (-t113 * t148 + t132) * MDP(18) + t159 * MDP(19) - t113 * t145 * MDP(20) + (t159 * pkin(8) + t113 * t126) * MDP(21) + (-t114 * t120 - t115 * t119) * MDP(23) + (t106 * t119 + t116 * t114) * MDP(27) + (t106 * t120 + t116 * t115) * MDP(28) + (t99 * t114 + t98 * t119) * MDP(29) + (t107 * t115 - t108 * t114 - t119 * t94 + t120 * t95) * MDP(30) + (-t99 * t115 - t98 * t120) * MDP(31) + (t107 * t95 + t108 * t94 + t98 * t99) * MDP(32) + (-pkin(7) * MDP(10) + MDP(7) + t169 * t108 - t170 * t107 + (t171 * pkin(8) - MDP(14)) * t148 + t152) * t149 + (-pkin(7) * MDP(9) + MDP(6) + (-t140 + t142) * MDP(12) + (-pkin(2) * t145 - t188) * MDP(16) + (-pkin(2) * t148 + pkin(7) * t145) * MDP(17) + t148 * t145 * MDP(11) + t156 * t126) * t146; MDP(8) + t140 * MDP(11) + (t180 * pkin(8) ^ 2 + t126 ^ 2) * MDP(21) + (t107 ^ 2 + t108 ^ 2 + t99 ^ 2) * MDP(32) + t180 * pkin(8) * t200 + (MDP(22) * t120 + 0.2e1 * t116 * MDP(28) + t107 * t193 + t119 * t195 + t99 * t192) * t120 + 0.2e1 * (pkin(2) * MDP(16) - t126 * MDP(18)) * t148 + 0.2e1 * (t116 * MDP(27) + t99 * MDP(29) - t108 * MDP(30)) * t119 + 0.2e1 * (-pkin(2) * MDP(17) - t126 * MDP(20) + t173) * t145; (-0.2e1 * t139 - t181) * MDP(18) + t112 * MDP(20) + (-pkin(3) * t110 + qJ(4) * t109) * MDP(21) - t96 * MDP(29) + (-t122 * t114 + t115 * t123) * MDP(30) - t97 * MDP(31) + (t94 * t122 + t123 * t95) * MDP(32) + (-MDP(15) - 0.2e1 * t179 + (-pkin(5) - t123) * MDP(29) + (-qJ(6) + t122) * MDP(31) + t153) * t149 + (t196 * MDP(19) + t158) * t146 + t199; t148 * MDP(14) + (qJ(4) * t148 - t190) * MDP(19) + (-t122 * t119 + t120 * t123) * MDP(30) + t161 * t108 + t163 * t107 + ((MDP(21) * qJ(4) - t171) * t148 + (-MDP(16) + t167) * t145) * pkin(8) - t152; 0.2e1 * pkin(3) * MDP(18) + 0.2e1 * t179 + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(21) + (t122 ^ 2 + t123 ^ 2) * MDP(32) + 0.2e1 * t175 + 0.2e1 * t174 + t123 * t194 + t122 * t192 + t172; MDP(19) * t184 + t110 * MDP(21) + (-t114 * t144 - t115 * t147) * MDP(30) + (t144 * t94 - t147 * t95) * MDP(32) + (t169 * t144 + t170 * t147 + MDP(18)) * t149; (-t119 * t144 - t120 * t147) * MDP(30) + (-t107 * t147 + t108 * t144) * MDP(32) + (pkin(8) * MDP(21) + MDP(19)) * t145; t161 * t144 - t163 * t147 + t167; MDP(21) + (t144 ^ 2 + t147 ^ 2) * MDP(32); t149 * MDP(26) + (t96 + 0.2e1 * t189) * MDP(29) + (-pkin(5) * t115 - qJ(6) * t114) * MDP(30) + (t97 + 0.2e1 * t187) * MDP(31) + (-pkin(5) * t95 + qJ(6) * t94) * MDP(32) + t201; (-pkin(5) * t120 - qJ(6) * t119) * MDP(30) + t162 * t108 - t164 * t107 + t154; (-0.2e1 * pkin(5) - t124) * MDP(29) + (-0.2e1 * qJ(6) + t125) * MDP(31) + (-t123 * pkin(5) + t122 * qJ(6)) * MDP(32) + t153; t162 * t144 + t164 * t147; MDP(26) + pkin(5) * t194 + 0.2e1 * qJ(6) * MDP(31) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(32); -t149 * MDP(29) + t115 * MDP(30) + t95 * MDP(32); t120 * MDP(30) + t107 * MDP(32); t165; -t147 * MDP(32); t166; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
