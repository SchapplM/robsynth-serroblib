% Calculate joint inertia matrix for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR13_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR13_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR13_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:35
% EndTime: 2019-12-31 21:46:37
% DurationCPUTime: 0.66s
% Computational Cost: add. (621->179), mult. (1417->252), div. (0->0), fcn. (1405->8), ass. (0->84)
t127 = cos(pkin(5));
t129 = sin(qJ(3));
t132 = cos(qJ(3));
t126 = sin(pkin(5));
t130 = sin(qJ(2));
t165 = t126 * t130;
t111 = t127 * t129 + t132 * t165;
t178 = 0.2e1 * t111;
t133 = cos(qJ(2));
t164 = t126 * t133;
t177 = 2 * MDP(16);
t176 = 2 * MDP(18);
t175 = 2 * MDP(19);
t174 = 2 * MDP(20);
t173 = 2 * MDP(27);
t172 = 2 * MDP(28);
t171 = pkin(3) + pkin(9);
t170 = pkin(4) + pkin(8);
t169 = pkin(1) * t130;
t168 = pkin(1) * t133;
t167 = (MDP(21) * pkin(3));
t117 = t170 * t132;
t166 = t117 * t132;
t163 = t127 * MDP(8);
t128 = sin(qJ(5));
t131 = cos(qJ(5));
t162 = t128 * t131;
t161 = t130 * MDP(6);
t119 = pkin(3) * t164;
t148 = pkin(7) * t164;
t107 = t148 + (pkin(8) + t169) * t127;
t108 = (-pkin(2) * t133 - pkin(8) * t130 - pkin(1)) * t126;
t97 = -t129 * t107 + t108 * t132;
t96 = t119 - t97;
t160 = t96 * MDP(21);
t98 = t132 * t107 + t129 * t108;
t159 = MDP(15) * t133;
t158 = MDP(21) * qJ(4);
t157 = MDP(21) * pkin(8) ^ 2;
t156 = MDP(27) * t128;
t155 = MDP(28) * t131;
t110 = -t127 * t132 + t129 * t165;
t101 = t110 * t131 + t128 * t164;
t154 = t101 * MDP(23);
t153 = t101 * MDP(25);
t102 = t110 * t128 - t131 * t164;
t152 = t102 * MDP(22);
t151 = t110 * MDP(12);
t150 = MDP(11) + MDP(26);
t149 = MDP(17) - MDP(20);
t147 = qJ(4) * t164;
t146 = MDP(23) * t162;
t145 = -qJ(4) * t129 - pkin(2);
t144 = MDP(19) - t167;
t143 = -t128 * MDP(24) - t131 * MDP(25);
t142 = t131 * MDP(27) - t128 * MDP(28);
t95 = t147 - t98;
t118 = pkin(7) * t165;
t106 = t118 + (-pkin(2) - t168) * t127;
t141 = MDP(12) + t143;
t140 = MDP(18) + t142;
t139 = -qJ(4) * t111 + t106;
t91 = pkin(4) * t111 + pkin(9) * t164 + t96;
t92 = t110 * t171 + t139;
t89 = -t128 * t92 + t131 * t91;
t90 = t128 * t91 + t131 * t92;
t138 = t102 * MDP(24) + t89 * MDP(27) - t90 * MDP(28) + t153;
t137 = (-MDP(27) * t171 + MDP(24)) * t131 + (MDP(28) * t171 - MDP(25)) * t128;
t136 = -pkin(3) * MDP(18) + MDP(13) + t137;
t125 = t132 ^ 2;
t124 = t131 ^ 2;
t123 = t129 ^ 2;
t122 = t128 ^ 2;
t121 = t126 ^ 2;
t116 = t170 * t129;
t115 = -pkin(3) * t132 + t145;
t114 = -t132 * t171 + t145;
t113 = t127 * t169 + t148;
t112 = t127 * t168 - t118;
t100 = t114 * t131 + t116 * t128;
t99 = -t114 * t128 + t116 * t131;
t94 = pkin(3) * t110 + t139;
t93 = -pkin(4) * t110 - t95;
t1 = [(t94 ^ 2 + t95 ^ 2 + t96 ^ 2) * MDP(21) + t121 * t130 ^ 2 * MDP(4) + MDP(1) + (0.2e1 * t126 * t161 + t163) * t127 + t150 * t111 ^ 2 + (MDP(24) * t178 + t152 + 0.2e1 * t154) * t102 + (0.2e1 * MDP(5) * t130 + t159) * t121 * t133 + (-t101 * t93 + t111 * t89) * t173 + (t110 * t95 + t111 * t96) * t176 + (t102 * t93 - t111 * t90) * t172 + 0.2e1 * (-t113 * t127 - t121 * t169) * MDP(10) + (t106 * t110 - t164 * t97) * t177 + (-t110 * t94 - t96 * t164) * t175 + 0.2e1 * (t106 * t111 + t98 * t164) * MDP(17) + (-t111 * t94 + t164 * t95) * t174 + 0.2e1 * (t112 * t127 + t121 * t168) * MDP(9) + (-t151 + t153) * t178 + 0.2e1 * (-t111 * MDP(13) + t110 * MDP(14) + MDP(7) * t127) * t164; t163 + t112 * MDP(9) - t113 * MDP(10) + (-t101 * t117 + t111 * t99) * MDP(27) + (-t100 * t111 + t102 * t117) * MDP(28) + (t133 * MDP(7) + t161) * t126 + (-t110 * MDP(16) - t111 * MDP(17)) * pkin(2) + (-t110 * MDP(19) - t111 * MDP(20) + t94 * MDP(21)) * t115 + (-MDP(13) * t164 - t151 + t106 * MDP(17) + t96 * MDP(18) - t94 * MDP(20) + t150 * t111 + (t111 * MDP(18) + t160 + (MDP(16) - MDP(19)) * t164) * pkin(8) + t138) * t129 + (-MDP(14) * t164 - t106 * MDP(16) - t95 * MDP(18) + t94 * MDP(19) + (-t102 * MDP(23) + t93 * MDP(27)) * t131 + (-t93 * MDP(28) - t152 - t154) * t128 + t141 * t111 + (-t110 * MDP(18) - t95 * MDP(21) + t149 * t164) * pkin(8)) * t132; pkin(2) * t132 * t177 + MDP(8) + (t115 * MDP(21) + t132 * t175) * t115 + (MDP(22) * t122 + 0.2e1 * t146 + t157) * t125 + (t150 + t157) * t123 + 0.2e1 * (-pkin(2) * MDP(17) - t115 * MDP(20) + t141 * t132) * t129 + (t129 * t99 + t131 * t166) * t173 + (-t100 * t129 - t128 * t166) * t172 + (t123 + t125) * pkin(8) * t176; -t126 * t159 + t97 * MDP(16) - t98 * MDP(17) + (0.2e1 * t119 - t97) * MDP(19) + (-0.2e1 * t147 + t98) * MDP(20) + (-pkin(3) * t96 - qJ(4) * t95) * MDP(21) + t131 * t152 + (t101 * t131 - t102 * t128) * MDP(23) + (-qJ(4) * t101 + t128 * t93) * MDP(27) + (qJ(4) * t102 + t131 * t93) * MDP(28) + (-qJ(4) * MDP(18) - MDP(14)) * t110 + t136 * t111; (t155 + t156) * t117 + t136 * t129 + (MDP(14) - MDP(22) * t162 + (t122 - t124) * MDP(23) + t140 * qJ(4)) * t132 + ((-t149 + t158) * t132 + (-MDP(16) + t144) * t129) * pkin(8); -0.2e1 * t146 + t124 * MDP(22) + MDP(15) + (-2 * MDP(19) + t167) * pkin(3) + (t174 + 0.2e1 * t155 + 0.2e1 * t156 + t158) * qJ(4); -MDP(19) * t164 + t111 * t140 + t160; (pkin(8) * MDP(21) + t140) * t129; t144; MDP(21); t111 * MDP(26) + t138; t129 * MDP(26) + t99 * MDP(27) - t100 * MDP(28) + t132 * t143; t137; t142; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
