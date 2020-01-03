% Calculate joint inertia matrix for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR13_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR13_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR13_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:33:54
% EndTime: 2019-12-31 20:33:56
% DurationCPUTime: 0.49s
% Computational Cost: add. (600->130), mult. (1241->195), div. (0->0), fcn. (1311->8), ass. (0->66)
t126 = sin(qJ(5));
t129 = cos(qJ(5));
t135 = (MDP(27) * t129 - MDP(28) * t126) * pkin(4);
t124 = sin(pkin(9));
t125 = cos(pkin(9));
t127 = sin(qJ(4));
t130 = cos(qJ(4));
t109 = t124 * t127 - t130 * t125;
t110 = t124 * t130 + t125 * t127;
t154 = pkin(7) + qJ(3);
t113 = t154 * t124;
t114 = t154 * t125;
t97 = -t130 * t113 - t114 * t127;
t87 = -pkin(8) * t110 + t97;
t98 = -t113 * t127 + t114 * t130;
t88 = -pkin(8) * t109 + t98;
t92 = t129 * t109 + t110 * t126;
t93 = -t109 * t126 + t110 * t129;
t140 = t93 * MDP(24) - t92 * MDP(25) + (-t126 * t88 + t129 * t87) * MDP(27) - (t126 * t87 + t129 * t88) * MDP(28);
t164 = t110 * MDP(17) - t109 * MDP(18) + t97 * MDP(20) - t98 * MDP(21) + t140;
t163 = (MDP(14) * qJ(3));
t128 = sin(qJ(2));
t103 = t110 * t128;
t104 = t109 * t128;
t85 = t129 * t103 - t104 * t126;
t86 = -t103 * t126 - t104 * t129;
t153 = t86 * MDP(24) - t85 * MDP(25);
t131 = cos(qJ(2));
t112 = -pkin(2) * t131 - qJ(3) * t128 - pkin(1);
t107 = t125 * t112;
t156 = pkin(6) * t124;
t94 = -pkin(7) * t125 * t128 + t107 + (-pkin(3) - t156) * t131;
t155 = pkin(6) * t131;
t100 = t124 * t112 + t125 * t155;
t151 = t124 * t128;
t96 = -pkin(7) * t151 + t100;
t81 = -t127 * t96 + t130 * t94;
t77 = -pkin(4) * t131 + pkin(8) * t104 + t81;
t82 = t127 * t94 + t130 * t96;
t80 = -pkin(8) * t103 + t82;
t72 = -t126 * t80 + t129 * t77;
t73 = t126 * t77 + t129 * t80;
t162 = t72 * MDP(27) - t73 * MDP(28) + t153;
t161 = 2 * MDP(13);
t160 = -2 * MDP(16);
t159 = 0.2e1 * MDP(21);
t158 = -2 * MDP(23);
t157 = 0.2e1 * MDP(28);
t152 = pkin(2) * MDP(14);
t148 = t92 * MDP(27);
t147 = t93 * MDP(22);
t111 = pkin(3) * t151 + t128 * pkin(6);
t145 = t109 * MDP(20);
t144 = t110 * MDP(15);
t143 = t124 * MDP(12);
t142 = t125 * MDP(11);
t141 = MDP(19) + MDP(26);
t118 = -pkin(3) * t125 - pkin(2);
t138 = t124 * MDP(11) + t125 * MDP(12);
t137 = -t104 * MDP(17) - t103 * MDP(18);
t134 = -t142 + t143 - t152;
t122 = t128 ^ 2;
t102 = pkin(4) * t109 + t118;
t99 = -t124 * t155 + t107;
t95 = pkin(4) * t103 + t111;
t1 = [MDP(1) + t122 * MDP(4) - 0.2e1 * pkin(1) * t128 * MDP(10) + (pkin(6) ^ 2 * t122 + t100 ^ 2 + t99 ^ 2) * MDP(14) + (t86 * MDP(22) + t85 * t158) * t86 + t141 * t131 ^ 2 - (-t104 * MDP(15) + t103 * t160) * t104 + 0.2e1 * (t128 * MDP(5) + pkin(1) * MDP(9) - t137 - t153) * t131 + 0.2e1 * (t122 * t156 - t131 * t99) * MDP(11) + 0.2e1 * (pkin(6) * t122 * t125 + t100 * t131) * MDP(12) + 0.2e1 * (t103 * t111 - t131 * t81) * MDP(20) + (-t104 * t111 + t131 * t82) * t159 + 0.2e1 * (-t131 * t72 + t85 * t95) * MDP(27) + (t131 * t73 + t86 * t95) * t157 + (-t100 * t124 - t125 * t99) * t128 * t161; -t104 * t144 + (-t103 * t110 + t104 * t109) * MDP(16) + (t103 * t118 + t109 * t111) * MDP(20) + (-t104 * t118 + t110 * t111) * MDP(21) + t86 * t147 + (-t85 * t93 - t86 * t92) * MDP(23) + (t102 * t85 + t92 * t95) * MDP(27) + (t102 * t86 + t95 * t93) * MDP(28) + (MDP(6) - t138 * pkin(2) + (-MDP(9) + t134) * pkin(6)) * t128 + (-pkin(6) * MDP(10) + t138 * qJ(3) + MDP(7) - t164) * t131 + (MDP(13) + t163) * (t100 * t125 - t124 * t99); 0.2e1 * t118 * t145 + 0.2e1 * t102 * t148 + MDP(8) + (0.2e1 * t142 - 0.2e1 * t143 + t152) * pkin(2) + (t102 * t157 + t92 * t158 + t147) * t93 + (t109 * t160 + t118 * t159 + t144) * t110 + (t161 + t163) * (t124 ^ 2 + t125 ^ 2) * qJ(3); MDP(20) * t103 - MDP(21) * t104 + t85 * MDP(27) + t86 * MDP(28) + (pkin(6) * MDP(14) + t138) * t128; t110 * MDP(21) + t93 * MDP(28) + t134 + t145 + t148; MDP(14); t81 * MDP(20) - t82 * MDP(21) + (-t141 - t135) * t131 + t137 + t162; t164; 0; 0.2e1 * t135 + t141; -t131 * MDP(26) + t162; t140; 0; MDP(26) + t135; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
