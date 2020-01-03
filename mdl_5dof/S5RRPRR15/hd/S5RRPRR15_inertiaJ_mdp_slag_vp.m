% Calculate joint inertia matrix for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR15_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR15_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR15_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:10
% EndTime: 2019-12-31 20:43:11
% DurationCPUTime: 0.34s
% Computational Cost: add. (310->112), mult. (569->153), div. (0->0), fcn. (511->6), ass. (0->62)
t140 = pkin(3) + pkin(6);
t104 = sin(qJ(4));
t107 = cos(qJ(4));
t109 = -pkin(2) - pkin(7);
t103 = sin(qJ(5));
t106 = cos(qJ(5));
t86 = t103 * t107 + t106 * t104;
t87 = -t103 * t104 + t106 * t107;
t135 = -pkin(8) + t109;
t89 = t135 * t104;
t90 = t135 * t107;
t117 = t87 * MDP(24) - t86 * MDP(25) + (-t103 * t89 + t106 * t90) * MDP(27) - (t103 * t90 + t106 * t89) * MDP(28);
t139 = (MDP(20) * t109 + MDP(17)) * t107 + (-MDP(21) * t109 - MDP(18)) * t104 + t117;
t138 = 0.2e1 * MDP(27);
t137 = 0.2e1 * MDP(28);
t105 = sin(qJ(2));
t136 = t105 * pkin(4);
t134 = t87 * MDP(27) - t86 * MDP(28);
t108 = cos(qJ(2));
t93 = t140 * t108;
t133 = MDP(14) * pkin(2);
t92 = t140 * t105;
t132 = t104 * t92;
t118 = -t105 * qJ(3) - pkin(1);
t85 = t109 * t108 + t118;
t119 = pkin(8) * t108 - t85;
t69 = -t119 * t107 + t132;
t131 = t106 * t69;
t130 = t104 * t107;
t129 = t104 * t108;
t128 = t107 * t108;
t78 = t103 * t129 - t106 * t128;
t76 = t78 * MDP(25);
t79 = t86 * t108;
t127 = t79 * MDP(22);
t77 = t79 * MDP(24);
t126 = t104 * MDP(20);
t125 = t107 * MDP(21);
t124 = pkin(6) ^ 2 * MDP(14);
t123 = MDP(19) + MDP(26);
t122 = 0.2e1 * t108;
t121 = t105 * MDP(26) + t76 - t77;
t120 = MDP(16) * t130;
t88 = t107 * t92;
t68 = t119 * t104 + t136 + t88;
t65 = -t103 * t69 + t106 * t68;
t116 = MDP(12) - t133;
t115 = -t104 * MDP(17) - t107 * MDP(18);
t114 = t107 * MDP(20) - t104 * MDP(21);
t113 = (MDP(27) * t106 - MDP(28) * t103) * pkin(4);
t112 = pkin(6) * MDP(14) + MDP(11) + t114;
t102 = t108 ^ 2;
t101 = t107 ^ 2;
t100 = t105 ^ 2;
t99 = t104 ^ 2;
t95 = t104 * pkin(4) + qJ(3);
t91 = -t108 * pkin(2) + t118;
t80 = pkin(4) * t128 + t93;
t71 = t107 * t85 + t132;
t70 = -t104 * t85 + t88;
t66 = t103 * t68 + t131;
t1 = [t91 ^ 2 * MDP(14) + MDP(1) - (0.2e1 * t78 * MDP(23) - t127) * t79 + (t91 * MDP(12) + pkin(1) * MDP(9)) * t122 + (t99 * MDP(15) + 0.2e1 * t120 + t124) * t102 + (MDP(4) + t123 + t124) * t100 + (-0.2e1 * pkin(1) * MDP(10) - 0.2e1 * t91 * MDP(13) - 0.2e1 * t77 + 0.2e1 * t76 + (MDP(5) + t115) * t122) * t105 + 0.2e1 * (t70 * t105 + t93 * t128) * MDP(20) + 0.2e1 * (-t71 * t105 - t93 * t129) * MDP(21) + (t65 * t105 - t80 * t78) * t138 + (-t66 * t105 - t80 * t79) * t137 + 0.2e1 * (t100 + t102) * MDP(11) * pkin(6); -t87 * t127 + (t87 * t78 + t79 * t86) * MDP(23) + (-t95 * t78 + t80 * t86) * MDP(27) + (-t95 * t79 + t80 * t87) * MDP(28) + (t125 + t126) * t93 + (MDP(7) - MDP(15) * t130 + (-t101 + t99) * MDP(16) + (-MDP(10) + MDP(13)) * pkin(6) + t112 * qJ(3)) * t108 + (-pkin(2) * MDP(11) + MDP(6) + (-MDP(9) + t116) * pkin(6) + t139) * t105; -0.2e1 * t120 + t95 * t86 * t138 + t101 * MDP(15) + MDP(8) + (-0.2e1 * MDP(12) + t133) * pkin(2) + (MDP(22) * t87 - 0.2e1 * t86 * MDP(23) + t95 * t137) * t87 + (MDP(14) * qJ(3) + 0.2e1 * MDP(13) + 0.2e1 * t125 + 0.2e1 * t126) * qJ(3); (t112 + t134) * t105; t116; MDP(14); t105 * MDP(19) + t70 * MDP(20) - t71 * MDP(21) + (t106 * t136 + t65) * MDP(27) + (-t131 + (-t68 - t136) * t103) * MDP(28) + t115 * t108 + t121; t139; t114 + t134; 0.2e1 * t113 + t123; MDP(27) * t65 - MDP(28) * t66 + t121; t117; t134; MDP(26) + t113; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
