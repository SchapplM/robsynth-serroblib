% Calculate joint inertia matrix for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR10_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:46
% EndTime: 2019-12-31 19:10:47
% DurationCPUTime: 0.34s
% Computational Cost: add. (471->105), mult. (923->152), div. (0->0), fcn. (1014->8), ass. (0->63)
t111 = sin(qJ(4));
t114 = cos(qJ(4));
t119 = -t111 * MDP(20) - t114 * MDP(21);
t110 = sin(qJ(5));
t113 = cos(qJ(5));
t94 = t110 * t111 - t113 * t114;
t95 = t110 * t114 + t113 * t111;
t141 = pkin(7) + pkin(8);
t98 = t141 * t111;
t99 = t141 * t114;
t122 = t95 * MDP(24) - t94 * MDP(25) + (-t110 * t99 - t113 * t98) * MDP(27) - (-t110 * t98 + t113 * t99) * MDP(28);
t145 = t111 * MDP(17) + t114 * MDP(18) + t119 * pkin(7) + t122;
t109 = cos(pkin(9));
t101 = -t109 * pkin(2) - pkin(1);
t144 = 0.2e1 * t101;
t143 = -2 * MDP(23);
t142 = 0.2e1 * MDP(28);
t108 = sin(pkin(9));
t112 = sin(qJ(3));
t139 = cos(qJ(3));
t92 = t112 * t108 - t139 * t109;
t140 = t92 * pkin(4);
t138 = pkin(1) * MDP(7);
t137 = pkin(6) + qJ(2);
t86 = t94 * MDP(27);
t136 = -t95 * MDP(28) - t86;
t93 = t139 * t108 + t112 * t109;
t135 = t111 * t93;
t96 = t137 * t108;
t97 = t137 * t109;
t80 = -t112 * t96 + t139 * t97;
t133 = t114 * t80;
t78 = t92 * pkin(3) - t93 * pkin(7) + t101;
t69 = t133 + (-pkin(8) * t93 + t78) * t111;
t134 = t113 * t69;
t132 = t114 * t93;
t131 = t108 * MDP(5);
t130 = t109 * MDP(4);
t129 = t111 * t114;
t75 = t95 * t93;
t73 = t75 * MDP(25);
t76 = t94 * t93;
t74 = t76 * MDP(24);
t128 = t93 * MDP(14);
t127 = t95 * MDP(22);
t125 = MDP(19) + MDP(26);
t124 = t92 * MDP(26) - t73 - t74;
t123 = MDP(16) * t129;
t70 = -t111 * t80 + t114 * t78;
t68 = -pkin(8) * t132 + t140 + t70;
t65 = -t110 * t69 + t113 * t68;
t121 = t114 * MDP(17) - t111 * MDP(18);
t120 = t114 * MDP(20) - t111 * MDP(21);
t79 = t112 * t97 + t139 * t96;
t118 = -MDP(13) - t120;
t117 = (MDP(27) * t113 - MDP(28) * t110) * pkin(4);
t107 = t114 ^ 2;
t106 = t111 ^ 2;
t103 = -t114 * pkin(4) - pkin(3);
t72 = pkin(4) * t135 + t79;
t71 = t111 * t78 + t133;
t66 = t110 * t68 + t134;
t1 = [t128 * t144 + MDP(1) + t125 * t92 ^ 2 - (-t76 * MDP(22) + t75 * t143) * t76 + (0.2e1 * t130 - 0.2e1 * t131 + t138) * pkin(1) + (t107 * MDP(15) + MDP(8) - 0.2e1 * t123) * t93 ^ 2 + (MDP(13) * t144 - 0.2e1 * t74 - 0.2e1 * t73 + 0.2e1 * (-MDP(9) + t121) * t93) * t92 + 0.2e1 * (t79 * t135 + t70 * t92) * MDP(20) + 0.2e1 * (t79 * t132 - t71 * t92) * MDP(21) + 0.2e1 * (t65 * t92 + t72 * t75) * MDP(27) + (-t66 * t92 - t72 * t76) * t142 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t108 ^ 2 + t109 ^ 2) * qJ(2); t128 - t130 + t131 - t138 + (-t118 + t136) * t92; MDP(7); -t80 * MDP(14) - t76 * t127 + (-t95 * t75 + t76 * t94) * MDP(23) + (t103 * t75 + t72 * t94) * MDP(27) + (-t103 * t76 + t72 * t95) * MDP(28) + t118 * t79 + (MDP(10) + MDP(15) * t129 + (-t106 + t107) * MDP(16) + t119 * pkin(3)) * t93 + (-MDP(11) + t145) * t92; 0; 0.2e1 * t123 + 0.2e1 * t103 * t86 + t106 * MDP(15) + MDP(12) + 0.2e1 * t120 * pkin(3) + (t103 * t142 + t94 * t143 + t127) * t95; t92 * MDP(19) + t70 * MDP(20) - t71 * MDP(21) + (t113 * t140 + t65) * MDP(27) + (-t134 + (-t68 - t140) * t110) * MDP(28) + t121 * t93 + t124; t120 + t136; t145; 0.2e1 * t117 + t125; t65 * MDP(27) - t66 * MDP(28) + t124; t136; t122; MDP(26) + t117; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
