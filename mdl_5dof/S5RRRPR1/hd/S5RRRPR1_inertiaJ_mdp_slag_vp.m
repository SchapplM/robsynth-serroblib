% Calculate joint inertia matrix for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:50
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR1_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:48:59
% EndTime: 2021-01-15 22:49:01
% DurationCPUTime: 0.31s
% Computational Cost: add. (650->109), mult. (1190->158), div. (0->0), fcn. (1315->8), ass. (0->58)
t135 = 2 * MDP(18);
t104 = sin(pkin(9));
t105 = cos(pkin(9));
t107 = sin(qJ(3));
t108 = sin(qJ(2));
t110 = cos(qJ(3));
t111 = cos(qJ(2));
t92 = t107 * t108 - t110 * t111;
t93 = t107 * t111 + t110 * t108;
t78 = t104 * t93 + t105 * t92;
t102 = -t111 * pkin(2) - pkin(1);
t83 = t92 * pkin(3) + t102;
t134 = 0.2e1 * t78 * pkin(4) + 0.2e1 * t83;
t133 = 0.2e1 * t102;
t132 = 0.2e1 * t111;
t131 = pkin(6) + pkin(7);
t130 = pkin(2) * t107;
t129 = t104 * pkin(3);
t106 = sin(qJ(5));
t103 = t105 * pkin(3);
t99 = t103 + pkin(4);
t128 = t106 * t99;
t109 = cos(qJ(5));
t79 = -t104 * t92 + t105 * t93;
t69 = t106 * t79 + t109 * t78;
t127 = t69 * MDP(27);
t101 = t110 * pkin(2) + pkin(3);
t95 = t105 * t101;
t85 = -t104 * t130 + t95;
t84 = pkin(4) + t85;
t87 = t104 * t101 + t105 * t130;
t73 = -t106 * t87 + t109 * t84;
t126 = t73 * MDP(27);
t74 = -t106 * t84 - t109 * t87;
t125 = t74 * MDP(28);
t124 = t78 * MDP(18);
t123 = t79 * MDP(19);
t122 = t83 * MDP(21);
t96 = t109 * t99;
t121 = (-t106 * t129 + t96) * MDP(27);
t120 = (-t109 * t129 - t128) * MDP(28);
t119 = t104 * MDP(19);
t118 = t110 * MDP(16);
t117 = MDP(15) + MDP(26);
t97 = t131 * t108;
t98 = t131 * t111;
t116 = -t107 * t98 - t110 * t97;
t75 = -t93 * qJ(4) + t116;
t114 = t107 * t97 - t110 * t98;
t76 = -t92 * qJ(4) - t114;
t65 = -t104 * t76 + t105 * t75;
t61 = -t79 * pkin(8) + t65;
t66 = t104 * t75 + t105 * t76;
t62 = -t78 * pkin(8) + t66;
t70 = -t106 * t78 + t109 * t79;
t115 = t70 * MDP(24) - t69 * MDP(25) + (-t106 * t62 + t109 * t61) * MDP(27) + (-t106 * t61 - t109 * t62) * MDP(28);
t113 = t93 * MDP(13) - t92 * MDP(14) + t116 * MDP(16) + t114 * MDP(17) + t65 * MDP(18) - t66 * MDP(19) + t115;
t1 = [MDP(1) + pkin(1) * MDP(9) * t132 + t92 * MDP(16) * t133 + 0.2e1 * (-t65 * t79 - t66 * t78) * MDP(20) + (t65 ^ 2 + t66 ^ 2) * MDP(21) + t127 * t134 + (t122 + 0.2e1 * t123 + 0.2e1 * t124) * t83 + (MDP(11) * t93 - 0.2e1 * t92 * MDP(12) + MDP(17) * t133) * t93 + (MDP(22) * t70 - 0.2e1 * t69 * MDP(23) + MDP(28) * t134) * t70 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t108 + MDP(5) * t132) * t108; t111 * MDP(7) + t108 * MDP(6) + (-t87 * t78 - t85 * t79) * MDP(20) + (t65 * t85 + t66 * t87) * MDP(21) + (-t111 * MDP(10) - t108 * MDP(9)) * pkin(6) + t113; MDP(8) + (t85 ^ 2 + t87 ^ 2) * MDP(21) + 0.2e1 * (-t107 * MDP(17) + t118) * pkin(2) + t85 * t135 - 0.2e1 * t87 * MDP(19) + 0.2e1 * t126 + 0.2e1 * t125 + t117; ((-t104 * t78 - t105 * t79) * MDP(20) + (t104 * t66 + t105 * t65) * MDP(21)) * pkin(3) + t113; (t103 + t95) * MDP(18) - t101 * t119 + (t73 + t96) * MDP(27) + (t74 - t128) * MDP(28) + (t105 * t85 * MDP(21) + (t87 * MDP(21) - t106 * MDP(27) - t109 * MDP(28) - MDP(19)) * t104) * pkin(3) + (t118 + (-MDP(18) * t104 - MDP(19) * t105 - MDP(17)) * t107) * pkin(2) + t117; 0.2e1 * t121 + 0.2e1 * t120 + t117 + (t105 * t135 - 0.2e1 * t119 + (t104 ^ 2 + t105 ^ 2) * MDP(21) * pkin(3)) * pkin(3); t70 * MDP(28) + t122 + t123 + t124 + t127; 0; 0; MDP(21); t115; MDP(26) + t125 + t126; MDP(26) + t120 + t121; 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
