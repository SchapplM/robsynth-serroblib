% Calculate joint inertia matrix for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:22
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP8_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:21:23
% EndTime: 2021-01-16 00:21:26
% DurationCPUTime: 0.60s
% Computational Cost: add. (654->158), mult. (1274->224), div. (0->0), fcn. (1238->6), ass. (0->62)
t111 = sin(qJ(3));
t114 = cos(qJ(3));
t110 = sin(qJ(4));
t113 = cos(qJ(4));
t137 = pkin(7) + pkin(8);
t96 = t137 * t111;
t97 = t137 * t114;
t79 = -t110 * t97 - t113 * t96;
t92 = t110 * t114 + t113 * t111;
t71 = -t92 * qJ(5) + t79;
t80 = -t110 * t96 + t113 * t97;
t129 = t113 * t114;
t91 = t110 * t111 - t129;
t72 = -t91 * qJ(5) + t80;
t121 = t92 * MDP(20) - t91 * MDP(21) + t79 * MDP(23) - t80 * MDP(24) + t71 * MDP(25) - t72 * MDP(26);
t145 = t121 - (t111 * MDP(16) + t114 * MDP(17)) * pkin(7) + t111 * MDP(13) + t114 * MDP(14);
t112 = sin(qJ(2));
t86 = t92 * t112;
t131 = t111 * t112;
t87 = -t110 * t131 + t112 * t129;
t144 = -t87 * MDP(20) + t86 * MDP(21);
t115 = cos(qJ(2));
t132 = pkin(8) * t112;
t135 = pkin(6) * t111;
t95 = -t115 * pkin(2) - t112 * pkin(7) - pkin(1);
t90 = t114 * t95;
t74 = -t114 * t132 + t90 + (-pkin(3) - t135) * t115;
t133 = pkin(6) * t115;
t123 = t114 * t133;
t76 = t123 + (t95 - t132) * t111;
t67 = -t110 * t76 + t113 * t74;
t117 = -t87 * qJ(5) + t67;
t143 = t67 * MDP(23) + MDP(25) * t117 - t144;
t141 = -2 * MDP(19);
t140 = 0.2e1 * MDP(24);
t139 = 0.2e1 * MDP(25);
t138 = 0.2e1 * MDP(26);
t136 = pkin(3) * t110;
t134 = pkin(6) * t114;
t130 = t111 * t114;
t127 = t87 * MDP(18);
t94 = pkin(3) * t131 + t112 * pkin(6);
t105 = t113 * pkin(3);
t102 = t105 + pkin(4);
t126 = t102 * MDP(28);
t125 = t113 * MDP(23);
t124 = MDP(15) + MDP(22);
t103 = -t114 * pkin(3) - pkin(2);
t122 = MDP(12) * t130;
t68 = t110 * t74 + t113 * t76;
t120 = t114 * MDP(13) - t111 * MDP(14);
t66 = -t86 * qJ(5) + t68;
t108 = t114 ^ 2;
t107 = t112 ^ 2;
t106 = t111 ^ 2;
t99 = t115 * t136;
t83 = t91 * pkin(4) + t103;
t82 = t111 * t95 + t123;
t81 = -t111 * t133 + t90;
t75 = t86 * pkin(4) + t94;
t65 = -t115 * pkin(4) + t117;
t1 = [MDP(1) - 0.2e1 * pkin(1) * t112 * MDP(10) + (t65 ^ 2 + t66 ^ 2 + t75 ^ 2) * MDP(28) + (t141 * t86 + t127) * t87 + t124 * t115 ^ 2 + (t108 * MDP(11) + MDP(4) - 0.2e1 * t122) * t107 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t120) * t112 + t144) * t115 + 0.2e1 * (t107 * t135 - t81 * t115) * MDP(16) + 0.2e1 * (t107 * t134 + t82 * t115) * MDP(17) + 0.2e1 * (-t67 * t115 + t94 * t86) * MDP(23) + (t68 * t115 + t94 * t87) * t140 + (-t65 * t115 + t75 * t86) * t139 + (t66 * t115 + t75 * t87) * t138 + 0.2e1 * (-t65 * t87 - t66 * t86) * MDP(27); t92 * t127 + (-t92 * t86 - t87 * t91) * MDP(19) + (t103 * t86 + t94 * t91) * MDP(23) + (t103 * t87 + t94 * t92) * MDP(24) + (t75 * t91 + t83 * t86) * MDP(25) + (t75 * t92 + t83 * t87) * MDP(26) + (-t65 * t92 - t66 * t91 - t71 * t87 - t72 * t86) * MDP(27) + (t65 * t71 + t66 * t72 + t75 * t83) * MDP(28) + (-pkin(6) * MDP(10) + MDP(7) - t145) * t115 + (MDP(6) - pkin(6) * MDP(9) + MDP(11) * t130 + (-t106 + t108) * MDP(12) + (-pkin(2) * t111 - t134) * MDP(16) + (-pkin(2) * t114 + t135) * MDP(17)) * t112; MDP(8) + t106 * MDP(11) + 0.2e1 * t122 + (t71 ^ 2 + t72 ^ 2 + t83 ^ 2) * MDP(28) + (MDP(18) * t92 - 0.2e1 * t71 * MDP(27) + t103 * t140 + t138 * t83 + t141 * t91) * t92 + 0.2e1 * (t114 * MDP(16) - t111 * MDP(17)) * pkin(2) + 0.2e1 * (t103 * MDP(23) + t83 * MDP(25) - t72 * MDP(27)) * t91; t81 * MDP(16) - t82 * MDP(17) + (-t68 + t99) * MDP(24) + (-t66 + t99) * MDP(26) + (-t102 * t87 - t136 * t86) * MDP(27) + (t65 * t102 + t136 * t66) * MDP(28) + t120 * t112 + (-pkin(3) * t125 + (-pkin(4) - t102) * MDP(25) - t124) * t115 + t143; (-t102 * t92 - t136 * t91) * MDP(27) + (t71 * t102 + t136 * t72) * MDP(28) + t145; (t139 + t126) * t102 + (0.2e1 * t125 + (MDP(28) * t136 - 0.2e1 * MDP(24) - 0.2e1 * MDP(26)) * t110) * pkin(3) + t124; -t115 * MDP(22) - t68 * MDP(24) - t66 * MDP(26) + (-0.2e1 * t115 * MDP(25) - t87 * MDP(27) + t65 * MDP(28)) * pkin(4) + t143; (-t92 * MDP(27) + t71 * MDP(28)) * pkin(4) + t121; MDP(22) + (0.2e1 * pkin(4) + t105) * MDP(25) + pkin(4) * t126 + (t125 + (-MDP(24) - MDP(26)) * t110) * pkin(3); MDP(22) + (MDP(28) * pkin(4) + t139) * pkin(4); t86 * MDP(25) + t87 * MDP(26) + t75 * MDP(28); t91 * MDP(25) + t92 * MDP(26) + t83 * MDP(28); 0; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
