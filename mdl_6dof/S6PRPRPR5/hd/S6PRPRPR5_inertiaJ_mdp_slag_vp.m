% Calculate joint inertia matrix for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRPR5_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:45:04
% EndTime: 2019-03-08 19:45:05
% DurationCPUTime: 0.33s
% Computational Cost: add. (384->110), mult. (789->152), div. (0->0), fcn. (867->10), ass. (0->63)
t137 = MDP(8) * qJ(3);
t94 = cos(pkin(11));
t123 = t94 * MDP(5);
t92 = sin(pkin(11));
t124 = t92 * MDP(6);
t133 = pkin(2) * MDP(8);
t100 = cos(qJ(4));
t97 = sin(qJ(4));
t82 = t100 * t92 + t97 * t94;
t87 = -t94 * pkin(3) - pkin(2);
t108 = -t82 * qJ(5) + t87;
t81 = -t100 * t94 + t97 * t92;
t73 = t81 * pkin(4) + t108;
t136 = t73 * MDP(19) + (MDP(14) - MDP(17)) * t81 - t123 + t124 - t133;
t135 = -0.2e1 * MDP(17);
t134 = pkin(4) + pkin(9);
t96 = sin(qJ(6));
t132 = t81 * t96;
t99 = cos(qJ(6));
t131 = t81 * t99;
t93 = sin(pkin(6));
t98 = sin(qJ(2));
t130 = t93 * t98;
t129 = t96 * t99;
t128 = pkin(8) + qJ(3);
t126 = MDP(19) * pkin(4);
t101 = cos(qJ(2));
t125 = t101 * t93;
t122 = t96 * MDP(25);
t121 = t99 * MDP(26);
t120 = qJ(5) * MDP(19);
t119 = MDP(8) + MDP(19);
t117 = -MDP(15) + MDP(18);
t116 = 0.2e1 * t82;
t115 = MDP(21) * t129;
t114 = MDP(17) - t126;
t83 = t128 * t92;
t84 = t128 * t94;
t75 = t100 * t84 - t97 * t83;
t74 = t100 * t83 + t97 * t84;
t112 = -MDP(14) + t114;
t111 = MDP(22) * t96 + MDP(23) * t99;
t110 = t99 * MDP(25) - t96 * MDP(26);
t109 = t121 + t122;
t107 = MDP(16) + t110;
t106 = -t109 - t117;
t104 = t99 * MDP(22) - t96 * MDP(23) - t110 * t134;
t95 = cos(pkin(6));
t91 = t99 ^ 2;
t90 = t96 ^ 2;
t85 = t93 ^ 2 * t101 ^ 2;
t78 = t94 * t130 + t95 * t92;
t77 = -t92 * t130 + t95 * t94;
t71 = t100 * t78 + t97 * t77;
t70 = -t100 * t77 + t97 * t78;
t69 = -t81 * pkin(5) + t75;
t68 = t82 * pkin(5) + t74;
t67 = t134 * t81 + t108;
t66 = t99 * t125 - t96 * t70;
t65 = t96 * t125 + t99 * t70;
t64 = t99 * t67 + t96 * t68;
t63 = -t96 * t67 + t99 * t68;
t1 = [MDP(1) + (t77 ^ 2 + t78 ^ 2 + t85) * MDP(8) + (t70 ^ 2 + t71 ^ 2 + t85) * MDP(19); (t70 * t82 - t71 * t81) * MDP(16) + (t70 * t74 + t71 * t75) * MDP(19) + (-t71 * t131 + t65 * t82) * MDP(25) + (t71 * t132 + t66 * t82) * MDP(26) + (-t98 * MDP(4) + (t117 * t82 + MDP(3) - t136) * t101) * t93 + (MDP(7) + t137) * (-t77 * t92 + t78 * t94); MDP(2) + (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) * MDP(19) + (t87 * MDP(15) - t73 * MDP(18)) * t116 + (MDP(9) + MDP(24)) * t82 ^ 2 + (t90 * MDP(20) + 0.2e1 * t115) * t81 ^ 2 + (0.2e1 * t123 - 0.2e1 * t124 + t133) * pkin(2) + (0.2e1 * t87 * MDP(14) + t73 * t135 + (-MDP(10) + t111) * t116) * t81 + 0.2e1 * (t74 * t82 - t75 * t81) * MDP(16) + 0.2e1 * (-t69 * t131 + t63 * t82) * MDP(25) + 0.2e1 * (t69 * t132 - t64 * t82) * MDP(26) + (0.2e1 * MDP(7) + t137) * (t92 ^ 2 + t94 ^ 2) * qJ(3); -t119 * t125; t106 * t82 + t136; t119; t112 * t70 + (-t106 + t120) * t71; (t117 + t120) * t75 + t112 * t74 + t109 * t69 + (-pkin(4) * MDP(16) + MDP(11) + t104) * t82 + (-MDP(12) + MDP(20) * t129 + (-t90 + t91) * MDP(21) - t107 * qJ(5)) * t81; 0; -0.2e1 * t115 + t91 * MDP(20) + MDP(13) + (t135 + t126) * pkin(4) + (0.2e1 * MDP(18) + t120 + 0.2e1 * t121 + 0.2e1 * t122) * qJ(5); t70 * MDP(19); t74 * MDP(19) + t107 * t82; 0; t114; MDP(19); t65 * MDP(25) + t66 * MDP(26); t82 * MDP(24) + t63 * MDP(25) - t64 * MDP(26) + t111 * t81; -t109; t104; t110; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
