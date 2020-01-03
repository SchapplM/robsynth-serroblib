% Calculate joint inertia matrix for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RRPPP1_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:32
% EndTime: 2019-12-31 19:24:34
% DurationCPUTime: 0.56s
% Computational Cost: add. (744->164), mult. (1659->243), div. (0->0), fcn. (1619->6), ass. (0->55)
t93 = sin(pkin(5));
t111 = qJ(3) * t93;
t96 = sin(qJ(2));
t97 = cos(qJ(2));
t83 = -t97 * pkin(2) - t96 * t111 - pkin(1);
t95 = cos(pkin(5));
t102 = qJ(3) * t95 + pkin(7);
t84 = t102 * t96;
t94 = cos(pkin(8));
t121 = (t83 * t93 - t84 * t95) * t94;
t120 = 0.2e1 * t97;
t92 = sin(pkin(8));
t119 = t92 * t93;
t118 = t92 * t95;
t117 = t93 * t94;
t116 = t93 * t97;
t115 = t94 * t95;
t114 = t95 * t97;
t113 = pkin(3) + qJ(5);
t85 = t102 * t97;
t81 = t92 * t85;
t112 = pkin(3) * t116 + t81;
t80 = pkin(2) * t118 + t94 * t111;
t77 = -t94 * t114 + t96 * t92;
t67 = t95 * t83 + t93 * t84;
t78 = t92 * t114 + t96 * t94;
t98 = -t78 * qJ(4) + t67;
t60 = t113 * t77 + t98;
t110 = t60 * MDP(22);
t62 = t77 * pkin(3) + t98;
t109 = t62 * MDP(18);
t108 = t67 * MDP(14);
t107 = MDP(15) + MDP(19);
t106 = MDP(16) - MDP(21);
t105 = MDP(18) + MDP(22);
t66 = -t84 * t118 + t83 * t119 + t94 * t85;
t104 = -pkin(2) * t94 - pkin(3);
t103 = -qJ(4) * t92 - pkin(2);
t101 = MDP(11) - t106;
t100 = MDP(12) - MDP(17) - MDP(20);
t71 = -t95 * qJ(4) - t80;
t63 = qJ(4) * t116 - t66;
t91 = t93 ^ 2;
t87 = t92 * t111;
t79 = pkin(2) * t115 - t87;
t73 = (-pkin(3) * t94 + t103) * t93;
t72 = t104 * t95 + t87;
t70 = (-t113 * t94 + t103) * t93;
t69 = pkin(4) * t117 - t71;
t68 = pkin(4) * t119 + t87 + (-qJ(5) + t104) * t95;
t65 = -t81 + t121;
t64 = t112 - t121;
t61 = -t77 * pkin(4) - t63;
t59 = t84 * t115 + t78 * pkin(4) + (qJ(5) * t97 - t83 * t94) * t93 + t112;
t1 = [MDP(1) + pkin(1) * MDP(9) * t120 + (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) * MDP(14) + (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) * MDP(18) + (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) * MDP(22) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t96 + MDP(5) * t120) * t96 + 0.2e1 * (-t65 * t116 + t67 * t77) * MDP(11) + 0.2e1 * (t66 * t116 + t67 * t78) * MDP(12) + 0.2e1 * (-t65 * t78 - t66 * t77) * MDP(13) + 0.2e1 * (t63 * t77 + t64 * t78) * MDP(15) + 0.2e1 * (-t64 * t116 - t62 * t77) * MDP(16) + 0.2e1 * (t63 * t116 - t62 * t78) * MDP(17) + 0.2e1 * (t59 * t78 - t61 * t77) * MDP(19) + 0.2e1 * (-t61 * t116 - t60 * t78) * MDP(20) + 0.2e1 * (t59 * t116 + t60 * t77) * MDP(21); t96 * MDP(6) + t97 * MDP(7) + (-t80 * t77 - t79 * t78) * MDP(13) + (t65 * t79 + t66 * t80) * MDP(14) + (t71 * t77 + t72 * t78) * MDP(15) + (t63 * t71 + t64 * t72) * MDP(18) + (t68 * t78 - t69 * t77) * MDP(19) + (t59 * t68 + t61 * t69) * MDP(22) + (-t77 * MDP(16) - t78 * MDP(17) + t109) * t73 + (-t78 * MDP(20) + t77 * MDP(21) + t110) * t70 + (t65 * MDP(11) - t66 * MDP(12) + t64 * MDP(16) - t63 * MDP(17) + t61 * MDP(20) - t59 * MDP(21)) * t95 + (-t97 * MDP(10) - t96 * MDP(9)) * pkin(7) + ((-pkin(2) * t77 - t67 * t94 - t79 * t97) * MDP(11) + (-pkin(2) * t78 + t67 * t92 + t80 * t97) * MDP(12) + (-t65 * t92 + t66 * t94) * MDP(13) - pkin(2) * t108 + (-t63 * t94 + t64 * t92) * MDP(15) + (t62 * t94 - t72 * t97) * MDP(16) + (-t62 * t92 + t71 * t97) * MDP(17) + (t59 * t92 + t61 * t94) * MDP(19) + (-t60 * t92 - t69 * t97) * MDP(20) + (-t60 * t94 + t68 * t97) * MDP(21)) * t93; MDP(8) + (t91 * pkin(2) ^ 2 + t79 ^ 2 + t80 ^ 2) * MDP(14) + (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) * MDP(18) + (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) * MDP(22) + 0.2e1 * (MDP(11) * t94 - MDP(12) * t92) * t91 * pkin(2) + 0.2e1 * (t79 * MDP(11) - t80 * MDP(12) + t72 * MDP(16) - t71 * MDP(17) + t69 * MDP(20) - t68 * MDP(21)) * t95 + 0.2e1 * ((-t79 * t92 + t80 * t94) * MDP(13) + (-t71 * t94 + t72 * t92) * MDP(15) + (t68 * t92 + t69 * t94) * MDP(19) + (t94 * MDP(16) - t92 * MDP(17)) * t73 + (-t92 * MDP(20) - t94 * MDP(21)) * t70) * t93; t100 * t78 + t101 * t77 + t108 + t109 + t110; t73 * MDP(18) + t70 * MDP(22) + (-pkin(2) * MDP(14) + t100 * t92 - t101 * t94) * t93; MDP(14) + t105; t64 * MDP(18) + t59 * MDP(22) - t106 * t116 + t107 * t78; t72 * MDP(18) + t68 * MDP(22) + t106 * t95 + t107 * t119; 0; t105; -t77 * MDP(19) - MDP(20) * t116 + t61 * MDP(22); MDP(19) * t117 + t95 * MDP(20) + t69 * MDP(22); 0; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
