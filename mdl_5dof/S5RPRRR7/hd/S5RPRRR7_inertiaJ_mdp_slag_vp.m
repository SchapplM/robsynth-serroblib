% Calculate joint inertia matrix for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRR7_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:04:10
% EndTime: 2019-12-31 19:04:11
% DurationCPUTime: 0.32s
% Computational Cost: add. (279->94), mult. (553->138), div. (0->0), fcn. (521->8), ass. (0->53)
t93 = sin(qJ(4));
t96 = cos(qJ(4));
t101 = MDP(17) * t93 + MDP(18) * t96;
t92 = sin(qJ(5));
t95 = cos(qJ(5));
t78 = t92 * t93 - t95 * t96;
t79 = t92 * t96 + t93 * t95;
t121 = pkin(7) + pkin(8);
t80 = t121 * t93;
t81 = t121 * t96;
t105 = t79 * MDP(21) - t78 * MDP(22) + (-t80 * t95 - t81 * t92) * MDP(24) - (-t80 * t92 + t81 * t95) * MDP(25);
t126 = t93 * MDP(14) + t96 * MDP(15) - t101 * pkin(7) + t105;
t99 = (MDP(24) * t95 - MDP(25) * t92) * pkin(4);
t94 = sin(qJ(3));
t71 = t79 * t94;
t72 = t78 * t94;
t114 = -t72 * MDP(21) - t71 * MDP(22);
t90 = sin(pkin(9));
t82 = pkin(1) * t90 + pkin(6);
t119 = t82 * t93;
t120 = pkin(8) * t94;
t91 = cos(pkin(9));
t83 = -pkin(1) * t91 - pkin(2);
t97 = cos(qJ(3));
t77 = -pkin(3) * t97 - pkin(7) * t94 + t83;
t73 = t96 * t77;
t59 = -t96 * t120 + t73 + (-pkin(4) - t119) * t97;
t117 = t82 * t97;
t107 = t96 * t117;
t60 = t107 + (t77 - t120) * t93;
t56 = t95 * t59 - t60 * t92;
t57 = t59 * t92 + t95 * t60;
t125 = t56 * MDP(24) - t57 * MDP(25) + t114;
t123 = -2 * MDP(20);
t122 = 0.2e1 * MDP(25);
t118 = t82 * t96;
t116 = t93 * t96;
t115 = -t71 * MDP(24) + t72 * MDP(25);
t111 = t78 * MDP(24);
t110 = t79 * MDP(19);
t109 = t94 * MDP(11);
t108 = MDP(16) + MDP(23);
t106 = MDP(13) * t116;
t104 = MDP(14) * t96 - MDP(15) * t93;
t102 = t96 * MDP(17) - t93 * MDP(18);
t88 = t96 ^ 2;
t87 = t94 ^ 2;
t86 = t93 ^ 2;
t85 = -pkin(4) * t96 - pkin(3);
t74 = (pkin(4) * t93 + t82) * t94;
t62 = t77 * t93 + t107;
t61 = -t93 * t117 + t73;
t1 = [0.2e1 * t83 * t109 + MDP(1) + (t90 ^ 2 + t91 ^ 2) * MDP(4) * pkin(1) ^ 2 + t108 * t97 ^ 2 - (-t72 * MDP(19) + t71 * t123) * t72 + (t88 * MDP(12) + MDP(5) - 0.2e1 * t106) * t87 + 0.2e1 * (-t83 * MDP(10) + (MDP(6) - t104) * t94 - t114) * t97 + 0.2e1 * (t87 * t119 - t61 * t97) * MDP(17) + 0.2e1 * (t87 * t118 + t62 * t97) * MDP(18) + 0.2e1 * (-t56 * t97 + t71 * t74) * MDP(24) + (t57 * t97 - t72 * t74) * t122; 0; MDP(4); -t72 * t110 + (-t71 * t79 + t72 * t78) * MDP(20) + (t71 * t85 + t74 * t78) * MDP(24) + (-t72 * t85 + t74 * t79) * MDP(25) + (-t82 * MDP(11) + MDP(8) - t126) * t97 + (MDP(7) - t82 * MDP(10) + MDP(12) * t116 + (-t86 + t88) * MDP(13) + (-pkin(3) * t93 - t118) * MDP(17) + (-pkin(3) * t96 + t119) * MDP(18)) * t94; -t109 + (-t79 * MDP(25) + MDP(10) + t102 - t111) * t97; 0.2e1 * t106 + 0.2e1 * t85 * t111 + MDP(12) * t86 + MDP(9) + 0.2e1 * t102 * pkin(3) + (t85 * t122 + t78 * t123 + t110) * t79; t61 * MDP(17) - t62 * MDP(18) + (-t108 - t99) * t97 + t104 * t94 + t125; -t101 * t94 + t115; t126; t108 + 0.2e1 * t99; -t97 * MDP(23) + t125; t115; t105; MDP(23) + t99; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
