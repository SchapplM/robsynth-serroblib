% Calculate joint inertia matrix for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RRRRP3_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:26
% EndTime: 2019-12-31 21:49:27
% DurationCPUTime: 0.30s
% Computational Cost: add. (341->99), mult. (549->121), div. (0->0), fcn. (378->6), ass. (0->49)
t92 = sin(qJ(4));
t90 = t92 ^ 2;
t95 = cos(qJ(4));
t118 = t95 ^ 2 + t90;
t94 = sin(qJ(2));
t127 = pkin(1) * t94;
t97 = cos(qJ(2));
t82 = t97 * pkin(1) + pkin(2);
t93 = sin(qJ(3));
t96 = cos(qJ(3));
t69 = -t96 * t127 - t93 * t82;
t67 = pkin(8) - t69;
t123 = t118 * t67;
t80 = t93 * pkin(2) + pkin(8);
t130 = t118 * t80;
t129 = t95 * MDP(15);
t100 = -0.2e1 * t92 * MDP(16) + 0.2e1 * t129;
t128 = 2 * MDP(18);
t126 = t96 * pkin(2);
t120 = t93 * t127 - t96 * t82;
t72 = -t95 * pkin(4) - t92 * qJ(5) - pkin(3);
t61 = t72 + t120;
t70 = t72 - t126;
t125 = -t61 - t70;
t124 = -t61 - t72;
t122 = -t70 - t72;
t119 = pkin(8) * t118;
t117 = pkin(3) * MDP(16);
t116 = t120 * MDP(8);
t115 = t69 * MDP(9);
t114 = t93 * MDP(9);
t113 = t97 * MDP(5);
t112 = MDP(20) * t92;
t111 = t92 * MDP(19);
t110 = t92 * MDP(12) + t95 * MDP(13) + (-t92 * pkin(4) + t95 * qJ(5)) * MDP(18);
t109 = 0.2e1 * t92 * t95 * MDP(11) + t90 * MDP(10) + MDP(7);
t107 = t118 * MDP(20);
t106 = -MDP(20) * pkin(4) - MDP(17);
t105 = MDP(4) + t109;
t103 = -0.2e1 * t95 * MDP(17) - 0.2e1 * t111;
t102 = pkin(3) * t129 + t109;
t101 = (t96 * MDP(8) - t114) * pkin(2);
t99 = (MDP(20) * qJ(5) - MDP(16) + MDP(19)) * t95 + (-MDP(15) + t106) * t92;
t84 = t92 * MDP(18);
t81 = -pkin(3) - t126;
t76 = t81 * t92;
t66 = -pkin(3) + t120;
t64 = t66 * t92;
t1 = [MDP(1) - t66 * t100 + t67 ^ 2 * t107 + (t61 * MDP(20) + t103) * t61 + 0.2e1 * (-t94 * MDP(6) + t113) * pkin(1) - 0.2e1 * t116 + 0.2e1 * t115 + t123 * t128 + t105; (-t120 + t126) * MDP(8) + (t76 + t64) * MDP(16) + (t130 + t123) * MDP(18) + (t130 * t67 + t61 * t70) * MDP(20) + (-pkin(2) - t82) * t114 + t125 * t111 + (t113 + (-MDP(9) * t96 - MDP(6)) * t94) * pkin(1) + ((-t66 - t81) * MDP(15) + t125 * MDP(17)) * t95 + t105; t130 * t128 - t81 * t100 + t80 ^ 2 * t107 + (t70 * MDP(20) + t103) * t70 + 0.2e1 * t101 + t105; -t116 + t115 + t64 * MDP(16) + (t119 + t123) * MDP(18) + (pkin(8) * t123 + t61 * t72) * MDP(20) + (-t66 * MDP(15) + t124 * MDP(17)) * t95 + (t124 * MDP(19) - t117) * t92 + t102; t76 * MDP(16) + (t119 + t130) * MDP(18) + (pkin(8) * t130 + t70 * t72) * MDP(20) + (-t81 * MDP(15) + t122 * MDP(17)) * t95 + (t122 * MDP(19) - t117) * t92 + t101 + t102; t119 * t128 + pkin(8) ^ 2 * t107 + (MDP(20) * t72 + t103) * t72 + pkin(3) * t100 + t109; t99 * t67 + t110; t99 * t80 + t110; t99 * pkin(8) + t110; MDP(14) + 0.2e1 * pkin(4) * MDP(17) + 0.2e1 * qJ(5) * MDP(19) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(20); t67 * t112 + t84; t80 * t112 + t84; pkin(8) * t112 + t84; t106; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
