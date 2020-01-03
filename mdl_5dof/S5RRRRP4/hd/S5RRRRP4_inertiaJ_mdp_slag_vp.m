% Calculate joint inertia matrix for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRRP4_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:14
% EndTime: 2019-12-31 21:51:15
% DurationCPUTime: 0.31s
% Computational Cost: add. (443->106), mult. (741->136), div. (0->0), fcn. (680->6), ass. (0->54)
t140 = -MDP(19) - MDP(21);
t139 = -MDP(20) + MDP(23);
t106 = sin(qJ(3));
t108 = cos(qJ(3));
t114 = t108 * MDP(12) - t106 * MDP(13);
t105 = sin(qJ(4));
t123 = t105 * t106;
t130 = cos(qJ(4));
t87 = -t130 * t108 + t123;
t120 = t130 * t106;
t88 = t105 * t108 + t120;
t103 = t105 * pkin(3);
t93 = t103 + qJ(5);
t96 = t130 * pkin(3) + pkin(4);
t138 = t108 * MDP(10) + t106 * MDP(9) + (-t93 * t87 - t96 * t88) * MDP(22);
t137 = -0.2e1 * t88;
t136 = 2 * MDP(21);
t135 = 0.2e1 * MDP(22);
t134 = 2 * MDP(23);
t133 = -pkin(8) - pkin(7);
t107 = sin(qJ(2));
t95 = t107 * pkin(1) + pkin(7);
t131 = -pkin(8) - t95;
t109 = cos(qJ(2));
t129 = t109 * pkin(1);
t104 = t108 * pkin(8);
t86 = t108 * t95 + t104;
t70 = t105 * t86 - t131 * t120;
t71 = t131 * t123 + t130 * t86;
t128 = t70 * t88 - t71 * t87;
t91 = t108 * pkin(7) + t104;
t79 = t105 * t91 - t133 * t120;
t80 = t133 * t123 + t130 * t91;
t127 = t79 * t88 - t80 * t87;
t99 = -t108 * pkin(3) - pkin(2);
t72 = t87 * pkin(4) - t88 * qJ(5) + t99;
t65 = t72 - t129;
t126 = t65 + t72;
t125 = t88 * MDP(16) - t87 * MDP(17);
t90 = t99 - t129;
t124 = t90 + t99;
t119 = pkin(4) * t136 + MDP(18);
t118 = t88 ^ 2 * MDP(14) + t87 * MDP(15) * t137 + MDP(4) + (MDP(7) * t106 + 0.2e1 * MDP(8) * t108) * t106;
t117 = t139 * t71 + t140 * t70 + t125;
t116 = t139 * t80 + t140 * t79 + t125;
t115 = MDP(23) * t137 + t87 * t136;
t113 = -t106 * MDP(12) - t108 * MDP(13);
t112 = (t109 * MDP(5) - t107 * MDP(6)) * pkin(1);
t111 = t130 * MDP(19) - t105 * MDP(20);
t110 = 0.2e1 * t87 * MDP(19) + 0.2e1 * t88 * MDP(20);
t98 = -pkin(2) - t129;
t84 = t88 * MDP(22);
t73 = (-pkin(4) * t88 - t87 * qJ(5)) * MDP(22);
t1 = [MDP(1) + t128 * t135 + (t70 ^ 2 + t71 ^ 2) * MDP(24) + t90 * t110 + (t65 * MDP(24) + t115) * t65 + t118 - 0.2e1 * t114 * t98 + 0.2e1 * t112; (t127 + t128) * MDP(22) + (t65 * t72 + t70 * t79 + t71 * t80) * MDP(24) + t112 + (t124 * MDP(20) - t126 * MDP(23)) * t88 + (t124 * MDP(19) + t126 * MDP(21)) * t87 + t118 + t114 * (pkin(2) - t98); t127 * t135 + (t79 ^ 2 + t80 ^ 2) * MDP(24) + t99 * t110 + (MDP(24) * t72 + t115) * t72 + 0.2e1 * t114 * pkin(2) + t118; (-t70 * t96 + t71 * t93) * MDP(24) + t113 * t95 + t117 + t138; (-t79 * t96 + t80 * t93) * MDP(24) + t113 * pkin(7) + t116 + t138; MDP(11) + MDP(18) + (t93 ^ 2 + t96 ^ 2) * MDP(24) + 0.2e1 * t111 * pkin(3) + t96 * t136 + t93 * t134; t73 + (-t70 * pkin(4) + t71 * qJ(5)) * MDP(24) + t117; t73 + (-t79 * pkin(4) + t80 * qJ(5)) * MDP(24) + t116; (0.2e1 * qJ(5) + t103) * MDP(23) + (t96 * pkin(4) + t93 * qJ(5)) * MDP(24) + (t130 * MDP(21) + t111) * pkin(3) + t119; qJ(5) * t134 + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(24) + t119; t70 * MDP(24) + t84; t79 * MDP(24) + t84; -t96 * MDP(24) - MDP(21); -MDP(24) * pkin(4) - MDP(21); MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
