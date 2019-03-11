% Calculate joint inertia matrix for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRRP2_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:01:18
% EndTime: 2019-03-09 02:01:19
% DurationCPUTime: 0.34s
% Computational Cost: add. (659->124), mult. (1150->176), div. (0->0), fcn. (1189->8), ass. (0->64)
t102 = sin(qJ(4));
t136 = cos(qJ(4));
t97 = sin(pkin(10));
t99 = cos(pkin(10));
t86 = t102 * t99 + t136 * t97;
t142 = 0.2e1 * t86;
t101 = sin(qJ(5));
t103 = cos(qJ(5));
t122 = MDP(22) - MDP(25);
t141 = MDP(14) + (MDP(21) + MDP(23)) * t103 - t122 * t101;
t100 = cos(pkin(9));
t92 = -t100 * pkin(1) - pkin(2);
t87 = -t99 * pkin(3) + t92;
t140 = 0.2e1 * t87;
t85 = t102 * t97 - t136 * t99;
t139 = pkin(8) * t85;
t138 = t85 * pkin(5);
t98 = sin(pkin(9));
t90 = t98 * pkin(1) + qJ(3);
t137 = pkin(7) + t90;
t95 = t101 ^ 2;
t135 = t95 * t86;
t76 = t85 * pkin(4) - t86 * pkin(8) + t87;
t81 = t137 * t97;
t82 = t137 * t99;
t78 = -t102 * t81 + t136 * t82;
t71 = t101 * t76 + t103 * t78;
t134 = t97 ^ 2 + t99 ^ 2;
t96 = t103 ^ 2;
t133 = t95 + t96;
t132 = t103 * t86;
t131 = t85 * qJ(6);
t130 = t92 * MDP(8);
t129 = t97 * MDP(6);
t128 = t99 * MDP(5);
t80 = t96 * t86;
t127 = (-t80 - t135) * MDP(24);
t126 = t86 * MDP(15);
t113 = -t103 * pkin(5) - t101 * qJ(6);
t88 = -pkin(4) + t113;
t125 = t88 * MDP(26);
t124 = t103 * MDP(17);
t121 = t134 * MDP(8);
t120 = t101 * t78 - t103 * t76;
t119 = t133 * MDP(26);
t118 = -MDP(26) * pkin(5) - MDP(23);
t117 = pkin(8) * MDP(26) + MDP(24);
t116 = -pkin(4) * t86 - t139;
t115 = -t86 * t88 + t139;
t114 = MDP(21) - t118;
t112 = pkin(5) * t101 - t103 * qJ(6);
t111 = MDP(26) * qJ(6) - t122;
t110 = -t120 * MDP(21) - t71 * MDP(22);
t77 = t102 * t82 + t136 * t81;
t72 = t112 * t86 + t77;
t109 = -t77 * MDP(21) - t72 * MDP(23);
t108 = t77 * MDP(22) - t72 * MDP(25);
t107 = t103 * MDP(18) - t101 * MDP(19);
t106 = -t114 * t101 + t111 * t103;
t84 = t86 ^ 2;
t83 = t85 ^ 2;
t69 = t120 - t138;
t68 = t131 + t71;
t1 = [MDP(1) + t126 * t140 + t83 * MDP(20) + (t68 ^ 2 + t69 ^ 2 + t72 ^ 2) * MDP(26) + (-0.2e1 * t128 + 0.2e1 * t129 + t130) * t92 + (t100 ^ 2 + t98 ^ 2) * MDP(4) * pkin(1) ^ 2 + (t96 * MDP(16) - 0.2e1 * t101 * t124 + MDP(9)) * t84 + (MDP(14) * t140 + (-MDP(10) + t107) * t142) * t85 + 0.2e1 * (-t69 * MDP(23) + t68 * MDP(25) + t110) * t85 + ((t69 * MDP(24) + t108) * t103 + (-t68 * MDP(24) - t109) * t101) * t142 + (0.2e1 * t134 * MDP(7) + t121 * t90) * t90; (t72 * t85 + (t69 * t101 + t68 * t103) * t86) * MDP(26); MDP(4) + t121 + (t133 * t84 + t83) * MDP(26); -t128 + t129 + t130 + t126 + t127 + (t68 * t101 - t69 * t103) * MDP(26) + t141 * t85; 0; MDP(8) + t119; t86 * MDP(11) - t85 * MDP(12) - t77 * MDP(14) - t78 * MDP(15) + (t80 - t135) * MDP(17) + t72 * t125 + (t85 * MDP(19) + t116 * MDP(22) + t115 * MDP(25) + t117 * t68 + t109) * t103 + (MDP(16) * t132 + t85 * MDP(18) + t116 * MDP(21) - t115 * MDP(23) + t117 * t69 + t108) * t101; -t127 + (pkin(8) * t119 - MDP(15)) * t86 + (t125 - t141) * t85; 0; MDP(13) + t95 * MDP(16) + (t133 * pkin(8) ^ 2 + t88 ^ 2) * MDP(26) + 0.2e1 * t133 * MDP(24) * pkin(8) + 0.2e1 * (pkin(4) * MDP(21) - t88 * MDP(23)) * t103 + 0.2e1 * (-pkin(4) * MDP(22) - t88 * MDP(25) + t124) * t101; t85 * MDP(20) + (-t120 + 0.2e1 * t138) * MDP(23) + (0.2e1 * t131 + t71) * MDP(25) + (-t69 * pkin(5) + t68 * qJ(6)) * MDP(26) + (t113 * MDP(24) + t107) * t86 + t110; t106 * t86; t111 * t101 + t114 * t103; t101 * MDP(18) + t103 * MDP(19) - t112 * MDP(24) + t106 * pkin(8); MDP(20) + 0.2e1 * pkin(5) * MDP(23) + 0.2e1 * qJ(6) * MDP(25) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(26); -t85 * MDP(23) + MDP(24) * t132 + t69 * MDP(26); t101 * t86 * MDP(26); -t103 * MDP(26); t117 * t101; t118; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
