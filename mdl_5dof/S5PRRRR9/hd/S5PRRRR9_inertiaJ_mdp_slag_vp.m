% Calculate joint inertia matrix for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR9_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:12
% EndTime: 2019-12-05 17:21:14
% DurationCPUTime: 0.41s
% Computational Cost: add. (287->108), mult. (643->165), div. (0->0), fcn. (663->10), ass. (0->60)
t102 = sin(qJ(4));
t106 = cos(qJ(4));
t101 = sin(qJ(5));
t105 = cos(qJ(5));
t87 = t101 * t102 - t105 * t106;
t88 = t101 * t106 + t105 * t102;
t133 = pkin(8) + pkin(9);
t91 = t133 * t102;
t92 = t133 * t106;
t115 = t88 * MDP(21) - t87 * MDP(22) + (-t101 * t92 - t105 * t91) * MDP(24) - (-t101 * t91 + t105 * t92) * MDP(25);
t138 = t115 - (MDP(17) * t102 + MDP(18) * t106) * pkin(8) + t102 * MDP(14) + t106 * MDP(15);
t109 = (MDP(24) * t105 - MDP(25) * t101) * pkin(4);
t103 = sin(qJ(3));
t80 = t88 * t103;
t81 = t87 * t103;
t128 = -t81 * MDP(21) - t80 * MDP(22);
t107 = cos(qJ(3));
t124 = t103 * t106;
t132 = pkin(7) * t102;
t90 = -t107 * pkin(3) - t103 * pkin(8) - pkin(2);
t86 = t106 * t90;
t68 = -pkin(9) * t124 + t86 + (-pkin(4) - t132) * t107;
t130 = pkin(7) * t107;
t117 = t106 * t130;
t71 = t117 + (-pkin(9) * t103 + t90) * t102;
t63 = -t101 * t71 + t105 * t68;
t64 = t101 * t68 + t105 * t71;
t137 = t63 * MDP(24) - t64 * MDP(25) + t128;
t135 = -2 * MDP(20);
t134 = 0.2e1 * MDP(25);
t131 = pkin(7) * t106;
t108 = cos(qJ(2));
t99 = sin(pkin(5));
t126 = t108 * t99;
t100 = cos(pkin(5));
t104 = sin(qJ(2));
t127 = t104 * t99;
t83 = t100 * t103 + t107 * t127;
t69 = -t83 * t102 - t106 * t126;
t70 = -t102 * t126 + t83 * t106;
t65 = -t101 * t70 + t105 * t69;
t66 = t101 * t69 + t105 * t70;
t129 = t65 * MDP(24) - t66 * MDP(25);
t125 = t102 * t106;
t121 = t87 * MDP(24);
t120 = t88 * MDP(19);
t119 = MDP(11) * t103;
t118 = MDP(16) + MDP(23);
t116 = MDP(13) * t125;
t114 = t106 * MDP(14) - t102 * MDP(15);
t112 = t106 * MDP(17) - t102 * MDP(18);
t97 = t106 ^ 2;
t96 = t103 ^ 2;
t95 = t102 ^ 2;
t94 = -t106 * pkin(4) - pkin(3);
t89 = (pkin(4) * t102 + pkin(7)) * t103;
t82 = -t100 * t107 + t103 * t127;
t77 = t102 * t90 + t117;
t76 = -t102 * t130 + t86;
t1 = [MDP(1); (t82 * t102 * t103 - t69 * t107) * MDP(17) + (t70 * t107 + t82 * t124) * MDP(18) + (-t65 * t107 + t82 * t80) * MDP(24) + (t66 * t107 - t82 * t81) * MDP(25) + (-t104 * MDP(4) + (MDP(10) * t107 + MDP(3) - t119) * t108) * t99; -0.2e1 * pkin(2) * t119 + MDP(2) + t118 * t107 ^ 2 - (-t81 * MDP(19) + t80 * t135) * t81 + (t97 * MDP(12) + MDP(5) - 0.2e1 * t116) * t96 + 0.2e1 * (pkin(2) * MDP(10) + (MDP(6) - t114) * t103 - t128) * t107 + 0.2e1 * (-t76 * t107 + t96 * t132) * MDP(17) + 0.2e1 * (t77 * t107 + t96 * t131) * MDP(18) + 0.2e1 * (-t63 * t107 + t89 * t80) * MDP(24) + (t64 * t107 - t89 * t81) * t134; -t83 * MDP(11) + (t88 * MDP(25) - MDP(10) - t112 + t121) * t82; -t81 * t120 + (-t88 * t80 + t81 * t87) * MDP(20) + (t94 * t80 + t89 * t87) * MDP(24) + (-t94 * t81 + t89 * t88) * MDP(25) + (-pkin(7) * MDP(11) + MDP(8) - t138) * t107 + (MDP(7) - pkin(7) * MDP(10) + MDP(12) * t125 + (-t95 + t97) * MDP(13) + (-pkin(3) * t102 - t131) * MDP(17) + (-pkin(3) * t106 + t132) * MDP(18)) * t103; 0.2e1 * t116 + 0.2e1 * t94 * t121 + t95 * MDP(12) + MDP(9) + 0.2e1 * t112 * pkin(3) + (t94 * t134 + t87 * t135 + t120) * t88; t69 * MDP(17) - t70 * MDP(18) + t129; t76 * MDP(17) - t77 * MDP(18) + (-t118 - t109) * t107 + t114 * t103 + t137; t138; 0.2e1 * t109 + t118; t129; -t107 * MDP(23) + t137; t115; MDP(23) + t109; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
