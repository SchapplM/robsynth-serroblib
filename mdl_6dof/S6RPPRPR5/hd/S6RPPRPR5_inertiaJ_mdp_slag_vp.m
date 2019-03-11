% Calculate joint inertia matrix for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRPR5_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:49:22
% EndTime: 2019-03-09 01:49:23
% DurationCPUTime: 0.35s
% Computational Cost: add. (334->118), mult. (570->167), div. (0->0), fcn. (514->6), ass. (0->57)
t90 = sin(pkin(9));
t91 = cos(pkin(9));
t119 = t90 ^ 2 + t91 ^ 2;
t110 = t119 * MDP(20);
t128 = qJ(5) * t110;
t126 = -2 * MDP(22);
t125 = 2 * MDP(27);
t98 = cos(qJ(4));
t124 = pkin(8) * t98;
t92 = -pkin(7) + qJ(2);
t123 = t90 * t92;
t96 = sin(qJ(4));
t122 = t92 * t96;
t93 = pkin(1) + qJ(3);
t121 = pkin(8) + qJ(5);
t95 = sin(qJ(6));
t97 = cos(qJ(6));
t78 = t95 * t90 - t97 * t91;
t73 = t78 * MDP(26);
t79 = t97 * t90 + t95 * t91;
t120 = t79 * MDP(27) + t73;
t80 = t96 * pkin(4) - t98 * qJ(5) + t93;
t68 = t91 * t122 + t90 * t80;
t88 = t96 ^ 2;
t89 = t98 ^ 2;
t118 = -t88 - t89;
t117 = pkin(4) * MDP(20);
t116 = MDP(20) * t96;
t72 = t78 * t98;
t115 = t72 * MDP(21);
t114 = t90 * MDP(18);
t113 = t91 * MDP(17);
t112 = t92 * MDP(20);
t111 = t119 * MDP(19);
t76 = t91 * t80;
t67 = -t90 * t122 + t76;
t109 = -t67 * t91 - t68 * t90;
t107 = -MDP(16) + t111;
t106 = -t113 + t114;
t105 = -t90 * MDP(17) - t91 * MDP(18);
t70 = t79 * t98;
t104 = -t72 * MDP(23) - t70 * MDP(24);
t103 = -t106 + t117;
t81 = t121 * t90;
t82 = t121 * t91;
t102 = t79 * MDP(23) - t78 * MDP(24) + (-t97 * t81 - t95 * t82) * MDP(26) - (-t95 * t81 + t97 * t82) * MDP(27);
t101 = -MDP(15) + t106 + t120;
t100 = (qJ(2) ^ 2);
t85 = -t91 * pkin(5) - pkin(4);
t77 = (pkin(5) * t90 - t92) * t98;
t71 = t78 * t96;
t69 = t79 * t96;
t64 = -t90 * t124 + t68;
t63 = -t91 * t124 + t76 + (pkin(5) - t123) * t96;
t62 = t95 * t63 + t97 * t64;
t61 = t97 * t63 - t95 * t64;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2 + t100) * MDP(6)) + (t93 ^ 2 + t100) * MDP(9) + t89 * MDP(10) + (t89 * t92 ^ 2 + t67 ^ 2 + t68 ^ 2) * MDP(20) + t88 * MDP(25) - (t70 * t126 - t115) * t72 + 0.2e1 * t93 * MDP(8) + 0.2e1 * (-t89 * t123 + t67 * t96) * MDP(17) + 0.2e1 * (-t89 * t92 * t91 - t68 * t96) * MDP(18) + 0.2e1 * (t61 * t96 + t77 * t70) * MDP(26) + (-t62 * t96 - t77 * t72) * t125 + 0.2e1 * (t93 * MDP(16) + t109 * MDP(19)) * t98 + (2 * (MDP(5) + MDP(7)) * qJ(2)) + 0.2e1 * (-t98 * MDP(11) + t93 * MDP(15) + t104) * t96; t109 * MDP(20) - (pkin(1) * MDP(6)) - t93 * MDP(9) + t101 * t96 + t107 * t98 + MDP(4) - MDP(8); MDP(6) + MDP(9) + t110; MDP(7) + qJ(2) * MDP(9) + t89 * t112 + (-t69 * t96 - t98 * t70) * MDP(26) + (t71 * t96 + t98 * t72) * MDP(27) + (t118 * MDP(18) + t68 * t116) * t91 + (t118 * MDP(17) - t67 * t116) * t90; 0; MDP(9) + (t119 * t88 + t89) * MDP(20); -t79 * t115 + (-t79 * t70 + t72 * t78) * MDP(22) + (t85 * t70 + t77 * t78) * MDP(26) + (-t85 * t72 + t77 * t79) * MDP(27) + (MDP(12) + t105 * pkin(4) + (MDP(15) + t103) * t92) * t98 + (-t92 * MDP(16) + t105 * qJ(5) - MDP(13) + t102) * t96 + (qJ(5) * MDP(20) + MDP(19)) * (-t67 * t90 + t68 * t91); 0; (t107 + t128) * t96 + (-t101 + t117) * t98; 0.2e1 * t85 * t73 + MDP(14) + (0.2e1 * t113 - 0.2e1 * t114 + t117) * pkin(4) + (MDP(21) * t79 + t85 * t125 + t78 * t126) * t79 + (0.2e1 * t111 + t128) * qJ(5); t70 * MDP(26) - t72 * MDP(27) + (-t105 - t112) * t98; 0; -t98 * MDP(20); -t103 + t120; MDP(20); t96 * MDP(25) + t61 * MDP(26) - t62 * MDP(27) + t104; t120; -t69 * MDP(26) + t71 * MDP(27); t102; 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
