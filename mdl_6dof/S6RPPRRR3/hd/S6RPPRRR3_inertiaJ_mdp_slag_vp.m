% Calculate joint inertia matrix for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRRR3_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:23:50
% EndTime: 2019-03-09 02:23:51
% DurationCPUTime: 0.34s
% Computational Cost: add. (343->114), mult. (638->162), div. (0->0), fcn. (606->8), ass. (0->60)
t114 = cos(qJ(4));
t109 = sin(qJ(6));
t110 = sin(qJ(5));
t112 = cos(qJ(6));
t113 = cos(qJ(5));
t93 = t109 * t113 + t112 * t110;
t85 = t93 * t114;
t92 = t109 * t110 - t112 * t113;
t87 = t92 * t114;
t147 = -t87 * MDP(24) - t85 * MDP(25);
t141 = pkin(8) + pkin(9);
t94 = t141 * t110;
t95 = t141 * t113;
t122 = t93 * MDP(24) - t92 * MDP(25) + (-t109 * t95 - t112 * t94) * MDP(27) - (-t109 * t94 + t112 * t95) * MDP(28);
t145 = t110 * MDP(20) + t113 * MDP(21);
t146 = t110 * MDP(17) + t113 * MDP(18) - pkin(8) * t145 + t122;
t143 = -2 * MDP(23);
t142 = 0.2e1 * MDP(28);
t140 = pkin(9) * t114;
t111 = sin(qJ(4));
t139 = t111 * pkin(5);
t84 = t93 * t111;
t86 = t92 * t111;
t138 = -t84 * MDP(27) + t86 * MDP(28);
t137 = -t85 * MDP(27) + t87 * MDP(28);
t108 = cos(pkin(10));
t99 = -t108 * pkin(1) - pkin(2);
t96 = -pkin(7) + t99;
t136 = t110 * t96;
t134 = t113 * t96;
t124 = t111 * t134;
t107 = sin(pkin(10));
t97 = t107 * pkin(1) + qJ(3);
t89 = t111 * pkin(4) - t114 * pkin(8) + t97;
t70 = t124 + (t89 - t140) * t110;
t135 = t112 * t70;
t133 = t110 * t113;
t132 = t92 * MDP(27);
t131 = t93 * MDP(22);
t127 = t114 * MDP(14);
t126 = MDP(19) + MDP(26);
t125 = t111 * MDP(26) + t147;
t123 = MDP(16) * t133;
t83 = t113 * t89;
t69 = -t113 * t140 + t83 + (pkin(5) - t136) * t111;
t66 = -t109 * t70 + t112 * t69;
t121 = t113 * MDP(17) - t110 * MDP(18);
t120 = t113 * MDP(20) - t110 * MDP(21);
t118 = (MDP(27) * t112 - MDP(28) * t109) * pkin(5);
t117 = -t93 * MDP(28) + MDP(13) + t120 - t132;
t106 = t114 ^ 2;
t105 = t113 ^ 2;
t104 = t111 ^ 2;
t103 = t110 ^ 2;
t101 = -t113 * pkin(5) - pkin(4);
t88 = (pkin(5) * t110 - t96) * t114;
t72 = t110 * t89 + t124;
t71 = -t111 * t136 + t83;
t67 = t109 * t69 + t135;
t1 = [MDP(1) + (t97 ^ 2 + t99 ^ 2) * MDP(7) - (-t87 * MDP(22) + t85 * t143) * t87 + (t107 ^ 2 + t108 ^ 2) * MDP(4) * pkin(1) ^ 2 + t126 * t104 + (t105 * MDP(15) + MDP(8) - 0.2e1 * t123) * t106 + 0.2e1 * ((-MDP(9) + t121) * t114 + t147) * t111 + 0.2e1 * t99 * MDP(5) + 0.2e1 * (-t106 * t136 + t71 * t111) * MDP(20) + 0.2e1 * (-t106 * t134 - t72 * t111) * MDP(21) + 0.2e1 * (t66 * t111 + t88 * t85) * MDP(27) + (-t67 * t111 - t88 * t87) * t142 + 0.2e1 * (t111 * MDP(13) + MDP(6) + t127) * t97; 0; MDP(4) + MDP(7); MDP(5) + t99 * MDP(7) + (-t84 * t111 - t114 * t85) * MDP(27) + (t86 * t111 + t114 * t87) * MDP(28) + t145 * (-t104 - t106); 0; MDP(7); -t87 * t131 + (-t93 * t85 + t87 * t92) * MDP(23) + (t101 * t85 + t88 * t92) * MDP(27) + (-t101 * t87 + t88 * t93) * MDP(28) + (-t96 * MDP(14) - MDP(11) + t146) * t111 + (MDP(10) + t96 * MDP(13) + MDP(15) * t133 + (-t103 + t105) * MDP(16) + (-pkin(4) * t110 + t134) * MDP(20) + (-pkin(4) * t113 - t136) * MDP(21)) * t114; -t117 * t111 - t127; -t111 * MDP(14) + t117 * t114; 0.2e1 * t123 + 0.2e1 * t101 * t132 + t103 * MDP(15) + MDP(12) + 0.2e1 * t120 * pkin(4) + (t101 * t142 + t92 * t143 + t131) * t93; t111 * MDP(19) + t71 * MDP(20) - t72 * MDP(21) + (t112 * t139 + t66) * MDP(27) + (-t135 + (-t69 - t139) * t109) * MDP(28) + t121 * t114 + t125; -t114 * t145 + t137; -t111 * t145 + t138; t146; 0.2e1 * t118 + t126; t66 * MDP(27) - t67 * MDP(28) + t125; t137; t138; t122; MDP(26) + t118; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
