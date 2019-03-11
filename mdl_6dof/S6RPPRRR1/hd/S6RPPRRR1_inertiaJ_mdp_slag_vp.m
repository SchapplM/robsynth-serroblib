% Calculate joint inertia matrix for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPPRRR1_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:40
% EndTime: 2019-03-09 02:18:41
% DurationCPUTime: 0.36s
% Computational Cost: add. (598->93), mult. (1086->133), div. (0->0), fcn. (1246->10), ass. (0->60)
t110 = sin(qJ(6));
t113 = cos(qJ(6));
t120 = t113 * MDP(28) - t110 * MDP(29);
t111 = sin(qJ(5));
t141 = cos(qJ(5));
t106 = sin(pkin(11));
t108 = cos(pkin(11));
t112 = sin(qJ(4));
t114 = cos(qJ(4));
t90 = t112 * t106 - t114 * t108;
t91 = t114 * t106 + t112 * t108;
t84 = -t111 * t90 + t141 * t91;
t81 = t84 * MDP(22);
t83 = t111 * t91 + t141 * t90;
t118 = -(MDP(21) + t120) * t83 - t81;
t134 = t90 * MDP(14);
t147 = -t91 * MDP(15) + t118 - t134;
t133 = t110 * MDP(25) + t113 * MDP(26);
t109 = cos(pkin(10));
t96 = -t109 * pkin(1) - pkin(2);
t92 = -t108 * pkin(3) + t96;
t85 = t90 * pkin(4) + t92;
t145 = 0.2e1 * t85;
t144 = 0.2e1 * t92;
t107 = sin(pkin(10));
t95 = t107 * pkin(1) + qJ(3);
t142 = pkin(7) + t95;
t88 = t142 * t106;
t89 = t142 * t108;
t127 = -t112 * t89 - t114 * t88;
t72 = -t91 * pkin(8) + t127;
t123 = t112 * t88 - t114 * t89;
t73 = -t90 * pkin(8) - t123;
t68 = t111 * t73 - t141 * t72;
t65 = t68 * t110;
t140 = t68 * t113;
t139 = t96 * MDP(8);
t138 = t106 * MDP(6);
t137 = t108 * MDP(5);
t136 = t110 * t113;
t135 = t83 * MDP(27);
t132 = t106 ^ 2 + t108 ^ 2;
t104 = t110 ^ 2;
t128 = MDP(24) * t136;
t129 = t104 * MDP(23) + MDP(20) + 0.2e1 * t128;
t126 = t132 * MDP(8);
t125 = -pkin(5) * t84 - pkin(9) * t83;
t97 = t111 * pkin(4) + pkin(9);
t98 = -t141 * pkin(4) - pkin(5);
t124 = -t83 * t97 + t84 * t98;
t121 = MDP(25) * t113 - MDP(26) * t110;
t119 = -MDP(28) * t110 - MDP(29) * t113;
t105 = t113 ^ 2;
t69 = t111 * t72 + t141 * t73;
t117 = -t68 * MDP(21) - t69 * MDP(22) + ((-t104 + t105) * MDP(24) + MDP(23) * t136 + MDP(18)) * t84 + (-MDP(19) + t133) * t83;
t116 = (t141 * MDP(21) - t111 * MDP(22)) * pkin(4);
t70 = t83 * pkin(5) - t84 * pkin(9) + t85;
t64 = t110 * t70 + t113 * t69;
t63 = -t110 * t69 + t113 * t70;
t1 = [t134 * t144 + t81 * t145 + MDP(1) + (-0.2e1 * t137 + 0.2e1 * t138 + t139) * t96 + (t107 ^ 2 + t109 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-0.2e1 * t90 * MDP(10) + MDP(15) * t144 + MDP(9) * t91) * t91 + (t105 * MDP(23) + MDP(16) - 0.2e1 * t128) * t84 ^ 2 + (MDP(21) * t145 + t135 + 0.2e1 * (-MDP(17) + t121) * t84) * t83 + 0.2e1 * (t63 * t83 + t84 * t65) * MDP(28) + 0.2e1 * (t84 * t140 - t64 * t83) * MDP(29) + (0.2e1 * t132 * MDP(7) + t126 * t95) * t95; 0; MDP(4) + t126; -t137 + t138 + t139 - t147; 0; MDP(8); t91 * MDP(11) - t90 * MDP(12) + t127 * MDP(14) + t123 * MDP(15) + (t124 * t110 - t140) * MDP(28) + (t124 * t113 + t65) * MDP(29) + t117; t147; 0; -0.2e1 * t120 * t98 + MDP(13) + 0.2e1 * t116 + t129; (t125 * t110 - t140) * MDP(28) + (t125 * t113 + t65) * MDP(29) + t117; t118; 0; t116 + t129 + t120 * (pkin(5) - t98); 0.2e1 * pkin(5) * t120 + t129; t63 * MDP(28) - t64 * MDP(29) + t121 * t84 + t135; t119 * t84; t120; t119 * t97 + t133; t119 * pkin(9) + t133; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
