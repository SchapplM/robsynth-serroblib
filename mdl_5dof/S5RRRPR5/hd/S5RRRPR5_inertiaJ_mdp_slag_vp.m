% Calculate joint inertia matrix for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR5_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:14:18
% EndTime: 2019-12-31 21:14:19
% DurationCPUTime: 0.34s
% Computational Cost: add. (551->100), mult. (1038->157), div. (0->0), fcn. (1141->8), ass. (0->53)
t105 = sin(qJ(5));
t108 = cos(qJ(5));
t130 = t105 * MDP(22) + t108 * MDP(23);
t118 = t108 * MDP(25) - t105 * MDP(26);
t110 = cos(qJ(2));
t97 = -t110 * pkin(2) - pkin(1);
t135 = 0.2e1 * t97;
t134 = pkin(6) + pkin(7);
t106 = sin(qJ(3));
t133 = pkin(2) * t106;
t103 = sin(pkin(9));
t104 = cos(pkin(9));
t109 = cos(qJ(3));
t96 = t109 * pkin(2) + pkin(3);
t84 = t103 * t96 + t104 * t133;
t107 = sin(qJ(2));
t90 = t134 * t107;
t91 = t134 * t110;
t121 = -t106 * t91 - t109 * t90;
t88 = t106 * t110 + t109 * t107;
t113 = -t88 * qJ(4) + t121;
t119 = t106 * t107 - t109 * t110;
t120 = t106 * t90 - t109 * t91;
t72 = -t119 * qJ(4) - t120;
t66 = t103 * t72 - t104 * t113;
t64 = t66 * t105;
t131 = t66 * t108;
t129 = t105 * t108;
t76 = t103 * t88 + t104 * t119;
t128 = t76 * MDP(24);
t101 = t105 ^ 2;
t124 = MDP(21) * t129;
t125 = t101 * MDP(20) + MDP(15) + 0.2e1 * t124;
t77 = -t103 * t119 + t104 * t88;
t83 = -t103 * t133 + t104 * t96;
t81 = -pkin(4) - t83;
t82 = pkin(8) + t84;
t123 = -t76 * t82 + t77 * t81;
t94 = t103 * pkin(3) + pkin(8);
t95 = -t104 * pkin(3) - pkin(4);
t122 = -t76 * t94 + t77 * t95;
t117 = -MDP(25) * t105 - MDP(26) * t108;
t116 = (t109 * MDP(16) - t106 * MDP(17)) * pkin(2);
t115 = (MDP(22) * t108 - MDP(23) * t105) * t77;
t102 = t108 ^ 2;
t114 = t88 * MDP(13) - t119 * MDP(14) + t121 * MDP(16) + t120 * MDP(17) + ((-t101 + t102) * MDP(21) + MDP(20) * t129) * t77 + t130 * t76;
t112 = 0.2e1 * t118;
t80 = t119 * pkin(3) + t97;
t68 = t103 * t113 + t104 * t72;
t65 = t76 * pkin(4) - t77 * pkin(8) + t80;
t63 = t105 * t65 + t108 * t68;
t62 = -t105 * t68 + t108 * t65;
t1 = [MDP(1) + t119 * MDP(16) * t135 + (t66 ^ 2 + t68 ^ 2 + t80 ^ 2) * MDP(19) + (t102 * MDP(20) - 0.2e1 * t124) * t77 ^ 2 + (MDP(11) * t88 - 0.2e1 * t119 * MDP(12) + MDP(17) * t135) * t88 + t128 * t76 + 0.2e1 * (t66 * t77 - t68 * t76) * MDP(18) + 0.2e1 * (t62 * t76 + t77 * t64) * MDP(25) + 0.2e1 * (t77 * t131 - t63 * t76) * MDP(26) + (MDP(4) * t107 + 0.2e1 * t110 * MDP(5)) * t107 + 0.2e1 * (-t107 * MDP(10) + t110 * MDP(9)) * pkin(1) + 0.2e1 * t115 * t76; t107 * MDP(6) + t110 * MDP(7) + (-t84 * t76 - t83 * t77) * MDP(18) + (-t66 * t83 + t68 * t84) * MDP(19) + (t123 * t105 - t131) * MDP(25) + (t123 * t108 + t64) * MDP(26) + (-t110 * MDP(10) - t107 * MDP(9)) * pkin(6) + t114; MDP(8) + (t83 ^ 2 + t84 ^ 2) * MDP(19) - t81 * t112 + 0.2e1 * t116 + t125; (t122 * t105 - t131) * MDP(25) + (t122 * t108 + t64) * MDP(26) + ((-t103 * t76 - t104 * t77) * MDP(18) + (t103 * t68 - t104 * t66) * MDP(19)) * pkin(3) + t114; (t103 * t84 + t104 * t83) * MDP(19) * pkin(3) + t116 + t125 - t118 * (t81 + t95); -t95 * t112 + (t103 ^ 2 + t104 ^ 2) * MDP(19) * pkin(3) ^ 2 + t125; t80 * MDP(19) + t118 * t76; 0; 0; MDP(19); t62 * MDP(25) - t63 * MDP(26) + t115 + t128; t117 * t82 + t130; t117 * t94 + t130; t118; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
