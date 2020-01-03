% Calculate joint inertia matrix for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR6_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:15:20
% EndTime: 2020-01-03 12:15:21
% DurationCPUTime: 0.24s
% Computational Cost: add. (467->92), mult. (815->113), div. (0->0), fcn. (859->8), ass. (0->59)
t117 = sin(qJ(5));
t121 = cos(qJ(5));
t118 = sin(qJ(4));
t122 = cos(qJ(4));
t120 = sin(qJ(2));
t109 = t120 * pkin(1) + pkin(7);
t119 = sin(qJ(3));
t98 = (-pkin(8) - t109) * t119;
t123 = cos(qJ(3));
t116 = t123 * pkin(8);
t99 = t123 * t109 + t116;
t138 = -t118 * t99 + t122 * t98;
t101 = t118 * t123 + t122 * t119;
t155 = t101 * pkin(9);
t74 = t138 - t155;
t135 = -t118 * t98 - t122 * t99;
t100 = t118 * t119 - t122 * t123;
t97 = t100 * pkin(9);
t75 = -t135 - t97;
t161 = (-t117 * t75 + t121 * t74) * MDP(26) + (-t117 * t74 - t121 * t75) * MDP(27);
t104 = (-pkin(7) - pkin(8)) * t119;
t105 = t123 * pkin(7) + t116;
t136 = t122 * t104 - t118 * t105;
t76 = t136 - t155;
t134 = -t118 * t104 - t122 * t105;
t77 = -t134 - t97;
t160 = (-t117 * t77 + t121 * t76) * MDP(26) + (-t117 * t76 - t121 * t77) * MDP(27);
t159 = t123 * MDP(10) + t119 * MDP(9);
t133 = t123 * MDP(12) - t119 * MDP(13);
t158 = t100 * MDP(19) + t101 * MDP(20);
t83 = t121 * t100 + t117 * t101;
t84 = -t117 * t100 + t121 * t101;
t157 = t83 * MDP(26) + t84 * MDP(27);
t156 = pkin(3) * t118;
t124 = cos(qJ(2));
t154 = t124 * pkin(1);
t152 = t84 * MDP(23) - t83 * MDP(24);
t110 = t122 * pkin(3) + pkin(4);
t106 = t121 * t110;
t148 = (-t117 * t156 + t106) * MDP(26);
t147 = (-t117 * t110 - t121 * t156) * MDP(27);
t143 = t117 * MDP(27);
t141 = t122 * MDP(19);
t139 = MDP(18) + MDP(25);
t112 = -t123 * pkin(3) - pkin(2);
t137 = t101 * MDP(16) - t100 * MDP(17) + t152;
t89 = t100 * pkin(4) + t112;
t132 = -MDP(12) * t119 - MDP(13) * t123;
t131 = (t124 * MDP(5) - t120 * MDP(6)) * pkin(1);
t130 = MDP(4) + (MDP(21) * t84 - 0.2e1 * MDP(22) * t83) * t84 + (MDP(7) * t119 + 0.2e1 * MDP(8) * t123) * t119 + (MDP(14) * t101 - 0.2e1 * MDP(15) * t100) * t101;
t129 = (t121 * MDP(26) - t143) * pkin(4);
t128 = t138 * MDP(19) + t135 * MDP(20) + t137 + t161;
t127 = t136 * MDP(19) + t134 * MDP(20) + t137 + t160;
t126 = 0.2e1 * t157;
t125 = 0.2e1 * t158;
t111 = -pkin(2) - t154;
t103 = t112 - t154;
t88 = t89 - t154;
t1 = [t103 * t125 - 0.2e1 * t111 * t133 + t88 * t126 + MDP(1) + t130 + 0.2e1 * t131; t130 + t131 + t133 * (pkin(2) - t111) + t157 * (t88 + t89) + t158 * (t103 + t112); 0.2e1 * pkin(2) * t133 + t112 * t125 + t89 * t126 + t130; t109 * t132 + t128 + t159; pkin(7) * t132 + t127 + t159; MDP(11) + 0.2e1 * (-t118 * MDP(20) + t141) * pkin(3) + 0.2e1 * t148 + 0.2e1 * t147 + t139; t128; t127; (t121 * pkin(4) + t106) * MDP(26) + (-pkin(4) - t110) * t143 + (t141 + (-MDP(26) * t117 - MDP(27) * t121 - MDP(20)) * t118) * pkin(3) + t139; 0.2e1 * t129 + t139; t152 + t161; t152 + t160; MDP(25) + t147 + t148; MDP(25) + t129; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
