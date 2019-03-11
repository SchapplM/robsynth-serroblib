% Calculate joint inertia matrix for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6PRPRRR6_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:48:48
% EndTime: 2019-03-08 20:48:49
% DurationCPUTime: 0.44s
% Computational Cost: add. (353->134), mult. (740->195), div. (0->0), fcn. (758->10), ass. (0->71)
t118 = sin(qJ(6));
t119 = sin(qJ(5));
t122 = cos(qJ(6));
t123 = cos(qJ(5));
t104 = t118 * t123 + t122 * t119;
t124 = cos(qJ(4));
t94 = t104 * t124;
t103 = t118 * t119 - t122 * t123;
t96 = t103 * t124;
t162 = -t96 * MDP(24) - t94 * MDP(25);
t157 = pkin(9) + pkin(10);
t106 = t157 * t119;
t107 = t157 * t123;
t133 = t104 * MDP(24) - t103 * MDP(25) + (-t122 * t106 - t118 * t107) * MDP(27) - (-t118 * t106 + t122 * t107) * MDP(28);
t160 = t119 * MDP(20) + t123 * MDP(21);
t161 = t119 * MDP(17) + t123 * MDP(18) - pkin(9) * t160 + t133;
t159 = -2 * MDP(23);
t158 = 0.2e1 * MDP(28);
t156 = (pkin(2) * MDP(7));
t120 = sin(qJ(4));
t155 = t120 * pkin(5);
t116 = sin(pkin(6));
t121 = sin(qJ(2));
t150 = t116 * t121;
t117 = cos(pkin(6));
t125 = cos(qJ(2));
t149 = t116 * t125;
t98 = t117 * t124 - t120 * t149;
t81 = -t98 * t119 + t123 * t150;
t82 = t119 * t150 + t98 * t123;
t76 = -t118 * t82 + t122 * t81;
t77 = t118 * t81 + t122 * t82;
t154 = t76 * MDP(27) - t77 * MDP(28);
t93 = t104 * t120;
t95 = t103 * t120;
t153 = -t93 * MDP(27) + t95 * MDP(28);
t105 = t120 * pkin(4) - t124 * pkin(9) + qJ(3);
t126 = -pkin(2) - pkin(8);
t145 = t123 * t126;
t136 = t120 * t145;
t80 = t136 + (-pkin(10) * t124 + t105) * t119;
t152 = t122 * t80;
t151 = qJ(3) * MDP(7);
t148 = t119 * t123;
t147 = t119 * t126;
t146 = t123 * t124;
t143 = t103 * MDP(27);
t142 = t104 * MDP(22);
t139 = t124 * MDP(14);
t138 = MDP(19) + MDP(26);
t137 = t120 * MDP(26) + t162;
t135 = MDP(16) * t148;
t134 = MDP(5) - t156;
t101 = t123 * t105;
t79 = -pkin(10) * t146 + t101 + (pkin(5) - t147) * t120;
t72 = -t118 * t80 + t122 * t79;
t132 = t123 * MDP(17) - t119 * MDP(18);
t131 = t123 * MDP(20) - t119 * MDP(21);
t129 = (MDP(27) * t122 - MDP(28) * t118) * pkin(5);
t128 = -t104 * MDP(28) + MDP(13) + t131 - t143;
t115 = t124 ^ 2;
t114 = t123 ^ 2;
t113 = t120 ^ 2;
t112 = t119 ^ 2;
t109 = -t123 * pkin(5) - pkin(4);
t102 = (pkin(5) * t119 - t126) * t124;
t97 = t117 * t120 + t124 * t149;
t88 = t119 * t105 + t136;
t87 = -t120 * t147 + t101;
t73 = t118 * t79 + t152;
t1 = [MDP(1) + (t117 ^ 2 + (t121 ^ 2 + t125 ^ 2) * t116 ^ 2) * MDP(7); (t97 * t119 * t124 + t81 * t120) * MDP(20) + (-t82 * t120 + t97 * t146) * MDP(21) + (t76 * t120 + t97 * t94) * MDP(27) + (-t77 * t120 - t97 * t96) * MDP(28) + ((MDP(3) - t134) * t125 + (t120 * MDP(13) - MDP(4) + MDP(6) + t139 + t151) * t121) * t116; MDP(2) - (-t96 * MDP(22) + t94 * t159) * t96 + t138 * t113 + ((-2 * MDP(5) + t156) * pkin(2)) + (0.2e1 * MDP(6) + 0.2e1 * t139 + t151) * qJ(3) + (t114 * MDP(15) + MDP(8) - 0.2e1 * t135) * t115 + 0.2e1 * (qJ(3) * MDP(13) + (-MDP(9) + t132) * t124 + t162) * t120 + 0.2e1 * (-t115 * t147 + t87 * t120) * MDP(20) + 0.2e1 * (-t115 * t145 - t88 * t120) * MDP(21) + 0.2e1 * (t102 * t94 + t72 * t120) * MDP(27) + (-t102 * t96 - t73 * t120) * t158; -MDP(7) * t149; (-t93 * t120 - t124 * t94) * MDP(27) + (t95 * t120 + t124 * t96) * MDP(28) + t134 + t160 * (-t113 - t115); MDP(7); -t98 * MDP(14) - t128 * t97; -t96 * t142 + (t96 * t103 - t104 * t94) * MDP(23) + (t102 * t103 + t109 * t94) * MDP(27) + (t102 * t104 - t109 * t96) * MDP(28) + (-t126 * MDP(14) - MDP(11) + t161) * t120 + (MDP(10) + t126 * MDP(13) + MDP(15) * t148 + (-t112 + t114) * MDP(16) + (-pkin(4) * t119 + t145) * MDP(20) + (-pkin(4) * t123 - t147) * MDP(21)) * t124; -t120 * MDP(14) + t128 * t124; 0.2e1 * t135 + 0.2e1 * t109 * t143 + t112 * MDP(15) + MDP(12) + 0.2e1 * t131 * pkin(4) + (t103 * t159 + t109 * t158 + t142) * t104; t81 * MDP(20) - t82 * MDP(21) + t154; t120 * MDP(19) + t87 * MDP(20) - t88 * MDP(21) + (t122 * t155 + t72) * MDP(27) + (-t152 + (-t79 - t155) * t118) * MDP(28) + t132 * t124 + t137; -t120 * t160 + t153; t161; 0.2e1 * t129 + t138; t154; t72 * MDP(27) - t73 * MDP(28) + t137; t153; t133; MDP(26) + t129; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
