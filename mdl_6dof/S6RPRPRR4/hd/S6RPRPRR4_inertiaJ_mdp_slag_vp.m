% Calculate joint inertia matrix for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR4_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:46:28
% EndTime: 2019-03-09 03:46:29
% DurationCPUTime: 0.40s
% Computational Cost: add. (393->125), mult. (690->168), div. (0->0), fcn. (634->8), ass. (0->73)
t115 = sin(pkin(10));
t106 = t115 * pkin(1) + pkin(7);
t162 = pkin(4) + t106;
t118 = sin(qJ(5));
t121 = cos(qJ(5));
t123 = -pkin(3) - pkin(8);
t117 = sin(qJ(6));
t120 = cos(qJ(6));
t100 = -t117 * t118 + t120 * t121;
t157 = -pkin(9) + t123;
t101 = t157 * t118;
t102 = t157 * t121;
t99 = t117 * t121 + t120 * t118;
t134 = t100 * MDP(25) - t99 * MDP(26) + (-t117 * t101 + t120 * t102) * MDP(28) - (t120 * t101 + t117 * t102) * MDP(29);
t161 = (MDP(21) * t123 + MDP(18)) * t121 + (-MDP(22) * t123 - MDP(19)) * t118 + t134;
t160 = -2 * MDP(24);
t159 = 0.2e1 * MDP(29);
t119 = sin(qJ(3));
t158 = t119 * pkin(5);
t122 = cos(qJ(3));
t149 = t121 * t122;
t150 = t118 * t122;
t88 = t117 * t150 - t120 * t149;
t89 = t99 * t122;
t156 = t88 * MDP(28) + t89 * MDP(29);
t155 = t100 * MDP(28) - t99 * MDP(29);
t154 = MDP(15) * pkin(3);
t97 = t162 * t119;
t153 = t118 * t97;
t116 = cos(pkin(10));
t107 = -t116 * pkin(1) - pkin(2);
t128 = -t119 * qJ(4) + t107;
t87 = t123 * t122 + t128;
t135 = pkin(9) * t122 - t87;
t75 = -t135 * t121 + t153;
t152 = t120 * t75;
t151 = t118 * t121;
t84 = t88 * MDP(26);
t148 = t89 * MDP(25);
t147 = t99 * MDP(28);
t98 = t162 * t122;
t112 = t119 ^ 2;
t114 = t122 ^ 2;
t146 = t112 + t114;
t145 = MDP(23) * t100;
t144 = qJ(4) * MDP(15);
t143 = t106 ^ 2 * MDP(15);
t142 = t118 * MDP(21);
t141 = t121 * MDP(22);
t140 = -MDP(11) + MDP(14);
t139 = MDP(20) + MDP(27);
t138 = 0.2e1 * t122;
t137 = t119 * MDP(27) - t148 + t84;
t136 = MDP(17) * t151;
t91 = t121 * t97;
t74 = t135 * t118 + t158 + t91;
t71 = -t117 * t75 + t120 * t74;
t133 = MDP(13) - t154;
t132 = MDP(10) - t133;
t131 = -t118 * MDP(18) - t121 * MDP(19);
t130 = t121 * MDP(21) - t118 * MDP(22);
t129 = t141 + t142;
t127 = (MDP(28) * t120 - MDP(29) * t117) * pkin(5);
t126 = t106 * MDP(15) + MDP(12) + t130;
t113 = t121 ^ 2;
t111 = t118 ^ 2;
t108 = t118 * pkin(5) + qJ(4);
t92 = -t122 * pkin(3) + t128;
t82 = pkin(5) * t149 + t98;
t77 = t121 * t87 + t153;
t76 = -t118 * t87 + t91;
t72 = t117 * t74 + t152;
t1 = [t92 ^ 2 * MDP(15) + MDP(1) + (t89 * MDP(23) + t88 * t160) * t89 + (t115 ^ 2 + t116 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-t107 * MDP(10) + t92 * MDP(13)) * t138 + (t111 * MDP(16) + 0.2e1 * t136 + t143) * t114 + (MDP(5) + t139 + t143) * t112 + (0.2e1 * t107 * MDP(11) - 0.2e1 * t92 * MDP(14) - 0.2e1 * t148 + 0.2e1 * t84 + (MDP(6) + t131) * t138) * t119 + 0.2e1 * (t76 * t119 + t98 * t149) * MDP(21) + 0.2e1 * (-t77 * t119 - t98 * t150) * MDP(22) + 0.2e1 * (t71 * t119 - t82 * t88) * MDP(28) + (-t72 * t119 - t82 * t89) * t159 + 0.2e1 * t146 * MDP(12) * t106; 0; t146 * MDP(15) + MDP(4); -t89 * t145 + (t100 * t88 + t89 * t99) * MDP(24) + (-t108 * t88 + t82 * t99) * MDP(28) + (t82 * t100 - t108 * t89) * MDP(29) + t129 * t98 + (MDP(8) - MDP(16) * t151 + (t111 - t113) * MDP(17) + t140 * t106 + t126 * qJ(4)) * t122 + (-pkin(3) * MDP(12) - t132 * t106 + MDP(7) + t161) * t119; t132 * t122 + (t100 * MDP(29) + t129 + t140 + t144 + t147) * t119; -0.2e1 * t136 + 0.2e1 * t108 * t147 + t113 * MDP(16) + MDP(9) + (-0.2e1 * MDP(13) + t154) * pkin(3) + (t108 * t159 + t99 * t160 + t145) * t100 + (0.2e1 * MDP(14) + 0.2e1 * t141 + 0.2e1 * t142 + t144) * qJ(4); (t126 + t155) * t119; -t122 * MDP(15); t133; MDP(15); t119 * MDP(20) + t76 * MDP(21) - t77 * MDP(22) + (t120 * t158 + t71) * MDP(28) + (-t152 + (-t74 - t158) * t117) * MDP(29) + t131 * t122 + t137; -t130 * t122 + t156; t161; t130 + t155; 0.2e1 * t127 + t139; t71 * MDP(28) - t72 * MDP(29) + t137; t156; t134; t155; MDP(27) + t127; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
