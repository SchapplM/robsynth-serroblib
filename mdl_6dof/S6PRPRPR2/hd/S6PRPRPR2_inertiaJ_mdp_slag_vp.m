% Calculate joint inertia matrix for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRPR2_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:32:55
% EndTime: 2019-03-08 19:32:56
% DurationCPUTime: 0.46s
% Computational Cost: add. (465->126), mult. (1040->199), div. (0->0), fcn. (1131->12), ass. (0->66)
t154 = qJ(5) * MDP(16);
t111 = sin(pkin(12));
t114 = cos(pkin(12));
t141 = t111 ^ 2 + t114 ^ 2;
t153 = t141 * t154;
t152 = MDP(15) + t154;
t150 = -2 * MDP(18);
t149 = 2 * MDP(23);
t118 = sin(qJ(4));
t148 = pkin(9) * t118;
t147 = pkin(9) + qJ(5);
t112 = sin(pkin(11));
t104 = t112 * pkin(2) + pkin(8);
t121 = cos(qJ(4));
t144 = t104 * t121;
t115 = cos(pkin(11));
t105 = -t115 * pkin(2) - pkin(3);
t96 = -t121 * pkin(4) - t118 * qJ(5) + t105;
t83 = t111 * t96 + t114 * t144;
t146 = pkin(4) * MDP(16);
t145 = t104 * t111;
t117 = sin(qJ(6));
t120 = cos(qJ(6));
t97 = t117 * t111 - t120 * t114;
t92 = t97 * t118;
t143 = t92 * MDP(17);
t142 = t97 * MDP(22);
t139 = t111 * MDP(14);
t138 = t114 * MDP(13);
t137 = t141 * MDP(15);
t116 = cos(pkin(6));
t113 = sin(pkin(6));
t119 = sin(qJ(2));
t122 = cos(qJ(2));
t90 = (t112 * t122 + t115 * t119) * t113;
t85 = t116 * t118 + t90 * t121;
t88 = (t112 * t119 - t115 * t122) * t113;
t78 = -t85 * t111 + t88 * t114;
t79 = t88 * t111 + t85 * t114;
t135 = -t78 * t111 + t79 * t114;
t94 = t114 * t96;
t82 = -t111 * t144 + t94;
t134 = -t82 * t111 + t83 * t114;
t98 = t120 * t111 + t117 * t114;
t91 = t98 * t118;
t133 = t92 * MDP(19) + t91 * MDP(20);
t132 = (-t117 * t79 + t120 * t78) * MDP(22) - (t117 * t78 + t120 * t79) * MDP(23);
t131 = t91 * MDP(22) - t92 * MDP(23);
t130 = t111 * MDP(13) + t114 * MDP(14);
t129 = -t138 + t139 - t146;
t128 = t104 * MDP(16) + t130;
t100 = t147 * t111;
t101 = t147 * t114;
t127 = t98 * MDP(19) - t97 * MDP(20) + (-t120 * t100 - t117 * t101) * MDP(22) - (-t117 * t100 + t120 * t101) * MDP(23);
t126 = t98 * MDP(23) + t129 + t142;
t125 = MDP(11) - t126;
t110 = t121 ^ 2;
t109 = t118 ^ 2;
t106 = -t114 * pkin(5) - pkin(4);
t95 = (pkin(5) * t111 + t104) * t118;
t84 = -t116 * t121 + t90 * t118;
t81 = -t111 * t148 + t83;
t80 = -t114 * t148 + t94 + (-pkin(5) - t145) * t121;
t77 = t117 * t80 + t120 * t81;
t76 = -t117 * t81 + t120 * t80;
t1 = [MDP(1) + (t116 ^ 2 + t88 ^ 2 + t90 ^ 2) * MDP(5) + (t78 ^ 2 + t79 ^ 2 + t84 ^ 2) * MDP(16); (t78 * t82 + t79 * t83) * MDP(16) + t131 * t84 + (t122 * MDP(3) - t119 * MDP(4)) * t113 + (t112 * t90 - t115 * t88) * MDP(5) * pkin(2) + (-t88 * MDP(11) - t78 * MDP(13) + t79 * MDP(14) - t132) * t121 + (t88 * MDP(12) + (-t111 * t79 - t114 * t78) * MDP(15) + t128 * t84) * t118; MDP(2) + t109 * MDP(6) + (t109 * t104 ^ 2 + t82 ^ 2 + t83 ^ 2) * MDP(16) + t110 * MDP(21) - (t91 * t150 - t143) * t92 + (t112 ^ 2 + t115 ^ 2) * MDP(5) * pkin(2) ^ 2 + 0.2e1 * (-t105 * MDP(11) + t118 * MDP(7) + t133) * t121 + 0.2e1 * (t109 * t145 - t82 * t121) * MDP(13) + 0.2e1 * (t109 * t104 * t114 + t83 * t121) * MDP(14) + 0.2e1 * (-t76 * t121 + t95 * t91) * MDP(22) + (t77 * t121 - t95 * t92) * t149 + 0.2e1 * (t105 * MDP(12) + (-t111 * t83 - t114 * t82) * MDP(15)) * t118; t116 * MDP(5) + (t135 * t118 - t84 * t121) * MDP(16); (t134 - t144) * t118 * MDP(16); MDP(5) + (t141 * t109 + t110) * MDP(16); -t85 * MDP(12) - t125 * t84 + t152 * t135; -t98 * t143 + (-t98 * t91 + t92 * t97) * MDP(18) + (t106 * t91 + t95 * t97) * MDP(22) + (-t106 * t92 + t95 * t98) * MDP(23) + (-t104 * MDP(12) + t130 * qJ(5) + MDP(9) - t127) * t121 + (MDP(8) - t130 * pkin(4) + (-MDP(11) + t129) * t104) * t118 + t152 * t134; (-MDP(12) + t137 + t153) * t118 + t125 * t121; 0.2e1 * t106 * t142 + MDP(10) + (0.2e1 * t138 - 0.2e1 * t139 + t146) * pkin(4) + (MDP(17) * t98 + t106 * t149 + t97 * t150) * t98 + (0.2e1 * t137 + t153) * qJ(5); t84 * MDP(16); t128 * t118 + t131; -t121 * MDP(16); t126; MDP(16); t132; -t121 * MDP(21) + t76 * MDP(22) - t77 * MDP(23) - t133; -t131; t127; 0; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
