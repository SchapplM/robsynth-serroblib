% Calculate joint inertia matrix for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR4_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:14:29
% EndTime: 2019-03-08 22:14:31
% DurationCPUTime: 0.43s
% Computational Cost: add. (377->128), mult. (748->187), div. (0->0), fcn. (770->10), ass. (0->70)
t111 = sin(qJ(6));
t115 = cos(qJ(6));
t163 = MDP(28) * t111 + MDP(29) * t115;
t162 = pkin(8) * MDP(15) + MDP(13);
t140 = t115 * MDP(28);
t126 = -t111 * MDP(29) + t140;
t124 = MDP(21) + t126;
t112 = sin(qJ(5));
t116 = cos(qJ(5));
t113 = sin(qJ(3));
t117 = cos(qJ(3));
t110 = sin(pkin(6));
t114 = sin(qJ(2));
t151 = t110 * t114;
t153 = cos(pkin(6));
t87 = t113 * t151 - t153 * t117;
t88 = t153 * t113 + t117 * t151;
t82 = t88 * t112 - t87 * t116;
t83 = t87 * t112 + t88 * t116;
t161 = -t83 * MDP(22) - t124 * t82;
t160 = -t112 * MDP(22) + t124 * t116;
t97 = -t117 * pkin(3) - t113 * qJ(4) - pkin(2);
t89 = t117 * pkin(4) - t97;
t159 = 0.2e1 * t89;
t158 = pkin(8) - pkin(9);
t119 = -pkin(3) - pkin(4);
t95 = t112 * qJ(4) - t116 * t119;
t93 = pkin(5) + t95;
t157 = pkin(5) + t93;
t92 = -t112 * t117 + t116 * t113;
t155 = t111 * t92;
t154 = t115 * t92;
t152 = MDP(22) * t92;
t118 = cos(qJ(2));
t150 = t110 * t118;
t91 = t112 * t113 + t116 * t117;
t148 = t91 * MDP(27);
t147 = t95 * MDP(21);
t96 = t116 * qJ(4) + t112 * t119;
t146 = t96 * MDP(22);
t107 = t113 ^ 2;
t145 = t117 ^ 2 + t107;
t141 = t115 * MDP(24);
t106 = t111 ^ 2;
t139 = t106 * MDP(23) + MDP(20);
t138 = -MDP(11) + MDP(14);
t137 = t111 * t141;
t136 = 0.2e1 * t137 + t139;
t135 = -pkin(3) * MDP(15) - MDP(12);
t134 = -pkin(5) * t92 - pkin(10) * t91;
t94 = -pkin(10) + t96;
t133 = -t91 * t94 + t92 * t93;
t132 = -MDP(10) + t135;
t129 = qJ(4) * MDP(15) + t138;
t98 = t158 * t113;
t99 = t158 * t117;
t85 = t112 * t99 - t116 * t98;
t128 = -t91 * MDP(26) + t85 * MDP(28);
t127 = MDP(25) * t115 - MDP(26) * t111;
t123 = 0.2e1 * t126;
t122 = -MDP(23) * t154 - t91 * MDP(25) - t85 * MDP(29);
t108 = t115 ^ 2;
t86 = t112 * t98 + t116 * t99;
t121 = -t91 * MDP(19) - t85 * MDP(21) - t86 * MDP(22) + (MDP(18) - (t106 - t108) * MDP(24)) * t92;
t81 = t91 * pkin(5) - t92 * pkin(10) + t89;
t80 = t111 * t150 + t83 * t115;
t79 = -t83 * t111 + t115 * t150;
t78 = t111 * t81 + t115 * t86;
t77 = -t111 * t86 + t115 * t81;
t1 = [MDP(1) + (t110 ^ 2 * t118 ^ 2 + t87 ^ 2 + t88 ^ 2) * MDP(15); (t82 * t155 + t79 * t91) * MDP(28) + (t82 * t154 - t80 * t91) * MDP(29) + (-t114 * MDP(4) + (-MDP(15) * t97 + MDP(21) * t91 + t152 + MDP(3) + (MDP(10) + MDP(12)) * t117 + t138 * t113) * t118) * t110 + t162 * (t87 * t113 + t88 * t117); MDP(2) + t107 * MDP(5) + (t145 * pkin(8) ^ 2 + t97 ^ 2) * MDP(15) + t152 * t159 + (MDP(21) * t159 + t148) * t91 + 0.2e1 * (t85 * t155 + t77 * t91) * MDP(28) + 0.2e1 * (t85 * t154 - t78 * t91) * MDP(29) + 0.2e1 * t145 * MDP(13) * pkin(8) + 0.2e1 * (pkin(2) * MDP(10) - t97 * MDP(12)) * t117 + 0.2e1 * (-pkin(2) * MDP(11) - t97 * MDP(14) + t117 * MDP(6)) * t113 + (0.2e1 * (-MDP(17) + t127) * t91 + (t108 * MDP(23) + MDP(16) - 0.2e1 * t137) * t92) * t92; t129 * t88 + t132 * t87 - t161; t113 * MDP(7) + t117 * MDP(8) + (-t113 * pkin(3) + t117 * qJ(4)) * MDP(13) + (t133 * MDP(29) + t128) * t115 + (t133 * MDP(28) + t122) * t111 + (t132 * t113 + t129 * t117) * pkin(8) - t121; MDP(9) + 0.2e1 * pkin(3) * MDP(12) + 0.2e1 * qJ(4) * MDP(14) + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(15) + t93 * t123 + 0.2e1 * t147 + 0.2e1 * t146 + t136; t87 * MDP(15); t162 * t113 + t163 * (-t112 * t91 - t116 * t92); t135 - t160; MDP(15); t161; (t134 * MDP(29) - t128) * t115 + (t134 * MDP(28) - t122) * t111 + t121; -t147 - t146 - t157 * t140 + (t157 * MDP(29) - 0.2e1 * t141) * t111 - t139; t160; pkin(5) * t123 + t136; t79 * MDP(28) - t80 * MDP(29); t77 * MDP(28) - t78 * MDP(29) + t127 * t92 + t148; (-t94 * MDP(29) - MDP(26)) * t115 + (-t94 * MDP(28) - MDP(25)) * t111; -t163 * t112; t111 * MDP(25) + t115 * MDP(26) - pkin(10) * t163; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
