% Calculate joint inertia matrix for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP1_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:38:51
% EndTime: 2021-01-16 02:38:54
% DurationCPUTime: 0.67s
% Computational Cost: add. (768->171), mult. (1546->247), div. (0->0), fcn. (1726->10), ass. (0->72)
t114 = cos(qJ(5));
t112 = sin(qJ(5));
t130 = MDP(22) + MDP(24);
t124 = t130 * t112;
t131 = MDP(21) + MDP(23);
t158 = t131 * t114 + MDP(12) - t124;
t109 = sin(pkin(11));
t111 = cos(pkin(11));
t113 = sin(qJ(3));
t151 = cos(qJ(3));
t98 = t109 * t151 + t111 * t113;
t157 = 0.2e1 * t98;
t156 = qJ(4) + pkin(8);
t110 = sin(pkin(6));
t150 = sin(qJ(2));
t129 = t110 * t150;
t143 = cos(pkin(6));
t117 = -t113 * t129 + t143 * t151;
t93 = t143 * t113 + t151 * t129;
t84 = t109 * t93 - t111 * t117;
t154 = t84 ^ 2;
t106 = -t151 * pkin(3) - pkin(2);
t153 = 0.2e1 * t106;
t97 = t109 * t113 - t111 * t151;
t152 = pkin(5) * t97;
t149 = MDP(15) * pkin(3);
t148 = MDP(25) * pkin(5);
t147 = MDP(26) * pkin(5);
t146 = t112 * t98;
t101 = t156 * t151;
t126 = t156 * t113;
t91 = t111 * t101 - t109 * t126;
t145 = t114 * t91;
t144 = t114 * t98;
t115 = cos(qJ(2));
t142 = t110 * t115;
t104 = pkin(3) * t109 + pkin(9);
t138 = qJ(6) + t104;
t95 = t138 * t114;
t141 = t95 * MDP(24);
t140 = t97 * MDP(20);
t139 = t98 * MDP(13);
t107 = t112 ^ 2;
t108 = t114 ^ 2;
t137 = t107 + t108;
t105 = -pkin(3) * t111 - pkin(4);
t100 = -pkin(5) * t114 + t105;
t136 = t100 * MDP(26);
t135 = t106 * MDP(15);
t134 = t112 * MDP(19);
t133 = t112 * MDP(24);
t132 = t114 * MDP(23);
t128 = t112 * t114 * MDP(17);
t127 = MDP(10) * t151;
t88 = t97 * pkin(4) - t98 * pkin(9) + t106;
t79 = -t112 * t91 + t114 * t88;
t89 = t101 * t109 + t111 * t126;
t125 = MDP(23) + t147;
t123 = -t104 * MDP(21) + MDP(18);
t122 = MDP(21) + t125;
t77 = -qJ(6) * t144 + t152 + t79;
t78 = t145 + (-qJ(6) * t98 + t88) * t112;
t121 = t112 * t78 + t114 * t77;
t86 = t109 * t117 + t111 * t93;
t81 = -t112 * t86 - t114 * t142;
t82 = -t112 * t142 + t114 * t86;
t120 = t112 * t82 + t114 * t81;
t119 = t112 * MDP(23) + t114 * MDP(24);
t118 = t79 * MDP(21) - (t112 * t88 + t145) * MDP(22) - t78 * MDP(24);
t94 = t138 * t112;
t83 = pkin(5) * t146 + t89;
t1 = [MDP(1) + (t110 ^ 2 * t115 ^ 2 + t86 ^ 2 + t154) * MDP(15) + (t81 ^ 2 + t82 ^ 2 + t154) * MDP(26); -t86 * t97 * MDP(14) + (t84 * t89 + t86 * t91) * MDP(15) + (t77 * t81 + t78 * t82 + t83 * t84) * MDP(26) + (t84 * MDP(14) - t120 * MDP(25)) * t98 + t130 * (t84 * t144 - t82 * t97) + t131 * (t84 * t146 + t81 * t97) + (-t150 * MDP(4) + (-MDP(11) * t113 - t97 * MDP(12) + MDP(3) + t127 - t135 - t139) * t115) * t110; MDP(2) + 0.2e1 * pkin(2) * t127 + t139 * t153 + (t106 ^ 2 + t89 ^ 2 + t91 ^ 2) * MDP(15) + (t77 ^ 2 + t78 ^ 2 + t83 ^ 2) * MDP(26) + (t108 * MDP(16) - 0.2e1 * t128) * t98 ^ 2 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t113 + 0.2e1 * t151 * MDP(6)) * t113 + (MDP(12) * t153 + t140 + (MDP(18) * t114 - t134) * t157) * t97 + 0.2e1 * (-t91 * MDP(14) + t77 * MDP(23) + t118) * t97 + (-t121 * MDP(25) + t119 * t83 + (t112 * MDP(21) + t114 * MDP(22) + MDP(14)) * t89) * t157; t117 * MDP(10) - t93 * MDP(11) + (-t112 * t81 + t114 * t82) * MDP(25) + (-t81 * t94 + t82 * t95) * MDP(26) + (t109 * t149 - MDP(13)) * t86 + (-t111 * t149 + t136 - t158) * t84; t113 * MDP(7) + t151 * MDP(8) - t89 * MDP(12) - t91 * MDP(13) + (t100 * t83 - t77 * t94 + t78 * t95) * MDP(26) + (-t107 + t108) * MDP(17) * t98 + (-t94 * MDP(23) - t141) * t97 + (-t113 * MDP(10) - t151 * MDP(11)) * pkin(8) + (t97 * MDP(19) - t89 * MDP(21) + (-t104 * t97 + t105 * t98) * MDP(22) - t83 * MDP(23) + t100 * t98 * MDP(24) + (t94 * t98 + t78) * MDP(25)) * t114 + (t89 * MDP(22) + t83 * MDP(24) - t77 * MDP(25) + t123 * t97 + (t114 * MDP(16) + t105 * MDP(21) + t100 * MDP(23) - t95 * MDP(25)) * t98) * t112 + ((-t109 * t97 - t111 * t98) * MDP(14) + (t109 * t91 - t111 * t89) * MDP(15)) * pkin(3); MDP(9) + t107 * MDP(16) + 0.2e1 * t128 + 0.2e1 * (t112 * t94 + t114 * t95) * MDP(25) + (t94 ^ 2 + t95 ^ 2) * MDP(26) + (t109 ^ 2 + t111 ^ 2) * MDP(15) * pkin(3) ^ 2 + (-0.2e1 * t132 + 0.2e1 * t133 + t136) * t100 + 0.2e1 * (-t114 * MDP(21) + t112 * MDP(22)) * t105 + 0.2e1 * (MDP(12) * t111 - MDP(13) * t109) * pkin(3); -MDP(15) * t142 + t120 * MDP(26); t135 + t121 * MDP(26) + (-t137 * MDP(25) + MDP(13)) * t98 + t158 * t97; (t112 * t95 - t114 * t94) * MDP(26); t137 * MDP(26) + MDP(15); t122 * t81 - t130 * t82; t140 + (t79 + 0.2e1 * t152) * MDP(23) + t77 * t147 + (-t134 + (-MDP(23) * qJ(6) + MDP(18) - t148) * t114) * t98 + t118; -t141 - t125 * t94 + (-t104 * MDP(22) + MDP(19)) * t114 + (t123 - t148) * t112; t122 * t114 - t124; MDP(20) + (0.2e1 * MDP(23) + t147) * pkin(5); t84 * MDP(26); MDP(26) * t83 + t119 * t98; -t132 + t133 + t136; 0; 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
