% Calculate joint inertia matrix for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRPR1_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:40:03
% EndTime: 2019-03-09 01:40:04
% DurationCPUTime: 0.38s
% Computational Cost: add. (670->110), mult. (1222->163), div. (0->0), fcn. (1336->10), ass. (0->66)
t115 = sin(qJ(6));
t117 = cos(qJ(6));
t112 = cos(pkin(11));
t110 = sin(pkin(10));
t113 = cos(pkin(10));
t116 = sin(qJ(4));
t148 = cos(qJ(4));
t95 = t148 * t110 + t116 * t113;
t143 = t112 * t95;
t109 = sin(pkin(11));
t93 = t116 * t110 - t148 * t113;
t114 = cos(pkin(9));
t103 = -t114 * pkin(1) - pkin(2);
t96 = -t113 * pkin(3) + t103;
t80 = t93 * pkin(4) - t95 * qJ(5) + t96;
t111 = sin(pkin(9));
t100 = t111 * pkin(1) + qJ(3);
t147 = pkin(7) + t100;
t88 = t147 * t110;
t89 = t147 * t113;
t84 = -t116 * t88 + t148 * t89;
t74 = -t109 * t84 + t112 * t80;
t72 = t93 * pkin(5) - pkin(8) * t143 + t74;
t144 = t109 * t95;
t75 = t109 * t80 + t112 * t84;
t73 = -pkin(8) * t144 + t75;
t94 = t117 * t109 + t115 * t112;
t81 = t94 * t95;
t92 = t115 * t109 - t117 * t112;
t82 = t92 * t95;
t154 = (-t115 * t73 + t117 * t72) * MDP(25) - (t115 * t72 + t117 * t73) * MDP(26) - t82 * MDP(22) - t81 * MDP(23);
t135 = t109 ^ 2 + t112 ^ 2;
t130 = t135 * MDP(19);
t152 = qJ(5) * t130;
t151 = t135 * MDP(18);
t104 = -t112 * pkin(5) - pkin(4);
t150 = 0.2e1 * t104;
t149 = -2 * MDP(21);
t146 = pkin(8) + qJ(5);
t145 = pkin(4) * MDP(19);
t142 = t103 * MDP(8);
t141 = t110 * MDP(6);
t140 = t113 * MDP(5);
t139 = t82 * MDP(20);
t138 = t95 * t151;
t137 = t92 * MDP(25);
t136 = t95 * MDP(15);
t134 = t110 ^ 2 + t113 ^ 2;
t133 = t109 * MDP(17);
t132 = t112 * MDP(16);
t131 = t134 * MDP(8);
t129 = t75 * t109 + t74 * t112;
t128 = -t74 * t109 + t75 * t112;
t125 = t81 * MDP(25) - t82 * MDP(26);
t124 = -t94 * MDP(26) - t137;
t123 = MDP(16) * t109 + MDP(17) * t112;
t83 = t116 * t89 + t148 * t88;
t97 = t146 * t109;
t98 = t146 * t112;
t122 = t94 * MDP(22) - t92 * MDP(23) + (-t115 * t98 - t117 * t97) * MDP(25) - (-t115 * t97 + t117 * t98) * MDP(26);
t121 = t124 + t132 - t133;
t120 = -t121 - t145;
t91 = t95 ^ 2;
t90 = t93 ^ 2;
t76 = pkin(5) * t144 + t83;
t1 = [MDP(1) + t91 * MDP(9) + 0.2e1 * t96 * t136 + (t74 ^ 2 + t75 ^ 2 + t83 ^ 2) * MDP(19) + t90 * MDP(24) - (t81 * t149 - t139) * t82 + (t111 ^ 2 + t114 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-0.2e1 * t140 + 0.2e1 * t141 + t142) * t103 + 0.2e1 * t125 * t76 + 0.2e1 * (-t129 * MDP(18) + t123 * t83) * t95 + (0.2e1 * t134 * MDP(7) + t131 * t100) * t100 + 0.2e1 * (-t95 * MDP(10) + t96 * MDP(14) + t74 * MDP(16) - t75 * MDP(17) + t154) * t93; (t128 * t95 + t83 * t93) * MDP(19); MDP(4) + t131 + (t135 * t91 + t90) * MDP(19); -t140 + t141 + t142 + t136 - t138 + t129 * MDP(19) + (MDP(14) + t121) * t93; 0; MDP(8) + t130; t95 * MDP(11) - t83 * MDP(14) - t84 * MDP(15) + (-pkin(4) * t144 - t83 * t112) * MDP(16) + (-pkin(4) * t143 + t83 * t109) * MDP(17) + t128 * MDP(18) + (-t83 * pkin(4) + t128 * qJ(5)) * MDP(19) - t94 * t139 + (-t94 * t81 + t82 * t92) * MDP(21) + (t104 * t81 + t76 * t92) * MDP(25) + (-t104 * t82 + t76 * t94) * MDP(26) + (-t123 * qJ(5) - MDP(12) + t122) * t93; t138 + (-MDP(15) + t152) * t95 + (-MDP(14) + t120) * t93; 0; t137 * t150 + MDP(13) + (0.2e1 * t132 - 0.2e1 * t133 + t145) * pkin(4) + (MDP(20) * t94 + MDP(26) * t150 + t92 * t149) * t94 + (0.2e1 * t151 + t152) * qJ(5); t83 * MDP(19) + t123 * t95 + t125; t93 * MDP(19); 0; t120; MDP(19); t93 * MDP(24) + t154; -t125; t124; t122; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
