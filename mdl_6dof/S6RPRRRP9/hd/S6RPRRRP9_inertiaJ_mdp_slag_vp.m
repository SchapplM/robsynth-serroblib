% Calculate joint inertia matrix for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRRP9_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:28:33
% EndTime: 2019-03-09 06:28:35
% DurationCPUTime: 0.59s
% Computational Cost: add. (643->165), mult. (1189->236), div. (0->0), fcn. (1174->6), ass. (0->67)
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t147 = pkin(8) + pkin(9);
t101 = t147 * t114;
t102 = t147 * t117;
t113 = sin(qJ(5));
t116 = cos(qJ(5));
t80 = -t116 * t101 - t102 * t113;
t81 = -t101 * t113 + t102 * t116;
t97 = t113 * t114 - t116 * t117;
t98 = t113 * t117 + t114 * t116;
t125 = t98 * MDP(23) - t97 * MDP(24) + t80 * MDP(26) - t81 * MDP(27);
t150 = MDP(19) * t114 + MDP(20) * t117;
t152 = t114 * MDP(16) + t117 * MDP(17) - pkin(8) * t150 + t125;
t118 = cos(qJ(3));
t90 = t98 * t118;
t137 = t117 * t118;
t139 = t114 * t118;
t92 = -t113 * t139 + t116 * t137;
t151 = t92 * MDP(23) - t90 * MDP(24);
t149 = -2 * MDP(22);
t148 = 2 * MDP(28);
t146 = (pkin(1) * MDP(6));
t145 = pkin(4) * t113;
t115 = sin(qJ(3));
t144 = t115 * pkin(4);
t89 = t98 * t115;
t91 = t97 * t115;
t143 = -t89 * MDP(26) + t91 * MDP(27);
t142 = MDP(29) * pkin(5);
t100 = pkin(3) * t115 - pkin(8) * t118 + qJ(2);
t119 = -pkin(1) - pkin(7);
t136 = t117 * t119;
t127 = t115 * t136;
t77 = t127 + (-pkin(9) * t118 + t100) * t114;
t141 = t116 * t77;
t140 = t114 * t117;
t138 = t114 * t119;
t107 = -pkin(4) * t117 - pkin(3);
t84 = pkin(5) * t97 + t107;
t135 = t84 * MDP(29);
t134 = t92 * MDP(21);
t130 = t116 * MDP(26);
t129 = MDP(18) + MDP(25);
t128 = t115 * MDP(25) + t151;
t126 = MDP(15) * t140;
t95 = t117 * t100;
t75 = -pkin(9) * t137 + t95 + (pkin(4) - t138) * t115;
t70 = -t113 * t77 + t116 * t75;
t96 = pkin(4) * t139 - t118 * t119;
t71 = t113 * t75 + t141;
t124 = t97 * MDP(26) + t98 * MDP(27);
t123 = t117 * MDP(16) - t114 * MDP(17);
t122 = MDP(19) * t117 - MDP(20) * t114;
t112 = t118 ^ 2;
t111 = t117 ^ 2;
t110 = t115 ^ 2;
t109 = t114 ^ 2;
t106 = pkin(4) * t116 + pkin(5);
t83 = t114 * t100 + t127;
t82 = -t115 * t138 + t95;
t76 = t90 * pkin(5) + t96;
t73 = -qJ(6) * t97 + t81;
t72 = -qJ(6) * t98 + t80;
t69 = -qJ(6) * t90 + t71;
t68 = pkin(5) * t115 - qJ(6) * t92 + t70;
t1 = [(t68 ^ 2 + t69 ^ 2 + t76 ^ 2) * MDP(29) + MDP(1) + (t149 * t90 + t134) * t92 + t129 * t110 + ((-2 * MDP(4) + t146) * pkin(1)) + (0.2e1 * t118 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + (t111 * MDP(14) + MDP(7) - 0.2e1 * t126) * t112 + 0.2e1 * (qJ(2) * MDP(12) + (-MDP(8) + t123) * t118 + t151) * t115 + 0.2e1 * (-t112 * t138 + t82 * t115) * MDP(19) + 0.2e1 * (-t112 * t136 - t83 * t115) * MDP(20) + 0.2e1 * (t115 * t70 + t90 * t96) * MDP(26) + 0.2e1 * (-t115 * t71 + t92 * t96) * MDP(27) + (-t68 * t92 - t69 * t90) * t148; MDP(4) - t146 + (-t115 * t89 - t118 * t90) * MDP(26) + (t115 * t91 - t118 * t92) * MDP(27) + (t89 * t92 + t90 * t91) * MDP(28) + (-t118 * t76 - t68 * t89 - t69 * t91) * MDP(29) + t150 * (-t110 - t112); MDP(6) + (t89 ^ 2 + t91 ^ 2 + t112) * MDP(29); t98 * t134 + (-t90 * t98 - t92 * t97) * MDP(22) + (t107 * t90 + t96 * t97) * MDP(26) + (t107 * t92 + t96 * t98) * MDP(27) + (-t68 * t98 - t69 * t97 - t72 * t92 - t73 * t90) * MDP(28) + (t68 * t72 + t69 * t73 + t76 * t84) * MDP(29) + (-t119 * MDP(13) - MDP(10) + t152) * t115 + (MDP(9) + t119 * MDP(12) + MDP(14) * t140 + (-t109 + t111) * MDP(15) + (-pkin(3) * t114 + t136) * MDP(19) + (-pkin(3) * t117 - t138) * MDP(20)) * t118; -t115 * MDP(13) + (t89 * t98 + t91 * t97) * MDP(28) + (-t72 * t89 - t73 * t91) * MDP(29) + (MDP(12) + t122 - t124 - t135) * t118; MDP(11) + t109 * MDP(14) + 0.2e1 * t126 - t73 * t97 * t148 + (t72 ^ 2 + t73 ^ 2 + t84 ^ 2) * MDP(29) + (MDP(21) * t98 - t148 * t72 + t149 * t97) * t98 + 0.2e1 * pkin(3) * t122 + 0.2e1 * t107 * t124; t115 * MDP(18) + t82 * MDP(19) - t83 * MDP(20) + (t116 * t144 + t70) * MDP(26) + (-t141 + (-t75 - t144) * t113) * MDP(27) + (-t106 * t92 - t145 * t90) * MDP(28) + (t106 * t68 + t145 * t69) * MDP(29) + t123 * t118 + t128; (-t106 * t89 - t145 * t91) * MDP(29) - t150 * t115 + t143; (-t106 * t98 - t145 * t97) * MDP(28) + (t106 * t72 + t145 * t73) * MDP(29) + t152; t106 ^ 2 * MDP(29) + (0.2e1 * t130 + (MDP(29) * t145 - 0.2e1 * MDP(27)) * t113) * pkin(4) + t129; t70 * MDP(26) - t71 * MDP(27) + (-t92 * MDP(28) + t68 * MDP(29)) * pkin(5) + t128; -t142 * t89 + t143; (-t98 * MDP(28) + MDP(29) * t72) * pkin(5) + t125; t106 * t142 + MDP(25) + (-MDP(27) * t113 + t130) * pkin(4); MDP(29) * pkin(5) ^ 2 + MDP(25); t76 * MDP(29); -t118 * MDP(29); t135; 0; 0; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
