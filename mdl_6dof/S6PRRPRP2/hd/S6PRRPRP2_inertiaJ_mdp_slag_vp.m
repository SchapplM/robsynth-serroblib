% Calculate joint inertia matrix for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPRP2_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:32:17
% EndTime: 2019-03-08 21:32:18
% DurationCPUTime: 0.49s
% Computational Cost: add. (747->152), mult. (1494->230), div. (0->0), fcn. (1671->10), ass. (0->66)
t111 = sin(pkin(11));
t113 = cos(pkin(11));
t115 = sin(qJ(3));
t153 = cos(qJ(3));
t101 = t111 * t153 + t113 * t115;
t157 = 0.2e1 * t101;
t156 = qJ(4) + pkin(8);
t106 = pkin(3) * t111 + pkin(9);
t132 = t106 * MDP(24) + MDP(22);
t114 = sin(qJ(5));
t117 = cos(qJ(5));
t137 = MDP(20) - MDP(23);
t138 = MDP(19) + MDP(21);
t155 = -t137 * t114 + t138 * t117;
t112 = sin(pkin(6));
t116 = sin(qJ(2));
t147 = t112 * t116;
t150 = cos(pkin(6));
t120 = -t115 * t147 + t150 * t153;
t97 = t150 * t115 + t153 * t147;
t86 = t111 * t97 - t113 * t120;
t154 = t86 ^ 2;
t100 = t111 * t115 - t113 * t153;
t152 = pkin(5) * t100;
t108 = -t153 * pkin(3) - pkin(2);
t91 = t100 * pkin(4) - t101 * pkin(9) + t108;
t103 = t156 * t153;
t135 = t156 * t115;
t95 = t113 * t103 - t111 * t135;
t81 = t114 * t91 + t117 * t95;
t151 = MDP(13) * pkin(3);
t149 = qJ(6) * t100;
t148 = t100 * t106;
t118 = cos(qJ(2));
t146 = t112 * t118;
t107 = -pkin(3) * t113 - pkin(4);
t130 = -pkin(5) * t117 - qJ(6) * t114;
t98 = t107 + t130;
t145 = t98 * MDP(24);
t109 = t114 ^ 2;
t110 = t117 ^ 2;
t144 = t109 + t110;
t143 = MDP(15) * t117;
t142 = MDP(22) * t101;
t141 = t100 * MDP(18);
t139 = t108 * MDP(13);
t136 = MDP(10) * t153;
t134 = t114 * t95 - t117 * t91;
t133 = -MDP(24) * pkin(5) - MDP(21);
t93 = t103 * t111 + t113 * t135;
t131 = MDP(19) - t133;
t129 = pkin(5) * t114 - qJ(6) * t117;
t78 = t149 + t81;
t79 = t134 - t152;
t128 = t114 * t78 - t117 * t79;
t88 = t111 * t120 + t113 * t97;
t84 = t114 * t88 + t117 * t146;
t85 = -t114 * t146 + t117 * t88;
t127 = t114 * t85 - t117 * t84;
t125 = -t101 * t98 + t148;
t124 = MDP(24) * qJ(6) - t137;
t123 = -t134 * MDP(19) - t81 * MDP(20);
t122 = t101 * t107 - t148;
t121 = t117 * MDP(16) - t114 * MDP(17);
t82 = t129 * t101 + t93;
t1 = [MDP(1) + (t112 ^ 2 * t118 ^ 2 + t88 ^ 2 + t154) * MDP(13) + (t84 ^ 2 + t85 ^ 2 + t154) * MDP(24); (t86 * t93 + t88 * t95) * MDP(13) + (t78 * t85 + t79 * t84 + t82 * t86) * MDP(24) + (-t88 * MDP(12) - t137 * t85 - t138 * t84) * t100 + (-t116 * MDP(4) + (-MDP(11) * t115 + MDP(3) + t136 - t139) * t118) * t112 + (-t127 * MDP(22) + (t138 * t114 + t137 * t117 + MDP(12)) * t86) * t101; MDP(2) + 0.2e1 * pkin(2) * t136 + (t108 ^ 2 + t93 ^ 2 + t95 ^ 2) * MDP(13) + (t78 ^ 2 + t79 ^ 2 + t82 ^ 2) * MDP(24) + (t110 * MDP(14) - 0.2e1 * t114 * t143) * t101 ^ 2 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t115 + 0.2e1 * t153 * MDP(6)) * t115 + (t121 * t157 + t141) * t100 + 0.2e1 * (-t95 * MDP(12) - t79 * MDP(21) + t78 * MDP(23) + t123) * t100 + (-t128 * MDP(22) + (t114 * MDP(21) - t117 * MDP(23)) * t82 + (t114 * MDP(19) + t117 * MDP(20) + MDP(12)) * t93) * t157; t120 * MDP(10) - t97 * MDP(11) + t88 * t111 * t151 + (-t113 * t151 + t145 - t155) * t86 + t132 * (t114 * t84 + t117 * t85); t82 * t145 + t115 * MDP(7) + t153 * MDP(8) + (-t109 + t110) * MDP(15) * t101 + (-t115 * MDP(10) - t153 * MDP(11)) * pkin(8) + (t100 * MDP(17) - t93 * MDP(19) + t122 * MDP(20) - t82 * MDP(21) + t125 * MDP(23) + t132 * t78) * t117 + (t117 * t101 * MDP(14) + t100 * MDP(16) + t122 * MDP(19) + t93 * MDP(20) - t125 * MDP(21) - t82 * MDP(23) + t132 * t79) * t114 + ((-t100 * t111 - t101 * t113) * MDP(12) + (t111 * t95 - t113 * t93) * MDP(13)) * pkin(3); MDP(9) + t109 * MDP(14) + (t144 * t106 ^ 2 + t98 ^ 2) * MDP(24) + (t111 ^ 2 + t113 ^ 2) * MDP(13) * pkin(3) ^ 2 + 0.2e1 * t144 * MDP(22) * t106 + 0.2e1 * (-t107 * MDP(19) - t98 * MDP(21)) * t117 + 0.2e1 * (t107 * MDP(20) - t98 * MDP(23) + t143) * t114; -MDP(13) * t146 + t127 * MDP(24); t128 * MDP(24) + t155 * t100 - t144 * t142 + t139; 0; t144 * MDP(24) + MDP(13); t124 * t85 - t131 * t84; t141 + (-t134 + 0.2e1 * t152) * MDP(21) + (0.2e1 * t149 + t81) * MDP(23) + (-pkin(5) * t79 + qJ(6) * t78) * MDP(24) + (t130 * MDP(22) + t121) * t101 + t123; t114 * MDP(16) + t117 * MDP(17) - t129 * MDP(22) + (-t131 * t114 + t124 * t117) * t106; t124 * t114 + t131 * t117; MDP(18) + 0.2e1 * pkin(5) * MDP(21) + 0.2e1 * qJ(6) * MDP(23) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(24); t84 * MDP(24); -t100 * MDP(21) + t79 * MDP(24) + t117 * t142; t132 * t114; -t117 * MDP(24); t133; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
