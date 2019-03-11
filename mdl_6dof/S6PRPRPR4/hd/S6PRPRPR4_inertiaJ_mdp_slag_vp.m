% Calculate joint inertia matrix for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRPR4_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:41:28
% EndTime: 2019-03-08 19:41:29
% DurationCPUTime: 0.44s
% Computational Cost: add. (694->135), mult. (1459->207), div. (0->0), fcn. (1698->12), ass. (0->76)
t169 = MDP(8) * qJ(3);
t116 = sin(pkin(12));
t119 = cos(pkin(12));
t151 = t116 ^ 2 + t119 ^ 2;
t168 = t151 * MDP(18);
t110 = -t119 * pkin(5) - pkin(4);
t167 = 0.2e1 * t110;
t120 = cos(pkin(11));
t111 = -t120 * pkin(3) - pkin(2);
t166 = 0.2e1 * t111;
t165 = -2 * MDP(21);
t164 = cos(qJ(4));
t163 = pkin(2) * MDP(8);
t162 = pkin(8) + qJ(3);
t161 = pkin(9) + qJ(5);
t117 = sin(pkin(11));
t123 = sin(qJ(4));
t101 = t123 * t117 - t164 * t120;
t103 = t164 * t117 + t123 * t120;
t92 = t101 * pkin(4) - t103 * qJ(5) + t111;
t105 = t162 * t117;
t107 = t162 * t120;
t97 = -t123 * t105 + t164 * t107;
t82 = t116 * t92 + t119 * t97;
t160 = pkin(4) * MDP(19);
t159 = t103 * t116;
t158 = t103 * t119;
t157 = t117 * MDP(6);
t118 = sin(pkin(6));
t124 = sin(qJ(2));
t156 = t118 * t124;
t126 = cos(qJ(2));
t155 = t118 * t126;
t154 = t120 * MDP(5);
t122 = sin(qJ(6));
t125 = cos(qJ(6));
t102 = t125 * t116 + t122 * t119;
t88 = t102 * t103;
t153 = t88 * MDP(23);
t100 = t122 * t116 - t125 * t119;
t89 = t100 * t103;
t152 = t89 * MDP(22);
t149 = MDP(20) * t102;
t148 = t100 * MDP(25);
t147 = t101 * MDP(24);
t146 = t116 * MDP(17);
t145 = t119 * MDP(16);
t81 = -t116 * t97 + t119 * t92;
t144 = t151 * MDP(19);
t143 = t82 * t116 + t81 * t119;
t142 = -t81 * t116 + t82 * t119;
t121 = cos(pkin(6));
t98 = -t117 * t156 + t121 * t120;
t99 = t121 * t117 + t120 * t156;
t87 = t123 * t98 + t164 * t99;
t83 = -t116 * t87 - t119 * t155;
t84 = -t116 * t155 + t119 * t87;
t141 = t84 * t116 + t83 * t119;
t79 = t101 * pkin(5) - pkin(9) * t158 + t81;
t80 = -pkin(9) * t159 + t82;
t138 = (-t122 * t80 + t125 * t79) * MDP(25) - (t122 * t79 + t125 * t80) * MDP(26);
t137 = (-t122 * t84 + t125 * t83) * MDP(25) - (t122 * t83 + t125 * t84) * MDP(26);
t136 = t88 * MDP(25) - t89 * MDP(26);
t135 = MDP(16) * t116 + MDP(17) * t119;
t134 = -t102 * MDP(26) - t148;
t95 = t164 * t105 + t123 * t107;
t133 = -t154 + t157 - t163;
t132 = t95 * MDP(19) + t136;
t104 = t161 * t116;
t106 = t161 * t119;
t131 = t102 * MDP(22) - t100 * MDP(23) + (-t125 * t104 - t122 * t106) * MDP(25) - (-t122 * t104 + t125 * t106) * MDP(26);
t130 = t134 + t145 - t146;
t129 = -t130 - t160;
t86 = t123 * t99 - t164 * t98;
t85 = pkin(5) * t159 + t95;
t1 = [MDP(1) + (t118 ^ 2 * t126 ^ 2 + t98 ^ 2 + t99 ^ 2) * MDP(8) + (t83 ^ 2 + t84 ^ 2 + t86 ^ 2) * MDP(19); (t83 * t81 + t84 * t82) * MDP(19) + t132 * t86 + (t83 * MDP(16) - t84 * MDP(17) + t137) * t101 + (-t141 * MDP(18) + t135 * t86) * t103 + (-t124 * MDP(4) + (-t101 * MDP(14) - t103 * MDP(15) + MDP(3) - t133) * t126) * t118 + (MDP(7) + t169) * (-t98 * t117 + t99 * t120); MDP(2) + (t81 ^ 2 + t82 ^ 2 + t95 ^ 2) * MDP(19) - (-t89 * MDP(20) + t88 * t165) * t89 + (MDP(15) * t166 + MDP(9) * t103) * t103 + (0.2e1 * t154 - 0.2e1 * t157 + t163) * pkin(2) + (-0.2e1 * t103 * MDP(10) + MDP(14) * t166 + t147 - 0.2e1 * t152 - 0.2e1 * t153) * t101 + 0.2e1 * t136 * t85 + 0.2e1 * (t81 * MDP(16) - t82 * MDP(17) + t138) * t101 + 0.2e1 * (-t143 * MDP(18) + t135 * t95) * t103 + (0.2e1 * MDP(7) + t169) * (t117 ^ 2 + t120 ^ 2) * qJ(3); t141 * MDP(19) - MDP(8) * t155; t143 * MDP(19) + (MDP(15) - t168) * t103 + (MDP(14) + t130) * t101 + t133; MDP(8) + t144; -t87 * MDP(15) + (-MDP(14) + t129) * t86 + (MDP(19) * qJ(5) + MDP(18)) * (-t83 * t116 + t84 * t119); t103 * MDP(11) - t95 * MDP(14) - t97 * MDP(15) + (-pkin(4) * t159 - t95 * t119) * MDP(16) + (-pkin(4) * t158 + t95 * t116) * MDP(17) + t142 * MDP(18) + (-t95 * pkin(4) + t142 * qJ(5)) * MDP(19) - t89 * t149 + (t89 * t100 - t102 * t88) * MDP(21) + (t85 * t100 + t110 * t88) * MDP(25) + (t85 * t102 - t110 * t89) * MDP(26) + (-t135 * qJ(5) - MDP(12) + t131) * t101; 0; t148 * t167 + MDP(13) + (0.2e1 * t145 - 0.2e1 * t146 + t160) * pkin(4) + (MDP(26) * t167 + t100 * t165 + t149) * t102 + (t144 * qJ(5) + 0.2e1 * t168) * qJ(5); t86 * MDP(19); t135 * t103 + t132; 0; t129; MDP(19); t137; t138 + t147 - t152 - t153; t134; t131; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
