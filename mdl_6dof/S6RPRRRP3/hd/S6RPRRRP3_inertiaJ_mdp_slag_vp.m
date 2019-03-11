% Calculate joint inertia matrix for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRRP3_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:05:18
% EndTime: 2019-03-09 06:05:21
% DurationCPUTime: 0.85s
% Computational Cost: add. (905->189), mult. (1640->262), div. (0->0), fcn. (1603->8), ass. (0->77)
t136 = sin(qJ(4));
t139 = cos(qJ(4));
t143 = MDP(17) * t136 + MDP(18) * t139;
t135 = sin(qJ(5));
t138 = cos(qJ(5));
t163 = t138 * t139;
t114 = t135 * t136 - t163;
t115 = t135 * t139 + t136 * t138;
t154 = MDP(25) - MDP(28);
t155 = MDP(24) + MDP(26);
t172 = -pkin(9) - pkin(8);
t118 = t172 * t139;
t153 = t172 * t136;
t98 = -t118 * t135 - t138 * t153;
t99 = -t138 * t118 + t135 * t153;
t147 = t115 * MDP(21) - t114 * MDP(22) - t154 * t99 - t155 * t98;
t184 = t136 * MDP(14) + t139 * MDP(15) - t143 * pkin(8) + t147;
t158 = t135 * MDP(25);
t181 = (MDP(24) * t138 - t158) * pkin(4);
t137 = sin(qJ(3));
t107 = t115 * t137;
t165 = t136 * t137;
t108 = -t135 * t165 + t137 * t163;
t160 = t108 * MDP(21) - t107 * MDP(22);
t134 = cos(pkin(10));
t123 = -pkin(1) * t134 - pkin(2);
t140 = cos(qJ(3));
t113 = -pkin(3) * t140 - pkin(8) * t137 + t123;
t109 = t139 * t113;
t133 = sin(pkin(10));
t121 = pkin(1) * t133 + pkin(7);
t168 = t121 * t136;
t171 = pkin(9) * t137;
t87 = -t139 * t171 + t109 + (-pkin(4) - t168) * t140;
t166 = t121 * t140;
t152 = t139 * t166;
t89 = t152 + (t113 - t171) * t136;
t82 = -t135 * t89 + t138 * t87;
t179 = t82 * MDP(24) + t160;
t177 = -2 * MDP(20);
t176 = 0.2e1 * MDP(25);
t175 = 2 * MDP(26);
t174 = 2 * MDP(27);
t173 = 2 * MDP(28);
t170 = t140 * pkin(5);
t83 = t135 * t87 + t138 * t89;
t169 = t107 * t115;
t167 = t121 * t139;
t164 = t136 * t139;
t162 = t140 * qJ(6);
t110 = pkin(4) * t165 + t121 * t137;
t159 = MDP(19) * t115;
t157 = t137 * MDP(11);
t156 = MDP(16) + MDP(23);
t127 = -pkin(4) * t139 - pkin(3);
t151 = MDP(13) * t164;
t150 = pkin(5) * t175 + MDP(23);
t149 = -t107 * t155 - t108 * t154;
t148 = -pkin(5) * t107 + qJ(6) * t108;
t146 = MDP(14) * t139 - MDP(15) * t136;
t144 = t139 * MDP(17) - t136 * MDP(18);
t132 = t140 ^ 2;
t131 = t139 ^ 2;
t130 = t137 ^ 2;
t129 = t136 ^ 2;
t128 = t135 * pkin(4);
t125 = pkin(4) * t138 + pkin(5);
t122 = t128 + qJ(6);
t106 = t108 ^ 2;
t93 = t113 * t136 + t152;
t92 = -t136 * t166 + t109;
t91 = t108 * t114;
t90 = pkin(5) * t114 - qJ(6) * t115 + t127;
t84 = -t148 + t110;
t81 = -t82 + t170;
t80 = -t162 + t83;
t1 = [t106 * MDP(19) + (t80 ^ 2 + t81 ^ 2 + t84 ^ 2) * MDP(29) + 0.2e1 * t123 * t157 + t108 * t107 * t177 + MDP(1) + (t133 ^ 2 + t134 ^ 2) * MDP(4) * pkin(1) ^ 2 + t156 * t132 + (MDP(12) * t131 + MDP(5) - 0.2e1 * t151) * t130 + 0.2e1 * (-t123 * MDP(10) + (MDP(6) - t146) * t137 - t160) * t140 + 0.2e1 * (t130 * t168 - t140 * t92) * MDP(17) + 0.2e1 * (t130 * t167 + t140 * t93) * MDP(18) + 0.2e1 * (t107 * t110 - t140 * t82) * MDP(24) + (t108 * t110 + t140 * t83) * t176 + (t107 * t84 + t140 * t81) * t175 + (-t107 * t80 + t108 * t81) * t174 + (-t108 * t84 - t140 * t80) * t173; (t107 * t81 + t108 * t80 - t140 * t84) * MDP(29); MDP(4) + (t107 ^ 2 + t106 + t132) * MDP(29); t108 * t159 + (-t91 - t169) * MDP(20) + (t107 * t127 + t110 * t114) * MDP(24) + (t108 * t127 + t110 * t115) * MDP(25) + (t107 * t90 + t114 * t84) * MDP(26) + (-t107 * t99 + t108 * t98 - t114 * t80 + t115 * t81) * MDP(27) + (-t108 * t90 - t115 * t84) * MDP(28) + (t80 * t99 + t81 * t98 + t84 * t90) * MDP(29) + (-t121 * MDP(11) + MDP(8) - t184) * t140 + (MDP(7) - t121 * MDP(10) + MDP(12) * t164 + (-t129 + t131) * MDP(13) + (-pkin(3) * t136 - t167) * MDP(17) + (-pkin(3) * t139 + t168) * MDP(18)) * t137; -t157 + (-t91 + t169) * MDP(27) + (t107 * t98 + t108 * t99) * MDP(29) + (-t90 * MDP(29) - t114 * t155 - t115 * t154 + MDP(10) + t144) * t140; MDP(9) + t129 * MDP(12) + 0.2e1 * t151 + (t90 ^ 2 + t98 ^ 2 + t99 ^ 2) * MDP(29) + (-0.2e1 * t90 * MDP(28) + t114 * t177 + t127 * t176 + t98 * t174 + t159) * t115 + 0.2e1 * t144 * pkin(3) + 0.2e1 * (MDP(24) * t127 + MDP(26) * t90 - MDP(27) * t99) * t114; t92 * MDP(17) - t93 * MDP(18) + t82 * MDP(26) + (-t107 * t122 - t108 * t125) * MDP(27) + (t122 * t80 - t125 * t81) * MDP(29) + t146 * t137 + ((-pkin(5) - t125) * MDP(26) + (-qJ(6) - t122) * MDP(28) - t181 - t156) * t140 - t154 * t83 + t179; (-t107 * t125 + t108 * t122) * MDP(29) - t143 * t137 + t149; (-t114 * t122 - t115 * t125) * MDP(27) + (t122 * t99 - t125 * t98) * MDP(29) + t184; (t122 ^ 2 + t125 ^ 2) * MDP(29) + 0.2e1 * t181 + t125 * t175 + t122 * t173 + t156; -t140 * MDP(23) - t83 * MDP(25) + (t82 - 0.2e1 * t170) * MDP(26) + (-pkin(5) * t108 - qJ(6) * t107) * MDP(27) + (-0.2e1 * t162 + t83) * MDP(28) + (-pkin(5) * t81 + qJ(6) * t80) * MDP(29) + t179; MDP(29) * t148 + t149; (-pkin(5) * t115 - qJ(6) * t114) * MDP(27) + (-pkin(5) * t98 + qJ(6) * t99) * MDP(29) + t147; (0.2e1 * qJ(6) + t128) * MDP(28) + (pkin(5) * t125 + qJ(6) * t122) * MDP(29) + (t138 * t155 - t158) * pkin(4) + t150; qJ(6) * t173 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(29) + t150; MDP(26) * t140 + MDP(27) * t108 + MDP(29) * t81; t107 * MDP(29); MDP(27) * t115 + MDP(29) * t98; -MDP(29) * t125 - MDP(26); -MDP(29) * pkin(5) - MDP(26); MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
