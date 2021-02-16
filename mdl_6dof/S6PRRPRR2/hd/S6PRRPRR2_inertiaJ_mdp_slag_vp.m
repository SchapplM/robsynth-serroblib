% Calculate joint inertia matrix for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:49
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR2_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:46:41
% EndTime: 2021-01-16 03:46:44
% DurationCPUTime: 0.57s
% Computational Cost: add. (708->153), mult. (1462->239), div. (0->0), fcn. (1686->12), ass. (0->81)
t184 = qJ(4) + pkin(8);
t138 = sin(pkin(12));
t132 = pkin(3) * t138 + pkin(9);
t142 = sin(qJ(5));
t145 = cos(qJ(5));
t154 = t142 * MDP(21) + t145 * MDP(22);
t176 = pkin(10) + t132;
t120 = t176 * t142;
t121 = t176 * t145;
t141 = sin(qJ(6));
t144 = cos(qJ(6));
t127 = t141 * t142 - t144 * t145;
t128 = t141 * t145 + t142 * t144;
t157 = t128 * MDP(25) - t127 * MDP(26) + (-t120 * t144 - t121 * t141) * MDP(28) - (-t120 * t141 + t121 * t144) * MDP(29);
t183 = MDP(18) * t142 + MDP(19) * t145 - t154 * t132 + t157;
t179 = cos(qJ(3));
t135 = -t179 * pkin(3) - pkin(2);
t182 = 0.2e1 * t135;
t181 = -2 * MDP(24);
t180 = 0.2e1 * MDP(29);
t178 = sin(qJ(2));
t140 = cos(pkin(12));
t143 = sin(qJ(3));
t124 = t138 * t143 - t140 * t179;
t177 = pkin(5) * t124;
t139 = sin(pkin(6));
t161 = t139 * t178;
t172 = cos(pkin(6));
t114 = t172 * t143 + t179 * t161;
t150 = -t143 * t161 + t172 * t179;
t101 = t140 * t114 + t138 * t150;
t146 = cos(qJ(2));
t168 = t139 * t146;
t94 = -t101 * t142 - t145 * t168;
t95 = t101 * t145 - t142 * t168;
t87 = -t141 * t95 + t144 * t94;
t88 = t141 * t94 + t144 * t95;
t175 = t87 * MDP(28) - t88 * MDP(29);
t174 = MDP(15) * pkin(3);
t125 = t138 * t179 + t140 * t143;
t109 = t124 * pkin(4) - t125 * pkin(9) + t135;
t130 = t184 * t179;
t158 = t184 * t143;
t112 = t140 * t130 - t138 * t158;
t171 = t112 * t145;
t91 = t171 + (-pkin(10) * t125 + t109) * t142;
t173 = t144 * t91;
t170 = t125 * t142;
t169 = t125 * t145;
t167 = t142 * t145;
t116 = t127 * MDP(28);
t166 = -t128 * MDP(29) - t116;
t102 = t128 * t125;
t97 = t102 * MDP(26);
t103 = t127 * t125;
t98 = t103 * MDP(25);
t165 = t125 * MDP(13);
t164 = t128 * MDP(23);
t163 = MDP(20) + MDP(27);
t162 = t124 * MDP(27) - t97 - t98;
t133 = -pkin(3) * t140 - pkin(4);
t160 = MDP(17) * t167;
t159 = MDP(10) * t179;
t92 = t145 * t109 - t112 * t142;
t90 = -pkin(10) * t169 + t177 + t92;
t83 = -t141 * t91 + t144 * t90;
t110 = t130 * t138 + t140 * t158;
t156 = t135 * MDP(15) + t165;
t155 = MDP(21) * t145 - MDP(22) * t142;
t153 = -MDP(12) - t155;
t152 = (MDP(28) * t144 - MDP(29) * t141) * pkin(5);
t151 = (MDP(18) * t145 - MDP(19) * t142) * t125;
t149 = -t153 + t166;
t137 = t145 ^ 2;
t136 = t142 ^ 2;
t129 = -pkin(5) * t145 + t133;
t99 = t114 * t138 - t140 * t150;
t96 = pkin(5) * t170 + t110;
t93 = t109 * t142 + t171;
t84 = t141 * t90 + t173;
t1 = [MDP(1) + (t139 ^ 2 * t146 ^ 2 + t101 ^ 2 + t99 ^ 2) * MDP(15); (-t101 * t124 + t125 * t99) * MDP(14) + (t101 * t112 + t110 * t99) * MDP(15) + (t124 * t94 + t99 * t170) * MDP(21) + (-t124 * t95 + t99 * t169) * MDP(22) + (t102 * t99 + t124 * t87) * MDP(28) + (-t103 * t99 - t124 * t88) * MDP(29) + (-t178 * MDP(4) + (-MDP(11) * t143 - t124 * MDP(12) + MDP(3) - t156 + t159) * t146) * t139; MDP(2) + 0.2e1 * pkin(2) * t159 + t165 * t182 + (t110 ^ 2 + t112 ^ 2 + t135 ^ 2) * MDP(15) + (t137 * MDP(16) - 0.2e1 * t160) * t125 ^ 2 + t163 * t124 ^ 2 - (-t103 * MDP(23) + t102 * t181) * t103 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t143 + 0.2e1 * t179 * MDP(6)) * t143 + (MDP(12) * t182 + 0.2e1 * t151 - 0.2e1 * t97 - 0.2e1 * t98) * t124 + 0.2e1 * (t110 * t125 - t112 * t124) * MDP(14) + 0.2e1 * (t110 * t170 + t124 * t92) * MDP(21) + 0.2e1 * (t110 * t169 - t124 * t93) * MDP(22) + 0.2e1 * (t102 * t96 + t124 * t83) * MDP(28) + (-t103 * t96 - t124 * t84) * t180; t150 * MDP(10) - t114 * MDP(11) + (t138 * t174 - MDP(13)) * t101 + (-t140 * t174 - t149) * t99; t143 * MDP(7) + t179 * MDP(8) - t112 * MDP(13) - t103 * t164 + (-t102 * t128 + t103 * t127) * MDP(24) + (t102 * t129 + t127 * t96) * MDP(28) + (-t103 * t129 + t128 * t96) * MDP(29) + (-t143 * MDP(10) - t179 * MDP(11)) * pkin(8) + t153 * t110 + (MDP(16) * t167 + (-t136 + t137) * MDP(17) + t154 * t133) * t125 + t183 * t124 + ((-t124 * t138 - t125 * t140) * MDP(14) + (-t110 * t140 + t112 * t138) * MDP(15)) * pkin(3); 0.2e1 * t160 + 0.2e1 * t129 * t116 + t136 * MDP(16) + MDP(9) + (t138 ^ 2 + t140 ^ 2) * MDP(15) * pkin(3) ^ 2 + (t127 * t181 + t129 * t180 + t164) * t128 - 0.2e1 * t155 * t133 + 0.2e1 * (MDP(12) * t140 - MDP(13) * t138) * pkin(3); -MDP(15) * t168; t149 * t124 + t156; 0; MDP(15); MDP(21) * t94 - MDP(22) * t95 + t175; t124 * MDP(20) + t92 * MDP(21) - t93 * MDP(22) + (t144 * t177 + t83) * MDP(28) + (-t173 + (-t90 - t177) * t141) * MDP(29) + t151 + t162; t183; t155 + t166; 0.2e1 * t152 + t163; t175; t83 * MDP(28) - t84 * MDP(29) + t162; t157; t166; MDP(27) + t152; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
