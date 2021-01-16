% Calculate joint inertia matrix for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR1_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:06:07
% EndTime: 2021-01-16 02:06:10
% DurationCPUTime: 0.56s
% Computational Cost: add. (784->154), mult. (1618->232), div. (0->0), fcn. (1866->12), ass. (0->79)
t179 = qJ(4) + pkin(8);
t169 = MDP(15) * pkin(3);
t126 = sin(pkin(12));
t129 = cos(pkin(12));
t162 = t126 ^ 2 + t129 ^ 2;
t178 = t162 * MDP(18);
t132 = sin(qJ(3));
t128 = sin(pkin(6));
t171 = sin(qJ(2));
t153 = t128 * t171;
t168 = cos(pkin(6));
t172 = cos(qJ(3));
t108 = t168 * t132 + t172 * t153;
t127 = sin(pkin(11));
t130 = cos(pkin(11));
t139 = -t132 * t153 + t168 * t172;
t95 = t108 * t127 - t130 * t139;
t177 = t95 ^ 2;
t117 = t179 * t172;
t151 = t179 * t132;
t104 = t117 * t127 + t130 * t151;
t176 = t104 ^ 2;
t122 = -pkin(3) * t130 - pkin(4);
t116 = -pkin(5) * t129 + t122;
t175 = 0.2e1 * t116;
t123 = -t172 * pkin(3) - pkin(2);
t174 = 0.2e1 * t123;
t173 = -2 * MDP(21);
t119 = pkin(3) * t127 + qJ(5);
t170 = pkin(9) + t119;
t111 = t127 * t132 - t130 * t172;
t113 = t127 * t172 + t130 * t132;
t102 = t111 * pkin(4) - t113 * qJ(5) + t123;
t106 = t130 * t117 - t127 * t151;
t89 = t126 * t102 + t129 * t106;
t131 = sin(qJ(6));
t133 = cos(qJ(6));
t112 = t126 * t131 - t133 * t129;
t94 = t112 * t113;
t167 = MDP(22) * t94;
t114 = t126 * t133 + t129 * t131;
t93 = t114 * t113;
t166 = MDP(23) * t93;
t165 = t113 * t126;
t164 = t113 * t129;
t134 = cos(qJ(2));
t163 = t128 * t134;
t161 = MDP(16) * t129;
t160 = MDP(17) * t126;
t159 = MDP(19) * t122;
t158 = MDP(20) * t114;
t157 = MDP(24) * t111;
t156 = MDP(25) * t112;
t155 = t113 * MDP(13);
t154 = t123 * MDP(15);
t152 = MDP(10) * t172;
t88 = t129 * t102 - t106 * t126;
t150 = t162 * MDP(19);
t149 = t126 * t89 + t129 * t88;
t148 = -t126 * t88 + t129 * t89;
t97 = t130 * t108 + t127 * t139;
t90 = -t126 * t97 - t129 * t163;
t91 = -t126 * t163 + t129 * t97;
t147 = t126 * t91 + t129 * t90;
t86 = pkin(5) * t111 - pkin(9) * t164 + t88;
t87 = -pkin(9) * t165 + t89;
t145 = (-t131 * t87 + t133 * t86) * MDP(25) - (t131 * t86 + t133 * t87) * MDP(26);
t144 = MDP(25) * (-t131 * t91 + t133 * t90) - MDP(26) * (t131 * t90 + t133 * t91);
t143 = t93 * MDP(25) - t94 * MDP(26);
t142 = t126 * MDP(16) + t129 * MDP(17);
t141 = -MDP(26) * t114 - t156;
t140 = MDP(14) + t142;
t138 = t141 - t160 + t161;
t109 = t170 * t126;
t110 = t170 * t129;
t137 = t114 * MDP(22) - t112 * MDP(23) + (-t109 * t133 - t110 * t131) * MDP(25) - (-t109 * t131 + t110 * t133) * MDP(26);
t136 = -t138 + t159;
t92 = pkin(5) * t165 + t104;
t1 = [MDP(1) + (t128 ^ 2 * t134 ^ 2 + t97 ^ 2 + t177) * MDP(15) + (t90 ^ 2 + t91 ^ 2 + t177) * MDP(19); t97 * t106 * MDP(15) + (t88 * t90 + t89 * t91) * MDP(19) + ((MDP(15) + MDP(19)) * t104 + t143) * t95 + (-MDP(14) * t97 + MDP(16) * t90 - MDP(17) * t91 + t144) * t111 + (-t147 * MDP(18) + t140 * t95) * t113 + (-t171 * MDP(4) + (-MDP(11) * t132 - t111 * MDP(12) + MDP(3) + t152 - t154 - t155) * t134) * t128; MDP(2) + 0.2e1 * pkin(2) * t152 + t155 * t174 + (t106 ^ 2 + t123 ^ 2 + t176) * MDP(15) + (t88 ^ 2 + t89 ^ 2 + t176) * MDP(19) - (-MDP(20) * t94 + t93 * t173) * t94 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t132 + 0.2e1 * t172 * MDP(6)) * t132 + (MDP(12) * t174 + t157 - 0.2e1 * t166 - 0.2e1 * t167) * t111 + 0.2e1 * t143 * t92 + 0.2e1 * (-t106 * MDP(14) + t88 * MDP(16) - t89 * MDP(17) + t145) * t111 + 0.2e1 * (-t149 * MDP(18) + t140 * t104) * t113; t139 * MDP(10) - t108 * MDP(11) + (t127 * t169 - MDP(13)) * t97 + (-t130 * t169 - MDP(12) + t136) * t95 + (MDP(19) * t119 + MDP(18)) * (-t126 * t90 + t129 * t91); t132 * MDP(7) + t172 * MDP(8) - t104 * MDP(12) - t106 * MDP(13) + (-t104 * t129 + t122 * t165) * MDP(16) + (t104 * t126 + t122 * t164) * MDP(17) + t148 * MDP(18) + (t104 * t122 + t148 * t119) * MDP(19) - t94 * t158 + (t112 * t94 - t114 * t93) * MDP(21) + (t112 * t92 + t116 * t93) * MDP(25) + (t114 * t92 - t116 * t94) * MDP(26) + (-t142 * t119 + t137) * t111 + (-t132 * MDP(10) - t172 * MDP(11)) * pkin(8) + ((-t111 * t127 - t113 * t130) * MDP(14) + (-t104 * t130 + t106 * t127) * MDP(15)) * pkin(3); t156 * t175 + MDP(9) + (t159 + 0.2e1 * t160 - 0.2e1 * t161) * t122 + (MDP(26) * t175 + t112 * t173 + t158) * t114 + (t150 * t119 + 0.2e1 * t178) * t119 + (0.2e1 * MDP(12) * t130 - 0.2e1 * MDP(13) * t127 + (t127 ^ 2 + t130 ^ 2) * t169) * pkin(3); -MDP(15) * t163 + t147 * MDP(19); t154 + t149 * MDP(19) + (MDP(13) - t178) * t113 + (MDP(12) + t138) * t111; 0; MDP(15) + t150; t95 * MDP(19); MDP(19) * t104 + t142 * t113 + t143; t136; 0; MDP(19); t144; t145 + t157 - t166 - t167; t137; t141; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
