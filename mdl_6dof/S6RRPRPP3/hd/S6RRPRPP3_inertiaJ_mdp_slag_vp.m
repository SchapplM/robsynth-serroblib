% Calculate joint inertia matrix for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPRPP3_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:57:11
% EndTime: 2019-03-09 09:57:14
% DurationCPUTime: 0.91s
% Computational Cost: add. (1052->209), mult. (1977->284), div. (0->0), fcn. (1978->6), ass. (0->75)
t128 = sin(pkin(9));
t129 = cos(pkin(9));
t132 = sin(qJ(4));
t171 = cos(qJ(4));
t179 = -t132 * t128 + t171 * t129;
t178 = (MDP(14) * qJ(3));
t177 = 2 * MDP(13);
t176 = -2 * MDP(16);
t175 = 2 * MDP(21);
t174 = 2 * MDP(22);
t173 = 2 * MDP(26);
t172 = 2 * MDP(28);
t170 = pkin(7) * t128;
t134 = cos(qJ(2));
t169 = pkin(7) * t134;
t112 = t128 * t171 + t132 * t129;
t133 = sin(qJ(2));
t106 = t112 * t133;
t168 = t106 * pkin(5);
t130 = pkin(4) + qJ(6);
t167 = pkin(8) + qJ(3);
t114 = -t134 * pkin(2) - t133 * qJ(3) - pkin(1);
t105 = t128 * t114 + t129 * t169;
t162 = t128 * t133;
t100 = -pkin(8) * t162 + t105;
t109 = t129 * t114;
t97 = -t129 * t133 * pkin(8) + t109 + (-pkin(3) - t170) * t134;
t166 = t132 * t100 - t171 * t97;
t89 = t171 * t100 + t132 * t97;
t165 = pkin(2) * MDP(14);
t164 = qJ(5) * t106;
t163 = qJ(5) * t179;
t160 = t134 * qJ(5);
t113 = pkin(3) * t162 + t133 * pkin(7);
t158 = t106 * MDP(18);
t107 = t179 * t133;
t157 = t107 * MDP(15);
t156 = t107 * MDP(17);
t155 = t128 * MDP(12);
t154 = t129 * MDP(11);
t153 = MDP(21) - MDP(24);
t152 = MDP(22) + MDP(26);
t151 = -MDP(23) + MDP(28);
t150 = MDP(24) + MDP(27);
t149 = MDP(25) + MDP(29);
t148 = -0.2e1 * t160 + t89;
t124 = t134 * pkin(4);
t86 = t124 + t166;
t121 = -t129 * pkin(3) - pkin(2);
t146 = MDP(20) + t151;
t145 = MDP(21) - t150;
t144 = -pkin(4) * MDP(25) + MDP(23);
t115 = t167 * t128;
t116 = t167 * t129;
t101 = t171 * t115 + t132 * t116;
t85 = t160 - t89;
t143 = -t107 * qJ(5) + t113;
t142 = -t107 * pkin(5) - t86;
t140 = MDP(11) * t128 + MDP(12) * t129;
t139 = -t112 * qJ(5) + t121;
t102 = -t132 * t115 + t171 * t116;
t138 = -t154 + t155 - t165;
t92 = t112 * pkin(5) + t101;
t93 = pkin(5) * t179 + t102;
t137 = -t112 * MDP(17) - MDP(18) * t179 - t93 * MDP(27) + t92 * MDP(28);
t135 = qJ(5) ^ 2;
t127 = t133 ^ 2;
t104 = -t128 * t169 + t109;
t96 = -pkin(4) * t179 + t139;
t91 = -t130 * t179 + t139;
t90 = t106 * pkin(4) + t143;
t87 = t106 * t130 + t143;
t84 = -t85 - t168;
t83 = t134 * qJ(6) - t142;
t1 = [t127 * MDP(4) + (t127 * pkin(7) ^ 2 + t104 ^ 2 + t105 ^ 2) * MDP(14) + (t85 ^ 2 + t86 ^ 2 + t90 ^ 2) * MDP(25) + (t83 ^ 2 + t84 ^ 2 + t87 ^ 2) * MDP(29) + MDP(1) + (t106 * t176 + t157) * t107 + (t134 * MDP(19) + 0.2e1 * pkin(1) * MDP(9) - 0.2e1 * t156 + 0.2e1 * t158) * t134 + 0.2e1 * (-t104 * t134 + t127 * t170) * MDP(11) + 0.2e1 * (t127 * pkin(7) * t129 + t105 * t134) * MDP(12) + 0.2e1 * (t113 * t106 + t134 * t166) * MDP(20) + (t113 * t107 + t89 * t134) * t175 + (t85 * t106 + t86 * t107) * t174 + 0.2e1 * (-t90 * t106 - t86 * t134) * MDP(23) + 0.2e1 * (-t90 * t107 + t85 * t134) * MDP(24) + (-t84 * t106 + t83 * t107) * t173 + 0.2e1 * (-t87 * t107 - t84 * t134) * MDP(27) + (t87 * t106 + t83 * t134) * t172 + (-0.2e1 * pkin(1) * MDP(10) + 0.2e1 * t134 * MDP(5) + (-t104 * t129 - t105 * t128) * t177) * t133; t112 * t157 + (-t112 * t106 + t107 * t179) * MDP(16) + (t121 * t106 - t113 * t179) * MDP(20) + (t121 * t107 + t113 * t112) * MDP(21) + (t101 * t107 - t102 * t106 + t86 * t112 - t179 * t85) * MDP(22) + (-t96 * t106 + t179 * t90) * MDP(23) + (-t96 * t107 - t90 * t112) * MDP(24) + (t86 * t101 - t85 * t102 + t90 * t96) * MDP(25) + (-t93 * t106 + t92 * t107 + t83 * t112 + t179 * t84) * MDP(26) + (-t91 * t107 - t87 * t112) * MDP(27) + (t91 * t106 - t179 * t87) * MDP(28) + (t83 * t92 + t84 * t93 + t87 * t91) * MDP(29) + (-pkin(7) * MDP(10) + MDP(7) + t153 * t102 + (MDP(20) - MDP(23)) * t101 + t140 * qJ(3) + t137) * t134 + (MDP(6) - t140 * pkin(2) + (-MDP(9) + t138) * pkin(7)) * t133 + (MDP(13) + t178) * (-t104 * t128 + t105 * t129); MDP(8) + (t101 ^ 2 + t102 ^ 2 + t96 ^ 2) * MDP(25) + (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) * MDP(29) + (0.2e1 * t154 - 0.2e1 * t155 + t165) * pkin(2) - 0.2e1 * (t121 * MDP(20) - t96 * MDP(23) + t91 * MDP(28)) * t179 + (MDP(15) * t112 - 0.2e1 * t96 * MDP(24) - 0.2e1 * t91 * MDP(27) + t121 * t175 - t176 * t179) * t112 + (t101 * t112 + t102 * t179) * t174 + (t92 * t112 + t179 * t93) * t173 + (t177 + t178) * (t128 ^ 2 + t129 ^ 2) * qJ(3); t90 * MDP(25) + t87 * MDP(29) + (pkin(7) * MDP(14) + t140) * t133 + t145 * t107 + t146 * t106; t96 * MDP(25) + t91 * MDP(29) + t112 * t145 - t146 * t179 + t138; MDP(14) + t149; t156 - t158 - t166 * MDP(20) - t89 * MDP(21) + (-pkin(4) * t107 - t164) * MDP(22) + (0.2e1 * t124 + t166) * MDP(23) + t148 * MDP(24) + (-t86 * pkin(4) - t85 * qJ(5)) * MDP(25) + (-t130 * t107 - t164) * MDP(26) + (t148 - t168) * MDP(27) + t142 * MDP(28) + (t84 * qJ(5) - t83 * t130) * MDP(29) + (-MDP(19) + (-qJ(6) - t130) * MDP(28)) * t134; (-pkin(4) * t112 + t163) * MDP(22) + (-t130 * t112 + t163) * MDP(26) + (t93 * qJ(5) - t92 * t130) * MDP(29) + (qJ(5) * MDP(25) - t153) * t102 + (-MDP(20) + t144) * t101 - t137; 0; MDP(19) - 0.2e1 * pkin(4) * MDP(23) + (pkin(4) ^ 2 + t135) * MDP(25) + t130 * t172 + (t130 ^ 2 + t135) * MDP(29) + 0.2e1 * t150 * qJ(5); t86 * MDP(25) + t83 * MDP(29) + t107 * t152 + t134 * t151; t101 * MDP(25) + t92 * MDP(29) + t112 * t152; 0; -t130 * MDP(29) - MDP(28) + t144; t149; -t106 * MDP(26) - t134 * MDP(27) + t84 * MDP(29); MDP(26) * t179 + t93 * MDP(29); 0; MDP(29) * qJ(5) + MDP(27); 0; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
