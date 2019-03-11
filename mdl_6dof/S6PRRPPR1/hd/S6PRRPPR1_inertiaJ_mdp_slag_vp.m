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
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPPR1_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:02:14
% EndTime: 2019-03-08 21:02:15
% DurationCPUTime: 0.48s
% Computational Cost: add. (752->142), mult. (1552->220), div. (0->0), fcn. (1798->12), ass. (0->79)
t178 = qJ(4) + pkin(8);
t125 = sin(pkin(12));
t128 = cos(pkin(12));
t160 = t125 ^ 2 + t128 ^ 2;
t177 = t160 * MDP(16);
t131 = sin(qJ(3));
t127 = sin(pkin(6));
t132 = sin(qJ(2));
t164 = t127 * t132;
t168 = cos(pkin(6));
t172 = cos(qJ(3));
t107 = t168 * t131 + t172 * t164;
t126 = sin(pkin(11));
t129 = cos(pkin(11));
t139 = -t131 * t164 + t168 * t172;
t94 = t126 * t107 - t129 * t139;
t176 = t94 ^ 2;
t116 = t178 * t172;
t151 = t178 * t131;
t103 = t126 * t116 + t129 * t151;
t175 = t103 ^ 2;
t170 = t129 * pkin(3);
t121 = -pkin(4) - t170;
t115 = -t128 * pkin(5) + t121;
t174 = 0.2e1 * t115;
t173 = -2 * MDP(19);
t171 = t126 * pkin(3);
t118 = qJ(5) + t171;
t169 = pkin(9) + t118;
t110 = t126 * t131 - t129 * t172;
t112 = t126 * t172 + t129 * t131;
t122 = -t172 * pkin(3) - pkin(2);
t101 = t110 * pkin(4) - t112 * qJ(5) + t122;
t105 = t129 * t116 - t126 * t151;
t88 = t125 * t101 + t128 * t105;
t96 = t129 * t107 + t126 * t139;
t167 = MDP(13) * t96;
t166 = t112 * t125;
t165 = t112 * t128;
t134 = cos(qJ(2));
t163 = t127 * t134;
t130 = sin(qJ(6));
t133 = cos(qJ(6));
t113 = t133 * t125 + t130 * t128;
t92 = t113 * t112;
t162 = t92 * MDP(21);
t111 = t130 * t125 - t133 * t128;
t93 = t111 * t112;
t161 = t93 * MDP(20);
t159 = MDP(18) * t113;
t158 = t110 * MDP(22);
t157 = t111 * MDP(23);
t156 = t121 * MDP(17);
t155 = t122 * MDP(13);
t154 = t125 * MDP(15);
t153 = t128 * MDP(14);
t152 = MDP(10) * t172;
t87 = t128 * t101 - t125 * t105;
t150 = t160 * MDP(17);
t149 = t88 * t125 + t87 * t128;
t148 = -t87 * t125 + t88 * t128;
t89 = -t125 * t96 - t128 * t163;
t90 = -t125 * t163 + t128 * t96;
t147 = t90 * t125 + t89 * t128;
t85 = t110 * pkin(5) - pkin(9) * t165 + t87;
t86 = -pkin(9) * t166 + t88;
t145 = (-t130 * t86 + t133 * t85) * MDP(23) - (t130 * t85 + t133 * t86) * MDP(24);
t144 = (-t130 * t90 + t133 * t89) * MDP(23) - (t130 * t89 + t133 * t90) * MDP(24);
t143 = t92 * MDP(23) - t93 * MDP(24);
t142 = t125 * MDP(14) + t128 * MDP(15);
t141 = -t113 * MDP(24) - t157;
t140 = MDP(12) + t142;
t108 = t169 * t125;
t109 = t169 * t128;
t138 = t113 * MDP(20) - t111 * MDP(21) + (-t133 * t108 - t130 * t109) * MDP(23) - (-t130 * t108 + t133 * t109) * MDP(24);
t137 = t141 + t153 - t154;
t136 = -t137 + t156;
t91 = pkin(5) * t166 + t103;
t1 = [MDP(1) + (t127 ^ 2 * t134 ^ 2 + t96 ^ 2 + t176) * MDP(13) + (t89 ^ 2 + t90 ^ 2 + t176) * MDP(17); t105 * t167 + (t89 * t87 + t90 * t88) * MDP(17) + ((MDP(13) + MDP(17)) * t103 + t143) * t94 + (-t96 * MDP(12) + t89 * MDP(14) - t90 * MDP(15) + t144) * t110 + (-t147 * MDP(16) + t140 * t94) * t112 + (-t132 * MDP(4) + (-MDP(11) * t131 + MDP(3) + t152 - t155) * t134) * t127; MDP(2) + 0.2e1 * pkin(2) * t152 + (t105 ^ 2 + t122 ^ 2 + t175) * MDP(13) + (t87 ^ 2 + t88 ^ 2 + t175) * MDP(17) - (-t93 * MDP(18) + t92 * t173) * t93 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t131 + 0.2e1 * t172 * MDP(6)) * t131 + (t158 - 0.2e1 * t161 - 0.2e1 * t162) * t110 + 0.2e1 * t143 * t91 + 0.2e1 * (-t105 * MDP(12) + t87 * MDP(14) - t88 * MDP(15) + t145) * t110 + 0.2e1 * (-t149 * MDP(16) + t140 * t103) * t112; t139 * MDP(10) - t107 * MDP(11) + t167 * t171 + (-MDP(13) * t170 + t136) * t94 + (MDP(17) * t118 + MDP(16)) * (-t89 * t125 + t90 * t128); t131 * MDP(7) + t172 * MDP(8) + (-t103 * t128 + t121 * t166) * MDP(14) + (t103 * t125 + t121 * t165) * MDP(15) + t148 * MDP(16) + (t103 * t121 + t148 * t118) * MDP(17) - t93 * t159 + (t93 * t111 - t113 * t92) * MDP(19) + (t91 * t111 + t115 * t92) * MDP(23) + (t91 * t113 - t115 * t93) * MDP(24) + (-t142 * t118 + t138) * t110 + (-t131 * MDP(10) - t172 * MDP(11)) * pkin(8) + ((-t110 * t126 - t112 * t129) * MDP(12) + (-t103 * t129 + t105 * t126) * MDP(13)) * pkin(3); t157 * t174 + MDP(9) + (t126 ^ 2 + t129 ^ 2) * MDP(13) * pkin(3) ^ 2 + (-0.2e1 * t153 + 0.2e1 * t154 + t156) * t121 + (MDP(24) * t174 + t111 * t173 + t159) * t113 + (t150 * t118 + 0.2e1 * t177) * t118; -MDP(13) * t163 + t147 * MDP(17); t149 * MDP(17) + t137 * t110 - t112 * t177 + t155; 0; MDP(13) + t150; t94 * MDP(17); t103 * MDP(17) + t142 * t112 + t143; t136; 0; MDP(17); t144; t145 + t158 - t161 - t162; t138; t141; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
