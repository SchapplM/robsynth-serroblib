% Calculate joint inertia matrix for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRPRRR4_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:38:57
% EndTime: 2019-03-08 20:38:58
% DurationCPUTime: 0.42s
% Computational Cost: add. (624->133), mult. (1309->202), div. (0->0), fcn. (1524->12), ass. (0->81)
t181 = MDP(8) * qJ(3);
t138 = sin(qJ(5));
t142 = cos(qJ(5));
t150 = -t138 * MDP(21) - t142 * MDP(22);
t137 = sin(qJ(6));
t141 = cos(qJ(6));
t119 = t137 * t138 - t141 * t142;
t120 = t137 * t142 + t141 * t138;
t176 = pkin(9) + pkin(10);
t123 = t176 * t138;
t124 = t176 * t142;
t154 = t120 * MDP(25) - t119 * MDP(26) + (-t141 * t123 - t137 * t124) * MDP(28) - (-t137 * t123 + t141 * t124) * MDP(29);
t180 = t138 * MDP(18) + t142 * MDP(19) + t150 * pkin(9) + t154;
t135 = cos(pkin(12));
t126 = -t135 * pkin(3) - pkin(2);
t179 = 0.2e1 * t126;
t178 = -2 * MDP(24);
t177 = 0.2e1 * MDP(29);
t175 = cos(qJ(4));
t174 = pkin(2) * MDP(8);
t133 = sin(pkin(12));
t139 = sin(qJ(4));
t117 = t139 * t133 - t175 * t135;
t173 = t117 * pkin(5);
t172 = pkin(8) + qJ(3);
t134 = sin(pkin(6));
t143 = cos(qJ(2));
t165 = t134 * t143;
t136 = cos(pkin(6));
t140 = sin(qJ(2));
t166 = t134 * t140;
t108 = -t133 * t166 + t136 * t135;
t109 = t136 * t133 + t135 * t166;
t97 = t139 * t108 + t175 * t109;
t91 = -t138 * t97 - t142 * t165;
t92 = -t138 * t165 + t142 * t97;
t84 = -t137 * t92 + t141 * t91;
t85 = t137 * t91 + t141 * t92;
t171 = t84 * MDP(28) - t85 * MDP(29);
t118 = t175 * t133 + t139 * t135;
t101 = t117 * pkin(4) - t118 * pkin(9) + t126;
t121 = t172 * t133;
t122 = t172 * t135;
t103 = -t139 * t121 + t175 * t122;
t162 = t142 * t103;
t88 = t162 + (-pkin(10) * t118 + t101) * t138;
t170 = t141 * t88;
t169 = t118 * t138;
t168 = t118 * t142;
t167 = t133 * MDP(6);
t164 = t135 * MDP(5);
t163 = t138 * t142;
t98 = t120 * t118;
t94 = t98 * MDP(26);
t99 = t119 * t118;
t95 = t99 * MDP(25);
t111 = t119 * MDP(28);
t161 = -t120 * MDP(29) - t111;
t159 = t118 * MDP(15);
t158 = t120 * MDP(23);
t157 = MDP(20) + MDP(27);
t156 = t117 * MDP(27) - t94 - t95;
t155 = MDP(17) * t163;
t89 = t142 * t101 - t138 * t103;
t87 = -pkin(10) * t168 + t173 + t89;
t80 = -t137 * t88 + t141 * t87;
t152 = t142 * MDP(18) - t138 * MDP(19);
t151 = t142 * MDP(21) - t138 * MDP(22);
t102 = t175 * t121 + t139 * t122;
t149 = -MDP(14) - t151;
t148 = (MDP(28) * t141 - MDP(29) * t137) * pkin(5);
t147 = t159 - t164 + t167 - t174;
t146 = -t149 + t161;
t132 = t142 ^ 2;
t131 = t138 ^ 2;
t128 = -t142 * pkin(5) - pkin(4);
t96 = -t175 * t108 + t139 * t109;
t93 = pkin(5) * t169 + t102;
t90 = t138 * t101 + t162;
t81 = t137 * t87 + t170;
t1 = [MDP(1) + (t134 ^ 2 * t143 ^ 2 + t108 ^ 2 + t109 ^ 2) * MDP(8); (t91 * t117 + t96 * t169) * MDP(21) + (-t92 * t117 + t96 * t168) * MDP(22) + (t84 * t117 + t96 * t98) * MDP(28) + (-t85 * t117 - t96 * t99) * MDP(29) + (-t140 * MDP(4) + (-t117 * MDP(14) + MDP(3) - t147) * t143) * t134 + (MDP(7) + t181) * (-t108 * t133 + t109 * t135); t159 * t179 + MDP(2) - (-t99 * MDP(23) + t98 * t178) * t99 + t157 * t117 ^ 2 + (0.2e1 * t164 - 0.2e1 * t167 + t174) * pkin(2) + (t132 * MDP(16) + MDP(9) - 0.2e1 * t155) * t118 ^ 2 + (MDP(14) * t179 - 0.2e1 * t95 - 0.2e1 * t94 + 0.2e1 * (-MDP(10) + t152) * t118) * t117 + 0.2e1 * (t102 * t169 + t89 * t117) * MDP(21) + 0.2e1 * (t102 * t168 - t90 * t117) * MDP(22) + 0.2e1 * (t80 * t117 + t93 * t98) * MDP(28) + (-t81 * t117 - t93 * t99) * t177 + (0.2e1 * MDP(7) + t181) * (t133 ^ 2 + t135 ^ 2) * qJ(3); -MDP(8) * t165; t146 * t117 + t147; MDP(8); -t97 * MDP(15) - t146 * t96; -t103 * MDP(15) - t99 * t158 + (t99 * t119 - t120 * t98) * MDP(24) + (t93 * t119 + t128 * t98) * MDP(28) + (t93 * t120 - t128 * t99) * MDP(29) + t149 * t102 + (MDP(11) + MDP(16) * t163 + (-t131 + t132) * MDP(17) + t150 * pkin(4)) * t118 + (-MDP(12) + t180) * t117; 0; 0.2e1 * t155 + 0.2e1 * t128 * t111 + t131 * MDP(16) + MDP(13) + 0.2e1 * t151 * pkin(4) + (t119 * t178 + t128 * t177 + t158) * t120; t91 * MDP(21) - t92 * MDP(22) + t171; t117 * MDP(20) + t89 * MDP(21) - t90 * MDP(22) + (t141 * t173 + t80) * MDP(28) + (-t170 + (-t87 - t173) * t137) * MDP(29) + t152 * t118 + t156; t151 + t161; t180; 0.2e1 * t148 + t157; t171; t80 * MDP(28) - t81 * MDP(29) + t156; t161; t154; MDP(27) + t148; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
