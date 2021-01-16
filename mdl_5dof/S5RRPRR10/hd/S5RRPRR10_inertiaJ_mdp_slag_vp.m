% Calculate joint inertia matrix for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR10_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:02:47
% EndTime: 2021-01-15 22:02:51
% DurationCPUTime: 0.56s
% Computational Cost: add. (778->162), mult. (1927->244), div. (0->0), fcn. (2062->10), ass. (0->82)
t131 = sin(qJ(5));
t134 = cos(qJ(5));
t142 = MDP(27) * t131 + MDP(28) * t134;
t139 = t131 * MDP(24) + t134 * MDP(25) - t142 * pkin(9);
t176 = 0.2e1 * MDP(27);
t175 = 0.2e1 * MDP(28);
t136 = cos(qJ(2));
t174 = pkin(1) * t136;
t173 = pkin(7) + qJ(3);
t127 = sin(pkin(10));
t129 = cos(pkin(10));
t128 = sin(pkin(5));
t163 = t128 * t136;
t133 = sin(qJ(2));
t164 = t128 * t133;
t112 = t127 * t164 - t129 * t163;
t113 = (t127 * t136 + t129 * t133) * t128;
t117 = (-pkin(2) * t136 - pkin(1)) * t128;
t100 = pkin(3) * t112 - pkin(8) * t113 + t117;
t132 = sin(qJ(4));
t135 = cos(qJ(4));
t130 = cos(pkin(5));
t119 = t130 * t174;
t107 = pkin(2) * t130 - t173 * t164 + t119;
t151 = pkin(1) * t130 * t133;
t110 = t173 * t163 + t151;
t97 = t127 * t107 + t129 * t110;
t95 = pkin(8) * t130 + t97;
t92 = t100 * t135 - t132 * t95;
t89 = -pkin(4) * t112 - t92;
t172 = t131 * t89;
t171 = t134 * t89;
t106 = t113 * t135 + t130 * t132;
t99 = t106 * t134 + t112 * t131;
t170 = MDP(22) * t99;
t105 = t113 * t132 - t130 * t135;
t169 = t105 * t132;
t121 = pkin(2) * t127 + pkin(8);
t168 = t121 * t131;
t167 = t121 * t134;
t166 = t121 * t135;
t123 = t128 ^ 2;
t165 = t123 * t133;
t162 = t130 * MDP(8);
t98 = t106 * t131 - t134 * t112;
t161 = t98 * MDP(25);
t160 = t99 * MDP(24);
t96 = t129 * t107 - t127 * t110;
t122 = -pkin(2) * t129 - pkin(3);
t159 = MDP(20) * t122;
t158 = MDP(21) * t132;
t157 = MDP(22) * t134;
t156 = MDP(26) * t135;
t155 = t105 * MDP(26);
t154 = t106 * MDP(15);
t153 = t106 * MDP(16);
t152 = t122 * MDP(21);
t150 = MDP(23) * t131 * t134;
t149 = -t121 * MDP(20) + MDP(17);
t148 = -t121 * MDP(21) + MDP(18);
t94 = -pkin(3) * t130 - t96;
t93 = t100 * t132 + t135 * t95;
t147 = MDP(11) * t129 - MDP(12) * t127;
t146 = MDP(24) * t134 - MDP(25) * t131;
t116 = -pkin(4) * t135 - pkin(9) * t132 + t122;
t103 = t116 * t134 - t131 * t166;
t104 = t116 * t131 + t134 * t166;
t144 = MDP(27) * t103 - MDP(28) * t104;
t143 = MDP(27) * t134 - MDP(28) * t131;
t141 = -MDP(16) + t146;
t140 = (MDP(6) * t133 + MDP(7) * t136) * t128;
t90 = pkin(9) * t112 + t93;
t91 = pkin(4) * t105 - pkin(9) * t106 + t94;
t87 = -t131 * t90 + t134 * t91;
t88 = t131 * t91 + t134 * t90;
t138 = t87 * MDP(27) - t88 * MDP(28) + t155 + t160 - t161;
t126 = t134 ^ 2;
t125 = t132 ^ 2;
t124 = t131 ^ 2;
t115 = pkin(7) * t163 + t151;
t114 = -pkin(7) * t164 + t119;
t1 = [t112 ^ 2 * MDP(19) + (t117 ^ 2 + t96 ^ 2 + t97 ^ 2) * MDP(14) + MDP(1) + (-0.2e1 * MDP(23) * t98 + t170) * t99 + (MDP(4) * t133 + 0.2e1 * MDP(5) * t136) * t165 + (0.2e1 * t112 * MDP(17) + t154) * t106 + (0.2e1 * t140 + t162) * t130 + (-0.2e1 * MDP(18) * t112 - 0.2e1 * t153 + t155 + 0.2e1 * t160 - 0.2e1 * t161) * t105 + (t105 * t87 + t89 * t98) * t176 + (-t105 * t88 + t89 * t99) * t175 + 0.2e1 * (t105 * t94 + t112 * t92) * MDP(20) + 0.2e1 * (t106 * t94 - t112 * t93) * MDP(21) + 0.2e1 * (-t112 * t97 - t113 * t96) * MDP(13) + 0.2e1 * (t112 * t117 + t130 * t96) * MDP(11) + 0.2e1 * (t113 * t117 - t130 * t97) * MDP(12) + 0.2e1 * (-pkin(1) * t165 - t115 * t130) * MDP(10) + 0.2e1 * (t114 * t130 + t123 * t174) * MDP(9); t162 + t114 * MDP(9) - t115 * MDP(10) + MDP(11) * t96 - MDP(12) * t97 + t106 * t152 + t140 + (t144 + t159) * t105 + ((-t112 * t127 - t113 * t129) * MDP(13) + (t127 * t97 + t129 * t96) * MDP(14) + t147 * t130) * pkin(2) + (-t94 * MDP(20) + t148 * t112 - t138 + t153) * t135 + (t154 + t94 * MDP(21) + t99 * t157 + (-t131 * t99 - t134 * t98) * MDP(23) + (t121 * t98 + t172) * MDP(27) + (t121 * t99 + t171) * MDP(28) + t149 * t112 + t141 * t105) * t132; MDP(8) + (t156 - 0.2e1 * t159) * t135 + (MDP(22) * t126 + MDP(15) - 0.2e1 * t150) * t125 + (-t103 * t135 + t125 * t168) * t176 + (t104 * t135 + t125 * t167) * t175 + 0.2e1 * (-t135 * t141 + t152) * t132 + (0.2e1 * t147 + (t127 ^ 2 + t129 ^ 2) * MDP(14) * pkin(2)) * pkin(2); t113 * MDP(12) + t117 * MDP(14) + (-t131 * t169 - t135 * t98) * MDP(27) + (-t134 * t169 - t135 * t99) * MDP(28) + (MDP(20) * t135 + MDP(11) - t158) * t112; 0; MDP(14); t106 * MDP(17) + t112 * MDP(19) + t92 * MDP(20) - t93 * MDP(21) + t131 * t170 + (-t131 * t98 + t134 * t99) * MDP(23) + (-pkin(4) * t98 - t171) * MDP(27) + (-pkin(4) * t99 + t172) * MDP(28) + (-MDP(18) + t139) * t105; (-t139 + t148) * t135 + (t131 * t157 + (-t124 + t126) * MDP(23) + (-pkin(4) * t131 - t167) * MDP(27) + (-pkin(4) * t134 + t168) * MDP(28) + t149) * t132; -t158 + (MDP(20) + t143) * t135; MDP(22) * t124 + 0.2e1 * pkin(4) * t143 + MDP(19) + 0.2e1 * t150; t138; t146 * t132 + t144 - t156; -t142 * t132; t139; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
