% Calculate joint inertia matrix for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP10_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:35:56
% EndTime: 2021-01-16 00:36:00
% DurationCPUTime: 0.84s
% Computational Cost: add. (1011->213), mult. (2279->303), div. (0->0), fcn. (2341->8), ass. (0->87)
t122 = sin(pkin(5));
t177 = 0.2e1 * t122;
t125 = sin(qJ(3));
t176 = -0.2e1 * t125;
t124 = sin(qJ(4));
t168 = qJ(5) + pkin(9);
t111 = t168 * t124;
t127 = cos(qJ(4));
t139 = MDP(23) * t124 + MDP(24) * t127;
t112 = t168 * t127;
t150 = t112 * MDP(26);
t175 = t124 * MDP(20) + t127 * MDP(21) - t111 * MDP(25) - t139 * pkin(9) - MDP(14) - t150;
t174 = 0.2e1 * MDP(25);
t173 = 2 * MDP(27);
t126 = sin(qJ(2));
t172 = pkin(1) * t126;
t129 = cos(qJ(2));
t171 = pkin(1) * t129;
t170 = pkin(8) * t124;
t128 = cos(qJ(3));
t169 = pkin(8) * t128;
t167 = MDP(28) * pkin(4);
t166 = pkin(2) * MDP(16);
t165 = pkin(2) * MDP(17);
t164 = pkin(8) * MDP(17);
t163 = qJ(5) * t125;
t162 = t122 * t126;
t161 = t122 * t129;
t123 = cos(pkin(5));
t160 = t123 * MDP(8);
t159 = t126 * MDP(6);
t105 = t123 * t125 + t128 * t162;
t95 = t105 * t124 + t127 * t161;
t158 = t95 * MDP(19);
t157 = t95 * MDP(21);
t96 = t105 * t127 - t124 * t161;
t156 = t96 * MDP(18);
t155 = t96 * MDP(20);
t154 = MDP(15) * t129;
t104 = -t123 * t128 + t125 * t162;
t153 = t104 * MDP(22);
t152 = t105 * MDP(12);
t151 = t105 * MDP(13);
t117 = -t127 * pkin(4) - pkin(3);
t149 = t117 * MDP(28);
t148 = t124 * MDP(18);
t147 = t124 * MDP(21);
t146 = t124 * MDP(26);
t145 = t127 * MDP(25);
t144 = t127 * t169;
t143 = pkin(7) * t161;
t142 = t124 * t127 * MDP(19);
t114 = pkin(7) * t162;
t101 = t114 + (-pkin(2) - t171) * t123;
t89 = t104 * pkin(3) - t105 * pkin(9) + t101;
t102 = t143 + (pkin(8) + t172) * t123;
t103 = (-pkin(2) * t129 - pkin(8) * t126 - pkin(1)) * t122;
t93 = t128 * t102 + t125 * t103;
t91 = -pkin(9) * t161 + t93;
t85 = -t124 * t91 + t127 * t89;
t92 = -t125 * t102 + t128 * t103;
t141 = pkin(8) * MDP(16) - MDP(13);
t140 = -pkin(4) * MDP(27) + MDP(20);
t90 = pkin(3) * t161 - t92;
t86 = t124 * t89 + t127 * t91;
t138 = t124 * MDP(25) + t127 * MDP(26);
t137 = -t96 * qJ(5) + t85;
t136 = t127 * MDP(20) - MDP(12) - t147;
t87 = t95 * pkin(4) + t90;
t135 = t95 * MDP(25) + t96 * MDP(26) + t87 * MDP(28);
t110 = -t128 * pkin(3) - t125 * pkin(9) - pkin(2);
t108 = t127 * t110;
t97 = t144 + (t110 - t163) * t124;
t134 = (-t124 * t169 + t108) * MDP(23) - (t124 * t110 + t144) * MDP(24) - t97 * MDP(26);
t133 = -t145 + t146 + t149;
t94 = -t127 * t163 + t108 + (-pkin(4) - t170) * t128;
t132 = -t94 * MDP(25) - t134;
t84 = -t95 * qJ(5) + t86;
t130 = t85 * MDP(23) - t86 * MDP(24) - t84 * MDP(26) + t153 + t155 - t157;
t121 = t127 ^ 2;
t119 = t124 ^ 2;
t118 = t122 ^ 2;
t109 = (pkin(4) * t124 + pkin(8)) * t125;
t107 = t123 * t172 + t143;
t106 = t123 * t171 - t114;
t83 = t104 * pkin(4) + t137;
t1 = [(t83 ^ 2 + t84 ^ 2 + t87 ^ 2) * MDP(28) + t105 ^ 2 * MDP(11) + MDP(1) + t118 * t126 ^ 2 * MDP(4) + (t156 - 0.2e1 * t158) * t96 + (t159 * t177 + t160) * t123 + ((MDP(7) * t123 - t151) * t177 + (0.2e1 * MDP(5) * t126 + t154) * t118) * t129 + (0.2e1 * MDP(14) * t161 - 0.2e1 * t152 + t153 + 0.2e1 * t155 - 0.2e1 * t157) * t104 + (-t83 * t96 - t84 * t95) * t173 + (t83 * t104 + t87 * t95) * t174 + 0.2e1 * (t85 * t104 + t90 * t95) * MDP(23) + 0.2e1 * (-t84 * t104 + t87 * t96) * MDP(26) + 0.2e1 * (-t86 * t104 + t90 * t96) * MDP(24) + 0.2e1 * (-t107 * t123 - t118 * t172) * MDP(10) + 0.2e1 * (t101 * t104 - t92 * t161) * MDP(16) + 0.2e1 * (t101 * t105 + t93 * t161) * MDP(17) + 0.2e1 * (t106 * t123 + t118 * t171) * MDP(9); t160 + t106 * MDP(9) - t107 * MDP(10) - t105 * t165 + (-t94 * t96 - t97 * t95) * MDP(27) + (t83 * t94 + t84 * t97) * MDP(28) + (t129 * MDP(7) + t159) * t122 + t135 * t109 + (-t132 - t166) * t104 + (t152 - t101 * MDP(16) - t83 * MDP(25) + (-MDP(14) + t164) * t161 - t130) * t128 + (t105 * MDP(11) + t101 * MDP(17) + (t95 * MDP(23) + t96 * MDP(24)) * pkin(8) + (-t96 * MDP(19) + t90 * MDP(23) + t87 * MDP(25) - t84 * MDP(27)) * t124 + (t90 * MDP(24) + t87 * MDP(26) - t83 * MDP(27) + t156 - t158) * t127 + t141 * t161 + t136 * t104) * t125; MDP(8) + t165 * t176 + (t109 ^ 2 + t94 ^ 2 + t97 ^ 2) * MDP(28) + 0.2e1 * ((-t124 * t97 - t127 * t94) * MDP(27) + t138 * t109) * t125 + (t128 * MDP(22) + t136 * t176 + 0.2e1 * t132 + 0.2e1 * t166) * t128 + (t121 * MDP(18) + 0.2e1 * t139 * pkin(8) + MDP(11) - 0.2e1 * t142) * t125 ^ 2; t151 - t122 * t154 + t92 * MDP(16) - t93 * MDP(17) + t96 * t148 + (-t124 * t95 + t96 * t127) * MDP(19) + (-pkin(3) * t95 - t90 * t127) * MDP(23) + (-pkin(3) * t96 + t90 * t124) * MDP(24) + (t117 * t95 - t87 * t127) * MDP(25) + (t117 * t96 + t87 * t124) * MDP(26) + (t111 * t96 - t112 * t95 - t83 * t124 + t84 * t127) * MDP(27) + (-t83 * t111 + t84 * t112 + t87 * t117) * MDP(28) + t175 * t104; (-t94 * t124 + t97 * t127) * MDP(27) + (-t94 * t111 + t97 * t112) * MDP(28) + t133 * t109 + (-t164 - t175) * t128 + (t127 * t148 + (-t119 + t121) * MDP(19) + (-pkin(3) * t124 - pkin(8) * t127) * MDP(23) + (-pkin(3) * t127 + t170) * MDP(24) + (t111 * t127 - t112 * t124) * MDP(27) + t138 * t117 - t141) * t125; MDP(15) + t119 * MDP(18) + 0.2e1 * t142 + (t111 * t124 + t112 * t127) * t173 + (t111 ^ 2 + t112 ^ 2) * MDP(28) + (-0.2e1 * t145 + 0.2e1 * t146 + t149) * t117 + 0.2e1 * (t127 * MDP(23) - t124 * MDP(24)) * pkin(3); t137 * MDP(25) + (-t96 * MDP(27) + t83 * MDP(28) + t104 * t174) * pkin(4) + t130; t94 * t167 + t108 * MDP(25) + (-MDP(22) + (-0.2e1 * pkin(4) - t170) * MDP(25)) * t128 + (-t147 + (-qJ(5) * MDP(25) + t140) * t127) * t125 + t134; -t150 + (-MDP(24) * pkin(9) + MDP(21)) * t127 - (MDP(25) + t167) * t111 + (-MDP(23) * pkin(9) + t140) * t124; MDP(22) + (t174 + t167) * pkin(4); t135; t109 * MDP(28) + t138 * t125; t133; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
