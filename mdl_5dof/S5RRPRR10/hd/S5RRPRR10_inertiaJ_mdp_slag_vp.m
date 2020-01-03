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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR10_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:26:41
% EndTime: 2019-12-31 20:26:43
% DurationCPUTime: 0.55s
% Computational Cost: add. (736->151), mult. (1821->229), div. (0->0), fcn. (1964->10), ass. (0->81)
t125 = sin(qJ(5));
t128 = cos(qJ(5));
t136 = MDP(25) * t125 + MDP(26) * t128;
t133 = t125 * MDP(22) + t128 * MDP(23) - t136 * pkin(9);
t169 = 0.2e1 * MDP(25);
t168 = 0.2e1 * MDP(26);
t130 = cos(qJ(2));
t167 = pkin(1) * t130;
t166 = pkin(7) + qJ(3);
t126 = sin(qJ(4));
t121 = sin(pkin(10));
t122 = sin(pkin(5));
t123 = cos(pkin(10));
t127 = sin(qJ(2));
t107 = (t121 * t130 + t123 * t127) * t122;
t124 = cos(pkin(5));
t129 = cos(qJ(4));
t99 = t107 * t126 - t124 * t129;
t165 = t126 * t99;
t157 = t122 * t130;
t158 = t122 * t127;
t106 = t121 * t158 - t123 * t157;
t113 = t124 * t167;
t101 = t124 * pkin(2) - t166 * t158 + t113;
t144 = t124 * t127 * pkin(1);
t104 = t166 * t157 + t144;
t92 = t121 * t101 + t123 * t104;
t90 = t124 * pkin(8) + t92;
t111 = (-pkin(2) * t130 - pkin(1)) * t122;
t95 = t106 * pkin(3) - t107 * pkin(8) + t111;
t87 = -t126 * t90 + t129 * t95;
t84 = -t106 * pkin(4) - t87;
t164 = t84 * t125;
t163 = t84 * t128;
t115 = t121 * pkin(2) + pkin(8);
t162 = t115 * t125;
t161 = t115 * t128;
t160 = t115 * t129;
t117 = t122 ^ 2;
t159 = t117 * t127;
t156 = t124 * MDP(8);
t100 = t107 * t129 + t124 * t126;
t93 = t100 * t125 - t106 * t128;
t155 = t93 * MDP(23);
t94 = t100 * t128 + t106 * t125;
t154 = t94 * MDP(20);
t153 = t94 * MDP(22);
t152 = t99 * MDP(24);
t151 = t100 * MDP(13);
t150 = t100 * MDP(14);
t116 = -t123 * pkin(2) - pkin(3);
t149 = t116 * MDP(18);
t148 = t116 * MDP(19);
t147 = t126 * MDP(19);
t146 = t128 * MDP(20);
t145 = t129 * MDP(24);
t143 = t125 * t128 * MDP(21);
t91 = t123 * t101 - t121 * t104;
t142 = -t115 * MDP(18) + MDP(15);
t141 = -t115 * MDP(19) + MDP(16);
t88 = t126 * t95 + t129 * t90;
t110 = -t129 * pkin(4) - t126 * pkin(9) + t116;
t97 = t128 * t110 - t125 * t160;
t98 = t125 * t110 + t128 * t160;
t140 = t97 * MDP(25) - t98 * MDP(26);
t139 = MDP(22) * t128 - MDP(23) * t125;
t137 = t128 * MDP(25) - t125 * MDP(26);
t89 = -t124 * pkin(3) - t91;
t135 = -MDP(14) + t139;
t134 = (t127 * MDP(6) + t130 * MDP(7)) * t122;
t85 = t106 * pkin(9) + t88;
t86 = t99 * pkin(4) - t100 * pkin(9) + t89;
t82 = -t125 * t85 + t128 * t86;
t83 = t125 * t86 + t128 * t85;
t132 = t82 * MDP(25) - t83 * MDP(26) + t152 + t153 - t155;
t120 = t128 ^ 2;
t119 = t126 ^ 2;
t118 = t125 ^ 2;
t109 = pkin(7) * t157 + t144;
t108 = -pkin(7) * t158 + t113;
t1 = [(t111 ^ 2 + t91 ^ 2 + t92 ^ 2) * MDP(12) + t106 ^ 2 * MDP(17) + MDP(1) + (-0.2e1 * t93 * MDP(21) + t154) * t94 + (MDP(4) * t127 + 0.2e1 * MDP(5) * t130) * t159 + (0.2e1 * t106 * MDP(15) + t151) * t100 + (0.2e1 * t134 + t156) * t124 + (-0.2e1 * t106 * MDP(16) - 0.2e1 * t150 + t152 + 0.2e1 * t153 - 0.2e1 * t155) * t99 + 0.2e1 * (t108 * t124 + t117 * t167) * MDP(9) + 0.2e1 * (-pkin(1) * t159 - t109 * t124) * MDP(10) + 0.2e1 * (-t92 * t106 - t91 * t107) * MDP(11) + 0.2e1 * (t87 * t106 + t89 * t99) * MDP(18) + 0.2e1 * (t89 * t100 - t88 * t106) * MDP(19) + (t82 * t99 + t84 * t93) * t169 + (-t83 * t99 + t84 * t94) * t168; t100 * t148 - t109 * MDP(10) + t156 + t108 * MDP(9) + t134 + (t140 + t149) * t99 + ((-t106 * t121 - t107 * t123) * MDP(11) + (t121 * t92 + t123 * t91) * MDP(12)) * pkin(2) + (-t89 * MDP(18) + t141 * t106 - t132 + t150) * t129 + (t151 + t89 * MDP(19) + t94 * t146 + (-t125 * t94 - t128 * t93) * MDP(21) + (t115 * t93 + t164) * MDP(25) + (t115 * t94 + t163) * MDP(26) + t142 * t106 + t135 * t99) * t126; MDP(8) + (t121 ^ 2 + t123 ^ 2) * MDP(12) * pkin(2) ^ 2 + (t145 - 0.2e1 * t149) * t129 + (t120 * MDP(20) + MDP(13) - 0.2e1 * t143) * t119 + (t119 * t162 - t97 * t129) * t169 + (t119 * t161 + t98 * t129) * t168 + 0.2e1 * (-t129 * t135 + t148) * t126; t111 * MDP(12) + (-t125 * t165 - t129 * t93) * MDP(25) + (-t128 * t165 - t94 * t129) * MDP(26) + (t129 * MDP(18) - t147) * t106; 0; MDP(12); t100 * MDP(15) + t106 * MDP(17) + t87 * MDP(18) - t88 * MDP(19) + t125 * t154 + (-t125 * t93 + t94 * t128) * MDP(21) + (-pkin(4) * t93 - t163) * MDP(25) + (-pkin(4) * t94 + t164) * MDP(26) + (-MDP(16) + t133) * t99; (-t133 + t141) * t129 + (t125 * t146 + (-t118 + t120) * MDP(21) + (-pkin(4) * t125 - t161) * MDP(25) + (-pkin(4) * t128 + t162) * MDP(26) + t142) * t126; -t147 + (MDP(18) + t137) * t129; t118 * MDP(20) + 0.2e1 * pkin(4) * t137 + MDP(17) + 0.2e1 * t143; t132; t139 * t126 + t140 - t145; -t136 * t126; t133; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
