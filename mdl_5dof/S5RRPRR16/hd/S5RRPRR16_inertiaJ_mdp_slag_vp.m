% Calculate joint inertia matrix for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR16_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR16_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR16_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:28
% EndTime: 2019-12-31 20:47:30
% DurationCPUTime: 0.69s
% Computational Cost: add. (494->161), mult. (1131->226), div. (0->0), fcn. (1109->8), ass. (0->72)
t114 = sin(pkin(5));
t121 = cos(qJ(2));
t152 = t114 * t121;
t117 = sin(qJ(4));
t120 = cos(qJ(4));
t118 = sin(qJ(2));
t153 = t114 * t118;
t103 = pkin(7) * t153;
t115 = cos(pkin(5));
t158 = pkin(1) * t121;
t136 = -pkin(2) - t158;
t87 = pkin(3) * t153 + t103 + (-pkin(8) + t136) * t115;
t122 = -pkin(2) - pkin(8);
t134 = -qJ(3) * t118 - pkin(1);
t91 = (t122 * t121 + t134) * t114;
t84 = -t117 * t91 + t120 * t87;
t85 = t117 * t87 + t120 * t91;
t98 = t115 * t120 - t117 * t152;
t164 = t98 * MDP(17) + t84 * MDP(20) - t85 * MDP(21);
t163 = 0.2e1 * t115;
t116 = sin(qJ(5));
t119 = cos(qJ(5));
t162 = t116 * MDP(27) + t119 * MDP(28);
t160 = 0.2e1 * MDP(27);
t159 = 0.2e1 * MDP(28);
t157 = pkin(2) * MDP(14);
t97 = t115 * t117 + t120 * t152;
t156 = t117 * t97;
t82 = -pkin(4) * t153 - t84;
t155 = t82 * t116;
t154 = t82 * t119;
t151 = t116 * t122;
t150 = t119 * t122;
t88 = t98 * t116 - t119 * t153;
t149 = t88 * MDP(25);
t89 = t116 * t153 + t98 * t119;
t148 = t89 * MDP(22);
t147 = t89 * MDP(24);
t146 = t97 * MDP(26);
t145 = t98 * MDP(16);
t100 = t115 * t118 * pkin(1) + pkin(7) * t152;
t142 = qJ(3) * MDP(21);
t140 = t117 * MDP(21);
t139 = t119 * MDP(22);
t137 = t122 * MDP(21);
t107 = t115 * qJ(3);
t94 = -t107 - t100;
t135 = t116 * t119 * MDP(23);
t133 = t122 * MDP(20) + MDP(17);
t90 = pkin(3) * t152 - t94;
t132 = -t100 * MDP(10) + (t115 * t158 - t103) * MDP(9);
t101 = t117 * pkin(4) - t120 * pkin(9) + qJ(3);
t92 = t119 * t101 - t117 * t151;
t93 = t116 * t101 + t117 * t150;
t130 = t92 * MDP(27) - t93 * MDP(28);
t129 = MDP(24) * t119 - MDP(25) * t116;
t128 = t119 * MDP(27) - t116 * MDP(28);
t126 = -MDP(16) + t129;
t125 = t116 * MDP(24) + t119 * MDP(25) - pkin(9) * t162;
t83 = pkin(9) * t153 + t85;
t86 = t97 * pkin(4) - t98 * pkin(9) + t90;
t80 = -t116 * t83 + t119 * t86;
t81 = t116 * t86 + t119 * t83;
t124 = t80 * MDP(27) - t81 * MDP(28) + t146 + t147 - t149;
t123 = -MDP(18) + t125;
t113 = t120 ^ 2;
t112 = t119 ^ 2;
t110 = t117 ^ 2;
t109 = t116 ^ 2;
t96 = t136 * t115 + t103;
t95 = (-pkin(2) * t121 + t134) * t114;
t1 = [t98 ^ 2 * MDP(15) + (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) * MDP(14) + MDP(1) + t115 ^ 2 * MDP(8) + (-0.2e1 * t88 * MDP(23) + t148) * t89 + (-0.2e1 * MDP(18) * t153 - 0.2e1 * t145 + t146 + 0.2e1 * t147 - 0.2e1 * t149) * t97 + (t80 * t97 + t82 * t88) * t160 + (-t81 * t97 + t82 * t89) * t159 + 0.2e1 * (t97 * MDP(20) + t98 * MDP(21)) * t90 + (t96 * MDP(12) - t94 * MDP(13) + t132) * t163 + ((t118 * MDP(6) + t121 * MDP(7)) * t163 + 0.2e1 * (-t94 * MDP(11) + t95 * MDP(12)) * t121 + ((MDP(19) + MDP(4)) * t118 ^ 2 + 0.2e1 * (-MDP(10) * t118 + MDP(9) * t121) * pkin(1)) * t114 + 0.2e1 * (t96 * MDP(11) - t95 * MDP(13) + MDP(5) * t152 + t164) * t118) * t114; t103 * MDP(12) + (0.2e1 * t107 + t100) * MDP(13) + (-t96 * pkin(2) - t94 * qJ(3)) * MDP(14) + t98 * t142 + (MDP(8) + (-0.2e1 * pkin(2) - t158) * MDP(12)) * t115 + (qJ(3) * MDP(20) + t130) * t97 + (t90 * MDP(20) + t124 - t145) * t117 + ((qJ(3) * MDP(11) + MDP(7)) * t121 + (-pkin(2) * MDP(11) + MDP(6) + (-MDP(18) - t137) * t117) * t118) * t114 + (t98 * MDP(15) + t90 * MDP(21) + t89 * t139 + (-t116 * t89 - t119 * t88) * MDP(23) + (-t122 * t88 + t155) * MDP(27) + (-t122 * t89 + t154) * MDP(28) + t133 * t153 + t126 * t97) * t120 + t132; t110 * MDP(26) + MDP(8) + (-0.2e1 * MDP(12) + t157) * pkin(2) + (MDP(14) * qJ(3) + 0.2e1 * t117 * MDP(20) + 0.2e1 * MDP(13)) * qJ(3) + (t112 * MDP(22) + MDP(15) - 0.2e1 * t135) * t113 + (-t113 * t151 + t92 * t117) * t160 + (-t113 * t150 - t93 * t117) * t159 + 0.2e1 * (t117 * t126 + t142) * t120; t115 * MDP(12) + t96 * MDP(14) + (-t116 * t156 - t120 * t88) * MDP(27) + (-t119 * t156 - t120 * t89) * MDP(28) + (t120 * MDP(20) + MDP(11) - t140) * t153; MDP(12) - t157 + t162 * (-t110 - t113); MDP(14); MDP(19) * t153 + t116 * t148 + (-t116 * t88 + t89 * t119) * MDP(23) + (-pkin(4) * t88 - t154) * MDP(27) + (-pkin(4) * t89 + t155) * MDP(28) + t123 * t97 + t164; (t123 - t137) * t117 + (t116 * t139 + (-t109 + t112) * MDP(23) + (-pkin(4) * t116 + t150) * MDP(27) + (-pkin(4) * t119 - t151) * MDP(28) + t133) * t120; -t140 + (MDP(20) + t128) * t120; t109 * MDP(22) + 0.2e1 * pkin(4) * t128 + MDP(19) + 0.2e1 * t135; t124; t117 * MDP(26) + t129 * t120 + t130; -t162 * t117; t125; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
