% Calculate joint inertia matrix for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP7_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:46
% EndTime: 2019-12-31 21:57:47
% DurationCPUTime: 0.49s
% Computational Cost: add. (634->136), mult. (1143->183), div. (0->0), fcn. (1121->6), ass. (0->71)
t126 = sin(qJ(3));
t114 = t126 * pkin(2) + pkin(8);
t125 = sin(qJ(4));
t123 = t125 ^ 2;
t128 = cos(qJ(4));
t124 = t128 ^ 2;
t157 = t123 + t124;
t159 = t157 * t114;
t176 = t125 * MDP(20) + t128 * MDP(21);
t127 = sin(qJ(2));
t169 = cos(qJ(3));
t170 = cos(qJ(2));
t103 = t126 * t127 - t169 * t170;
t175 = 0.2e1 * t103;
t174 = pkin(7) + pkin(6);
t173 = t125 * MDP(24);
t116 = -t170 * pkin(2) - pkin(1);
t172 = 0.2e1 * t116;
t171 = 2 * MDP(26);
t168 = pkin(8) * t103;
t167 = t103 * pkin(4);
t150 = t169 * pkin(2);
t115 = -t150 - pkin(3);
t166 = pkin(3) - t115;
t104 = t126 * t170 + t169 * t127;
t89 = t103 * pkin(3) - t104 * pkin(8) + t116;
t107 = t174 * t127;
t108 = t174 * t170;
t95 = -t126 * t107 + t169 * t108;
t84 = t125 * t89 + t128 * t95;
t165 = t103 * qJ(5);
t164 = t103 * t114;
t163 = t125 * t128;
t146 = t125 * t95 - t128 * t89;
t82 = t146 - t167;
t162 = t82 * MDP(28);
t141 = pkin(4) * t125 - t128 * qJ(5);
t94 = t169 * t107 + t126 * t108;
t85 = t141 * t104 + t94;
t161 = t85 * MDP(27);
t142 = -t128 * pkin(4) - t125 * qJ(5);
t106 = -pkin(3) + t142;
t101 = -t150 + t106;
t160 = -t101 - t106;
t158 = t157 * pkin(8);
t156 = MDP(28) * t106;
t155 = MDP(28) * t114;
t154 = MDP(28) * t125;
t153 = t101 * MDP(28);
t152 = t103 * MDP(22);
t151 = 0.2e1 * t170;
t149 = -t141 * MDP(26) + t176;
t148 = MDP(19) * t163;
t147 = t123 * MDP(18) + MDP(15) + 0.2e1 * t148;
t145 = -MDP(28) * pkin(4) - MDP(25);
t144 = t157 * MDP(28);
t143 = -pkin(3) * t104 - t168;
t140 = -t104 * t106 + t168;
t139 = -MDP(23) * t146 - t84 * MDP(24);
t138 = -t94 * MDP(23) - t85 * MDP(25);
t137 = t101 * t104 - t164;
t136 = t104 * t115 - t164;
t135 = t128 * MDP(20) - t125 * MDP(21);
t134 = t128 * MDP(23) - t173;
t133 = -0.2e1 * t128 * MDP(25) - 0.2e1 * t125 * MDP(27);
t132 = (t169 * MDP(16) - t126 * MDP(17)) * pkin(2);
t81 = t165 + t84;
t131 = (t82 * t125 + t81 * t128) * MDP(26) - t95 * MDP(17) + (t173 - MDP(16)) * t94 + (MDP(13) + (-t123 + t124) * MDP(19) + MDP(18) * t163) * t104 + (-MDP(14) + t176) * t103;
t130 = (MDP(28) * qJ(5) - MDP(24) + MDP(27)) * t128 + (-MDP(23) + t145) * t125;
t118 = t125 * MDP(26);
t1 = [MDP(1) + pkin(1) * MDP(9) * t151 + (t81 ^ 2 + t82 ^ 2 + t85 ^ 2) * MDP(28) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t127 + MDP(5) * t151) * t127 + (MDP(16) * t172 + t152) * t103 + (-t82 * MDP(25) + t81 * MDP(27) + t139) * t175 + (MDP(17) * t172 + (-MDP(12) + t135) * t175 + 0.2e1 * (t94 * MDP(24) + t82 * MDP(26) - t161) * t128 + 0.2e1 * (-t81 * MDP(26) - t138) * t125 + (t124 * MDP(18) + MDP(11) - 0.2e1 * t148) * t104) * t104; (-t170 * MDP(10) - t127 * MDP(9)) * pkin(6) + (t136 * MDP(23) + t137 * MDP(25) + t82 * t155 - t161) * t125 + (t136 * MDP(24) - t137 * MDP(27) + t81 * t155 + t138) * t128 + t131 + t170 * MDP(7) + t127 * MDP(6) + t85 * t153; MDP(8) + t159 * t171 + t114 ^ 2 * t144 + (t133 + t153) * t101 + t147 - 0.2e1 * t134 * t115 + 0.2e1 * t132; t85 * t156 + (pkin(8) * t81 * MDP(28) + t143 * MDP(24) + t140 * MDP(27) + t138) * t128 + (t143 * MDP(23) - t140 * MDP(25) + pkin(8) * t162 - t161) * t125 + t131; (t158 + t159) * MDP(26) + (pkin(8) * t159 + t101 * t106) * MDP(28) + t132 + (t166 * MDP(23) + t160 * MDP(25)) * t128 + (-t166 * MDP(24) + t160 * MDP(27)) * t125 + t147; t158 * t171 + pkin(8) ^ 2 * t144 + (t133 + t156) * t106 + 0.2e1 * t134 * pkin(3) + t147; t152 + (-t146 + 0.2e1 * t167) * MDP(25) + (0.2e1 * t165 + t84) * MDP(27) + (-t82 * pkin(4) + t81 * qJ(5)) * MDP(28) + (t142 * MDP(26) + t135) * t104 + t139; t130 * t114 + t149; t130 * pkin(8) + t149; MDP(22) + 0.2e1 * pkin(4) * MDP(25) + 0.2e1 * qJ(5) * MDP(27) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(28); t128 * t104 * MDP(26) - t103 * MDP(25) + t162; t114 * t154 + t118; pkin(8) * t154 + t118; t145; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
