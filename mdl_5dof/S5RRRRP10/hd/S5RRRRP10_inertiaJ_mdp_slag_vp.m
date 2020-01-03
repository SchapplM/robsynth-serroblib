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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP10_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:12:23
% EndTime: 2019-12-31 22:12:25
% DurationCPUTime: 0.63s
% Computational Cost: add. (745->181), mult. (1702->267), div. (0->0), fcn. (1735->8), ass. (0->83)
t118 = sin(pkin(5));
t167 = 0.2e1 * t118;
t120 = sin(qJ(4));
t123 = cos(qJ(4));
t166 = t120 * MDP(20) + t123 * MDP(21) - (MDP(23) * t120 + MDP(24) * t123) * pkin(9) - MDP(14);
t165 = 0.2e1 * MDP(23);
t164 = 0.2e1 * MDP(24);
t163 = 2 * MDP(25);
t122 = sin(qJ(2));
t162 = pkin(1) * t122;
t125 = cos(qJ(2));
t161 = pkin(1) * t125;
t160 = pkin(8) * t120;
t159 = pkin(8) * t123;
t124 = cos(qJ(3));
t158 = pkin(8) * t124;
t157 = -qJ(5) - pkin(9);
t156 = MDP(26) * pkin(4);
t155 = pkin(2) * MDP(16);
t154 = pkin(2) * MDP(17);
t153 = pkin(8) * MDP(17);
t148 = t118 * t125;
t121 = sin(qJ(3));
t119 = cos(pkin(5));
t134 = pkin(7) * t148;
t98 = t134 + (pkin(8) + t162) * t119;
t99 = (-pkin(2) * t125 - pkin(8) * t122 - pkin(1)) * t118;
t88 = -t121 * t98 + t124 * t99;
t86 = pkin(3) * t148 - t88;
t152 = t86 * t120;
t151 = t86 * t123;
t150 = qJ(5) * t121;
t149 = t118 * t122;
t147 = t119 * MDP(8);
t146 = t122 * MDP(6);
t101 = t119 * t121 + t124 * t149;
t91 = t101 * t120 + t123 * t148;
t145 = t91 * MDP(21);
t92 = t101 * t123 - t120 * t148;
t144 = t92 * MDP(18);
t143 = t92 * MDP(20);
t142 = MDP(15) * t125;
t141 = MDP(18) * t123;
t100 = -t119 * t124 + t121 * t149;
t140 = t100 * MDP(22);
t139 = t101 * MDP(12);
t138 = t101 * MDP(13);
t137 = t120 * MDP(21);
t136 = t124 * MDP(22);
t135 = t123 * t158;
t133 = t120 * t123 * MDP(19);
t110 = pkin(7) * t149;
t97 = t110 + (-pkin(2) - t161) * t119;
t85 = t100 * pkin(3) - t101 * pkin(9) + t97;
t89 = t121 * t99 + t124 * t98;
t87 = -pkin(9) * t148 + t89;
t81 = -t120 * t87 + t123 * t85;
t132 = pkin(8) * MDP(16) - MDP(13);
t131 = -(MDP(25) * pkin(4)) + MDP(20);
t82 = t120 * t85 + t123 * t87;
t106 = -t124 * pkin(3) - t121 * pkin(9) - pkin(2);
t104 = t123 * t106;
t95 = -t120 * t158 + t104;
t96 = t120 * t106 + t135;
t130 = t95 * MDP(23) - t96 * MDP(24);
t128 = t123 * MDP(20) - MDP(12) - t137;
t126 = t81 * MDP(23) - t82 * MDP(24) + t140 + t143 - t145;
t117 = t123 ^ 2;
t116 = t121 ^ 2;
t115 = t120 ^ 2;
t114 = t118 ^ 2;
t113 = -t123 * pkin(4) - pkin(3);
t108 = t157 * t123;
t107 = t157 * t120;
t105 = (pkin(4) * t120 + pkin(8)) * t121;
t103 = t119 * t162 + t134;
t102 = t119 * t161 - t110;
t93 = t135 + (t106 - t150) * t120;
t90 = -t123 * t150 + t104 + (-pkin(4) - t160) * t124;
t83 = t91 * pkin(4) + t86;
t80 = -t91 * qJ(5) + t82;
t79 = t100 * pkin(4) - t92 * qJ(5) + t81;
t1 = [t114 * t122 ^ 2 * MDP(4) + t101 ^ 2 * MDP(11) + (t79 ^ 2 + t80 ^ 2 + t83 ^ 2) * MDP(26) + MDP(1) + (-0.2e1 * t91 * MDP(19) + t144) * t92 + (t146 * t167 + t147) * t119 + ((MDP(7) * t119 - t138) * t167 + (0.2e1 * MDP(5) * t122 + t142) * t114) * t125 + (0.2e1 * MDP(14) * t148 - 0.2e1 * t139 + t140 + 0.2e1 * t143 - 0.2e1 * t145) * t100 + 0.2e1 * (t102 * t119 + t114 * t161) * MDP(9) + 0.2e1 * (-t103 * t119 - t114 * t162) * MDP(10) + 0.2e1 * (t97 * t100 - t88 * t148) * MDP(16) + 0.2e1 * (t97 * t101 + t89 * t148) * MDP(17) + (t81 * t100 + t86 * t91) * t165 + (-t82 * t100 + t86 * t92) * t164 + (-t79 * t92 - t80 * t91) * t163; t147 + t102 * MDP(9) - t103 * MDP(10) - t101 * t154 + (-t90 * t92 - t93 * t91) * MDP(25) + (t83 * t105 + t79 * t90 + t80 * t93) * MDP(26) + (t125 * MDP(7) + t146) * t118 + (t130 - t155) * t100 + (t139 - t97 * MDP(16) + (-MDP(14) + t153) * t148 - t126) * t124 + (t101 * MDP(11) + t97 * MDP(17) + t92 * t141 + (-t120 * t92 - t123 * t91) * MDP(19) + (pkin(8) * t91 + t152) * MDP(23) + (pkin(8) * t92 + t151) * MDP(24) + (-t120 * t80 - t123 * t79) * MDP(25) + t132 * t148 + t128 * t100) * t121; MDP(8) + (t105 ^ 2 + t90 ^ 2 + t93 ^ 2) * MDP(26) + (t136 + 0.2e1 * t155) * t124 + (t117 * MDP(18) + MDP(11) - 0.2e1 * t133) * t116 + (t116 * t160 - t95 * t124) * t165 + (t116 * t159 + t96 * t124) * t164 + (-0.2e1 * t154 + (-t120 * t93 - t123 * t90) * t163 - 0.2e1 * t128 * t124) * t121; t138 - t118 * t142 + t88 * MDP(16) - t89 * MDP(17) + t120 * t144 + (-t120 * t91 + t92 * t123) * MDP(19) + (-pkin(3) * t91 - t151) * MDP(23) + (-pkin(3) * t92 + t152) * MDP(24) + (-t107 * t92 + t108 * t91 - t79 * t120 + t80 * t123) * MDP(25) + (t79 * t107 - t80 * t108 + t83 * t113) * MDP(26) + t166 * t100; (-t90 * t120 + t93 * t123) * MDP(25) + (t105 * t113 + t90 * t107 - t93 * t108) * MDP(26) + (-t153 - t166) * t124 + (t120 * t141 + (-t115 + t117) * MDP(19) + (-pkin(3) * t120 - t159) * MDP(23) + (-pkin(3) * t123 + t160) * MDP(24) + (-t107 * t123 + t108 * t120) * MDP(25) - t132) * t121; MDP(15) + t115 * MDP(18) + 0.2e1 * t133 + (-t107 * t120 - t108 * t123) * t163 + (t107 ^ 2 + t108 ^ 2 + t113 ^ 2) * MDP(26) + 0.2e1 * (t123 * MDP(23) - t120 * MDP(24)) * pkin(3); (-t92 * MDP(25) + t79 * MDP(26)) * pkin(4) + t126; t90 * t156 - t136 + (t131 * t123 - t137) * t121 + t130; t107 * t156 + (-MDP(24) * pkin(9) + MDP(21)) * t123 + (-MDP(23) * pkin(9) + t131) * t120; MDP(26) * (pkin(4) ^ 2) + MDP(22); t83 * MDP(26); t105 * MDP(26); t113 * MDP(26); 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
