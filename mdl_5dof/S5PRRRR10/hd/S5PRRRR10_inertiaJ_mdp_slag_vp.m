% Calculate joint inertia matrix for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR10_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:08
% EndTime: 2019-12-05 17:26:10
% DurationCPUTime: 0.54s
% Computational Cost: add. (482->155), mult. (1242->241), div. (0->0), fcn. (1360->12), ass. (0->81)
t114 = sin(pkin(6));
t164 = 0.2e1 * t114;
t118 = sin(qJ(5));
t122 = cos(qJ(5));
t127 = -(MDP(24) * t118 + MDP(25) * t122) * pkin(10) + t118 * MDP(21) + t122 * MDP(22);
t163 = 0.2e1 * MDP(24);
t162 = 0.2e1 * MDP(25);
t120 = sin(qJ(3));
t161 = pkin(2) * t120;
t124 = cos(qJ(3));
t160 = pkin(2) * t124;
t159 = pkin(9) * t118;
t158 = pkin(9) * t122;
t123 = cos(qJ(4));
t157 = pkin(9) * t123;
t156 = pkin(3) * MDP(17);
t155 = pkin(3) * MDP(18);
t119 = sin(qJ(4));
t115 = sin(pkin(5));
t116 = cos(pkin(6));
t117 = cos(pkin(5));
t125 = cos(qJ(2));
t100 = -t115 * t125 * t114 + t117 * t116;
t121 = sin(qJ(2));
t148 = t116 * t125;
t151 = t114 * t120;
t91 = t117 * t151 + (t120 * t148 + t121 * t124) * t115;
t83 = -t100 * t123 + t91 * t119;
t154 = t119 * t83;
t150 = t114 * t124;
t137 = pkin(8) * t150;
t98 = t137 + (pkin(9) + t161) * t116;
t99 = (-pkin(3) * t124 - pkin(9) * t120 - pkin(2)) * t114;
t88 = -t119 * t98 + t123 * t99;
t86 = pkin(4) * t150 - t88;
t153 = t86 * t118;
t152 = t86 * t122;
t149 = t116 * MDP(9);
t147 = t120 * MDP(7);
t102 = t119 * t116 + t123 * t151;
t92 = t118 * t102 + t122 * t150;
t146 = t92 * MDP(22);
t93 = t122 * t102 - t118 * t150;
t145 = t93 * MDP(19);
t144 = t93 * MDP(21);
t143 = MDP(16) * t124;
t101 = -t123 * t116 + t119 * t151;
t142 = t101 * MDP(23);
t141 = t102 * MDP(13);
t140 = t102 * MDP(14);
t139 = t122 * MDP(19);
t138 = t123 * MDP(23);
t136 = t118 * t122 * MDP(20);
t135 = pkin(9) * MDP(17) - MDP(14);
t134 = pkin(9) * MDP(18) - MDP(15);
t89 = t119 * t99 + t123 * t98;
t106 = -t123 * pkin(4) - t119 * pkin(10) - pkin(3);
t94 = t122 * t106 - t118 * t157;
t95 = t118 * t106 + t122 * t157;
t133 = t94 * MDP(24) - t95 * MDP(25);
t132 = MDP(21) * t122 - MDP(22) * t118;
t130 = t122 * MDP(24) - t118 * MDP(25);
t108 = pkin(8) * t151;
t97 = t108 + (-pkin(3) - t160) * t116;
t128 = -MDP(13) + t132;
t85 = t101 * pkin(4) - t102 * pkin(10) + t97;
t87 = -pkin(10) * t150 + t89;
t79 = -t118 * t87 + t122 * t85;
t80 = t118 * t85 + t122 * t87;
t126 = t79 * MDP(24) - t80 * MDP(25) + t142 + t144 - t146;
t113 = t122 ^ 2;
t112 = t119 ^ 2;
t111 = t118 ^ 2;
t110 = t114 ^ 2;
t104 = t116 * t161 + t137;
t103 = t116 * t160 - t108;
t90 = -t117 * t150 + (t120 * t121 - t124 * t148) * t115;
t84 = t100 * t119 + t91 * t123;
t82 = t90 * t118 + t84 * t122;
t81 = -t84 * t118 + t90 * t122;
t1 = [MDP(1); (-t100 * t150 - t90 * t116) * MDP(10) + (t100 * t151 - t91 * t116) * MDP(11) + (t90 * t101 + t83 * t150) * MDP(17) + (t90 * t102 + t84 * t150) * MDP(18) + (t81 * t101 + t83 * t92) * MDP(24) + (-t82 * t101 + t83 * t93) * MDP(25) + (t125 * MDP(3) - t121 * MDP(4)) * t115; t110 * t120 ^ 2 * MDP(5) + t102 ^ 2 * MDP(12) + MDP(2) + (-0.2e1 * t92 * MDP(20) + t145) * t93 + (t147 * t164 + t149) * t116 + ((MDP(8) * t116 - t140) * t164 + (0.2e1 * MDP(6) * t120 + t143) * t110) * t124 + (0.2e1 * MDP(15) * t150 - 0.2e1 * t141 + t142 + 0.2e1 * t144 - 0.2e1 * t146) * t101 + 0.2e1 * (t103 * t116 + t110 * t160) * MDP(10) + 0.2e1 * (-t104 * t116 - t110 * t161) * MDP(11) + 0.2e1 * (t97 * t101 - t88 * t150) * MDP(17) + 0.2e1 * (t97 * t102 + t89 * t150) * MDP(18) + (t79 * t101 + t86 * t92) * t163 + (-t80 * t101 + t86 * t93) * t162; -t91 * MDP(11) + (t118 * t154 - t81 * t123) * MDP(24) + (t122 * t154 + t82 * t123) * MDP(25) + (-t123 * MDP(17) + t119 * MDP(18) - MDP(10)) * t90; -t102 * t155 + t103 * MDP(10) - t104 * MDP(11) + t149 + (t124 * MDP(8) + t147) * t114 + (t133 - t156) * t101 + (-t97 * MDP(17) + t134 * t150 - t126 + t141) * t123 + (t102 * MDP(12) + t97 * MDP(18) + t93 * t139 + (-t118 * t93 - t122 * t92) * MDP(20) + (pkin(9) * t92 + t153) * MDP(24) + (pkin(9) * t93 + t152) * MDP(25) + t135 * t150 + t128 * t101) * t119; MDP(9) + (t138 + 0.2e1 * t156) * t123 + (t113 * MDP(19) + MDP(12) - 0.2e1 * t136) * t112 + (t112 * t159 - t94 * t123) * t163 + (t112 * t158 + t95 * t123) * t162 + 0.2e1 * (-t123 * t128 - t155) * t119; -t84 * MDP(18) + (-MDP(17) - t130) * t83; t140 - t114 * t143 + t88 * MDP(17) - t89 * MDP(18) + t118 * t145 + (-t118 * t92 + t93 * t122) * MDP(20) + (-pkin(4) * t92 - t152) * MDP(24) + (-pkin(4) * t93 + t153) * MDP(25) + (-MDP(15) + t127) * t101; (-t127 - t134) * t123 + (t118 * t139 + (-t111 + t113) * MDP(20) + (-pkin(4) * t118 - t158) * MDP(24) + (-pkin(4) * t122 + t159) * MDP(25) - t135) * t119; t111 * MDP(19) + 0.2e1 * pkin(4) * t130 + MDP(16) + 0.2e1 * t136; t81 * MDP(24) - t82 * MDP(25); t126; t132 * t119 + t133 - t138; t127; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
