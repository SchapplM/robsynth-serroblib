% Calculate joint inertia matrix for
% S5RRRRP6
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
%   see S5RRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:11
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP6_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:10:26
% EndTime: 2021-01-16 00:10:29
% DurationCPUTime: 0.56s
% Computational Cost: add. (634->153), mult. (1143->202), div. (0->0), fcn. (1139->6), ass. (0->69)
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t133 = t114 * MDP(20) + t117 * MDP(21);
t115 = sin(qJ(3));
t116 = sin(qJ(2));
t118 = cos(qJ(2));
t149 = cos(qJ(3));
t95 = t115 * t116 - t118 * t149;
t157 = 0.2e1 * t95;
t151 = pkin(7) + pkin(6);
t100 = t151 * t116;
t101 = t151 * t118;
t84 = -t115 * t100 + t101 * t149;
t143 = t117 * t84;
t96 = t115 * t118 + t116 * t149;
t144 = qJ(5) * t96;
t107 = -t118 * pkin(2) - pkin(1);
t79 = t95 * pkin(3) - t96 * pkin(8) + t107;
t72 = t143 + (t79 - t144) * t114;
t156 = t72 * MDP(27);
t126 = -t117 * MDP(25) + t114 * MDP(26);
t155 = t114 * MDP(24);
t154 = 0.2e1 * t107;
t153 = 0.2e1 * t118;
t152 = 0.2e1 * MDP(27);
t150 = t95 * pkin(4);
t148 = t117 * pkin(4);
t105 = -pkin(2) * t149 - pkin(3);
t147 = pkin(3) - t105;
t146 = -qJ(5) - pkin(8);
t145 = MDP(28) * pkin(4);
t106 = -pkin(3) - t148;
t97 = t105 - t148;
t142 = t106 + t97;
t141 = MDP(24) * t95;
t140 = MDP(25) * t96;
t139 = t114 * t117;
t104 = t115 * pkin(2) + pkin(8);
t134 = qJ(5) + t104;
t91 = t134 * t114;
t138 = t91 * MDP(25);
t137 = t95 * MDP(22);
t136 = t97 * MDP(28);
t98 = t146 * t114;
t135 = t98 * MDP(25);
t132 = t106 * MDP(28);
t131 = t114 * MDP(21);
t130 = t114 * MDP(27);
t128 = MDP(19) * t139;
t112 = t114 ^ 2;
t127 = t112 * MDP(18) + MDP(15) + 0.2e1 * t128;
t73 = -t114 * t84 + t117 * t79;
t83 = t100 * t149 + t115 * t101;
t76 = t114 * t96 * pkin(4) + t83;
t125 = -t83 * MDP(23) - t76 * MDP(25);
t124 = t117 * MDP(23) - t155;
t123 = -t114 * MDP(23) - t117 * MDP(24);
t122 = 0.2e1 * t126;
t121 = t73 * MDP(23) - (t114 * t79 + t143) * MDP(24) - t72 * MDP(26);
t120 = (MDP(16) * t149 - t115 * MDP(17)) * pkin(2);
t113 = t117 ^ 2;
t119 = -t84 * MDP(17) + t117 * t156 + (t155 - MDP(16)) * t83 + ((-t112 + t113) * MDP(19) + MDP(18) * t139 + MDP(13)) * t96 + (-MDP(14) + t133) * t95;
t99 = t146 * t117;
t94 = t99 * t117;
t92 = t134 * t117;
t86 = t92 * t117;
t75 = t76 * t114;
t70 = -t117 * t144 + t150 + t73;
t1 = [MDP(1) + pkin(1) * MDP(9) * t153 + (t70 ^ 2 + t72 ^ 2 + t76 ^ 2) * MDP(28) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t116 + MDP(5) * t153) * t116 + (MDP(16) * t154 + t137) * t95 + (t70 * MDP(25) + t121) * t157 + (MDP(17) * t154 + (t117 * MDP(20) - MDP(12) - t131) * t157 + 0.2e1 * (t83 * MDP(24) + t76 * MDP(26) - t70 * MDP(27)) * t117 + 0.2e1 * (-t125 - t156) * t114 + (t113 * MDP(18) + MDP(11) - 0.2e1 * t128) * t96) * t96; (-t104 * t141 + (t105 * MDP(24) + t97 * MDP(26) + t91 * MDP(27)) * t96 + t125) * t117 + (-t92 * t95 + t75) * MDP(26) + ((-t104 * t95 + t105 * t96) * MDP(23) + t97 * t140 + (-t92 * t96 - t70) * MDP(27)) * t114 + (-t118 * MDP(10) - t116 * MDP(9)) * pkin(6) + t119 + t116 * MDP(6) + t118 * MDP(7) + (-t70 * t91 + t72 * t92 + t76 * t97) * MDP(28) - t95 * t138; MDP(8) + (t91 * t114 + t86) * t152 + (t91 ^ 2 + t92 ^ 2) * MDP(28) + (t122 + t136) * t97 + t127 - 0.2e1 * t105 * t124 + 0.2e1 * t120; t119 + (t76 * t106 + t70 * t98 - t72 * t99) * MDP(28) + (-pkin(8) * t141 + (-pkin(3) * MDP(24) + t106 * MDP(26) - t98 * MDP(27)) * t96 + t125) * t117 + t95 * t135 + ((-pkin(3) * t96 - pkin(8) * t95) * MDP(23) + t106 * t140 + (t96 * t99 - t70) * MDP(27)) * t114 + (t99 * t95 + t75) * MDP(26); (t86 - t94) * MDP(27) + (t97 * t106 - t91 * t98 - t92 * t99) * MDP(28) + t120 + (MDP(23) * t147 - MDP(25) * t142) * t117 + (-t147 * MDP(24) + t142 * MDP(26) + (t91 - t98) * MDP(27)) * t114 + t127; (-t98 * t114 - t94) * t152 + (t98 ^ 2 + t99 ^ 2) * MDP(28) + (t122 + t132) * t106 + 0.2e1 * t124 * pkin(3) + t127; t137 + (t73 + 0.2e1 * t150) * MDP(25) + t70 * t145 + (-t131 + (-MDP(25) * qJ(5) - MDP(27) * pkin(4) + MDP(20)) * t117) * t96 + t121; -t138 - t92 * MDP(26) + t123 * t104 + (-t91 * MDP(28) - t130) * pkin(4) + t133; t135 + t99 * MDP(26) + t123 * pkin(8) + (MDP(28) * t98 - t130) * pkin(4) + t133; MDP(22) + (0.2e1 * MDP(25) + t145) * pkin(4); t76 * MDP(28) + (t114 * MDP(25) + t117 * MDP(26)) * t96; t126 + t136; t126 + t132; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
