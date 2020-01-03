% Calculate joint inertia matrix for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR7_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:22:29
% EndTime: 2019-12-31 22:22:30
% DurationCPUTime: 0.39s
% Computational Cost: add. (629->120), mult. (1173->162), div. (0->0), fcn. (1292->8), ass. (0->68)
t117 = sin(qJ(5));
t121 = cos(qJ(5));
t143 = t117 * MDP(27) + t121 * MDP(28);
t142 = MDP(30) * t121;
t126 = -0.2e1 * MDP(31) * t117 + 0.2e1 * t142;
t124 = cos(qJ(2));
t110 = -t124 * pkin(2) - pkin(1);
t119 = sin(qJ(3));
t120 = sin(qJ(2));
t123 = cos(qJ(3));
t99 = t119 * t120 - t123 * t124;
t91 = t99 * pkin(3) + t110;
t155 = 0.2e1 * t91;
t154 = 0.2e1 * t110;
t153 = 0.2e1 * t124;
t152 = pkin(6) + pkin(7);
t151 = pkin(2) * t119;
t150 = pkin(4) * t117;
t122 = cos(qJ(4));
t149 = t122 * pkin(3);
t118 = sin(qJ(4));
t100 = t119 * t124 + t123 * t120;
t102 = t152 * t120;
t103 = t152 * t124;
t137 = -t123 * t102 - t119 * t103;
t79 = -t100 * pkin(8) + t137;
t132 = t119 * t102 - t123 * t103;
t80 = -t99 * pkin(8) - t132;
t75 = t118 * t80 - t122 * t79;
t71 = t75 * t117;
t148 = t75 * t121;
t147 = t117 * t121;
t86 = t118 * t100 + t122 * t99;
t146 = t86 * MDP(29);
t109 = t123 * pkin(2) + pkin(3);
t104 = t122 * t109;
t94 = -t118 * t151 + t104;
t145 = t94 * MDP(23);
t95 = -t118 * t109 - t122 * t151;
t144 = t95 * MDP(24);
t141 = t118 * MDP(24);
t140 = t123 * MDP(16);
t139 = MDP(26) * t147;
t115 = t117 ^ 2;
t138 = t115 * MDP(25) + MDP(22) + 0.2e1 * t139;
t136 = MDP(15) + t138;
t87 = t122 * t100 - t118 * t99;
t135 = -pkin(4) * t87 - pkin(9) * t86;
t92 = -pkin(4) - t94;
t93 = pkin(9) - t95;
t134 = -t86 * t93 + t87 * t92;
t107 = t118 * pkin(3) + pkin(9);
t108 = -pkin(4) - t149;
t133 = -t107 * t86 + t108 * t87;
t131 = MDP(27) * t121 - MDP(28) * t117;
t129 = -MDP(30) * t117 - MDP(31) * t121;
t128 = (t122 * MDP(23) - t141) * pkin(3);
t116 = t121 ^ 2;
t76 = t118 * t79 + t122 * t80;
t127 = -t75 * MDP(23) - t76 * MDP(24) + ((-t115 + t116) * MDP(26) + MDP(25) * t147 + MDP(20)) * t87 + (-MDP(21) + t143) * t86;
t125 = t100 * MDP(13) - t99 * MDP(14) + t137 * MDP(16) + t132 * MDP(17) + t127;
t114 = pkin(4) * t121;
t105 = t108 * t117;
t90 = t92 * t117;
t74 = t86 * pkin(4) - t87 * pkin(9) + t91;
t70 = t117 * t74 + t121 * t76;
t69 = -t117 * t76 + t121 * t74;
t1 = [t99 * MDP(16) * t154 + t87 * MDP(24) * t155 + pkin(1) * MDP(9) * t153 + MDP(1) + (MDP(23) * t155 + t146 + 0.2e1 * (-MDP(19) + t131) * t87) * t86 + 0.2e1 * (t69 * t86 + t87 * t71) * MDP(30) + 0.2e1 * (t87 * t148 - t70 * t86) * MDP(31) + (t116 * MDP(25) + MDP(18) - 0.2e1 * t139) * t87 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t120 + MDP(5) * t153) * t120 + (MDP(11) * t100 - 0.2e1 * t99 * MDP(12) + MDP(17) * t154) * t100; t125 + t120 * MDP(6) + t124 * MDP(7) + (t134 * t117 - t148) * MDP(30) + (t134 * t121 + t71) * MDP(31) + (-t124 * MDP(10) - t120 * MDP(9)) * pkin(6); MDP(8) - t92 * t126 + 0.2e1 * (-t119 * MDP(17) + t140) * pkin(2) + 0.2e1 * t145 + 0.2e1 * t144 + t136; t125 + (t133 * t117 - t148) * MDP(30) + (t133 * t121 + t71) * MDP(31); (t104 + t149) * MDP(23) + (t105 + t90) * MDP(31) + (-t108 - t92) * t142 + (-pkin(3) - t109) * t141 + (t140 + (-MDP(23) * t118 - MDP(24) * t122 - MDP(17)) * t119) * pkin(2) + t136; -t108 * t126 + 0.2e1 * t128 + t136; (t135 * t117 - t148) * MDP(30) + (t135 * t121 + t71) * MDP(31) + t127; t145 + t144 + (-t92 * t121 + t114) * MDP(30) + (t90 - t150) * MDP(31) + t138; (-t108 * t121 + t114) * MDP(30) + (t105 - t150) * MDP(31) + t128 + t138; pkin(4) * t126 + t138; t69 * MDP(30) - t70 * MDP(31) + t131 * t87 + t146; t129 * t93 + t143; t129 * t107 + t143; t129 * pkin(9) + t143; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
