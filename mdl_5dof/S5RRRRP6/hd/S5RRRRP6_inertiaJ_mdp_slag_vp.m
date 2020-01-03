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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP6_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:35
% EndTime: 2019-12-31 21:54:36
% DurationCPUTime: 0.35s
% Computational Cost: add. (475->120), mult. (868->166), div. (0->0), fcn. (859->6), ass. (0->57)
t105 = sin(qJ(4));
t108 = cos(qJ(4));
t123 = t105 * MDP(20) + t108 * MDP(21);
t132 = 2 * MDP(25);
t135 = t105 * MDP(24);
t110 = cos(qJ(2));
t98 = -t110 * pkin(2) - pkin(1);
t134 = 0.2e1 * t98;
t133 = 0.2e1 * t110;
t131 = -pkin(7) - pkin(6);
t109 = cos(qJ(3));
t129 = t109 * pkin(2);
t96 = -pkin(3) - t129;
t130 = pkin(3) - t96;
t106 = sin(qJ(3));
t107 = sin(qJ(2));
t91 = t131 * t107;
t92 = t131 * t110;
t75 = t106 * t91 - t109 * t92;
t128 = t108 * t75;
t87 = t106 * t110 + t109 * t107;
t127 = MDP(25) * t87;
t126 = t105 * t108;
t102 = t108 * qJ(5);
t74 = -t106 * t92 - t109 * t91;
t125 = t74 * MDP(23);
t86 = t106 * t107 - t109 * t110;
t124 = t86 * MDP(22);
t122 = MDP(21) * t105;
t121 = t105 * MDP(25);
t120 = t108 * MDP(23);
t103 = t105 ^ 2;
t118 = MDP(19) * t126;
t119 = t103 * MDP(18) + MDP(15) + 0.2e1 * t118;
t97 = -t108 * pkin(4) - pkin(3);
t70 = t86 * pkin(3) - t87 * pkin(8) + t98;
t65 = -t105 * t75 + t108 * t70;
t117 = -pkin(3) * t87 - pkin(8) * t86;
t95 = t106 * pkin(2) + pkin(8);
t116 = -t86 * t95 + t87 * t96;
t115 = t65 * MDP(23) - (t105 * t70 + t128) * MDP(24);
t114 = t120 - t135;
t113 = t105 * MDP(23) + t108 * MDP(24);
t112 = (t109 * MDP(16) - t106 * MDP(17)) * pkin(2);
t104 = t108 ^ 2;
t64 = t128 + (-qJ(5) * t87 + t70) * t105;
t111 = MDP(25) * t108 * t64 - t75 * MDP(17) + (t135 - MDP(16)) * t74 + ((-t103 + t104) * MDP(19) + MDP(18) * t126 + MDP(13)) * t87 + (-MDP(14) + t123) * t86;
t90 = t108 * pkin(8) + t102;
t89 = (-qJ(5) - pkin(8)) * t105;
t88 = t97 - t129;
t85 = t90 * t108;
t83 = t108 * t95 + t102;
t82 = (-qJ(5) - t95) * t105;
t77 = t83 * t108;
t67 = t105 * t87 * pkin(4) + t74;
t62 = t86 * pkin(4) - t87 * t102 + t65;
t1 = [MDP(1) + pkin(1) * MDP(9) * t133 + (t62 ^ 2 + t64 ^ 2 + t67 ^ 2) * MDP(26) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t107 + MDP(5) * t133) * t107 + (MDP(17) * t134 + (-t105 * t64 - t108 * t62) * t132 + 0.2e1 * t113 * t74 + (t104 * MDP(18) + MDP(11) - 0.2e1 * t118) * t87) * t87 + (MDP(16) * t134 + t124 + 0.2e1 * (MDP(20) * t108 - MDP(12) - t122) * t87 + 0.2e1 * t115) * t86; (-t110 * MDP(10) - t107 * MDP(9)) * pkin(6) + t107 * MDP(6) + t110 * MDP(7) + (t116 * MDP(23) + (-t83 * t87 - t62) * MDP(25)) * t105 + (t116 * MDP(24) - t82 * t127 - t125) * t108 + (t62 * t82 + t64 * t83 + t67 * t88) * MDP(26) + t111; MDP(8) + (-t82 * t105 + t77) * t132 + (t82 ^ 2 + t83 ^ 2 + t88 ^ 2) * MDP(26) + t119 - 0.2e1 * t114 * t96 + 0.2e1 * t112; (t62 * t89 + t64 * t90 + t67 * t97) * MDP(26) + (t117 * MDP(24) - t89 * t127 - t125) * t108 + (t117 * MDP(23) + (-t87 * t90 - t62) * MDP(25)) * t105 + t111; (t77 + t85) * MDP(25) + (t82 * t89 + t83 * t90 + t88 * t97) * MDP(26) + t130 * t120 + t112 + (-t130 * MDP(24) + (-t82 - t89) * MDP(25)) * t105 + t119; (-t89 * t105 + t85) * t132 + (t89 ^ 2 + t90 ^ 2 + t97 ^ 2) * MDP(26) + 0.2e1 * t114 * pkin(3) + t119; t62 * pkin(4) * MDP(26) + t124 + (-t122 + (-MDP(25) * pkin(4) + MDP(20)) * t108) * t87 + t115; -t113 * t95 + (t82 * MDP(26) - t121) * pkin(4) + t123; -t113 * pkin(8) + (MDP(26) * t89 - t121) * pkin(4) + t123; MDP(26) * pkin(4) ^ 2 + MDP(22); t67 * MDP(26); t88 * MDP(26); t97 * MDP(26); 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
