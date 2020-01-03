% Calculate joint inertia matrix for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR10_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:42
% EndTime: 2019-12-31 19:44:44
% DurationCPUTime: 0.41s
% Computational Cost: add. (331->119), mult. (649->161), div. (0->0), fcn. (572->6), ass. (0->59)
t140 = pkin(6) * MDP(14);
t115 = MDP(12) - MDP(17);
t116 = MDP(11) + MDP(15);
t96 = sin(pkin(8));
t97 = cos(pkin(8));
t139 = t115 * t97 + t116 * t96;
t100 = cos(qJ(5));
t98 = sin(qJ(5));
t80 = t100 * t96 - t98 * t97;
t138 = MDP(14) + MDP(18);
t114 = t96 * qJ(4) + pkin(2);
t134 = pkin(3) + pkin(4);
t77 = t134 * t97 + t114;
t137 = 0.2e1 * t77;
t99 = sin(qJ(2));
t136 = 0.2e1 * t99;
t135 = -2 * MDP(20);
t133 = pkin(7) * t99;
t101 = cos(qJ(2));
t132 = pkin(6) * t101;
t130 = -pkin(7) + qJ(3);
t83 = -t101 * pkin(2) - t99 * qJ(3) - pkin(1);
t73 = t97 * t132 + t96 * t83;
t128 = pkin(2) * MDP(14);
t75 = t80 * t99;
t126 = t75 * MDP(22);
t79 = t100 * t97 + t98 * t96;
t76 = t79 * t99;
t125 = t76 * MDP(19);
t124 = t76 * MDP(21);
t123 = t79 * MDP(24);
t81 = -t97 * pkin(3) - t114;
t122 = t81 * MDP(18);
t121 = t96 * MDP(12);
t120 = t96 * MDP(17);
t119 = t97 * MDP(11);
t118 = t97 * MDP(15);
t117 = t101 * MDP(23);
t113 = qJ(4) * t97 - pkin(6);
t88 = t96 * t132;
t72 = t97 * t83 - t88;
t70 = -t101 * qJ(4) + t73;
t92 = t101 * pkin(3);
t71 = -t72 + t92;
t110 = t70 * t97 + t71 * t96;
t109 = -t72 * t96 + t73 * t97;
t108 = t96 * MDP(11) + t97 * MDP(12);
t107 = t96 * MDP(15) - t97 * MDP(17);
t65 = t101 * pkin(4) + t88 + t92 + (-t83 - t133) * t97;
t66 = t96 * t133 + t70;
t106 = (t100 * t65 - t98 * t66) * MDP(24) - (t100 * t66 + t98 * t65) * MDP(25);
t105 = -t75 * MDP(24) + t76 * MDP(25);
t104 = t100 * MDP(24) - t98 * MDP(25);
t84 = t130 * t96;
t85 = t130 * t97;
t103 = t80 * MDP(21) - t79 * MDP(22) + (t100 * t84 - t98 * t85) * MDP(24) - (t100 * t85 + t98 * t84) * MDP(25);
t74 = (pkin(3) * t96 - t113) * t99;
t69 = (-t134 * t96 + t113) * t99;
t1 = [MDP(1) + (t72 ^ 2 + t73 ^ 2) * MDP(14) + (t70 ^ 2 + t71 ^ 2 + t74 ^ 2) * MDP(18) + (-t75 * t135 + t125) * t76 + (MDP(5) * t136 + 0.2e1 * pkin(1) * MDP(9) + t117 + 0.2e1 * t124 + 0.2e1 * t126) * t101 + 0.2e1 * t105 * t69 + 0.2e1 * (-t72 * MDP(11) + t73 * MDP(12) + t71 * MDP(15) - t70 * MDP(17) + t106) * t101 + ((-t72 * t97 - t73 * t96) * MDP(13) + (-t70 * t96 + t71 * t97) * MDP(16) + t107 * t74) * t136 + (-0.2e1 * pkin(1) * MDP(10) + (MDP(4) + (0.2e1 * t108 + t140) * pkin(6)) * t99) * t99; t109 * MDP(13) + t110 * MDP(16) + t80 * t125 + (t80 * t75 - t76 * t79) * MDP(20) + (t69 * t79 - t77 * t75) * MDP(24) + (t69 * t80 + t77 * t76) * MDP(25) + (-t118 - t120 + t122) * t74 + (MDP(6) + t107 * t81 - t108 * pkin(2) + (-MDP(9) - t119 + t121 - t128) * pkin(6)) * t99 + (-pkin(6) * MDP(10) + MDP(7) + t103) * t101 + (t109 * MDP(14) + t110 * MDP(18) + t139 * t101) * qJ(3); MDP(8) + t123 * t137 + (-0.2e1 * t118 - 0.2e1 * t120 + t122) * t81 + (0.2e1 * t119 - 0.2e1 * t121 + t128) * pkin(2) + (MDP(19) * t80 + MDP(25) * t137 + t79 * t135) * t80 + (t138 * qJ(3) + 0.2e1 * MDP(13) + 0.2e1 * MDP(16)) * (t96 ^ 2 + t97 ^ 2) * qJ(3); t74 * MDP(18) + (t139 + t140) * t99 - t105; -t80 * MDP(25) + t115 * t96 - t116 * t97 + t122 - t123 - t128; t138; t97 * t99 * MDP(16) + t71 * MDP(18) + (MDP(15) + t104) * t101; (MDP(18) * qJ(3) + MDP(16)) * t96; 0; MDP(18); t106 + t117 + t124 + t126; t103; 0; t104; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
