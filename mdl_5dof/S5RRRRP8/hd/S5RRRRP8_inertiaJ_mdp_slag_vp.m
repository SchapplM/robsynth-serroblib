% Calculate joint inertia matrix for
% S5RRRRP8
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
%   see S5RRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP8_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:10
% EndTime: 2019-12-31 22:02:11
% DurationCPUTime: 0.46s
% Computational Cost: add. (479->127), mult. (946->191), div. (0->0), fcn. (909->6), ass. (0->57)
t102 = sin(qJ(3));
t105 = cos(qJ(3));
t101 = sin(qJ(4));
t104 = cos(qJ(4));
t128 = pkin(7) + pkin(8);
t89 = t128 * t102;
t90 = t128 * t105;
t72 = -t101 * t90 - t104 * t89;
t73 = -t101 * t89 + t104 * t90;
t118 = t104 * t105;
t84 = t101 * t102 - t118;
t85 = t101 * t105 + t102 * t104;
t111 = t85 * MDP(20) - t84 * MDP(21) + t72 * MDP(23) - t73 * MDP(24);
t134 = t111 - (MDP(16) * t102 + MDP(17) * t105) * pkin(7) + t102 * MDP(13) + t105 * MDP(14);
t115 = t104 * MDP(23);
t133 = (-MDP(24) * t101 + t115) * pkin(3);
t103 = sin(qJ(2));
t79 = t85 * t103;
t120 = t102 * t103;
t80 = -t101 * t120 + t103 * t118;
t122 = t80 * MDP(20) - t79 * MDP(21);
t106 = cos(qJ(2));
t123 = pkin(8) * t103;
t126 = pkin(6) * t102;
t88 = -pkin(2) * t106 - t103 * pkin(7) - pkin(1);
t83 = t105 * t88;
t67 = -t105 * t123 + t83 + (-pkin(3) - t126) * t106;
t124 = pkin(6) * t106;
t113 = t105 * t124;
t69 = t113 + (t88 - t123) * t102;
t62 = -t101 * t69 + t104 * t67;
t63 = t101 * t67 + t104 * t69;
t132 = t62 * MDP(23) - t63 * MDP(24) + t122;
t130 = -2 * MDP(19);
t129 = 2 * MDP(25);
t127 = pkin(3) * t101;
t125 = pkin(6) * t105;
t87 = pkin(3) * t120 + t103 * pkin(6);
t121 = MDP(18) * t80;
t119 = t102 * t105;
t114 = MDP(15) + MDP(22);
t95 = -pkin(3) * t105 - pkin(2);
t112 = MDP(12) * t119;
t110 = MDP(13) * t105 - MDP(14) * t102;
t99 = t105 ^ 2;
t98 = t103 ^ 2;
t97 = t102 ^ 2;
t94 = pkin(3) * t104 + pkin(4);
t76 = pkin(4) * t84 + t95;
t75 = t102 * t88 + t113;
t74 = -t102 * t124 + t83;
t68 = pkin(4) * t79 + t87;
t65 = -qJ(5) * t84 + t73;
t64 = -qJ(5) * t85 + t72;
t61 = -qJ(5) * t79 + t63;
t60 = -pkin(4) * t106 - t80 * qJ(5) + t62;
t1 = [MDP(1) - 0.2e1 * pkin(1) * t103 * MDP(10) + (t60 ^ 2 + t61 ^ 2 + t68 ^ 2) * MDP(26) + (t79 * t130 + t121) * t80 + t114 * t106 ^ 2 + (MDP(11) * t99 + MDP(4) - 0.2e1 * t112) * t98 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t110) * t103 - t122) * t106 + 0.2e1 * (-t106 * t74 + t98 * t126) * MDP(16) + 0.2e1 * (t106 * t75 + t98 * t125) * MDP(17) + 0.2e1 * (-t106 * t62 + t87 * t79) * MDP(23) + 0.2e1 * (t106 * t63 + t87 * t80) * MDP(24) + (-t60 * t80 - t61 * t79) * t129; t85 * t121 + (-t79 * t85 - t80 * t84) * MDP(19) + (t95 * t79 + t87 * t84) * MDP(23) + (t95 * t80 + t87 * t85) * MDP(24) + (-t60 * t85 - t61 * t84 - t64 * t80 - t65 * t79) * MDP(25) + (t60 * t64 + t61 * t65 + t68 * t76) * MDP(26) + (-pkin(6) * MDP(10) + MDP(7) - t134) * t106 + (MDP(6) - pkin(6) * MDP(9) + MDP(11) * t119 + (-t97 + t99) * MDP(12) + (-pkin(2) * t102 - t125) * MDP(16) + (-pkin(2) * t105 + t126) * MDP(17)) * t103; MDP(8) + t97 * MDP(11) + 0.2e1 * t112 - t65 * t84 * t129 + (t64 ^ 2 + t65 ^ 2 + t76 ^ 2) * MDP(26) + (MDP(18) * t85 - t64 * t129 + t84 * t130) * t85 + 0.2e1 * (MDP(23) * t84 + MDP(24) * t85) * t95 + 0.2e1 * (MDP(16) * t105 - MDP(17) * t102) * pkin(2); t74 * MDP(16) - t75 * MDP(17) + (-t79 * t127 - t80 * t94) * MDP(25) + (t61 * t127 + t60 * t94) * MDP(26) + (-t114 - t133) * t106 + t110 * t103 + t132; (-t84 * t127 - t85 * t94) * MDP(25) + (t65 * t127 + t64 * t94) * MDP(26) + t134; t94 ^ 2 * MDP(26) + (0.2e1 * t115 + (MDP(26) * t127 - 0.2e1 * MDP(24)) * t101) * pkin(3) + t114; -t106 * MDP(22) + (-MDP(25) * t80 + MDP(26) * t60) * pkin(4) + t132; (-MDP(25) * t85 + MDP(26) * t64) * pkin(4) + t111; MDP(26) * pkin(4) * t94 + MDP(22) + t133; MDP(26) * pkin(4) ^ 2 + MDP(22); t68 * MDP(26); t76 * MDP(26); 0; 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
