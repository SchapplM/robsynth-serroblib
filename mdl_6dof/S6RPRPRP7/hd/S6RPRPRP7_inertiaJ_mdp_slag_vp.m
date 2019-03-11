% Calculate joint inertia matrix for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPRP7_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:22:47
% EndTime: 2019-03-09 03:22:49
% DurationCPUTime: 0.34s
% Computational Cost: add. (546->119), mult. (897->170), div. (0->0), fcn. (924->6), ass. (0->57)
t98 = -pkin(1) - pkin(7);
t129 = -qJ(4) + t98;
t128 = sin(qJ(3));
t93 = sin(pkin(9));
t94 = cos(pkin(9));
t97 = cos(qJ(3));
t78 = -t93 * t128 + t94 * t97;
t76 = t78 ^ 2;
t127 = (pkin(1) * MDP(6));
t95 = sin(qJ(5));
t96 = cos(qJ(5));
t126 = t95 * t96;
t112 = t129 * t97;
t82 = t129 * t128;
t72 = t93 * t112 + t94 * t82;
t125 = t96 * t72;
t91 = t95 ^ 2;
t92 = t96 ^ 2;
t124 = t91 + t92;
t123 = MDP(15) * pkin(3);
t122 = MDP(24) * pkin(5);
t121 = qJ(6) * t78;
t89 = t128 * pkin(3) + qJ(2);
t87 = t93 * pkin(3) + pkin(8);
t120 = qJ(6) + t87;
t119 = MDP(19) * t95;
t70 = -t94 * t112 + t93 * t82;
t67 = t95 * t78 * pkin(5) + t70;
t118 = t67 * MDP(24);
t88 = -t94 * pkin(3) - pkin(4);
t83 = -t96 * pkin(5) + t88;
t117 = t83 * MDP(24);
t116 = t95 * MDP(21);
t115 = t95 * MDP(22);
t114 = t96 * MDP(22);
t113 = MDP(17) * t126;
t79 = -t94 * t128 - t93 * t97;
t69 = -t79 * pkin(4) - t78 * pkin(8) + t89;
t65 = t96 * t69 - t95 * t72;
t111 = t128 * MDP(13);
t110 = t124 * MDP(23);
t109 = -MDP(23) * pkin(5) + MDP(18);
t108 = MDP(21) + t122;
t63 = -t79 * pkin(5) - t96 * t121 + t65;
t64 = t125 + (t69 - t121) * t95;
t107 = t63 * t96 + t64 * t95;
t106 = t63 * t95 - t64 * t96;
t73 = t120 * t95;
t74 = t120 * t96;
t105 = -t73 * t96 + t74 * t95;
t104 = -t73 * t95 - t74 * t96;
t103 = t65 * MDP(21) - (t95 * t69 + t125) * MDP(22);
t102 = -t96 * MDP(21) + t115;
t101 = t114 + t116;
t100 = MDP(14) + t101;
t77 = t79 ^ 2;
t1 = [MDP(1) + (t70 ^ 2 + t72 ^ 2 + t89 ^ 2) * MDP(15) + t77 * MDP(20) + (t63 ^ 2 + t64 ^ 2 + t67 ^ 2) * MDP(24) + (MDP(7) * t97 - 0.2e1 * t128 * MDP(8)) * t97 + (t92 * MDP(16) - 0.2e1 * t113) * t76 + ((-2 * MDP(4) + t127) * pkin(1)) + (0.2e1 * t128 * MDP(12) + 0.2e1 * t97 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + 0.2e1 * (-t107 * MDP(23) + t100 * t70) * t78 + 0.2e1 * ((-MDP(18) * t96 + t119) * t78 + t72 * MDP(14) - t103) * t79; -t77 * t116 - t127 + MDP(4) + (-t72 * MDP(15) + t106 * MDP(24) + (-MDP(14) - t114) * t79) * t79 + (-t70 * MDP(15) - t100 * t78 - t118) * t78; MDP(6) + (t77 + t76) * MDP(15) + (t124 * t77 + t76) * MDP(24); -t128 * MDP(10) - t98 * t111 - t106 * MDP(23) + (-t63 * t73 + t64 * t74 + t67 * t83) * MDP(24) + (t98 * MDP(12) + MDP(9)) * t97 + t102 * t70 + (-t95 * MDP(18) - t96 * MDP(19) + t101 * t87) * t79 + (MDP(16) * t126 + (-t91 + t92) * MDP(17) - t105 * MDP(23) + t101 * t88) * t78 + ((-t78 * t94 + t79 * t93) * MDP(14) + (-t70 * t94 + t72 * t93) * MDP(15)) * pkin(3); t97 * MDP(12) - t111 + (t94 * t123 - t102 - t117) * t78 + (t104 * MDP(24) - t93 * t123 - t110) * t79; MDP(11) + t91 * MDP(16) + 0.2e1 * t113 - 0.2e1 * t104 * MDP(23) + (t73 ^ 2 + t74 ^ 2 + t83 ^ 2) * MDP(24) + (t93 ^ 2 + t94 ^ 2) * MDP(15) * pkin(3) ^ 2 + 0.2e1 * t102 * t88; t89 * MDP(15) + t107 * MDP(24) + t102 * t79 - t78 * t110; 0; t105 * MDP(24); t124 * MDP(24) + MDP(15); t63 * t122 - t79 * MDP(20) + (t109 * t96 - t119) * t78 + t103; (t108 * t95 + t114) * t79; -t73 * t122 + (-MDP(22) * t87 + MDP(19)) * t96 + (-MDP(21) * t87 + t109) * t95; t108 * t96 - t115; MDP(24) * pkin(5) ^ 2 + MDP(20); t118; -t78 * MDP(24); t117; 0; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
