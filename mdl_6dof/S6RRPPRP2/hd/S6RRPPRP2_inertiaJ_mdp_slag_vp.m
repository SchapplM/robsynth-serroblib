% Calculate joint inertia matrix for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RRPPRP2_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:31:58
% EndTime: 2019-03-09 08:31:59
% DurationCPUTime: 0.37s
% Computational Cost: add. (599->129), mult. (1028->180), div. (0->0), fcn. (1062->6), ass. (0->61)
t102 = sin(pkin(9));
t103 = cos(pkin(9));
t105 = sin(qJ(2));
t107 = cos(qJ(2));
t86 = t102 * t105 - t103 * t107;
t137 = 0.2e1 * t86;
t104 = sin(qJ(5));
t136 = 0.2e1 * t104;
t106 = cos(qJ(5));
t119 = t106 * MDP(23);
t135 = t119 + MDP(15);
t87 = t102 * t107 + t103 * t105;
t99 = -t107 * pkin(2) - pkin(1);
t109 = -t87 * qJ(4) + t99;
t77 = t86 * pkin(3) + t109;
t134 = -0.2e1 * t77;
t132 = 0.2e1 * t107;
t95 = t102 * pkin(2) + qJ(4);
t131 = t95 * t86;
t130 = -qJ(3) - pkin(7);
t129 = MDP(25) * pkin(5);
t100 = t104 ^ 2;
t128 = t100 * t86;
t90 = t130 * t105;
t91 = t130 * t107;
t78 = -t102 * t91 - t103 * t90;
t74 = t87 * pkin(4) + t78;
t127 = t104 * t74;
t98 = -t103 * pkin(2) - pkin(3);
t94 = -pkin(8) + t98;
t126 = -qJ(6) + t94;
t125 = t86 * MDP(14);
t124 = t87 * MDP(21);
t101 = t106 ^ 2;
t93 = -t100 - t101;
t123 = -t93 * MDP(25) + MDP(16);
t122 = MDP(20) * t106;
t121 = t104 * MDP(23);
t120 = t106 * MDP(18);
t118 = MDP(11) + MDP(13);
t80 = t102 * t90 - t103 * t91;
t117 = t78 ^ 2 + t80 ^ 2;
t72 = (pkin(3) + pkin(8)) * t86 + t109;
t116 = qJ(6) * t86 + t72;
t115 = -MDP(24) * pkin(5) + MDP(19);
t114 = MDP(22) + t129;
t113 = -t87 * t94 + t131;
t73 = t106 * t74;
t67 = t87 * pkin(5) - t116 * t104 + t73;
t68 = t116 * t106 + t127;
t112 = -t67 * t104 + t68 * t106;
t111 = (-t104 * t72 + t73) * MDP(22) - (t106 * t72 + t127) * MDP(23);
t110 = -t106 * MDP(22) + t121;
t89 = t104 * pkin(5) + t95;
t84 = t126 * t106;
t83 = t126 * t104;
t82 = t101 * t86;
t76 = t83 * t104 + t84 * t106;
t75 = -t86 * pkin(4) + t80;
t71 = (-pkin(5) * t106 - pkin(4)) * t86 + t80;
t1 = [MDP(1) + pkin(1) * MDP(9) * t132 + (t99 ^ 2 + t117) * MDP(12) + t125 * t134 + (t77 ^ 2 + t117) * MDP(16) + (t67 ^ 2 + t68 ^ 2 + t71 ^ 2) * MDP(25) + (t100 * MDP(17) + t120 * t136) * t86 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t105 + MDP(5) * t132) * t105 + (MDP(15) * t134 + t124 + (MDP(19) * t104 + t122) * t137) * t87 + 0.2e1 * (t118 * t78 + t111) * t87 + (t112 * MDP(24) + t110 * t75 - t118 * t80) * t137; t105 * MDP(6) + t107 * MDP(7) + (t98 * t87 - t131) * MDP(13) + t78 * MDP(14) + t80 * MDP(15) + (t78 * t98 + t80 * t95) * MDP(16) + (t82 - t128) * MDP(18) + (t67 * t84 + t68 * t83 + t71 * t89) * MDP(25) + (-t107 * MDP(10) - t105 * MDP(9)) * pkin(7) + (t87 * MDP(19) - t113 * MDP(22) + t75 * MDP(23) + (t83 * t86 - t67) * MDP(24)) * t106 + (t86 * t106 * MDP(17) - t87 * MDP(20) + t75 * MDP(22) + t113 * MDP(23) + (-t84 * t86 - t68) * MDP(24)) * t104 + ((-t102 * t86 - t103 * t87) * MDP(11) + (t102 * t80 - t103 * t78) * MDP(12)) * pkin(2); MDP(8) + (t95 ^ 2 + t98 ^ 2) * MDP(16) + t101 * MDP(17) + (t83 ^ 2 + t84 ^ 2 + t89 ^ 2) * MDP(25) + (t102 ^ 2 + t103 ^ 2) * MDP(12) * pkin(2) ^ 2 + (t95 * MDP(22) - t120) * t136 + 0.2e1 * t98 * MDP(14) - 0.2e1 * t76 * MDP(24) + 0.2e1 * t135 * t95; t99 * MDP(12) - t125 + t77 * MDP(16) + (t82 + t128) * MDP(24) + t112 * MDP(25) + (-t104 * MDP(22) - t135) * t87; (-t84 * t104 + t83 * t106) * MDP(25); MDP(12) + t123; t78 * MDP(16) + (t68 * t104 + t67 * t106) * MDP(25) + (MDP(13) - t110) * t87; t98 * MDP(16) + t93 * MDP(24) + t76 * MDP(25) + MDP(14); 0; t123; t67 * t129 + t124 + (t115 * t104 + t122) * t86 + t111; t84 * t129 + (-t94 * MDP(23) - MDP(20)) * t104 + (t94 * MDP(22) + t115) * t106; -t114 * t104 - t119; t114 * t106 - t121; MDP(25) * pkin(5) ^ 2 + MDP(21); t71 * MDP(25); t89 * MDP(25); 0; 0; 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
