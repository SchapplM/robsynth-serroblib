% Calculate joint inertia matrix for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRPRP10_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:32:33
% EndTime: 2019-03-09 03:32:34
% DurationCPUTime: 0.39s
% Computational Cost: add. (415->143), mult. (626->181), div. (0->0), fcn. (468->4), ass. (0->62)
t110 = MDP(24) - MDP(27);
t111 = MDP(23) + MDP(25);
t91 = sin(qJ(5));
t93 = cos(qJ(5));
t97 = -t110 * t91 + t111 * t93;
t131 = -MDP(14) - t97;
t94 = cos(qJ(3));
t128 = 0.2e1 * t94;
t85 = t91 ^ 2;
t87 = t93 ^ 2;
t81 = t85 + t87;
t130 = t81 * MDP(28);
t106 = -MDP(28) * pkin(5) - MDP(25);
t129 = (MDP(28) * qJ(6) - t110) * t91 + (MDP(23) - t106) * t93;
t127 = -2 * MDP(15);
t126 = 0.2e1 * MDP(25);
t125 = t94 * pkin(5);
t124 = (pkin(1) * MDP(6));
t120 = t94 * qJ(6);
t92 = sin(qJ(3));
t75 = t92 * pkin(3) - t94 * qJ(4) + qJ(2);
t66 = t92 * pkin(8) + t75;
t96 = -pkin(1) - pkin(7);
t77 = (pkin(4) - t96) * t94;
t63 = t93 * t66 + t91 * t77;
t60 = t120 + t63;
t123 = t60 * t91;
t122 = t91 * t93;
t86 = t92 ^ 2;
t88 = t94 ^ 2;
t82 = t86 + t88;
t121 = pkin(3) * MDP(17);
t95 = -pkin(3) - pkin(8);
t119 = MDP(28) * t95;
t108 = t91 * t66 - t93 * t77;
t61 = t108 - t125;
t118 = t61 * MDP(28);
t104 = -t91 * pkin(5) + t93 * qJ(6);
t74 = qJ(4) - t104;
t117 = t74 * MDP(28);
t116 = t96 ^ 2 * MDP(17);
t115 = t93 * MDP(20);
t114 = t93 * MDP(28);
t113 = t96 * MDP(17);
t112 = -MDP(13) + MDP(16);
t109 = MDP(19) * t122;
t107 = -MDP(15) + t121;
t78 = pkin(5) * t93 + t91 * qJ(6);
t59 = -t61 * t93 + t123;
t101 = t91 * MDP(20) + t93 * MDP(21);
t100 = -t108 * MDP(23) - t63 * MDP(24);
t99 = -t93 * MDP(23) + t91 * MDP(24);
t98 = -t93 * MDP(25) - t91 * MDP(27);
t83 = t92 * t96;
t80 = t93 * t95 * t94;
t79 = t94 * pkin(3) + t92 * qJ(4);
t76 = -t92 * pkin(4) + t83;
t73 = t82 * t96;
t72 = t81 * t95;
t69 = t81 * t94;
t64 = t83 + (-pkin(4) - t78) * t92;
t1 = [MDP(1) + (t60 ^ 2 + t61 ^ 2 + t64 ^ 2) * MDP(28) + (-0.2e1 * t94 * MDP(16) + MDP(17) * t75) * t75 + (MDP(22) + MDP(7) + t116) * t88 + (t85 * MDP(18) + 0.2e1 * t109 + t116) * t86 + ((-2 * MDP(4) + t124) * pkin(1)) + (MDP(13) * t128 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + (0.2e1 * qJ(2) * MDP(12) + t75 * t127 + (-MDP(8) + t101) * t128) * t92 - 0.2e1 * t73 * MDP(14) + (-t61 * MDP(25) + t60 * MDP(27) + t100) * t128 + 0.2e1 * ((t60 * t93 + t61 * t91) * MDP(26) + t99 * t76 + t98 * t64) * t92; MDP(4) - t124 + t73 * MDP(17) + (-t59 * t94 + t64 * t92) * MDP(28) + t131 * t82; MDP(6) + t82 * MDP(17) + (t81 * t88 + t86) * MDP(28); -t79 * MDP(14) + (t76 * t91 + t80) * MDP(23) + (t64 * t91 + t80) * MDP(25) - t59 * MDP(26) + (t95 * t123 + t64 * t74) * MDP(28) + (t76 * MDP(24) - t64 * MDP(27) - t95 * t118) * t93 + (t115 + MDP(9) + (MDP(12) + t107) * t96 + (-t110 * t95 - MDP(21)) * t91) * t94 + (-MDP(10) + MDP(18) * t122 + (-t85 + t87) * MDP(19) + t112 * t96 + t98 * t74 + (t99 + t113) * qJ(4)) * t92; t79 * MDP(17) + t69 * MDP(26) + (-t81 * t119 + MDP(12) - MDP(15)) * t94 + (t110 * t93 + t111 * t91 + t112 + t117) * t92; -0.2e1 * t109 - 0.2e1 * t72 * MDP(26) + t87 * MDP(18) + MDP(11) + t95 ^ 2 * t130 + (-0.2e1 * t93 * MDP(27) + t91 * t126 + t117) * t74 + (t127 + t121) * pkin(3) + (MDP(17) * qJ(4) + 0.2e1 * t91 * MDP(23) + 0.2e1 * t93 * MDP(24) + 0.2e1 * MDP(16)) * qJ(4); t59 * MDP(28) + (-t113 - t131) * t94; -t94 * MDP(17) - t69 * MDP(28); -t81 * MDP(26) + t72 * MDP(28) - t107; MDP(17) + t130; t94 * MDP(22) + (-t108 + 0.2e1 * t125) * MDP(25) + (0.2e1 * t120 + t63) * MDP(27) + (-t61 * pkin(5) + t60 * qJ(6)) * MDP(28) + (t104 * MDP(26) + t101) * t92 + t100; -t129 * t94; -t91 * MDP(21) - t78 * MDP(26) + t129 * t95 + t115; t78 * MDP(28) + t97; MDP(22) + pkin(5) * t126 + 0.2e1 * qJ(6) * MDP(27) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(28); t91 * t92 * MDP(26) - t94 * MDP(25) + t118; t94 * t114; (MDP(26) - t119) * t93; -t114; t106; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
