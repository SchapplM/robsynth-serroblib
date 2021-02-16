% Calculate joint inertia matrix for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:51
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRRP1_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:50:14
% EndTime: 2021-01-16 00:50:16
% DurationCPUTime: 0.48s
% Computational Cost: add. (549->153), mult. (1358->230), div. (0->0), fcn. (1549->12), ass. (0->65)
t102 = sin(qJ(4));
t137 = 0.2e1 * t102;
t101 = sin(qJ(5));
t136 = pkin(9) * t101;
t105 = cos(qJ(4));
t135 = pkin(9) * t105;
t95 = sin(pkin(12));
t97 = sin(pkin(6));
t134 = t95 * t97;
t98 = cos(pkin(12));
t133 = t97 * t98;
t132 = -qJ(6) - pkin(10);
t131 = MDP(23) * pkin(5);
t103 = sin(qJ(3));
t96 = sin(pkin(7));
t130 = t103 * t96;
t106 = cos(qJ(3));
t129 = t106 * t96;
t128 = t101 * t102;
t104 = cos(qJ(5));
t127 = t102 * t104;
t88 = t132 * t104;
t126 = t88 * MDP(21);
t90 = -t104 * pkin(5) - pkin(4);
t125 = t90 * MDP(23);
t124 = t101 * MDP(16);
t123 = t101 * MDP(21);
t122 = t102 * MDP(12);
t121 = t104 * MDP(20);
t120 = MDP(18) + MDP(20);
t119 = MDP(19) + MDP(21);
t118 = t104 * t135;
t117 = t101 * t104 * MDP(14);
t116 = -MDP(22) * pkin(5) + MDP(15);
t115 = MDP(20) + t131;
t114 = -t105 * MDP(11) - MDP(4);
t113 = MDP(18) + t115;
t112 = MDP(18) * t101 + MDP(19) * t104;
t111 = t101 * MDP(20) + t104 * MDP(21);
t100 = cos(pkin(6));
t99 = cos(pkin(7));
t110 = t100 * t96 + t99 * t133;
t86 = -t105 * pkin(4) - t102 * pkin(10) - pkin(3);
t78 = t118 + (-qJ(6) * t102 + t86) * t101;
t84 = t104 * t86;
t109 = (-t101 * t135 + t84) * MDP(18) - (t101 * t86 + t118) * MDP(19) - t78 * MDP(21);
t108 = -t121 + t123 + t125;
t107 = t119 * t101 - t120 * t104 - MDP(11) + t125;
t94 = t104 ^ 2;
t92 = t101 ^ 2;
t87 = t132 * t101;
t85 = (pkin(5) * t101 + pkin(9)) * t102;
t83 = t102 * t99 + t105 * t130;
t82 = t102 * t130 - t105 * t99;
t81 = t100 * t99 - t96 * t133;
t77 = -t101 * t129 + t104 * t83;
t76 = -t101 * t83 - t104 * t129;
t75 = t110 * t103 + t106 * t134;
t74 = t103 * t134 - t110 * t106;
t73 = -qJ(6) * t127 + t84 + (-pkin(5) - t136) * t105;
t70 = t102 * t81 + t105 * t75;
t69 = t102 * t75 - t105 * t81;
t68 = t74 * t101 + t104 * t70;
t67 = -t101 * t70 + t74 * t104;
t1 = [MDP(1) + (t100 ^ 2 + (t95 ^ 2 + t98 ^ 2) * t97 ^ 2) * MDP(2) + (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) * MDP(23); t100 * MDP(2) + (t67 * t76 + t68 * t77 + t69 * t82) * MDP(23); MDP(2) + (t76 ^ 2 + t77 ^ 2 + t82 ^ 2) * MDP(23); -t75 * MDP(5) + (t67 * t73 + t68 * t78 + t69 * t85) * MDP(23) + t114 * t74 + t119 * (t68 * t105 + t69 * t127) + t120 * (-t67 * t105 + t69 * t128) + (t74 * MDP(12) + (-t101 * t68 - t104 * t67) * MDP(22)) * t102; (t76 * t73 + t77 * t78 + t82 * t85) * MDP(23) + t119 * (t77 * t105 + t82 * t127) + t120 * (-t76 * t105 + t82 * t128) + (-t101 * t77 - t104 * t76) * MDP(22) * t102 + (-t103 * MDP(5) + (-t114 - t122) * t106) * t96; MDP(3) - 0.2e1 * pkin(3) * t122 + (t73 ^ 2 + t78 ^ 2 + t85 ^ 2) * MDP(23) + (0.2e1 * pkin(3) * MDP(11) + t105 * MDP(17) + (-t104 * MDP(15) + MDP(7) + t124) * t137) * t105 + 0.2e1 * (-t73 * MDP(20) - t109) * t105 + ((-t101 * t78 - t104 * t73) * MDP(22) + t111 * t85) * t137 + (t94 * MDP(13) + 0.2e1 * t112 * pkin(9) + MDP(6) - 0.2e1 * t117) * t102 ^ 2; -t70 * MDP(12) + (-t67 * t101 + t68 * t104) * MDP(22) + (t67 * t87 - t68 * t88) * MDP(23) + t107 * t69; -t83 * MDP(12) + (-t76 * t101 + t77 * t104) * MDP(22) + (t76 * t87 - t77 * t88) * MDP(23) + t107 * t82; (-t73 * t101 + t78 * t104) * MDP(22) + (t73 * t87 - t78 * t88) * MDP(23) + t108 * t85 + (-pkin(9) * MDP(12) - t101 * MDP(15) - t104 * MDP(16) - t87 * MDP(20) + t112 * pkin(10) + MDP(9) - t126) * t105 + (MDP(8) - pkin(9) * MDP(11) + (-t92 + t94) * MDP(14) + (-pkin(9) * MDP(18) - pkin(4) * MDP(19) + t90 * MDP(21) - t87 * MDP(22)) * t104 + (t104 * MDP(13) - pkin(4) * MDP(18) + pkin(9) * MDP(19) + t90 * MDP(20) + t88 * MDP(22)) * t101) * t102; MDP(10) + t92 * MDP(13) + 0.2e1 * t117 + 0.2e1 * (-t87 * t101 - t88 * t104) * MDP(22) + (t87 ^ 2 + t88 ^ 2) * MDP(23) + (-0.2e1 * t121 + 0.2e1 * t123 + t125) * t90 + 0.2e1 * (t104 * MDP(18) - t101 * MDP(19)) * pkin(4); t113 * t67 - t119 * t68; t113 * t76 - t119 * t77; t73 * t131 + t84 * MDP(20) + (-MDP(17) + (-0.2e1 * pkin(5) - t136) * MDP(20)) * t105 + (-t124 + (-MDP(20) * qJ(6) + t116) * t104) * t102 + t109; t126 + t115 * t87 + (-MDP(19) * pkin(10) + MDP(16)) * t104 + (-MDP(18) * pkin(10) + t116) * t101; MDP(17) + (0.2e1 * MDP(20) + t131) * pkin(5); t69 * MDP(23); t82 * MDP(23); t85 * MDP(23) + t111 * t102; t108; 0; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
