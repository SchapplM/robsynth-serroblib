% Calculate joint inertia matrix for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRRP5_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:51:52
% EndTime: 2021-01-16 01:51:55
% DurationCPUTime: 0.46s
% Computational Cost: add. (375->145), mult. (737->195), div. (0->0), fcn. (691->8), ass. (0->61)
t108 = MDP(21) + MDP(23);
t90 = cos(qJ(5));
t130 = t108 * t90;
t109 = MDP(20) + MDP(22);
t79 = -t90 * pkin(5) - pkin(4);
t114 = t79 * MDP(25);
t87 = sin(qJ(5));
t129 = t108 * t87 - t109 * t90 - MDP(13) + t114;
t128 = (pkin(2) * MDP(7));
t85 = sin(pkin(6));
t92 = cos(qJ(2));
t125 = t85 * t92;
t86 = cos(pkin(6));
t88 = sin(qJ(4));
t91 = cos(qJ(4));
t70 = t91 * t125 + t86 * t88;
t127 = t70 * t91;
t89 = sin(qJ(2));
t126 = t85 * t89;
t124 = t87 * t90;
t93 = -pkin(2) - pkin(8);
t123 = t87 * t93;
t122 = t90 * t93;
t121 = -qJ(6) - pkin(9);
t81 = t87 ^ 2;
t83 = t90 ^ 2;
t120 = t81 + t83;
t118 = MDP(25) * pkin(5);
t117 = qJ(6) * t91;
t116 = qJ(3) * MDP(7);
t78 = t121 * t90;
t115 = t78 * MDP(23);
t113 = t87 * MDP(18);
t112 = t87 * MDP(23);
t111 = t90 * MDP(22);
t110 = t91 * MDP(14);
t107 = t88 * t122;
t106 = MDP(16) * t124;
t105 = MDP(5) - t128;
t104 = -MDP(24) * pkin(5) + MDP(17);
t103 = MDP(22) + t118;
t76 = t88 * pkin(4) - t91 * pkin(9) + qJ(3);
t72 = t90 * t76;
t64 = -t90 * t117 + t72 + (pkin(5) - t123) * t88;
t65 = t107 + (t76 - t117) * t87;
t101 = -t64 * t87 + t65 * t90;
t71 = -t88 * t125 + t86 * t91;
t66 = t90 * t126 - t71 * t87;
t67 = t87 * t126 + t71 * t90;
t100 = -t66 * t87 + t67 * t90;
t77 = t121 * t87;
t99 = -t77 * t87 - t78 * t90;
t98 = MDP(20) + t103;
t97 = -MDP(20) * t87 - MDP(21) * t90;
t96 = t87 * MDP(22) + t90 * MDP(23);
t95 = (-t88 * t123 + t72) * MDP(20) - (t87 * t76 + t107) * MDP(21) - t65 * MDP(23);
t94 = -t111 + t112 + t114;
t84 = t91 ^ 2;
t82 = t88 ^ 2;
t73 = (pkin(5) * t87 - t93) * t91;
t1 = [MDP(1) + (t86 ^ 2 + (t89 ^ 2 + t92 ^ 2) * t85 ^ 2) * MDP(7) + (t66 ^ 2 + t67 ^ 2 + t70 ^ 2) * MDP(25); (t66 * t64 + t67 * t65 + t70 * t73) * MDP(25) + (-t66 * t90 - t67 * t87) * MDP(24) * t91 + t108 * (t90 * t127 - t67 * t88) + t109 * (t87 * t127 + t66 * t88) + ((MDP(3) - t105) * t92 + (t88 * MDP(13) - MDP(4) + MDP(6) + t110 + t116) * t89) * t85; MDP(2) + t82 * MDP(19) + (t64 ^ 2 + t65 ^ 2 + t73 ^ 2) * MDP(25) + ((-2 * MDP(5) + t128) * pkin(2)) + (0.2e1 * MDP(6) + 0.2e1 * t110 + t116) * qJ(3) + 0.2e1 * ((-t64 * t90 - t65 * t87) * MDP(24) + t96 * t73) * t91 + (t83 * MDP(15) + 0.2e1 * t97 * t93 + MDP(8) - 0.2e1 * t106) * t84 + 0.2e1 * (qJ(3) * MDP(13) + (t90 * MDP(17) - MDP(9) - t113) * t91 + t64 * MDP(22) + t95) * t88; -MDP(7) * t125 + (t100 * t88 - t127) * MDP(25); (t101 * t88 - t73 * t91) * MDP(25) + t105 + (t109 * t87 + t130) * (-t82 - t84); MDP(7) + (t120 * t82 + t84) * MDP(25); -t71 * MDP(14) + t100 * MDP(24) + (t66 * t77 - t67 * t78) * MDP(25) + t129 * t70; t101 * MDP(24) + (t64 * t77 - t65 * t78) * MDP(25) + t94 * t73 + (-t93 * MDP(14) + t87 * MDP(17) + t90 * MDP(18) + t77 * MDP(22) + t97 * pkin(9) - MDP(11) + t115) * t88 + (MDP(10) + t93 * MDP(13) + MDP(15) * t124 + (-t81 + t83) * MDP(16) + (-pkin(4) * t87 + t122) * MDP(20) + (-pkin(4) * t90 - t123) * MDP(21) + (-t77 * t90 + t78 * t87) * MDP(24) + t96 * t79) * t91; (t120 * MDP(24) + t99 * MDP(25) - MDP(14)) * t88 - t129 * t91; MDP(12) + t81 * MDP(15) + 0.2e1 * t106 + 0.2e1 * t99 * MDP(24) + (t77 ^ 2 + t78 ^ 2) * MDP(25) + (-0.2e1 * t111 + 0.2e1 * t112 + t114) * t79 + 0.2e1 * (t90 * MDP(20) - t87 * MDP(21)) * pkin(4); -t108 * t67 + t98 * t66; t64 * t118 + t72 * MDP(22) + (MDP(19) + (0.2e1 * pkin(5) - t123) * MDP(22)) * t88 + (-t113 + (-MDP(22) * qJ(6) + t104) * t90) * t91 + t95; (-t98 * t87 - t130) * t88; t115 + (-MDP(21) * pkin(9) + MDP(18)) * t90 + t103 * t77 + (-MDP(20) * pkin(9) + t104) * t87; MDP(19) + (0.2e1 * MDP(22) + t118) * pkin(5); t70 * MDP(25); t73 * MDP(25) + t96 * t91; -t91 * MDP(25); t94; 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
