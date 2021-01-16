% Calculate joint inertia matrix for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRPR1_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:06:21
% EndTime: 2021-01-16 01:06:23
% DurationCPUTime: 0.33s
% Computational Cost: add. (445->117), mult. (951->194), div. (0->0), fcn. (1077->12), ass. (0->58)
t101 = sin(qJ(2));
t103 = cos(qJ(2));
t94 = sin(pkin(11));
t95 = sin(pkin(6));
t97 = cos(pkin(11));
t76 = (t101 * t94 - t103 * t97) * t95;
t129 = t76 ^ 2;
t127 = cos(qJ(4));
t90 = -t97 * pkin(2) - pkin(3);
t85 = -t127 * pkin(4) + t90;
t128 = 0.2e1 * t85;
t126 = MDP(5) * pkin(2);
t100 = sin(qJ(4));
t93 = sin(pkin(12));
t96 = cos(pkin(12));
t84 = t96 * t100 + t93 * t127;
t99 = sin(qJ(6));
t125 = t84 * t99;
t124 = MDP(16) * pkin(4);
t102 = cos(qJ(6));
t123 = t102 * t84;
t122 = MDP(18) * t99;
t82 = t100 * t93 - t96 * t127;
t121 = MDP(21) * t82;
t120 = t84 * MDP(14);
t119 = t99 * MDP(23);
t118 = pkin(2) * t94 + pkin(8);
t117 = t127 * MDP(11);
t116 = t93 * t124 - MDP(14);
t88 = pkin(4) * t93 + pkin(9);
t89 = -pkin(4) * t96 - pkin(5);
t115 = -t82 * t88 + t84 * t89;
t114 = t85 * MDP(16) + t120;
t113 = t102 * MDP(22) - t119;
t112 = (-qJ(5) - t118) * t100;
t111 = t127 * t118;
t78 = (t101 * t97 + t103 * t94) * t95;
t98 = cos(pkin(6));
t110 = -t100 * t78 + t98 * t127;
t109 = MDP(13) + t113;
t108 = -t100 * MDP(12) + t117;
t107 = (MDP(19) * t102 - MDP(20) * t99) * t84;
t106 = -t96 * t124 - t109;
t92 = t102 ^ 2;
t91 = t99 ^ 2;
t81 = t84 ^ 2;
t80 = t127 * qJ(5) + t111;
t75 = t98 * t100 + t127 * t78;
t73 = t93 * t112 + t96 * t80;
t71 = -t96 * t112 + t80 * t93;
t70 = t82 * pkin(5) - t84 * pkin(9) + t85;
t69 = t93 * t110 + t96 * t75;
t67 = -t96 * t110 + t75 * t93;
t66 = t102 * t73 + t70 * t99;
t65 = t102 * t70 - t73 * t99;
t64 = t102 * t69 + t76 * t99;
t63 = t102 * t76 - t69 * t99;
t1 = [MDP(1) + (t78 ^ 2 + t98 ^ 2 + t129) * MDP(5) + (t67 ^ 2 + t69 ^ 2 + t129) * MDP(16); t78 * t94 * t126 + (t67 * t84 - t69 * t82) * MDP(15) + (t67 * t71 + t69 * t73) * MDP(16) + (t67 * t125 + t63 * t82) * MDP(22) + (t67 * t123 - t64 * t82) * MDP(23) + (MDP(3) * t103 - MDP(4) * t101) * t95 + (t82 * MDP(13) - t97 * t126 - t108 + t114) * t76; MDP(2) - 0.2e1 * t90 * t117 + t120 * t128 + (t71 ^ 2 + t73 ^ 2 + t85 ^ 2) * MDP(16) + (t92 * MDP(17) - 0.2e1 * t102 * t122) * t81 + (t94 ^ 2 + t97 ^ 2) * MDP(5) * pkin(2) ^ 2 + (0.2e1 * t90 * MDP(12) + MDP(6) * t100 + 0.2e1 * t127 * MDP(7)) * t100 + (MDP(13) * t128 + 0.2e1 * t107 + t121) * t82 + 0.2e1 * (t71 * t84 - t73 * t82) * MDP(15) + 0.2e1 * (t71 * t125 + t65 * t82) * MDP(22) + 0.2e1 * (t71 * t123 - t66 * t82) * MDP(23); t98 * MDP(5) + (t67 * t82 + t69 * t84) * MDP(16); (t71 * t82 + t73 * t84) * MDP(16); MDP(5) + (t82 ^ 2 + t81) * MDP(16); t110 * MDP(11) - t75 * MDP(12) + t106 * t67 + t116 * t69; t127 * MDP(9) - MDP(12) * t111 - t71 * MDP(13) - t73 * MDP(14) + (-t91 + t92) * MDP(18) * t84 + (-t118 * MDP(11) + MDP(8)) * t100 + (t82 * MDP(19) + t115 * MDP(22) + t71 * MDP(23)) * t99 + (MDP(17) * t125 + t82 * MDP(20) - t71 * MDP(22) + t115 * MDP(23)) * t102 + ((-t82 * t93 - t84 * t96) * MDP(15) + (-t71 * t96 + t73 * t93) * MDP(16)) * pkin(4); t106 * t82 + t116 * t84 + t108; 0.2e1 * t89 * t119 + t91 * MDP(17) + MDP(10) + (t93 ^ 2 + t96 ^ 2) * MDP(16) * pkin(4) ^ 2 + 0.2e1 * (-MDP(22) * t89 + t122) * t102 + 0.2e1 * (MDP(13) * t96 - MDP(14) * t93) * pkin(4); t76 * MDP(16); t109 * t82 + t114; 0; 0; MDP(16); MDP(22) * t63 - MDP(23) * t64; t65 * MDP(22) - t66 * MDP(23) + t107 + t121; (-MDP(22) * t99 - MDP(23) * t102) * t84; (-t88 * MDP(22) + MDP(19)) * t99 + (-t88 * MDP(23) + MDP(20)) * t102; t113; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
