% Calculate joint inertia matrix for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRRP1_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:24:20
% EndTime: 2021-01-16 01:24:23
% DurationCPUTime: 0.48s
% Computational Cost: add. (447->141), mult. (959->207), div. (0->0), fcn. (985->10), ass. (0->60)
t113 = MDP(19) + MDP(21);
t114 = MDP(18) + MDP(20);
t96 = cos(qJ(5));
t84 = -t96 * pkin(5) - pkin(4);
t119 = t84 * MDP(23);
t93 = sin(qJ(5));
t128 = t113 * t93 - t114 * t96 - MDP(11) + t119;
t89 = sin(pkin(11));
t90 = sin(pkin(6));
t91 = cos(pkin(11));
t95 = sin(qJ(2));
t98 = cos(qJ(2));
t74 = (t89 * t98 + t91 * t95) * t90;
t92 = cos(pkin(6));
t94 = sin(qJ(4));
t97 = cos(qJ(4));
t70 = t94 * t74 - t92 * t97;
t127 = t70 * t94;
t81 = t89 * pkin(2) + pkin(8);
t126 = t81 * t93;
t125 = t81 * t97;
t124 = t93 * t96;
t123 = -qJ(6) - pkin(9);
t122 = MDP(23) * pkin(5);
t121 = qJ(6) * t94;
t79 = t123 * t96;
t120 = t79 * MDP(21);
t118 = t93 * MDP(16);
t117 = t93 * MDP(21);
t116 = t96 * MDP(20);
t115 = t97 * MDP(11);
t112 = t96 * t125;
t111 = MDP(14) * t124;
t82 = -t91 * pkin(2) - pkin(3);
t110 = -MDP(22) * pkin(5) + MDP(15);
t109 = MDP(20) + t122;
t71 = t97 * t74 + t92 * t94;
t72 = (t89 * t95 - t91 * t98) * t90;
t64 = -t93 * t71 + t96 * t72;
t65 = t96 * t71 + t93 * t72;
t107 = -t64 * t93 + t65 * t96;
t77 = -t97 * pkin(4) - t94 * pkin(9) + t82;
t75 = t96 * t77;
t66 = -t96 * t121 + t75 + (-pkin(5) - t126) * t97;
t67 = t112 + (t77 - t121) * t93;
t106 = -t66 * t93 + t67 * t96;
t78 = t123 * t93;
t105 = -t78 * t93 - t79 * t96;
t104 = MDP(18) + t109;
t103 = MDP(18) * t93 + MDP(19) * t96;
t102 = t93 * MDP(20) + t96 * MDP(21);
t101 = (-t93 * t125 + t75) * MDP(18) - (t93 * t77 + t112) * MDP(19) - t67 * MDP(21);
t100 = -t116 + t117 + t119;
t88 = t97 ^ 2;
t87 = t96 ^ 2;
t86 = t94 ^ 2;
t85 = t93 ^ 2;
t83 = t87 * t94;
t76 = (pkin(5) * t93 + t81) * t94;
t1 = [MDP(1) + (t72 ^ 2 + t74 ^ 2 + t92 ^ 2) * MDP(5) + (t64 ^ 2 + t65 ^ 2 + t70 ^ 2) * MDP(23); -t72 * t115 + (t64 * t66 + t65 * t67 + t70 * t76) * MDP(23) + (t98 * MDP(3) - t95 * MDP(4)) * t90 + t114 * (t93 * t127 - t64 * t97) + t113 * (t96 * t127 + t65 * t97) + (-t72 * t91 + t74 * t89) * MDP(5) * pkin(2) + (t72 * MDP(12) + (-t64 * t96 - t65 * t93) * MDP(22)) * t94; MDP(2) - 0.2e1 * t82 * t115 + t88 * MDP(17) + (t66 ^ 2 + t67 ^ 2 + t76 ^ 2) * MDP(23) + (t89 ^ 2 + t91 ^ 2) * MDP(5) * pkin(2) ^ 2 + 0.2e1 * (-t66 * MDP(20) - t101) * t97 + (t87 * MDP(13) + 0.2e1 * t103 * t81 + MDP(6) - 0.2e1 * t111) * t86 + 0.2e1 * (t82 * MDP(12) + (-t96 * MDP(15) + MDP(7) + t118) * t97 + (-t66 * t96 - t67 * t93) * MDP(22) + t102 * t76) * t94; t92 * MDP(5) + (t107 * t94 - t70 * t97) * MDP(23); (t106 * t94 - t76 * t97) * MDP(23); MDP(5) + (t88 + (t85 + t87) * t86) * MDP(23); -t71 * MDP(12) + t107 * MDP(22) + (t64 * t78 - t65 * t79) * MDP(23) + t128 * t70; t83 * MDP(14) + t106 * MDP(22) + (t66 * t78 - t67 * t79) * MDP(23) + t100 * t76 + (-t81 * MDP(12) - t93 * MDP(15) - t96 * MDP(16) - t78 * MDP(20) + t103 * pkin(9) + MDP(9) - t120) * t97 + (MDP(8) - t81 * MDP(11) + MDP(13) * t124 - t85 * MDP(14) + (-pkin(4) * t93 - t81 * t96) * MDP(18) + (-pkin(4) * t96 + t126) * MDP(19) + (-t78 * t96 + t79 * t93) * MDP(22) + t102 * t84) * t94; t83 * MDP(22) + (t85 * MDP(22) + t105 * MDP(23) - MDP(12)) * t94 - t128 * t97; MDP(10) + t85 * MDP(13) + 0.2e1 * t111 + 0.2e1 * t105 * MDP(22) + (t78 ^ 2 + t79 ^ 2) * MDP(23) + (-0.2e1 * t116 + 0.2e1 * t117 + t119) * t84 + 0.2e1 * (t96 * MDP(18) - t93 * MDP(19)) * pkin(4); t104 * t64 - t113 * t65; t66 * t122 + t75 * MDP(20) + (-MDP(17) + (-0.2e1 * pkin(5) - t126) * MDP(20)) * t97 + (-t118 + (-MDP(20) * qJ(6) + t110) * t96) * t94 + t101; (-t104 * t93 - t113 * t96) * t94; t120 + (-MDP(19) * pkin(9) + MDP(16)) * t96 + t109 * t78 + (-MDP(18) * pkin(9) + t110) * t93; MDP(17) + (0.2e1 * MDP(20) + t122) * pkin(5); t70 * MDP(23); t76 * MDP(23) + t102 * t94; -t97 * MDP(23); t100; 0; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
