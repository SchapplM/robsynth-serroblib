% Calculate joint inertia matrix for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP7_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:33
% EndTime: 2019-12-31 20:01:34
% DurationCPUTime: 0.35s
% Computational Cost: add. (497->114), mult. (901->169), div. (0->0), fcn. (907->6), ass. (0->47)
t107 = cos(qJ(2));
t80 = sin(pkin(8));
t81 = cos(pkin(8));
t83 = sin(qJ(2));
t70 = t80 * t107 + t81 * t83;
t110 = 0.2e1 * t70;
t109 = qJ(3) + pkin(6);
t69 = -t81 * t107 + t80 * t83;
t108 = pkin(4) * t69;
t75 = pkin(2) * t80 + pkin(7);
t106 = t69 * t75;
t77 = -t107 * pkin(2) - pkin(1);
t62 = t69 * pkin(3) - t70 * pkin(7) + t77;
t72 = t109 * t107;
t98 = t109 * t83;
t66 = t81 * t72 - t80 * t98;
t82 = sin(qJ(4));
t84 = cos(qJ(4));
t58 = t82 * t62 + t84 * t66;
t78 = t82 ^ 2;
t79 = t84 ^ 2;
t105 = t78 + t79;
t104 = qJ(5) * t69;
t103 = MDP(14) * t84;
t102 = MDP(21) * t70;
t101 = t69 * MDP(17);
t100 = 0.2e1 * t107;
t99 = -MDP(19) + MDP(22);
t76 = -pkin(2) * t81 - pkin(3);
t97 = -t84 * t62 + t66 * t82;
t64 = t72 * t80 + t81 * t98;
t96 = -MDP(23) * pkin(4) - MDP(20);
t95 = t75 * MDP(23) + MDP(21);
t94 = -pkin(4) * t84 - qJ(5) * t82;
t93 = pkin(4) * t82 - qJ(5) * t84;
t55 = t104 + t58;
t56 = t97 - t108;
t92 = t55 * t82 - t56 * t84;
t67 = t76 + t94;
t91 = t67 * t70 - t106;
t90 = t70 * t76 - t106;
t89 = MDP(18) - t96;
t88 = MDP(23) * qJ(5) + t99;
t87 = t84 * MDP(15) - t82 * MDP(16);
t86 = -MDP(18) * t97 - t58 * MDP(19);
t59 = t93 * t70 + t64;
t1 = [MDP(1) + pkin(1) * MDP(9) * t100 + (t64 ^ 2 + t66 ^ 2 + t77 ^ 2) * MDP(12) + (t55 ^ 2 + t56 ^ 2 + t59 ^ 2) * MDP(23) + (t79 * MDP(13) - 0.2e1 * t82 * t103) * t70 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t83 + MDP(5) * t100) * t83 + (t87 * t110 + t101) * t69 + 0.2e1 * (-t66 * MDP(11) - t56 * MDP(20) + t55 * MDP(22) + t86) * t69 + (-t92 * MDP(21) + (t82 * MDP(20) - t84 * MDP(22)) * t59 + (t82 * MDP(18) + t84 * MDP(19) + MDP(11)) * t64) * t110; t59 * t67 * MDP(23) + t83 * MDP(6) + t107 * MDP(7) + (-t78 + t79) * MDP(14) * t70 + (-t107 * MDP(10) - t83 * MDP(9)) * pkin(6) + (t69 * MDP(16) - t64 * MDP(18) + t90 * MDP(19) - t59 * MDP(20) - t91 * MDP(22) + t95 * t55) * t84 + (t84 * t70 * MDP(13) + t69 * MDP(15) + t90 * MDP(18) + t64 * MDP(19) + t91 * MDP(20) - t59 * MDP(22) + t95 * t56) * t82 + ((-t69 * t80 - t70 * t81) * MDP(11) + (-t64 * t81 + t66 * t80) * MDP(12)) * pkin(2); MDP(8) + t78 * MDP(13) + (t105 * t75 ^ 2 + t67 ^ 2) * MDP(23) + (t80 ^ 2 + t81 ^ 2) * MDP(12) * pkin(2) ^ 2 + 0.2e1 * t105 * MDP(21) * t75 + 0.2e1 * (-t76 * MDP(18) - t67 * MDP(20)) * t84 + 0.2e1 * (t76 * MDP(19) - t67 * MDP(22) + t103) * t82; t77 * MDP(12) + t92 * MDP(23) - t105 * t102 + ((MDP(18) + MDP(20)) * t84 + t99 * t82) * t69; 0; t105 * MDP(23) + MDP(12); t101 + (-t97 + 0.2e1 * t108) * MDP(20) + (0.2e1 * t104 + t58) * MDP(22) + (-pkin(4) * t56 + qJ(5) * t55) * MDP(23) + (t94 * MDP(21) + t87) * t70 + t86; t82 * MDP(15) + t84 * MDP(16) - t93 * MDP(21) + (-t89 * t82 + t88 * t84) * t75; t88 * t82 + t89 * t84; MDP(17) + 0.2e1 * pkin(4) * MDP(20) + 0.2e1 * qJ(5) * MDP(22) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(23); -t69 * MDP(20) + t56 * MDP(23) + t84 * t102; t95 * t82; -t84 * MDP(23); t96; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
