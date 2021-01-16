% Calculate joint inertia matrix for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP6_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:29:59
% EndTime: 2021-01-15 20:30:02
% DurationCPUTime: 0.39s
% Computational Cost: add. (511->123), mult. (938->172), div. (0->0), fcn. (944->6), ass. (0->52)
t79 = sin(pkin(8));
t80 = cos(pkin(8));
t82 = sin(qJ(2));
t84 = cos(qJ(2));
t69 = t79 * t84 + t80 * t82;
t112 = 0.2e1 * t69;
t76 = -pkin(2) * t84 - pkin(1);
t111 = 0.2e1 * t76;
t110 = 0.2e1 * t84;
t68 = t79 * t82 - t80 * t84;
t109 = pkin(4) * t68;
t107 = -qJ(3) - pkin(6);
t72 = t107 * t84;
t95 = t107 * t82;
t64 = -t80 * t72 + t79 * t95;
t83 = cos(qJ(4));
t108 = t64 * t83;
t81 = sin(qJ(4));
t77 = t81 ^ 2;
t78 = t83 ^ 2;
t106 = t77 + t78;
t105 = MDP(25) * pkin(4);
t104 = qJ(5) * t69;
t74 = pkin(2) * t79 + pkin(7);
t103 = qJ(5) + t74;
t75 = -pkin(2) * t80 - pkin(3);
t71 = -pkin(4) * t83 + t75;
t102 = MDP(25) * t71;
t66 = t103 * t83;
t101 = t66 * MDP(23);
t100 = t68 * MDP(19);
t99 = t81 * MDP(18);
t98 = t81 * MDP(23);
t97 = t83 * MDP(22);
t96 = t81 * t83 * MDP(16);
t61 = pkin(3) * t68 - pkin(7) * t69 + t76;
t57 = t83 * t61 - t64 * t81;
t62 = -t72 * t79 - t80 * t95;
t94 = -MDP(24) * pkin(4) + MDP(17);
t93 = MDP(22) + t105;
t92 = (-MDP(21) - MDP(23)) * t81;
t55 = -t83 * t104 + t109 + t57;
t56 = t108 + (t61 - t104) * t81;
t91 = t55 * t83 + t56 * t81;
t90 = -t83 * MDP(20) + t81 * MDP(21);
t89 = -t81 * MDP(20) - t83 * MDP(21);
t88 = t81 * MDP(22) + t83 * MDP(23);
t87 = t57 * MDP(20) - (t61 * t81 + t108) * MDP(21) - t56 * MDP(23);
t86 = -t97 + t98 + t102;
t65 = t103 * t81;
t59 = pkin(4) * t69 * t81 + t62;
t1 = [MDP(1) + pkin(1) * MDP(9) * t110 + (t62 ^ 2 + t64 ^ 2 + t76 ^ 2) * MDP(14) + (t55 ^ 2 + t56 ^ 2 + t59 ^ 2) * MDP(25) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t82 + MDP(5) * t110) * t82 + (MDP(11) * t111 + t100 + (MDP(17) * t83 - t99) * t112) * t68 + 0.2e1 * (-t64 * MDP(13) + t55 * MDP(22) + t87) * t68 + (-t91 * MDP(24) + t88 * t59 + (MDP(13) - t89) * t62) * t112 + (MDP(12) * t111 + (t78 * MDP(15) - 0.2e1 * t96) * t69) * t69; t82 * MDP(6) + t84 * MDP(7) - t64 * MDP(12) + (-t55 * t81 + t56 * t83) * MDP(24) + (-t55 * t65 + t56 * t66) * MDP(25) + t86 * t59 + (-MDP(11) + t90) * t62 + (t81 * MDP(17) + t83 * MDP(18) - t65 * MDP(22) + t89 * t74 - t101) * t68 + (-MDP(10) * t84 - MDP(9) * t82) * pkin(6) + ((-t77 + t78) * MDP(16) + (MDP(21) * t75 + MDP(23) * t71 + MDP(24) * t65) * t83 + (MDP(15) * t83 + MDP(20) * t75 + MDP(22) * t71 - MDP(24) * t66) * t81) * t69 + ((-t68 * t79 - t69 * t80) * MDP(13) + (-t62 * t80 + t64 * t79) * MDP(14)) * pkin(2); MDP(8) + t77 * MDP(15) + 0.2e1 * t96 + 0.2e1 * (t65 * t81 + t66 * t83) * MDP(24) + (t65 ^ 2 + t66 ^ 2) * MDP(25) + (t79 ^ 2 + t80 ^ 2) * MDP(14) * pkin(2) ^ 2 + (-0.2e1 * t97 + 0.2e1 * t98 + t102) * t71 + 0.2e1 * t90 * t75 + 0.2e1 * (MDP(11) * t80 - MDP(12) * t79) * pkin(2); t76 * MDP(14) + t91 * MDP(25) + (-t106 * MDP(24) + MDP(12)) * t69 + (MDP(11) + (MDP(20) + MDP(22)) * t83 + t92) * t68; (-t65 * t83 + t66 * t81) * MDP(25); t106 * MDP(25) + MDP(14); t100 + (t57 + 0.2e1 * t109) * MDP(22) + t55 * t105 + (-t99 + (-MDP(22) * qJ(5) + t94) * t83) * t69 + t87; -t101 + (-MDP(21) * t74 + MDP(18)) * t83 - t93 * t65 + (-MDP(20) * t74 + t94) * t81; t92 + (MDP(20) + t93) * t83; MDP(19) + (0.2e1 * MDP(22) + t105) * pkin(4); MDP(25) * t59 + t88 * t69; t86; 0; 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
