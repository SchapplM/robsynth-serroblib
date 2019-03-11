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
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRRP5_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:16:44
% EndTime: 2019-03-08 20:16:45
% DurationCPUTime: 0.29s
% Computational Cost: add. (288->125), mult. (574->180), div. (0->0), fcn. (535->8), ass. (0->58)
t110 = 2 * MDP(22);
t109 = (pkin(2) * MDP(7));
t72 = sin(pkin(6));
t76 = sin(qJ(2));
t108 = t72 * t76;
t79 = cos(qJ(2));
t107 = t72 * t79;
t74 = sin(qJ(5));
t77 = cos(qJ(5));
t106 = t74 * t77;
t80 = -pkin(2) - pkin(8);
t105 = t74 * t80;
t104 = t77 * t80;
t103 = -qJ(6) - pkin(9);
t68 = t74 ^ 2;
t70 = t77 ^ 2;
t102 = t68 + t70;
t75 = sin(qJ(4));
t69 = t75 ^ 2;
t78 = cos(qJ(4));
t71 = t78 ^ 2;
t101 = -t69 - t71;
t100 = MDP(23) * pkin(5);
t99 = qJ(6) * t78;
t98 = MDP(13) * t75;
t97 = MDP(18) * t74;
t96 = MDP(23) * t75;
t95 = MDP(7) * qJ(3);
t73 = cos(pkin(6));
t60 = -t75 * t107 + t73 * t78;
t56 = t74 * t108 + t60 * t77;
t94 = t56 * MDP(21);
t66 = -pkin(5) * t77 - pkin(4);
t93 = t66 * MDP(23);
t92 = t77 * MDP(21);
t91 = t78 * MDP(23);
t90 = t75 * t104;
t89 = MDP(16) * t106;
t88 = MDP(5) - t109;
t87 = -MDP(22) * pkin(5) + MDP(17);
t86 = MDP(20) + t100;
t55 = t77 * t108 - t60 * t74;
t85 = -t55 * t74 + t56 * t77;
t64 = t103 * t74;
t65 = t103 * t77;
t84 = -t64 * t74 - t65 * t77;
t83 = MDP(20) * t77 - MDP(21) * t74;
t82 = MDP(20) * t74 + t92;
t81 = MDP(13) + t83 - t93;
t63 = pkin(4) * t75 - pkin(9) * t78 + qJ(3);
t62 = (pkin(5) * t74 - t80) * t78;
t61 = t77 * t63;
t59 = t78 * t107 + t73 * t75;
t58 = t74 * t63 + t90;
t57 = -t75 * t105 + t61;
t54 = t90 + (t63 - t99) * t74;
t53 = -t77 * t99 + t61 + (pkin(5) - t105) * t75;
t1 = [MDP(1) + (t73 ^ 2 + (t76 ^ 2 + t79 ^ 2) * t72 ^ 2) * MDP(7) + (t55 ^ 2 + t56 ^ 2 + t59 ^ 2) * MDP(23); (t53 * t55 + t54 * t56 + t59 * t62) * MDP(23) + (MDP(20) * t55 - t94) * t75 + ((-t55 * t77 - t56 * t74) * MDP(22) + t82 * t59) * t78 + ((MDP(3) - t88) * t79 + (MDP(14) * t78 - MDP(4) + MDP(6) + t95 + t98) * t76) * t72; MDP(2) + t69 * MDP(19) + (t53 ^ 2 + t54 ^ 2 + t62 ^ 2) * MDP(23) + ((-2 * MDP(5) + t109) * pkin(2)) + (0.2e1 * MDP(6) + t95 + 0.2e1 * t98) * qJ(3) + (MDP(15) * t70 + MDP(8) - 0.2e1 * t89) * t71 + 0.2e1 * (-t71 * t105 + t57 * t75) * MDP(20) + 0.2e1 * (-t71 * t104 - t58 * t75) * MDP(21) + (0.2e1 * qJ(3) * MDP(14) + (-t53 * t77 - t54 * t74) * t110 + 0.2e1 * (MDP(17) * t77 - MDP(9) - t97) * t75) * t78; -MDP(7) * t107 + (-t59 * t78 + t85 * t75) * MDP(23); -t62 * t91 + (t101 * MDP(21) + t54 * t96) * t77 + (t101 * MDP(20) - t53 * t96) * t74 + t88; MDP(7) + (t102 * t69 + t71) * MDP(23); -t60 * MDP(14) + t85 * MDP(22) + (t55 * t64 - t56 * t65) * MDP(23) - t81 * t59; (-t53 * t74 + t54 * t77) * MDP(22) + (t53 * t64 - t54 * t65 + t62 * t66) * MDP(23) + (-t80 * MDP(14) + t74 * MDP(17) + t77 * MDP(18) - t82 * pkin(9) - MDP(11)) * t75 + (MDP(10) + t80 * MDP(13) + MDP(15) * t106 + (-t68 + t70) * MDP(16) + (-pkin(4) * t74 + t104) * MDP(20) + (-pkin(4) * t77 - t105) * MDP(21) + (-t64 * t77 + t65 * t74) * MDP(22)) * t78; t81 * t78 + (t102 * MDP(22) + t84 * MDP(23) - MDP(14)) * t75; MDP(12) + t68 * MDP(15) + 0.2e1 * t89 + t84 * t110 + (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) * MDP(23) + 0.2e1 * t83 * pkin(4); t86 * t55 - t94; t53 * t100 + MDP(19) * t75 + MDP(20) * t57 - MDP(21) * t58 + (t87 * t77 - t97) * t78; (-t86 * t74 - t92) * t75; t64 * t100 + (-MDP(21) * pkin(9) + MDP(18)) * t77 + (-MDP(20) * pkin(9) + t87) * t74; MDP(23) * pkin(5) ^ 2 + MDP(19); t59 * MDP(23); t62 * MDP(23); -t91; t93; 0; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
