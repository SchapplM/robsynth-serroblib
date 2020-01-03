% Calculate joint inertia matrix for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRRR8_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:06
% EndTime: 2019-12-31 19:06:07
% DurationCPUTime: 0.14s
% Computational Cost: add. (209->67), mult. (328->85), div. (0->0), fcn. (301->6), ass. (0->45)
t78 = sin(qJ(5));
t79 = sin(qJ(4));
t81 = cos(qJ(5));
t82 = cos(qJ(4));
t63 = t78 * t79 - t81 * t82;
t64 = t78 * t82 + t81 * t79;
t113 = -t64 * MDP(19) + t63 * MDP(20);
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t95 = t82 * MDP(15);
t90 = -t79 * MDP(16) + t95;
t98 = t63 * MDP(22);
t112 = (-t64 * MDP(23) + MDP(8) + t90 - t98) * t83 - t80 * MDP(9);
t108 = t82 * pkin(4);
t84 = -pkin(1) - pkin(2);
t68 = t80 * qJ(2) - t83 * t84;
t66 = pkin(3) + t68;
t58 = t66 + t108;
t111 = -0.2e1 * t58;
t74 = -pkin(3) - t108;
t110 = 0.2e1 * t74;
t109 = pkin(7) + pkin(8);
t107 = pkin(3) + t66;
t69 = t83 * qJ(2) + t80 * t84;
t67 = -pkin(7) + t69;
t106 = pkin(8) - t67;
t105 = (-MDP(22) * t64 + MDP(23) * t63) * t80;
t104 = t58 - t74;
t103 = t68 * MDP(8);
t102 = t69 * MDP(9);
t100 = t79 ^ 2 * MDP(10) + MDP(7);
t96 = t82 * MDP(11);
t94 = 0.2e1 * t79 * t96 + t100;
t70 = t109 * t79;
t71 = t109 * t82;
t93 = (-t81 * t70 - t78 * t71) * MDP(22) + (t78 * t70 - t81 * t71) * MDP(23) - t113;
t54 = t106 * t79;
t55 = t106 * t82;
t92 = (t81 * t54 + t78 * t55) * MDP(22) + (-t78 * t54 + t81 * t55) * MDP(23) + t113;
t91 = -t79 * MDP(12) - t82 * MDP(13);
t89 = -MDP(15) * t79 - MDP(16) * t82;
t88 = MDP(17) * t64 - 0.2e1 * t63 * MDP(18);
t87 = (MDP(22) * t81 - MDP(23) * t78) * pkin(4);
t86 = 0.2e1 * t90;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + t98 * t111 + t66 * t86 + 0.2e1 * t103 + 0.2e1 * t102 + (MDP(23) * t111 + t88) * t64 + t94; -pkin(1) * MDP(6) - MDP(4) - t112; MDP(6); -t103 - t102 - t107 * t95 + t104 * t98 + (t107 * MDP(16) - 0.2e1 * t96) * t79 + (t104 * MDP(23) - t88) * t64 - t100; t112; t98 * t110 + pkin(3) * t86 + (MDP(23) * t110 + t88) * t64 + t94; t89 * t67 + t91 + t92; t89 * t80 + t105; t89 * pkin(7) - t91 + t93; MDP(14) + MDP(21) + 0.2e1 * t87; t92; t105; t93; MDP(21) + t87; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
