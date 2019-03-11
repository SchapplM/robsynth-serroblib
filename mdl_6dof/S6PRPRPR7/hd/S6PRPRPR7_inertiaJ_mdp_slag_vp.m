% Calculate joint inertia matrix for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRPR7_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:53:58
% EndTime: 2019-03-08 19:53:58
% DurationCPUTime: 0.21s
% Computational Cost: add. (204->103), mult. (392->135), div. (0->0), fcn. (335->8), ass. (0->60)
t84 = cos(qJ(4));
t118 = 0.2e1 * t84;
t80 = sin(qJ(6));
t83 = cos(qJ(6));
t92 = t83 * MDP(24) - t80 * MDP(25);
t117 = -2 * MDP(16);
t87 = -pkin(2) - pkin(8);
t116 = pkin(5) - t87;
t115 = (pkin(2) * MDP(7));
t78 = sin(pkin(6));
t82 = sin(qJ(2));
t114 = t78 * t82;
t85 = cos(qJ(2));
t113 = t78 * t85;
t81 = sin(qJ(4));
t112 = t80 * t81;
t111 = t80 * t83;
t110 = t81 * t83;
t74 = t81 ^ 2;
t76 = t84 ^ 2;
t69 = t74 + t76;
t109 = (MDP(18) * pkin(4));
t108 = qJ(3) * MDP(7);
t64 = t81 * pkin(4) - t84 * qJ(5) + qJ(3);
t107 = t64 * MDP(18);
t106 = t87 ^ 2 * MDP(18);
t105 = t80 * MDP(24);
t102 = t83 * MDP(25);
t101 = t87 * MDP(18);
t100 = qJ(5) * MDP(18);
t99 = MDP(13) - MDP(16);
t98 = MDP(14) - MDP(17);
t97 = MDP(20) * t111;
t96 = MDP(5) - t115;
t95 = MDP(16) - t109;
t79 = cos(pkin(6));
t59 = t84 * t113 + t79 * t81;
t60 = -t81 * t113 + t79 * t84;
t54 = t59 * t84 - t60 * t81;
t94 = MDP(13) - t95;
t93 = MDP(21) * t80 + MDP(22) * t83;
t91 = t102 + t105;
t90 = t91 - t98;
t89 = -t92 + t101;
t86 = -pkin(4) - pkin(9);
t88 = (MDP(24) * t86 + MDP(21)) * t83 + (-MDP(25) * t86 - MDP(22)) * t80;
t75 = t83 ^ 2;
t73 = t80 ^ 2;
t72 = t78 ^ 2;
t70 = t72 * t82 ^ 2;
t67 = pkin(4) * t84 + qJ(5) * t81;
t66 = t116 * t84;
t65 = t116 * t81;
t63 = t69 * t87;
t62 = t81 * pkin(9) + t64;
t58 = t83 * t114 + t59 * t80;
t57 = -t80 * t114 + t59 * t83;
t56 = t83 * t62 + t80 * t66;
t55 = -t80 * t62 + t83 * t66;
t1 = [MDP(1) + (t72 * t85 ^ 2 + t79 ^ 2 + t70) * MDP(7) + (t59 ^ 2 + t60 ^ 2 + t70) * MDP(18); (-t60 * t110 + t57 * t84) * MDP(24) + (t60 * t112 - t58 * t84) * MDP(25) + ((MDP(3) - t96) * t85 + (t99 * t81 + t98 * t84 - MDP(4) + MDP(6) + t107 + t108) * t82) * t78 + (MDP(15) - t101) * t54; MDP(2) + (-0.2e1 * t84 * MDP(17) + t107) * t64 + (MDP(23) + MDP(8) + t106) * t76 + (t73 * MDP(19) + t106 + 0.2e1 * t97) * t74 + ((-2 * MDP(5) + t115) * pkin(2)) + (MDP(14) * t118 + 0.2e1 * MDP(6) + t108) * qJ(3) + (0.2e1 * qJ(3) * MDP(13) + t64 * t117 + (-MDP(9) + t93) * t118) * t81 - 0.2e1 * t63 * MDP(15) + 0.2e1 * (t65 * t110 + t55 * t84) * MDP(24) + 0.2e1 * (-t65 * t112 - t56 * t84) * MDP(25); -t54 * MDP(18) - MDP(7) * t113; t63 * MDP(18) + t96 + (-MDP(15) - t92) * t69; t69 * MDP(18) + MDP(7); -t94 * t59 + (t90 + t100) * t60; -t67 * MDP(15) - t91 * t65 + (t94 * t87 + MDP(10) + t88) * t84 + (-MDP(11) + MDP(19) * t111 + (-t73 + t75) * MDP(20) - t98 * t87 + t89 * qJ(5)) * t81; t67 * MDP(18) + t90 * t81 + t99 * t84; -0.2e1 * t97 + t75 * MDP(19) + MDP(12) + (t117 + t109) * pkin(4) + (0.2e1 * MDP(17) + t100 + 0.2e1 * t102 + 0.2e1 * t105) * qJ(5); t59 * MDP(18); (MDP(15) - t89) * t84; -t84 * MDP(18); t95; MDP(18); t57 * MDP(24) - t58 * MDP(25); t84 * MDP(23) + t55 * MDP(24) - t56 * MDP(25) + t93 * t81; -t92 * t84; t88; t92; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
