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
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S6PPRRRP1_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:54:26
% EndTime: 2019-03-08 18:54:26
% DurationCPUTime: 0.31s
% Computational Cost: add. (426->130), mult. (1067->207), div. (0->0), fcn. (1217->12), ass. (0->61)
t115 = 2 * MDP(20);
t86 = sin(qJ(5));
t114 = pkin(9) * t86;
t89 = cos(qJ(5));
t113 = pkin(9) * t89;
t90 = cos(qJ(4));
t112 = pkin(9) * t90;
t81 = sin(pkin(7));
t88 = sin(qJ(3));
t111 = t81 * t88;
t91 = cos(qJ(3));
t110 = t81 * t91;
t83 = cos(pkin(12));
t84 = cos(pkin(7));
t109 = t83 * t84;
t108 = t86 * t89;
t107 = -qJ(6) - pkin(10);
t106 = MDP(21) * pkin(5);
t87 = sin(qJ(4));
t105 = qJ(6) * t87;
t104 = MDP(17) * t90;
t80 = sin(pkin(12));
t82 = sin(pkin(6));
t85 = cos(pkin(6));
t61 = t85 * t111 + (t88 * t109 + t80 * t91) * t82;
t67 = -t81 * t82 * t83 + t84 * t85;
t58 = t61 * t90 + t67 * t87;
t60 = -t85 * t110 + (-t91 * t109 + t80 * t88) * t82;
t56 = t58 * t89 + t60 * t86;
t103 = t56 * MDP(19);
t69 = t90 * t111 + t84 * t87;
t63 = -t86 * t110 + t89 * t69;
t102 = t63 * MDP(19);
t75 = -pkin(5) * t89 - pkin(4);
t101 = t75 * MDP(21);
t100 = t86 * MDP(16);
t99 = t89 * t112;
t98 = MDP(14) * t108;
t97 = -MDP(20) * pkin(5) + MDP(15);
t96 = MDP(18) + t106;
t95 = t89 * MDP(18) - t86 * MDP(19);
t94 = t86 * MDP(18) + t89 * MDP(19);
t93 = t90 * MDP(11) - t87 * MDP(12) + MDP(4);
t92 = -MDP(11) - t95 + t101;
t79 = t89 ^ 2;
t78 = t87 ^ 2;
t77 = t86 ^ 2;
t74 = t107 * t89;
t73 = t107 * t86;
t72 = -pkin(4) * t90 - pkin(10) * t87 - pkin(3);
t71 = (pkin(5) * t86 + pkin(9)) * t87;
t70 = t89 * t72;
t68 = t87 * t111 - t84 * t90;
t66 = t72 * t86 + t99;
t65 = -t86 * t112 + t70;
t64 = t99 + (t72 - t105) * t86;
t62 = -t89 * t110 - t86 * t69;
t59 = -t89 * t105 + t70 + (-pkin(5) - t114) * t90;
t57 = t61 * t87 - t67 * t90;
t55 = -t58 * t86 + t60 * t89;
t1 = [MDP(1) + (t85 ^ 2 + (t80 ^ 2 + t83 ^ 2) * t82 ^ 2) * MDP(2) + (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) * MDP(21); t85 * MDP(2) + (t55 * t62 + t56 * t63 + t57 * t68) * MDP(21); MDP(2) + (t62 ^ 2 + t63 ^ 2 + t68 ^ 2) * MDP(21); -t61 * MDP(5) + (t55 * t59 + t56 * t64 + t57 * t71) * MDP(21) + (-t55 * MDP(18) + t103) * t90 + ((-t55 * t89 - t56 * t86) * MDP(20) + t94 * t57) * t87 - t93 * t60; (t59 * t62 + t63 * t64 + t68 * t71) * MDP(21) + (-t62 * MDP(18) + t102) * t90 + ((-t62 * t89 - t63 * t86) * MDP(20) + t94 * t68) * t87 + (-t88 * MDP(5) + t93 * t91) * t81; MDP(3) + (t59 ^ 2 + t64 ^ 2 + t71 ^ 2) * MDP(21) + (0.2e1 * pkin(3) * MDP(11) + t104) * t90 + (t79 * MDP(13) + MDP(6) - 0.2e1 * t98) * t78 + 0.2e1 * (t78 * t114 - t65 * t90) * MDP(18) + 0.2e1 * (t78 * t113 + t66 * t90) * MDP(19) + (-0.2e1 * pkin(3) * MDP(12) + (-t59 * t89 - t64 * t86) * t115 + 0.2e1 * (-t89 * MDP(15) + MDP(7) + t100) * t90) * t87; -t58 * MDP(12) + (-t55 * t86 + t56 * t89) * MDP(20) + (t55 * t73 - t56 * t74) * MDP(21) + t92 * t57; -t69 * MDP(12) + (-t62 * t86 + t63 * t89) * MDP(20) + (t62 * t73 - t63 * t74) * MDP(21) + t92 * t68; (-t59 * t86 + t64 * t89) * MDP(20) + (t59 * t73 - t64 * t74 + t71 * t75) * MDP(21) + (-pkin(9) * MDP(12) - t86 * MDP(15) - t89 * MDP(16) + t94 * pkin(10) + MDP(9)) * t90 + (MDP(8) - pkin(9) * MDP(11) + MDP(13) * t108 + (-t77 + t79) * MDP(14) + (-pkin(4) * t86 - t113) * MDP(18) + (-pkin(4) * t89 + t114) * MDP(19) + (-t73 * t89 + t74 * t86) * MDP(20)) * t87; MDP(10) + t77 * MDP(13) + 0.2e1 * t98 + (-t73 * t86 - t74 * t89) * t115 + (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) * MDP(21) + 0.2e1 * t95 * pkin(4); t96 * t55 - t103; t96 * t62 - t102; t59 * t106 - t104 + t65 * MDP(18) - t66 * MDP(19) + (t97 * t89 - t100) * t87; t73 * t106 + (-MDP(19) * pkin(10) + MDP(16)) * t89 + (-MDP(18) * pkin(10) + t97) * t86; MDP(21) * pkin(5) ^ 2 + MDP(17); t57 * MDP(21); t68 * MDP(21); t71 * MDP(21); t101; 0; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
