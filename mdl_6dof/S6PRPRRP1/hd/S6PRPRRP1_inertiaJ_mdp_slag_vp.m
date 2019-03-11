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
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S6PRPRRP1_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:58:30
% EndTime: 2019-03-08 19:58:31
% DurationCPUTime: 0.29s
% Computational Cost: add. (348->120), mult. (758->188), div. (0->0), fcn. (777->10), ass. (0->56)
t114 = 2 * MDP(20);
t82 = sin(pkin(11));
t74 = pkin(2) * t82 + pkin(8);
t86 = sin(qJ(5));
t113 = t74 * t86;
t89 = cos(qJ(5));
t112 = t74 * t89;
t90 = cos(qJ(4));
t111 = t74 * t90;
t110 = t86 * t89;
t109 = -qJ(6) - pkin(9);
t108 = MDP(21) * pkin(5);
t87 = sin(qJ(4));
t107 = qJ(6) * t87;
t83 = sin(pkin(6));
t84 = cos(pkin(11));
t88 = sin(qJ(2));
t91 = cos(qJ(2));
t67 = (t82 * t91 + t84 * t88) * t83;
t85 = cos(pkin(6));
t64 = t67 * t90 + t85 * t87;
t65 = (t82 * t88 - t84 * t91) * t83;
t58 = t64 * t89 + t65 * t86;
t106 = MDP(19) * t58;
t105 = MDP(19) * t89;
t77 = -pkin(5) * t89 - pkin(4);
t104 = t77 * MDP(21);
t103 = t86 * MDP(16);
t102 = t89 * t111;
t101 = MDP(14) * t110;
t75 = -pkin(2) * t84 - pkin(3);
t100 = -MDP(20) * pkin(5) + MDP(15);
t99 = MDP(18) + t108;
t57 = -t64 * t86 + t65 * t89;
t98 = -t57 * t86 + t58 * t89;
t70 = -pkin(4) * t90 - pkin(9) * t87 + t75;
t68 = t89 * t70;
t59 = -t89 * t107 + t68 + (-pkin(5) - t113) * t90;
t60 = t102 + (t70 - t107) * t86;
t97 = -t59 * t86 + t60 * t89;
t71 = t109 * t86;
t72 = t109 * t89;
t96 = -t71 * t86 - t72 * t89;
t95 = MDP(18) * t89 - MDP(19) * t86;
t94 = MDP(18) * t86 + t105;
t93 = MDP(11) + t95 - t104;
t81 = t90 ^ 2;
t80 = t89 ^ 2;
t79 = t87 ^ 2;
t78 = t86 ^ 2;
t76 = t80 * t87;
t69 = (pkin(5) * t86 + t74) * t87;
t63 = t67 * t87 - t85 * t90;
t62 = t70 * t86 + t102;
t61 = -t86 * t111 + t68;
t1 = [MDP(1) + (t65 ^ 2 + t67 ^ 2 + t85 ^ 2) * MDP(5) + (t57 ^ 2 + t58 ^ 2 + t63 ^ 2) * MDP(21); (t57 * t59 + t58 * t60 + t63 * t69) * MDP(21) + (MDP(3) * t91 - MDP(4) * t88) * t83 + (-t65 * t84 + t67 * t82) * MDP(5) * pkin(2) + (-MDP(11) * t65 - MDP(18) * t57 + t106) * t90 + (t65 * MDP(12) + (-t57 * t89 - t58 * t86) * MDP(20) + t94 * t63) * t87; MDP(2) - 0.2e1 * t75 * t90 * MDP(11) + t81 * MDP(17) + (t59 ^ 2 + t60 ^ 2 + t69 ^ 2) * MDP(21) + (t82 ^ 2 + t84 ^ 2) * MDP(5) * pkin(2) ^ 2 + (t80 * MDP(13) + MDP(6) - 0.2e1 * t101) * t79 + 0.2e1 * (t79 * t113 - t61 * t90) * MDP(18) + 0.2e1 * (t79 * t112 + t62 * t90) * MDP(19) + (0.2e1 * t75 * MDP(12) + (-t59 * t89 - t60 * t86) * t114 + 0.2e1 * (-t89 * MDP(15) + MDP(7) + t103) * t90) * t87; t85 * MDP(5) + (-t63 * t90 + t98 * t87) * MDP(21); (-t69 * t90 + t97 * t87) * MDP(21); MDP(5) + (t81 + (t78 + t80) * t79) * MDP(21); -t64 * MDP(12) + t98 * MDP(20) + (t57 * t71 - t58 * t72) * MDP(21) - t93 * t63; t76 * MDP(14) + t97 * MDP(20) + (t59 * t71 - t60 * t72 + t69 * t77) * MDP(21) + (-t74 * MDP(12) - t86 * MDP(15) - t89 * MDP(16) + t94 * pkin(9) + MDP(9)) * t90 + (MDP(8) - t74 * MDP(11) + MDP(13) * t110 - t78 * MDP(14) + (-pkin(4) * t86 - t112) * MDP(18) + (-pkin(4) * t89 + t113) * MDP(19) + (-t71 * t89 + t72 * t86) * MDP(20)) * t87; t76 * MDP(20) + (t78 * MDP(20) + t96 * MDP(21) - MDP(12)) * t87 + t93 * t90; MDP(10) + t78 * MDP(13) + 0.2e1 * t101 + t96 * t114 + (t71 ^ 2 + t72 ^ 2 + t77 ^ 2) * MDP(21) + 0.2e1 * t95 * pkin(4); t99 * t57 - t106; t59 * t108 - MDP(17) * t90 + t61 * MDP(18) - t62 * MDP(19) + (t100 * t89 - t103) * t87; (-t99 * t86 - t105) * t87; t71 * t108 + (-MDP(19) * pkin(9) + MDP(16)) * t89 + (-MDP(18) * pkin(9) + t100) * t86; MDP(21) * pkin(5) ^ 2 + MDP(17); t63 * MDP(21); t69 * MDP(21); -t90 * MDP(21); t104; 0; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
