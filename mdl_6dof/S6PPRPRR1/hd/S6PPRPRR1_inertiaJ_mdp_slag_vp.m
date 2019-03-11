% Calculate joint inertia matrix for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRPRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S6PPRPRR1_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:43:39
% EndTime: 2019-03-08 18:43:40
% DurationCPUTime: 0.28s
% Computational Cost: add. (287->95), mult. (781->166), div. (0->0), fcn. (916->14), ass. (0->57)
t92 = sin(qJ(6));
t95 = cos(qJ(6));
t101 = MDP(19) * t92 + MDP(20) * t95;
t119 = t92 * MDP(16) + t95 * MDP(17) - t101 * pkin(10);
t118 = MDP(6) * pkin(3);
t84 = sin(pkin(13));
t78 = pkin(3) * t84 + pkin(9);
t117 = t78 * t92;
t116 = t78 * t95;
t96 = cos(qJ(5));
t115 = t78 * t96;
t86 = sin(pkin(7));
t94 = sin(qJ(3));
t114 = t86 * t94;
t97 = cos(qJ(3));
t113 = t86 * t97;
t89 = cos(pkin(12));
t90 = cos(pkin(7));
t112 = t89 * t90;
t93 = sin(qJ(5));
t111 = t92 * t93;
t110 = t92 * t95;
t109 = t93 * t95;
t108 = t93 * MDP(13);
t107 = t96 * MDP(18);
t106 = MDP(15) * t110;
t88 = cos(pkin(13));
t79 = -pkin(3) * t88 - pkin(4);
t105 = -MDP(12) * t96 + t108;
t104 = MDP(16) * t95 - MDP(17) * t92;
t102 = MDP(19) * t95 - MDP(20) * t92;
t100 = MDP(12) + t102;
t85 = sin(pkin(12));
t87 = sin(pkin(6));
t91 = cos(pkin(6));
t99 = t91 * t113 + (t97 * t112 - t85 * t94) * t87;
t83 = t95 ^ 2;
t82 = t93 ^ 2;
t81 = t92 ^ 2;
t76 = -pkin(5) * t96 - pkin(10) * t93 + t79;
t75 = -t86 * t87 * t89 + t90 * t91;
t74 = (t84 * t97 + t88 * t94) * t86;
t72 = -t88 * t113 + t84 * t114;
t71 = t74 * t96 + t90 * t93;
t70 = t74 * t93 - t90 * t96;
t69 = t95 * t115 + t76 * t92;
t68 = -t92 * t115 + t76 * t95;
t67 = t91 * t114 + (t94 * t112 + t85 * t97) * t87;
t65 = t71 * t95 + t72 * t92;
t64 = -t71 * t92 + t72 * t95;
t63 = t88 * t67 + t84 * t99;
t61 = t67 * t84 - t88 * t99;
t60 = t63 * t96 + t75 * t93;
t59 = t63 * t93 - t75 * t96;
t58 = t60 * t95 + t61 * t92;
t57 = -t60 * t92 + t61 * t95;
t1 = [MDP(1) + (t91 ^ 2 + (t85 ^ 2 + t89 ^ 2) * t87 ^ 2) * MDP(2) + (t61 ^ 2 + t63 ^ 2 + t75 ^ 2) * MDP(6); t91 * MDP(2) + (t61 * t72 + t63 * t74 + t75 * t90) * MDP(6); MDP(2) + (t72 ^ 2 + t74 ^ 2 + t90 ^ 2) * MDP(6); t99 * MDP(4) - t67 * MDP(5) + (t59 * t111 - t57 * t96) * MDP(19) + (t59 * t109 + t58 * t96) * MDP(20) + t105 * t61 + (-t61 * t88 + t63 * t84) * t118; (t70 * t111 - t64 * t96) * MDP(19) + (t70 * t109 + t65 * t96) * MDP(20) + (MDP(4) * t97 - MDP(5) * t94) * t86 + t105 * t72 + (-t72 * t88 + t74 * t84) * t118; MDP(3) + (t84 ^ 2 + t88 ^ 2) * MDP(6) * pkin(3) ^ 2 + (-0.2e1 * MDP(12) * t79 + t107) * t96 + (MDP(14) * t83 + MDP(7) - 0.2e1 * t106) * t82 + 0.2e1 * (t82 * t117 - t68 * t96) * MDP(19) + 0.2e1 * (t82 * t116 + t69 * t96) * MDP(20) + 0.2e1 * (t79 * MDP(13) + (MDP(8) - t104) * t96) * t93; t75 * MDP(6); t90 * MDP(6); 0; MDP(6); -t60 * MDP(13) - t100 * t59; -t71 * MDP(13) - t100 * t70; (-t78 * MDP(13) + MDP(10) - t119) * t96 + (MDP(9) - t78 * MDP(12) + MDP(14) * t110 + (-t81 + t83) * MDP(15) + (-pkin(5) * t92 - t116) * MDP(19) + (-pkin(5) * t95 + t117) * MDP(20)) * t93; t100 * t96 - t108; t81 * MDP(14) + 0.2e1 * pkin(5) * t102 + MDP(11) + 0.2e1 * t106; MDP(19) * t57 - MDP(20) * t58; MDP(19) * t64 - MDP(20) * t65; t68 * MDP(19) - t69 * MDP(20) + t104 * t93 - t107; -t101 * t93; t119; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
