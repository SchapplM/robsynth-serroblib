% Calculate joint inertia matrix for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPPRRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_inertiaJ_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S6PPPRRR1_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:40:41
% EndTime: 2019-03-08 18:40:42
% DurationCPUTime: 0.32s
% Computational Cost: add. (360->99), mult. (1026->181), div. (0->0), fcn. (1256->16), ass. (0->65)
t101 = sin(qJ(6));
t104 = cos(qJ(6));
t129 = -(MDP(19) * t101 + MDP(20) * t104) * pkin(11) + t101 * MDP(16) + t104 * MDP(17);
t128 = pkin(10) * t101;
t127 = pkin(10) * t104;
t105 = cos(qJ(5));
t126 = pkin(10) * t105;
t91 = sin(pkin(14));
t94 = sin(pkin(7));
t125 = t91 * t94;
t96 = cos(pkin(14));
t124 = t94 * t96;
t97 = cos(pkin(13));
t99 = cos(pkin(7));
t123 = t97 * t99;
t100 = cos(pkin(6));
t122 = t100 * t94;
t103 = sin(qJ(4));
t93 = sin(pkin(8));
t121 = t103 * t93;
t106 = cos(qJ(4));
t120 = t106 * t93;
t102 = sin(qJ(5));
t119 = t101 * t102;
t118 = t101 * t104;
t117 = t102 * t104;
t116 = t105 * MDP(18);
t115 = MDP(15) * t118;
t92 = sin(pkin(13));
t95 = sin(pkin(6));
t73 = t96 * t122 + (t96 * t123 - t91 * t92) * t95;
t82 = -t95 * t97 * t94 + t100 * t99;
t98 = cos(pkin(8));
t114 = t73 * t98 + t82 * t93;
t113 = MDP(16) * t104 - MDP(17) * t101;
t111 = t104 * MDP(19) - t101 * MDP(20);
t109 = t98 * t124 + t93 * t99;
t108 = t105 * MDP(12) - t102 * MDP(13) + MDP(5);
t107 = -MDP(12) - t111;
t90 = t104 ^ 2;
t89 = t102 ^ 2;
t88 = t101 ^ 2;
t85 = -t105 * pkin(5) - t102 * pkin(11) - pkin(4);
t84 = t102 * t98 + t105 * t121;
t83 = t102 * t121 - t105 * t98;
t81 = -t93 * t124 + t98 * t99;
t80 = t101 * t85 + t104 * t126;
t79 = -t101 * t126 + t104 * t85;
t78 = -t101 * t120 + t104 * t84;
t77 = -t101 * t84 - t104 * t120;
t76 = t109 * t103 + t106 * t125;
t75 = t103 * t125 - t109 * t106;
t74 = t95 * t92 * t96 + (t95 * t123 + t122) * t91;
t72 = t102 * t81 + t105 * t76;
t71 = t102 * t76 - t105 * t81;
t70 = -t73 * t93 + t82 * t98;
t69 = t101 * t75 + t104 * t72;
t68 = -t101 * t72 + t104 * t75;
t67 = t114 * t103 + t74 * t106;
t66 = t74 * t103 - t114 * t106;
t65 = t70 * t102 + t67 * t105;
t64 = t67 * t102 - t70 * t105;
t63 = t66 * t101 + t65 * t104;
t62 = -t65 * t101 + t66 * t104;
t1 = [MDP(1) + (t100 ^ 2 + (t92 ^ 2 + t97 ^ 2) * t95 ^ 2) * MDP(2) + (t73 ^ 2 + t74 ^ 2 + t82 ^ 2) * MDP(3); t100 * MDP(2) + (t82 * t99 + (t73 * t96 + t74 * t91) * t94) * MDP(3); MDP(2) + (t99 ^ 2 + (t91 ^ 2 + t96 ^ 2) * t94 ^ 2) * MDP(3); t82 * MDP(3); t99 * MDP(3); MDP(3); -t67 * MDP(6) + (-t62 * t105 + t64 * t119) * MDP(19) + (t63 * t105 + t64 * t117) * MDP(20) - t108 * t66; -t76 * MDP(6) + (-t68 * t105 + t71 * t119) * MDP(19) + (t69 * t105 + t71 * t117) * MDP(20) - t108 * t75; (-t77 * t105 + t83 * t119) * MDP(19) + (t78 * t105 + t83 * t117) * MDP(20) + (-t103 * MDP(6) + t108 * t106) * t93; MDP(4) + (0.2e1 * pkin(4) * MDP(12) + t116) * t105 + (t90 * MDP(14) + MDP(7) - 0.2e1 * t115) * t89 + 0.2e1 * (-t79 * t105 + t89 * t128) * MDP(19) + 0.2e1 * (t80 * t105 + t89 * t127) * MDP(20) + 0.2e1 * (-pkin(4) * MDP(13) + (MDP(8) - t113) * t105) * t102; -t65 * MDP(13) + t107 * t64; -t72 * MDP(13) + t107 * t71; -t84 * MDP(13) + t107 * t83; (-pkin(10) * MDP(13) + MDP(10) - t129) * t105 + (MDP(9) - pkin(10) * MDP(12) + MDP(14) * t118 + (-t88 + t90) * MDP(15) + (-pkin(5) * t101 - t127) * MDP(19) + (-pkin(5) * t104 + t128) * MDP(20)) * t102; t88 * MDP(14) + 0.2e1 * pkin(5) * t111 + MDP(11) + 0.2e1 * t115; t62 * MDP(19) - t63 * MDP(20); t68 * MDP(19) - t69 * MDP(20); t77 * MDP(19) - t78 * MDP(20); t79 * MDP(19) - t80 * MDP(20) + t113 * t102 - t116; t129; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
