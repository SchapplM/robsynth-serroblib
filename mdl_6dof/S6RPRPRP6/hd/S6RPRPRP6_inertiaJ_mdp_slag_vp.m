% Calculate joint inertia matrix for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRPRP6_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:19:46
% EndTime: 2019-03-09 03:19:48
% DurationCPUTime: 0.34s
% Computational Cost: add. (559->126), mult. (976->164), div. (0->0), fcn. (1018->6), ass. (0->58)
t134 = -2 * MDP(16);
t133 = pkin(3) + pkin(8);
t132 = pkin(1) * MDP(7);
t100 = sin(qJ(3));
t102 = cos(qJ(3));
t130 = pkin(7) + qJ(2);
t97 = sin(pkin(9));
t86 = t130 * t97;
t98 = cos(pkin(9));
t87 = t130 * t98;
t76 = t100 * t87 + t102 * t86;
t83 = t100 * t98 + t102 * t97;
t72 = t83 * pkin(4) + t76;
t99 = sin(qJ(5));
t131 = t99 * t72;
t128 = MDP(27) * pkin(5);
t127 = (pkin(3) * MDP(18));
t82 = t100 * t97 - t102 * t98;
t126 = qJ(4) * t82;
t125 = t97 * MDP(5);
t124 = t98 * MDP(4);
t123 = t99 * MDP(24);
t122 = t99 * MDP(25);
t95 = t99 ^ 2;
t101 = cos(qJ(5));
t96 = t101 ^ 2;
t89 = -t95 - t96;
t121 = -t89 * MDP(27) + MDP(18);
t120 = -qJ(6) - t133;
t119 = MDP(18) * qJ(4);
t118 = MDP(22) * t101;
t117 = t101 * MDP(25);
t116 = MDP(17) - MDP(14);
t92 = -t98 * pkin(2) - pkin(1);
t114 = t101 * t99 * MDP(20);
t105 = -t83 * qJ(4) + t92;
t70 = t133 * t82 + t105;
t113 = qJ(6) * t82 + t70;
t112 = MDP(16) - t127;
t111 = -MDP(26) * pkin(5) + MDP(21);
t110 = MDP(24) + t128;
t71 = t101 * t72;
t65 = t83 * pkin(5) - t113 * t99 + t71;
t66 = t113 * t101 + t131;
t109 = t66 * t101 - t65 * t99;
t108 = t133 * t83 + t126;
t77 = -t100 * t86 + t102 * t87;
t107 = (-t99 * t70 + t71) * MDP(24) - (t101 * t70 + t131) * MDP(25);
t106 = -t101 * MDP(24) + t122;
t91 = t99 * pkin(5) + qJ(4);
t85 = t120 * t101;
t84 = t120 * t99;
t79 = t96 * t82;
t75 = t85 * t101 + t84 * t99;
t74 = t82 * pkin(3) + t105;
t73 = -t82 * pkin(4) + t77;
t69 = (-pkin(5) * t101 - pkin(4)) * t82 + t77;
t1 = [MDP(1) + (t74 ^ 2 + t76 ^ 2 + t77 ^ 2) * MDP(18) + (t65 ^ 2 + t66 ^ 2 + t69 ^ 2) * MDP(27) + (MDP(8) + MDP(23)) * t83 ^ 2 + (t95 * MDP(19) + 0.2e1 * t114) * t82 ^ 2 + (0.2e1 * t124 - 0.2e1 * t125 + t132) * pkin(1) + (0.2e1 * t92 * MDP(13) + t74 * t134) * t82 + 0.2e1 * (-t77 * MDP(15) + t109 * MDP(26) + t106 * t73) * t82 + (MDP(7) * qJ(2) + 2 * MDP(6)) * (t97 ^ 2 + t98 ^ 2) * qJ(2) + 0.2e1 * (t92 * MDP(14) - t74 * MDP(17) + (MDP(21) * t99 - MDP(9) + t118) * t82 + t76 * MDP(15) + t107) * t83; -t124 + t125 - t132 + t74 * MDP(18) + t79 * MDP(26) + t109 * MDP(27) + (t95 * MDP(26) + MDP(13) - MDP(16)) * t82 + (-t116 - t117 - t123) * t83; MDP(7) + t121; t83 * MDP(10) - t82 * MDP(11) + (-pkin(3) * t83 - t126) * MDP(15) + (-t95 * t82 + t79) * MDP(20) + (t65 * t85 + t66 * t84 + t69 * t91) * MDP(27) + (t116 + t119) * t77 + (-MDP(13) + t112) * t76 + (-t83 * MDP(22) + t73 * MDP(24) + t108 * MDP(25) + (-t82 * t85 - t66) * MDP(26)) * t99 + (t99 * t82 * MDP(19) + t83 * MDP(21) - t108 * MDP(24) + t73 * MDP(25) + (t82 * t84 - t65) * MDP(26)) * t101; (t101 * t84 - t99 * t85) * MDP(27); MDP(12) + t96 * MDP(19) - 0.2e1 * t114 - 0.2e1 * t75 * MDP(26) + (t84 ^ 2 + t85 ^ 2 + t91 ^ 2) * MDP(27) + (t134 + t127) * pkin(3) + (0.2e1 * MDP(17) + 0.2e1 * t117 + t119 + 0.2e1 * t123) * qJ(4); t76 * MDP(18) + (t65 * t101 + t66 * t99) * MDP(27) + (MDP(15) - t106) * t83; 0; t89 * MDP(26) + t75 * MDP(27) + t112; t121; t65 * t128 + t83 * MDP(23) + (t111 * t99 + t118) * t82 + t107; -t110 * t99 - t117; t85 * t128 + (MDP(25) * t133 - MDP(22)) * t99 + (-MDP(24) * t133 + t111) * t101; t110 * t101 - t122; MDP(27) * pkin(5) ^ 2 + MDP(23); t69 * MDP(27); 0; t91 * MDP(27); 0; 0; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
