% Calculate joint inertia matrix for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPRP2_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:06:14
% EndTime: 2019-03-09 03:06:16
% DurationCPUTime: 0.42s
% Computational Cost: add. (724->136), mult. (1244->200), div. (0->0), fcn. (1284->8), ass. (0->61)
t103 = sin(qJ(5));
t105 = cos(qJ(5));
t125 = MDP(20) - MDP(23);
t141 = (MDP(19) + MDP(21)) * t105 - t125 * t103;
t101 = cos(pkin(10));
t104 = sin(qJ(3));
t138 = cos(qJ(3));
t99 = sin(pkin(10));
t87 = -t101 * t138 + t99 * t104;
t85 = t87 ^ 2;
t139 = t87 * pkin(5);
t93 = t99 * pkin(3) + pkin(8);
t137 = t87 * t93;
t89 = t101 * t104 + t99 * t138;
t97 = t103 ^ 2;
t136 = t97 * t89;
t102 = cos(pkin(9));
t96 = -t102 * pkin(1) - pkin(2);
t90 = -t138 * pkin(3) + t96;
t76 = t87 * pkin(4) - t89 * pkin(8) + t90;
t100 = sin(pkin(9));
t94 = t100 * pkin(1) + pkin(7);
t120 = (-qJ(4) - t94) * t104;
t123 = t138 * t94;
t84 = t138 * qJ(4) + t123;
t79 = t101 * t84 + t99 * t120;
t71 = t103 * t76 + t105 * t79;
t98 = t105 ^ 2;
t135 = t97 + t98;
t134 = MDP(13) * pkin(3);
t133 = t105 * t89;
t132 = t87 * qJ(6);
t82 = t98 * t89;
t131 = (-t82 - t136) * MDP(22);
t113 = -t105 * pkin(5) - t103 * qJ(6);
t95 = -t101 * pkin(3) - pkin(4);
t83 = t113 + t95;
t130 = t83 * MDP(24);
t129 = t87 * MDP(17);
t128 = t105 * MDP(15);
t127 = t105 * MDP(16);
t124 = 0.2e1 * t103;
t77 = -t101 * t120 + t99 * t84;
t122 = t138 * MDP(10);
t121 = t103 * t79 - t105 * t76;
t119 = t135 * MDP(24);
t118 = -MDP(24) * pkin(5) - MDP(21);
t117 = t93 * MDP(24) + MDP(22);
t116 = t83 * t89 - t137;
t115 = t89 * t95 - t137;
t114 = MDP(19) - t118;
t112 = pkin(5) * t103 - t105 * qJ(6);
t68 = t132 + t71;
t69 = t121 - t139;
t111 = t68 * t103 - t69 * t105;
t110 = MDP(24) * qJ(6) - t125;
t109 = -t121 * MDP(19) - t71 * MDP(20);
t108 = -t114 * t103 + t110 * t105;
t86 = t89 ^ 2;
t72 = t112 * t89 + t77;
t1 = [MDP(1) - 0.2e1 * t96 * t122 + (t77 ^ 2 + t79 ^ 2 + t90 ^ 2) * MDP(13) + t98 * t86 * MDP(14) + t85 * MDP(18) + (t68 ^ 2 + t69 ^ 2 + t72 ^ 2) * MDP(24) + (t100 ^ 2 + t102 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-t86 * t128 - t89 * t129) * t124 + (0.2e1 * t96 * MDP(11) + MDP(5) * t104 + 0.2e1 * t138 * MDP(6)) * t104 + 0.2e1 * (-t79 * MDP(12) - t69 * MDP(21) + t68 * MDP(23) + t109) * t87 + 0.2e1 * (t87 * t127 - t111 * MDP(22) + (t103 * MDP(21) - t105 * MDP(23)) * t72 + (t103 * MDP(19) + t105 * MDP(20) + MDP(12)) * t77) * t89; (t77 * t87 + t79 * t89) * MDP(13) + (t72 * t87 + (t69 * t103 + t68 * t105) * t89) * MDP(24); MDP(4) + (t86 + t85) * MDP(13) + (t135 * t86 + t85) * MDP(24); t138 * MDP(8) - MDP(11) * t123 + (t82 - t136) * MDP(15) + t72 * t130 + (-t94 * MDP(10) + MDP(7)) * t104 + (-t77 * MDP(19) + t115 * MDP(20) - t72 * MDP(21) - t116 * MDP(23) + t117 * t68 + t129) * t105 + (MDP(14) * t133 + t87 * MDP(16) + t115 * MDP(19) + t77 * MDP(20) + t116 * MDP(21) - t72 * MDP(23) + t117 * t69) * t103 + ((-t101 * t89 - t87 * t99) * MDP(12) + (-t101 * t77 + t79 * t99) * MDP(13)) * pkin(3); t122 - t104 * MDP(11) - t131 + (t93 * t119 + t99 * t134) * t89 + (-t101 * t134 + t130 - t141) * t87; MDP(9) + t97 * MDP(14) + (t135 * t93 ^ 2 + t83 ^ 2) * MDP(24) + 0.2e1 * t135 * MDP(22) * t93 + (t101 ^ 2 + t99 ^ 2) * MDP(13) * pkin(3) ^ 2 + 0.2e1 * (-t95 * MDP(19) - t83 * MDP(21)) * t105 + (t95 * MDP(20) - t83 * MDP(23) + t128) * t124; t90 * MDP(13) + t111 * MDP(24) + t141 * t87 + t131; 0; 0; MDP(13) + t119; t87 * MDP(18) + (-t121 + 0.2e1 * t139) * MDP(21) + (0.2e1 * t132 + t71) * MDP(23) + (-t69 * pkin(5) + t68 * qJ(6)) * MDP(24) + (-t103 * MDP(17) + t113 * MDP(22) + t127) * t89 + t109; t108 * t89; t103 * MDP(16) + t105 * MDP(17) - t112 * MDP(22) + t108 * t93; t110 * t103 + t114 * t105; MDP(18) + 0.2e1 * pkin(5) * MDP(21) + 0.2e1 * qJ(6) * MDP(23) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(24); -t87 * MDP(21) + MDP(22) * t133 + t69 * MDP(24); t103 * t89 * MDP(24); t117 * t103; -t105 * MDP(24); t118; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
