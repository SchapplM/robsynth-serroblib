% Calculate joint inertia matrix for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRRP6_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:20:57
% EndTime: 2019-03-08 20:20:58
% DurationCPUTime: 0.44s
% Computational Cost: add. (380->136), mult. (745->190), div. (0->0), fcn. (689->8), ass. (0->59)
t131 = pkin(9) * MDP(25);
t140 = MDP(23) + t131;
t122 = MDP(20) + MDP(22);
t121 = MDP(21) - MDP(24);
t95 = sin(qJ(5));
t98 = cos(qJ(5));
t139 = -t121 * t98 - t122 * t95;
t113 = -t98 * pkin(5) - t95 * qJ(6);
t85 = -pkin(4) + t113;
t125 = t85 * MDP(25);
t138 = t121 * t95 - t122 * t98 - MDP(13) + t125;
t137 = 0.2e1 * pkin(5);
t136 = (pkin(2) * MDP(7));
t100 = cos(qJ(2));
t93 = sin(pkin(6));
t130 = t100 * t93;
t94 = cos(pkin(6));
t96 = sin(qJ(4));
t99 = cos(qJ(4));
t78 = t99 * t130 + t94 * t96;
t135 = t78 * t99;
t97 = sin(qJ(2));
t134 = t93 * t97;
t101 = -pkin(2) - pkin(8);
t128 = t101 * t98;
t84 = t96 * pkin(4) - t99 * pkin(9) + qJ(3);
t75 = t96 * t128 + t95 * t84;
t89 = t95 ^ 2;
t91 = t98 ^ 2;
t133 = t89 + t91;
t129 = t101 * t95;
t127 = MDP(23) * t99;
t126 = qJ(3) * MDP(7);
t124 = t98 * MDP(16);
t123 = t99 * MDP(14);
t119 = MDP(5) - t136;
t118 = t133 * MDP(23);
t117 = -MDP(25) * pkin(5) - MDP(22);
t114 = 0.2e1 * qJ(6) * MDP(24) + MDP(19);
t112 = -pkin(5) * t95 + t98 * qJ(6);
t79 = -t96 * t130 + t94 * t99;
t70 = -t98 * t134 + t79 * t95;
t71 = t95 * t134 + t79 * t98;
t111 = t70 * t95 + t71 * t98;
t72 = t96 * qJ(6) + t75;
t81 = t98 * t84;
t73 = -t81 + (-pkin(5) + t129) * t96;
t110 = t72 * t98 + t73 * t95;
t109 = -MDP(20) + t117;
t108 = MDP(25) * qJ(6) - t121;
t107 = t98 * MDP(17) - t95 * MDP(18);
t106 = t95 * MDP(17) + t98 * MDP(18);
t105 = (-t96 * t129 + t81) * MDP(20) - t75 * MDP(21);
t104 = t95 * MDP(22) - t98 * MDP(24);
t103 = t108 * t98 + t109 * t95;
t92 = t99 ^ 2;
t90 = t96 ^ 2;
t77 = (-t101 - t112) * t99;
t1 = [MDP(1) + (t94 ^ 2 + (t100 ^ 2 + t97 ^ 2) * t93 ^ 2) * MDP(7) + (t70 ^ 2 + t71 ^ 2 + t78 ^ 2) * MDP(25); (t70 * t73 + t71 * t72 + t78 * t77) * MDP(25) + (t70 * t98 - t71 * t95) * t127 - t121 * (-t98 * t135 + t71 * t96) + ((MDP(3) - t119) * t100 + (t96 * MDP(13) - MDP(4) + MDP(6) + t123 + t126) * t97) * t93 + t122 * (t95 * t135 - t70 * t96); MDP(2) + t90 * MDP(19) + (t72 ^ 2 + t73 ^ 2 + t77 ^ 2) * MDP(25) + ((-2 * MDP(5) + t136) * pkin(2)) + (0.2e1 * MDP(6) + 0.2e1 * t123 + t126) * qJ(3) + 0.2e1 * ((-t72 * t95 + t73 * t98) * MDP(23) + t104 * t77) * t99 + (t91 * MDP(15) - 0.2e1 * t95 * t124 + MDP(8) + 0.2e1 * (-t95 * MDP(20) - t98 * MDP(21)) * t101) * t92 + 0.2e1 * (qJ(3) * MDP(13) + (-MDP(9) + t107) * t99 - t73 * MDP(22) + t72 * MDP(24) + t105) * t96; -MDP(7) * t130 + (t111 * t96 - t135) * MDP(25); (t110 * t96 - t77 * t99) * MDP(25) + t119 + t139 * (t90 + t92); MDP(7) + (t133 * t90 + t92) * MDP(25); -t79 * MDP(14) + t140 * t111 + t138 * t78; (-t98 * MDP(22) - t95 * MDP(24) + t125) * t77 + (-t101 * MDP(14) + t139 * pkin(9) - MDP(11) + t106) * t96 + (MDP(10) + t101 * MDP(13) + t98 * t95 * MDP(15) + (-t89 + t91) * MDP(16) + (-pkin(4) * t95 + t128) * MDP(20) + (-pkin(4) * t98 - t129) * MDP(21) + t104 * t85) * t99 + t140 * t110; (t133 * t131 - MDP(14) + t118) * t96 - t138 * t99; MDP(12) + t89 * MDP(15) + (t133 * pkin(9) ^ 2 + t85 ^ 2) * MDP(25) + 0.2e1 * pkin(9) * t118 + 0.2e1 * (pkin(4) * MDP(20) - t85 * MDP(22)) * t98 + 0.2e1 * (-pkin(4) * MDP(21) - t85 * MDP(24) + t124) * t95; t108 * t71 + t109 * t70; t81 * MDP(22) + t75 * MDP(24) + (-t73 * pkin(5) + t72 * qJ(6)) * MDP(25) + ((t137 - t129) * MDP(22) + t114) * t96 + (t113 * MDP(23) + t107) * t99 + t105; t103 * t96; t112 * MDP(23) + t103 * pkin(9) + t106; MDP(22) * t137 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(25) + t114; t70 * MDP(25); -t96 * MDP(22) + t73 * MDP(25) + t98 * t127; t95 * t96 * MDP(25); t140 * t95; t117; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
