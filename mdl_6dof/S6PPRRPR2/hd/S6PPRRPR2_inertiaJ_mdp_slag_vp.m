% Calculate joint inertia matrix for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRPR2_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:50:56
% EndTime: 2019-03-08 18:50:57
% DurationCPUTime: 0.32s
% Computational Cost: add. (303->107), mult. (761->166), div. (0->0), fcn. (844->12), ass. (0->64)
t125 = MDP(16) * pkin(9);
t135 = MDP(13) + t125;
t117 = MDP(12) - MDP(15);
t96 = sin(qJ(4));
t115 = -qJ(5) * t96 - pkin(3);
t99 = cos(qJ(4));
t81 = -pkin(4) * t99 + t115;
t121 = t81 * MDP(16);
t134 = t117 * t96 - (MDP(11) - MDP(14)) * t99 - MDP(4) + t121;
t133 = 0.2e1 * t99;
t101 = -pkin(4) - pkin(10);
t132 = pkin(5) + pkin(9);
t90 = sin(pkin(7));
t97 = sin(qJ(3));
t131 = t90 * t97;
t92 = cos(pkin(12));
t93 = cos(pkin(7));
t130 = t92 * t93;
t95 = sin(qJ(6));
t98 = cos(qJ(6));
t129 = t95 * t98;
t128 = t95 * t99;
t127 = t98 * t99;
t126 = MDP(16) * pkin(4);
t100 = cos(qJ(3));
t124 = t100 * t90;
t123 = MDP(22) * t95;
t122 = MDP(23) * t98;
t120 = qJ(5) * MDP(16);
t119 = pkin(9) ^ 2 * MDP(16);
t116 = MDP(18) * t129;
t114 = MDP(14) - t126;
t110 = -MDP(11) + t114;
t109 = -t117 + t120;
t108 = -MDP(19) * t95 - MDP(20) * t98;
t107 = MDP(22) * t98 - t95 * MDP(23);
t106 = t122 + t123;
t105 = MDP(13) + t107;
t104 = t106 + t109;
t103 = t98 * MDP(19) - t95 * MDP(20) + t107 * t101;
t94 = cos(pkin(6));
t91 = sin(pkin(6));
t89 = sin(pkin(12));
t88 = t99 ^ 2;
t87 = t98 ^ 2;
t86 = t96 ^ 2;
t85 = t95 ^ 2;
t83 = t132 * t99;
t82 = t132 * t96;
t80 = t101 * t99 + t115;
t78 = t99 * t131 + t93 * t96;
t77 = t96 * t131 - t93 * t99;
t76 = -t90 * t91 * t92 + t93 * t94;
t75 = t98 * t124 - t77 * t95;
t74 = t95 * t124 + t77 * t98;
t73 = t80 * t98 + t82 * t95;
t72 = -t80 * t95 + t82 * t98;
t71 = t94 * t131 + (t100 * t89 + t97 * t130) * t91;
t70 = t89 * t91 * t97 + (-t91 * t130 - t90 * t94) * t100;
t68 = t71 * t99 + t76 * t96;
t67 = t71 * t96 - t76 * t99;
t66 = t67 * t95 + t70 * t98;
t65 = t67 * t98 - t70 * t95;
t1 = [MDP(1) + (t94 ^ 2 + (t89 ^ 2 + t92 ^ 2) * t91 ^ 2) * MDP(2) + (t67 ^ 2 + t68 ^ 2 + t70 ^ 2) * MDP(16); t94 * MDP(2) + (-t70 * t124 + t67 * t77 + t68 * t78) * MDP(16); MDP(2) + (t100 ^ 2 * t90 ^ 2 + t77 ^ 2 + t78 ^ 2) * MDP(16); -t71 * MDP(5) + (t68 * t127 + t65 * t96) * MDP(22) + (-t68 * t128 - t66 * t96) * MDP(23) + t134 * t70 + t135 * (t67 * t96 + t68 * t99); (t78 * t127 + t74 * t96) * MDP(22) + (-t78 * t128 + t75 * t96) * MDP(23) + (-t97 * MDP(5) - t134 * t100) * t90 + t135 * (t77 * t96 + t78 * t99); pkin(3) * MDP(11) * t133 + MDP(3) + (MDP(14) * t133 + t121) * t81 + (t85 * MDP(17) + 0.2e1 * t116 + t119) * t88 + (MDP(21) + MDP(6) + t119) * t86 + 0.2e1 * (-pkin(3) * MDP(12) - t81 * MDP(15) + (MDP(7) + t108) * t99) * t96 + 0.2e1 * (t83 * t127 + t72 * t96) * MDP(22) + 0.2e1 * (-t83 * t128 - t73 * t96) * MDP(23) + 0.2e1 * (t86 + t88) * MDP(13) * pkin(9); t104 * t68 + t110 * t67; t104 * t78 + t110 * t77; t106 * t83 + (-pkin(4) * MDP(13) + MDP(8) + t103) * t96 + (MDP(9) - MDP(17) * t129 + (t85 - t87) * MDP(18) + t105 * qJ(5)) * t99 + (t109 * t99 + t110 * t96) * pkin(9); -0.2e1 * t116 + t87 * MDP(17) + MDP(10) + (-0.2e1 * MDP(14) + t126) * pkin(4) + (0.2e1 * MDP(15) + t120 + 0.2e1 * t122 + 0.2e1 * t123) * qJ(5); t67 * MDP(16); t77 * MDP(16); (t105 + t125) * t96; t114; MDP(16); MDP(22) * t65 - MDP(23) * t66; MDP(22) * t74 + MDP(23) * t75; MDP(21) * t96 + t72 * MDP(22) - t73 * MDP(23) + t108 * t99; t103; t107; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
