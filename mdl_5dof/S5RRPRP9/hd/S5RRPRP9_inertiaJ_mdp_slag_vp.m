% Calculate joint inertia matrix for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP9_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:30
% EndTime: 2019-12-31 20:07:31
% DurationCPUTime: 0.47s
% Computational Cost: add. (557->140), mult. (1082->203), div. (0->0), fcn. (1041->6), ass. (0->56)
t125 = cos(qJ(4));
t96 = sin(pkin(8));
t97 = cos(pkin(8));
t98 = sin(qJ(4));
t134 = t125 * t97 - t98 * t96;
t133 = (MDP(14) * qJ(3));
t132 = 2 * MDP(13);
t131 = -2 * MDP(16);
t130 = 2 * MDP(21);
t129 = 2 * MDP(22);
t128 = 2 * MDP(23);
t127 = 2 * MDP(24);
t126 = pkin(6) * t96;
t100 = cos(qJ(2));
t124 = pkin(6) * t100;
t123 = t100 * pkin(4);
t99 = sin(qJ(2));
t122 = t96 * t99;
t120 = pkin(7) + qJ(3);
t85 = -t100 * pkin(2) - t99 * qJ(3) - pkin(1);
t80 = t97 * t85;
t70 = -t97 * t99 * pkin(7) + t80 + (-pkin(3) - t126) * t100;
t76 = t97 * t124 + t96 * t85;
t72 = -pkin(7) * t122 + t76;
t65 = t125 * t72 + t98 * t70;
t84 = pkin(3) * t122 + t99 * pkin(6);
t118 = pkin(2) * MDP(14);
t117 = t100 * qJ(5);
t83 = t125 * t96 + t98 * t97;
t77 = t83 * t99;
t116 = t77 * MDP(18);
t78 = t134 * t99;
t115 = t78 * MDP(15);
t114 = t78 * MDP(17);
t113 = t96 * MDP(12);
t112 = t97 * MDP(11);
t111 = t100 * MDP(19);
t110 = MDP(20) + MDP(22);
t109 = MDP(21) - MDP(24);
t91 = -t97 * pkin(3) - pkin(2);
t107 = t120 * t96;
t64 = t125 * t70 - t98 * t72;
t106 = -MDP(25) * pkin(4) - MDP(22);
t104 = t96 * MDP(11) + t97 * MDP(12);
t103 = t83 * MDP(17) + MDP(18) * t134;
t102 = -t112 + t113 - t118;
t95 = t99 ^ 2;
t86 = t120 * t97;
t75 = -t96 * t124 + t80;
t74 = -t98 * t107 + t125 * t86;
t73 = t125 * t107 + t98 * t86;
t69 = -pkin(4) * t134 - t83 * qJ(5) + t91;
t66 = t77 * pkin(4) - t78 * qJ(5) + t84;
t63 = -t64 + t123;
t62 = -t117 + t65;
t1 = [MDP(1) + t95 * MDP(4) + (t95 * pkin(6) ^ 2 + t75 ^ 2 + t76 ^ 2) * MDP(14) + (t62 ^ 2 + t63 ^ 2 + t66 ^ 2) * MDP(25) + (t77 * t131 + t115) * t78 + (0.2e1 * pkin(1) * MDP(9) + t111 - 0.2e1 * t114 + 0.2e1 * t116) * t100 + 0.2e1 * (-t75 * t100 + t95 * t126) * MDP(11) + 0.2e1 * (t95 * pkin(6) * t97 + t76 * t100) * MDP(12) + 0.2e1 * (-t64 * t100 + t84 * t77) * MDP(20) + (t65 * t100 + t84 * t78) * t130 + (t63 * t100 + t66 * t77) * t129 + (-t62 * t77 + t63 * t78) * t128 + (-t62 * t100 - t66 * t78) * t127 + (-0.2e1 * pkin(1) * MDP(10) + 0.2e1 * t100 * MDP(5) + (-t75 * t97 - t76 * t96) * t132) * t99; t83 * t115 + (t134 * t78 - t83 * t77) * MDP(16) + (-t134 * t84 + t91 * t77) * MDP(20) + (t91 * t78 + t84 * t83) * MDP(21) + (-t134 * t66 + t69 * t77) * MDP(22) + (t134 * t62 + t63 * t83 + t73 * t78 - t74 * t77) * MDP(23) + (-t66 * t83 - t69 * t78) * MDP(24) + (t62 * t74 + t63 * t73 + t66 * t69) * MDP(25) + (MDP(6) - t104 * pkin(2) + (-MDP(9) + t102) * pkin(6)) * t99 + (-pkin(6) * MDP(10) + t104 * qJ(3) + t109 * t74 + t110 * t73 + MDP(7) - t103) * t100 + (MDP(13) + t133) * (-t75 * t96 + t76 * t97); MDP(8) + (t69 ^ 2 + t73 ^ 2 + t74 ^ 2) * MDP(25) + (0.2e1 * t112 - 0.2e1 * t113 + t118) * pkin(2) - 0.2e1 * (t91 * MDP(20) + t69 * MDP(22) - t74 * MDP(23)) * t134 + (MDP(15) * t83 - 0.2e1 * t69 * MDP(24) + t73 * t128 + t91 * t130 - t131 * t134) * t83 + (t132 + t133) * (t96 ^ 2 + t97 ^ 2) * qJ(3); t66 * MDP(25) + t109 * t78 + t110 * t77 + (pkin(6) * MDP(14) + t104) * t99; t69 * MDP(25) + t109 * t83 - t110 * t134 + t102; MDP(14) + MDP(25); t114 - t116 - t111 + t64 * MDP(20) - t65 * MDP(21) + (t64 - 0.2e1 * t123) * MDP(22) + (-pkin(4) * t78 - t77 * qJ(5)) * MDP(23) + (-0.2e1 * t117 + t65) * MDP(24) + (-t63 * pkin(4) + t62 * qJ(5)) * MDP(25); (-pkin(4) * t83 + qJ(5) * t134) * MDP(23) + (MDP(25) * qJ(5) - t109) * t74 + (-MDP(20) + t106) * t73 + t103; 0; MDP(19) + pkin(4) * t129 + qJ(5) * t127 + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(25); t100 * MDP(22) + t78 * MDP(23) + t63 * MDP(25); t83 * MDP(23) + t73 * MDP(25); 0; t106; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
