% Calculate joint inertia matrix for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRRP2_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:03:07
% EndTime: 2019-03-08 20:03:09
% DurationCPUTime: 0.46s
% Computational Cost: add. (457->135), mult. (979->208), div. (0->0), fcn. (998->10), ass. (0->58)
t140 = pkin(9) * MDP(23) + MDP(21);
t129 = MDP(18) + MDP(20);
t104 = sin(qJ(5));
t107 = cos(qJ(5));
t128 = MDP(19) - MDP(22);
t121 = -t107 * pkin(5) - t104 * qJ(6);
t90 = -pkin(4) + t121;
t139 = t90 * MDP(23) + t128 * t104 - t129 * t107 - MDP(11);
t108 = cos(qJ(4));
t105 = sin(qJ(4));
t102 = cos(pkin(11));
t94 = -t102 * pkin(2) - pkin(3);
t88 = -t108 * pkin(4) - t105 * pkin(9) + t94;
t100 = sin(pkin(11));
t93 = t100 * pkin(2) + pkin(8);
t78 = t107 * t108 * t93 + t104 * t88;
t96 = t104 ^ 2;
t98 = t107 ^ 2;
t138 = t96 + t98;
t136 = t104 * t93;
t135 = t104 * t105;
t134 = t104 * t108;
t133 = t105 * t107;
t132 = t104 * MDP(15);
t131 = t107 * MDP(14);
t130 = t108 * MDP(11);
t127 = t138 * pkin(9);
t125 = -MDP(23) * pkin(5) - MDP(20);
t123 = -0.2e1 * qJ(6) * MDP(22) - MDP(17);
t122 = -MDP(18) + t125;
t120 = -pkin(5) * t104 + t107 * qJ(6);
t103 = cos(pkin(6));
t101 = sin(pkin(6));
t106 = sin(qJ(2));
t109 = cos(qJ(2));
t85 = (t100 * t109 + t102 * t106) * t101;
t80 = t103 * t105 + t85 * t108;
t83 = (t100 * t106 - t102 * t109) * t101;
t72 = t80 * t104 - t83 * t107;
t73 = t83 * t104 + t80 * t107;
t119 = t72 * t104 + t73 * t107;
t75 = -t108 * qJ(6) + t78;
t87 = t107 * t88;
t76 = -t87 + (pkin(5) + t136) * t108;
t118 = t76 * t104 + t75 * t107;
t117 = pkin(4) * MDP(18) - t90 * MDP(20);
t116 = -pkin(4) * MDP(19) - t90 * MDP(22);
t115 = MDP(23) * qJ(6) - t128;
t114 = (-t93 * t134 + t87) * MDP(18) - t78 * MDP(19);
t113 = t107 * MDP(15) - t104 * MDP(16);
t112 = t122 * t104 + t115 * t107;
t99 = t108 ^ 2;
t97 = t105 ^ 2;
t95 = t98 * t105;
t92 = pkin(9) * t134;
t81 = (-t120 + t93) * t105;
t79 = -t103 * t108 + t85 * t105;
t1 = [MDP(1) + (t103 ^ 2 + t83 ^ 2 + t85 ^ 2) * MDP(5) + (t72 ^ 2 + t73 ^ 2 + t79 ^ 2) * MDP(23); -t83 * t130 + (t72 * t76 + t73 * t75 + t79 * t81) * MDP(23) + t128 * (t73 * t108 + t79 * t133) + (t109 * MDP(3) - t106 * MDP(4)) * t101 + (t100 * t85 - t102 * t83) * MDP(5) * pkin(2) + (t83 * MDP(12) + (-t104 * t73 + t107 * t72) * MDP(21)) * t105 + t129 * (t72 * t108 + t79 * t135); MDP(2) - 0.2e1 * t94 * t130 + t99 * MDP(17) + (t75 ^ 2 + t76 ^ 2 + t81 ^ 2) * MDP(23) + (t100 ^ 2 + t102 ^ 2) * MDP(5) * pkin(2) ^ 2 + 0.2e1 * (t76 * MDP(20) - t75 * MDP(22) - t114) * t108 + (t98 * MDP(13) - 0.2e1 * t104 * t131 + MDP(6) + 0.2e1 * (t104 * MDP(18) + t107 * MDP(19)) * t93) * t97 + 0.2e1 * (t94 * MDP(12) + (MDP(7) - t113) * t108 + (-t104 * t75 + t107 * t76) * MDP(21) + (t104 * MDP(20) - t107 * MDP(22)) * t81) * t105; t103 * MDP(5) + (t119 * t105 - t79 * t108) * MDP(23); (t118 * t105 - t81 * t108) * MDP(23); MDP(5) + (t138 * t97 + t99) * MDP(23); -t80 * MDP(12) + t140 * t119 + t139 * t79; t95 * MDP(14) + t92 * MDP(18) + (-t81 * t107 + t92) * MDP(20) + t118 * MDP(21) - t81 * t104 * MDP(22) + (t118 * pkin(9) + t81 * t90) * MDP(23) + (-t93 * MDP(12) - t132 + MDP(9) + (t128 * pkin(9) - MDP(16)) * t107) * t108 + (-t93 * MDP(11) - t96 * MDP(14) + MDP(8) + (-t93 * MDP(18) + t116) * t107 + (t107 * MDP(13) + t93 * MDP(19) - t117) * t104) * t105; t95 * MDP(21) + (t96 * MDP(21) + MDP(23) * t127 - MDP(12)) * t105 - t139 * t108; MDP(10) + t96 * MDP(13) + (t138 * pkin(9) ^ 2 + t90 ^ 2) * MDP(23) + 0.2e1 * MDP(21) * t127 + 0.2e1 * t117 * t107 + 0.2e1 * (t116 + t131) * t104; t115 * t73 + t122 * t72; t87 * MDP(20) + t78 * MDP(22) + (-t76 * pkin(5) + t75 * qJ(6)) * MDP(23) + ((-0.2e1 * pkin(5) - t136) * MDP(20) + t123) * t108 + (t121 * MDP(21) + t113) * t105 + t114; t112 * t105; t107 * MDP(16) + t120 * MDP(21) + t112 * pkin(9) + t132; 0.2e1 * pkin(5) * MDP(20) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(23) - t123; t72 * MDP(23); t108 * MDP(20) + MDP(21) * t133 + t76 * MDP(23); MDP(23) * t135; t140 * t104; t125; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
