% Calculate Gravitation load on the joints for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:37:02
% EndTime: 2019-03-09 09:37:03
% DurationCPUTime: 0.46s
% Computational Cost: add. (233->82), mult. (325->115), div. (0->0), fcn. (301->10), ass. (0->46)
t125 = MDP(10) - MDP(13);
t124 = MDP(9) - MDP(12) + MDP(17);
t101 = cos(qJ(1));
t99 = sin(qJ(1));
t79 = g(1) * t101 + g(2) * t99;
t98 = sin(qJ(2));
t123 = t79 * t98;
t100 = cos(qJ(2));
t89 = t98 * qJ(3);
t109 = t100 * pkin(2) + t89;
t76 = g(3) * t98 + t79 * t100;
t122 = g(1) * t99;
t118 = g(3) * t100;
t95 = pkin(10) + qJ(5);
t88 = qJ(6) + t95;
t84 = sin(t88);
t117 = t99 * t84;
t85 = cos(t88);
t116 = t99 * t85;
t86 = sin(t95);
t115 = t99 * t86;
t87 = cos(t95);
t114 = t99 * t87;
t96 = sin(pkin(10));
t113 = t99 * t96;
t97 = cos(pkin(10));
t112 = t99 * t97;
t108 = t101 * t98;
t67 = t85 * t108 - t117;
t68 = t84 * t108 + t116;
t69 = t101 * t84 + t98 * t116;
t70 = t101 * t85 - t98 * t117;
t110 = (-g(1) * t67 - g(2) * t69 + t85 * t118) * MDP(31) + (g(1) * t68 - g(2) * t70 - t84 * t118) * MDP(32);
t107 = qJ(4) * t100;
t106 = t100 * t101;
t105 = pkin(2) * t106 + t99 * pkin(7) + (pkin(1) + t89) * t101;
t103 = -pkin(1) - t109;
t92 = t101 * pkin(7);
t82 = qJ(3) * t106;
t80 = t99 * t100 * qJ(3);
t75 = -t118 + t123;
t74 = t101 * t87 - t98 * t115;
t73 = t101 * t86 + t98 * t114;
t72 = t86 * t108 + t114;
t71 = t87 * t108 - t115;
t1 = [(-g(1) * t92 - g(2) * t105 - t103 * t122) * MDP(14) + (-g(1) * (t101 * t97 - t98 * t113) - g(2) * (t96 * t108 + t112)) * MDP(15) + (-g(1) * (-t101 * t96 - t98 * t112) - g(2) * (t97 * t108 - t113)) * MDP(16) + (-g(1) * (pkin(3) * t101 + t92) - g(2) * (qJ(4) * t106 + t105) + (-g(1) * (t103 - t107) - g(2) * pkin(3)) * t99) * MDP(18) + (-g(1) * t74 - g(2) * t72) * MDP(24) + (g(1) * t73 - g(2) * t71) * MDP(25) + (-g(1) * t70 - g(2) * t68) * MDP(31) + (g(1) * t69 - g(2) * t67) * MDP(32) + (MDP(3) - MDP(11)) * t79 + (t124 * t100 - t125 * t98 + MDP(2)) * (-g(2) * t101 + t122); (-g(1) * (-pkin(2) * t108 + t82) - g(2) * (-pkin(2) * t98 * t99 + t80) - g(3) * t109) * MDP(14) + (-g(1) * t82 - g(2) * t80 - g(3) * (t107 + t109) + (pkin(2) + qJ(4)) * t123) * MDP(18) + t124 * t75 + (-t96 * MDP(15) - t97 * MDP(16) - t86 * MDP(24) - t87 * MDP(25) - t84 * MDP(31) - t85 * MDP(32) + t125) * t76; (-MDP(14) - MDP(18)) * t75; -t76 * MDP(18); (-g(1) * t71 - g(2) * t73 + t87 * t118) * MDP(24) + (g(1) * t72 - g(2) * t74 - t86 * t118) * MDP(25) + t110; t110;];
taug  = t1;
