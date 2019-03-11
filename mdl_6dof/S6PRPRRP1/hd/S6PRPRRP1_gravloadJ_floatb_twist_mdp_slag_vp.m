% Calculate Gravitation load on the joints for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S6PRPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:58:30
% EndTime: 2019-03-08 19:58:32
% DurationCPUTime: 0.68s
% Computational Cost: add. (328->87), mult. (830->149), div. (0->0), fcn. (1025->12), ass. (0->51)
t141 = MDP(12) - MDP(20);
t100 = sin(pkin(10));
t102 = cos(pkin(10));
t107 = sin(qJ(2));
t103 = cos(pkin(6));
t110 = cos(qJ(2));
t125 = t103 * t110;
t140 = -t100 * t125 - t102 * t107;
t105 = sin(qJ(5));
t133 = cos(pkin(11));
t119 = t110 * t133;
t99 = sin(pkin(11));
t117 = -t107 * t99 + t119;
t91 = -t107 * t133 - t110 * t99;
t89 = t91 * t103;
t76 = -t100 * t117 + t102 * t89;
t136 = t105 * t76;
t79 = -t100 * t89 - t102 * t117;
t135 = t105 * t79;
t101 = sin(pkin(6));
t88 = t91 * t101;
t134 = t105 * t88;
t132 = t100 * t107;
t106 = sin(qJ(4));
t131 = t101 * t106;
t130 = t101 * t107;
t109 = cos(qJ(4));
t129 = t101 * t109;
t128 = t101 * t110;
t126 = t103 * t107;
t124 = t105 * t109;
t108 = cos(qJ(5));
t123 = t108 * t109;
t122 = MDP(21) + MDP(5);
t120 = t102 * t125;
t72 = t102 * t129 - t106 * t76;
t74 = -t100 * t129 - t106 * t79;
t82 = -t103 * t109 - t106 * t88;
t116 = g(1) * t74 + g(2) * t72 + g(3) * t82;
t113 = t103 * t117;
t112 = -g(1) * t140 - g(3) * t128;
t104 = -qJ(6) - pkin(9);
t98 = pkin(5) * t108 + pkin(4);
t92 = pkin(2) * t120;
t87 = -t101 * t119 + t99 * t130;
t83 = t103 * t106 - t109 * t88;
t80 = -t100 * t113 + t102 * t91;
t77 = t100 * t91 + t102 * t113;
t75 = t100 * t131 - t109 * t79;
t73 = -t102 * t131 - t109 * t76;
t1 = [(-MDP(1) - t122) * g(3); (-g(2) * (t120 - t132) + t112) * MDP(3) + (-g(1) * (t100 * t126 - t102 * t110) - g(2) * (-t100 * t110 - t102 * t126) + g(3) * t130) * MDP(4) + (-g(2) * t92 + (g(2) * t132 + t112) * pkin(2)) * MDP(5) + (-g(1) * (t80 * t123 - t135) - g(2) * (t77 * t123 - t136) - g(3) * (-t87 * t123 - t134)) * MDP(18) + (-g(1) * (-t108 * t79 - t80 * t124) - g(2) * (-t108 * t76 - t77 * t124) - g(3) * (-t108 * t88 + t87 * t124)) * MDP(19) + (-g(1) * (pkin(2) * t140 - pkin(5) * t135 - t79 * pkin(8)) - g(2) * (-pkin(2) * t132 - pkin(5) * t136 - pkin(8) * t76 + t92) - g(3) * (pkin(2) * t128 - pkin(5) * t134 - t88 * pkin(8))) * MDP(21) + (-t109 * MDP(11) + (t104 * t106 - t109 * t98 - pkin(3)) * MDP(21) + t141 * t106) * (g(1) * t80 + g(2) * t77 - g(3) * t87); t122 * (-g(3) * t103 + (-g(1) * t100 + g(2) * t102) * t101); (-g(1) * (-t104 * t75 - t74 * t98) - g(2) * (-t104 * t73 - t72 * t98) - g(3) * (-t104 * t83 - t82 * t98)) * MDP(21) + t141 * (g(1) * t75 + g(2) * t73 + g(3) * t83) + (MDP(18) * t108 - MDP(19) * t105 + MDP(11)) * t116; (-g(1) * (t105 * t80 - t108 * t75) - g(2) * (t105 * t77 - t108 * t73) - g(3) * (-t105 * t87 - t108 * t83)) * MDP(19) + (pkin(5) * MDP(21) + MDP(18)) * (-g(1) * (-t105 * t75 - t108 * t80) - g(2) * (-t105 * t73 - t108 * t77) - g(3) * (-t105 * t83 + t108 * t87)); -t116 * MDP(21);];
taug  = t1;
