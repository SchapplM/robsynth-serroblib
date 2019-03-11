% Calculate Gravitation load on the joints for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:06:32
% EndTime: 2019-03-08 21:06:33
% DurationCPUTime: 0.64s
% Computational Cost: add. (295->94), mult. (548->156), div. (0->0), fcn. (623->12), ass. (0->52)
t117 = cos(pkin(6));
t93 = sin(pkin(6));
t97 = sin(qJ(2));
t125 = t93 * t97;
t96 = sin(qJ(3));
t99 = cos(qJ(3));
t132 = t117 * t99 - t96 * t125;
t124 = t93 * t99;
t100 = cos(qJ(2));
t92 = sin(pkin(10));
t112 = t92 * t117;
t116 = cos(pkin(10));
t82 = t116 * t100 - t97 * t112;
t131 = t92 * t124 - t82 * t96;
t130 = g(3) * t93;
t91 = qJ(3) + pkin(11);
t89 = sin(t91);
t95 = sin(qJ(6));
t128 = t89 * t95;
t98 = cos(qJ(6));
t127 = t89 * t98;
t126 = t92 * t93;
t94 = -qJ(4) - pkin(8);
t123 = t94 * t97;
t106 = t117 * t116;
t79 = -t100 * t106 + t92 * t97;
t80 = t92 * t100 + t97 * t106;
t88 = pkin(3) * t99 + pkin(2);
t122 = -t79 * t88 - t80 * t94;
t81 = t100 * t112 + t116 * t97;
t121 = -t81 * t88 - t82 * t94;
t120 = t100 * t93;
t119 = t100 * t95;
t118 = t100 * t98;
t115 = MDP(13) + MDP(17);
t111 = t93 * t116;
t109 = t131 * pkin(3);
t90 = cos(t91);
t108 = pkin(4) * t90 + qJ(5) * t89;
t107 = t132 * pkin(3);
t104 = -t99 * t111 - t80 * t96;
t102 = -g(1) * t81 - g(2) * t79 + g(3) * t120;
t101 = t104 * pkin(3);
t83 = t88 * t120;
t76 = t117 * t89 + t90 * t125;
t75 = -t117 * t90 + t89 * t125;
t72 = t89 * t126 + t82 * t90;
t71 = -t90 * t126 + t82 * t89;
t70 = -t89 * t111 + t80 * t90;
t69 = t90 * t111 + t80 * t89;
t66 = -g(1) * t71 - g(2) * t69 - g(3) * t75;
t1 = [(-MDP(1) - t115) * g(3); (-g(1) * t121 - g(2) * t122 - g(3) * (-t93 * t123 + t83)) * MDP(13) + (-g(1) * (-t108 * t81 + t121) - g(2) * (-t108 * t79 + t122) - g(3) * t83 - (t108 * t100 - t123) * t130) * MDP(17) + (-g(1) * (-t81 * t128 + t82 * t98) - g(2) * (-t79 * t128 + t80 * t98) - (t89 * t119 + t97 * t98) * t130) * MDP(23) + (-g(1) * (-t81 * t127 - t82 * t95) - g(2) * (-t79 * t127 - t80 * t95) - (t89 * t118 - t95 * t97) * t130) * MDP(24) + (-MDP(10) * t99 + MDP(11) * t96 + t90 * MDP(15) - t89 * MDP(16) - MDP(3)) * t102 + (MDP(4) - MDP(12) - MDP(14)) * (g(1) * t82 + g(2) * t80 + g(3) * t125); (-g(1) * t131 - g(2) * t104 - g(3) * t132) * MDP(10) + (-g(1) * (-t96 * t126 - t82 * t99) - g(2) * (t96 * t111 - t80 * t99) - g(3) * (-t117 * t96 - t97 * t124)) * MDP(11) + (-g(1) * t109 - g(2) * t101 - g(3) * t107) * MDP(13) + t66 * MDP(15) + (-g(1) * (-pkin(4) * t71 + qJ(5) * t72 + t109) - g(2) * (-t69 * pkin(4) + t70 * qJ(5) + t101) - g(3) * (-pkin(4) * t75 + qJ(5) * t76 + t107)) * MDP(17) + (MDP(23) * t95 + MDP(24) * t98 + MDP(16)) * (-g(1) * t72 - g(2) * t70 - g(3) * t76); t115 * t102; t66 * MDP(17); (-g(1) * (t71 * t98 - t81 * t95) - g(2) * (t69 * t98 - t79 * t95) - g(3) * (t93 * t119 + t75 * t98)) * MDP(23) + (-g(1) * (-t71 * t95 - t81 * t98) - g(2) * (-t69 * t95 - t79 * t98) - g(3) * (t93 * t118 - t75 * t95)) * MDP(24);];
taug  = t1;
