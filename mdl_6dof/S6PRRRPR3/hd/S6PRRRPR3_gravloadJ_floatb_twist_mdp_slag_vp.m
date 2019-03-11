% Calculate Gravitation load on the joints for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:14:02
% EndTime: 2019-03-08 23:14:04
% DurationCPUTime: 0.60s
% Computational Cost: add. (411->82), mult. (688->134), div. (0->0), fcn. (803->12), ass. (0->42)
t126 = MDP(17) - MDP(20);
t125 = MDP(18) - MDP(21);
t93 = sin(pkin(6));
t122 = g(3) * t93;
t91 = qJ(3) + qJ(4);
t89 = sin(t91);
t94 = sin(qJ(6));
t121 = t89 * t94;
t97 = cos(qJ(6));
t120 = t89 * t97;
t92 = sin(pkin(11));
t119 = t92 * t93;
t96 = sin(qJ(2));
t118 = t93 * t96;
t98 = cos(qJ(3));
t117 = t93 * t98;
t99 = cos(qJ(2));
t116 = t94 * t99;
t115 = t97 * t99;
t114 = cos(pkin(6));
t113 = cos(pkin(11));
t112 = t92 * t114;
t111 = t93 * t113;
t109 = t114 * t113;
t80 = t96 * t109 + t92 * t99;
t82 = -t96 * t112 + t113 * t99;
t110 = g(1) * t82 + g(2) * t80;
t90 = cos(t91);
t72 = t90 * t111 + t80 * t89;
t74 = -t90 * t119 + t82 * t89;
t77 = -t114 * t90 + t89 * t118;
t106 = g(1) * t74 + g(2) * t72 + g(3) * t77;
t73 = -t89 * t111 + t80 * t90;
t75 = t89 * t119 + t82 * t90;
t78 = t114 * t89 + t90 * t118;
t108 = (-t94 * MDP(28) - t97 * MDP(29) + t125) * (g(1) * t75 + g(2) * t73 + g(3) * t78) + t126 * t106;
t102 = -g(1) * (-t74 * pkin(4) + qJ(5) * t75) - g(2) * (-t72 * pkin(4) + t73 * qJ(5)) - g(3) * (-t77 * pkin(4) + t78 * qJ(5));
t95 = sin(qJ(3));
t101 = -g(1) * (t92 * t117 - t82 * t95) - g(2) * (-t98 * t111 - t80 * t95) - g(3) * (t114 * t98 - t95 * t118);
t81 = t99 * t112 + t113 * t96;
t79 = -t99 * t109 + t92 * t96;
t1 = [(-MDP(1) - MDP(22)) * g(3); (t96 * t122 + t110) * (-pkin(9) - pkin(8)) * MDP(22) + (-g(1) * (-t81 * t121 + t82 * t97) - g(2) * (-t79 * t121 + t80 * t97) - (t89 * t116 + t96 * t97) * t122) * MDP(28) + (-g(1) * (-t81 * t120 - t82 * t94) - g(2) * (-t79 * t120 - t80 * t94) - (t89 * t115 - t94 * t96) * t122) * MDP(29) + (MDP(4) - MDP(19)) * (g(3) * t118 + t110) + (-t126 * t90 + t125 * t89 - MDP(3) - t98 * MDP(10) + t95 * MDP(11) + (-pkin(3) * t98 - pkin(4) * t90 - qJ(5) * t89 - pkin(2)) * MDP(22)) * (-g(1) * t81 - g(2) * t79 + t99 * t122); t101 * MDP(10) + (-g(1) * (-t95 * t119 - t82 * t98) - g(2) * (t95 * t111 - t80 * t98) - g(3) * (-t114 * t95 - t96 * t117)) * MDP(11) + (t101 * pkin(3) + t102) * MDP(22) + t108; t102 * MDP(22) + t108; -t106 * MDP(22); (-g(1) * (t74 * t97 - t81 * t94) - g(2) * (t72 * t97 - t79 * t94) - g(3) * (t93 * t116 + t77 * t97)) * MDP(28) + (-g(1) * (-t74 * t94 - t81 * t97) - g(2) * (-t72 * t94 - t79 * t97) - g(3) * (t93 * t115 - t77 * t94)) * MDP(29);];
taug  = t1;
