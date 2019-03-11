% Calculate Gravitation load on the joints for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:59:18
% EndTime: 2019-03-08 23:59:20
% DurationCPUTime: 0.45s
% Computational Cost: add. (412->87), mult. (685->142), div. (0->0), fcn. (795->12), ass. (0->46)
t130 = MDP(18) - MDP(26);
t101 = sin(qJ(2));
t104 = cos(qJ(2));
t119 = cos(pkin(11));
t120 = cos(pkin(6));
t113 = t120 * t119;
t96 = sin(pkin(11));
t82 = t101 * t96 - t104 * t113;
t129 = g(2) * t82;
t83 = t101 * t113 + t96 * t104;
t128 = g(2) * t83;
t97 = sin(pkin(6));
t127 = g(3) * t97;
t95 = qJ(3) + qJ(4);
t94 = cos(t95);
t99 = sin(qJ(5));
t126 = t94 * t99;
t125 = t96 * t97;
t124 = t101 * t97;
t102 = cos(qJ(5));
t123 = t102 * t94;
t103 = cos(qJ(3));
t122 = t103 * t97;
t121 = t104 * t99;
t118 = t102 * t104;
t117 = pkin(5) * t99 + pkin(8) + pkin(9);
t116 = t96 * t120;
t115 = t97 * t119;
t93 = sin(t95);
t74 = t94 * t115 + t83 * t93;
t85 = -t101 * t116 + t119 * t104;
t76 = -t94 * t125 + t85 * t93;
t80 = -t120 * t94 + t93 * t124;
t111 = g(1) * t76 + g(2) * t74 + g(3) * t80;
t75 = -t93 * t115 + t83 * t94;
t77 = t93 * t125 + t85 * t94;
t81 = t120 * t93 + t94 * t124;
t114 = t130 * (g(1) * t77 + g(2) * t75 + g(3) * t81) + (t102 * MDP(24) - t99 * MDP(25) + MDP(17)) * t111;
t91 = pkin(5) * t102 + pkin(4);
t98 = -qJ(6) - pkin(10);
t112 = pkin(3) * t103 + t91 * t94 - t93 * t98 + pkin(2);
t108 = -g(1) * (-t76 * t91 - t77 * t98) - g(2) * (-t74 * t91 - t75 * t98) - g(3) * (-t80 * t91 - t81 * t98);
t100 = sin(qJ(3));
t106 = -g(1) * (-t100 * t85 + t96 * t122) - g(2) * (-t83 * t100 - t103 * t115) - g(3) * (-t100 * t124 + t120 * t103);
t84 = t119 * t101 + t104 * t116;
t1 = [(-MDP(1) - MDP(27)) * g(3); (g(1) * t85 + g(3) * t124 + t128) * MDP(4) + (-g(1) * (-t84 * t123 + t85 * t99) - g(2) * (-t82 * t123 + t83 * t99) - (t101 * t99 + t94 * t118) * t127) * MDP(24) + (-g(1) * (t102 * t85 + t84 * t126) - g(2) * (t102 * t83 + t82 * t126) - (t101 * t102 - t94 * t121) * t127) * MDP(25) + (-g(1) * (-t112 * t84 + t117 * t85) - t117 * t128 + t112 * t129 - (t117 * t101 + t112 * t104) * t127) * MDP(27) + (-t103 * MDP(10) + t100 * MDP(11) - t94 * MDP(17) + t130 * t93 - MDP(3)) * (-g(1) * t84 + t104 * t127 - t129); t106 * MDP(10) + (-g(1) * (-t100 * t125 - t103 * t85) - g(2) * (t100 * t115 - t83 * t103) - g(3) * (-t120 * t100 - t101 * t122)) * MDP(11) + (t106 * pkin(3) + t108) * MDP(27) + t114; t108 * MDP(27) + t114; (-g(1) * (-t102 * t77 - t84 * t99) - g(2) * (-t102 * t75 - t82 * t99) - g(3) * (-t102 * t81 + t97 * t121)) * MDP(25) + (pkin(5) * MDP(27) + MDP(24)) * (-g(1) * (t102 * t84 - t77 * t99) - g(2) * (t102 * t82 - t75 * t99) - g(3) * (-t97 * t118 - t81 * t99)); -t111 * MDP(27);];
taug  = t1;
