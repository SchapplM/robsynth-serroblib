% Calculate Gravitation load on the joints for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:48:24
% EndTime: 2019-03-09 02:48:25
% DurationCPUTime: 0.54s
% Computational Cost: add. (282->79), mult. (394->113), div. (0->0), fcn. (390->10), ass. (0->45)
t100 = cos(qJ(1));
t98 = sin(qJ(1));
t80 = g(1) * t100 + g(2) * t98;
t91 = pkin(9) + qJ(3);
t88 = sin(t91);
t124 = t80 * t88;
t89 = cos(t91);
t68 = -g(3) * t89 + t124;
t127 = MDP(15) + MDP(19);
t126 = -MDP(16) + MDP(21);
t125 = MDP(14) - MDP(17) - MDP(20);
t115 = t89 * pkin(3) + t88 * qJ(4);
t121 = g(3) * t88;
t119 = pkin(3) * t100;
t92 = sin(pkin(10));
t117 = t98 * t92;
t94 = cos(pkin(10));
t116 = t98 * t94;
t114 = t100 * t92;
t113 = t100 * t94;
t96 = -pkin(7) - qJ(2);
t112 = t100 * t96;
t111 = qJ(4) * t100;
t110 = -MDP(18) - MDP(22);
t95 = cos(pkin(9));
t87 = pkin(2) * t95 + pkin(1);
t109 = t100 * t87 + t88 * t111 + t89 * t119;
t79 = g(1) * t98 - g(2) * t100;
t107 = pkin(4) * t94 + qJ(5) * t92;
t71 = t117 * t89 + t113;
t72 = t116 * t89 - t114;
t97 = sin(qJ(6));
t99 = cos(qJ(6));
t106 = t71 * t99 - t72 * t97;
t105 = t71 * t97 + t72 * t99;
t104 = t92 * t99 - t94 * t97;
t103 = t92 * t97 + t94 * t99;
t101 = (-g(1) * (-t87 - t115) + g(2) * t96) * t98;
t77 = t89 * t111;
t75 = t98 * t89 * qJ(4);
t74 = t113 * t89 + t117;
t73 = t114 * t89 - t116;
t64 = t73 * t97 + t74 * t99;
t63 = t73 * t99 - t74 * t97;
t1 = [(-g(1) * (-t98 * pkin(1) + qJ(2) * t100) - g(2) * (pkin(1) * t100 + t98 * qJ(2))) * MDP(7) + (g(1) * t112 - g(2) * t109 + t101) * MDP(18) + (-g(1) * (-t72 * pkin(4) - t71 * qJ(5) - t112) - g(2) * (pkin(4) * t74 + qJ(5) * t73 + t109) + t101) * MDP(22) + (g(1) * t105 - g(2) * t64) * MDP(28) + (g(1) * t106 - g(2) * t63) * MDP(29) + (MDP(3) - MDP(6)) * t80 + t127 * (g(1) * t72 - g(2) * t74) + t126 * (g(1) * t71 - g(2) * t73) + (-t125 * t88 + t89 * MDP(13) + MDP(4) * t95 - MDP(5) * sin(pkin(9)) + MDP(2)) * t79; (-MDP(7) + t110) * t79; (-g(1) * (-t119 * t88 + t77) - g(2) * (-pkin(3) * t88 * t98 + t75) - g(3) * t115) * MDP(18) + (-g(1) * t77 - g(2) * t75 - g(3) * (t107 * t89 + t115) + (pkin(3) + t107) * t124) * MDP(22) + t125 * (t80 * t89 + t121) + (t103 * MDP(28) + t104 * MDP(29) + t126 * t92 + t127 * t94 + MDP(13)) * t68; t110 * t68; (-g(1) * t73 - g(2) * t71 - t121 * t92) * MDP(22); (-g(1) * t63 - g(2) * t106) * MDP(28) + (g(1) * t64 + g(2) * t105) * MDP(29) + (-MDP(28) * t104 + MDP(29) * t103) * t121;];
taug  = t1;
