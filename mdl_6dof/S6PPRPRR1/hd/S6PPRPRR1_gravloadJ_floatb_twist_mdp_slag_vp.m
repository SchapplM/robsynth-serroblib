% Calculate Gravitation load on the joints for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRPRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S6PPRPRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:43:39
% EndTime: 2019-03-08 18:43:39
% DurationCPUTime: 0.39s
% Computational Cost: add. (406->78), mult. (1126->148), div. (0->0), fcn. (1459->16), ass. (0->50)
t101 = sin(pkin(11));
t102 = sin(pkin(7));
t105 = cos(pkin(12));
t107 = cos(pkin(7));
t108 = cos(pkin(6));
t103 = sin(pkin(6));
t131 = t103 * t107;
t106 = cos(pkin(11));
t132 = t103 * t106;
t133 = t102 * t103;
t100 = sin(pkin(12));
t130 = t106 * t108;
t92 = -t101 * t100 + t105 * t130;
t134 = t101 * t108;
t94 = -t106 * t100 - t105 * t134;
t136 = g(1) * (t101 * t133 + t107 * t94) - g(2) * (t102 * t132 - t107 * t92) + g(3) * (t102 * t108 + t105 * t131);
t135 = t101 * t103;
t109 = sin(qJ(6));
t113 = cos(qJ(5));
t129 = t109 * t113;
t112 = cos(qJ(6));
t128 = t112 * t113;
t127 = MDP(2) + MDP(6);
t104 = cos(pkin(13));
t111 = sin(qJ(3));
t114 = cos(qJ(3));
t99 = sin(pkin(13));
t126 = t114 * t104 - t111 * t99;
t125 = t111 * t104 + t114 * t99;
t93 = t100 * t130 + t101 * t105;
t95 = -t100 * t134 + t106 * t105;
t119 = g(3) * t100 * t103 + g(1) * t95 + g(2) * t93;
t88 = t125 * t102;
t90 = t125 * t107;
t118 = t126 * t95 + t88 * t135 + t94 * t90;
t117 = t126 * t93 - t88 * t132 + t92 * t90;
t116 = t108 * t88 + (t100 * t126 + t105 * t90) * t103;
t110 = sin(qJ(5));
t91 = -t105 * t133 + t108 * t107;
t89 = t126 * t107;
t87 = t126 * t102;
t84 = t101 * t131 - t94 * t102;
t83 = -t92 * t102 - t106 * t131;
t81 = t108 * t87 + (-t100 * t125 + t105 * t89) * t103;
t79 = t91 * t110 + t113 * t116;
t76 = -t125 * t95 + t87 * t135 + t94 * t89;
t73 = -t125 * t93 - t87 * t132 + t92 * t89;
t71 = t84 * t110 + t113 * t118;
t69 = t83 * t110 + t113 * t117;
t1 = [(-MDP(1) - t127) * g(3); t127 * (-g(3) * t108 + (-g(1) * t101 + g(2) * t106) * t103); (t136 * t111 + t119 * t114) * MDP(5) + (-g(1) * (t109 * t118 + t76 * t128) - g(2) * (t109 * t117 + t73 * t128) - g(3) * (t109 * t116 + t81 * t128)) * MDP(19) + (-g(1) * (t112 * t118 - t76 * t129) - g(2) * (t112 * t117 - t73 * t129) - g(3) * (t112 * t116 - t81 * t129)) * MDP(20) + (-t113 * MDP(12) + MDP(13) * t110) * (g(1) * t76 + g(2) * t73 + g(3) * t81) + (pkin(3) * MDP(6) + MDP(4)) * (t119 * t111 - t136 * t114); (-g(1) * t84 - g(2) * t83 - g(3) * t91) * MDP(6); (g(1) * t71 + g(2) * t69 + g(3) * t79) * MDP(13) + (-MDP(19) * t112 + MDP(20) * t109 - MDP(12)) * (g(1) * (-t110 * t118 + t84 * t113) + g(2) * (-t110 * t117 + t83 * t113) + g(3) * (-t110 * t116 + t91 * t113)); (-g(1) * (-t71 * t109 - t76 * t112) - g(2) * (-t69 * t109 - t73 * t112) - g(3) * (-t79 * t109 - t81 * t112)) * MDP(19) + (-g(1) * (t76 * t109 - t71 * t112) - g(2) * (t73 * t109 - t69 * t112) - g(3) * (t81 * t109 - t79 * t112)) * MDP(20);];
taug  = t1;
