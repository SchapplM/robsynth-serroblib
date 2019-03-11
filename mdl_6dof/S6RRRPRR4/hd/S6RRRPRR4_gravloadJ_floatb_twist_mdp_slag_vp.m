% Calculate Gravitation load on the joints for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:17:46
% EndTime: 2019-03-09 18:17:47
% DurationCPUTime: 0.30s
% Computational Cost: add. (371->79), mult. (363->115), div. (0->0), fcn. (336->12), ass. (0->44)
t106 = cos(qJ(2));
t101 = qJ(2) + qJ(3);
t97 = sin(t101);
t98 = cos(t101);
t121 = t98 * pkin(3) + t97 * qJ(4);
t127 = t106 * pkin(2) + t121;
t126 = MDP(17) - MDP(20);
t105 = sin(qJ(1));
t107 = cos(qJ(1));
t113 = g(1) * t107 + g(2) * t105;
t76 = -g(3) * t98 + t113 * t97;
t125 = pkin(3) * t97;
t124 = g(3) * t97;
t120 = t105 * t98;
t100 = pkin(11) + qJ(5);
t96 = qJ(6) + t100;
t90 = sin(t96);
t91 = cos(t96);
t78 = t107 * t91 + t120 * t90;
t79 = t107 * t90 - t120 * t91;
t119 = t107 * t98;
t80 = t105 * t91 - t119 * t90;
t81 = t105 * t90 + t119 * t91;
t122 = (-g(1) * t80 + g(2) * t78 + t90 * t124) * MDP(34) + (g(1) * t81 - g(2) * t79 + t124 * t91) * MDP(35);
t102 = sin(pkin(11));
t118 = t105 * t102;
t103 = cos(pkin(11));
t117 = t105 * t103;
t116 = t107 * t102;
t115 = t107 * t103;
t104 = sin(qJ(2));
t114 = -pkin(2) * t104 - t125;
t111 = pkin(1) + t127;
t94 = sin(t100);
t95 = cos(t100);
t110 = t126 * (t113 * t98 + t124) + (MDP(18) * t103 - MDP(19) * t102 + MDP(27) * t95 - MDP(28) * t94 + MDP(34) * t91 - MDP(35) * t90 + MDP(16)) * t76;
t108 = -pkin(8) - pkin(7);
t88 = qJ(4) * t119;
t87 = qJ(4) * t120;
t85 = t105 * t94 + t119 * t95;
t84 = t105 * t95 - t119 * t94;
t83 = t107 * t94 - t120 * t95;
t82 = t107 * t95 + t120 * t94;
t1 = [t113 * MDP(3) + (-g(1) * (-t117 * t98 + t116) - g(2) * (t115 * t98 + t118)) * MDP(18) + (-g(1) * (t118 * t98 + t115) - g(2) * (-t116 * t98 + t117)) * MDP(19) + ((g(1) * t108 - g(2) * t111) * t107 + (g(1) * t111 + g(2) * t108) * t105) * MDP(21) + (-g(1) * t83 - g(2) * t85) * MDP(27) + (-g(1) * t82 - g(2) * t84) * MDP(28) + (-g(1) * t79 - g(2) * t81) * MDP(34) + (-g(1) * t78 - g(2) * t80) * MDP(35) + (-t104 * MDP(10) + MDP(16) * t98 + t106 * MDP(9) - t126 * t97 + MDP(2)) * (g(1) * t105 - g(2) * t107); (-g(3) * t106 + t104 * t113) * MDP(9) + (g(3) * t104 + t106 * t113) * MDP(10) + (-g(1) * (t107 * t114 + t88) - g(2) * (t105 * t114 + t87) - g(3) * t127) * MDP(21) + t110; (-g(1) * (-t107 * t125 + t88) - g(2) * (-t105 * t125 + t87) - g(3) * t121) * MDP(21) + t110; -t76 * MDP(21); (-g(1) * t84 + g(2) * t82 + t124 * t94) * MDP(27) + (g(1) * t85 - g(2) * t83 + t124 * t95) * MDP(28) + t122; t122;];
taug  = t1;
