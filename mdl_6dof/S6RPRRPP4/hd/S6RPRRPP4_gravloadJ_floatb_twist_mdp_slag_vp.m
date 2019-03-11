% Calculate Gravitation load on the joints for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:40:47
% EndTime: 2019-03-09 04:40:48
% DurationCPUTime: 0.48s
% Computational Cost: add. (343->94), mult. (382->127), div. (0->0), fcn. (352->10), ass. (0->48)
t93 = pkin(9) + qJ(3);
t88 = sin(t93);
t97 = -qJ(5) - pkin(8);
t122 = t88 * t97;
t101 = cos(qJ(4));
t87 = pkin(4) * t101 + pkin(3);
t90 = cos(t93);
t110 = t90 * t87 - t122;
t127 = MDP(14) - MDP(22) - MDP(25);
t100 = sin(qJ(1));
t102 = cos(qJ(1));
t108 = g(1) * t102 + g(2) * t100;
t63 = -g(3) * t90 + t108 * t88;
t126 = g(3) * t88;
t121 = t90 * t97;
t94 = qJ(4) + pkin(10);
t89 = sin(t94);
t120 = t100 * t89;
t91 = cos(t94);
t119 = t100 * t91;
t99 = sin(qJ(4));
t118 = t100 * t99;
t117 = t102 * t90;
t116 = t102 * t91;
t115 = t102 * t99;
t114 = t100 * t101;
t113 = t101 * t102;
t112 = -MDP(23) - MDP(27);
t111 = t90 * t115;
t109 = pkin(5) * t91 + qJ(6) * t89;
t78 = g(1) * t100 - g(2) * t102;
t70 = t90 * t118 + t113;
t65 = t90 * t120 + t116;
t67 = t89 * t117 - t119;
t106 = g(1) * t67 + g(2) * t65 + t89 * t126;
t96 = cos(pkin(9));
t86 = pkin(2) * t96 + pkin(1);
t98 = -pkin(7) - qJ(2);
t105 = pkin(4) * t118 - t100 * t98 + t87 * t117 + (-t122 + t86) * t102;
t64 = t108 * t90 + t126;
t104 = pkin(4) * t115 - t102 * t98 + (-t86 - t110) * t100;
t83 = pkin(4) * t114;
t73 = t90 * t113 + t118;
t72 = -t111 + t114;
t71 = -t90 * t114 + t115;
t68 = t90 * t116 + t120;
t66 = -t102 * t89 + t90 * t119;
t1 = [(-g(1) * (-t100 * pkin(1) + t102 * qJ(2)) - g(2) * (pkin(1) * t102 + t100 * qJ(2))) * MDP(7) + (-g(1) * t71 - g(2) * t73) * MDP(20) + (-g(1) * t70 - g(2) * t72) * MDP(21) + (-g(1) * t104 - g(2) * t105) * MDP(23) + (g(1) * t66 - g(2) * t68) * MDP(24) + (g(1) * t65 - g(2) * t67) * MDP(26) + (-g(1) * (-t66 * pkin(5) - t65 * qJ(6) + t104) - g(2) * (t68 * pkin(5) + t67 * qJ(6) + t105)) * MDP(27) + (MDP(3) - MDP(6)) * t108 + (MDP(2) + t90 * MDP(13) + MDP(4) * t96 - MDP(5) * sin(pkin(9)) - t127 * t88) * t78; (-MDP(7) + t112) * t78; (-g(3) * t110 + t108 * (t87 * t88 + t121)) * MDP(23) + (-g(3) * (t109 * t90 + t110) + t108 * (t121 - (-t109 - t87) * t88)) * MDP(27) + t127 * t64 + (t101 * MDP(20) - t99 * MDP(21) + t91 * MDP(24) + t89 * MDP(26) + MDP(13)) * t63; (-g(1) * t72 + g(2) * t70 + t99 * t126) * MDP(20) + (g(1) * t73 - g(2) * t71 + t101 * t126) * MDP(21) + (-g(1) * t83 + (g(2) * t113 + t64 * t99) * pkin(4)) * MDP(23) + t106 * MDP(24) + (-g(1) * t68 - g(2) * t66 - t91 * t126) * MDP(26) + (-g(1) * (-pkin(4) * t111 - t67 * pkin(5) + t68 * qJ(6) + t83) - g(2) * (-t70 * pkin(4) - t65 * pkin(5) + t66 * qJ(6)) - (-pkin(4) * t99 - pkin(5) * t89 + qJ(6) * t91) * t126) * MDP(27); t112 * t63; -t106 * MDP(27);];
taug  = t1;
