% Calculate Gravitation load on the joints for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:23:13
% EndTime: 2019-03-09 22:23:14
% DurationCPUTime: 0.28s
% Computational Cost: add. (354->79), mult. (367->117), div. (0->0), fcn. (355->10), ass. (0->46)
t92 = sin(qJ(1));
t95 = cos(qJ(1));
t101 = g(1) * t95 + g(2) * t92;
t120 = MDP(10) - MDP(25);
t91 = sin(qJ(2));
t114 = g(3) * t91;
t89 = qJ(3) + qJ(4);
t86 = cos(t89);
t107 = t95 * t86;
t94 = cos(qJ(2));
t111 = t92 * t94;
t85 = sin(t89);
t67 = t111 * t85 + t107;
t108 = t95 * t85;
t69 = -t108 * t94 + t92 * t86;
t119 = -g(1) * t69 + g(2) * t67 + t85 * t114;
t71 = -g(3) * t94 + t101 * t91;
t90 = sin(qJ(3));
t79 = t90 * pkin(3) + pkin(4) * t85;
t112 = pkin(7) + t79;
t84 = pkin(11) + qJ(6) + t89;
t81 = sin(t84);
t110 = t95 * t81;
t82 = cos(t84);
t109 = t95 * t82;
t106 = t95 * t90;
t93 = cos(qJ(3));
t105 = t95 * t93;
t63 = t111 * t81 + t109;
t64 = -t111 * t82 + t110;
t65 = -t110 * t94 + t92 * t82;
t66 = t109 * t94 + t92 * t81;
t104 = (-g(1) * t65 + g(2) * t63 + t81 * t114) * MDP(32) + (g(1) * t66 - g(2) * t64 + t114 * t82) * MDP(33);
t80 = t93 * pkin(3) + pkin(4) * t86;
t68 = -t111 * t86 + t108;
t70 = t107 * t94 + t92 * t85;
t102 = t119 * MDP(23) + (g(1) * t70 - g(2) * t68 + t114 * t86) * MDP(24) + t104;
t78 = pkin(2) + t80;
t88 = -qJ(5) - pkin(9) - pkin(8);
t99 = t94 * t78 - t91 * t88;
t97 = pkin(1) + t99;
t76 = t105 * t94 + t92 * t90;
t75 = -t106 * t94 + t92 * t93;
t74 = -t111 * t93 + t106;
t73 = t111 * t90 + t105;
t1 = [t101 * MDP(3) + (-g(1) * t74 - g(2) * t76) * MDP(16) + (-g(1) * t73 - g(2) * t75) * MDP(17) + (-g(1) * t68 - g(2) * t70) * MDP(23) + (-g(1) * t67 - g(2) * t69) * MDP(24) + ((-g(1) * t112 - g(2) * t97) * t95 + (g(1) * t97 - g(2) * t112) * t92) * MDP(26) + (-g(1) * t64 - g(2) * t66) * MDP(32) + (-g(1) * t63 - g(2) * t65) * MDP(33) + (t94 * MDP(9) - t120 * t91 + MDP(2)) * (g(1) * t92 - g(2) * t95); (-g(3) * t99 + t101 * (t78 * t91 + t88 * t94)) * MDP(26) + t120 * (t101 * t94 + t114) + (t93 * MDP(16) - t90 * MDP(17) + t86 * MDP(23) - t85 * MDP(24) + t82 * MDP(32) - t81 * MDP(33) + MDP(9)) * t71; (-g(1) * t75 + g(2) * t73 + t114 * t90) * MDP(16) + (g(1) * t76 - g(2) * t74 + t114 * t93) * MDP(17) + (-g(1) * (-t95 * t94 * t79 + t92 * t80) - g(2) * (-t111 * t79 - t95 * t80) + t79 * t114) * MDP(26) + t102; t119 * pkin(4) * MDP(26) + t102; -t71 * MDP(26); t104;];
taug  = t1;
