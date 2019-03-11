% Calculate Gravitation load on the joints for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:43:59
% EndTime: 2019-03-08 21:44:01
% DurationCPUTime: 0.61s
% Computational Cost: add. (288->89), mult. (719->140), div. (0->0), fcn. (833->10), ass. (0->50)
t148 = -MDP(10) + MDP(13) - MDP(23);
t147 = MDP(11) - MDP(14);
t111 = cos(qJ(3));
t146 = pkin(3) * t111 + pkin(2);
t108 = sin(qJ(3));
t145 = -qJ(4) * t108 - t146;
t144 = MDP(24) * pkin(5) + MDP(21);
t105 = sin(pkin(6));
t141 = g(3) * t105;
t110 = cos(qJ(5));
t103 = pkin(5) * t110 + pkin(4);
t140 = pkin(8) + t103;
t138 = cos(pkin(6));
t137 = cos(pkin(10));
t109 = sin(qJ(2));
t135 = t105 * t109;
t134 = t105 * t111;
t112 = cos(qJ(2));
t133 = t105 * t112;
t107 = sin(qJ(5));
t132 = t107 * t108;
t131 = t108 * t110;
t130 = t108 * t112;
t129 = t110 * t112;
t128 = MDP(15) + MDP(24);
t104 = sin(pkin(10));
t120 = t138 * t137;
t88 = t104 * t109 - t112 * t120;
t127 = t145 * t88;
t122 = t104 * t138;
t90 = t109 * t137 + t112 * t122;
t126 = t145 * t90;
t121 = t105 * t137;
t89 = t104 * t112 + t109 * t120;
t77 = t89 * t108 + t111 * t121;
t78 = -t108 * t121 + t111 * t89;
t125 = -t77 * pkin(3) + qJ(4) * t78;
t91 = -t109 * t122 + t112 * t137;
t79 = -t104 * t134 + t108 * t91;
t80 = t104 * t105 * t108 + t111 * t91;
t124 = -t79 * pkin(3) + qJ(4) * t80;
t92 = t108 * t135 - t111 * t138;
t93 = t108 * t138 + t109 * t134;
t123 = -t92 * pkin(3) + qJ(4) * t93;
t119 = g(3) * (t105 * qJ(4) * t130 + pkin(8) * t135 + t146 * t133);
t106 = -qJ(6) - pkin(9);
t118 = pkin(5) * t132 - t106 * t111;
t117 = g(1) * t79 + g(2) * t77 + g(3) * t92;
t116 = g(1) * t80 + g(2) * t78 + g(3) * t93;
t1 = [(-MDP(1) - t128) * g(3); (-g(1) * (pkin(8) * t91 + t126) - g(2) * (pkin(8) * t89 + t127) - t119) * MDP(15) + (-g(1) * (t110 * t91 - t132 * t90) - g(2) * (t110 * t89 - t132 * t88) - (t107 * t130 + t109 * t110) * t141) * MDP(21) + (-g(1) * (-t107 * t91 - t131 * t90) - g(2) * (-t107 * t89 - t131 * t88) - (-t107 * t109 + t108 * t129) * t141) * MDP(22) + (-g(1) * (-t118 * t90 + t140 * t91 + t126) - g(2) * (-t118 * t88 + t140 * t89 + t127) - t119 - (t103 * t109 + t112 * t118) * t141) * MDP(24) + (MDP(4) - MDP(12)) * (g(1) * t91 + g(2) * t89 + g(3) * t135) + (t147 * t108 + t148 * t111 - MDP(3)) * (-g(1) * t90 - g(2) * t88 + g(3) * t133); (-g(1) * t124 - g(2) * t125 - g(3) * t123) * MDP(15) + (-g(1) * (t106 * t79 + t124) - g(2) * (t106 * t77 + t125) - g(3) * (t106 * t92 + t123)) * MDP(24) - t148 * t117 + (-t110 * MDP(22) - t144 * t107 + t147) * t116; -t128 * t117; (-g(1) * (-t107 * t79 - t110 * t90) - g(2) * (-t107 * t77 - t110 * t88) - g(3) * (t105 * t129 - t92 * t107)) * MDP(22) + t144 * (-g(1) * (-t107 * t90 + t110 * t79) - g(2) * (-t107 * t88 + t110 * t77) - g(3) * (t107 * t133 + t92 * t110)); -t116 * MDP(24);];
taug  = t1;
