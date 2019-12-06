% Calculate Gravitation load on the joints for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:36
% EndTime: 2019-12-05 17:00:38
% DurationCPUTime: 0.48s
% Computational Cost: add. (281->75), mult. (729->125), div. (0->0), fcn. (879->10), ass. (0->46)
t142 = MDP(11) - MDP(20);
t140 = MDP(17) + MDP(19);
t139 = MDP(18) - MDP(21);
t106 = sin(pkin(9));
t110 = sin(qJ(2));
t113 = cos(qJ(2));
t132 = cos(pkin(9));
t133 = cos(pkin(5));
t122 = t133 * t132;
t95 = t106 * t110 - t113 * t122;
t124 = t106 * t133;
t97 = t132 * t110 + t113 * t124;
t141 = g(1) * t97 + g(2) * t95;
t107 = sin(pkin(5));
t131 = t107 * t110;
t112 = cos(qJ(3));
t130 = t107 * t112;
t129 = t107 * t113;
t108 = sin(qJ(4));
t128 = t108 * t112;
t111 = cos(qJ(4));
t127 = t111 * t112;
t126 = t111 * t113;
t125 = t108 * t129;
t123 = t107 * t132;
t109 = sin(qJ(3));
t121 = pkin(3) * t112 + pkin(8) * t109 + pkin(2);
t96 = t106 * t113 + t110 * t122;
t84 = -t109 * t123 + t96 * t112;
t74 = t84 * t108 - t95 * t111;
t98 = -t110 * t124 + t132 * t113;
t86 = t106 * t107 * t109 + t98 * t112;
t76 = t86 * t108 - t97 * t111;
t100 = t133 * t109 + t110 * t130;
t87 = t100 * t108 + t107 * t126;
t119 = g(1) * t76 + g(2) * t74 + g(3) * t87;
t90 = (t108 * t110 + t112 * t126) * t107;
t89 = -t111 * t131 + t112 * t125;
t88 = t100 * t111 - t125;
t82 = t98 * t108 - t97 * t127;
t81 = -t98 * t111 - t97 * t128;
t80 = t96 * t108 - t95 * t127;
t79 = -t96 * t111 - t95 * t128;
t77 = t97 * t108 + t86 * t111;
t75 = t95 * t108 + t84 * t111;
t1 = [(-MDP(1) - MDP(22)) * g(3); (g(1) * t98 + g(2) * t96 + g(3) * t131) * MDP(4) + (-g(1) * (t82 * pkin(4) + t98 * pkin(7) + t81 * qJ(5)) - g(2) * (t80 * pkin(4) + t96 * pkin(7) + t79 * qJ(5)) + t141 * t121 + (-t90 * pkin(4) - t89 * qJ(5) - (pkin(7) * t110 + t121 * t113) * t107) * g(3)) * MDP(22) + t140 * (-g(1) * t82 - g(2) * t80 - g(3) * t90) + t139 * (g(1) * t81 + g(2) * t79 + g(3) * t89) + (-t112 * MDP(10) + t142 * t109 - MDP(3)) * (g(3) * t129 - t141); (-pkin(8) * MDP(22) + t142) * (g(1) * t86 + g(2) * t84 + g(3) * t100) + (-MDP(10) - t140 * t111 + t139 * t108 - MDP(22) * (pkin(4) * t111 + qJ(5) * t108 + pkin(3))) * (g(3) * (-t109 * t131 + t133 * t112) + g(2) * (-t96 * t109 - t112 * t123) + g(1) * (t106 * t130 - t98 * t109)); (-g(1) * (-t76 * pkin(4) + t77 * qJ(5)) - g(2) * (-t74 * pkin(4) + t75 * qJ(5)) - g(3) * (-t87 * pkin(4) + t88 * qJ(5))) * MDP(22) + t140 * t119 + t139 * (g(1) * t77 + g(2) * t75 + g(3) * t88); -t119 * MDP(22);];
taug = t1;
