% Calculate Gravitation load on the joints for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:54:25
% EndTime: 2019-03-08 18:54:27
% DurationCPUTime: 0.47s
% Computational Cost: add. (467->79), mult. (1279->134), div. (0->0), fcn. (1627->14), ass. (0->52)
t154 = MDP(12) - MDP(20);
t137 = sin(pkin(12));
t138 = sin(pkin(11));
t124 = t138 * t137;
t141 = cos(pkin(12));
t142 = cos(pkin(11));
t131 = t142 * t141;
t144 = cos(pkin(6));
t115 = -t144 * t131 + t124;
t139 = sin(pkin(7));
t140 = sin(pkin(6));
t128 = t140 * t139;
t143 = cos(pkin(7));
t153 = t115 * t143 + t142 * t128;
t125 = t138 * t141;
t129 = t142 * t137;
t116 = t144 * t125 + t129;
t127 = t140 * t138;
t152 = t116 * t143 - t139 * t127;
t151 = t141 * t143 * t140 + t144 * t139;
t145 = cos(qJ(3));
t104 = sin(qJ(5));
t108 = cos(qJ(4));
t136 = t104 * t108;
t107 = cos(qJ(5));
t135 = t107 * t108;
t134 = MDP(2) + MDP(21);
t130 = t142 * t140;
t126 = t140 * t137;
t105 = sin(qJ(4));
t109 = t115 * t139 - t143 * t130;
t106 = sin(qJ(3));
t96 = t144 * t129 + t125;
t83 = -t153 * t106 + t96 * t145;
t78 = t83 * t105 - t109 * t108;
t110 = t116 * t139 + t143 * t127;
t97 = -t144 * t124 + t131;
t85 = -t152 * t106 + t97 * t145;
t80 = t85 * t105 - t110 * t108;
t114 = -t141 * t128 + t144 * t143;
t91 = t151 * t106 + t145 * t126;
t86 = t91 * t105 - t114 * t108;
t122 = g(1) * t80 + g(2) * t78 + g(3) * t86;
t103 = -qJ(6) - pkin(10);
t102 = t107 * pkin(5) + pkin(4);
t90 = t106 * t126 - t151 * t145;
t87 = t114 * t105 + t91 * t108;
t84 = t97 * t106 + t152 * t145;
t82 = t96 * t106 + t153 * t145;
t81 = t110 * t105 + t85 * t108;
t79 = t109 * t105 + t83 * t108;
t1 = [(-MDP(1) - t134) * g(3); t134 * (-g(1) * t127 + g(2) * t130 - g(3) * t144); (-g(1) * (t85 * t104 - t84 * t135) - g(2) * (t83 * t104 - t82 * t135) - g(3) * (t91 * t104 - t90 * t135)) * MDP(18) + (-g(1) * (t85 * t107 + t84 * t136) - g(2) * (t83 * t107 + t82 * t136) - g(3) * (t91 * t107 + t90 * t136)) * MDP(19) + (MDP(5) - (pkin(5) * t104 + pkin(9)) * MDP(21)) * (g(1) * t85 + g(2) * t83 + g(3) * t91) + (-t154 * t105 + MDP(4) + t108 * MDP(11) + (t102 * t108 - t103 * t105 + pkin(3)) * MDP(21)) * (g(1) * t84 + g(2) * t82 + g(3) * t90); (-g(1) * (-t80 * t102 - t81 * t103) - g(2) * (-t78 * t102 - t79 * t103) - g(3) * (-t86 * t102 - t87 * t103)) * MDP(21) + t154 * (g(1) * t81 + g(2) * t79 + g(3) * t87) + (MDP(18) * t107 - MDP(19) * t104 + MDP(11)) * t122; (-g(1) * (-t84 * t104 - t81 * t107) - g(2) * (-t82 * t104 - t79 * t107) - g(3) * (-t90 * t104 - t87 * t107)) * MDP(19) + (pkin(5) * MDP(21) + MDP(18)) * (-g(1) * (-t81 * t104 + t84 * t107) - g(2) * (-t79 * t104 + t82 * t107) - g(3) * (-t87 * t104 + t90 * t107)); -t122 * MDP(21);];
taug  = t1;
