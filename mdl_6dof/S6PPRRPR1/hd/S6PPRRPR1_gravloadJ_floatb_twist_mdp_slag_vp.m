% Calculate Gravitation load on the joints for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:47:29
% EndTime: 2019-03-08 18:47:30
% DurationCPUTime: 0.48s
% Computational Cost: add. (532->90), mult. (1405->155), div. (0->0), fcn. (1797->16), ass. (0->55)
t139 = sin(pkin(12));
t140 = sin(pkin(11));
t125 = t140 * t139;
t143 = cos(pkin(12));
t144 = cos(pkin(11));
t132 = t144 * t143;
t146 = cos(pkin(6));
t115 = -t146 * t132 + t125;
t141 = sin(pkin(7));
t142 = sin(pkin(6));
t129 = t142 * t141;
t145 = cos(pkin(7));
t154 = t115 * t145 + t144 * t129;
t126 = t140 * t143;
t130 = t144 * t139;
t116 = t146 * t126 + t130;
t128 = t142 * t140;
t153 = t116 * t145 - t141 * t128;
t152 = t143 * t145 * t142 + t146 * t141;
t151 = MDP(12) - MDP(15);
t147 = cos(qJ(3));
t104 = pkin(13) + qJ(6);
t102 = sin(t104);
t109 = cos(qJ(4));
t138 = t102 * t109;
t103 = cos(t104);
t137 = t103 * t109;
t105 = sin(pkin(13));
t136 = t105 * t109;
t106 = cos(pkin(13));
t135 = t106 * t109;
t134 = MDP(16) + MDP(2);
t131 = t144 * t142;
t127 = t142 * t139;
t107 = sin(qJ(4));
t110 = t115 * t141 - t145 * t131;
t108 = sin(qJ(3));
t96 = t146 * t130 + t126;
t83 = -t154 * t108 + t96 * t147;
t78 = t107 * t83 - t110 * t109;
t111 = t116 * t141 + t145 * t128;
t97 = -t146 * t125 + t132;
t85 = -t153 * t108 + t97 * t147;
t80 = t107 * t85 - t111 * t109;
t114 = -t143 * t129 + t146 * t145;
t91 = t152 * t108 + t147 * t127;
t86 = t107 * t91 - t114 * t109;
t123 = g(1) * t80 + g(2) * t78 + g(3) * t86;
t90 = t108 * t127 - t152 * t147;
t87 = t114 * t107 + t91 * t109;
t84 = t108 * t97 + t153 * t147;
t82 = t96 * t108 + t154 * t147;
t81 = t111 * t107 + t85 * t109;
t79 = t110 * t107 + t83 * t109;
t1 = [(-MDP(1) - t134) * g(3); t134 * (-g(1) * t128 + g(2) * t131 - g(3) * t146); (-g(1) * (t85 * t105 - t84 * t135) - g(2) * (t83 * t105 - t82 * t135) - g(3) * (t91 * t105 - t90 * t135)) * MDP(13) + (-g(1) * (t85 * t106 + t84 * t136) - g(2) * (t83 * t106 + t82 * t136) - g(3) * (t91 * t106 + t90 * t136)) * MDP(14) + (-g(1) * (t85 * t102 - t84 * t137) - g(2) * (t83 * t102 - t82 * t137) - g(3) * (t91 * t102 - t90 * t137)) * MDP(22) + (-g(1) * (t85 * t103 + t84 * t138) - g(2) * (t83 * t103 + t82 * t138) - g(3) * (t91 * t103 + t90 * t138)) * MDP(23) + (-pkin(9) * MDP(16) + MDP(5)) * (g(1) * t85 + g(2) * t83 + g(3) * t91) + (-t151 * t107 + MDP(16) * (pkin(4) * t109 + qJ(5) * t107 + pkin(3)) + t109 * MDP(11) + MDP(4)) * (g(1) * t84 + g(2) * t82 + g(3) * t90); (-g(1) * (-pkin(4) * t80 + qJ(5) * t81) - g(2) * (-pkin(4) * t78 + qJ(5) * t79) - g(3) * (-pkin(4) * t86 + qJ(5) * t87)) * MDP(16) + t151 * (g(1) * t81 + g(2) * t79 + g(3) * t87) + (MDP(13) * t106 - MDP(14) * t105 + MDP(22) * t103 - MDP(23) * t102 + MDP(11)) * t123; -t123 * MDP(16); (-g(1) * (-t102 * t81 + t103 * t84) - g(2) * (-t102 * t79 + t103 * t82) - g(3) * (-t102 * t87 + t103 * t90)) * MDP(22) + (-g(1) * (-t102 * t84 - t103 * t81) - g(2) * (-t102 * t82 - t103 * t79) - g(3) * (-t102 * t90 - t103 * t87)) * MDP(23);];
taug  = t1;
