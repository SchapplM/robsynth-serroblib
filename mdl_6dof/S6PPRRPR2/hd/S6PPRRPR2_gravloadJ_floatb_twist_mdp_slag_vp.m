% Calculate Gravitation load on the joints for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:50:55
% EndTime: 2019-03-08 18:50:56
% DurationCPUTime: 0.51s
% Computational Cost: add. (455->77), mult. (1266->131), div. (0->0), fcn. (1615->14), ass. (0->51)
t149 = MDP(11) - MDP(14);
t133 = sin(pkin(12));
t134 = sin(pkin(11));
t121 = t134 * t133;
t137 = cos(pkin(12));
t138 = cos(pkin(11));
t128 = t138 * t137;
t140 = cos(pkin(6));
t111 = -t140 * t128 + t121;
t135 = sin(pkin(7));
t136 = sin(pkin(6));
t125 = t136 * t135;
t139 = cos(pkin(7));
t148 = t111 * t139 + t138 * t125;
t122 = t134 * t137;
t126 = t138 * t133;
t112 = t140 * t122 + t126;
t124 = t136 * t134;
t147 = t112 * t139 - t135 * t124;
t146 = t137 * t139 * t136 + t140 * t135;
t145 = MDP(12) - MDP(15);
t141 = cos(qJ(3));
t101 = sin(qJ(6));
t102 = sin(qJ(4));
t132 = t101 * t102;
t104 = cos(qJ(6));
t131 = t102 * t104;
t130 = MDP(16) + MDP(2);
t127 = t138 * t136;
t123 = t136 * t133;
t105 = cos(qJ(4));
t106 = t111 * t135 - t139 * t127;
t103 = sin(qJ(3));
t95 = t140 * t126 + t122;
t82 = -t148 * t103 + t95 * t141;
t77 = t82 * t102 - t105 * t106;
t107 = t112 * t135 + t139 * t124;
t96 = -t140 * t121 + t128;
t84 = -t147 * t103 + t96 * t141;
t79 = t84 * t102 - t105 * t107;
t110 = -t137 * t125 + t140 * t139;
t90 = t146 * t103 + t141 * t123;
t85 = t90 * t102 - t105 * t110;
t119 = g(1) * t79 + g(2) * t77 + g(3) * t85;
t89 = t103 * t123 - t146 * t141;
t86 = t102 * t110 + t90 * t105;
t83 = t96 * t103 + t147 * t141;
t81 = t95 * t103 + t148 * t141;
t80 = t102 * t107 + t84 * t105;
t78 = t102 * t106 + t82 * t105;
t1 = [(-MDP(1) - t130) * g(3); t130 * (-g(1) * t124 + g(2) * t127 - g(3) * t140); (-g(1) * (t84 * t104 - t83 * t132) - g(2) * (t82 * t104 - t81 * t132) - g(3) * (t90 * t104 - t89 * t132)) * MDP(22) + (-g(1) * (-t84 * t101 - t83 * t131) - g(2) * (-t82 * t101 - t81 * t131) - g(3) * (-t90 * t101 - t89 * t131)) * MDP(23) + (-pkin(9) * MDP(16) - MDP(13) + MDP(5)) * (g(1) * t84 + g(2) * t82 + g(3) * t90) + (MDP(4) + t149 * t105 - t145 * t102 + MDP(16) * (pkin(4) * t105 + qJ(5) * t102 + pkin(3))) * (g(1) * t83 + g(2) * t81 + g(3) * t89); (-g(1) * (-t79 * pkin(4) + t80 * qJ(5)) - g(2) * (-t77 * pkin(4) + t78 * qJ(5)) - g(3) * (-t85 * pkin(4) + t86 * qJ(5))) * MDP(16) + (-MDP(22) * t101 - MDP(23) * t104 + t145) * (g(1) * t80 + g(2) * t78 + g(3) * t86) + t149 * t119; -t119 * MDP(16); (-g(1) * (-t83 * t101 + t79 * t104) - g(2) * (-t81 * t101 + t77 * t104) - g(3) * (-t89 * t101 + t85 * t104)) * MDP(22) + (-g(1) * (-t79 * t101 - t83 * t104) - g(2) * (-t77 * t101 - t81 * t104) - g(3) * (-t85 * t101 - t89 * t104)) * MDP(23);];
taug  = t1;
