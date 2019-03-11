% Calculate Gravitation load on the joints for
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:46:49
% EndTime: 2019-03-09 13:46:51
% DurationCPUTime: 0.68s
% Computational Cost: add. (494->111), mult. (1138->195), div. (0->0), fcn. (1449->14), ass. (0->58)
t121 = qJ(5) + qJ(6);
t119 = sin(t121);
t120 = cos(t121);
t126 = sin(qJ(4));
t130 = cos(qJ(4));
t122 = sin(pkin(12));
t127 = sin(qJ(2));
t131 = cos(qJ(2));
t156 = cos(pkin(12));
t114 = -t131 * t122 - t127 * t156;
t124 = cos(pkin(6));
t107 = t114 * t124;
t128 = sin(qJ(1));
t132 = cos(qJ(1));
t136 = -t127 * t122 + t131 * t156;
t140 = -t107 * t132 + t128 * t136;
t123 = sin(pkin(6));
t152 = t123 * t132;
t92 = -t126 * t152 + t130 * t140;
t133 = t124 * t136;
t97 = t128 * t114 + t132 * t133;
t168 = t119 * t92 + t120 * t97;
t167 = -t119 * t97 + t120 * t92;
t125 = sin(qJ(5));
t129 = cos(qJ(5));
t166 = t125 * t92 + t129 * t97;
t165 = -t125 * t97 + t129 * t92;
t162 = g(3) * t123;
t106 = t114 * t123;
t103 = -t106 * t130 + t124 * t126;
t105 = t136 * t123;
t100 = t114 * t132 - t128 * t133;
t139 = t128 * t107 + t132 * t136;
t153 = t123 * t128;
t95 = t126 * t153 + t130 * t139;
t87 = -t100 * t120 - t119 * t95;
t88 = -t100 * t119 + t120 * t95;
t161 = (-g(1) * t87 + g(2) * t168 - g(3) * (-t103 * t119 - t105 * t120)) * MDP(32) + (g(1) * t88 + g(2) * t167 - g(3) * (-t103 * t120 + t105 * t119)) * MDP(33);
t155 = t119 * t130;
t154 = t120 * t130;
t151 = t125 * t130;
t149 = t127 * t132;
t148 = t128 * t127;
t147 = t128 * t131;
t146 = t129 * t130;
t145 = t131 * t132;
t141 = g(1) * t128 - g(2) * t132;
t138 = t126 * t140 + t130 * t152;
t137 = t124 * t145 - t148;
t111 = -t124 * t147 - t149;
t118 = pkin(2) * t131 + pkin(1);
t112 = -t124 * t148 + t145;
t110 = -t124 * t149 - t147;
t108 = pkin(2) * t124 * t127 + (-pkin(8) - qJ(3)) * t123;
t94 = -t126 * t139 + t130 * t153;
t90 = -t100 * t125 + t129 * t95;
t89 = -t100 * t129 - t125 * t95;
t1 = [t141 * MDP(2) + (-g(1) * t110 - g(2) * t112) * MDP(9) + (g(1) * t137 - g(2) * t111) * MDP(10) + (-g(1) * (-t108 * t132 - t128 * t118) - g(2) * (-t128 * t108 + t118 * t132)) * MDP(12) + (g(1) * t92 - g(2) * t95) * MDP(18) + (-g(1) * t138 - g(2) * t94) * MDP(19) + (g(1) * t165 - g(2) * t90) * MDP(25) + (-g(1) * t166 - g(2) * t89) * MDP(26) + (g(1) * t167 - g(2) * t88) * MDP(32) + (-g(1) * t168 - g(2) * t87) * MDP(33) + (-t123 * MDP(11) + MDP(3)) * (g(1) * t132 + g(2) * t128); (g(1) * t112 - g(2) * t110 + t127 * t162) * MDP(10) + (-g(1) * (t100 * t146 + t125 * t139) - g(2) * (t125 * t140 + t97 * t146) - g(3) * (t105 * t146 - t106 * t125)) * MDP(25) + (-g(1) * (-t100 * t151 + t129 * t139) - g(2) * (t129 * t140 - t97 * t151) - g(3) * (-t105 * t151 - t106 * t129)) * MDP(26) + (-g(1) * (t100 * t154 + t119 * t139) - g(2) * (t119 * t140 + t97 * t154) - g(3) * (t105 * t154 - t106 * t119)) * MDP(32) + (-g(1) * (-t100 * t155 + t120 * t139) - g(2) * (t120 * t140 - t97 * t155) - g(3) * (-t105 * t155 - t106 * t120)) * MDP(33) + (-t130 * MDP(18) + MDP(19) * t126) * (g(1) * t100 + g(2) * t97 + g(3) * t105) + (MDP(12) * pkin(2) + MDP(9)) * (-g(1) * t111 - g(2) * t137 - t131 * t162); (-g(3) * t124 - t141 * t123) * MDP(12); (g(1) * t95 + g(2) * t92 + g(3) * t103) * MDP(19) + (-MDP(25) * t129 + MDP(26) * t125 - MDP(32) * t120 + MDP(33) * t119 - MDP(18)) * (g(1) * t94 - g(2) * t138 + g(3) * (t106 * t126 + t124 * t130)); (-g(1) * t89 + g(2) * t166 - g(3) * (-t103 * t125 - t105 * t129)) * MDP(25) + (g(1) * t90 + g(2) * t165 - g(3) * (-t103 * t129 + t105 * t125)) * MDP(26) + t161; t161;];
taug  = t1;
