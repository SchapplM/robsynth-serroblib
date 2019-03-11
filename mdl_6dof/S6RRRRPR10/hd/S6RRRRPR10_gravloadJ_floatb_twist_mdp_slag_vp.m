% Calculate Gravitation load on the joints for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:08:33
% EndTime: 2019-03-09 23:08:35
% DurationCPUTime: 0.91s
% Computational Cost: add. (543->112), mult. (920->179), div. (0->0), fcn. (1079->12), ass. (0->54)
t172 = MDP(23) - MDP(26);
t171 = MDP(10) - MDP(25);
t169 = MDP(24) - MDP(27);
t124 = sin(qJ(2));
t125 = sin(qJ(1));
t128 = cos(qJ(2));
t129 = cos(qJ(1));
t158 = cos(pkin(6));
t145 = t129 * t158;
t107 = t124 * t125 - t128 * t145;
t122 = sin(qJ(6));
t126 = cos(qJ(6));
t108 = t124 * t145 + t125 * t128;
t120 = qJ(3) + qJ(4);
t118 = sin(t120);
t119 = cos(t120);
t121 = sin(pkin(6));
t150 = t121 * t129;
t96 = t108 * t118 + t119 * t150;
t168 = t107 * t126 + t122 * t96;
t167 = -t107 * t122 + t126 * t96;
t146 = t125 * t158;
t110 = -t124 * t146 + t128 * t129;
t123 = sin(qJ(3));
t127 = cos(qJ(3));
t151 = t121 * t127;
t102 = -t110 * t123 + t125 * t151;
t137 = t108 * t123 + t127 * t150;
t153 = t121 * t124;
t166 = g(2) * t137 - g(3) * (-t123 * t153 + t127 * t158) - g(1) * t102;
t165 = g(1) * t129 + g(2) * t125;
t159 = g(3) * t121;
t155 = t118 * t122;
t154 = t118 * t126;
t152 = t121 * t125;
t149 = t122 * t128;
t148 = t126 * t128;
t97 = t108 * t119 - t118 * t150;
t144 = t108 * t127 - t123 * t150;
t140 = g(1) * t110 + g(2) * t108;
t101 = t110 * t119 + t118 * t152;
t106 = t118 * t158 + t119 * t153;
t100 = t110 * t118 - t119 * t152;
t105 = t118 * t153 - t119 * t158;
t136 = g(1) * t100 + g(2) * t96 + g(3) * t105;
t139 = (-MDP(34) * t122 - MDP(35) * t126 + t169) * (g(1) * t101 + g(2) * t97 + g(3) * t106) + t172 * t136;
t132 = -g(1) * (-t100 * pkin(4) + qJ(5) * t101) - g(2) * (-t96 * pkin(4) + qJ(5) * t97) - g(3) * (-t105 * pkin(4) + t106 * qJ(5));
t130 = -pkin(10) - pkin(9);
t117 = pkin(3) * t127 + pkin(2);
t109 = t129 * t124 + t128 * t146;
t103 = t110 * t127 + t123 * t152;
t93 = t100 * t122 + t109 * t126;
t92 = t100 * t126 - t109 * t122;
t1 = [(g(1) * t125 - g(2) * t129) * MDP(2) + t165 * MDP(3) + (g(1) * t108 - g(2) * t110) * MDP(9) + (g(1) * t144 - g(2) * t103) * MDP(16) + (-g(1) * t137 - g(2) * t102) * MDP(17) + (-g(1) * (-t125 * pkin(1) - pkin(4) * t97 - qJ(5) * t96 + t107 * t130 - t108 * t117) - g(2) * (t129 * pkin(1) + t101 * pkin(4) + t100 * qJ(5) - t109 * t130 + t110 * t117) - t165 * t121 * (pkin(3) * t123 + pkin(8))) * MDP(28) + (g(1) * t168 - g(2) * t93) * MDP(34) + (g(1) * t167 - g(2) * t92) * MDP(35) + t169 * (-g(1) * t96 + g(2) * t100) - t172 * (-g(1) * t97 + g(2) * t101) - t171 * (g(1) * t107 - g(2) * t109); (t124 * t159 + t140) * MDP(28) * t130 + (-g(1) * (-t109 * t155 + t110 * t126) - g(2) * (-t107 * t155 + t108 * t126) - (t118 * t149 + t124 * t126) * t159) * MDP(34) + (-g(1) * (-t109 * t154 - t110 * t122) - g(2) * (-t107 * t154 - t108 * t122) - (t118 * t148 - t122 * t124) * t159) * MDP(35) + t171 * (g(3) * t153 + t140) + (-t172 * t119 + t169 * t118 - MDP(9) - t127 * MDP(16) + t123 * MDP(17) + (-pkin(4) * t119 - qJ(5) * t118 - t117) * MDP(28)) * (-g(1) * t109 - g(2) * t107 + t128 * t159); t166 * MDP(16) + (g(1) * t103 + g(2) * t144 - g(3) * (-t123 * t158 - t124 * t151)) * MDP(17) + (t166 * pkin(3) + t132) * MDP(28) + t139; MDP(28) * t132 + t139; -t136 * MDP(28); (-g(1) * t92 - g(2) * t167 - g(3) * (t105 * t126 + t121 * t149)) * MDP(34) + (g(1) * t93 + g(2) * t168 - g(3) * (-t105 * t122 + t121 * t148)) * MDP(35);];
taug  = t1;
