% Calculate Gravitation load on the joints for
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:55:46
% EndTime: 2019-03-09 22:55:49
% DurationCPUTime: 0.93s
% Computational Cost: add. (642->136), mult. (1025->222), div. (0->0), fcn. (1212->14), ass. (0->59)
t170 = MDP(24) - MDP(27);
t128 = sin(qJ(2));
t129 = sin(qJ(1));
t131 = cos(qJ(2));
t159 = cos(pkin(6));
t165 = cos(qJ(1));
t142 = t159 * t165;
t107 = t128 * t129 - t131 * t142;
t122 = pkin(12) + qJ(6);
t118 = sin(t122);
t119 = cos(t122);
t108 = t128 * t142 + t129 * t131;
t123 = qJ(3) + qJ(4);
t120 = sin(t123);
t121 = cos(t123);
t125 = sin(pkin(6));
t147 = t125 * t165;
t97 = t108 * t121 - t120 * t147;
t169 = -t107 * t119 + t118 * t97;
t168 = t107 * t118 + t119 * t97;
t145 = t129 * t159;
t110 = -t128 * t145 + t165 * t131;
t127 = sin(qJ(3));
t130 = cos(qJ(3));
t149 = t125 * t130;
t102 = -t110 * t127 + t129 * t149;
t136 = t108 * t127 + t130 * t147;
t151 = t125 * t128;
t167 = g(2) * t136 - g(3) * (-t127 * t151 + t159 * t130) - g(1) * t102;
t109 = t165 * t128 + t131 * t145;
t166 = -g(1) * t109 - g(2) * t107;
t161 = g(2) * t129;
t160 = g(3) * t125;
t156 = t118 * t121;
t155 = t119 * t121;
t124 = sin(pkin(12));
t154 = t121 * t124;
t126 = cos(pkin(12));
t153 = t121 * t126;
t152 = t121 * t131;
t150 = t125 * t129;
t148 = t125 * t131;
t146 = t127 * t165;
t144 = t108 * t130 - t125 * t146;
t141 = g(1) * t110 + g(2) * t108;
t101 = t110 * t121 + t120 * t150;
t106 = t159 * t120 + t121 * t151;
t100 = t110 * t120 - t121 * t150;
t105 = t120 * t151 - t159 * t121;
t96 = t108 * t120 + t121 * t147;
t138 = g(1) * t100 + g(2) * t96 + g(3) * t105;
t139 = t170 * (g(1) * t101 + g(2) * t97 + g(3) * t106) + (MDP(25) * t126 - MDP(26) * t124 + MDP(34) * t119 - MDP(35) * t118 + MDP(23)) * t138;
t134 = -g(1) * (-t100 * pkin(4) + qJ(5) * t101) - g(2) * (-t96 * pkin(4) + t97 * qJ(5)) - g(3) * (-t105 * pkin(4) + t106 * qJ(5));
t132 = -pkin(10) - pkin(9);
t117 = pkin(3) * t130 + pkin(2);
t103 = t110 * t130 + t127 * t150;
t92 = t101 * t119 + t109 * t118;
t91 = -t101 * t118 + t109 * t119;
t1 = [(g(1) * t129 - g(2) * t165) * MDP(2) + (g(1) * t165 + t161) * MDP(3) + (g(1) * t108 - g(2) * t110) * MDP(9) + (-g(1) * t107 + g(2) * t109) * MDP(10) + (g(1) * t144 - g(2) * t103) * MDP(16) + (-g(1) * t136 - g(2) * t102) * MDP(17) + (g(1) * t97 - g(2) * t101) * MDP(23) + (-g(1) * (-t107 * t124 - t126 * t97) - g(2) * (t101 * t126 + t109 * t124)) * MDP(25) + (-g(1) * (-t107 * t126 + t124 * t97) - g(2) * (-t101 * t124 + t109 * t126)) * MDP(26) + (-g(1) * (-t129 * pkin(1) - pkin(4) * t97 - qJ(5) * t96 + t107 * t132 - t108 * t117) - g(2) * (t165 * pkin(1) + t101 * pkin(4) + t100 * qJ(5) - t109 * t132 + t110 * t117) + (-g(1) * (pkin(3) * t146 + t165 * pkin(8)) - (pkin(3) * t127 + pkin(8)) * t161) * t125) * MDP(28) + (g(1) * t168 - g(2) * t92) * MDP(34) + (-g(1) * t169 - g(2) * t91) * MDP(35) + t170 * (-g(1) * t96 + g(2) * t100); (g(3) * t151 + t141) * MDP(10) + (-g(1) * (-t109 * t153 + t110 * t124) - g(2) * (-t107 * t153 + t108 * t124) - (t124 * t128 + t126 * t152) * t160) * MDP(25) + (-g(1) * (t109 * t154 + t110 * t126) - g(2) * (t107 * t154 + t108 * t126) - (-t124 * t152 + t126 * t128) * t160) * MDP(26) + ((t128 * t160 + t141) * t132 + (-t131 * t160 - t166) * (pkin(4) * t121 + qJ(5) * t120 + t117)) * MDP(28) + (-g(1) * (-t109 * t155 + t110 * t118) - g(2) * (-t107 * t155 + t108 * t118) - (t118 * t128 + t119 * t152) * t160) * MDP(34) + (-g(1) * (t109 * t156 + t110 * t119) - g(2) * (t107 * t156 + t108 * t119) - (-t118 * t152 + t119 * t128) * t160) * MDP(35) + (-t130 * MDP(16) + t127 * MDP(17) - t121 * MDP(23) + t120 * t170 - MDP(9)) * (g(3) * t148 + t166); t167 * MDP(16) + (g(1) * t103 + g(2) * t144 - g(3) * (-t159 * t127 - t128 * t149)) * MDP(17) + (t167 * pkin(3) + t134) * MDP(28) + t139; t134 * MDP(28) + t139; -t138 * MDP(28); (-g(1) * t91 + g(2) * t169 - g(3) * (-t106 * t118 - t119 * t148)) * MDP(34) + (g(1) * t92 + g(2) * t168 - g(3) * (-t106 * t119 + t118 * t148)) * MDP(35);];
taug  = t1;
