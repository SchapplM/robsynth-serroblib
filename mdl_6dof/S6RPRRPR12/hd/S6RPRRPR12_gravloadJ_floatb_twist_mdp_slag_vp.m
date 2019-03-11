% Calculate Gravitation load on the joints for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:52:59
% EndTime: 2019-03-09 05:53:02
% DurationCPUTime: 1.12s
% Computational Cost: add. (682->114), mult. (1865->187), div. (0->0), fcn. (2377->14), ass. (0->56)
t147 = cos(qJ(1));
t183 = sin(pkin(12));
t187 = cos(pkin(6));
t168 = t187 * t183;
t185 = cos(pkin(12));
t191 = sin(qJ(1));
t131 = t147 * t168 + t185 * t191;
t144 = sin(qJ(3));
t192 = cos(qJ(3));
t169 = t187 * t185;
t157 = -t147 * t169 + t183 * t191;
t141 = sin(pkin(6));
t184 = sin(pkin(7));
t177 = t141 * t184;
t186 = cos(pkin(7));
t202 = t147 * t177 + t157 * t186;
t118 = -t131 * t192 + t144 * t202;
t143 = sin(qJ(4));
t146 = cos(qJ(4));
t178 = t141 * t186;
t193 = -t147 * t178 + t157 * t184;
t109 = t118 * t143 + t146 * t193;
t115 = t131 * t144 + t192 * t202;
t142 = sin(qJ(6));
t145 = cos(qJ(6));
t206 = t109 * t142 - t115 * t145;
t205 = t109 * t145 + t115 * t142;
t110 = t118 * t146 - t143 * t193;
t201 = MDP(20) - MDP(23);
t153 = t147 * t183 + t169 * t191;
t197 = t153 * t186 - t177 * t191;
t196 = MDP(14) - MDP(22);
t194 = MDP(21) - MDP(24);
t182 = qJ(2) * t141;
t181 = t142 * t143;
t180 = t143 * t145;
t179 = pkin(1) * t147 + t182 * t191;
t174 = -pkin(1) * t191 + t147 * t182;
t167 = t186 * t185;
t166 = t184 * t187;
t163 = g(1) * t191 - g(2) * t147;
t132 = t147 * t185 - t168 * t191;
t120 = t132 * t192 - t144 * t197;
t148 = -t153 * t184 - t178 * t191;
t111 = t120 * t143 + t146 * t148;
t125 = t144 * t166 + (t144 * t167 + t183 * t192) * t141;
t152 = -t177 * t185 + t186 * t187;
t113 = t125 * t143 - t146 * t152;
t161 = g(1) * t111 - g(2) * t109 + g(3) * t113;
t124 = -t192 * t166 + (t144 * t183 - t167 * t192) * t141;
t119 = t132 * t144 + t192 * t197;
t114 = t125 * t146 + t143 * t152;
t112 = t120 * t146 - t143 * t148;
t106 = t111 * t142 + t119 * t145;
t105 = t111 * t145 - t119 * t142;
t1 = [t163 * MDP(2) + (g(1) * t131 - g(2) * t132) * MDP(4) + (-g(1) * t157 + g(2) * t153) * MDP(5) + (-g(1) * t174 - g(2) * t179) * MDP(7) + (-g(1) * t118 - g(2) * t120) * MDP(13) + (-g(1) * (-t131 * pkin(2) + t118 * pkin(3) + t110 * pkin(4) - pkin(10) * t115 + t109 * qJ(5) + t174) - g(2) * (t132 * pkin(2) + t120 * pkin(3) + t112 * pkin(4) + t119 * pkin(10) + t111 * qJ(5) + t179) + (g(1) * t193 + g(2) * t148) * pkin(9)) * MDP(25) + (-g(1) * t206 - g(2) * t106) * MDP(31) + (-g(1) * t205 - g(2) * t105) * MDP(32) + t194 * (g(1) * t109 + g(2) * t111) - t201 * (g(1) * t110 + g(2) * t112) + t196 * (-g(1) * t115 + g(2) * t119) + (MDP(6) * t141 - MDP(3)) * (-g(1) * t147 - g(2) * t191); (MDP(25) + MDP(7)) * (-g(3) * t187 - t141 * t163); (-g(1) * (-t119 * t181 + t120 * t145) - g(2) * (-t115 * t181 - t118 * t145) - g(3) * (-t124 * t181 + t125 * t145)) * MDP(31) + (-g(1) * (-t119 * t180 - t120 * t142) - g(2) * (-t115 * t180 + t118 * t142) - g(3) * (-t124 * t180 - t125 * t142)) * MDP(32) + (-MDP(25) * pkin(10) + t196) * (g(1) * t120 - g(2) * t118 + g(3) * t125) + (MDP(13) + MDP(25) * (pkin(4) * t146 + qJ(5) * t143 + pkin(3)) + t201 * t146 - t194 * t143) * (g(1) * t119 + g(2) * t115 + g(3) * t124); (-g(1) * (-pkin(4) * t111 + qJ(5) * t112) - g(2) * (pkin(4) * t109 - qJ(5) * t110) - g(3) * (-pkin(4) * t113 + qJ(5) * t114)) * MDP(25) + (-MDP(31) * t142 - MDP(32) * t145 + t194) * (g(1) * t112 - g(2) * t110 + g(3) * t114) + t201 * t161; -t161 * MDP(25); (-g(1) * t105 + g(2) * t205 - g(3) * (t113 * t145 - t124 * t142)) * MDP(31) + (g(1) * t106 - g(2) * t206 - g(3) * (-t113 * t142 - t124 * t145)) * MDP(32);];
taug  = t1;
