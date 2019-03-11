% Calculate Gravitation load on the joints for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:44:44
% EndTime: 2019-03-09 05:44:48
% DurationCPUTime: 1.36s
% Computational Cost: add. (793->136), mult. (2076->226), div. (0->0), fcn. (2657->16), ass. (0->59)
t184 = sin(pkin(12));
t189 = cos(pkin(6));
t171 = t189 * t184;
t187 = cos(pkin(12));
t193 = sin(qJ(1));
t195 = cos(qJ(1));
t131 = t195 * t171 + t193 * t187;
t147 = sin(qJ(3));
t194 = cos(qJ(3));
t173 = t189 * t187;
t158 = -t195 * t173 + t193 * t184;
t185 = sin(pkin(7));
t186 = sin(pkin(6));
t169 = t186 * t185;
t188 = cos(pkin(7));
t203 = t158 * t188 + t195 * t169;
t118 = -t131 * t194 + t147 * t203;
t146 = sin(qJ(4));
t148 = cos(qJ(4));
t170 = t188 * t186;
t196 = t158 * t185 - t195 * t170;
t110 = t118 * t148 - t146 * t196;
t115 = t131 * t147 + t194 * t203;
t143 = pkin(13) + qJ(6);
t140 = sin(t143);
t141 = cos(t143);
t207 = t110 * t140 + t115 * t141;
t206 = t110 * t141 - t115 * t140;
t109 = t118 * t146 + t148 * t196;
t154 = t193 * t173 + t195 * t184;
t199 = t154 * t188 - t193 * t169;
t198 = t187 * t170 + t189 * t185;
t197 = MDP(21) - MDP(24);
t183 = t140 * t148;
t182 = t141 * t148;
t144 = sin(pkin(13));
t181 = t144 * t148;
t145 = cos(pkin(13));
t180 = t145 * t148;
t176 = t186 * t193;
t179 = t195 * pkin(1) + qJ(2) * t176;
t177 = t195 * t186;
t178 = -t193 * pkin(1) + qJ(2) * t177;
t168 = t186 * t184;
t132 = -t193 * t171 + t195 * t187;
t120 = t132 * t194 - t147 * t199;
t149 = -t154 * t185 - t193 * t170;
t111 = t120 * t146 + t149 * t148;
t125 = t147 * t198 + t194 * t168;
t153 = -t187 * t169 + t189 * t188;
t113 = t125 * t146 - t153 * t148;
t166 = g(1) * t111 - g(2) * t109 + g(3) * t113;
t124 = t147 * t168 - t194 * t198;
t119 = t132 * t147 + t194 * t199;
t114 = t125 * t148 + t153 * t146;
t112 = t120 * t148 - t149 * t146;
t105 = t112 * t141 + t119 * t140;
t104 = -t112 * t140 + t119 * t141;
t1 = [(g(1) * t193 - g(2) * t195) * MDP(2) + (g(1) * t195 + g(2) * t193) * MDP(3) + (g(1) * t131 - g(2) * t132) * MDP(4) + (-g(1) * t158 + g(2) * t154) * MDP(5) + (-g(1) * t177 - g(2) * t176) * MDP(6) + (-g(1) * t178 - g(2) * t179) * MDP(7) + (-g(1) * t118 - g(2) * t120) * MDP(13) + (-g(1) * t115 + g(2) * t119) * MDP(14) + (-g(1) * t110 - g(2) * t112) * MDP(20) + (-g(1) * (t110 * t145 - t115 * t144) - g(2) * (t112 * t145 + t119 * t144)) * MDP(22) + (-g(1) * (-t110 * t144 - t115 * t145) - g(2) * (-t112 * t144 + t119 * t145)) * MDP(23) + (-g(1) * (-t131 * pkin(2) + t118 * pkin(3) + t110 * pkin(4) - pkin(10) * t115 + t109 * qJ(5) + t178) - g(2) * (t132 * pkin(2) + t120 * pkin(3) + t112 * pkin(4) + t119 * pkin(10) + t111 * qJ(5) + t179) + (g(1) * t196 + g(2) * t149) * pkin(9)) * MDP(25) + (-g(1) * t206 - g(2) * t105) * MDP(31) + (g(1) * t207 - g(2) * t104) * MDP(32) + t197 * (g(1) * t109 + g(2) * t111); (MDP(25) + MDP(7)) * (-g(1) * t176 + g(2) * t177 - g(3) * t189); (-g(1) * (-t119 * t180 + t120 * t144) - g(2) * (-t115 * t180 - t118 * t144) - g(3) * (-t124 * t180 + t125 * t144)) * MDP(22) + (-g(1) * (t119 * t181 + t120 * t145) - g(2) * (t115 * t181 - t118 * t145) - g(3) * (t124 * t181 + t125 * t145)) * MDP(23) + (-g(1) * (-t119 * t182 + t120 * t140) - g(2) * (-t115 * t182 - t118 * t140) - g(3) * (-t124 * t182 + t125 * t140)) * MDP(31) + (-g(1) * (t119 * t183 + t120 * t141) - g(2) * (t115 * t183 - t118 * t141) - g(3) * (t124 * t183 + t125 * t141)) * MDP(32) + (-pkin(10) * MDP(25) + MDP(14)) * (g(1) * t120 - g(2) * t118 + g(3) * t125) + (MDP(25) * (pkin(4) * t148 + qJ(5) * t146 + pkin(3)) + t148 * MDP(20) + MDP(13) - t197 * t146) * (g(1) * t119 + g(2) * t115 + g(3) * t124); (-g(1) * (-pkin(4) * t111 + qJ(5) * t112) - g(2) * (pkin(4) * t109 - qJ(5) * t110) - g(3) * (-pkin(4) * t113 + qJ(5) * t114)) * MDP(25) + t197 * (g(1) * t112 - g(2) * t110 + g(3) * t114) + (MDP(22) * t145 - MDP(23) * t144 + MDP(31) * t141 - MDP(32) * t140 + MDP(20)) * t166; -t166 * MDP(25); (-g(1) * t104 - g(2) * t207 - g(3) * (-t114 * t140 + t124 * t141)) * MDP(31) + (g(1) * t105 - g(2) * t206 - g(3) * (-t114 * t141 - t124 * t140)) * MDP(32);];
taug  = t1;
