% Calculate Gravitation load on the joints for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:14:12
% EndTime: 2019-03-09 12:14:15
% DurationCPUTime: 0.95s
% Computational Cost: add. (722->127), mult. (1829->199), div. (0->0), fcn. (2329->12), ass. (0->66)
t231 = MDP(19) - MDP(28);
t208 = MDP(25) + MDP(27);
t230 = MDP(26) - MDP(29);
t178 = cos(pkin(6));
t176 = sin(pkin(11));
t181 = sin(qJ(2));
t185 = cos(qJ(2));
t221 = cos(pkin(11));
t194 = -t185 * t176 - t181 * t221;
t160 = t194 * t178;
t182 = sin(qJ(1));
t186 = cos(qJ(1));
t193 = -t181 * t176 + t185 * t221;
t148 = -t160 * t186 + t182 * t193;
t180 = sin(qJ(4));
t184 = cos(qJ(4));
t177 = sin(pkin(6));
t216 = t177 * t186;
t134 = t148 * t184 - t180 * t216;
t159 = t193 * t178;
t147 = t186 * t159 + t182 * t194;
t179 = sin(qJ(5));
t183 = cos(qJ(5));
t121 = t134 * t179 + t147 * t183;
t122 = t134 * t183 - t147 * t179;
t211 = t182 * t185;
t213 = t181 * t186;
t165 = -t178 * t211 - t213;
t217 = t177 * t185;
t232 = -g(1) * t165 - g(3) * t217;
t150 = -t159 * t182 + t186 * t194;
t157 = t193 * t177;
t229 = -g(1) * t150 - g(2) * t147 - g(3) * t157;
t218 = t177 * t182;
t215 = t179 * t184;
t212 = t182 * t181;
t210 = t183 * t184;
t209 = t185 * t186;
t206 = t178 * t209;
t203 = -t148 * t180 - t184 * t216;
t161 = pkin(2) * t178 * t181 + (-pkin(8) - qJ(3)) * t177;
t175 = pkin(2) * t185 + pkin(1);
t202 = -t161 * t182 + t186 * t175;
t198 = g(1) * t182 - g(2) * t186;
t149 = -t182 * t160 - t186 * t193;
t197 = -t161 * t186 - t182 * t175;
t138 = -t149 * t184 + t180 * t218;
t125 = t138 * t179 + t150 * t183;
t158 = t194 * t177;
t153 = -t158 * t184 + t178 * t180;
t131 = t153 * t179 + t157 * t183;
t192 = g(1) * t125 + g(2) * t121 + g(3) * t131;
t168 = pkin(2) * t206;
t166 = -t178 * t212 + t209;
t164 = -t178 * t213 - t211;
t163 = -t206 + t212;
t140 = t157 * t210 - t158 * t179;
t139 = t157 * t215 + t158 * t183;
t137 = -t149 * t180 - t184 * t218;
t132 = t153 * t183 - t157 * t179;
t130 = -t149 * t179 + t150 * t210;
t129 = t149 * t183 + t150 * t215;
t128 = t147 * t210 + t148 * t179;
t127 = t147 * t215 - t148 * t183;
t126 = t138 * t183 - t150 * t179;
t1 = [t198 * MDP(2) + (-g(1) * t164 - g(2) * t166) * MDP(9) + (-g(1) * t163 - g(2) * t165) * MDP(10) + (-g(1) * t197 - g(2) * t202) * MDP(12) + (g(1) * t134 - g(2) * t138) * MDP(18) + (-g(1) * (-t148 * pkin(3) - pkin(4) * t134 - pkin(5) * t122 + t147 * pkin(9) + pkin(10) * t203 - qJ(6) * t121 + t197) - g(2) * (-pkin(3) * t149 + pkin(4) * t138 + pkin(5) * t126 - pkin(9) * t150 + pkin(10) * t137 + qJ(6) * t125 + t202)) * MDP(30) + t208 * (g(1) * t122 - g(2) * t126) + t230 * (-g(1) * t121 + g(2) * t125) + t231 * (g(1) * t203 + g(2) * t137) + (-MDP(11) * t177 + MDP(3)) * (g(1) * t186 + g(2) * t182); (g(2) * t163 + t232) * MDP(9) + (g(3) * t177 * t181 + g(1) * t166 - g(2) * t164) * MDP(10) + (-g(2) * t168 + (g(2) * t212 + t232) * pkin(2)) * MDP(12) + (-g(1) * (pkin(2) * t165 + t130 * pkin(5) - t149 * pkin(9) + t129 * qJ(6)) - g(2) * (-pkin(2) * t212 + t128 * pkin(5) + pkin(9) * t148 + t127 * qJ(6) + t168) - g(3) * (pkin(2) * t217 + t140 * pkin(5) - t158 * pkin(9) + t139 * qJ(6)) + t229 * (pkin(4) * t184 + pkin(10) * t180 + pkin(3))) * MDP(30) + t230 * (g(1) * t129 + g(2) * t127 + g(3) * t139) + t208 * (-g(1) * t130 - g(2) * t128 - g(3) * t140) + (t184 * MDP(18) - t231 * t180) * t229; (MDP(12) + MDP(30)) * (-g(3) * t178 - t198 * t177); (-pkin(10) * MDP(30) + t231) * (g(1) * t138 + g(2) * t134 + g(3) * t153) + (MDP(18) + MDP(30) * (pkin(5) * t183 + qJ(6) * t179 + pkin(4)) + t208 * t183 - t230 * t179) * (-g(3) * (t158 * t180 + t178 * t184) - g(2) * t203 + g(1) * t137); (-g(1) * (-pkin(5) * t125 + qJ(6) * t126) - g(2) * (-pkin(5) * t121 + qJ(6) * t122) - g(3) * (-pkin(5) * t131 + qJ(6) * t132)) * MDP(30) + t208 * t192 + t230 * (g(1) * t126 + g(2) * t122 + g(3) * t132); -t192 * MDP(30);];
taug  = t1;
