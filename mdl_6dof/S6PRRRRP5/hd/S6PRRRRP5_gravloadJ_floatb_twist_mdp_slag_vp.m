% Calculate Gravitation load on the joints for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:25:55
% EndTime: 2019-03-09 00:25:57
% DurationCPUTime: 0.81s
% Computational Cost: add. (661->135), mult. (1810->229), div. (0->0), fcn. (2295->14), ass. (0->70)
t194 = MDP(18) - MDP(26);
t183 = cos(pkin(12));
t185 = cos(pkin(6));
t169 = t185 * t183;
t180 = sin(pkin(12));
t190 = sin(qJ(2));
t192 = cos(qJ(2));
t153 = -t192 * t169 + t180 * t190;
t181 = sin(pkin(7));
t182 = sin(pkin(6));
t165 = t182 * t181;
t184 = cos(pkin(7));
t197 = t153 * t184 + t183 * t165;
t167 = t185 * t180;
t152 = t192 * t167 + t183 * t190;
t164 = t182 * t180;
t196 = t152 * t184 - t181 * t164;
t166 = t184 * t182;
t195 = t192 * t166 + t185 * t181;
t191 = cos(qJ(3));
t141 = sin(qJ(5));
t145 = cos(qJ(4));
t179 = t141 * t145;
t144 = cos(qJ(5));
t178 = t144 * t145;
t177 = pkin(5) * t141 + pkin(10);
t176 = pkin(9) * t181;
t142 = sin(qJ(4));
t175 = t142 * t181;
t143 = sin(qJ(3));
t174 = t143 * t184;
t173 = t145 * t181;
t172 = t184 * t191;
t171 = t192 * t182;
t170 = t182 * t190;
t131 = t190 * t169 + t180 * t192;
t109 = t131 * t191 - t197 * t143;
t146 = t153 * t181 - t183 * t166;
t100 = t109 * t142 - t146 * t145;
t132 = -t190 * t167 + t183 * t192;
t111 = t132 * t191 - t196 * t143;
t147 = t152 * t181 + t184 * t164;
t102 = t111 * t142 - t147 * t145;
t123 = t195 * t143 + t191 * t170;
t151 = -t192 * t165 + t185 * t184;
t112 = t123 * t142 - t151 * t145;
t162 = g(1) * t102 + g(2) * t100 + g(3) * t112;
t157 = t190 * t166;
t156 = t190 * t165;
t140 = -qJ(6) - pkin(11);
t139 = pkin(5) * t144 + pkin(4);
t129 = -t143 * t157 + t191 * t171;
t128 = t143 * t171 + t191 * t157;
t122 = t143 * t170 - t195 * t191;
t119 = t129 * t145 + t142 * t156;
t118 = t129 * t142 - t145 * t156;
t117 = -t132 * t174 - t152 * t191;
t116 = t132 * t172 - t152 * t143;
t115 = -t131 * t174 - t153 * t191;
t114 = t131 * t172 - t153 * t143;
t113 = t123 * t145 + t151 * t142;
t110 = t132 * t143 + t196 * t191;
t108 = t131 * t143 + t197 * t191;
t107 = t117 * t145 + t132 * t175;
t106 = t117 * t142 - t132 * t173;
t105 = t115 * t145 + t131 * t175;
t104 = t115 * t142 - t131 * t173;
t103 = t111 * t145 + t147 * t142;
t101 = t109 * t145 + t146 * t142;
t1 = [(-MDP(1) - MDP(27)) * g(3); (g(1) * t152 + g(2) * t153 - g(3) * t171) * MDP(3) + (g(1) * t132 + g(2) * t131 + g(3) * t170) * MDP(4) + (-g(1) * t117 - g(2) * t115 - g(3) * t129) * MDP(10) + (g(1) * t116 + g(2) * t114 + g(3) * t128) * MDP(11) + (-g(1) * t107 - g(2) * t105 - g(3) * t119) * MDP(17) + (-g(1) * (t107 * t144 + t116 * t141) - g(2) * (t105 * t144 + t114 * t141) - g(3) * (t119 * t144 + t128 * t141)) * MDP(24) + (-g(1) * (-t107 * t141 + t116 * t144) - g(2) * (-t105 * t141 + t114 * t144) - g(3) * (-t119 * t141 + t128 * t144)) * MDP(25) + (-g(1) * (-t152 * pkin(2) + t117 * pkin(3) - t106 * t140 + t107 * t139 + t177 * t116 + t132 * t176) - g(2) * (-t153 * pkin(2) + t115 * pkin(3) - t104 * t140 + t105 * t139 + t177 * t114 + t131 * t176) - g(3) * (pkin(2) * t171 + t129 * pkin(3) + pkin(9) * t156 - t118 * t140 + t119 * t139 + t177 * t128)) * MDP(27) + t194 * (g(1) * t106 + g(2) * t104 + g(3) * t118); (-g(1) * (-t110 * t178 + t111 * t141) - g(2) * (-t108 * t178 + t109 * t141) - g(3) * (-t122 * t178 + t123 * t141)) * MDP(24) + (-g(1) * (t110 * t179 + t111 * t144) - g(2) * (t108 * t179 + t109 * t144) - g(3) * (t122 * t179 + t123 * t144)) * MDP(25) + (-t177 * MDP(27) + MDP(11)) * (g(1) * t111 + g(2) * t109 + g(3) * t123) + (-t194 * t142 + MDP(10) + t145 * MDP(17) + (t139 * t145 - t140 * t142 + pkin(3)) * MDP(27)) * (g(1) * t110 + g(2) * t108 + g(3) * t122); (-g(1) * (-t102 * t139 - t103 * t140) - g(2) * (-t100 * t139 - t101 * t140) - g(3) * (-t112 * t139 - t113 * t140)) * MDP(27) + t194 * (g(1) * t103 + g(2) * t101 + g(3) * t113) + (MDP(24) * t144 - MDP(25) * t141 + MDP(17)) * t162; (-g(1) * (-t103 * t144 - t110 * t141) - g(2) * (-t101 * t144 - t108 * t141) - g(3) * (-t113 * t144 - t122 * t141)) * MDP(25) + (pkin(5) * MDP(27) + MDP(24)) * (-g(1) * (-t103 * t141 + t110 * t144) - g(2) * (-t101 * t141 + t108 * t144) - g(3) * (-t113 * t141 + t122 * t144)); -t162 * MDP(27);];
taug  = t1;
