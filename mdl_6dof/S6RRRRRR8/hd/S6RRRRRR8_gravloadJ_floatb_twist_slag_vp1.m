% Calculate Gravitation load on the joints for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:48:19
% EndTime: 2019-03-10 04:48:27
% DurationCPUTime: 2.42s
% Computational Cost: add. (1609->278), mult. (3776->418), div. (0->0), fcn. (4776->16), ass. (0->116)
t107 = sin(qJ(6));
t110 = cos(qJ(6));
t104 = qJ(4) + qJ(5);
t101 = sin(t104);
t102 = cos(t104);
t105 = sin(pkin(7));
t112 = cos(qJ(1));
t106 = sin(pkin(6));
t161 = cos(pkin(7));
t145 = t106 * t161;
t162 = cos(pkin(6));
t183 = cos(qJ(2));
t135 = t162 * t183;
t180 = sin(qJ(2));
t181 = sin(qJ(1));
t83 = -t112 * t135 + t180 * t181;
t200 = -t83 * t105 + t112 * t145;
t109 = sin(qJ(3));
t144 = t109 * t161;
t154 = t106 * t112;
t151 = t105 * t154;
t182 = cos(qJ(3));
t134 = t162 * t180;
t84 = t112 * t134 + t181 * t183;
t41 = -t109 * t151 - t144 * t83 + t182 * t84;
t14 = -t101 * t200 + t102 * t41;
t133 = t161 * t182;
t40 = t84 * t109 + t133 * t83 + t151 * t182;
t207 = t107 * t14 - t110 * t40;
t206 = -t107 * t40 - t110 * t14;
t108 = sin(qJ(4));
t111 = cos(qJ(4));
t205 = t108 * t41 + t111 * t200;
t169 = t108 * t200;
t204 = -t111 * t41 + t169;
t185 = pkin(13) + rSges(7,3);
t123 = t112 * t180 + t135 * t181;
t117 = t123 * t105 + t181 * t145;
t148 = t106 * t181;
t195 = -t105 * t148 + t123 * t161;
t85 = t112 * t183 - t134 * t181;
t45 = -t195 * t109 + t85 * t182;
t19 = -t45 * t108 + t117 * t111;
t146 = t105 * t162;
t65 = t109 * t146 + (t144 * t183 + t180 * t182) * t106;
t149 = t106 * t183;
t82 = -t105 * t149 + t161 * t162;
t199 = -t108 * t65 + t111 * t82;
t197 = t110 * rSges(7,1) - t107 * rSges(7,2) + pkin(5);
t196 = -t101 * t41 - t102 * t200;
t194 = g(1) * t85 + g(2) * t84;
t44 = t109 * t85 + t195 * t182;
t147 = t106 * t180;
t64 = t109 * t147 - t133 * t149 - t146 * t182;
t193 = g(1) * t44 + g(2) * t40 + g(3) * t64;
t192 = rSges(6,1) * t196 - t14 * rSges(6,2);
t186 = pkin(11) + rSges(5,3);
t17 = t101 * t45 - t102 * t117;
t18 = t101 * t117 + t45 * t102;
t184 = -t17 * rSges(6,1) - t18 * rSges(6,2);
t179 = t102 * pkin(5);
t178 = t105 * pkin(10);
t34 = -t101 * t65 + t102 * t82;
t35 = t101 * t82 + t102 * t65;
t177 = t34 * rSges(6,1) - t35 * rSges(6,2);
t100 = pkin(4) * t111 + pkin(3);
t113 = -pkin(12) - pkin(11);
t176 = -t40 * t100 - t41 * t113;
t175 = -t44 * t100 - t45 * t113;
t174 = -t64 * t100 - t65 * t113;
t140 = t105 * t147;
t173 = pkin(2) * t149 + pkin(10) * t140;
t163 = t112 * pkin(1) + pkin(9) * t148;
t160 = t101 * t105;
t159 = t102 * t105;
t158 = t102 * t107;
t157 = t102 * t110;
t156 = t105 * t108;
t155 = t105 * t111;
t50 = -t109 * t83 + t133 * t84;
t51 = -t144 * t84 - t182 * t83;
t78 = t83 * pkin(2);
t153 = t51 * t100 - t50 * t113 - t78;
t52 = -t109 * t123 + t133 * t85;
t53 = -t123 * t182 - t144 * t85;
t80 = t123 * pkin(2);
t152 = t53 * t100 - t52 * t113 - t80;
t143 = t205 * pkin(4);
t142 = t19 * pkin(4);
t141 = t199 * pkin(4);
t138 = -pkin(1) * t181 + pkin(9) * t154;
t130 = t108 * t140;
t132 = t161 * t180;
t75 = (t109 * t183 + t132 * t182) * t106;
t76 = (-t109 * t132 + t182 * t183) * t106;
t137 = pkin(4) * t130 + t76 * t100 - t75 * t113 + t173;
t131 = -rSges(6,1) * t102 + rSges(6,2) * t101;
t127 = t185 * t14 + t196 * t197;
t126 = -t197 * t17 + t185 * t18;
t125 = t185 * t35 + t197 * t34;
t124 = -t84 * pkin(2) + t200 * pkin(10) + t138;
t122 = pkin(4) * t169 - t100 * t41 + t113 * t40 + t124;
t120 = t194 * (pkin(4) * t108 + pkin(10)) * t105;
t118 = t85 * pkin(2) + t117 * pkin(10) + t163;
t115 = t117 * t108;
t116 = pkin(4) * t115 + t45 * t100 - t44 * t113 + t118;
t55 = t101 * t140 + t76 * t102;
t54 = t101 * t76 - t102 * t140;
t24 = t102 * t53 + t160 * t85;
t23 = t101 * t53 - t159 * t85;
t22 = t102 * t51 + t160 * t84;
t21 = t101 * t51 - t159 * t84;
t20 = t45 * t111 + t115;
t2 = t107 * t44 + t110 * t18;
t1 = -t107 * t18 + t110 * t44;
t3 = [-m(2) * (g(1) * (-rSges(2,1) * t181 - t112 * rSges(2,2)) + g(2) * (t112 * rSges(2,1) - rSges(2,2) * t181)) - m(3) * (g(1) * (-t84 * rSges(3,1) + t83 * rSges(3,2) + rSges(3,3) * t154 + t138) + g(2) * (t85 * rSges(3,1) - rSges(3,2) * t123 + rSges(3,3) * t148 + t163)) - m(4) * (g(1) * (-rSges(4,1) * t41 + rSges(4,2) * t40 + rSges(4,3) * t200 + t124) + g(2) * (t45 * rSges(4,1) - t44 * rSges(4,2) + rSges(4,3) * t117 + t118)) - m(5) * (g(1) * (t204 * rSges(5,1) + t205 * rSges(5,2) - t41 * pkin(3) - t186 * t40 + t124) + g(2) * (t20 * rSges(5,1) + t19 * rSges(5,2) + t45 * pkin(3) + t186 * t44 + t118)) - m(6) * (g(1) * (-rSges(6,1) * t14 - rSges(6,2) * t196 - rSges(6,3) * t40 + t122) + g(2) * (t18 * rSges(6,1) - t17 * rSges(6,2) + t44 * rSges(6,3) + t116)) - m(7) * (g(1) * (t206 * rSges(7,1) + t207 * rSges(7,2) - t14 * pkin(5) + t185 * t196 + t122) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t18 * pkin(5) + t17 * t185 + t116)) -m(3) * (g(1) * (-rSges(3,1) * t123 - t85 * rSges(3,2)) + g(2) * (-rSges(3,1) * t83 - rSges(3,2) * t84) + g(3) * (rSges(3,1) * t183 - rSges(3,2) * t180) * t106) - m(4) * (g(1) * (rSges(4,1) * t53 - rSges(4,2) * t52 - t80) + g(2) * (rSges(4,1) * t51 - rSges(4,2) * t50 - t78) + g(3) * (t76 * rSges(4,1) - t75 * rSges(4,2) + t173) + (rSges(4,3) * g(3) * t147 + t194 * (rSges(4,3) + pkin(10))) * t105) - m(5) * (g(1) * (t53 * pkin(3) - t80 + t85 * t178 + (t111 * t53 + t156 * t85) * rSges(5,1) + (-t108 * t53 + t155 * t85) * rSges(5,2) + t186 * t52) + g(2) * (t51 * pkin(3) - t78 + t84 * t178 + (t111 * t51 + t156 * t84) * rSges(5,1) + (-t108 * t51 + t155 * t84) * rSges(5,2) + t186 * t50) + g(3) * (t76 * pkin(3) + (t76 * t111 + t130) * rSges(5,1) + (-t76 * t108 + t111 * t140) * rSges(5,2) + t186 * t75 + t173)) - m(6) * (g(1) * (rSges(6,1) * t24 - rSges(6,2) * t23 + rSges(6,3) * t52 + t152) + g(2) * (rSges(6,1) * t22 - rSges(6,2) * t21 + rSges(6,3) * t50 + t153) + g(3) * (rSges(6,1) * t55 - rSges(6,2) * t54 + rSges(6,3) * t75 + t137) + t120) - m(7) * (g(1) * (t24 * pkin(5) + (t107 * t52 + t110 * t24) * rSges(7,1) + (-t107 * t24 + t110 * t52) * rSges(7,2) + t152 + t185 * t23) + g(2) * (t22 * pkin(5) + (t107 * t50 + t110 * t22) * rSges(7,1) + (-t107 * t22 + t110 * t50) * rSges(7,2) + t153 + t185 * t21) + g(3) * (t55 * pkin(5) + (t107 * t75 + t110 * t55) * rSges(7,1) + (-t107 * t55 + t110 * t75) * rSges(7,2) + t185 * t54 + t137) + t120) -m(4) * (g(1) * (-rSges(4,1) * t44 - rSges(4,2) * t45) + g(2) * (-rSges(4,1) * t40 - rSges(4,2) * t41) + g(3) * (-rSges(4,1) * t64 - rSges(4,2) * t65)) - m(5) * ((g(1) * t45 + g(2) * t41 + g(3) * t65) * t186 + t193 * (-t111 * rSges(5,1) + t108 * rSges(5,2) - pkin(3))) - m(6) * (g(1) * (rSges(6,3) * t45 + t131 * t44 + t175) + g(2) * (rSges(6,3) * t41 + t131 * t40 + t176) + g(3) * (rSges(6,3) * t65 + t131 * t64 + t174)) + (-g(1) * (-t44 * t179 + (t107 * t45 - t157 * t44) * rSges(7,1) + (t110 * t45 + t158 * t44) * rSges(7,2) + t175) - g(2) * (-t40 * t179 + (t107 * t41 - t157 * t40) * rSges(7,1) + (t110 * t41 + t158 * t40) * rSges(7,2) + t176) - g(3) * (-t64 * t179 + (t107 * t65 - t157 * t64) * rSges(7,1) + (t110 * t65 + t158 * t64) * rSges(7,2) + t174) + t193 * t101 * t185) * m(7), -m(5) * (g(1) * (rSges(5,1) * t19 - rSges(5,2) * t20) + g(2) * (-rSges(5,1) * t205 + t204 * rSges(5,2)) + g(3) * (t199 * rSges(5,1) + (-t108 * t82 - t111 * t65) * rSges(5,2))) - m(6) * (g(1) * (t142 + t184) + g(2) * (-t143 + t192) + g(3) * (t141 + t177)) - m(7) * (g(1) * (t126 + t142) + g(2) * (t127 - t143) + g(3) * (t125 + t141)) -m(6) * (g(1) * t184 + g(2) * t192 + g(3) * t177) - m(7) * (g(1) * t126 + g(2) * t127 + g(3) * t125) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t207 * rSges(7,1) + t206 * rSges(7,2)) + g(3) * ((-t107 * t35 + t110 * t64) * rSges(7,1) + (-t107 * t64 - t110 * t35) * rSges(7,2)))];
taug  = t3(:);
