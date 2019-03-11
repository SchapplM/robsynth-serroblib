% Calculate Gravitation load on the joints for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR15_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR15_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:34:14
% EndTime: 2019-03-10 00:34:21
% DurationCPUTime: 2.11s
% Computational Cost: add. (1449->240), mult. (3907->346), div. (0->0), fcn. (4968->14), ass. (0->105)
t149 = cos(pkin(6));
t159 = cos(qJ(2));
t133 = t149 * t159;
t156 = sin(qJ(2));
t157 = sin(qJ(1));
t90 = cos(qJ(1));
t116 = -t133 * t90 + t157 * t156;
t146 = sin(pkin(7));
t147 = sin(pkin(6));
t124 = t147 * t146;
t148 = cos(pkin(7));
t176 = t116 * t148 + t90 * t124;
t158 = cos(qJ(3));
t132 = t149 * t156;
t69 = t132 * t90 + t157 * t159;
t87 = sin(qJ(3));
t34 = -t158 * t69 + t176 * t87;
t86 = sin(qJ(4));
t89 = cos(qJ(4));
t125 = t148 * t147;
t99 = t116 * t146 - t90 * t125;
t11 = t34 * t89 - t86 * t99;
t10 = t34 * t86 + t89 * t99;
t85 = sin(qJ(6));
t88 = cos(qJ(6));
t175 = -rSges(7,1) * t85 - rSges(7,2) * t88;
t103 = t133 * t157 + t156 * t90;
t91 = t103 * t146 + t157 * t125;
t172 = -pkin(4) * t89 - qJ(5) * t86;
t164 = pkin(11) + pkin(5);
t171 = rSges(7,1) * t88 - rSges(7,2) * t85 + t164;
t31 = t158 * t176 + t69 * t87;
t160 = pkin(12) + rSges(7,3);
t169 = t103 * t148 - t157 * t124;
t168 = t159 * t125 + t146 * t149;
t166 = t175 * t86;
t162 = pkin(11) + rSges(6,1);
t161 = pkin(11) + rSges(5,3);
t113 = t156 * t124;
t130 = t159 * t147;
t153 = pkin(2) * t130 + pkin(10) * t113;
t129 = t147 * t157;
t152 = t90 * pkin(1) + pkin(9) * t129;
t150 = qJ(5) + rSges(6,3);
t25 = t31 * pkin(3);
t145 = t172 * t31 - t25;
t70 = -t132 * t157 + t159 * t90;
t35 = t169 * t158 + t70 * t87;
t27 = t35 * pkin(3);
t144 = t172 * t35 - t27;
t128 = t147 * t156;
t52 = t128 * t87 - t168 * t158;
t51 = t52 * pkin(3);
t143 = t172 * t52 - t51;
t114 = t156 * t125;
t63 = -t114 * t87 + t130 * t158;
t142 = t63 * pkin(3) + t153;
t141 = t146 * pkin(10);
t140 = t86 * t146;
t139 = t87 * t148;
t138 = t89 * t146;
t137 = t90 * t147;
t45 = t113 * t86 + t63 * t89;
t136 = t45 * pkin(4) + t142;
t135 = -pkin(1) * t157 + pkin(9) * t137;
t131 = t148 * t158;
t127 = -rSges(5,1) * t89 + rSges(5,2) * t86;
t126 = rSges(6,2) * t89 - rSges(6,3) * t86;
t122 = qJ(5) - t175;
t40 = -t116 * t158 - t139 * t69;
t65 = t116 * pkin(2);
t118 = t40 * pkin(3) + t141 * t69 - t65;
t42 = -t103 * t158 - t139 * t70;
t67 = t103 * pkin(2);
t117 = t42 * pkin(3) + t141 * t70 - t67;
t17 = t140 * t69 + t40 * t89;
t110 = t17 * pkin(4) + t118;
t19 = t140 * t70 + t42 * t89;
t109 = t19 * pkin(4) + t117;
t108 = rSges(4,3) * t146 + t141;
t102 = -t124 * t159 + t148 * t149;
t97 = -t69 * pkin(2) - t99 * pkin(10) + t135;
t96 = t34 * pkin(3) + t97;
t95 = t70 * pkin(2) + t91 * pkin(10) + t152;
t94 = t11 * pkin(4) + t96;
t36 = t70 * t158 - t169 * t87;
t93 = t36 * pkin(3) + t95;
t13 = t36 * t89 + t86 * t91;
t92 = t13 * pkin(4) + t93;
t62 = t114 * t158 + t130 * t87;
t53 = t158 * t128 + t168 * t87;
t44 = -t113 * t89 + t63 * t86;
t41 = -t103 * t87 + t131 * t70;
t39 = -t116 * t87 + t131 * t69;
t30 = t102 * t86 + t53 * t89;
t29 = -t102 * t89 + t53 * t86;
t24 = t29 * pkin(4);
t18 = -t138 * t70 + t42 * t86;
t16 = -t138 * t69 + t40 * t86;
t12 = t36 * t86 - t89 * t91;
t6 = t12 * pkin(4);
t4 = t10 * pkin(4);
t3 = t12 * t85 + t35 * t88;
t2 = t12 * t88 - t35 * t85;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t157 - t90 * rSges(2,2)) + g(2) * (t90 * rSges(2,1) - rSges(2,2) * t157)) - m(3) * (g(1) * (-t69 * rSges(3,1) + rSges(3,2) * t116 + rSges(3,3) * t137 + t135) + g(2) * (t70 * rSges(3,1) - rSges(3,2) * t103 + rSges(3,3) * t129 + t152)) - m(4) * (g(1) * (t34 * rSges(4,1) + rSges(4,2) * t31 - rSges(4,3) * t99 + t97) + g(2) * (t36 * rSges(4,1) - t35 * rSges(4,2) + rSges(4,3) * t91 + t95)) - m(5) * (g(1) * (t11 * rSges(5,1) - t10 * rSges(5,2) - t161 * t31 + t96) + g(2) * (t13 * rSges(5,1) - t12 * rSges(5,2) + t161 * t35 + t93)) - m(6) * (g(1) * (-t11 * rSges(6,2) + t150 * t10 - t162 * t31 + t94) + g(2) * (-t13 * rSges(6,2) + t150 * t12 + t162 * t35 + t92)) - m(7) * (g(1) * (t10 * t122 + t11 * t160 - t171 * t31 + t94) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t12 * qJ(5) + t160 * t13 + t164 * t35 + t92)) -m(3) * (g(1) * (-rSges(3,1) * t103 - t70 * rSges(3,2)) + g(2) * (-rSges(3,1) * t116 - t69 * rSges(3,2)) + g(3) * (rSges(3,1) * t130 - rSges(3,2) * t128)) - m(4) * (g(1) * (t42 * rSges(4,1) - t41 * rSges(4,2) + t108 * t70 - t67) + g(2) * (t40 * rSges(4,1) - t39 * rSges(4,2) + t108 * t69 - t65) + g(3) * (t63 * rSges(4,1) - t62 * rSges(4,2) + rSges(4,3) * t113 + t153)) - m(5) * (g(1) * (t19 * rSges(5,1) - t18 * rSges(5,2) + t161 * t41 + t117) + g(2) * (t17 * rSges(5,1) - t16 * rSges(5,2) + t161 * t39 + t118) + g(3) * (rSges(5,1) * t45 - rSges(5,2) * t44 + t161 * t62 + t142)) - m(6) * (g(1) * (-t19 * rSges(6,2) + t150 * t18 + t162 * t41 + t109) + g(2) * (-t17 * rSges(6,2) + t150 * t16 + t162 * t39 + t110) + g(3) * (-rSges(6,2) * t45 + t150 * t44 + t162 * t62 + t136)) - m(7) * (g(1) * (t122 * t18 + t160 * t19 + t171 * t41 + t109) + g(2) * (t122 * t16 + t160 * t17 + t171 * t39 + t110) + g(3) * (t122 * t44 + t160 * t45 + t171 * t62 + t136)) -m(4) * (g(1) * (-rSges(4,1) * t35 - rSges(4,2) * t36) + g(2) * (-rSges(4,1) * t31 + rSges(4,2) * t34) + g(3) * (-rSges(4,1) * t52 - rSges(4,2) * t53)) - m(5) * (g(1) * (t127 * t35 + t161 * t36 - t27) + g(2) * (t127 * t31 - t161 * t34 - t25) + g(3) * (t127 * t52 + t161 * t53 - t51)) - m(6) * (g(1) * (t126 * t35 + t162 * t36 + t144) + g(2) * (t126 * t31 - t162 * t34 + t145) + g(3) * (t126 * t52 + t162 * t53 + t143)) + (-g(1) * (t166 * t35 + t171 * t36 + t144) - g(2) * (t166 * t31 - t171 * t34 + t145) - g(3) * (t166 * t52 + t171 * t53 + t143) - (-g(1) * t35 - g(2) * t31 - g(3) * t52) * t89 * t160) * m(7), -m(5) * (g(1) * (-rSges(5,1) * t12 - rSges(5,2) * t13) + g(2) * (rSges(5,1) * t10 + rSges(5,2) * t11) + g(3) * (-rSges(5,1) * t29 - rSges(5,2) * t30)) - m(6) * (g(1) * (rSges(6,2) * t12 + t13 * t150 - t6) + g(2) * (-rSges(6,2) * t10 - t11 * t150 + t4) + g(3) * (rSges(6,2) * t29 + t150 * t30 - t24)) + (-g(1) * (-t160 * t12 - t6) - g(2) * (t10 * t160 + t4) - g(3) * (-t160 * t29 - t24) - (g(1) * t13 - g(2) * t11 + g(3) * t30) * t122) * m(7) (-m(6) - m(7)) * (g(1) * t12 - g(2) * t10 + g(3) * t29) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * ((-t10 * t88 - t31 * t85) * rSges(7,1) + (t10 * t85 - t31 * t88) * rSges(7,2)) + g(3) * ((t29 * t88 - t52 * t85) * rSges(7,1) + (-t29 * t85 - t52 * t88) * rSges(7,2)))];
taug  = t1(:);
