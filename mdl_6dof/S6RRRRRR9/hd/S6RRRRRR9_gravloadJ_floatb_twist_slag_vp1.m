% Calculate Gravitation load on the joints for
% S6RRRRRR9
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
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:16:45
% EndTime: 2019-03-10 05:16:48
% DurationCPUTime: 2.81s
% Computational Cost: add. (1694->276), mult. (4393->413), div. (0->0), fcn. (5605->16), ass. (0->114)
t155 = cos(qJ(3));
t141 = cos(pkin(6));
t156 = cos(qJ(2));
t126 = t141 * t156;
t153 = sin(qJ(2));
t154 = sin(qJ(1));
t157 = cos(qJ(1));
t100 = -t126 * t157 + t154 * t153;
t138 = sin(pkin(7));
t139 = sin(pkin(6));
t115 = t139 * t138;
t140 = cos(pkin(7));
t177 = t100 * t140 + t157 * t115;
t125 = t141 * t153;
t60 = t125 * t157 + t154 * t156;
t82 = sin(qJ(3));
t30 = -t155 * t60 + t177 * t82;
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t116 = t140 * t139;
t90 = -t100 * t138 + t157 * t116;
t14 = t30 * t84 + t81 * t90;
t27 = t177 * t155 + t60 * t82;
t80 = sin(qJ(5));
t83 = cos(qJ(5));
t118 = t14 * t80 + t27 * t83;
t180 = t14 * t83 - t27 * t80;
t176 = t30 * t81 - t84 * t90;
t94 = t126 * t154 + t153 * t157;
t86 = t154 * t116 + t94 * t138;
t173 = t154 * t115 - t94 * t140;
t171 = t156 * t116 + t138 * t141;
t61 = -t125 * t154 + t156 * t157;
t32 = t61 * t155 + t173 * t82;
t16 = t32 * t84 + t81 * t86;
t120 = t139 * t153;
t45 = t155 * t120 + t171 * t82;
t59 = -t115 * t156 + t140 * t141;
t26 = t45 * t84 + t59 * t81;
t170 = g(1) * t16 - g(2) * t14 + g(3) * t26;
t15 = t32 * t81 - t84 * t86;
t25 = -t45 * t81 + t59 * t84;
t169 = -g(1) * t15 + g(2) * t176 + g(3) * t25;
t31 = -t173 * t155 + t61 * t82;
t44 = t120 * t82 - t171 * t155;
t168 = (-g(1) * t31 - g(2) * t27 - g(3) * t44) * t81;
t160 = t84 * pkin(4);
t159 = pkin(11) + rSges(5,3);
t158 = pkin(12) + rSges(6,3);
t152 = t30 * t80;
t151 = t32 * t80;
t150 = t45 * t80;
t79 = qJ(5) + qJ(6);
t76 = sin(t79);
t149 = t76 * t84;
t77 = cos(t79);
t148 = t77 * t84;
t147 = t80 * t84;
t146 = t83 * t84;
t75 = pkin(5) * t83 + pkin(4);
t145 = t84 * t75;
t105 = t153 * t115;
t122 = t156 * t139;
t144 = pkin(2) * t122 + pkin(10) * t105;
t121 = t139 * t154;
t143 = t157 * pkin(1) + pkin(9) * t121;
t142 = pkin(13) + pkin(12) + rSges(7,3);
t108 = t153 * t116;
t54 = -t108 * t82 + t122 * t155;
t137 = t54 * pkin(3) + t144;
t136 = t80 * pkin(5) + pkin(11);
t21 = t27 * pkin(3);
t135 = -pkin(11) * t30 - t21;
t23 = t31 * pkin(3);
t134 = t32 * pkin(11) - t23;
t43 = t44 * pkin(3);
t133 = t45 * pkin(11) - t43;
t132 = t138 * pkin(10);
t131 = t81 * t138;
t130 = t82 * t140;
t129 = t84 * t138;
t123 = t157 * t139;
t128 = -pkin(1) * t154 + pkin(9) * t123;
t124 = t140 * t155;
t119 = -rSges(5,1) * t84 + rSges(5,2) * t81;
t7 = -t16 * t80 + t31 * t83;
t117 = -t26 * t80 + t44 * t83;
t112 = t77 * rSges(7,1) - t76 * rSges(7,2) + t75;
t36 = -t100 * t155 - t130 * t60;
t55 = t100 * pkin(2);
t111 = t36 * pkin(3) + t132 * t60 - t55;
t38 = -t130 * t61 - t155 * t94;
t57 = t94 * pkin(2);
t110 = t38 * pkin(3) + t132 * t61 - t57;
t102 = t76 * rSges(7,1) + t77 * rSges(7,2) + t136;
t101 = rSges(4,3) * t138 + t132;
t5 = -t16 * t76 + t31 * t77;
t6 = t16 * t77 + t31 * t76;
t96 = m(7) * (g(1) * (t5 * rSges(7,1) - t6 * rSges(7,2)) + g(2) * ((t14 * t76 + t27 * t77) * rSges(7,1) + (t14 * t77 - t27 * t76) * rSges(7,2)) + g(3) * ((-t26 * t76 + t44 * t77) * rSges(7,1) + (-t26 * t77 - t44 * t76) * rSges(7,2)));
t91 = -t60 * pkin(2) + t90 * pkin(10) + t128;
t89 = t30 * pkin(3) + t91;
t88 = t61 * pkin(2) + t86 * pkin(10) + t143;
t87 = t32 * pkin(3) + t88;
t53 = t108 * t155 + t122 * t82;
t40 = t105 * t81 + t54 * t84;
t39 = -t105 * t84 + t54 * t81;
t37 = t124 * t61 - t82 * t94;
t35 = -t100 * t82 + t124 * t60;
t20 = t131 * t61 + t38 * t84;
t19 = -t129 * t61 + t38 * t81;
t18 = t131 * t60 + t36 * t84;
t17 = -t129 * t60 + t36 * t81;
t8 = t16 * t83 + t31 * t80;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t154 - rSges(2,2) * t157) + g(2) * (rSges(2,1) * t157 - rSges(2,2) * t154)) - m(3) * (g(1) * (-t60 * rSges(3,1) + rSges(3,2) * t100 + rSges(3,3) * t123 + t128) + g(2) * (t61 * rSges(3,1) - rSges(3,2) * t94 + rSges(3,3) * t121 + t143)) - m(4) * (g(1) * (t30 * rSges(4,1) + rSges(4,2) * t27 + rSges(4,3) * t90 + t91) + g(2) * (t32 * rSges(4,1) - t31 * rSges(4,2) + rSges(4,3) * t86 + t88)) - m(5) * (g(1) * (t14 * rSges(5,1) - rSges(5,2) * t176 - t159 * t27 + t89) + g(2) * (t16 * rSges(5,1) - t15 * rSges(5,2) + t159 * t31 + t87)) - m(6) * (g(1) * (t180 * rSges(6,1) - t118 * rSges(6,2) + t14 * pkin(4) - t27 * pkin(11) + t158 * t176 + t89) + g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t16 * pkin(4) + t31 * pkin(11) + t15 * t158 + t87)) - m(7) * (g(1) * (-t102 * t27 + t112 * t14 + t142 * t176 + t89) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t136 * t31 + t142 * t15 + t16 * t75 + t87)) -m(3) * (g(1) * (-rSges(3,1) * t94 - t61 * rSges(3,2)) + g(2) * (-rSges(3,1) * t100 - t60 * rSges(3,2)) + g(3) * (rSges(3,1) * t122 - rSges(3,2) * t120)) - m(4) * (g(1) * (t38 * rSges(4,1) - t37 * rSges(4,2) + t101 * t61 - t57) + g(2) * (t36 * rSges(4,1) - t35 * rSges(4,2) + t101 * t60 - t55) + g(3) * (t54 * rSges(4,1) - t53 * rSges(4,2) + rSges(4,3) * t105 + t144)) - m(5) * (g(1) * (t20 * rSges(5,1) - t19 * rSges(5,2) + t159 * t37 + t110) + g(2) * (t18 * rSges(5,1) - t17 * rSges(5,2) + t159 * t35 + t111) + g(3) * (rSges(5,1) * t40 - rSges(5,2) * t39 + t159 * t53 + t137)) - m(6) * (g(1) * (t20 * pkin(4) + t37 * pkin(11) + (t20 * t83 + t37 * t80) * rSges(6,1) + (-t20 * t80 + t37 * t83) * rSges(6,2) + t158 * t19 + t110) + g(2) * (t18 * pkin(4) + t35 * pkin(11) + (t18 * t83 + t35 * t80) * rSges(6,1) + (-t18 * t80 + t35 * t83) * rSges(6,2) + t158 * t17 + t111) + g(3) * (t40 * pkin(4) + t53 * pkin(11) + (t40 * t83 + t53 * t80) * rSges(6,1) + (-t40 * t80 + t53 * t83) * rSges(6,2) + t158 * t39 + t137)) - m(7) * (g(1) * (t102 * t37 + t112 * t20 + t142 * t19 + t110) + g(2) * (t102 * t35 + t112 * t18 + t142 * t17 + t111) + g(3) * (t102 * t53 + t112 * t40 + t142 * t39 + t137)) -m(4) * (g(1) * (-rSges(4,1) * t31 - rSges(4,2) * t32) + g(2) * (-rSges(4,1) * t27 + rSges(4,2) * t30) + g(3) * (-rSges(4,1) * t44 - rSges(4,2) * t45)) - m(5) * (g(1) * (t119 * t31 + t159 * t32 - t23) + g(2) * (t119 * t27 - t159 * t30 - t21) + g(3) * (t119 * t44 + t159 * t45 - t43)) + (-g(1) * (-t31 * t145 + pkin(5) * t151 + (-t148 * t31 + t32 * t76) * rSges(7,1) + (t149 * t31 + t32 * t77) * rSges(7,2) + t134) - g(2) * (-t27 * t145 - pkin(5) * t152 + (-t148 * t27 - t30 * t76) * rSges(7,1) + (t149 * t27 - t30 * t77) * rSges(7,2) + t135) - g(3) * (-t44 * t145 + pkin(5) * t150 + (-t148 * t44 + t45 * t76) * rSges(7,1) + (t149 * t44 + t45 * t77) * rSges(7,2) + t133) - t142 * t168) * m(7) + (-g(1) * (-t31 * t160 + (-t146 * t31 + t151) * rSges(6,1) + (t147 * t31 + t32 * t83) * rSges(6,2) + t134) - g(2) * (-t27 * t160 + (-t146 * t27 - t152) * rSges(6,1) + (t147 * t27 - t30 * t83) * rSges(6,2) + t135) - g(3) * (-t44 * t160 + (-t146 * t44 + t150) * rSges(6,1) + (t147 * t44 + t45 * t83) * rSges(6,2) + t133) - t158 * t168) * m(6), -m(5) * (g(1) * (-rSges(5,1) * t15 - rSges(5,2) * t16) + g(2) * (rSges(5,1) * t176 + rSges(5,2) * t14) + g(3) * (rSges(5,1) * t25 - rSges(5,2) * t26)) - m(6) * (t170 * t158 + t169 * (t83 * rSges(6,1) - t80 * rSges(6,2) + pkin(4))) - m(7) * (t169 * t112 + t170 * t142) -m(6) * (g(1) * (rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (t118 * rSges(6,1) + t180 * rSges(6,2)) + g(3) * (t117 * rSges(6,1) + (-t26 * t83 - t44 * t80) * rSges(6,2))) - t96 - m(7) * (g(1) * t7 + g(2) * t118 + g(3) * t117) * pkin(5), -t96];
taug  = t1(:);
