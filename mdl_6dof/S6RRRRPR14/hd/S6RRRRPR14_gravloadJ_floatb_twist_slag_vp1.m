% Calculate Gravitation load on the joints for
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR14_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR14_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:08:45
% EndTime: 2019-03-10 00:08:51
% DurationCPUTime: 2.55s
% Computational Cost: add. (1554->266), mult. (4034->395), div. (0->0), fcn. (5136->16), ass. (0->108)
t130 = sin(pkin(7));
t131 = sin(pkin(6));
t109 = t131 * t130;
t132 = cos(pkin(7));
t150 = cos(qJ(1));
t133 = cos(pkin(6));
t149 = cos(qJ(2));
t118 = t133 * t149;
t146 = sin(qJ(2));
t147 = sin(qJ(1));
t94 = -t118 * t150 + t147 * t146;
t166 = t150 * t109 + t94 * t132;
t148 = cos(qJ(3));
t117 = t133 * t146;
t53 = t117 * t150 + t147 * t149;
t77 = sin(qJ(3));
t23 = -t148 * t53 + t166 * t77;
t76 = sin(qJ(4));
t78 = cos(qJ(4));
t110 = t132 * t131;
t84 = -t150 * t110 + t94 * t130;
t7 = t23 * t78 - t76 * t84;
t6 = t23 * t76 + t78 * t84;
t89 = t118 * t147 + t146 * t150;
t79 = t147 * t110 + t89 * t130;
t20 = t148 * t166 + t53 * t77;
t165 = -t147 * t109 + t89 * t132;
t163 = t149 * t110 + t130 * t133;
t112 = t131 * t146;
t39 = t148 * t112 + t163 * t77;
t88 = -t109 * t149 + t132 * t133;
t19 = t39 * t78 + t76 * t88;
t54 = -t117 * t147 + t149 * t150;
t25 = t54 * t148 - t165 * t77;
t9 = t25 * t78 + t76 * t79;
t162 = g(1) * t9 - g(2) * t7 + g(3) * t19;
t18 = t39 * t76 - t78 * t88;
t8 = t25 * t76 - t78 * t79;
t161 = g(1) * t8 - g(2) * t6 + g(3) * t18;
t24 = t165 * t148 + t54 * t77;
t38 = t112 * t77 - t163 * t148;
t160 = (-g(1) * t24 - g(2) * t20 - g(3) * t38) * t76;
t152 = t78 * pkin(4);
t151 = pkin(11) + rSges(5,3);
t73 = sin(pkin(13));
t145 = t23 * t73;
t144 = t25 * t73;
t143 = t39 * t73;
t72 = pkin(13) + qJ(6);
t69 = sin(t72);
t142 = t69 * t78;
t70 = cos(t72);
t141 = t70 * t78;
t140 = t73 * t78;
t74 = cos(pkin(13));
t139 = t74 * t78;
t67 = pkin(5) * t74 + pkin(4);
t138 = t78 * t67;
t114 = t149 * t131;
t99 = t146 * t109;
t137 = pkin(2) * t114 + pkin(10) * t99;
t113 = t131 * t147;
t136 = t150 * pkin(1) + pkin(9) * t113;
t135 = pkin(12) + qJ(5) + rSges(7,3);
t134 = qJ(5) + rSges(6,3);
t102 = t146 * t110;
t47 = -t102 * t77 + t114 * t148;
t129 = t47 * pkin(3) + t137;
t128 = t73 * pkin(5) + pkin(11);
t14 = t20 * pkin(3);
t127 = -pkin(11) * t23 - t14;
t16 = t24 * pkin(3);
t126 = t25 * pkin(11) - t16;
t37 = t38 * pkin(3);
t125 = t39 * pkin(11) - t37;
t124 = t130 * pkin(10);
t123 = t76 * t130;
t122 = t77 * t132;
t121 = t78 * t130;
t115 = t150 * t131;
t120 = -pkin(1) * t147 + pkin(9) * t115;
t116 = t132 * t148;
t111 = -rSges(5,1) * t78 + rSges(5,2) * t76;
t106 = t70 * rSges(7,1) - t69 * rSges(7,2) + t67;
t29 = -t122 * t53 - t148 * t94;
t49 = t94 * pkin(2);
t105 = t29 * pkin(3) + t124 * t53 - t49;
t31 = -t122 * t54 - t148 * t89;
t51 = t89 * pkin(2);
t104 = t31 * pkin(3) + t124 * t54 - t51;
t96 = t69 * rSges(7,1) + t70 * rSges(7,2) + t128;
t95 = rSges(4,3) * t130 + t124;
t85 = -t53 * pkin(2) - t84 * pkin(10) + t120;
t82 = t23 * pkin(3) + t85;
t81 = t54 * pkin(2) + t79 * pkin(10) + t136;
t80 = t25 * pkin(3) + t81;
t46 = t102 * t148 + t114 * t77;
t33 = t47 * t78 + t76 * t99;
t32 = t47 * t76 - t78 * t99;
t30 = t116 * t54 - t77 * t89;
t28 = t116 * t53 - t77 * t94;
t13 = t123 * t54 + t31 * t78;
t12 = -t121 * t54 + t31 * t76;
t11 = t123 * t53 + t29 * t78;
t10 = -t121 * t53 + t29 * t76;
t3 = t24 * t69 + t70 * t9;
t2 = t24 * t70 - t69 * t9;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t147 - rSges(2,2) * t150) + g(2) * (rSges(2,1) * t150 - rSges(2,2) * t147)) - m(3) * (g(1) * (-t53 * rSges(3,1) + rSges(3,2) * t94 + rSges(3,3) * t115 + t120) + g(2) * (t54 * rSges(3,1) - rSges(3,2) * t89 + rSges(3,3) * t113 + t136)) - m(4) * (g(1) * (t23 * rSges(4,1) + rSges(4,2) * t20 - rSges(4,3) * t84 + t85) + g(2) * (t25 * rSges(4,1) - t24 * rSges(4,2) + rSges(4,3) * t79 + t81)) - m(5) * (g(1) * (t7 * rSges(5,1) - t6 * rSges(5,2) - t151 * t20 + t82) + g(2) * (t9 * rSges(5,1) - t8 * rSges(5,2) + t151 * t24 + t80)) - m(6) * (g(1) * (t7 * pkin(4) - t20 * pkin(11) + (-t20 * t73 + t7 * t74) * rSges(6,1) + (-t20 * t74 - t7 * t73) * rSges(6,2) + t134 * t6 + t82) + g(2) * (t24 * pkin(11) + t9 * pkin(4) + (t24 * t73 + t74 * t9) * rSges(6,1) + (t24 * t74 - t73 * t9) * rSges(6,2) + t134 * t8 + t80)) - m(7) * (g(1) * (t106 * t7 + t135 * t6 - t20 * t96 + t82) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t128 * t24 + t135 * t8 + t9 * t67 + t80)) -m(3) * (g(1) * (-rSges(3,1) * t89 - t54 * rSges(3,2)) + g(2) * (-rSges(3,1) * t94 - t53 * rSges(3,2)) + g(3) * (rSges(3,1) * t114 - rSges(3,2) * t112)) - m(4) * (g(1) * (t31 * rSges(4,1) - t30 * rSges(4,2) + t54 * t95 - t51) + g(2) * (t29 * rSges(4,1) - t28 * rSges(4,2) + t53 * t95 - t49) + g(3) * (t47 * rSges(4,1) - t46 * rSges(4,2) + rSges(4,3) * t99 + t137)) - m(5) * (g(1) * (t13 * rSges(5,1) - t12 * rSges(5,2) + t151 * t30 + t104) + g(2) * (t11 * rSges(5,1) - t10 * rSges(5,2) + t151 * t28 + t105) + g(3) * (rSges(5,1) * t33 - rSges(5,2) * t32 + t151 * t46 + t129)) - m(6) * (g(1) * (t13 * pkin(4) + t30 * pkin(11) + (t13 * t74 + t30 * t73) * rSges(6,1) + (-t13 * t73 + t30 * t74) * rSges(6,2) + t134 * t12 + t104) + g(2) * (t11 * pkin(4) + t28 * pkin(11) + (t11 * t74 + t28 * t73) * rSges(6,1) + (-t11 * t73 + t28 * t74) * rSges(6,2) + t134 * t10 + t105) + g(3) * (t33 * pkin(4) + t46 * pkin(11) + (t33 * t74 + t46 * t73) * rSges(6,1) + (-t33 * t73 + t46 * t74) * rSges(6,2) + t134 * t32 + t129)) - m(7) * (g(1) * (t106 * t13 + t12 * t135 + t30 * t96 + t104) + g(2) * (t10 * t135 + t106 * t11 + t28 * t96 + t105) + g(3) * (t106 * t33 + t135 * t32 + t46 * t96 + t129)) -m(4) * (g(1) * (-rSges(4,1) * t24 - rSges(4,2) * t25) + g(2) * (-rSges(4,1) * t20 + rSges(4,2) * t23) + g(3) * (-rSges(4,1) * t38 - rSges(4,2) * t39)) - m(5) * (g(1) * (t111 * t24 + t151 * t25 - t16) + g(2) * (t111 * t20 - t151 * t23 - t14) + g(3) * (t111 * t38 + t151 * t39 - t37)) + (-g(1) * (-t24 * t138 + pkin(5) * t144 + (-t141 * t24 + t25 * t69) * rSges(7,1) + (t142 * t24 + t25 * t70) * rSges(7,2) + t126) - g(2) * (-t20 * t138 - pkin(5) * t145 + (-t141 * t20 - t23 * t69) * rSges(7,1) + (t142 * t20 - t23 * t70) * rSges(7,2) + t127) - g(3) * (-t38 * t138 + pkin(5) * t143 + (-t141 * t38 + t39 * t69) * rSges(7,1) + (t142 * t38 + t39 * t70) * rSges(7,2) + t125) - t135 * t160) * m(7) + (-g(1) * (-t24 * t152 + (-t139 * t24 + t144) * rSges(6,1) + (t140 * t24 + t25 * t74) * rSges(6,2) + t126) - g(2) * (-t20 * t152 + (-t139 * t20 - t145) * rSges(6,1) + (t20 * t140 - t23 * t74) * rSges(6,2) + t127) - g(3) * (-t38 * t152 + (-t139 * t38 + t143) * rSges(6,1) + (t140 * t38 + t39 * t74) * rSges(6,2) + t125) - t134 * t160) * m(6), -m(5) * (g(1) * (-rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (rSges(5,1) * t6 + rSges(5,2) * t7) + g(3) * (-rSges(5,1) * t18 - rSges(5,2) * t19)) - m(6) * (t162 * t134 + t161 * (-t74 * rSges(6,1) + t73 * rSges(6,2) - pkin(4))) - m(7) * (-t161 * t106 + t162 * t135) (-m(6) - m(7)) * t161, -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * ((t20 * t70 + t69 * t7) * rSges(7,1) + (-t20 * t69 + t7 * t70) * rSges(7,2)) + g(3) * ((-t19 * t69 + t38 * t70) * rSges(7,1) + (-t19 * t70 - t38 * t69) * rSges(7,2)))];
taug  = t1(:);
