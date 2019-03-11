% Calculate Gravitation load on the joints for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:50:27
% EndTime: 2019-03-09 07:50:33
% DurationCPUTime: 2.90s
% Computational Cost: add. (2463->243), mult. (6934->369), div. (0->0), fcn. (9078->18), ass. (0->110)
t151 = cos(qJ(4));
t139 = cos(pkin(8));
t137 = sin(pkin(6));
t140 = cos(pkin(7));
t118 = t140 * t137;
t136 = sin(pkin(7));
t153 = cos(qJ(1));
t138 = cos(pkin(14));
t141 = cos(pkin(6));
t121 = t141 * t138;
t135 = sin(pkin(14));
t150 = sin(qJ(1));
t71 = -t153 * t121 + t150 * t135;
t108 = -t153 * t118 + t71 * t136;
t80 = sin(pkin(8));
t147 = t108 * t80;
t117 = t137 * t136;
t152 = cos(qJ(3));
t113 = t152 * t117;
t133 = t71 * t140;
t149 = sin(qJ(3));
t119 = t141 * t135;
t72 = t153 * t119 + t150 * t138;
t54 = t153 * t113 + t152 * t133 + t72 * t149;
t165 = t139 * t54 - t147;
t112 = t149 * t117;
t53 = -t153 * t112 - t149 * t133 + t72 * t152;
t83 = sin(qJ(4));
t22 = -t53 * t151 + t165 * t83;
t38 = t108 * t139 + t54 * t80;
t82 = sin(qJ(5));
t85 = cos(qJ(5));
t6 = t22 * t85 - t38 * t82;
t81 = sin(qJ(6));
t170 = t6 * t81;
t84 = cos(qJ(6));
t169 = t6 * t84;
t168 = t22 * t82 + t38 * t85;
t164 = t53 * t83;
t106 = -t150 * t121 - t153 * t135;
t101 = t106 * t140;
t107 = -t150 * t119 + t153 * t138;
t91 = -t152 * t101 + t107 * t149 - t150 * t113;
t96 = -t106 * t136 + t150 * t118;
t86 = t96 * t139 + t91 * t80;
t162 = t91 * t139 - t96 * t80;
t104 = -t138 * t117 + t141 * t140;
t116 = t137 * t135;
t159 = t138 * t118 + t141 * t136;
t95 = t149 * t116 - t152 * t159;
t161 = t104 * t80 - t95 * t139;
t158 = pkin(11) * t80;
t157 = t85 * pkin(5);
t155 = rSges(6,3) + pkin(12);
t154 = rSges(7,3) + pkin(13);
t146 = t80 * t82;
t145 = t80 * t85;
t144 = t81 * t85;
t143 = t84 * t85;
t126 = t137 * t150;
t142 = t153 * pkin(1) + qJ(2) * t126;
t131 = t83 * t139;
t127 = t153 * t137;
t129 = -t150 * pkin(1) + qJ(2) * t127;
t128 = t139 * t151;
t125 = -rSges(6,1) * t85 + rSges(6,2) * t82;
t28 = -t53 * t131 - t54 * t151;
t49 = t54 * pkin(3);
t124 = t28 * pkin(4) + t53 * t158 - t49;
t56 = t149 * t101 + t107 * t152 + t150 * t112;
t30 = -t56 * t131 - t91 * t151;
t51 = t91 * pkin(3);
t123 = t30 * pkin(4) + t56 * t158 - t51;
t65 = t152 * t116 + t149 * t159;
t42 = -t65 * t131 - t95 * t151;
t64 = t95 * pkin(3);
t122 = t42 * pkin(4) + t65 * t158 - t64;
t115 = rSges(7,1) * t84 - rSges(7,2) * t81 + pkin(5);
t109 = -t72 * pkin(2) - pkin(10) * t108 + t129;
t102 = -t53 * pkin(3) - pkin(11) * t38 + t109;
t99 = t22 * pkin(4) + t102;
t98 = t107 * pkin(2) + pkin(10) * t96 + t142;
t88 = t56 * pkin(3) + pkin(11) * t86 + t98;
t24 = t56 * t151 - t162 * t83;
t87 = t24 * pkin(4) + t88;
t48 = t104 * t139 + t95 * t80;
t41 = t65 * t128 - t95 * t83;
t35 = t65 * t151 + t161 * t83;
t34 = -t151 * t161 + t65 * t83;
t33 = t34 * pkin(4);
t32 = t65 * t146 + t42 * t85;
t31 = -t65 * t145 + t42 * t82;
t29 = t56 * t128 - t91 * t83;
t27 = t53 * t128 - t54 * t83;
t23 = t151 * t162 + t56 * t83;
t21 = -t54 * t128 + t151 * t147 - t164;
t19 = t151 * t165 + t164;
t17 = t23 * pkin(4);
t15 = t19 * pkin(4);
t14 = t35 * t85 + t48 * t82;
t13 = -t35 * t82 + t48 * t85;
t12 = t56 * t146 + t30 * t85;
t11 = -t56 * t145 + t30 * t82;
t10 = t53 * t146 + t28 * t85;
t9 = -t53 * t145 + t28 * t82;
t8 = t24 * t85 + t86 * t82;
t7 = t24 * t82 - t86 * t85;
t2 = t23 * t81 + t8 * t84;
t1 = t23 * t84 - t8 * t81;
t3 = [-m(2) * (g(1) * (-t150 * rSges(2,1) - t153 * rSges(2,2)) + g(2) * (t153 * rSges(2,1) - t150 * rSges(2,2))) - m(3) * (g(1) * (-t72 * rSges(3,1) + t71 * rSges(3,2) + rSges(3,3) * t127 + t129) + g(2) * (t107 * rSges(3,1) + t106 * rSges(3,2) + rSges(3,3) * t126 + t142)) - m(4) * (g(1) * (-rSges(4,1) * t53 + t54 * rSges(4,2) - rSges(4,3) * t108 + t109) + g(2) * (t56 * rSges(4,1) - t91 * rSges(4,2) + t96 * rSges(4,3) + t98)) - m(5) * (g(1) * (t22 * rSges(5,1) - t21 * rSges(5,2) - rSges(5,3) * t38 + t102) + g(2) * (t24 * rSges(5,1) - t23 * rSges(5,2) + t86 * rSges(5,3) + t88)) - m(6) * (g(1) * (t6 * rSges(6,1) - rSges(6,2) * t168 + t155 * t21 + t99) + g(2) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t155 * t23 + t87)) - m(7) * (g(1) * (t6 * pkin(5) + t21 * pkin(12) + (t21 * t81 + t169) * rSges(7,1) + (t21 * t84 - t170) * rSges(7,2) + t154 * t168 + t99) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t8 * pkin(5) + t23 * pkin(12) + t154 * t7 + t87)) (-m(3) - m(4) - m(5) - m(6) - m(7)) * (g(1) * t126 - g(2) * t127 + g(3) * t141) -m(4) * (g(1) * (-t91 * rSges(4,1) - t56 * rSges(4,2)) + g(2) * (-t54 * rSges(4,1) - t53 * rSges(4,2)) + g(3) * (-t95 * rSges(4,1) - t65 * rSges(4,2))) - m(6) * (g(1) * (rSges(6,1) * t12 - rSges(6,2) * t11 + t155 * t29 + t123) + g(2) * (rSges(6,1) * t10 - rSges(6,2) * t9 + t155 * t27 + t124) + g(3) * (rSges(6,1) * t32 - rSges(6,2) * t31 + t155 * t41 + t122)) - m(7) * (g(1) * (t12 * pkin(5) + t29 * pkin(12) + (t12 * t84 + t29 * t81) * rSges(7,1) + (-t12 * t81 + t29 * t84) * rSges(7,2) + t154 * t11 + t123) + g(2) * (t10 * pkin(5) + t27 * pkin(12) + (t10 * t84 + t27 * t81) * rSges(7,1) + (-t10 * t81 + t27 * t84) * rSges(7,2) + t154 * t9 + t124) + g(3) * (t32 * pkin(5) + t41 * pkin(12) + (t32 * t84 + t41 * t81) * rSges(7,1) + (-t32 * t81 + t41 * t84) * rSges(7,2) + t154 * t31 + t122)) + (-g(1) * (rSges(5,1) * t30 - rSges(5,2) * t29 - t51) - g(2) * (rSges(5,1) * t28 - rSges(5,2) * t27 - t49) - g(3) * (rSges(5,1) * t42 - rSges(5,2) * t41 - t64) - (g(1) * t56 + g(2) * t53 + g(3) * t65) * t80 * (rSges(5,3) + pkin(11))) * m(5), -m(5) * (g(1) * (-rSges(5,1) * t23 - rSges(5,2) * t24) + g(2) * (-rSges(5,1) * t19 + rSges(5,2) * t22) + g(3) * (-rSges(5,1) * t34 - rSges(5,2) * t35)) - m(6) * (g(1) * (t125 * t23 + t155 * t24 - t17) + g(2) * (t125 * t19 - t155 * t22 - t15) + g(3) * (t125 * t34 + t155 * t35 - t33)) + (-g(1) * (-t23 * t157 - t17 + t24 * pkin(12) + (-t23 * t143 + t24 * t81) * rSges(7,1) + (t23 * t144 + t24 * t84) * rSges(7,2)) - g(2) * (-t19 * t157 - t15 - t22 * pkin(12) + (-t19 * t143 - t22 * t81) * rSges(7,1) + (t19 * t144 - t22 * t84) * rSges(7,2)) - g(3) * (-t34 * t157 - t33 + t35 * pkin(12) + (-t34 * t143 + t35 * t81) * rSges(7,1) + (t34 * t144 + t35 * t84) * rSges(7,2)) - (-g(1) * t23 - g(2) * t19 - g(3) * t34) * t82 * t154) * m(7), -m(6) * (g(1) * (-rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (rSges(6,1) * t168 + rSges(6,2) * t6) + g(3) * (rSges(6,1) * t13 - rSges(6,2) * t14)) - m(7) * (g(1) * (-t115 * t7 + t154 * t8) + (t115 * t13 + t154 * t14) * g(3) + (t115 * t168 - t154 * t6) * g(2)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * ((t19 * t84 + t170) * rSges(7,1) + (-t19 * t81 + t169) * rSges(7,2)) + g(3) * ((-t14 * t81 + t34 * t84) * rSges(7,1) + (-t14 * t84 - t34 * t81) * rSges(7,2)))];
taug  = t3(:);
