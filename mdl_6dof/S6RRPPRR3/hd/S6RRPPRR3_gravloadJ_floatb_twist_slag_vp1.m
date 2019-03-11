% Calculate Gravitation load on the joints for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:54:34
% EndTime: 2019-03-09 08:54:37
% DurationCPUTime: 1.17s
% Computational Cost: add. (787->175), mult. (1615->259), div. (0->0), fcn. (1998->14), ass. (0->77)
t65 = sin(qJ(1));
t67 = cos(qJ(2));
t101 = t65 * t67;
t61 = cos(pkin(6));
t64 = sin(qJ(2));
t68 = cos(qJ(1));
t97 = t68 * t64;
t36 = -t61 * t101 - t97;
t121 = t36 * pkin(2);
t58 = sin(pkin(11));
t93 = cos(pkin(11));
t39 = -t67 * t58 - t64 * t93;
t73 = -t64 * t58 + t67 * t93;
t69 = t61 * t73;
t22 = t39 * t68 - t65 * t69;
t95 = t39 * t61;
t23 = t65 * t95 + t68 * t73;
t60 = cos(pkin(12));
t52 = pkin(4) * t60 + pkin(3);
t62 = -pkin(9) - qJ(4);
t124 = t22 * t52 - t23 * t62 + t121;
t102 = t65 * t64;
t99 = t67 * t68;
t34 = -t61 * t99 + t102;
t19 = t65 * t39 + t68 * t69;
t59 = sin(pkin(6));
t106 = t59 * t68;
t18 = -t65 * t73 + t68 * t95;
t56 = pkin(12) + qJ(5);
t54 = sin(t56);
t55 = cos(t56);
t5 = -t54 * t106 - t18 * t55;
t63 = sin(qJ(6));
t66 = cos(qJ(6));
t123 = t19 * t66 + t5 * t63;
t122 = t19 * t63 - t5 * t66;
t94 = qJ(4) + rSges(5,3);
t31 = t73 * t59;
t120 = g(1) * t22 + g(3) * t31;
t119 = -g(2) * t19 - t120;
t118 = pkin(2) * t67;
t114 = t55 * pkin(5);
t112 = rSges(7,3) + pkin(10);
t109 = t55 * t63;
t108 = t55 * t66;
t107 = t59 * t65;
t53 = pkin(1) + t118;
t103 = t65 * t53;
t92 = -m(5) - m(6) - m(7);
t57 = sin(pkin(12));
t91 = t57 * t107;
t90 = t57 * t106;
t32 = t39 * t59;
t50 = t59 * t118;
t87 = t31 * t52 + t32 * t62 + t50;
t86 = g(2) * t112;
t33 = pkin(2) * t61 * t64 + (-pkin(8) - qJ(3)) * t59;
t85 = rSges(4,3) * t59 - t33;
t84 = -t55 * t106 + t18 * t54;
t49 = t68 * t53;
t83 = -t65 * t33 + t49;
t80 = t34 * pkin(2);
t79 = rSges(6,1) * t55 - rSges(6,2) * t54;
t78 = -t68 * t33 - t103;
t76 = t66 * rSges(7,1) - t63 * rSges(7,2) + pkin(5);
t75 = t18 * t62 + t19 * t52 - t80;
t74 = pkin(4) * t91 + t22 * t62 + t23 * t52 + t83;
t70 = pkin(4) * t90 + t18 * t52 - t19 * t62 + t78;
t37 = -t102 * t61 + t99;
t35 = -t61 * t97 - t101;
t25 = -t32 * t55 + t54 * t61;
t24 = t32 * t54 + t55 * t61;
t9 = t107 * t54 + t23 * t55;
t8 = -t107 * t55 + t23 * t54;
t2 = -t22 * t63 + t66 * t9;
t1 = -t22 * t66 - t63 * t9;
t3 = [-m(2) * (g(1) * (-t65 * rSges(2,1) - rSges(2,2) * t68) + g(2) * (rSges(2,1) * t68 - t65 * rSges(2,2))) - m(3) * (g(1) * (t35 * rSges(3,1) + t34 * rSges(3,2) - t65 * pkin(1)) + g(2) * (t37 * rSges(3,1) + t36 * rSges(3,2) + pkin(1) * t68) + (g(1) * t68 + g(2) * t65) * t59 * (rSges(3,3) + pkin(8))) - m(4) * (g(1) * (rSges(4,1) * t18 - t19 * rSges(4,2) + t68 * t85 - t103) + g(2) * (rSges(4,1) * t23 + rSges(4,2) * t22 + t65 * t85 + t49)) - m(5) * (g(1) * (t18 * pkin(3) + (t18 * t60 + t90) * rSges(5,1) + (t106 * t60 - t18 * t57) * rSges(5,2) + t94 * t19 + t78) + g(2) * (t23 * pkin(3) + (t23 * t60 + t91) * rSges(5,1) + (t107 * t60 - t23 * t57) * rSges(5,2) - t94 * t22 + t83)) - m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t84 + t19 * rSges(6,3) + t70) + g(2) * (rSges(6,1) * t9 - rSges(6,2) * t8 - rSges(6,3) * t22 + t74)) - m(7) * (g(1) * (t122 * rSges(7,1) + t123 * rSges(7,2) - t5 * pkin(5) + t112 * t84 + t70) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t9 + t112 * t8 + t74)) -m(3) * (g(1) * (rSges(3,1) * t36 - rSges(3,2) * t37) + g(2) * (-rSges(3,1) * t34 + rSges(3,2) * t35) + g(3) * (rSges(3,1) * t67 - rSges(3,2) * t64) * t59) - m(4) * (g(1) * (t22 * rSges(4,1) - rSges(4,2) * t23 + t121) + g(2) * (rSges(4,1) * t19 + rSges(4,2) * t18 - t80) + g(3) * (rSges(4,1) * t31 + rSges(4,2) * t32 + t50)) - m(6) * (g(1) * (rSges(6,3) * t23 + t22 * t79 + t124) + g(2) * (-rSges(6,3) * t18 + t19 * t79 + t75) + g(3) * (-rSges(6,3) * t32 + t31 * t79 + t87)) - m(7) * (g(1) * (t22 * t114 + (t108 * t22 + t23 * t63) * rSges(7,1) + (-t109 * t22 + t23 * t66) * rSges(7,2) + t124) + g(2) * (t19 * t114 + (t108 * t19 - t18 * t63) * rSges(7,1) + (-t109 * t19 - t18 * t66) * rSges(7,2) + t75) + g(3) * (t31 * t114 + (t108 * t31 - t32 * t63) * rSges(7,1) + (-t109 * t31 - t32 * t66) * rSges(7,2) + t87) + (t112 * t120 + t19 * t86) * t54) + (-g(1) * (t23 * t94 + t121) - g(2) * (-t18 * t94 - t80) - g(3) * (-t32 * t94 + t50) + t119 * (t60 * rSges(5,1) - t57 * rSges(5,2) + pkin(3))) * m(5) (-m(4) + t92) * (g(3) * t61 + (g(1) * t65 - g(2) * t68) * t59) t92 * t119, -m(6) * (g(1) * (-rSges(6,1) * t8 - rSges(6,2) * t9) + g(2) * (rSges(6,1) * t84 - rSges(6,2) * t5) + g(3) * (rSges(6,1) * t24 - rSges(6,2) * t25)) - m(7) * (g(1) * (t112 * t9 - t76 * t8) + t5 * t86 + g(2) * t76 * t84 + (t112 * t25 + t76 * t24) * g(3)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t123 * rSges(7,1) + t122 * rSges(7,2)) + g(3) * ((-t25 * t63 - t31 * t66) * rSges(7,1) + (-t25 * t66 + t31 * t63) * rSges(7,2)))];
taug  = t3(:);
