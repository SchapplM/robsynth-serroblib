% Calculate Gravitation load on the joints for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:16:24
% EndTime: 2019-03-08 22:16:26
% DurationCPUTime: 0.82s
% Computational Cost: add. (600->118), mult. (1111->176), div. (0->0), fcn. (1311->14), ass. (0->64)
t43 = sin(qJ(2));
t45 = cos(qJ(2));
t71 = cos(pkin(11));
t72 = cos(pkin(6));
t61 = t72 * t71;
t70 = sin(pkin(11));
t15 = t43 * t61 + t70 * t45;
t60 = t72 * t70;
t17 = -t43 * t60 + t71 * t45;
t95 = g(1) * t17 + g(2) * t15;
t14 = t70 * t43 - t45 * t61;
t16 = t71 * t43 + t45 * t60;
t94 = g(1) * t16 + g(2) * t14;
t42 = sin(qJ(3));
t44 = cos(qJ(3));
t39 = sin(pkin(6));
t66 = t39 * t70;
t10 = t17 * t42 - t44 * t66;
t77 = t39 * t43;
t18 = t42 * t77 - t72 * t44;
t67 = t39 * t71;
t8 = t15 * t42 + t44 * t67;
t93 = g(1) * t10 + g(2) * t8 + g(3) * t18;
t11 = t17 * t44 + t42 * t66;
t19 = t72 * t42 + t44 * t77;
t9 = t15 * t44 - t42 * t67;
t92 = g(1) * t11 + g(2) * t9 + g(3) * t19;
t37 = pkin(12) + qJ(5);
t34 = qJ(6) + t37;
t29 = sin(t34);
t30 = cos(t34);
t40 = cos(pkin(12));
t31 = t40 * pkin(4) + pkin(3);
t33 = cos(t37);
t55 = rSges(7,1) * t30 - rSges(7,2) * t29 + pkin(5) * t33 + t31;
t41 = -pkin(9) - qJ(4);
t74 = rSges(7,3) + pkin(10) - t41;
t91 = t74 * t42 + t55 * t44;
t32 = sin(t37);
t56 = rSges(6,1) * t33 - rSges(6,2) * t32 + t31;
t75 = rSges(6,3) - t41;
t90 = t75 * t42 + t56 * t44;
t38 = sin(pkin(12));
t59 = rSges(5,1) * t40 - rSges(5,2) * t38 + pkin(3);
t73 = rSges(5,3) + qJ(4);
t89 = t73 * t42 + t59 * t44;
t80 = g(3) * t39;
t79 = t38 * pkin(4);
t78 = rSges(4,3) + pkin(8);
t76 = t39 * t45;
t69 = -m(5) - m(6) - m(7);
t68 = g(3) * (pkin(2) * t76 + pkin(8) * t77);
t65 = t14 * t33 - t32 * t9;
t64 = rSges(4,1) * t44 - rSges(4,2) * t42;
t63 = t38 * rSges(5,1) + t40 * rSges(5,2);
t62 = -t11 * t32 + t16 * t33;
t57 = -t19 * t32 - t33 * t76;
t54 = t29 * rSges(7,1) + t30 * rSges(7,2) + pkin(5) * t32 + t79;
t52 = t32 * rSges(6,1) + t33 * rSges(6,2) + t79;
t12 = t14 * pkin(2);
t13 = t16 * pkin(2);
t50 = -g(1) * t13 - g(2) * t12 + t68;
t49 = m(7) * (g(1) * ((-t11 * t29 + t16 * t30) * rSges(7,1) + (-t11 * t30 - t16 * t29) * rSges(7,2)) + g(2) * ((t14 * t30 - t29 * t9) * rSges(7,1) + (-t14 * t29 - t30 * t9) * rSges(7,2)) + g(3) * ((-t19 * t29 - t30 * t76) * rSges(7,1) + (-t19 * t30 + t29 * t76) * rSges(7,2)));
t1 = [(-m(2) - m(3) - m(4) + t69) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t16 - rSges(3,2) * t17) + g(2) * (-rSges(3,1) * t14 - rSges(3,2) * t15) + (rSges(3,1) * t45 - rSges(3,2) * t43) * t80) - m(4) * (g(1) * (-t64 * t16 + t78 * t17 - t13) + g(2) * (-t64 * t14 + t78 * t15 - t12) + t68 + (rSges(4,3) * t43 + t64 * t45) * t80) - m(5) * ((t63 * t43 + t45 * t89) * t80 + t50 + t95 * (pkin(8) + t63) - t94 * t89) - m(6) * ((t52 * t43 + t45 * t90) * t80 + t50 + t95 * (pkin(8) + t52) - t94 * t90) - m(7) * ((t54 * t43 + t45 * t91) * t80 + t50 + t95 * (pkin(8) + t54) - t94 * t91) -m(4) * (g(1) * (-rSges(4,1) * t10 - rSges(4,2) * t11) + g(2) * (-rSges(4,1) * t8 - rSges(4,2) * t9) + g(3) * (-rSges(4,1) * t18 - rSges(4,2) * t19)) - m(5) * (-t59 * t93 + t73 * t92) - m(6) * (-t56 * t93 + t75 * t92) - m(7) * (-t93 * t55 + t74 * t92) t69 * t93, -m(6) * (g(1) * (t62 * rSges(6,1) + (-t11 * t33 - t16 * t32) * rSges(6,2)) + g(2) * (t65 * rSges(6,1) + (-t14 * t32 - t33 * t9) * rSges(6,2)) + g(3) * (t57 * rSges(6,1) + (-t19 * t33 + t32 * t76) * rSges(6,2))) - t49 - m(7) * (g(1) * t62 + g(2) * t65 + g(3) * t57) * pkin(5), -t49];
taug  = t1(:);
