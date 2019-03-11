% Calculate Gravitation load on the joints for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:58:52
% EndTime: 2019-03-09 09:58:55
% DurationCPUTime: 0.78s
% Computational Cost: add. (349->154), mult. (557->201), div. (0->0), fcn. (534->8), ass. (0->71)
t84 = rSges(7,1) + pkin(5);
t81 = rSges(5,3) + pkin(8);
t43 = cos(qJ(1));
t40 = sin(qJ(1));
t87 = g(2) * t40;
t52 = g(1) * t43 + t87;
t63 = rSges(7,3) + qJ(6);
t91 = -m(6) - m(7);
t38 = sin(qJ(4));
t90 = pkin(4) * t38;
t89 = g(1) * t40;
t42 = cos(qJ(2));
t86 = g(3) * t42;
t41 = cos(qJ(4));
t85 = t41 * pkin(4);
t33 = t42 * pkin(2);
t82 = -rSges(7,2) - pkin(2);
t80 = -rSges(6,3) - pkin(2);
t79 = rSges(4,2) * t42;
t78 = t38 * t43;
t39 = sin(qJ(2));
t77 = t39 * t40;
t76 = t39 * t43;
t36 = qJ(4) + pkin(9);
t30 = cos(t36);
t75 = t40 * t30;
t74 = t40 * t38;
t73 = t40 * t41;
t72 = t40 * t42;
t71 = t41 * t43;
t70 = t42 * t43;
t37 = -qJ(5) - pkin(8);
t69 = rSges(7,2) - t37;
t68 = rSges(6,3) - t37;
t26 = pkin(4) * t78;
t58 = t39 * t73;
t67 = pkin(4) * t58 + t26;
t31 = t39 * qJ(3);
t66 = t31 + t33;
t65 = t43 * pkin(1) + t40 * pkin(7);
t64 = qJ(3) * t42;
t62 = -pkin(2) - t81;
t61 = pkin(4) * t74;
t60 = t38 * t76;
t59 = t39 * t71;
t28 = pkin(3) + t85;
t34 = t43 * pkin(7);
t57 = t43 * t28 + t37 * t72 + t34;
t56 = -pkin(1) - t31;
t55 = pkin(2) * t70 + t43 * t31 + t65;
t54 = g(1) * t62;
t53 = pkin(4) * t59 - t61;
t51 = rSges(3,1) * t42 - rSges(3,2) * t39;
t49 = rSges(5,1) * t38 + rSges(5,2) * t41;
t29 = sin(t36);
t48 = rSges(6,1) * t29 + rSges(6,2) * t30;
t47 = pkin(4) * t60 + t40 * t28 + t55;
t46 = (-qJ(3) - t90) * t39 - pkin(1);
t45 = t29 * t84 - t30 * t63;
t22 = t40 * t64;
t24 = t43 * t64;
t44 = g(1) * (t42 * t26 + t37 * t76 + t24) + g(2) * (t37 * t77 + t42 * t61 + t22) + g(3) * (t39 * t90 + t66);
t10 = -t39 * t74 + t71;
t9 = t58 + t78;
t8 = t60 + t73;
t7 = t59 - t74;
t4 = -t29 * t77 + t30 * t43;
t3 = t29 * t43 + t39 * t75;
t2 = t29 * t76 + t75;
t1 = t29 * t40 - t30 * t76;
t5 = [-m(2) * (g(1) * (-t40 * rSges(2,1) - rSges(2,2) * t43) + g(2) * (rSges(2,1) * t43 - t40 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t43 + t34) + g(2) * (rSges(3,1) * t70 - rSges(3,2) * t76 + t65) + (g(1) * (-pkin(1) - t51) + g(2) * rSges(3,3)) * t40) - m(4) * (g(1) * (rSges(4,1) * t43 + t34) + g(2) * (-rSges(4,2) * t70 + rSges(4,3) * t76 + t55) + (g(1) * (-rSges(4,3) * t39 - t33 + t56 + t79) + g(2) * rSges(4,1)) * t40) - m(5) * (g(1) * (t10 * rSges(5,1) - t9 * rSges(5,2) + pkin(3) * t43 + t34) + g(2) * (t8 * rSges(5,1) + t7 * rSges(5,2) + t81 * t70 + t55) + (g(2) * pkin(3) + g(1) * t56 + t42 * t54) * t40) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t57) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t68 * t70 + t47) + (t42 * t80 + t46) * t89) - m(7) * (g(1) * (t63 * t3 + t84 * t4 + t57) + g(2) * (t1 * t63 + t2 * t84 + t69 * t70 + t47) + (t42 * t82 + t46) * t89) -m(3) * (g(3) * t51 + t52 * (-rSges(3,1) * t39 - rSges(3,2) * t42)) - m(4) * (g(1) * (rSges(4,3) * t70 + t24) + g(2) * (rSges(4,3) * t72 + t22) + g(3) * (t66 - t79) + (g(3) * rSges(4,3) + t52 * (rSges(4,2) - pkin(2))) * t39) - m(5) * (g(1) * t24 + g(2) * t22 + g(3) * t66 + (g(3) * t81 + t52 * t49) * t42 + (g(3) * t49 + t43 * t54 + t62 * t87) * t39) - m(6) * ((g(3) * t68 + t52 * t48) * t42 + (g(3) * t48 + t52 * t80) * t39 + t44) - m(7) * ((g(3) * t45 + t52 * t82) * t39 + (g(3) * t69 + t52 * t45) * t42 + t44) (-m(4) - m(5) + t91) * (t39 * t52 - t86) -m(5) * (g(1) * (rSges(5,1) * t7 - rSges(5,2) * t8) + g(2) * (rSges(5,1) * t9 + rSges(5,2) * t10)) - m(6) * (g(1) * (-rSges(6,1) * t1 - rSges(6,2) * t2 + t53) + g(2) * (rSges(6,1) * t3 + rSges(6,2) * t4 + t67)) - m(7) * (g(1) * (-t1 * t84 + t2 * t63 + t53) + g(2) * (t3 * t84 - t4 * t63 + t67)) + (-m(5) * (-rSges(5,1) * t41 + rSges(5,2) * t38) - m(6) * (-rSges(6,1) * t30 + rSges(6,2) * t29 - t85) - m(7) * (-t63 * t29 - t84 * t30 - t85)) * t86, t91 * (g(3) * t39 + t42 * t52) -m(7) * (g(1) * t1 - g(2) * t3 + t30 * t86)];
taug  = t5(:);
