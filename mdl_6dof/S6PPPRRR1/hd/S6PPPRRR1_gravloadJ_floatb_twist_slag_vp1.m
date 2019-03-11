% Calculate Gravitation load on the joints for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPPRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:10
% EndTime: 2019-03-08 18:39:12
% DurationCPUTime: 0.54s
% Computational Cost: add. (1155->110), mult. (3262->192), div. (0->0), fcn. (4281->18), ass. (0->76)
t70 = sin(pkin(13));
t71 = sin(pkin(12));
t57 = t71 * t70;
t76 = cos(pkin(13));
t77 = cos(pkin(12));
t65 = t77 * t76;
t80 = cos(pkin(6));
t49 = t80 * t65 - t57;
t79 = cos(pkin(7));
t46 = t49 * t79;
t58 = t71 * t76;
t63 = t77 * t70;
t50 = t80 * t63 + t58;
t73 = sin(pkin(7));
t74 = sin(pkin(6));
t62 = t74 * t73;
t75 = cos(pkin(14));
t54 = t75 * t62;
t69 = sin(pkin(14));
t37 = -t75 * t46 + t50 * t69 + t77 * t54;
t64 = t77 * t74;
t43 = -t49 * t73 - t79 * t64;
t72 = sin(pkin(8));
t78 = cos(pkin(8));
t89 = t37 * t78 - t43 * t72;
t51 = -t80 * t58 - t63;
t47 = t51 * t79;
t52 = -t80 * t57 + t65;
t38 = -t75 * t47 + t52 * t69 - t71 * t54;
t61 = t74 * t71;
t44 = -t51 * t73 + t79 * t61;
t88 = t38 * t78 - t44 * t72;
t55 = t79 * t76 * t74;
t60 = t74 * t70;
t42 = t69 * t60 + (-t73 * t80 - t55) * t75;
t48 = -t76 * t62 + t80 * t79;
t87 = t42 * t78 - t48 * t72;
t34 = cos(qJ(5));
t86 = t34 * pkin(5);
t85 = rSges(6,3) + pkin(10);
t84 = rSges(7,3) + pkin(11);
t83 = cos(qJ(4));
t30 = sin(qJ(6));
t82 = t30 * t34;
t33 = cos(qJ(6));
t81 = t33 * t34;
t68 = -m(4) - m(5) - m(6) - m(7);
t67 = -m(3) + t68;
t31 = sin(qJ(5));
t66 = -rSges(6,1) * t34 + rSges(6,2) * t31;
t59 = t73 * t69;
t56 = rSges(7,1) * t33 - rSges(7,2) * t30 + pkin(5);
t53 = t74 * t59;
t32 = sin(qJ(4));
t27 = t69 * t55 + t80 * t59 + t75 * t60;
t23 = t42 * t72 + t48 * t78;
t22 = t69 * t47 + t52 * t75 + t71 * t53;
t21 = t69 * t46 + t50 * t75 - t77 * t53;
t17 = t38 * t72 + t44 * t78;
t16 = t37 * t72 + t43 * t78;
t15 = t27 * t83 - t87 * t32;
t14 = t27 * t32 + t87 * t83;
t13 = t14 * pkin(4);
t12 = t22 * t83 - t88 * t32;
t11 = t22 * t32 + t88 * t83;
t10 = t21 * t83 - t89 * t32;
t9 = t21 * t32 + t89 * t83;
t8 = t15 * t34 + t23 * t31;
t7 = -t15 * t31 + t23 * t34;
t6 = t11 * pkin(4);
t5 = t9 * pkin(4);
t4 = t12 * t34 + t17 * t31;
t3 = -t12 * t31 + t17 * t34;
t2 = t10 * t34 + t16 * t31;
t1 = -t10 * t31 + t16 * t34;
t18 = [(-m(2) + t67) * g(3), t67 * (g(1) * t61 - g(2) * t64 + g(3) * t80) t68 * (g(1) * t44 + g(2) * t43 + g(3) * t48) -m(5) * (g(1) * (-t11 * rSges(5,1) - t12 * rSges(5,2)) + g(2) * (-t9 * rSges(5,1) - t10 * rSges(5,2)) + g(3) * (-t14 * rSges(5,1) - t15 * rSges(5,2))) - m(6) * (g(1) * (t66 * t11 + t85 * t12 - t6) + g(2) * (t85 * t10 + t66 * t9 - t5) + g(3) * (t66 * t14 + t85 * t15 - t13)) + (-g(1) * (-t11 * t86 - t6 + t12 * pkin(10) + (-t11 * t81 + t12 * t30) * rSges(7,1) + (t11 * t82 + t12 * t33) * rSges(7,2)) - g(2) * (-t9 * t86 - t5 + t10 * pkin(10) + (t10 * t30 - t9 * t81) * rSges(7,1) + (t10 * t33 + t9 * t82) * rSges(7,2)) - g(3) * (-t14 * t86 - t13 + t15 * pkin(10) + (-t14 * t81 + t15 * t30) * rSges(7,1) + (t14 * t82 + t15 * t33) * rSges(7,2)) - (-g(1) * t11 - g(2) * t9 - g(3) * t14) * t31 * t84) * m(7), -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (t1 * rSges(6,1) - t2 * rSges(6,2)) + g(3) * (t7 * rSges(6,1) - t8 * rSges(6,2))) - m(7) * (g(3) * (t56 * t7 + t84 * t8) + (t56 * t1 + t84 * t2) * g(2) + (t56 * t3 + t84 * t4) * g(1)) -m(7) * (g(1) * ((t11 * t33 - t4 * t30) * rSges(7,1) + (-t11 * t30 - t4 * t33) * rSges(7,2)) + g(2) * ((-t2 * t30 + t9 * t33) * rSges(7,1) + (-t2 * t33 - t9 * t30) * rSges(7,2)) + g(3) * ((t14 * t33 - t8 * t30) * rSges(7,1) + (-t14 * t30 - t8 * t33) * rSges(7,2)))];
taug  = t18(:);
