% Calculate Gravitation load on the joints for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:18:22
% EndTime: 2019-03-08 21:18:23
% DurationCPUTime: 0.93s
% Computational Cost: add. (443->132), mult. (1005->188), div. (0->0), fcn. (1180->12), ass. (0->62)
t35 = sin(pkin(11));
t37 = cos(pkin(11));
t92 = -rSges(6,1) * t35 - rSges(6,2) * t37;
t41 = cos(qJ(3));
t39 = sin(qJ(3));
t64 = qJ(4) * t39;
t89 = -pkin(3) * t41 - t64;
t66 = rSges(7,3) + pkin(9) + qJ(5);
t40 = sin(qJ(2));
t42 = cos(qJ(2));
t60 = cos(pkin(10));
t61 = cos(pkin(6));
t50 = t61 * t60;
t59 = sin(pkin(10));
t16 = t59 * t40 - t42 * t50;
t49 = t61 * t59;
t18 = t60 * t40 + t42 * t49;
t88 = g(1) * t18 + g(2) * t16;
t87 = rSges(6,3) + qJ(5);
t36 = sin(pkin(6));
t71 = t36 * t40;
t21 = t61 * t39 + t41 * t71;
t17 = t40 * t50 + t59 * t42;
t55 = t36 * t60;
t6 = t17 * t41 - t39 * t55;
t19 = -t40 * t49 + t60 * t42;
t54 = t36 * t59;
t8 = t19 * t41 + t39 * t54;
t86 = g(1) * t8 + g(2) * t6 + g(3) * t21;
t85 = t92 * t39;
t34 = pkin(11) + qJ(6);
t32 = sin(t34);
t33 = cos(t34);
t45 = rSges(7,1) * t32 + rSges(7,2) * t33 + pkin(5) * t35;
t84 = t45 * t39 + t66 * t41;
t83 = t37 * rSges(6,1) - t35 * rSges(6,2) + pkin(4) + pkin(8);
t80 = -m(6) - m(7);
t75 = g(3) * t36;
t74 = rSges(5,1) + pkin(8);
t73 = rSges(4,3) + pkin(8);
t70 = t36 * t42;
t68 = t39 * t42;
t67 = t41 * t42;
t65 = pkin(2) * t70 + pkin(8) * t71;
t63 = rSges(5,3) + qJ(4);
t58 = -m(5) + t80;
t13 = t16 * pkin(2);
t57 = t89 * t16 - t13;
t14 = t18 * pkin(2);
t56 = t89 * t18 - t14;
t53 = rSges(4,1) * t41 - rSges(4,2) * t39;
t52 = rSges(5,2) * t41 - rSges(5,3) * t39;
t51 = g(3) * (t36 * pkin(3) * t67 + t64 * t70 + t65);
t48 = t33 * rSges(7,1) - t32 * rSges(7,2) + pkin(5) * t37 + pkin(4);
t46 = pkin(8) + t48;
t20 = t39 * t71 - t61 * t41;
t15 = t20 * pkin(3);
t7 = t19 * t39 - t41 * t54;
t5 = t17 * t39 + t41 * t55;
t4 = t7 * pkin(3);
t3 = t5 * pkin(3);
t1 = [(-m(2) - m(3) - m(4) + t58) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t18 - rSges(3,2) * t19) + g(2) * (-rSges(3,1) * t16 - rSges(3,2) * t17) + (rSges(3,1) * t42 - rSges(3,2) * t40) * t75) - m(4) * (g(1) * (-t53 * t18 + t73 * t19 - t14) + g(2) * (-t53 * t16 + t73 * t17 - t13) + g(3) * t65 + (rSges(4,3) * t40 + t53 * t42) * t75) - m(5) * (g(1) * (t52 * t18 + t74 * t19 + t56) + g(2) * (t52 * t16 + t74 * t17 + t57) + t51 + (rSges(5,1) * t40 - t52 * t42) * t75) - m(6) * (g(1) * (t85 * t18 + t83 * t19 + t56) + g(2) * (t85 * t16 + t83 * t17 + t57) + t51 + (t40 * pkin(4) + (t35 * t68 + t37 * t40) * rSges(6,1) + (-t35 * t40 + t37 * t68) * rSges(6,2)) * t75 + (-t88 * t41 + t67 * t75) * t87) - m(7) * (t51 + (t48 * t40 + t84 * t42) * t75 - t88 * t84 + (t46 * t17 + t57) * g(2) + (t46 * t19 + t56) * g(1)) -m(4) * (g(1) * (-rSges(4,1) * t7 - rSges(4,2) * t8) + g(2) * (-rSges(4,1) * t5 - rSges(4,2) * t6) + g(3) * (-rSges(4,1) * t20 - rSges(4,2) * t21)) - m(5) * (g(1) * (rSges(5,2) * t7 + t63 * t8 - t4) + g(2) * (rSges(5,2) * t5 + t63 * t6 - t3) + g(3) * (rSges(5,2) * t20 + t63 * t21 - t15)) + (-g(1) * (-t66 * t7 - t4) - g(2) * (-t66 * t5 - t3) - g(3) * (-t66 * t20 - t15) - t86 * (qJ(4) + t45)) * m(7) + (-g(1) * (-t87 * t7 - t4) - g(2) * (-t87 * t5 - t3) - g(3) * (-t87 * t20 - t15) - t86 * (qJ(4) - t92)) * m(6), t58 * (g(1) * t7 + g(2) * t5 + g(3) * t20) t80 * t86, -m(7) * (g(1) * ((-t18 * t32 + t33 * t7) * rSges(7,1) + (-t18 * t33 - t32 * t7) * rSges(7,2)) + g(2) * ((-t16 * t32 + t33 * t5) * rSges(7,1) + (-t16 * t33 - t32 * t5) * rSges(7,2)) + g(3) * ((t20 * t33 + t32 * t70) * rSges(7,1) + (-t20 * t32 + t33 * t70) * rSges(7,2)))];
taug  = t1(:);
