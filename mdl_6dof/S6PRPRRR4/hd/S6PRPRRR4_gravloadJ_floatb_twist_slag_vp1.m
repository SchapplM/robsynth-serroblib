% Calculate Gravitation load on the joints for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:36:07
% EndTime: 2019-03-08 20:36:09
% DurationCPUTime: 0.69s
% Computational Cost: add. (585->116), mult. (924->175), div. (0->0), fcn. (1079->14), ass. (0->60)
t36 = sin(pkin(6));
t42 = cos(qJ(2));
t67 = t36 * t42;
t86 = g(3) * t67;
t39 = sin(qJ(5));
t41 = cos(qJ(5));
t85 = t39 * rSges(6,1) + t41 * rSges(6,2);
t40 = sin(qJ(2));
t35 = sin(pkin(11));
t60 = cos(pkin(6));
t57 = t35 * t60;
t59 = cos(pkin(11));
t20 = t40 * t59 + t42 * t57;
t52 = t60 * t59;
t18 = t35 * t40 - t42 * t52;
t76 = g(2) * t18;
t84 = g(1) * t20 + t76;
t21 = -t40 * t57 + t42 * t59;
t32 = pkin(12) + qJ(4);
t28 = sin(t32);
t29 = cos(t32);
t69 = t35 * t36;
t11 = t21 * t29 + t28 * t69;
t68 = t36 * t40;
t15 = t28 * t60 + t29 * t68;
t19 = t35 * t42 + t40 * t52;
t56 = t36 * t59;
t9 = t19 * t29 - t28 * t56;
t83 = g(1) * t11 + g(2) * t9 + g(3) * t15;
t10 = -t21 * t28 + t29 * t69;
t14 = -t28 * t68 + t29 * t60;
t8 = -t19 * t28 - t29 * t56;
t82 = g(1) * t10 + g(2) * t8 + g(3) * t14;
t33 = qJ(5) + qJ(6);
t30 = sin(t33);
t31 = cos(t33);
t48 = rSges(7,1) * t31 - rSges(7,2) * t30 + pkin(5) * t41 + pkin(4);
t64 = rSges(7,3) + pkin(10) + pkin(9);
t81 = t28 * t64 + t29 * t48;
t50 = rSges(6,1) * t41 - rSges(6,2) * t39 + pkin(4);
t70 = rSges(6,3) + pkin(9);
t80 = t28 * t70 + t29 * t50;
t75 = g(2) * t19;
t37 = cos(pkin(12));
t26 = pkin(3) * t37 + pkin(2);
t72 = t26 * t86;
t71 = g(3) * t36;
t38 = -pkin(8) - qJ(3);
t63 = -t18 * t26 - t19 * t38;
t62 = -t20 * t26 - t21 * t38;
t61 = rSges(4,3) + qJ(3);
t58 = -m(4) - m(5) - m(6) - m(7);
t55 = t18 * t41 - t39 * t9;
t54 = rSges(5,1) * t29 - rSges(5,2) * t28;
t53 = -t11 * t39 + t20 * t41;
t51 = rSges(4,1) * t37 - rSges(4,2) * sin(pkin(12)) + pkin(2);
t49 = -t15 * t39 - t41 * t67;
t47 = t30 * rSges(7,1) + t31 * rSges(7,2) + t39 * pkin(5);
t45 = m(7) * (g(1) * ((-t11 * t30 + t20 * t31) * rSges(7,1) + (-t11 * t31 - t20 * t30) * rSges(7,2)) + g(2) * ((t18 * t31 - t30 * t9) * rSges(7,1) + (-t18 * t30 - t31 * t9) * rSges(7,2)) + g(3) * ((-t15 * t30 - t31 * t67) * rSges(7,1) + (-t15 * t31 + t30 * t67) * rSges(7,2)));
t1 = [(-m(2) - m(3) + t58) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t20 - rSges(3,2) * t21) + g(2) * (-rSges(3,1) * t18 - rSges(3,2) * t19) + (rSges(3,1) * t42 - rSges(3,2) * t40) * t71) - m(4) * (g(1) * (-t20 * t51 + t21 * t61) + t61 * t75 - t51 * t76 + (t40 * t61 + t42 * t51) * t71) - m(5) * (g(1) * (rSges(5,3) * t21 - t20 * t54 + t62) + g(2) * (rSges(5,3) * t19 - t18 * t54 + t63) + t72 + (t54 * t42 + (rSges(5,3) - t38) * t40) * t71) - m(6) * (g(1) * (t85 * t21 + t62) + g(2) * (t85 * t19 + t63) + t72 + ((-t38 + t85) * t40 + t80 * t42) * t71 - t84 * t80) - m(7) * (g(2) * t63 + t72 + t47 * t75 + ((-t38 + t47) * t40 + t81 * t42) * t71 - t84 * t81 + (t21 * t47 + t62) * g(1)) t58 * (t84 - t86) -m(5) * (g(1) * (rSges(5,1) * t10 - rSges(5,2) * t11) + g(2) * (rSges(5,1) * t8 - rSges(5,2) * t9) + g(3) * (rSges(5,1) * t14 - rSges(5,2) * t15)) - m(6) * (t82 * t50 + t83 * t70) - m(7) * (t82 * t48 + t83 * t64) -m(6) * (g(1) * (t53 * rSges(6,1) + (-t11 * t41 - t20 * t39) * rSges(6,2)) + g(2) * (t55 * rSges(6,1) + (-t18 * t39 - t41 * t9) * rSges(6,2)) + g(3) * (t49 * rSges(6,1) + (-t15 * t41 + t39 * t67) * rSges(6,2))) - t45 - m(7) * (g(1) * t53 + g(2) * t55 + g(3) * t49) * pkin(5), -t45];
taug  = t1(:);
