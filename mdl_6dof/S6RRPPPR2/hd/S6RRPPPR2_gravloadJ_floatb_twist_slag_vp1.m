% Calculate Gravitation load on the joints for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:09:58
% EndTime: 2019-03-09 08:10:00
% DurationCPUTime: 0.63s
% Computational Cost: add. (350->116), mult. (390->151), div. (0->0), fcn. (351->10), ass. (0->63)
t23 = qJ(2) + pkin(9);
t18 = sin(t23);
t13 = t18 * qJ(4);
t20 = cos(t23);
t14 = t20 * pkin(3);
t73 = t13 + t14;
t26 = -qJ(3) - pkin(7);
t72 = pkin(4) - t26;
t57 = rSges(7,3) + pkin(8) + qJ(5);
t31 = cos(qJ(1));
t29 = sin(qJ(1));
t68 = g(2) * t29;
t43 = g(1) * t31 + t68;
t54 = rSges(6,3) + qJ(5);
t71 = -m(6) - m(7);
t28 = sin(qJ(2));
t70 = pkin(2) * t28;
t67 = g(3) * t20;
t66 = rSges(3,3) + pkin(7);
t25 = cos(pkin(10));
t65 = rSges(6,2) * t25;
t22 = pkin(10) + qJ(6);
t17 = sin(t22);
t64 = t17 * t31;
t24 = sin(pkin(10));
t63 = t18 * t24;
t19 = cos(t22);
t62 = t19 * t31;
t61 = t29 * t17;
t60 = t29 * t19;
t59 = rSges(5,1) - t26;
t58 = rSges(4,3) - t26;
t56 = pkin(5) * t25 + t72;
t55 = qJ(4) * t20;
t53 = -m(5) + t71;
t52 = pkin(5) * t63;
t51 = -pkin(3) - t57;
t30 = cos(qJ(2));
t21 = t30 * pkin(2);
t16 = t21 + pkin(1);
t12 = t31 * t16;
t50 = t73 * t31 + t12;
t49 = -pkin(3) - t54;
t48 = t21 + t73;
t47 = -t16 - t13;
t46 = g(1) * t51;
t45 = g(2) * t50;
t44 = g(1) * t49;
t42 = rSges(3,1) * t30 - rSges(3,2) * t28;
t40 = rSges(4,1) * t20 - rSges(4,2) * t18;
t39 = rSges(6,1) * t24 + t65;
t38 = rSges(5,2) * t20 - rSges(5,3) * t18;
t37 = pkin(1) + t42;
t36 = rSges(6,1) * t25 - rSges(6,2) * t24 + t72;
t34 = rSges(7,1) * t17 + rSges(7,2) * t19 + pkin(5) * t24;
t7 = t29 * t55;
t9 = t31 * t55;
t32 = g(1) * (-t31 * t70 + t9) + g(2) * (-t29 * t70 + t7) + g(3) * t48;
t6 = -t18 * t61 + t62;
t5 = t18 * t60 + t64;
t4 = t18 * t64 + t60;
t3 = t18 * t62 - t61;
t1 = [-m(2) * (g(1) * (-t29 * rSges(2,1) - rSges(2,2) * t31) + g(2) * (rSges(2,1) * t31 - t29 * rSges(2,2))) - m(3) * ((g(1) * t66 + g(2) * t37) * t31 + (-g(1) * t37 + g(2) * t66) * t29) - m(4) * (g(2) * t12 + (g(1) * t58 + g(2) * t40) * t31 + (g(1) * (-t16 - t40) + g(2) * t58) * t29) - m(5) * (t45 + (g(1) * t59 - g(2) * t38) * t31 + (g(1) * (t38 + t47 - t14) + g(2) * t59) * t29) - m(6) * (t45 + (g(1) * t36 + g(2) * (rSges(6,1) * t63 + t18 * t65 + t54 * t20)) * t31 + (g(2) * t36 + t20 * t44 + (-t16 + (-qJ(4) - t39) * t18) * g(1)) * t29) - m(7) * (g(1) * (t6 * rSges(7,1) - t5 * rSges(7,2)) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t50) + (g(1) * t56 + g(2) * (t57 * t20 + t52)) * t31 + (g(1) * (t47 - t52) + g(2) * t56 + t20 * t46) * t29) -m(3) * (g(3) * t42 + t43 * (-rSges(3,1) * t28 - rSges(3,2) * t30)) - m(4) * (g(3) * (t21 + t40) + t43 * (-rSges(4,1) * t18 - rSges(4,2) * t20 - t70)) - m(5) * (g(1) * t9 + g(2) * t7 + g(3) * (-t38 + t48) + t43 * (rSges(5,3) * t20 - t70 + (rSges(5,2) - pkin(3)) * t18)) - m(6) * ((g(3) * t54 + t43 * t39) * t20 + (g(3) * t39 + t31 * t44 + t49 * t68) * t18 + t32) - m(7) * ((g(3) * t57 + t43 * t34) * t20 + (g(3) * t34 + t31 * t46 + t51 * t68) * t18 + t32) (-m(4) + t53) * (g(1) * t29 - g(2) * t31) t53 * (t43 * t18 - t67) t71 * (g(3) * t18 + t43 * t20) -m(7) * (g(1) * (rSges(7,1) * t3 - rSges(7,2) * t4) + g(2) * (rSges(7,1) * t5 + rSges(7,2) * t6) + (-rSges(7,1) * t19 + rSges(7,2) * t17) * t67)];
taug  = t1(:);
