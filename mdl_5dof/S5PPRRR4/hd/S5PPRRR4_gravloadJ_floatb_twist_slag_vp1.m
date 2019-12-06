% Calculate Gravitation load on the joints for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:36
% EndTime: 2019-12-05 15:18:39
% DurationCPUTime: 0.47s
% Computational Cost: add. (436->92), mult. (1199->159), div. (0->0), fcn. (1520->14), ass. (0->56)
t52 = sin(pkin(11));
t53 = sin(pkin(10));
t41 = t53 * t52;
t56 = cos(pkin(11));
t57 = cos(pkin(10));
t48 = t57 * t56;
t59 = cos(pkin(5));
t35 = -t59 * t48 + t41;
t54 = sin(pkin(6));
t55 = sin(pkin(5));
t45 = t55 * t54;
t58 = cos(pkin(6));
t68 = t35 * t58 + t57 * t45;
t42 = t53 * t56;
t46 = t57 * t52;
t36 = t59 * t42 + t46;
t44 = t55 * t53;
t67 = t36 * t58 - t54 * t44;
t66 = t56 * t58 * t55 + t59 * t54;
t32 = cos(qJ(4));
t65 = t32 * pkin(4);
t64 = rSges(5,3) + pkin(8);
t63 = rSges(6,3) + pkin(9);
t62 = cos(qJ(3));
t28 = sin(qJ(5));
t61 = t28 * t32;
t31 = cos(qJ(5));
t60 = t31 * t32;
t51 = -m(3) - m(4) - m(5) - m(6);
t29 = sin(qJ(4));
t50 = -rSges(5,1) * t32 + rSges(5,2) * t29;
t47 = t57 * t55;
t43 = t55 * t52;
t40 = rSges(6,1) * t31 - rSges(6,2) * t28 + pkin(4);
t30 = sin(qJ(3));
t23 = -t59 * t41 + t48;
t22 = t59 * t46 + t42;
t21 = -t56 * t45 + t59 * t58;
t17 = t36 * t54 + t58 * t44;
t16 = t35 * t54 - t58 * t47;
t15 = t30 * t66 + t62 * t43;
t14 = t30 * t43 - t62 * t66;
t13 = t14 * pkin(3);
t12 = t15 * t32 + t21 * t29;
t11 = -t15 * t29 + t21 * t32;
t10 = t23 * t62 - t30 * t67;
t9 = t23 * t30 + t62 * t67;
t8 = t22 * t62 - t30 * t68;
t7 = t22 * t30 + t62 * t68;
t6 = t9 * pkin(3);
t5 = t7 * pkin(3);
t4 = t10 * t32 + t17 * t29;
t3 = -t10 * t29 + t17 * t32;
t2 = t16 * t29 + t8 * t32;
t1 = t16 * t32 - t8 * t29;
t18 = [(-m(2) + t51) * g(3), t51 * (g(1) * t44 - g(2) * t47 + g(3) * t59), -m(4) * (g(1) * (-t9 * rSges(4,1) - t10 * rSges(4,2)) + g(2) * (-t7 * rSges(4,1) - t8 * rSges(4,2)) + g(3) * (-t14 * rSges(4,1) - t15 * rSges(4,2))) - m(5) * (g(1) * (t64 * t10 + t50 * t9 - t6) + g(2) * (t50 * t7 + t64 * t8 - t5) + g(3) * (t50 * t14 + t64 * t15 - t13)) + (-g(1) * (-t9 * t65 - t6 + t10 * pkin(8) + (t10 * t28 - t9 * t60) * rSges(6,1) + (t10 * t31 + t9 * t61) * rSges(6,2)) - g(2) * (-t7 * t65 - t5 + t8 * pkin(8) + (t8 * t28 - t7 * t60) * rSges(6,1) + (t8 * t31 + t7 * t61) * rSges(6,2)) - g(3) * (-t14 * t65 - t13 + t15 * pkin(8) + (-t14 * t60 + t15 * t28) * rSges(6,1) + (t14 * t61 + t15 * t31) * rSges(6,2)) - (-g(1) * t9 - g(2) * t7 - g(3) * t14) * t29 * t63) * m(6), -m(5) * (g(1) * (t3 * rSges(5,1) - t4 * rSges(5,2)) + g(2) * (t1 * rSges(5,1) - t2 * rSges(5,2)) + g(3) * (t11 * rSges(5,1) - t12 * rSges(5,2))) - m(6) * (g(1) * (t40 * t3 + t63 * t4) + (t40 * t11 + t63 * t12) * g(3) + (t40 * t1 + t63 * t2) * g(2)), -m(6) * (g(1) * ((-t4 * t28 + t9 * t31) * rSges(6,1) + (-t9 * t28 - t4 * t31) * rSges(6,2)) + g(2) * ((-t2 * t28 + t7 * t31) * rSges(6,1) + (-t2 * t31 - t7 * t28) * rSges(6,2)) + g(3) * ((-t12 * t28 + t14 * t31) * rSges(6,1) + (-t12 * t31 - t14 * t28) * rSges(6,2)))];
taug = t18(:);
