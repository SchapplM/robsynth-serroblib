% Calculate Gravitation load on the joints for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:03
% EndTime: 2019-12-31 20:20:04
% DurationCPUTime: 0.48s
% Computational Cost: add. (287->100), mult. (332->136), div. (0->0), fcn. (306->10), ass. (0->55)
t22 = qJ(2) + pkin(9);
t17 = sin(t22);
t55 = rSges(5,3) + pkin(7);
t66 = t55 * t17;
t45 = rSges(6,3) + pkin(8) + pkin(7);
t65 = t45 * t17;
t27 = sin(qJ(1));
t30 = cos(qJ(1));
t64 = g(1) * t30 + g(2) * t27;
t18 = cos(t22);
t23 = qJ(4) + qJ(5);
t20 = cos(t23);
t49 = t30 * t20;
t19 = sin(t23);
t54 = t27 * t19;
t5 = t18 * t54 + t49;
t50 = t30 * t19;
t53 = t27 * t20;
t6 = -t18 * t53 + t50;
t63 = -t5 * rSges(6,1) + t6 * rSges(6,2);
t7 = -t18 * t50 + t53;
t8 = t18 * t49 + t54;
t62 = t7 * rSges(6,1) - t8 * rSges(6,2);
t26 = sin(qJ(2));
t61 = pkin(2) * t26;
t25 = sin(qJ(4));
t60 = pkin(4) * t25;
t57 = g(3) * t17;
t56 = rSges(3,3) + pkin(6);
t52 = t27 * t25;
t28 = cos(qJ(4));
t51 = t27 * t28;
t48 = t30 * t25;
t47 = t30 * t28;
t24 = -qJ(3) - pkin(6);
t46 = rSges(4,3) - t24;
t44 = -t24 + t60;
t29 = cos(qJ(2));
t43 = t29 * rSges(3,1) - t26 * rSges(3,2);
t41 = t18 * rSges(4,1) - t17 * rSges(4,2);
t40 = -rSges(6,1) * t19 - rSges(6,2) * t20;
t39 = pkin(1) + t43;
t38 = rSges(5,1) * t28 - rSges(5,2) * t25 + pkin(3);
t11 = -t18 * t48 + t51;
t9 = t18 * t52 + t47;
t15 = t28 * pkin(4) + pkin(3);
t37 = rSges(6,1) * t20 - rSges(6,2) * t19 + t15;
t36 = t18 * pkin(3) + t66;
t34 = t18 * t15 + t65;
t21 = t29 * pkin(2);
t16 = t21 + pkin(1);
t14 = t30 * t16;
t12 = t18 * t47 + t52;
t10 = -t18 * t51 + t48;
t1 = [-m(2) * (g(1) * (-t27 * rSges(2,1) - t30 * rSges(2,2)) + g(2) * (t30 * rSges(2,1) - t27 * rSges(2,2))) - m(3) * ((g(1) * t56 + g(2) * t39) * t30 + (-g(1) * t39 + g(2) * t56) * t27) - m(4) * (g(2) * t14 + (g(1) * t46 + g(2) * t41) * t30 + (g(1) * (-t16 - t41) + g(2) * t46) * t27) - m(5) * (g(1) * (t10 * rSges(5,1) + t9 * rSges(5,2)) + g(2) * (t12 * rSges(5,1) + t11 * rSges(5,2) + t14) + (-g(1) * t24 + g(2) * t36) * t30 + (g(1) * (-t16 - t36) - g(2) * t24) * t27) - m(6) * (g(1) * (t6 * rSges(6,1) + t5 * rSges(6,2)) + g(2) * (t8 * rSges(6,1) + t7 * rSges(6,2) + t14) + (g(1) * t44 + g(2) * t34) * t30 + (g(1) * (-t16 - t34) + g(2) * t44) * t27), -m(3) * (g(3) * t43 + t64 * (-rSges(3,1) * t26 - rSges(3,2) * t29)) - m(4) * (g(3) * (t21 + t41) + t64 * (-rSges(4,1) * t17 - rSges(4,2) * t18 - t61)) - m(5) * (g(3) * (t38 * t18 + t21 + t66) + t64 * (-t38 * t17 + t55 * t18 - t61)) - m(6) * (g(3) * (t37 * t18 + t21 + t65) + t64 * (-t37 * t17 + t45 * t18 - t61)), (-m(4) - m(5) - m(6)) * (g(1) * t27 - g(2) * t30), -m(5) * (g(1) * (t11 * rSges(5,1) - t12 * rSges(5,2)) + g(2) * (-t9 * rSges(5,1) + t10 * rSges(5,2))) - m(6) * (g(1) * (t11 * pkin(4) + t62) + g(2) * (-pkin(4) * t9 + t63)) + (-m(5) * (-rSges(5,1) * t25 - rSges(5,2) * t28) - m(6) * (t40 - t60)) * t57, -m(6) * (g(1) * t62 + g(2) * t63 + t40 * t57)];
taug = t1(:);
