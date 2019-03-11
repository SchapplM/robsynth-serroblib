% Calculate Gravitation load on the joints for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:29:41
% EndTime: 2019-03-09 08:29:43
% DurationCPUTime: 0.65s
% Computational Cost: add. (344->119), mult. (420->161), div. (0->0), fcn. (382->8), ass. (0->57)
t19 = qJ(2) + pkin(9);
t16 = sin(t19);
t12 = t16 * qJ(4);
t17 = cos(t19);
t13 = t17 * pkin(3);
t66 = t12 + t13;
t65 = rSges(7,1) + pkin(5);
t21 = -qJ(3) - pkin(7);
t57 = pkin(4) - t21;
t50 = rSges(7,3) + qJ(6) + pkin(8);
t27 = cos(qJ(1));
t24 = sin(qJ(1));
t61 = g(2) * t24;
t38 = g(1) * t27 + t61;
t23 = sin(qJ(2));
t64 = pkin(2) * t23;
t25 = cos(qJ(5));
t63 = pkin(5) * t25;
t60 = g(3) * t17;
t59 = rSges(3,3) + pkin(7);
t58 = rSges(6,3) + pkin(8);
t22 = sin(qJ(5));
t56 = t22 * t27;
t55 = t24 * t22;
t54 = t24 * t25;
t53 = t25 * t27;
t52 = rSges(5,1) - t21;
t51 = rSges(4,3) - t21;
t49 = t63 + t57;
t48 = qJ(4) * t17;
t47 = -m(5) - m(6) - m(7);
t46 = -pkin(3) - t58;
t45 = pkin(5) * t16 * t22;
t44 = -pkin(3) - t50;
t26 = cos(qJ(2));
t18 = t26 * pkin(2);
t15 = t18 + pkin(1);
t11 = t27 * t15;
t43 = t66 * t27 + t11;
t42 = t18 + t66;
t41 = -t15 - t12;
t40 = g(1) * t46;
t39 = g(1) * t44;
t37 = rSges(3,1) * t26 - rSges(3,2) * t23;
t35 = rSges(4,1) * t17 - rSges(4,2) * t16;
t34 = rSges(6,1) * t22 + rSges(6,2) * t25;
t33 = rSges(5,2) * t17 - rSges(5,3) * t16;
t32 = pkin(1) + t37;
t2 = t16 * t53 - t55;
t4 = t16 * t54 + t56;
t31 = rSges(7,2) * t25 + t65 * t22;
t6 = t24 * t48;
t8 = t27 * t48;
t28 = g(1) * (-t27 * t64 + t8) + g(2) * (-t24 * t64 + t6) + g(3) * t42;
t5 = -t16 * t55 + t53;
t3 = t16 * t56 + t54;
t1 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - rSges(2,2) * t27) + g(2) * (rSges(2,1) * t27 - t24 * rSges(2,2))) - m(3) * ((g(1) * t59 + g(2) * t32) * t27 + (-g(1) * t32 + g(2) * t59) * t24) - m(4) * (g(2) * t11 + (g(1) * t51 + g(2) * t35) * t27 + (g(1) * (-t15 - t35) + g(2) * t51) * t24) - m(5) * (g(2) * t43 + (g(1) * t52 - g(2) * t33) * t27 + (g(1) * (t33 + t41 - t13) + g(2) * t52) * t24) - m(6) * (g(1) * (t5 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (t3 * rSges(6,1) + t2 * rSges(6,2) + t43) + (g(2) * t58 * t17 + g(1) * t57) * t27 + (g(1) * t41 + g(2) * t57 + t17 * t40) * t24) - m(7) * (g(1) * (t5 * rSges(7,1) - t4 * rSges(7,2)) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t43) + (g(1) * t49 + g(2) * (t50 * t17 + t45)) * t27 + (g(1) * (t41 - t45) + g(2) * t49 + t17 * t39) * t24) -m(3) * (g(3) * t37 + t38 * (-rSges(3,1) * t23 - rSges(3,2) * t26)) - m(4) * (g(3) * (t18 + t35) + t38 * (-rSges(4,1) * t16 - rSges(4,2) * t17 - t64)) - m(5) * (g(1) * t8 + g(2) * t6 + g(3) * (-t33 + t42) + t38 * (rSges(5,3) * t17 - t64 + (rSges(5,2) - pkin(3)) * t16)) - m(6) * ((g(3) * t58 + t38 * t34) * t17 + (g(3) * t34 + t27 * t40 + t46 * t61) * t16 + t28) - m(7) * ((g(3) * t50 + t38 * t31) * t17 + (g(3) * t31 + t27 * t39 + t44 * t61) * t16 + t28) (-m(4) + t47) * (g(1) * t24 - g(2) * t27) t47 * (t16 * t38 - t60) -m(6) * (g(1) * (rSges(6,1) * t2 - rSges(6,2) * t3) + g(2) * (rSges(6,1) * t4 + rSges(6,2) * t5)) - m(7) * (g(1) * (-t3 * rSges(7,2) + t65 * t2) + g(2) * (t5 * rSges(7,2) + t65 * t4)) + (-m(6) * (-rSges(6,1) * t25 + rSges(6,2) * t22) - m(7) * (-rSges(7,1) * t25 + rSges(7,2) * t22 - t63)) * t60, -m(7) * (g(3) * t16 + t17 * t38)];
taug  = t1(:);
