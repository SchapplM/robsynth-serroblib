% Calculate Gravitation load on the joints for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:36:59
% EndTime: 2019-03-09 03:37:00
% DurationCPUTime: 0.62s
% Computational Cost: add. (435->114), mult. (354->148), div. (0->0), fcn. (320->12), ass. (0->62)
t25 = qJ(3) + pkin(11);
t17 = sin(t25);
t62 = rSges(6,3) + pkin(8);
t74 = t62 * t17;
t52 = rSges(7,3) + pkin(9) + pkin(8);
t73 = t52 * t17;
t26 = qJ(1) + pkin(10);
t18 = sin(t26);
t20 = cos(t26);
t72 = g(1) * t20 + g(2) * t18;
t19 = cos(t25);
t27 = qJ(5) + qJ(6);
t22 = cos(t27);
t56 = t20 * t22;
t21 = sin(t27);
t61 = t18 * t21;
t5 = t19 * t61 + t56;
t57 = t20 * t21;
t60 = t18 * t22;
t6 = -t19 * t60 + t57;
t71 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t7 = -t19 * t57 + t60;
t8 = t19 * t56 + t61;
t70 = t7 * rSges(7,1) - t8 * rSges(7,2);
t30 = sin(qJ(3));
t69 = pkin(3) * t30;
t29 = sin(qJ(5));
t68 = pkin(5) * t29;
t65 = g(3) * t17;
t31 = sin(qJ(1));
t64 = t31 * pkin(1);
t63 = rSges(4,3) + pkin(7);
t59 = t18 * t29;
t32 = cos(qJ(5));
t58 = t18 * t32;
t55 = t20 * t29;
t54 = t20 * t32;
t28 = -qJ(4) - pkin(7);
t53 = rSges(5,3) - t28;
t33 = cos(qJ(3));
t23 = t33 * pkin(3);
t16 = t23 + pkin(2);
t34 = cos(qJ(1));
t24 = t34 * pkin(1);
t51 = t20 * t16 + t24;
t50 = -m(5) - m(6) - m(7);
t49 = g(1) * t64;
t48 = -t28 + t68;
t47 = t33 * rSges(4,1) - t30 * rSges(4,2);
t45 = t19 * rSges(5,1) - t17 * rSges(5,2);
t44 = -rSges(7,1) * t21 - rSges(7,2) * t22;
t43 = pkin(2) + t47;
t42 = rSges(6,1) * t32 - rSges(6,2) * t29 + pkin(4);
t11 = -t19 * t55 + t58;
t9 = t19 * t59 + t54;
t15 = t32 * pkin(5) + pkin(4);
t41 = rSges(7,1) * t22 - rSges(7,2) * t21 + t15;
t40 = t19 * pkin(4) + t74;
t38 = t19 * t15 + t73;
t12 = t19 * t54 + t59;
t10 = -t19 * t58 + t55;
t1 = [-m(2) * (g(1) * (-t31 * rSges(2,1) - t34 * rSges(2,2)) + g(2) * (t34 * rSges(2,1) - t31 * rSges(2,2))) - m(3) * (g(1) * (-t18 * rSges(3,1) - t20 * rSges(3,2) - t64) + g(2) * (t20 * rSges(3,1) - t18 * rSges(3,2) + t24)) - m(4) * (-t49 + g(2) * t24 + (g(1) * t63 + g(2) * t43) * t20 + (-g(1) * t43 + g(2) * t63) * t18) - m(5) * (-t49 + g(2) * t51 + (g(1) * t53 + g(2) * t45) * t20 + (g(1) * (-t16 - t45) + g(2) * t53) * t18) - m(6) * (g(1) * (t10 * rSges(6,1) + t9 * rSges(6,2) - t64) + g(2) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t51) + (-g(1) * t28 + g(2) * t40) * t20 + (g(1) * (-t16 - t40) - g(2) * t28) * t18) - m(7) * (g(1) * (t6 * rSges(7,1) + t5 * rSges(7,2) - t64) + g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) + t51) + (g(1) * t48 + g(2) * t38) * t20 + (g(1) * (-t16 - t38) + g(2) * t48) * t18) (-m(3) - m(4) + t50) * g(3), -m(4) * (g(3) * t47 + t72 * (-rSges(4,1) * t30 - rSges(4,2) * t33)) - m(5) * (g(3) * (t23 + t45) + t72 * (-rSges(5,1) * t17 - rSges(5,2) * t19 - t69)) - m(6) * (g(3) * (t42 * t19 + t23 + t74) + t72 * (-t42 * t17 + t62 * t19 - t69)) - m(7) * (g(3) * (t41 * t19 + t23 + t73) + t72 * (-t41 * t17 + t52 * t19 - t69)) t50 * (g(1) * t18 - g(2) * t20) -m(6) * (g(1) * (t11 * rSges(6,1) - t12 * rSges(6,2)) + g(2) * (-t9 * rSges(6,1) + t10 * rSges(6,2))) - m(7) * (g(1) * (t11 * pkin(5) + t70) + g(2) * (-t9 * pkin(5) + t71)) + (-m(6) * (-rSges(6,1) * t29 - rSges(6,2) * t32) - m(7) * (t44 - t68)) * t65, -m(7) * (g(1) * t70 + g(2) * t71 + t44 * t65)];
taug  = t1(:);
