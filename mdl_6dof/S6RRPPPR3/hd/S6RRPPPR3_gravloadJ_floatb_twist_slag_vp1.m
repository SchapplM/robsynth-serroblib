% Calculate Gravitation load on the joints for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:13:53
% EndTime: 2019-03-09 08:13:54
% DurationCPUTime: 0.64s
% Computational Cost: add. (244->128), mult. (423->168), div. (0->0), fcn. (384->8), ass. (0->60)
t56 = rSges(7,3) + pkin(8) + qJ(5);
t27 = sin(qJ(1));
t68 = g(2) * t27;
t29 = cos(qJ(1));
t69 = g(1) * t29;
t40 = t68 + t69;
t71 = -m(6) - m(7);
t70 = -pkin(2) - pkin(3);
t28 = cos(qJ(2));
t67 = g(3) * t28;
t19 = t28 * pkin(2);
t65 = rSges(4,1) * t28;
t26 = sin(qJ(2));
t64 = rSges(5,1) * t26;
t63 = rSges(5,2) * t28;
t24 = cos(pkin(9));
t13 = pkin(5) * t24 + pkin(4);
t62 = t13 * t26;
t61 = t26 * t29;
t22 = pkin(9) + qJ(6);
t14 = sin(t22);
t60 = t27 * t14;
t15 = cos(t22);
t59 = t27 * t15;
t58 = t27 * t28;
t57 = t28 * t29;
t16 = t26 * qJ(3);
t55 = t16 + t19;
t54 = t29 * pkin(1) + t27 * pkin(7);
t53 = qJ(3) * t28;
t52 = -rSges(5,3) - qJ(4);
t51 = rSges(6,3) + qJ(5);
t50 = -m(5) + t71;
t49 = rSges(5,2) + t70;
t48 = t28 * pkin(3) + t55;
t47 = -t56 + t70;
t46 = -pkin(1) - t16;
t23 = sin(pkin(9));
t45 = -pkin(5) * t23 - qJ(4);
t44 = -t51 + t70;
t43 = pkin(2) * t57 + t29 * t16 + t54;
t42 = g(1) * t49;
t41 = pkin(3) * t57 + t43;
t39 = g(1) * t47;
t38 = g(1) * t44;
t37 = rSges(3,1) * t28 - rSges(3,2) * t26;
t35 = rSges(6,1) * t24 - rSges(6,2) * t23 + pkin(4);
t34 = rSges(7,1) * t15 - rSges(7,2) * t14 + t13;
t33 = -t23 * rSges(6,1) - t24 * rSges(6,2) - qJ(4);
t32 = g(2) * t35;
t20 = t29 * pkin(7);
t31 = g(1) * t20 + g(2) * t41;
t10 = t29 * t53;
t8 = t27 * t53;
t30 = g(1) * t10 + g(2) * t8 + g(3) * t48;
t4 = t15 * t61 - t60;
t3 = -t14 * t61 - t59;
t2 = -t14 * t29 - t26 * t59;
t1 = -t15 * t29 + t26 * t60;
t5 = [-m(2) * (g(1) * (-t27 * rSges(2,1) - rSges(2,2) * t29) + g(2) * (rSges(2,1) * t29 - t27 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t29 + t20) + g(2) * (rSges(3,1) * t57 - rSges(3,2) * t61 + t54) + (g(1) * (-pkin(1) - t37) + g(2) * rSges(3,3)) * t27) - m(4) * (g(1) * (rSges(4,2) * t29 + t20) + g(2) * (rSges(4,1) * t57 + rSges(4,3) * t61 + t43) + (g(1) * (-rSges(4,3) * t26 - t19 + t46 - t65) + g(2) * rSges(4,2)) * t27) - m(5) * ((g(1) * t52 + g(2) * (-t63 + t64)) * t29 + (g(1) * (t46 - t64) + g(2) * t52 + t28 * t42) * t27 + t31) - m(6) * ((g(2) * t51 * t28 + g(1) * t33 + t26 * t32) * t29 + (g(2) * t33 + t28 * t38 + (-pkin(1) + (-qJ(3) - t35) * t26) * g(1)) * t27 + t31) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t20) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t41) + (g(1) * t45 + g(2) * (t56 * t28 + t62)) * t29 + (g(1) * (t46 - t62) + g(2) * t45 + t28 * t39) * t27) -m(3) * (g(3) * t37 + t40 * (-rSges(3,1) * t26 - rSges(3,2) * t28)) - m(4) * (g(1) * (rSges(4,3) * t57 + t10) + g(2) * (rSges(4,3) * t58 + t8) + g(3) * (t55 + t65) + (g(3) * rSges(4,3) + t40 * (-rSges(4,1) - pkin(2))) * t26) - m(5) * (g(1) * (rSges(5,1) * t57 + t10) + g(2) * (rSges(5,1) * t58 + t8) + g(3) * (t48 - t63) + (g(3) * rSges(5,1) + t29 * t42 + t49 * t68) * t26) - m(6) * ((g(3) * t51 + t27 * t32 + t35 * t69) * t28 + (g(3) * t35 + t29 * t38 + t44 * t68) * t26 + t30) - m(7) * ((g(3) * t56 + t40 * t34) * t28 + (g(3) * t34 + t29 * t39 + t47 * t68) * t26 + t30) (-m(4) + t50) * (t40 * t26 - t67) t50 * (-g(1) * t27 + g(2) * t29) t71 * (g(3) * t26 + t40 * t28) -m(7) * (g(1) * (rSges(7,1) * t3 - rSges(7,2) * t4) + g(2) * (-rSges(7,1) * t1 + rSges(7,2) * t2) + (rSges(7,1) * t14 + rSges(7,2) * t15) * t67)];
taug  = t5(:);
