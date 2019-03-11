% Calculate Gravitation load on the joints for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:04:30
% EndTime: 2019-03-09 03:04:32
% DurationCPUTime: 0.64s
% Computational Cost: add. (423->112), mult. (372->152), div. (0->0), fcn. (347->10), ass. (0->51)
t57 = rSges(7,1) + pkin(5);
t22 = qJ(3) + pkin(10);
t16 = sin(t22);
t18 = cos(t22);
t66 = t18 * pkin(4) + t16 * pkin(8);
t23 = qJ(1) + pkin(9);
t17 = sin(t23);
t19 = cos(t23);
t65 = g(1) * t19 + g(2) * t17;
t46 = rSges(7,3) + qJ(6);
t25 = sin(qJ(5));
t28 = cos(qJ(5));
t64 = t46 * t25 + t57 * t28;
t26 = sin(qJ(3));
t63 = pkin(3) * t26;
t24 = -qJ(4) - pkin(7);
t60 = g(2) * t24;
t59 = g(3) * t16;
t27 = sin(qJ(1));
t58 = t27 * pkin(1);
t56 = rSges(4,3) + pkin(7);
t55 = t16 * t19;
t54 = t17 * t18;
t53 = t18 * t19;
t52 = t18 * t25;
t51 = t18 * t28;
t50 = t19 * t25;
t49 = t19 * t28;
t48 = rSges(5,3) - t24;
t29 = cos(qJ(3));
t20 = t29 * pkin(3);
t15 = t20 + pkin(2);
t30 = cos(qJ(1));
t21 = t30 * pkin(1);
t47 = t19 * t15 + t21;
t45 = -m(5) - m(6) - m(7);
t44 = g(1) * t58;
t43 = t20 + t66;
t42 = pkin(4) * t53 + pkin(8) * t55 + t47;
t41 = pkin(8) * t54 - t17 * t63;
t40 = pkin(8) * t53 - t19 * t63;
t39 = -t19 * t24 - t58;
t38 = rSges(4,1) * t29 - rSges(4,2) * t26;
t36 = rSges(5,1) * t18 - rSges(5,2) * t16;
t35 = -t15 - t66;
t34 = pkin(2) + t38;
t4 = t17 * t25 + t18 * t49;
t3 = -t17 * t28 + t18 * t50;
t2 = t17 * t51 - t50;
t1 = t17 * t52 + t49;
t5 = [-m(2) * (g(1) * (-t27 * rSges(2,1) - rSges(2,2) * t30) + g(2) * (rSges(2,1) * t30 - t27 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t17 - rSges(3,2) * t19 - t58) + g(2) * (rSges(3,1) * t19 - rSges(3,2) * t17 + t21)) - m(4) * (-t44 + g(2) * t21 + (g(1) * t56 + g(2) * t34) * t19 + (-g(1) * t34 + g(2) * t56) * t17) - m(5) * (-t44 + g(2) * t47 + (g(1) * t48 + g(2) * t36) * t19 + (g(1) * (-t15 - t36) + g(2) * t48) * t17) - m(6) * (g(1) * (-rSges(6,1) * t2 + rSges(6,2) * t1 + t39) + g(2) * (rSges(6,1) * t4 - rSges(6,2) * t3 + rSges(6,3) * t55 + t42) + (g(1) * (-t16 * rSges(6,3) + t35) - t60) * t17) - m(7) * (g(1) * (-t46 * t1 - t57 * t2 + t39) + g(2) * (rSges(7,2) * t55 + t46 * t3 + t57 * t4 + t42) + (g(1) * (-t16 * rSges(7,2) + t35) - t60) * t17) (-m(3) - m(4) + t45) * g(3), -m(4) * (g(3) * t38 + t65 * (-rSges(4,1) * t26 - rSges(4,2) * t29)) - m(5) * (g(3) * (t20 + t36) + t65 * (-rSges(5,1) * t16 - rSges(5,2) * t18 - t63)) - m(6) * (g(1) * (rSges(6,3) * t53 + t40) + g(2) * (rSges(6,3) * t54 + t41) + g(3) * (rSges(6,1) * t51 - rSges(6,2) * t52 + t43) + (g(3) * rSges(6,3) + t65 * (-rSges(6,1) * t28 + rSges(6,2) * t25 - pkin(4))) * t16) - m(7) * (g(1) * t40 + g(2) * t41 + g(3) * t43 + (t65 * rSges(7,2) + g(3) * t64) * t18 + (g(3) * rSges(7,2) + t65 * (-pkin(4) - t64)) * t16) t45 * (g(1) * t17 - g(2) * t19) -m(6) * (g(1) * (-rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2)) - m(7) * (g(1) * (-t57 * t3 + t46 * t4) + g(2) * (-t57 * t1 + t46 * t2)) + (-m(6) * (-rSges(6,1) * t25 - rSges(6,2) * t28) - m(7) * (-t57 * t25 + t46 * t28)) * t59, -m(7) * (g(1) * t3 + g(2) * t1 + t25 * t59)];
taug  = t5(:);
