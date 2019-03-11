% Calculate Gravitation load on the joints for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:21:12
% EndTime: 2019-03-09 03:21:14
% DurationCPUTime: 0.59s
% Computational Cost: add. (268->114), mult. (354->146), div. (0->0), fcn. (319->8), ass. (0->50)
t17 = qJ(3) + pkin(9);
t13 = cos(t17);
t44 = rSges(6,3) + pkin(8);
t60 = t44 * t13;
t59 = rSges(7,1) + pkin(5);
t39 = rSges(7,3) + qJ(6) + pkin(8);
t25 = cos(qJ(1));
t49 = g(2) * t25;
t22 = sin(qJ(1));
t52 = g(1) * t22;
t34 = -t49 + t52;
t58 = g(1) * t25 + g(2) * t22;
t12 = sin(t17);
t23 = cos(qJ(5));
t11 = t23 * pkin(5) + pkin(4);
t20 = sin(qJ(5));
t27 = rSges(7,1) * t23 - rSges(7,2) * t20 + t11;
t57 = t39 * t12 + t27 * t13;
t28 = rSges(6,1) * t23 - rSges(6,2) * t20 + pkin(4);
t56 = t44 * t12 + t28 * t13;
t24 = cos(qJ(3));
t54 = pkin(3) * t24;
t8 = t22 * t54;
t55 = g(1) * t8;
t53 = pkin(5) * t20;
t48 = g(3) * t12;
t47 = t12 * pkin(4);
t21 = sin(qJ(3));
t46 = t21 * pkin(3);
t45 = rSges(4,3) + pkin(7);
t43 = t22 * t20;
t42 = t22 * t23;
t41 = t25 * t20;
t40 = t25 * t23;
t38 = t25 * pkin(1) + t22 * qJ(2);
t37 = -m(5) - m(6) - m(7);
t15 = t25 * qJ(2);
t19 = -qJ(4) - pkin(7);
t36 = t22 * t19 + t25 * t46 + t15;
t35 = t22 * t46 + t38;
t33 = t39 * t13;
t31 = t21 * rSges(4,1) + t24 * rSges(4,2);
t30 = rSges(5,1) * t13 - rSges(5,2) * t12;
t29 = t12 * rSges(5,1) + t13 * rSges(5,2);
t3 = t12 * t41 + t42;
t1 = -t12 * t43 + t40;
t26 = t12 * t11 - t33;
t4 = t12 * t40 - t43;
t2 = t12 * t42 + t41;
t5 = [-m(2) * (g(1) * (-t22 * rSges(2,1) - t25 * rSges(2,2)) + g(2) * (t25 * rSges(2,1) - t22 * rSges(2,2))) - m(3) * (g(1) * (t25 * rSges(3,3) + t15 + (rSges(3,2) - pkin(1)) * t22) + g(2) * (-t25 * rSges(3,2) + t22 * rSges(3,3) + t38)) - m(4) * (g(1) * t15 + g(2) * t38 + (g(1) * t31 + g(2) * t45) * t25 + (g(1) * (-pkin(1) - t45) + g(2) * t31) * t22) - m(5) * (g(1) * t36 + g(2) * t35 + (g(1) * t29 + g(2) * (rSges(5,3) - t19)) * t25 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t29) * t22) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) - t22 * pkin(1) + t25 * t47 + t36) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) - t25 * t19 + t22 * t47 + t35) - t58 * t60) - m(7) * (g(1) * (t4 * rSges(7,1) - t3 * rSges(7,2) + t36) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t35) + (g(1) * t26 + g(2) * (-t19 + t53)) * t25 + (g(1) * (-pkin(1) - t53) + g(2) * t26) * t22) (-m(3) - m(4) + t37) * t34, -m(4) * (-g(3) * t31 + t34 * (rSges(4,1) * t24 - rSges(4,2) * t21)) - m(5) * (g(1) * (t30 * t22 + t8) + g(3) * (-t29 - t46) + (-t30 - t54) * t49) - m(6) * (t55 + g(3) * (-t46 + t60) - t28 * t48 + t56 * t52 + (-t54 - t56) * t49) - m(7) * (t55 + g(3) * (t33 - t46) - t27 * t48 + t57 * t52 + (-t54 - t57) * t49) t37 * t58, -m(6) * (g(1) * (t1 * rSges(6,1) - t2 * rSges(6,2)) + g(2) * (t3 * rSges(6,1) + t4 * rSges(6,2))) - m(7) * (g(1) * (-t2 * rSges(7,2) + t59 * t1) + g(2) * (t4 * rSges(7,2) + t59 * t3)) + (-m(6) * (-rSges(6,1) * t20 - rSges(6,2) * t23) - m(7) * (-rSges(7,1) * t20 - rSges(7,2) * t23 - t53)) * g(3) * t13, -m(7) * (-t34 * t13 + t48)];
taug  = t5(:);
