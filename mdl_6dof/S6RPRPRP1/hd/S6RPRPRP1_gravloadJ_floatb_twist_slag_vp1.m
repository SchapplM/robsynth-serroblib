% Calculate Gravitation load on the joints for
% S6RPRPRP1
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
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:01:21
% EndTime: 2019-03-09 03:01:23
% DurationCPUTime: 0.52s
% Computational Cost: add. (385->107), mult. (336->137), div. (0->0), fcn. (301->10), ass. (0->48)
t39 = rSges(7,3) + qJ(6) + pkin(8);
t15 = qJ(3) + pkin(10);
t9 = sin(t15);
t56 = t39 * t9;
t46 = rSges(6,3) + pkin(8);
t55 = t46 * t9;
t54 = rSges(7,1) + pkin(5);
t16 = qJ(1) + pkin(9);
t10 = sin(t16);
t12 = cos(t16);
t53 = g(1) * t12 + g(2) * t10;
t20 = sin(qJ(3));
t52 = pkin(3) * t20;
t19 = sin(qJ(5));
t51 = pkin(5) * t19;
t21 = sin(qJ(1));
t48 = t21 * pkin(1);
t47 = rSges(4,3) + pkin(7);
t24 = cos(qJ(1));
t14 = t24 * pkin(1);
t23 = cos(qJ(3));
t13 = t23 * pkin(3);
t8 = t13 + pkin(2);
t45 = t12 * t8 + t14;
t44 = t10 * t19;
t22 = cos(qJ(5));
t43 = t10 * t22;
t42 = t12 * t19;
t41 = t12 * t22;
t18 = -qJ(4) - pkin(7);
t40 = rSges(5,3) - t18;
t38 = -m(5) - m(6) - m(7);
t37 = g(1) * t48;
t36 = -t18 + t51;
t11 = cos(t15);
t35 = t11 * rSges(5,1) - t9 * rSges(5,2);
t34 = t23 * rSges(4,1) - t20 * rSges(4,2);
t32 = pkin(2) + t34;
t31 = rSges(6,1) * t22 - rSges(6,2) * t19 + pkin(4);
t7 = t22 * pkin(5) + pkin(4);
t30 = rSges(7,1) * t22 - rSges(7,2) * t19 + t7;
t3 = -t11 * t42 + t43;
t1 = t11 * t44 + t41;
t29 = t11 * pkin(4) + t55;
t28 = t11 * t7 + t56;
t4 = t11 * t41 + t44;
t2 = -t11 * t43 + t42;
t5 = [-m(2) * (g(1) * (-t21 * rSges(2,1) - t24 * rSges(2,2)) + g(2) * (t24 * rSges(2,1) - t21 * rSges(2,2))) - m(3) * (g(1) * (-t10 * rSges(3,1) - t12 * rSges(3,2) - t48) + g(2) * (t12 * rSges(3,1) - t10 * rSges(3,2) + t14)) - m(4) * (-t37 + g(2) * t14 + (g(1) * t47 + g(2) * t32) * t12 + (-g(1) * t32 + g(2) * t47) * t10) - m(5) * (-t37 + g(2) * t45 + (g(1) * t40 + g(2) * t35) * t12 + (g(1) * (-t35 - t8) + g(2) * t40) * t10) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) - t48) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t45) + (-g(1) * t18 + g(2) * t29) * t12 + (g(1) * (-t29 - t8) - g(2) * t18) * t10) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t48) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t45) + (g(1) * t36 + g(2) * t28) * t12 + (g(1) * (-t28 - t8) + g(2) * t36) * t10) (-m(3) - m(4) + t38) * g(3), -m(4) * (g(3) * t34 + t53 * (-rSges(4,1) * t20 - rSges(4,2) * t23)) - m(5) * (g(3) * (t13 + t35) + t53 * (-rSges(5,1) * t9 - rSges(5,2) * t11 - t52)) - m(6) * (g(3) * (t31 * t11 + t13 + t55) + t53 * (t46 * t11 - t31 * t9 - t52)) - m(7) * (g(3) * (t30 * t11 + t13 + t56) + t53 * (t39 * t11 - t30 * t9 - t52)) t38 * (g(1) * t10 - g(2) * t12) -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) + t2 * rSges(6,2))) - m(7) * (g(1) * (-t4 * rSges(7,2) + t54 * t3) + g(2) * (t2 * rSges(7,2) - t54 * t1)) + (-m(6) * (-rSges(6,1) * t19 - rSges(6,2) * t22) - m(7) * (-rSges(7,1) * t19 - rSges(7,2) * t22 - t51)) * g(3) * t9, -m(7) * (-g(3) * t11 + t53 * t9)];
taug  = t5(:);
