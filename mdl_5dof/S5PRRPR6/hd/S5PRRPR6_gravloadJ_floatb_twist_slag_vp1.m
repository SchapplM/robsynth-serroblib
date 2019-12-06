% Calculate Gravitation load on the joints for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:07
% EndTime: 2019-12-05 16:30:09
% DurationCPUTime: 0.50s
% Computational Cost: add. (311->85), mult. (682->132), div. (0->0), fcn. (797->12), ass. (0->50)
t30 = sin(qJ(2));
t32 = cos(qJ(2));
t48 = sin(pkin(9));
t50 = cos(pkin(5));
t41 = t50 * t48;
t49 = cos(pkin(9));
t11 = -t30 * t41 + t49 * t32;
t42 = t50 * t49;
t9 = t30 * t42 + t48 * t32;
t71 = g(1) * t11 + g(2) * t9;
t10 = t49 * t30 + t32 * t41;
t8 = t48 * t30 - t32 * t42;
t70 = g(1) * t10 + g(2) * t8;
t29 = sin(qJ(3));
t31 = cos(qJ(3));
t26 = sin(pkin(5));
t54 = t26 * t30;
t12 = t29 * t54 - t50 * t31;
t46 = t26 * t49;
t2 = t9 * t29 + t31 * t46;
t45 = t26 * t48;
t4 = t11 * t29 - t31 * t45;
t69 = g(1) * t4 + g(2) * t2 + g(3) * t12;
t13 = t50 * t29 + t31 * t54;
t3 = -t29 * t46 + t9 * t31;
t5 = t11 * t31 + t29 * t45;
t68 = g(1) * t5 + g(2) * t3 + g(3) * t13;
t24 = pkin(10) + qJ(5);
t22 = sin(t24);
t23 = cos(t24);
t27 = cos(pkin(10));
t38 = t23 * rSges(6,1) - t22 * rSges(6,2) + t27 * pkin(4) + pkin(3);
t52 = rSges(6,3) + pkin(8) + qJ(4);
t67 = t52 * t29 + t38 * t31;
t25 = sin(pkin(10));
t40 = rSges(5,1) * t27 - rSges(5,2) * t25 + pkin(3);
t51 = rSges(5,3) + qJ(4);
t66 = t51 * t29 + t40 * t31;
t61 = -m(5) - m(6);
t56 = g(3) * t26;
t55 = rSges(4,3) + pkin(7);
t53 = t26 * t32;
t47 = g(3) * (pkin(2) * t53 + pkin(7) * t54);
t44 = rSges(4,1) * t31 - rSges(4,2) * t29;
t43 = t25 * rSges(5,1) + t27 * rSges(5,2);
t37 = t22 * rSges(6,1) + t23 * rSges(6,2) + t25 * pkin(4);
t6 = t8 * pkin(2);
t7 = t10 * pkin(2);
t35 = -g(1) * t7 - g(2) * t6 + t47;
t1 = [(-m(2) - m(3) - m(4) + t61) * g(3), -m(3) * (g(1) * (-t10 * rSges(3,1) - t11 * rSges(3,2)) + g(2) * (-t8 * rSges(3,1) - t9 * rSges(3,2)) + (rSges(3,1) * t32 - rSges(3,2) * t30) * t56) - m(4) * (g(1) * (-t44 * t10 + t55 * t11 - t7) + g(2) * (-t44 * t8 + t55 * t9 - t6) + t47 + (rSges(4,3) * t30 + t44 * t32) * t56) - m(5) * ((t43 * t30 + t66 * t32) * t56 + t35 + t71 * (pkin(7) + t43) - t70 * t66) - m(6) * ((t37 * t30 + t67 * t32) * t56 + t35 + t71 * (pkin(7) + t37) - t70 * t67), -m(4) * (g(1) * (-t4 * rSges(4,1) - t5 * rSges(4,2)) + g(2) * (-t2 * rSges(4,1) - t3 * rSges(4,2)) + g(3) * (-t12 * rSges(4,1) - t13 * rSges(4,2))) - m(5) * (-t69 * t40 + t68 * t51) - m(6) * (-t69 * t38 + t68 * t52), t61 * t69, -m(6) * (g(1) * ((t10 * t23 - t5 * t22) * rSges(6,1) + (-t10 * t22 - t5 * t23) * rSges(6,2)) + g(2) * ((-t3 * t22 + t8 * t23) * rSges(6,1) + (-t8 * t22 - t3 * t23) * rSges(6,2)) + g(3) * ((-t13 * t22 - t23 * t53) * rSges(6,1) + (-t13 * t23 + t22 * t53) * rSges(6,2)))];
taug = t1(:);
