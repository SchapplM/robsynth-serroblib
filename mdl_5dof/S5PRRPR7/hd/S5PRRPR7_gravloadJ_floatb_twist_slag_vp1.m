% Calculate Gravitation load on the joints for
% S5PRRPR7
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:35:01
% EndTime: 2019-12-05 16:35:05
% DurationCPUTime: 0.74s
% Computational Cost: add. (347->133), mult. (889->215), div. (0->0), fcn. (1073->12), ass. (0->67)
t46 = sin(qJ(3));
t49 = cos(qJ(3));
t88 = -pkin(3) * t49 - qJ(4) * t46;
t47 = sin(qJ(2));
t50 = cos(qJ(2));
t63 = cos(pkin(9));
t64 = cos(pkin(5));
t52 = t64 * t63;
t62 = sin(pkin(9));
t27 = t47 * t52 + t62 * t50;
t43 = sin(pkin(5));
t59 = t43 * t63;
t11 = t27 * t46 + t49 * t59;
t51 = t64 * t62;
t29 = -t47 * t51 + t63 * t50;
t58 = t43 * t62;
t13 = t29 * t46 - t49 * t58;
t75 = t43 * t47;
t30 = t46 * t75 - t64 * t49;
t87 = g(1) * t13 + g(2) * t11 + g(3) * t30;
t86 = -m(5) - m(6);
t81 = g(3) * t43;
t44 = cos(pkin(10));
t80 = t44 * pkin(4);
t79 = rSges(4,3) + pkin(7);
t78 = pkin(8) + rSges(6,3);
t77 = rSges(5,3) * t46;
t42 = sin(pkin(10));
t76 = t42 * t49;
t74 = t43 * t50;
t45 = sin(qJ(5));
t73 = t44 * t45;
t48 = cos(qJ(5));
t72 = t44 * t48;
t71 = t44 * t49;
t70 = t45 * t46;
t69 = t46 * t48;
t68 = t49 * t50;
t67 = pkin(2) * t74 + pkin(7) * t75;
t65 = rSges(5,3) + qJ(4);
t61 = t46 * t74;
t60 = t42 * t74;
t57 = t43 * pkin(3) * t68 + qJ(4) * t61 + t67;
t56 = rSges(4,1) * t49 - rSges(4,2) * t46;
t55 = -rSges(5,1) * t44 + rSges(5,2) * t42;
t26 = t62 * t47 - t50 * t52;
t23 = t26 * pkin(2);
t54 = t27 * pkin(7) + t88 * t26 - t23;
t28 = t63 * t47 + t50 * t51;
t24 = t28 * pkin(2);
t53 = t29 * pkin(7) + t88 * t28 - t24;
t31 = t64 * t46 + t49 * t75;
t25 = t30 * pkin(3);
t16 = (t42 * t47 + t44 * t68) * t43;
t15 = -t44 * t75 + t49 * t60;
t14 = t29 * t49 + t46 * t58;
t12 = t27 * t49 - t46 * t59;
t10 = t31 * t44 - t60;
t9 = t13 * pkin(3);
t8 = t11 * pkin(3);
t7 = -t28 * t71 + t29 * t42;
t6 = -t28 * t76 - t29 * t44;
t5 = -t26 * t71 + t27 * t42;
t4 = -t26 * t76 - t27 * t44;
t3 = t14 * t44 + t28 * t42;
t2 = t12 * t44 + t26 * t42;
t1 = [(-m(2) - m(3) - m(4) + t86) * g(3), -m(3) * (g(1) * (-t28 * rSges(3,1) - t29 * rSges(3,2)) + g(2) * (-t26 * rSges(3,1) - t27 * rSges(3,2)) + (rSges(3,1) * t50 - rSges(3,2) * t47) * t81) - m(4) * (g(1) * (-t56 * t28 + t79 * t29 - t24) + g(2) * (-t56 * t26 + t79 * t27 - t23) + g(3) * t67 + (rSges(4,3) * t47 + t56 * t50) * t81) - m(5) * (g(1) * (t7 * rSges(5,1) - t6 * rSges(5,2) - t28 * t77 + t53) + g(2) * (t5 * rSges(5,1) - t4 * rSges(5,2) - t26 * t77 + t54) + g(3) * (t16 * rSges(5,1) - t15 * rSges(5,2) + rSges(5,3) * t61 + t57)) - m(6) * (g(1) * (t7 * pkin(4) + (-t28 * t70 + t7 * t48) * rSges(6,1) + (-t28 * t69 - t7 * t45) * rSges(6,2) + t78 * t6 + t53) + g(2) * (t5 * pkin(4) + (-t26 * t70 + t5 * t48) * rSges(6,1) + (-t26 * t69 - t5 * t45) * rSges(6,2) + t78 * t4 + t54) + g(3) * (t16 * pkin(4) + (t16 * t48 + t45 * t61) * rSges(6,1) + (-t16 * t45 + t48 * t61) * rSges(6,2) + t78 * t15 + t57)), -m(4) * (g(1) * (-t13 * rSges(4,1) - t14 * rSges(4,2)) + g(2) * (-t11 * rSges(4,1) - t12 * rSges(4,2)) + g(3) * (-t30 * rSges(4,1) - t31 * rSges(4,2))) - m(5) * (g(1) * (t55 * t13 + t65 * t14 - t9) + g(2) * (t55 * t11 + t65 * t12 - t8) + g(3) * (t55 * t30 + t65 * t31 - t25)) + (-g(1) * (-t13 * t80 - t9 + t14 * qJ(4) + (-t13 * t72 + t14 * t45) * rSges(6,1) + (t13 * t73 + t14 * t48) * rSges(6,2)) - g(2) * (-t11 * t80 - t8 + t12 * qJ(4) + (-t11 * t72 + t12 * t45) * rSges(6,1) + (t11 * t73 + t12 * t48) * rSges(6,2)) - g(3) * (-t30 * t80 - t25 + t31 * qJ(4) + (-t30 * t72 + t31 * t45) * rSges(6,1) + (t30 * t73 + t31 * t48) * rSges(6,2)) + t87 * t42 * t78) * m(6), t86 * t87, -m(6) * (g(1) * ((t13 * t48 - t3 * t45) * rSges(6,1) + (-t13 * t45 - t3 * t48) * rSges(6,2)) + g(2) * ((t11 * t48 - t2 * t45) * rSges(6,1) + (-t11 * t45 - t2 * t48) * rSges(6,2)) + g(3) * ((-t10 * t45 + t30 * t48) * rSges(6,1) + (-t10 * t48 - t30 * t45) * rSges(6,2)))];
taug = t1(:);
