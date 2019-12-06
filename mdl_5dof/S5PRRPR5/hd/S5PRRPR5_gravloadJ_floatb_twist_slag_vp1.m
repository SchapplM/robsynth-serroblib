% Calculate Gravitation load on the joints for
% S5PRRPR5
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:19
% EndTime: 2019-12-05 16:25:23
% DurationCPUTime: 0.64s
% Computational Cost: add. (317->101), mult. (580->156), div. (0->0), fcn. (662->12), ass. (0->54)
t30 = sin(qJ(3));
t33 = cos(qJ(3));
t50 = cos(pkin(5));
t27 = sin(pkin(5));
t31 = sin(qJ(2));
t55 = t27 * t31;
t74 = -t30 * t55 + t50 * t33;
t34 = cos(qJ(2));
t26 = sin(pkin(9));
t46 = t26 * t50;
t49 = cos(pkin(9));
t15 = -t31 * t46 + t49 * t34;
t54 = t27 * t33;
t73 = -t15 * t30 + t26 * t54;
t39 = t50 * t49;
t13 = t26 * t34 + t31 * t39;
t45 = t27 * t49;
t36 = -t13 * t30 - t33 * t45;
t72 = pkin(3) * t36;
t53 = t27 * t34;
t71 = g(3) * t53;
t61 = rSges(6,3) + pkin(8);
t29 = sin(qJ(5));
t32 = cos(qJ(5));
t70 = t29 * rSges(6,1) + t32 * rSges(6,2);
t14 = t49 * t31 + t34 * t46;
t12 = t26 * t31 - t34 * t39;
t65 = g(2) * t12;
t69 = g(1) * t14 + t65;
t25 = qJ(3) + pkin(10);
t23 = sin(t25);
t24 = cos(t25);
t37 = t32 * rSges(6,1) - t29 * rSges(6,2) + pkin(4);
t68 = t61 * t23 + t37 * t24;
t67 = -m(5) - m(6);
t22 = t33 * pkin(3) + pkin(2);
t64 = t22 * t71;
t63 = g(3) * t27;
t62 = rSges(4,3) + pkin(7);
t28 = -qJ(4) - pkin(7);
t60 = -t12 * t22 - t13 * t28;
t59 = -t14 * t22 - t15 * t28;
t56 = t26 * t27;
t43 = t73 * pkin(3);
t41 = rSges(5,1) * t24 - rSges(5,2) * t23;
t40 = t74 * pkin(3);
t38 = rSges(4,1) * t33 - rSges(4,2) * t30 + pkin(2);
t9 = t50 * t23 + t24 * t55;
t8 = -t23 * t55 + t50 * t24;
t5 = t15 * t24 + t23 * t56;
t4 = -t15 * t23 + t24 * t56;
t3 = t13 * t24 - t23 * t45;
t2 = -t13 * t23 - t24 * t45;
t1 = [(-m(2) - m(3) - m(4) + t67) * g(3), -m(3) * (g(1) * (-t14 * rSges(3,1) - t15 * rSges(3,2)) + g(2) * (-t12 * rSges(3,1) - t13 * rSges(3,2)) + (rSges(3,1) * t34 - rSges(3,2) * t31) * t63) - m(4) * (g(1) * (-t38 * t14 + t62 * t15) + g(2) * t62 * t13 - t38 * t65 + (t62 * t31 + t38 * t34) * t63) - m(5) * (g(1) * (t15 * rSges(5,3) - t41 * t14 + t59) + g(2) * (t13 * rSges(5,3) - t41 * t12 + t60) + t64 + (t41 * t34 + (rSges(5,3) - t28) * t31) * t63) - m(6) * (g(1) * (t70 * t15 + t59) + g(2) * (t70 * t13 + t60) + t64 + ((-t28 + t70) * t31 + t68 * t34) * t63 - t69 * t68), -m(4) * (g(1) * (t73 * rSges(4,1) + (-t15 * t33 - t30 * t56) * rSges(4,2)) + g(2) * (t36 * rSges(4,1) + (-t13 * t33 + t30 * t45) * rSges(4,2)) + g(3) * (t74 * rSges(4,1) + (-t50 * t30 - t31 * t54) * rSges(4,2))) - m(5) * (g(1) * (t4 * rSges(5,1) - t5 * rSges(5,2) + t43) + g(2) * (t2 * rSges(5,1) - t3 * rSges(5,2) + t72) + g(3) * (t8 * rSges(5,1) - t9 * rSges(5,2) + t40)) + (-g(1) * (t61 * t5 + t43) - g(2) * (t61 * t3 + t72) - g(3) * (t61 * t9 + t40) - (g(1) * t4 + g(2) * t2 + g(3) * t8) * t37) * m(6), t67 * (t69 - t71), -m(6) * (g(1) * ((t14 * t32 - t5 * t29) * rSges(6,1) + (-t14 * t29 - t5 * t32) * rSges(6,2)) + g(2) * ((t12 * t32 - t3 * t29) * rSges(6,1) + (-t12 * t29 - t3 * t32) * rSges(6,2)) + g(3) * ((-t9 * t29 - t32 * t53) * rSges(6,1) + (t29 * t53 - t9 * t32) * rSges(6,2)))];
taug = t1(:);
