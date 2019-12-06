% Calculate Gravitation load on the joints for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:40
% EndTime: 2019-12-05 15:49:44
% DurationCPUTime: 0.64s
% Computational Cost: add. (332->99), mult. (839->166), div. (0->0), fcn. (1035->12), ass. (0->50)
t29 = sin(pkin(10));
t36 = sin(qJ(2));
t39 = cos(qJ(2));
t53 = cos(pkin(10));
t23 = -t39 * t29 - t36 * t53;
t30 = sin(pkin(9));
t32 = cos(pkin(9));
t33 = cos(pkin(5));
t42 = -t36 * t29 + t39 * t53;
t40 = t42 * t33;
t11 = t32 * t23 - t30 * t40;
t58 = t33 * t39;
t43 = -t30 * t58 - t32 * t36;
t41 = t43 * pkin(2);
t68 = t11 * pkin(3) + t41;
t67 = -t30 * t36 + t32 * t58;
t38 = cos(qJ(4));
t66 = t38 * pkin(4);
t65 = rSges(5,3) + pkin(7);
t64 = rSges(6,3) + pkin(8);
t31 = sin(pkin(5));
t35 = sin(qJ(4));
t62 = t31 * t35;
t61 = t31 * t38;
t59 = t33 * t36;
t34 = sin(qJ(5));
t57 = t34 * t38;
t37 = cos(qJ(5));
t55 = t37 * t38;
t19 = t42 * t31;
t27 = t31 * t39 * pkin(2);
t54 = t19 * pkin(3) + t27;
t52 = -m(4) - m(5) - m(6);
t49 = g(2) * t64;
t47 = t67 * pkin(2);
t46 = rSges(5,1) * t38 - rSges(5,2) * t35;
t21 = t23 * t33;
t9 = -t32 * t21 + t30 * t42;
t10 = -t30 * t21 - t32 * t42;
t8 = t30 * t23 + t32 * t40;
t45 = t8 * pkin(3) + t47;
t44 = rSges(6,1) * t37 - rSges(6,2) * t34 + pkin(4);
t20 = t23 * t31;
t14 = -t20 * t38 + t33 * t35;
t13 = t20 * t35 + t33 * t38;
t4 = -t10 * t38 + t30 * t62;
t3 = t10 * t35 + t30 * t61;
t2 = -t32 * t62 + t9 * t38;
t1 = -t32 * t61 - t9 * t35;
t5 = [(-m(2) - m(3) + t52) * g(3), -m(3) * (g(1) * (t43 * rSges(3,1) + (t30 * t59 - t32 * t39) * rSges(3,2)) + g(2) * (t67 * rSges(3,1) + (-t30 * t39 - t32 * t59) * rSges(3,2)) + g(3) * (rSges(3,1) * t39 - rSges(3,2) * t36) * t31) - m(4) * (g(1) * (t11 * rSges(4,1) + t10 * rSges(4,2) + t41) + g(2) * (t8 * rSges(4,1) - rSges(4,2) * t9 + t47) + g(3) * (t19 * rSges(4,1) + t20 * rSges(4,2) + t27)) - m(5) * (g(1) * (-t65 * t10 + t46 * t11 + t68) + g(2) * (t46 * t8 + t65 * t9 + t45) + g(3) * (t46 * t19 - t65 * t20 + t54)) - m(6) * (g(1) * (t11 * t66 - t10 * pkin(7) + (-t10 * t34 + t11 * t55) * rSges(6,1) + (-t10 * t37 - t11 * t57) * rSges(6,2) + t68) + g(2) * (t8 * t66 + t9 * pkin(7) + (t34 * t9 + t8 * t55) * rSges(6,1) + (t37 * t9 - t8 * t57) * rSges(6,2) + t45) + g(3) * (t19 * t66 - t20 * pkin(7) + (t19 * t55 - t20 * t34) * rSges(6,1) + (-t19 * t57 - t20 * t37) * rSges(6,2) + t54) + (t8 * t49 + (g(1) * t11 + g(3) * t19) * t64) * t35), t52 * (g(3) * t33 + (g(1) * t30 - g(2) * t32) * t31), -m(5) * (g(1) * (t3 * rSges(5,1) - t4 * rSges(5,2)) + g(2) * (t1 * rSges(5,1) - t2 * rSges(5,2)) + g(3) * (t13 * rSges(5,1) - t14 * rSges(5,2))) - m(6) * (g(1) * (t44 * t3 + t64 * t4) + t2 * t49 + g(2) * t44 * t1 + (t44 * t13 + t64 * t14) * g(3)), -m(6) * (g(1) * ((-t11 * t37 - t4 * t34) * rSges(6,1) + (t11 * t34 - t4 * t37) * rSges(6,2)) + g(2) * ((-t2 * t34 - t8 * t37) * rSges(6,1) + (-t2 * t37 + t8 * t34) * rSges(6,2)) + g(3) * ((-t14 * t34 - t19 * t37) * rSges(6,1) + (-t14 * t37 + t19 * t34) * rSges(6,2)))];
taug = t5(:);
