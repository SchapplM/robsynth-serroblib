% Calculate Gravitation load on the joints for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:18
% EndTime: 2019-12-05 18:42:19
% DurationCPUTime: 0.32s
% Computational Cost: add. (290->75), mult. (218->91), div. (0->0), fcn. (173->10), ass. (0->48)
t31 = qJ(3) + pkin(9);
t24 = sin(t31);
t25 = cos(t31);
t36 = cos(qJ(3));
t29 = t36 * pkin(3);
t73 = t25 * rSges(5,1) - t24 * rSges(5,2) + t29;
t72 = rSges(4,3) + pkin(7);
t33 = -qJ(4) - pkin(7);
t71 = rSges(5,3) - t33;
t70 = rSges(6,3) + pkin(8) - t33;
t69 = -pkin(2) - t73;
t55 = t36 * rSges(4,1);
t68 = -pkin(2) - t55;
t34 = sin(qJ(3));
t64 = t34 * pkin(3);
t66 = m(4) * (rSges(4,1) * t34 + rSges(4,2) * t36) + m(5) * (rSges(5,1) * t24 + rSges(5,2) * t25 + t64);
t32 = qJ(1) + qJ(2);
t28 = cos(t32);
t65 = g(3) * t28;
t35 = sin(qJ(1));
t63 = t35 * pkin(1);
t37 = cos(qJ(1));
t62 = t37 * pkin(1);
t26 = qJ(5) + t31;
t21 = cos(t26);
t61 = rSges(6,2) * t21;
t20 = sin(t26);
t60 = t20 * rSges(6,2);
t27 = sin(t32);
t59 = t20 * t27;
t14 = t21 * rSges(6,1);
t56 = t34 * rSges(4,2);
t54 = pkin(4) * t25 + t29;
t53 = t27 * t56 + t72 * t28;
t52 = g(2) * (rSges(6,1) * t59 + t27 * t61);
t51 = -pkin(2) - t54 - t14;
t50 = -t28 * rSges(3,1) + t27 * rSges(3,2);
t48 = t14 - t60;
t47 = -t27 * rSges(3,1) - t28 * rSges(3,2);
t45 = -rSges(6,1) * t20 - t61;
t44 = (t56 + t68) * t28;
t42 = (t51 + t60) * t28 - t70 * t27;
t41 = rSges(6,2) * t59 + t51 * t27 + t70 * t28;
t40 = -t71 * t27 + t69 * t28;
t39 = t69 * t27 + t71 * t28;
t38 = (-g(2) * t72 + g(3) * t68) * t27;
t7 = -pkin(4) * t24 - t64;
t1 = [-m(2) * (g(2) * (-t37 * rSges(2,1) + t35 * rSges(2,2)) + g(3) * (-t35 * rSges(2,1) - t37 * rSges(2,2))) - m(3) * (g(2) * (t50 - t62) + g(3) * (t47 - t63)) - m(4) * (g(2) * (t44 - t62) + g(3) * (t53 - t63) + t38) - m(5) * (g(2) * (t40 - t62) + g(3) * (t39 - t63)) - m(6) * (g(2) * (t42 - t62) + g(3) * (t41 - t63)), -m(3) * (g(2) * t50 + g(3) * t47) - m(4) * (g(2) * t44 + g(3) * t53 + t38) - m(5) * (g(2) * t40 + g(3) * t39) - m(6) * (g(2) * t42 + g(3) * t41), -m(6) * t52 + (-m(4) * (t55 - t56) - m(5) * t73 - m(6) * (t48 + t54)) * g(1) + (m(6) * t7 - t66) * g(2) * t27 + (-m(6) * (t45 + t7) + t66) * t65, (-m(5) - m(6)) * (g(2) * t28 + g(3) * t27), -m(6) * (g(1) * t48 + t45 * t65 + t52)];
taug = t1(:);
