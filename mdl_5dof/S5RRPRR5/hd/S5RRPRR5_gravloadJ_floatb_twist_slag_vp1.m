% Calculate Gravitation load on the joints for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:33:37
% EndTime: 2019-12-05 18:33:38
% DurationCPUTime: 0.27s
% Computational Cost: add. (276->71), mult. (199->89), div. (0->0), fcn. (157->10), ass. (0->46)
t34 = -pkin(7) - qJ(3);
t66 = rSges(5,3) - t34;
t65 = rSges(6,3) + pkin(8) - t34;
t64 = rSges(4,3) + qJ(3);
t33 = cos(pkin(9));
t63 = -rSges(4,1) * t33 - pkin(2);
t31 = qJ(1) + qJ(2);
t26 = sin(t31);
t30 = pkin(9) + qJ(4);
t25 = qJ(5) + t30;
t20 = sin(t25);
t52 = t20 * t26;
t21 = cos(t25);
t53 = rSges(6,2) * t21;
t62 = rSges(6,1) * t52 + t26 * t53;
t27 = cos(t31);
t61 = g(3) * t27;
t35 = sin(qJ(1));
t60 = t35 * pkin(1);
t36 = cos(qJ(1));
t59 = t36 * pkin(1);
t22 = t33 * pkin(3) + pkin(2);
t24 = cos(t30);
t57 = rSges(5,1) * t24;
t56 = rSges(4,2) * sin(pkin(9));
t23 = sin(t30);
t55 = rSges(5,2) * t23;
t54 = rSges(6,2) * t20;
t13 = t21 * rSges(6,1);
t51 = t23 * t26;
t50 = t26 * t56 + t64 * t27;
t15 = pkin(4) * t24;
t49 = -t15 - t22 - t13;
t48 = -t27 * rSges(3,1) + t26 * rSges(3,2);
t47 = -t22 - t57;
t46 = t13 - t54;
t45 = -rSges(3,1) * t26 - rSges(3,2) * t27;
t44 = rSges(5,1) * t23 + rSges(5,2) * t24;
t43 = -rSges(6,1) * t20 - t53;
t42 = (t56 + t63) * t27;
t41 = (t49 + t54) * t27 - t65 * t26;
t40 = rSges(6,2) * t52 + t49 * t26 + t65 * t27;
t39 = (t47 + t55) * t27 - t66 * t26;
t38 = rSges(5,2) * t51 + t47 * t26 + t66 * t27;
t37 = (-g(2) * t64 + g(3) * t63) * t26;
t1 = [-m(2) * (g(2) * (-rSges(2,1) * t36 + t35 * rSges(2,2)) + g(3) * (-t35 * rSges(2,1) - rSges(2,2) * t36)) - m(3) * (g(2) * (t48 - t59) + g(3) * (t45 - t60)) - m(4) * (g(2) * (t42 - t59) + g(3) * (t50 - t60) + t37) - m(5) * (g(2) * (t39 - t59) + g(3) * (t38 - t60)) - m(6) * (g(2) * (t41 - t59) + g(3) * (t40 - t60)), -m(3) * (g(2) * t48 + g(3) * t45) - m(4) * (g(2) * t42 + g(3) * t50 + t37) - m(5) * (g(2) * t39 + g(3) * t38) - m(6) * (g(2) * t41 + g(3) * t40), (-m(4) - m(5) - m(6)) * (g(2) * t27 + g(3) * t26), -m(5) * (g(1) * (-t55 + t57) + g(2) * t44 * t26) - m(6) * (g(1) * (t15 + t46) + g(2) * (pkin(4) * t51 + t62)) + (m(5) * t44 - m(6) * (-pkin(4) * t23 + t43)) * t61, -m(6) * (g(1) * t46 + g(2) * t62 + t43 * t61)];
taug = t1(:);
