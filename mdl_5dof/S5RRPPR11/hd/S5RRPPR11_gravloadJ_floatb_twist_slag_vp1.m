% Calculate Gravitation load on the joints for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:26
% EndTime: 2019-12-31 19:46:28
% DurationCPUTime: 0.58s
% Computational Cost: add. (183->102), mult. (319->142), div. (0->0), fcn. (291->8), ass. (0->50)
t45 = rSges(6,3) + pkin(7) + qJ(4);
t26 = cos(qJ(1));
t24 = sin(qJ(1));
t57 = g(2) * t24;
t33 = g(1) * t26 + t57;
t41 = rSges(5,3) + qJ(4);
t59 = -m(5) - m(6);
t25 = cos(qJ(2));
t56 = g(3) * t25;
t16 = t25 * pkin(2);
t21 = cos(pkin(8));
t54 = rSges(5,2) * t21;
t20 = sin(pkin(8));
t23 = sin(qJ(2));
t53 = t20 * t23;
t52 = t23 * t26;
t19 = pkin(8) + qJ(5);
t12 = sin(t19);
t51 = t24 * t12;
t13 = cos(t19);
t50 = t24 * t13;
t49 = t25 * rSges(4,2);
t48 = t26 * t12;
t47 = t26 * t13;
t46 = t26 * t25;
t14 = t23 * qJ(3);
t44 = t14 + t16;
t43 = t26 * pkin(1) + t24 * pkin(6);
t42 = qJ(3) * t25;
t40 = pkin(4) * t53;
t39 = -pkin(2) - t45;
t38 = -pkin(2) - t41;
t37 = -pkin(1) - t14;
t36 = pkin(2) * t46 + t26 * t14 + t43;
t35 = g(1) * t39;
t34 = g(1) * t38;
t32 = t25 * rSges(3,1) - t23 * rSges(3,2);
t30 = rSges(5,1) * t20 + t54;
t29 = t21 * rSges(5,1) - t20 * rSges(5,2) + pkin(3);
t28 = rSges(6,1) * t12 + rSges(6,2) * t13 + pkin(4) * t20;
t7 = t24 * t42;
t9 = t26 * t42;
t27 = g(1) * t9 + g(2) * t7 + g(3) * t44;
t17 = t26 * pkin(6);
t11 = t21 * pkin(4) + pkin(3);
t4 = -t23 * t51 + t47;
t3 = t23 * t50 + t48;
t2 = t23 * t48 + t50;
t1 = t23 * t47 - t51;
t5 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - t26 * rSges(2,2)) + g(2) * (t26 * rSges(2,1) - t24 * rSges(2,2))) - m(3) * (g(1) * (t26 * rSges(3,3) + t17) + g(2) * (rSges(3,1) * t46 - rSges(3,2) * t52 + t43) + (g(1) * (-pkin(1) - t32) + g(2) * rSges(3,3)) * t24) - m(4) * (g(1) * (t26 * rSges(4,1) + t17) + g(2) * (-rSges(4,2) * t46 + rSges(4,3) * t52 + t36) + (g(1) * (-t23 * rSges(4,3) - t16 + t37 + t49) + g(2) * rSges(4,1)) * t24) - m(5) * (g(1) * t17 + g(2) * t36 + (g(1) * t29 + g(2) * (rSges(5,1) * t53 + t23 * t54 + t41 * t25)) * t26 + (g(2) * t29 + t25 * t34 + (-pkin(1) + (-qJ(3) - t30) * t23) * g(1)) * t24) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t17) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t36) + (g(1) * t11 + g(2) * (t45 * t25 + t40)) * t26 + (g(1) * (t37 - t40) + g(2) * t11 + t25 * t35) * t24), -m(3) * (g(3) * t32 + t33 * (-rSges(3,1) * t23 - rSges(3,2) * t25)) - m(4) * (g(1) * (rSges(4,3) * t46 + t9) + g(2) * (t24 * t25 * rSges(4,3) + t7) + g(3) * (t44 - t49) + (g(3) * rSges(4,3) + t33 * (rSges(4,2) - pkin(2))) * t23) - m(5) * ((g(3) * t41 + t33 * t30) * t25 + (g(3) * t30 + t26 * t34 + t38 * t57) * t23 + t27) - m(6) * ((g(3) * t45 + t33 * t28) * t25 + (g(3) * t28 + t26 * t35 + t39 * t57) * t23 + t27), (-m(4) + t59) * (t33 * t23 - t56), t59 * (g(3) * t23 + t33 * t25), -m(6) * (g(1) * (t1 * rSges(6,1) - t2 * rSges(6,2)) + g(2) * (t3 * rSges(6,1) + t4 * rSges(6,2)) + (-rSges(6,1) * t13 + rSges(6,2) * t12) * t56)];
taug = t5(:);
