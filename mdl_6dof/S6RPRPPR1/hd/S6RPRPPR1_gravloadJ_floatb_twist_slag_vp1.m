% Calculate Gravitation load on the joints for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:38:00
% EndTime: 2019-03-09 02:38:00
% DurationCPUTime: 0.50s
% Computational Cost: add. (381->97), mult. (306->123), div. (0->0), fcn. (270->12), ass. (0->49)
t19 = qJ(3) + pkin(10);
t11 = sin(t19);
t14 = cos(t19);
t21 = sin(pkin(11));
t22 = cos(pkin(11));
t37 = rSges(6,1) * t22 - rSges(6,2) * t21 + pkin(4);
t46 = rSges(6,3) + qJ(5);
t59 = t46 * t11 + t37 * t14;
t47 = rSges(7,3) + pkin(8) + qJ(5);
t60 = t47 * t11;
t20 = qJ(1) + pkin(9);
t12 = sin(t20);
t15 = cos(t20);
t58 = g(1) * t15 + g(2) * t12;
t57 = -m(6) - m(7);
t25 = sin(qJ(3));
t56 = pkin(3) * t25;
t26 = sin(qJ(1));
t53 = t26 * pkin(1);
t52 = rSges(4,3) + pkin(7);
t28 = cos(qJ(1));
t17 = t28 * pkin(1);
t27 = cos(qJ(3));
t16 = t27 * pkin(3);
t9 = t16 + pkin(2);
t51 = t15 * t9 + t17;
t50 = t12 * t14;
t49 = t15 * t14;
t23 = -qJ(4) - pkin(7);
t48 = rSges(5,3) - t23;
t45 = -m(5) + t57;
t44 = g(1) * t53;
t43 = pkin(5) * t21 - t23;
t41 = t27 * rSges(4,1) - t25 * rSges(4,2);
t39 = t14 * rSges(5,1) - t11 * rSges(5,2);
t38 = pkin(2) + t41;
t18 = pkin(11) + qJ(6);
t10 = sin(t18);
t13 = cos(t18);
t8 = t22 * pkin(5) + pkin(4);
t36 = rSges(7,1) * t13 - rSges(7,2) * t10 + t8;
t35 = t21 * rSges(6,1) + t22 * rSges(6,2) - t23;
t34 = g(2) * t51 - t44;
t32 = t14 * t8 + t60;
t5 = t12 * t10 + t13 * t49;
t4 = -t10 * t49 + t12 * t13;
t3 = t15 * t10 - t13 * t50;
t2 = t10 * t50 + t15 * t13;
t1 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - t28 * rSges(2,2)) + g(2) * (t28 * rSges(2,1) - t26 * rSges(2,2))) - m(3) * (g(1) * (-t12 * rSges(3,1) - t15 * rSges(3,2) - t53) + g(2) * (t15 * rSges(3,1) - t12 * rSges(3,2) + t17)) - m(4) * (-t44 + g(2) * t17 + (g(1) * t52 + g(2) * t38) * t15 + (-g(1) * t38 + g(2) * t52) * t12) - m(5) * ((g(1) * t48 + g(2) * t39) * t15 + (g(1) * (-t39 - t9) + g(2) * t48) * t12 + t34) - m(6) * ((g(1) * t35 + t59 * g(2)) * t15 + (g(2) * t35 + (-t9 - t59) * g(1)) * t12 + t34) - m(7) * (g(1) * (t3 * rSges(7,1) + t2 * rSges(7,2) - t53) + g(2) * (t5 * rSges(7,1) + t4 * rSges(7,2) + t51) + (g(1) * t43 + g(2) * t32) * t15 + (g(1) * (-t32 - t9) + g(2) * t43) * t12) (-m(3) - m(4) + t45) * g(3), -m(4) * (g(3) * t41 + t58 * (-rSges(4,1) * t25 - rSges(4,2) * t27)) - m(5) * (g(3) * (t16 + t39) + t58 * (-rSges(5,1) * t11 - rSges(5,2) * t14 - t56)) - m(6) * (g(3) * (t16 + t59) + t58 * (-t37 * t11 + t46 * t14 - t56)) - m(7) * (g(3) * (t36 * t14 + t16 + t60) + t58 * (-t36 * t11 + t47 * t14 - t56)) t45 * (g(1) * t12 - g(2) * t15) t57 * (-g(3) * t14 + t58 * t11) -m(7) * (g(1) * (t4 * rSges(7,1) - t5 * rSges(7,2)) + g(2) * (-t2 * rSges(7,1) + t3 * rSges(7,2)) + g(3) * (-rSges(7,1) * t10 - rSges(7,2) * t13) * t11)];
taug  = t1(:);
