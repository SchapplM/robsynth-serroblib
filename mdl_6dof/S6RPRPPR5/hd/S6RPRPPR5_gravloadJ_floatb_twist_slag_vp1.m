% Calculate Gravitation load on the joints for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:50:05
% EndTime: 2019-03-09 02:50:07
% DurationCPUTime: 0.54s
% Computational Cost: add. (334->111), mult. (359->148), div. (0->0), fcn. (323->10), ass. (0->60)
t28 = -pkin(7) - qJ(2);
t70 = pkin(4) - t28;
t54 = rSges(7,3) + pkin(8) + qJ(5);
t30 = cos(qJ(1));
t29 = sin(qJ(1));
t67 = g(2) * t29;
t39 = g(1) * t30 + t67;
t49 = rSges(6,3) + qJ(5);
t69 = -m(6) - m(7);
t22 = pkin(9) + qJ(3);
t20 = cos(t22);
t66 = g(3) * t20;
t14 = t20 * pkin(3);
t25 = cos(pkin(10));
t64 = rSges(6,2) * t25;
t21 = pkin(10) + qJ(6);
t17 = sin(t21);
t63 = t17 * t30;
t18 = sin(t22);
t23 = sin(pkin(10));
t62 = t18 * t23;
t19 = cos(t21);
t61 = t19 * t30;
t60 = t20 * rSges(5,2);
t59 = t20 * t30;
t58 = t29 * t17;
t57 = t29 * t19;
t56 = rSges(5,1) - t28;
t55 = rSges(4,3) - t28;
t13 = t18 * qJ(4);
t53 = t13 + t14;
t52 = pkin(5) * t25 + t70;
t51 = qJ(4) * t20;
t50 = rSges(3,3) + qJ(2);
t48 = -m(5) + t69;
t47 = pkin(5) * t62;
t46 = -pkin(3) - t54;
t26 = cos(pkin(9));
t16 = pkin(2) * t26 + pkin(1);
t10 = t30 * t16;
t45 = pkin(3) * t59 + t30 * t13 + t10;
t44 = -pkin(3) - t49;
t43 = -t16 - t13;
t42 = g(1) * t46;
t41 = g(2) * t45;
t40 = g(1) * t44;
t38 = rSges(4,1) * t20 - rSges(4,2) * t18;
t36 = rSges(6,1) * t23 + t64;
t35 = t18 * rSges(5,3) - t60;
t34 = rSges(3,1) * t26 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t33 = rSges(6,1) * t25 - rSges(6,2) * t23 + t70;
t32 = rSges(7,1) * t17 + rSges(7,2) * t19 + pkin(5) * t23;
t7 = t29 * t51;
t9 = t30 * t51;
t31 = g(1) * t9 + g(2) * t7 + g(3) * t53;
t6 = -t18 * t58 + t61;
t5 = t18 * t57 + t63;
t4 = t18 * t63 + t57;
t3 = t18 * t61 - t58;
t1 = [-m(2) * (g(1) * (-t29 * rSges(2,1) - rSges(2,2) * t30) + g(2) * (rSges(2,1) * t30 - t29 * rSges(2,2))) - m(3) * ((g(1) * t50 + g(2) * t34) * t30 + (-g(1) * t34 + g(2) * t50) * t29) - m(4) * (g(2) * t10 + (g(1) * t55 + g(2) * t38) * t30 + (g(1) * (-t16 - t38) + g(2) * t55) * t29) - m(5) * (t41 + (g(1) * t56 + g(2) * t35) * t30 + (g(1) * (-t35 + t43 - t14) + g(2) * t56) * t29) - m(6) * (t41 + (g(1) * t33 + g(2) * (rSges(6,1) * t62 + t18 * t64 + t49 * t20)) * t30 + (g(2) * t33 + t20 * t40 + (-t16 + (-qJ(4) - t36) * t18) * g(1)) * t29) - m(7) * (g(1) * (t6 * rSges(7,1) - t5 * rSges(7,2)) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t45) + (g(1) * t52 + g(2) * (t54 * t20 + t47)) * t30 + (g(1) * (t43 - t47) + g(2) * t52 + t20 * t42) * t29) (-m(3) - m(4) + t48) * (g(1) * t29 - g(2) * t30) -m(4) * (g(3) * t38 + t39 * (-rSges(4,1) * t18 - rSges(4,2) * t20)) - m(5) * (g(1) * (rSges(5,3) * t59 + t9) + g(2) * (t29 * t20 * rSges(5,3) + t7) + g(3) * (t53 - t60) + (g(3) * rSges(5,3) + t39 * (rSges(5,2) - pkin(3))) * t18) - m(6) * ((g(3) * t49 + t39 * t36) * t20 + (g(3) * t36 + t30 * t40 + t44 * t67) * t18 + t31) - m(7) * ((g(3) * t54 + t39 * t32) * t20 + (g(3) * t32 + t30 * t42 + t46 * t67) * t18 + t31) t48 * (t39 * t18 - t66) t69 * (g(3) * t18 + t39 * t20) -m(7) * (g(1) * (rSges(7,1) * t3 - rSges(7,2) * t4) + g(2) * (rSges(7,1) * t5 + rSges(7,2) * t6) + (-rSges(7,1) * t19 + rSges(7,2) * t17) * t66)];
taug  = t1(:);
