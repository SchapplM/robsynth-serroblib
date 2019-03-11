% Calculate Gravitation load on the joints for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:04:42
% EndTime: 2019-03-09 05:04:44
% DurationCPUTime: 0.95s
% Computational Cost: add. (506->156), mult. (621->205), div. (0->0), fcn. (650->10), ass. (0->57)
t32 = qJ(1) + pkin(10);
t27 = sin(t32);
t28 = cos(t32);
t34 = sin(qJ(4));
t38 = cos(qJ(4));
t39 = cos(qJ(3));
t62 = t38 * t39;
t10 = t27 * t62 - t28 * t34;
t33 = sin(qJ(6));
t37 = cos(qJ(6));
t63 = t34 * t39;
t9 = t27 * t63 + t28 * t38;
t51 = t10 * t37 + t33 * t9;
t79 = -t10 * t33 + t9 * t37;
t80 = t79 * rSges(7,1) - t51 * rSges(7,2);
t72 = g(2) * t27;
t77 = g(1) * t28 + t72;
t76 = -m(6) - m(7);
t75 = -pkin(4) - pkin(5);
t74 = g(1) * t27;
t35 = sin(qJ(3));
t71 = g(3) * t35;
t36 = sin(qJ(1));
t70 = t36 * pkin(1);
t30 = t39 * pkin(3);
t69 = -rSges(6,1) - pkin(4);
t68 = -pkin(9) - rSges(7,3);
t66 = t27 * t39;
t65 = t28 * t35;
t64 = t28 * t39;
t61 = t35 * pkin(8) + t30;
t60 = rSges(6,3) + qJ(5);
t40 = cos(qJ(1));
t31 = t40 * pkin(1);
t59 = t28 * pkin(2) + t27 * pkin(7) + t31;
t58 = -pkin(2) - t30;
t57 = t28 * t68;
t56 = t28 * pkin(7) - t70;
t54 = pkin(4) * t62 + qJ(5) * t63 + t61;
t53 = pkin(3) * t64 + pkin(8) * t65 + t59;
t11 = -t27 * t38 + t28 * t63;
t12 = t27 * t34 + t28 * t62;
t2 = t11 * t37 - t12 * t33;
t3 = t11 * t33 + t12 * t37;
t52 = rSges(7,1) * t2 - rSges(7,2) * t3;
t50 = rSges(4,1) * t39 - rSges(4,2) * t35;
t46 = t33 * t34 + t37 * t38;
t47 = t33 * t38 - t34 * t37;
t48 = (-rSges(7,1) * t47 - rSges(7,2) * t46) * t35;
t45 = t12 * pkin(4) + t53;
t43 = -t10 * pkin(4) - t9 * qJ(5) + t56;
t21 = t35 * t38 * qJ(5);
t17 = pkin(8) * t64;
t15 = pkin(8) * t66;
t7 = t11 * pkin(4);
t5 = t9 * pkin(4);
t1 = [-m(2) * (g(1) * (-t36 * rSges(2,1) - rSges(2,2) * t40) + g(2) * (rSges(2,1) * t40 - t36 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t27 - rSges(3,2) * t28 - t70) + g(2) * (rSges(3,1) * t28 - rSges(3,2) * t27 + t31)) - m(4) * (g(1) * (rSges(4,3) * t28 + t56) + g(2) * (rSges(4,1) * t64 - rSges(4,2) * t65 + t59) + (g(1) * (-pkin(2) - t50) + g(2) * rSges(4,3)) * t27) - m(5) * (g(1) * (-rSges(5,1) * t10 + rSges(5,2) * t9 + t56) + g(2) * (rSges(5,1) * t12 - rSges(5,2) * t11 + rSges(5,3) * t65 + t53) + ((-rSges(5,3) - pkin(8)) * t35 + t58) * t74) - m(6) * (g(1) * (-rSges(6,1) * t10 - rSges(6,3) * t9 + t43) + g(2) * (rSges(6,1) * t12 + rSges(6,2) * t65 + t11 * t60 + t45) + ((-rSges(6,2) - pkin(8)) * t35 + t58) * t74) - m(7) * (g(1) * (-rSges(7,1) * t51 - rSges(7,2) * t79 - t10 * pkin(5) + t43) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + pkin(5) * t12 + qJ(5) * t11 + t35 * t57 + t45) + ((-pkin(8) - t68) * t35 + t58) * t74) (-m(3) - m(4) - m(5) + t76) * g(3), -m(4) * (g(3) * t50 + t77 * (-rSges(4,1) * t35 - rSges(4,2) * t39)) - m(5) * (g(1) * (rSges(5,3) * t64 + t17) + g(2) * (rSges(5,3) * t66 + t15) + g(3) * (rSges(5,1) * t62 - rSges(5,2) * t63 + t61) + (g(3) * rSges(5,3) + t77 * (-rSges(5,1) * t38 + rSges(5,2) * t34 - pkin(3))) * t35) - m(6) * (g(1) * (rSges(6,2) * t64 + t17) + g(2) * (rSges(6,2) * t66 + t15) + g(3) * (rSges(6,1) * t62 + rSges(6,3) * t63 + t54) + (g(3) * rSges(6,2) + t77 * (-t34 * t60 + t38 * t69 - pkin(3))) * t35) - m(7) * (g(1) * t17 + g(2) * t15 + g(3) * t54 + (g(3) * (rSges(7,1) * t46 - rSges(7,2) * t47 + t38 * pkin(5)) + g(1) * t57 + t68 * t72) * t39 + (g(3) * t68 + t77 * (-pkin(3) + (-t33 * rSges(7,1) - rSges(7,2) * t37 - qJ(5)) * t34 + (-t37 * rSges(7,1) + t33 * rSges(7,2) + t75) * t38)) * t35) -m(5) * (g(1) * (-rSges(5,1) * t11 - rSges(5,2) * t12) + g(2) * (-rSges(5,1) * t9 - rSges(5,2) * t10)) - m(6) * (g(1) * (-rSges(6,1) * t11 + t12 * t60 - t7) + g(2) * (-rSges(6,1) * t9 + t10 * t60 - t5) + g(3) * t21) - m(7) * (g(1) * (-t11 * pkin(5) + t12 * qJ(5) - t52 - t7) + g(2) * (-t9 * pkin(5) + t10 * qJ(5) - t5 - t80) + g(3) * (t21 - t48)) + ((m(5) * rSges(5,2) - m(6) * rSges(6,3)) * t38 + (m(5) * rSges(5,1) - m(6) * t69 - m(7) * t75) * t34) * t71, t76 * (g(1) * t11 + g(2) * t9 + t34 * t71) -m(7) * (g(1) * t52 + g(2) * t80 + g(3) * t48)];
taug  = t1(:);
