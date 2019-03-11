% Calculate Gravitation load on the joints for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:46
% EndTime: 2019-03-09 04:38:48
% DurationCPUTime: 0.74s
% Computational Cost: add. (460->129), mult. (477->171), div. (0->0), fcn. (457->10), ass. (0->59)
t26 = pkin(9) + qJ(3);
t24 = cos(t26);
t33 = sin(qJ(1));
t34 = cos(qJ(4));
t57 = t33 * t34;
t32 = sin(qJ(4));
t35 = cos(qJ(1));
t61 = t32 * t35;
t8 = -t24 * t61 + t57;
t65 = rSges(7,1) + pkin(5);
t64 = rSges(5,3) + pkin(8);
t74 = g(1) * t35 + g(2) * t33;
t50 = rSges(7,3) + qJ(6);
t27 = qJ(4) + pkin(10);
t23 = sin(t27);
t25 = cos(t27);
t73 = t50 * t23 + t65 * t25;
t72 = -m(6) - m(7);
t71 = pkin(4) * t32;
t31 = -pkin(7) - qJ(2);
t69 = g(2) * t31;
t21 = pkin(4) * t34 + pkin(3);
t11 = t24 * t21;
t67 = g(3) * t11;
t22 = sin(t26);
t66 = g(3) * t22;
t63 = t22 * t35;
t62 = t25 * t35;
t60 = t33 * t23;
t59 = t33 * t25;
t58 = t33 * t32;
t56 = t34 * t35;
t55 = t35 * t23;
t30 = -qJ(5) - pkin(8);
t54 = rSges(7,2) - t30;
t53 = rSges(4,3) - t31;
t52 = rSges(6,3) - t30;
t51 = rSges(3,3) + qJ(2);
t29 = cos(pkin(9));
t20 = pkin(2) * t29 + pkin(1);
t48 = -t20 - t11;
t47 = t33 * t22 * t30 + pkin(4) * t61 - t31 * t35;
t46 = rSges(4,1) * t24 - rSges(4,2) * t22;
t44 = rSges(6,1) * t25 - rSges(6,2) * t23;
t43 = t8 * pkin(4);
t42 = rSges(3,1) * t29 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t41 = rSges(5,1) * t34 - rSges(5,2) * t32 + pkin(3);
t6 = t24 * t58 + t56;
t13 = t35 * t20;
t39 = pkin(4) * t58 + t35 * t11 - t30 * t63 + t13;
t38 = pkin(3) * t24 + t64 * t22;
t37 = t6 * pkin(4);
t9 = t24 * t56 + t58;
t7 = -t24 * t57 + t61;
t5 = t24 * t62 + t60;
t4 = t24 * t55 - t59;
t3 = t24 * t59 - t55;
t2 = t24 * t60 + t62;
t1 = [-m(2) * (g(1) * (-t33 * rSges(2,1) - rSges(2,2) * t35) + g(2) * (rSges(2,1) * t35 - t33 * rSges(2,2))) - m(3) * ((g(1) * t51 + g(2) * t42) * t35 + (-g(1) * t42 + g(2) * t51) * t33) - m(4) * (g(2) * t13 + (g(1) * t53 + g(2) * t46) * t35 + (g(1) * (-t20 - t46) + g(2) * t53) * t33) - m(5) * (g(1) * (t7 * rSges(5,1) + t6 * rSges(5,2)) + g(2) * (t9 * rSges(5,1) + t8 * rSges(5,2) + t13) + (-g(1) * t31 + g(2) * t38) * t35 + (g(1) * (-t20 - t38) - t69) * t33) - m(6) * (g(1) * (-t3 * rSges(6,1) + t2 * rSges(6,2) + t47) + g(2) * (t5 * rSges(6,1) - t4 * rSges(6,2) + rSges(6,3) * t63 + t39) + (g(1) * (-rSges(6,3) * t22 + t48) - t69) * t33) - m(7) * (g(1) * (-t50 * t2 - t65 * t3 + t47) + g(2) * (rSges(7,2) * t63 + t50 * t4 + t65 * t5 + t39) + (g(1) * (-rSges(7,2) * t22 + t48) - t69) * t33) (-m(3) - m(4) - m(5) + t72) * (g(1) * t33 - g(2) * t35) -m(4) * (g(3) * t46 + t74 * (-rSges(4,1) * t22 - rSges(4,2) * t24)) - m(5) * ((g(3) * t41 + t74 * t64) * t24 + (g(3) * t64 - t74 * t41) * t22) - m(6) * (t67 + (g(3) * t44 + t74 * t52) * t24 + (g(3) * t52 + t74 * (-t21 - t44)) * t22) - m(7) * (t67 + (g(3) * t73 + t74 * t54) * t24 + (g(3) * t54 + t74 * (-t21 - t73)) * t22) -m(5) * (g(1) * (rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (-rSges(5,1) * t6 + rSges(5,2) * t7)) - m(6) * (g(1) * (-t4 * rSges(6,1) - t5 * rSges(6,2) + t43) + g(2) * (-t2 * rSges(6,1) - t3 * rSges(6,2) - t37)) - m(7) * (g(1) * (-t65 * t4 + t50 * t5 + t43) + g(2) * (-t65 * t2 + t50 * t3 - t37)) + (-m(5) * (-rSges(5,1) * t32 - rSges(5,2) * t34) - m(6) * (-rSges(6,1) * t23 - rSges(6,2) * t25 - t71) - m(7) * (-t65 * t23 + t50 * t25 - t71)) * t66, t72 * (-g(3) * t24 + t74 * t22) -m(7) * (g(1) * t4 + g(2) * t2 + t23 * t66)];
taug  = t1(:);
