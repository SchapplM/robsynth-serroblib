% Calculate Gravitation load on the joints for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:57:46
% EndTime: 2019-03-09 03:57:48
% DurationCPUTime: 0.61s
% Computational Cost: add. (312->121), mult. (372->158), div. (0->0), fcn. (338->10), ass. (0->64)
t27 = qJ(3) + pkin(10);
t21 = cos(t27);
t59 = rSges(6,3) + pkin(8);
t78 = t21 * t59;
t50 = rSges(7,3) + pkin(9) + pkin(8);
t35 = cos(qJ(1));
t64 = g(2) * t35;
t32 = sin(qJ(1));
t67 = g(1) * t32;
t77 = -t64 + t67;
t76 = g(1) * t35 + g(2) * t32;
t20 = sin(t27);
t33 = cos(qJ(5));
t19 = pkin(5) * t33 + pkin(4);
t28 = qJ(5) + qJ(6);
t22 = sin(t28);
t23 = cos(t28);
t38 = rSges(7,1) * t23 - rSges(7,2) * t22 + t19;
t75 = t50 * t20 + t38 * t21;
t30 = sin(qJ(5));
t39 = rSges(6,1) * t33 - rSges(6,2) * t30 + pkin(4);
t74 = t59 * t20 + t39 * t21;
t55 = t23 * t35;
t58 = t22 * t32;
t5 = -t20 * t58 + t55;
t56 = t23 * t32;
t57 = t22 * t35;
t6 = t20 * t56 + t57;
t73 = t5 * rSges(7,1) - t6 * rSges(7,2);
t7 = t20 * t57 + t56;
t8 = t20 * t55 - t58;
t72 = t7 * rSges(7,1) + t8 * rSges(7,2);
t34 = cos(qJ(3));
t71 = pkin(3) * t34;
t70 = pkin(4) * t20;
t69 = pkin(5) * t30;
t16 = t32 * t71;
t68 = g(1) * t16;
t63 = g(3) * t20;
t62 = g(3) * t21;
t31 = sin(qJ(3));
t61 = t31 * pkin(3);
t60 = rSges(4,3) + pkin(7);
t54 = t30 * t32;
t53 = t30 * t35;
t52 = t32 * t33;
t51 = t33 * t35;
t49 = t35 * pkin(1) + t32 * qJ(2);
t48 = -m(5) - m(6) - m(7);
t47 = t32 * t61 + t49;
t25 = t35 * qJ(2);
t29 = -qJ(4) - pkin(7);
t46 = t32 * t29 + t35 * t61 + t25;
t45 = t50 * t21;
t43 = rSges(4,1) * t31 + rSges(4,2) * t34;
t42 = rSges(5,1) * t21 - rSges(5,2) * t20;
t41 = rSges(5,1) * t20 + rSges(5,2) * t21;
t40 = -rSges(7,1) * t22 - rSges(7,2) * t23;
t11 = t20 * t53 + t52;
t9 = -t20 * t54 + t51;
t37 = t20 * t19 - t45;
t12 = t20 * t51 - t54;
t10 = t20 * t52 + t53;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t32 - rSges(2,2) * t35) + g(2) * (rSges(2,1) * t35 - rSges(2,2) * t32)) - m(3) * (g(1) * (rSges(3,3) * t35 + t25 + (rSges(3,2) - pkin(1)) * t32) + g(2) * (-rSges(3,2) * t35 + rSges(3,3) * t32 + t49)) - m(4) * (g(1) * t25 + g(2) * t49 + (g(1) * t43 + g(2) * t60) * t35 + (g(1) * (-pkin(1) - t60) + g(2) * t43) * t32) - m(5) * (g(1) * t46 + g(2) * t47 + (g(1) * t41 + g(2) * (rSges(5,3) - t29)) * t35 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t41) * t32) - m(6) * (g(1) * (rSges(6,1) * t12 - rSges(6,2) * t11 - t32 * pkin(1) + t35 * t70 + t46) + g(2) * (rSges(6,1) * t10 + rSges(6,2) * t9 - t35 * t29 + t32 * t70 + t47) - t76 * t78) - m(7) * (g(1) * (t8 * rSges(7,1) - t7 * rSges(7,2) + t46) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t47) + (g(1) * t37 + g(2) * (-t29 + t69)) * t35 + (g(1) * (-pkin(1) - t69) + g(2) * t37) * t32) (-m(3) - m(4) + t48) * t77, -m(4) * (-g(3) * t43 + t77 * (rSges(4,1) * t34 - rSges(4,2) * t31)) - m(5) * (g(1) * (t42 * t32 + t16) + g(3) * (-t41 - t61) + (-t42 - t71) * t64) - m(6) * (t68 + g(3) * (-t61 + t78) - t39 * t63 + t74 * t67 + (-t71 - t74) * t64) - m(7) * (t68 + g(3) * (t45 - t61) - t38 * t63 + t75 * t67 + (-t71 - t75) * t64) t48 * t76, -m(6) * (g(1) * (rSges(6,1) * t9 - rSges(6,2) * t10) + g(2) * (rSges(6,1) * t11 + rSges(6,2) * t12)) - m(7) * (g(1) * (t9 * pkin(5) + t73) + g(2) * (t11 * pkin(5) + t72)) + (-m(6) * (-rSges(6,1) * t30 - rSges(6,2) * t33) - m(7) * (t40 - t69)) * t62, -m(7) * (g(1) * t73 + g(2) * t72 + t40 * t62)];
taug  = t1(:);
