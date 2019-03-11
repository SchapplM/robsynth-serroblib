% Calculate Gravitation load on the joints for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:12:28
% EndTime: 2019-03-09 13:12:29
% DurationCPUTime: 0.61s
% Computational Cost: add. (569->116), mult. (412->145), div. (0->0), fcn. (354->12), ass. (0->68)
t98 = rSges(7,3) + pkin(10);
t39 = qJ(2) + pkin(11);
t35 = qJ(4) + t39;
t31 = qJ(5) + t35;
t26 = sin(t31);
t27 = cos(t31);
t41 = sin(qJ(6));
t85 = rSges(7,2) * t41;
t97 = t26 * t85 + t27 * t98;
t95 = t27 * rSges(6,1) - t26 * rSges(6,2);
t29 = sin(t35);
t30 = cos(t35);
t68 = t30 * rSges(5,1) - t29 * rSges(5,2);
t33 = sin(t39);
t34 = cos(t39);
t45 = cos(qJ(2));
t37 = t45 * pkin(2);
t94 = t34 * rSges(4,1) - t33 * rSges(4,2) + t37;
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t93 = g(1) * t46 + g(2) * t43;
t53 = t27 * pkin(5) + t26 * t98;
t44 = cos(qJ(6));
t86 = rSges(7,1) * t44;
t92 = (-pkin(5) - t86) * t26;
t91 = pkin(4) * t29;
t42 = sin(qJ(2));
t88 = t42 * pkin(2);
t87 = rSges(3,3) + pkin(7);
t80 = t43 * t41;
t79 = t43 * t44;
t78 = t46 * t41;
t77 = t46 * t44;
t40 = -qJ(3) - pkin(7);
t76 = rSges(4,3) - t40;
t38 = -pkin(8) + t40;
t75 = rSges(5,3) - t38;
t36 = -pkin(9) + t38;
t74 = rSges(6,3) - t36;
t73 = pkin(3) * t34 + t37;
t12 = pkin(1) + t73;
t71 = t97 * t43;
t70 = t97 * t46;
t25 = pkin(4) * t30;
t66 = t25 + t95;
t18 = -pkin(3) * t33 - t88;
t65 = t45 * rSges(3,1) - t42 * rSges(3,2);
t62 = -rSges(5,1) * t29 - rSges(5,2) * t30;
t60 = -rSges(6,1) * t26 - rSges(6,2) * t27;
t59 = pkin(1) + t65;
t58 = pkin(1) + t94;
t57 = t12 + t68;
t55 = t60 * t43;
t54 = t60 * t46;
t52 = t53 + (-t85 + t86) * t27;
t50 = t25 + t52;
t49 = g(1) * t70 + g(2) * t71;
t48 = t93 * t92;
t9 = t18 - t91;
t8 = t27 * t77 + t80;
t7 = -t27 * t78 + t79;
t6 = -t27 * t79 + t78;
t5 = t27 * t80 + t77;
t4 = t25 + t12;
t3 = t46 * t9;
t2 = t43 * t9;
t1 = t46 * t4;
t10 = [-m(2) * (g(1) * (-t43 * rSges(2,1) - t46 * rSges(2,2)) + g(2) * (t46 * rSges(2,1) - t43 * rSges(2,2))) - m(3) * ((g(1) * t87 + g(2) * t59) * t46 + (-g(1) * t59 + g(2) * t87) * t43) - m(4) * ((g(1) * t76 + g(2) * t58) * t46 + (-g(1) * t58 + g(2) * t76) * t43) - m(5) * ((g(1) * t75 + g(2) * t57) * t46 + (-g(1) * t57 + g(2) * t75) * t43) - m(6) * (g(2) * t1 + (g(1) * t74 + g(2) * t95) * t46 + (g(1) * (-t4 - t95) + g(2) * t74) * t43) - m(7) * (g(1) * (t6 * rSges(7,1) + t5 * rSges(7,2)) + g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) + t1) + (-g(1) * t36 + g(2) * t53) * t46 + (g(1) * (-t4 - t53) - g(2) * t36) * t43) -m(3) * (g(3) * t65 + t93 * (-rSges(3,1) * t42 - rSges(3,2) * t45)) - m(4) * (g(3) * t94 + t93 * (-rSges(4,1) * t33 - rSges(4,2) * t34 - t88)) - m(5) * (g(3) * (t68 + t73) + t93 * (t18 + t62)) - m(6) * (g(1) * (t3 + t54) + g(2) * (t2 + t55) + g(3) * (t66 + t73)) - m(7) * (g(1) * (t3 + t70) + g(2) * (t2 + t71) + g(3) * (t50 + t73) + t48) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t43 - g(2) * t46) -m(7) * t49 + (-m(5) * t68 - m(6) * t66 - m(7) * t50) * g(3) + t93 * (-m(5) * t62 - m(6) * (t60 - t91) - m(7) * (-t91 + t92)) -m(6) * (g(1) * t54 + g(2) * t55 + g(3) * t95) - m(7) * (g(3) * t52 + t48 + t49) -m(7) * (g(1) * (t7 * rSges(7,1) - t8 * rSges(7,2)) + g(2) * (-t5 * rSges(7,1) + t6 * rSges(7,2)) + g(3) * (-rSges(7,1) * t41 - rSges(7,2) * t44) * t26)];
taug  = t10(:);
