% Calculate Gravitation load on the joints for
% S6RRPRRR3
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
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:22:00
% EndTime: 2019-03-09 13:22:03
% DurationCPUTime: 0.75s
% Computational Cost: add. (551->150), mult. (529->196), div. (0->0), fcn. (499->12), ass. (0->78)
t38 = qJ(2) + pkin(11);
t31 = sin(t38);
t82 = rSges(5,3) + pkin(8);
t95 = t82 * t31;
t48 = -pkin(9) - pkin(8);
t68 = rSges(7,3) + pkin(10) - t48;
t94 = t68 * t31;
t69 = rSges(6,3) - t48;
t93 = t69 * t31;
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t92 = g(1) * t47 + g(2) * t44;
t40 = qJ(4) + qJ(5);
t35 = qJ(6) + t40;
t27 = sin(t35);
t28 = cos(t35);
t32 = cos(t38);
t80 = t32 * t44;
t5 = t27 * t80 + t28 * t47;
t6 = t27 * t47 - t28 * t80;
t91 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t79 = t32 * t47;
t7 = -t27 * t79 + t28 * t44;
t8 = t27 * t44 + t28 * t79;
t90 = t7 * rSges(7,1) - t8 * rSges(7,2);
t43 = sin(qJ(2));
t89 = pkin(2) * t43;
t42 = sin(qJ(4));
t88 = pkin(4) * t42;
t33 = sin(t40);
t87 = pkin(5) * t33;
t84 = g(3) * t31;
t83 = rSges(3,3) + pkin(7);
t34 = cos(t40);
t75 = t34 * t47;
t78 = t33 * t44;
t13 = t32 * t78 + t75;
t76 = t34 * t44;
t77 = t33 * t47;
t14 = -t32 * t76 + t77;
t81 = -t13 * rSges(6,1) + t14 * rSges(6,2);
t74 = t42 * t44;
t73 = t42 * t47;
t45 = cos(qJ(4));
t72 = t44 * t45;
t71 = t45 * t47;
t41 = -qJ(3) - pkin(7);
t70 = rSges(4,3) - t41;
t15 = -t32 * t77 + t76;
t16 = t32 * t75 + t78;
t67 = t15 * rSges(6,1) - t16 * rSges(6,2);
t22 = t87 + t88;
t66 = t22 - t41;
t36 = t45 * pkin(4);
t23 = pkin(5) * t34 + t36;
t65 = -t41 + t88;
t46 = cos(qJ(2));
t64 = rSges(3,1) * t46 - rSges(3,2) * t43;
t62 = rSges(4,1) * t32 - rSges(4,2) * t31;
t61 = -rSges(6,1) * t33 - rSges(6,2) * t34;
t60 = -rSges(7,1) * t27 - rSges(7,2) * t28;
t59 = pkin(1) + t64;
t58 = rSges(5,1) * t45 - rSges(5,2) * t42 + pkin(3);
t19 = -t32 * t73 + t72;
t17 = t32 * t74 + t71;
t29 = t36 + pkin(3);
t57 = rSges(6,1) * t34 - rSges(6,2) * t33 + t29;
t21 = pkin(3) + t23;
t56 = rSges(7,1) * t28 - rSges(7,2) * t27 + t21;
t55 = t32 * pkin(3) + t95;
t53 = t32 * t29 + t93;
t52 = t32 * t21 + t94;
t37 = t46 * pkin(2);
t30 = t37 + pkin(1);
t25 = t47 * t30;
t20 = t32 * t71 + t74;
t18 = -t32 * t72 + t73;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t44 - rSges(2,2) * t47) + g(2) * (rSges(2,1) * t47 - rSges(2,2) * t44)) - m(3) * ((g(1) * t83 + g(2) * t59) * t47 + (-g(1) * t59 + g(2) * t83) * t44) - m(4) * (g(2) * t25 + (g(1) * t70 + g(2) * t62) * t47 + (g(1) * (-t30 - t62) + g(2) * t70) * t44) - m(5) * (g(1) * (rSges(5,1) * t18 + rSges(5,2) * t17) + g(2) * (rSges(5,1) * t20 + rSges(5,2) * t19 + t25) + (-g(1) * t41 + g(2) * t55) * t47 + (g(1) * (-t30 - t55) - g(2) * t41) * t44) - m(6) * (g(1) * (t14 * rSges(6,1) + t13 * rSges(6,2)) + g(2) * (t16 * rSges(6,1) + t15 * rSges(6,2) + t25) + (g(1) * t65 + g(2) * t53) * t47 + (g(1) * (-t30 - t53) + g(2) * t65) * t44) - m(7) * (g(1) * (rSges(7,1) * t6 + rSges(7,2) * t5) + g(2) * (rSges(7,1) * t8 + rSges(7,2) * t7 + t25) + (g(1) * t66 + g(2) * t52) * t47 + (g(1) * (-t30 - t52) + g(2) * t66) * t44) -m(3) * (g(3) * t64 + t92 * (-rSges(3,1) * t43 - rSges(3,2) * t46)) - m(4) * (g(3) * (t37 + t62) + t92 * (-rSges(4,1) * t31 - rSges(4,2) * t32 - t89)) - m(5) * (g(3) * (t58 * t32 + t37 + t95) + t92 * (-t58 * t31 + t82 * t32 - t89)) - m(6) * (g(3) * (t57 * t32 + t37 + t93) + t92 * (-t57 * t31 + t69 * t32 - t89)) - m(7) * (g(3) * (t56 * t32 + t37 + t94) + t92 * (-t56 * t31 + t68 * t32 - t89)) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t44 - g(2) * t47) -m(5) * (g(1) * (rSges(5,1) * t19 - rSges(5,2) * t20) + g(2) * (-rSges(5,1) * t17 + rSges(5,2) * t18)) - m(6) * (g(1) * (t19 * pkin(4) + t67) + g(2) * (-t17 * pkin(4) + t81)) - m(7) * (g(1) * (-t22 * t79 + t23 * t44 + t90) + g(2) * (-t22 * t80 - t23 * t47 + t91)) + (-m(5) * (-rSges(5,1) * t42 - rSges(5,2) * t45) - m(6) * (t61 - t88) - m(7) * (-t22 + t60)) * t84, -m(6) * (g(1) * t67 + g(2) * t81) - m(7) * (g(1) * (t15 * pkin(5) + t90) + g(2) * (-t13 * pkin(5) + t91)) + (-m(6) * t61 - m(7) * (t60 - t87)) * t84, -m(7) * (g(1) * t90 + g(2) * t91 + t60 * t84)];
taug  = t1(:);
