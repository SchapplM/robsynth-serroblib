% Calculate Gravitation load on the joints for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:16
% EndTime: 2019-03-08 20:09:18
% DurationCPUTime: 0.59s
% Computational Cost: add. (580->132), mult. (993->202), div. (0->0), fcn. (1176->12), ass. (0->72)
t93 = rSges(7,2) + pkin(9);
t50 = pkin(11) + qJ(4);
t49 = cos(t50);
t92 = pkin(4) * t49;
t52 = sin(pkin(10));
t57 = sin(qJ(2));
t59 = cos(qJ(2));
t72 = cos(pkin(10));
t73 = cos(pkin(6));
t62 = t73 * t72;
t35 = t52 * t57 - t59 * t62;
t91 = g(2) * t35;
t53 = sin(pkin(6));
t90 = g(3) * t53;
t89 = rSges(7,1) + pkin(5);
t88 = rSges(6,3) + pkin(9);
t48 = sin(t50);
t87 = t35 * t48;
t68 = t52 * t73;
t37 = t72 * t57 + t59 * t68;
t86 = t37 * t48;
t85 = t48 * t59;
t56 = sin(qJ(5));
t84 = t49 * t56;
t58 = cos(qJ(5));
t83 = t49 * t58;
t82 = t52 * t53;
t81 = t53 * t57;
t80 = t53 * t59;
t55 = -pkin(8) - qJ(3);
t79 = t55 * t57;
t78 = t58 * t59;
t36 = t52 * t59 + t57 * t62;
t54 = cos(pkin(11));
t47 = t54 * pkin(3) + pkin(2);
t77 = -t35 * t47 - t36 * t55;
t38 = -t57 * t68 + t72 * t59;
t76 = -t37 * t47 - t38 * t55;
t75 = rSges(4,3) + qJ(3);
t74 = rSges(7,3) + qJ(6);
t71 = t56 * t80;
t39 = t47 * t80;
t70 = t39 + (pkin(9) * t48 + t92) * t80;
t69 = -m(4) - m(5) - m(6) - m(7);
t67 = t53 * t72;
t66 = -pkin(9) * t87 - t35 * t92 + t77;
t65 = -pkin(9) * t86 - t37 * t92 + t76;
t64 = rSges(5,1) * t49 - rSges(5,2) * t48;
t63 = rSges(6,1) * t58 - rSges(6,2) * t56;
t61 = rSges(4,1) * t54 - rSges(4,2) * sin(pkin(11)) + pkin(2);
t28 = t73 * t48 + t49 * t81;
t27 = -t48 * t81 + t73 * t49;
t26 = t27 * pkin(4);
t19 = (t49 * t78 + t56 * t57) * t53;
t18 = t49 * t71 - t58 * t81;
t17 = t28 * t58 - t71;
t16 = t28 * t56 + t53 * t78;
t15 = t38 * t49 + t48 * t82;
t14 = -t38 * t48 + t49 * t82;
t13 = t36 * t49 - t48 * t67;
t12 = -t36 * t48 - t49 * t67;
t11 = t14 * pkin(4);
t10 = t12 * pkin(4);
t8 = -t37 * t83 + t38 * t56;
t7 = -t37 * t84 - t38 * t58;
t6 = -t35 * t83 + t36 * t56;
t5 = -t35 * t84 - t36 * t58;
t4 = t15 * t58 + t37 * t56;
t3 = t15 * t56 - t37 * t58;
t2 = t13 * t58 + t35 * t56;
t1 = t13 * t56 - t35 * t58;
t9 = [(-m(2) - m(3) + t69) * g(3), -m(3) * (g(1) * (-t37 * rSges(3,1) - t38 * rSges(3,2)) + g(2) * (-t35 * rSges(3,1) - t36 * rSges(3,2)) + (rSges(3,1) * t59 - rSges(3,2) * t57) * t90) - m(4) * (g(1) * (-t61 * t37 + t75 * t38) + g(2) * t75 * t36 - t61 * t91 + (t75 * t57 + t61 * t59) * t90) - m(5) * (g(1) * (t38 * rSges(5,3) - t64 * t37 + t76) + g(2) * (t36 * rSges(5,3) - t64 * t35 + t77) + g(3) * t39 + (t64 * t59 + (rSges(5,3) - t55) * t57) * t90) - m(6) * (g(1) * (t8 * rSges(6,1) - t7 * rSges(6,2) - rSges(6,3) * t86 + t65) + g(2) * (t6 * rSges(6,1) - t5 * rSges(6,2) - rSges(6,3) * t87 + t66) + g(3) * (t19 * rSges(6,1) - t18 * rSges(6,2) + (rSges(6,3) * t85 - t79) * t53 + t70)) - m(7) * (g(1) * (-rSges(7,2) * t86 + t74 * t7 + t89 * t8 + t65) + g(2) * (-rSges(7,2) * t87 + t74 * t5 + t89 * t6 + t66) + g(3) * ((rSges(7,2) * t85 - t79) * t53 + t89 * t19 + t74 * t18 + t70)) t69 * (g(1) * t37 - g(3) * t80 + t91) -m(5) * (g(1) * (t14 * rSges(5,1) - t15 * rSges(5,2)) + g(2) * (t12 * rSges(5,1) - t13 * rSges(5,2)) + g(3) * (t27 * rSges(5,1) - t28 * rSges(5,2))) - m(6) * (g(1) * (t63 * t14 + t88 * t15 + t11) + g(2) * (t63 * t12 + t88 * t13 + t10) + g(3) * (t63 * t27 + t88 * t28 + t26)) + (-g(1) * (t93 * t15 + t11) - g(2) * (t93 * t13 + t10) - g(3) * (t93 * t28 + t26) - (g(1) * t14 + g(2) * t12 + g(3) * t27) * (t74 * t56 + t89 * t58)) * m(7), -m(6) * (g(1) * (-t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) - t2 * rSges(6,2)) + g(3) * (-t16 * rSges(6,1) - t17 * rSges(6,2))) - m(7) * (g(1) * (-t89 * t3 + t74 * t4) + g(2) * (-t89 * t1 + t74 * t2) + g(3) * (-t89 * t16 + t74 * t17)) -m(7) * (g(1) * t3 + g(2) * t1 + g(3) * t16)];
taug  = t9(:);
