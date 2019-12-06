% Calculate Gravitation load on the joints for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:16
% EndTime: 2019-12-05 17:23:18
% DurationCPUTime: 0.86s
% Computational Cost: add. (664->165), mult. (1827->270), div. (0->0), fcn. (2307->14), ass. (0->77)
t50 = sin(pkin(6));
t55 = sin(qJ(2));
t79 = cos(pkin(11));
t81 = cos(pkin(5));
t64 = t81 * t79;
t78 = sin(pkin(11));
t89 = cos(qJ(2));
t60 = t78 * t55 - t89 * t64;
t51 = sin(pkin(5));
t74 = t51 * t79;
t80 = cos(pkin(6));
t96 = t50 * t74 + t60 * t80;
t63 = t81 * t78;
t61 = t79 * t55 + t89 * t63;
t73 = t51 * t78;
t95 = -t50 * t73 + t61 * t80;
t94 = pkin(8) * t50;
t57 = cos(qJ(4));
t93 = t57 * pkin(4);
t91 = rSges(5,3) + pkin(9);
t90 = rSges(6,3) + pkin(10);
t88 = cos(qJ(3));
t53 = sin(qJ(4));
t87 = t50 * t53;
t86 = t50 * t57;
t85 = t51 * t55;
t52 = sin(qJ(5));
t84 = t52 * t57;
t56 = cos(qJ(5));
t83 = t56 * t57;
t75 = t51 * t89;
t77 = t50 * t85;
t82 = pkin(2) * t75 + pkin(8) * t77;
t54 = sin(qJ(3));
t72 = t54 * t80;
t36 = (-t55 * t72 + t88 * t89) * t51;
t76 = t36 * pkin(3) + t82;
t71 = t81 * t50;
t68 = t80 * t88;
t67 = -rSges(5,1) * t57 + rSges(5,2) * t53;
t40 = t55 * t64 + t78 * t89;
t20 = -t40 * t72 - t60 * t88;
t37 = t60 * pkin(2);
t66 = t20 * pkin(3) + t40 * t94 - t37;
t41 = -t55 * t63 + t79 * t89;
t22 = -t41 * t72 - t61 * t88;
t38 = t61 * pkin(2);
t65 = t22 * pkin(3) + t41 * t94 - t38;
t62 = rSges(6,1) * t56 - rSges(6,2) * t52 + pkin(4);
t39 = -t50 * t75 + t81 * t80;
t35 = (t89 * t54 + t55 * t68) * t51;
t29 = t61 * t50 + t80 * t73;
t28 = t60 * t50 - t80 * t74;
t27 = t54 * t71 + (t88 * t55 + t89 * t72) * t51;
t26 = t54 * t85 - t68 * t75 - t88 * t71;
t25 = t26 * pkin(3);
t24 = t36 * t57 + t53 * t77;
t23 = t36 * t53 - t57 * t77;
t21 = t41 * t68 - t61 * t54;
t19 = t40 * t68 - t60 * t54;
t16 = t27 * t57 + t39 * t53;
t15 = -t27 * t53 + t39 * t57;
t14 = t41 * t88 - t95 * t54;
t13 = t41 * t54 + t95 * t88;
t12 = t40 * t88 - t96 * t54;
t11 = t40 * t54 + t96 * t88;
t10 = t13 * pkin(3);
t9 = t11 * pkin(3);
t8 = t22 * t57 + t41 * t87;
t7 = t22 * t53 - t41 * t86;
t6 = t20 * t57 + t40 * t87;
t5 = t20 * t53 - t40 * t86;
t4 = t14 * t57 + t29 * t53;
t3 = -t14 * t53 + t29 * t57;
t2 = t12 * t57 + t28 * t53;
t1 = -t12 * t53 + t28 * t57;
t17 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(1) * (-t61 * rSges(3,1) - t41 * rSges(3,2)) + g(2) * (-t60 * rSges(3,1) - t40 * rSges(3,2)) + g(3) * (t89 * rSges(3,1) - rSges(3,2) * t55) * t51) - m(4) * (g(1) * (t22 * rSges(4,1) - t21 * rSges(4,2) - t38) + g(2) * (t20 * rSges(4,1) - t19 * rSges(4,2) - t37) + g(3) * (t36 * rSges(4,1) - t35 * rSges(4,2) + t82) + (g(3) * rSges(4,3) * t85 + (g(1) * t41 + g(2) * t40) * (rSges(4,3) + pkin(8))) * t50) - m(5) * (g(1) * (t8 * rSges(5,1) - t7 * rSges(5,2) + t91 * t21 + t65) + g(2) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t91 * t19 + t66) + g(3) * (t24 * rSges(5,1) - t23 * rSges(5,2) + t91 * t35 + t76)) - m(6) * (g(1) * (t8 * pkin(4) + t21 * pkin(9) + (t21 * t52 + t8 * t56) * rSges(6,1) + (t21 * t56 - t8 * t52) * rSges(6,2) + t90 * t7 + t65) + g(2) * (t6 * pkin(4) + t19 * pkin(9) + (t19 * t52 + t6 * t56) * rSges(6,1) + (t19 * t56 - t6 * t52) * rSges(6,2) + t90 * t5 + t66) + g(3) * (t24 * pkin(4) + t35 * pkin(9) + (t24 * t56 + t35 * t52) * rSges(6,1) + (-t24 * t52 + t35 * t56) * rSges(6,2) + t90 * t23 + t76)), -m(4) * (g(1) * (-t13 * rSges(4,1) - t14 * rSges(4,2)) + g(2) * (-t11 * rSges(4,1) - t12 * rSges(4,2)) + g(3) * (-t26 * rSges(4,1) - t27 * rSges(4,2))) - m(5) * (g(1) * (t67 * t13 + t91 * t14 - t10) + g(2) * (t67 * t11 + t91 * t12 - t9) + g(3) * (t67 * t26 + t91 * t27 - t25)) + (-g(1) * (-t13 * t93 - t10 + t14 * pkin(9) + (-t13 * t83 + t14 * t52) * rSges(6,1) + (t13 * t84 + t14 * t56) * rSges(6,2)) - g(2) * (-t11 * t93 - t9 + t12 * pkin(9) + (-t11 * t83 + t12 * t52) * rSges(6,1) + (t11 * t84 + t12 * t56) * rSges(6,2)) - g(3) * (-t26 * t93 - t25 + t27 * pkin(9) + (-t26 * t83 + t27 * t52) * rSges(6,1) + (t26 * t84 + t27 * t56) * rSges(6,2)) - (-g(1) * t13 - g(2) * t11 - g(3) * t26) * t53 * t90) * m(6), -m(5) * (g(1) * (t3 * rSges(5,1) - t4 * rSges(5,2)) + g(2) * (t1 * rSges(5,1) - t2 * rSges(5,2)) + g(3) * (t15 * rSges(5,1) - t16 * rSges(5,2))) - m(6) * (g(1) * (t62 * t3 + t90 * t4) + (t62 * t15 + t90 * t16) * g(3) + (t62 * t1 + t90 * t2) * g(2)), -m(6) * (g(1) * ((t13 * t56 - t4 * t52) * rSges(6,1) + (-t13 * t52 - t4 * t56) * rSges(6,2)) + g(2) * ((t11 * t56 - t2 * t52) * rSges(6,1) + (-t11 * t52 - t2 * t56) * rSges(6,2)) + g(3) * ((-t16 * t52 + t26 * t56) * rSges(6,1) + (-t16 * t56 - t26 * t52) * rSges(6,2)))];
taug = t17(:);
