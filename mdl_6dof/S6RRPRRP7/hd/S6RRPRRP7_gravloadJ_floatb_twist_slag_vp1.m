% Calculate Gravitation load on the joints for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:16:47
% EndTime: 2019-03-09 12:16:50
% DurationCPUTime: 0.96s
% Computational Cost: add. (410->160), mult. (944->201), div. (0->0), fcn. (1040->8), ass. (0->66)
t41 = sin(qJ(2));
t44 = cos(qJ(2));
t81 = sin(qJ(4));
t82 = cos(qJ(4));
t21 = t41 * t82 - t44 * t81;
t42 = sin(qJ(1));
t14 = t21 * t42;
t45 = cos(qJ(1));
t67 = t45 * t81;
t68 = t45 * t82;
t17 = -t41 * t68 + t44 * t67;
t20 = t41 * t81 + t44 * t82;
t40 = sin(qJ(5));
t43 = cos(qJ(5));
t72 = rSges(7,3) + qJ(6);
t84 = rSges(7,1) + pkin(5);
t109 = (t72 * t40 + t84 * t43) * m(7) * (g(1) * t17 - g(2) * t14 + g(3) * t20);
t54 = rSges(6,1) * t43 - rSges(6,2) * t40;
t87 = t17 * pkin(4);
t108 = -t54 * t17 - t87;
t90 = t14 * pkin(4);
t107 = t54 * t14 + t90;
t34 = t41 * qJ(3);
t75 = t44 * pkin(2) + t34;
t15 = t20 * t42;
t89 = t15 * pkin(9);
t104 = t15 * rSges(7,2) + t89 + t90;
t16 = -t41 * t67 - t44 * t68;
t88 = t16 * pkin(9);
t103 = -t16 * rSges(7,2) - t87 - t88;
t102 = g(1) * t45 + g(2) * t42;
t101 = t102 * t41;
t83 = rSges(6,3) + pkin(9);
t99 = t54 * t20 - t83 * t21;
t91 = g(3) * t21;
t86 = t20 * pkin(4);
t36 = t44 * pkin(3);
t78 = t41 * t45;
t77 = t44 * rSges(4,1);
t76 = t44 * t45;
t74 = t45 * pkin(1) + t42 * pkin(7);
t73 = qJ(3) * t44;
t71 = t36 + t75;
t38 = t45 * pkin(7);
t66 = -pkin(8) * t45 + t38;
t1 = t15 * t40 - t45 * t43;
t65 = t86 + t71;
t64 = pkin(2) * t76 + t45 * t34 + t74;
t63 = pkin(3) * t76 + t64;
t60 = (rSges(7,2) + pkin(9)) * t21;
t59 = t44 * rSges(3,1) - t41 * rSges(3,2);
t57 = t14 * rSges(5,1) - t15 * rSges(5,2);
t56 = -t17 * rSges(5,1) + t16 * rSges(5,2);
t55 = -t20 * rSges(5,1) - t21 * rSges(5,2);
t2 = t15 * t43 + t40 * t45;
t53 = -pkin(1) - t75;
t52 = -t15 * pkin(4) + t14 * pkin(9) + t66;
t51 = -t16 * pkin(4) + pkin(9) * t17 + t63;
t48 = g(1) * (t53 - t36);
t47 = (-pkin(2) - pkin(3)) * t101;
t46 = (-g(2) * pkin(8) + t48) * t42;
t28 = t45 * t73;
t26 = t42 * t73;
t6 = -t16 * t43 - t42 * t40;
t5 = -t16 * t40 + t42 * t43;
t3 = [-m(2) * (g(1) * (-t42 * rSges(2,1) - rSges(2,2) * t45) + g(2) * (rSges(2,1) * t45 - t42 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t45 + t38) + g(2) * (rSges(3,1) * t76 - rSges(3,2) * t78 + t74) + (g(1) * (-pkin(1) - t59) + g(2) * rSges(3,3)) * t42) - m(4) * (g(1) * (rSges(4,2) * t45 + t38) + g(2) * (rSges(4,1) * t76 + rSges(4,3) * t78 + t64) + (g(1) * (-t41 * rSges(4,3) + t53 - t77) + g(2) * rSges(4,2)) * t42) - m(5) * (g(1) * (-t15 * rSges(5,1) - t14 * rSges(5,2) - rSges(5,3) * t45 + t66) + g(2) * (-rSges(5,1) * t16 - rSges(5,2) * t17 + t63) + (t48 + g(2) * (-rSges(5,3) - pkin(8))) * t42) - m(6) * (g(1) * (-rSges(6,1) * t2 + rSges(6,2) * t1 + t14 * rSges(6,3) + t52) + g(2) * (rSges(6,1) * t6 - rSges(6,2) * t5 + rSges(6,3) * t17 + t51) + t46) - m(7) * (g(1) * (t14 * rSges(7,2) - t1 * t72 - t2 * t84 + t52) + g(2) * (rSges(7,2) * t17 + t72 * t5 + t6 * t84 + t51) + t46) -m(3) * (g(3) * t59 + t102 * (-rSges(3,1) * t41 - rSges(3,2) * t44)) - m(4) * (g(1) * (rSges(4,3) * t76 + t28) + g(2) * (t42 * t44 * rSges(4,3) + t26) + g(3) * (t75 + t77) + (g(3) * rSges(4,3) + t102 * (-rSges(4,1) - pkin(2))) * t41) - m(5) * (g(1) * (t28 - t56) + g(2) * (t26 - t57) + g(3) * (-t55 + t71) + t47) - m(6) * (g(1) * (t16 * rSges(6,3) - t108 + t28 + t88) + g(2) * (-t15 * rSges(6,3) - t107 + t26 - t89) + g(3) * (t65 + t99) + t47) - m(7) * (g(1) * (t28 - t103) + g(2) * (t26 - t104) + g(3) * (-t60 + t65) + t47) - t109 (-m(4) - m(5) - m(6) - m(7)) * (-g(3) * t44 + t101) -m(5) * (g(1) * t56 + g(2) * t57 + g(3) * t55) - m(6) * (g(1) * (-t83 * t16 + t108) + g(2) * (t83 * t15 + t107) + g(3) * (-t86 - t99)) - m(7) * (g(1) * t103 + g(2) * t104 + g(3) * (-t86 + t60)) + t109, -m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2)) - m(7) * (g(1) * (-t84 * t5 + t72 * t6) + g(2) * (-t84 * t1 + t72 * t2)) + (-m(6) * (-rSges(6,1) * t40 - rSges(6,2) * t43) - m(7) * (-t84 * t40 + t72 * t43)) * t91, -m(7) * (g(1) * t5 + g(2) * t1 + t40 * t91)];
taug  = t3(:);
