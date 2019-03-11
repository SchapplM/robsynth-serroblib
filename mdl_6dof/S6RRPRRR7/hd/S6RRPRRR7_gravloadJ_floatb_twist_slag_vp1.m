% Calculate Gravitation load on the joints for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:55:45
% EndTime: 2019-03-09 13:55:48
% DurationCPUTime: 1.25s
% Computational Cost: add. (432->156), mult. (869->204), div. (0->0), fcn. (940->10), ass. (0->65)
t39 = cos(qJ(5));
t26 = pkin(5) * t39 + pkin(4);
t35 = qJ(5) + qJ(6);
t27 = sin(t35);
t28 = cos(t35);
t37 = sin(qJ(2));
t40 = cos(qJ(2));
t74 = sin(qJ(4));
t75 = cos(qJ(4));
t15 = t37 * t75 - t40 * t74;
t38 = sin(qJ(1));
t10 = t15 * t38;
t41 = cos(qJ(1));
t61 = t41 * t74;
t62 = t41 * t75;
t13 = -t37 * t62 + t40 * t61;
t14 = t37 * t74 + t40 * t75;
t91 = -g(1) * t13 + g(2) * t10 - g(3) * t14;
t99 = (rSges(7,1) * t28 - rSges(7,2) * t27 + t26) * t91;
t36 = sin(qJ(5));
t98 = (rSges(6,1) * t39 - rSges(6,2) * t36 + pkin(4)) * t91;
t78 = g(3) * t15;
t97 = -pkin(9) - rSges(6,3);
t29 = t37 * qJ(3);
t69 = t40 * pkin(2) + t29;
t96 = -pkin(10) - pkin(9) - rSges(7,3);
t95 = g(1) * t41 + g(2) * t38;
t92 = t95 * t37;
t11 = t14 * t38;
t12 = -t37 * t61 - t40 * t62;
t90 = g(1) * t12 - g(2) * t11 - t78;
t50 = t11 * t27 - t28 * t41;
t51 = -t11 * t28 - t27 * t41;
t88 = -t50 * rSges(7,1) + t51 * rSges(7,2);
t5 = t12 * t27 - t28 * t38;
t6 = -t12 * t28 - t27 * t38;
t87 = t5 * rSges(7,1) - t6 * rSges(7,2);
t86 = pkin(5) * t36;
t31 = t40 * pkin(3);
t73 = rSges(4,1) * t40;
t72 = t36 * t41;
t71 = t37 * t41;
t70 = t40 * t41;
t68 = t41 * pkin(1) + t38 * pkin(7);
t66 = qJ(3) * t40;
t65 = t31 + t69;
t33 = t41 * pkin(7);
t60 = -t41 * pkin(8) + t33;
t59 = pkin(2) * t70 + t41 * t29 + t68;
t58 = pkin(3) * t70 + t59;
t57 = rSges(3,1) * t40 - rSges(3,2) * t37;
t55 = rSges(5,1) * t10 - rSges(5,2) * t11;
t54 = -rSges(5,1) * t13 + rSges(5,2) * t12;
t53 = -rSges(5,1) * t14 - rSges(5,2) * t15;
t52 = -t27 * rSges(7,1) - t28 * rSges(7,2);
t49 = -t11 * t39 - t72;
t48 = t11 * t36 - t39 * t41;
t7 = t12 * t36 - t38 * t39;
t46 = -pkin(1) - t69;
t44 = g(1) * (t46 - t31);
t43 = (-pkin(2) - pkin(3)) * t92;
t22 = t41 * t66;
t20 = t38 * t66;
t8 = -t12 * t39 - t36 * t38;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t38 - rSges(2,2) * t41) + g(2) * (rSges(2,1) * t41 - rSges(2,2) * t38)) - m(3) * (g(1) * (rSges(3,3) * t41 + t33) + g(2) * (rSges(3,1) * t70 - rSges(3,2) * t71 + t68) + (g(1) * (-pkin(1) - t57) + g(2) * rSges(3,3)) * t38) - m(4) * (g(1) * (rSges(4,2) * t41 + t33) + g(2) * (rSges(4,1) * t70 + rSges(4,3) * t71 + t59) + (g(1) * (-rSges(4,3) * t37 + t46 - t73) + g(2) * rSges(4,2)) * t38) - m(5) * (g(1) * (-rSges(5,1) * t11 - rSges(5,2) * t10 - rSges(5,3) * t41 + t60) + g(2) * (-rSges(5,1) * t12 - rSges(5,2) * t13 + t58) + (t44 + g(2) * (-rSges(5,3) - pkin(8))) * t38) - m(6) * (g(1) * (t49 * rSges(6,1) + t48 * rSges(6,2) - t11 * pkin(4) - t10 * t97 + t60) + g(2) * (rSges(6,1) * t8 + rSges(6,2) * t7 - pkin(4) * t12 - t13 * t97 + t58) + (-g(2) * pkin(8) + t44) * t38) - m(7) * (g(1) * (t51 * rSges(7,1) + t50 * rSges(7,2) - pkin(5) * t72 - t10 * t96 - t11 * t26 + t60) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) - t12 * t26 - t96 * t13 + t58) + (t44 + g(2) * (-pkin(8) - t86)) * t38) -m(3) * (g(3) * t57 + t95 * (-rSges(3,1) * t37 - rSges(3,2) * t40)) - m(4) * (g(1) * (rSges(4,3) * t70 + t22) + g(2) * (rSges(4,3) * t38 * t40 + t20) + g(3) * (t69 + t73) + (g(3) * rSges(4,3) + t95 * (-rSges(4,1) - pkin(2))) * t37) - m(5) * (g(1) * (t22 - t54) + g(2) * (t20 - t55) + g(3) * (-t53 + t65) + t43) + (-g(1) * (-t12 * t96 + t22) - g(2) * (t11 * t96 + t20) - g(3) * (t15 * t96 + t65) - t43 + t99) * m(7) + (-g(1) * (-t12 * t97 + t22) - g(2) * (t11 * t97 + t20) - g(3) * (t15 * t97 + t65) - t43 + t98) * m(6) (-m(4) - m(5) - m(6) - m(7)) * (-g(3) * t40 + t92) -m(5) * (g(1) * t54 + g(2) * t55 + g(3) * t53) - m(6) * (t90 * t97 + t98) - m(7) * (t90 * t96 + t99) -m(6) * (g(1) * (rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (-t48 * rSges(6,1) + t49 * rSges(6,2))) - m(7) * (g(1) * (t7 * pkin(5) + t87) + g(2) * (-t48 * pkin(5) + t88)) + (-m(6) * (-t36 * rSges(6,1) - t39 * rSges(6,2)) - m(7) * (t52 - t86)) * t78, -m(7) * (g(1) * t87 + g(2) * t88 + t52 * t78)];
taug  = t1(:);
