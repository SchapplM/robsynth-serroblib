% Calculate Gravitation load on the joints for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:26
% EndTime: 2019-03-08 22:30:28
% DurationCPUTime: 1.02s
% Computational Cost: add. (504->146), mult. (1133->215), div. (0->0), fcn. (1339->12), ass. (0->66)
t41 = sin(qJ(5));
t44 = cos(qJ(5));
t101 = -rSges(6,1) * t41 - rSges(6,2) * t44;
t45 = cos(qJ(3));
t42 = sin(qJ(3));
t72 = qJ(4) * t42;
t98 = -pkin(3) * t45 - t72;
t97 = rSges(6,3) + pkin(9);
t74 = rSges(7,3) + pkin(10) + pkin(9);
t43 = sin(qJ(2));
t46 = cos(qJ(2));
t69 = cos(pkin(11));
t70 = cos(pkin(6));
t57 = t70 * t69;
t68 = sin(pkin(11));
t21 = t68 * t43 - t46 * t57;
t56 = t70 * t68;
t23 = t69 * t43 + t46 * t56;
t96 = g(1) * t23 + g(2) * t21;
t22 = t43 * t57 + t68 * t46;
t40 = sin(pkin(6));
t64 = t40 * t69;
t11 = t22 * t45 - t42 * t64;
t24 = -t43 * t56 + t69 * t46;
t63 = t40 * t68;
t13 = t24 * t45 + t42 * t63;
t81 = t40 * t43;
t26 = t70 * t42 + t45 * t81;
t95 = g(1) * t13 + g(2) * t11 + g(3) * t26;
t94 = t101 * t42;
t39 = qJ(5) + qJ(6);
t37 = sin(t39);
t38 = cos(t39);
t51 = rSges(7,1) * t37 + rSges(7,2) * t38 + pkin(5) * t41;
t93 = t51 * t42 + t74 * t45;
t92 = t44 * rSges(6,1) - t41 * rSges(6,2) + pkin(4) + pkin(8);
t85 = g(3) * t40;
t84 = rSges(5,1) + pkin(8);
t83 = rSges(4,3) + pkin(8);
t80 = t40 * t46;
t78 = t41 * t46;
t76 = t44 * t46;
t75 = t45 * t46;
t73 = pkin(2) * t80 + pkin(8) * t81;
t71 = rSges(5,3) + qJ(4);
t67 = -m(5) - m(6) - m(7);
t18 = t21 * pkin(2);
t66 = t98 * t21 - t18;
t19 = t23 * pkin(2);
t65 = t98 * t23 - t19;
t62 = rSges(4,1) * t45 - rSges(4,2) * t42;
t61 = rSges(5,2) * t45 - rSges(5,3) * t42;
t10 = t22 * t42 + t45 * t64;
t60 = t10 * t44 - t21 * t41;
t12 = t24 * t42 - t45 * t63;
t59 = t12 * t44 - t23 * t41;
t58 = g(3) * (t40 * pkin(3) * t75 + t72 * t80 + t73);
t25 = t42 * t81 - t70 * t45;
t55 = t25 * t44 + t40 * t78;
t54 = t38 * rSges(7,1) - t37 * rSges(7,2) + pkin(5) * t44 + pkin(4);
t52 = pkin(8) + t54;
t49 = m(7) * (g(1) * ((t12 * t38 - t23 * t37) * rSges(7,1) + (-t12 * t37 - t23 * t38) * rSges(7,2)) + g(2) * ((t10 * t38 - t21 * t37) * rSges(7,1) + (-t10 * t37 - t21 * t38) * rSges(7,2)) + g(3) * ((t25 * t38 + t37 * t80) * rSges(7,1) + (-t25 * t37 + t38 * t80) * rSges(7,2)));
t20 = t25 * pkin(3);
t9 = t12 * pkin(3);
t8 = t10 * pkin(3);
t1 = [(-m(2) - m(3) - m(4) + t67) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t23 - rSges(3,2) * t24) + g(2) * (-rSges(3,1) * t21 - rSges(3,2) * t22) + (rSges(3,1) * t46 - rSges(3,2) * t43) * t85) - m(4) * (g(1) * (-t62 * t23 + t83 * t24 - t19) + g(2) * (-t62 * t21 + t83 * t22 - t18) + g(3) * t73 + (rSges(4,3) * t43 + t46 * t62) * t85) - m(5) * (g(1) * (t61 * t23 + t84 * t24 + t65) + g(2) * (t61 * t21 + t84 * t22 + t66) + t58 + (rSges(5,1) * t43 - t46 * t61) * t85) - m(6) * (g(1) * (t23 * t94 + t24 * t92 + t65) + g(2) * (t21 * t94 + t22 * t92 + t66) + t58 + (t43 * pkin(4) + (t42 * t78 + t43 * t44) * rSges(6,1) + (-t41 * t43 + t42 * t76) * rSges(6,2)) * t85 + (-t45 * t96 + t75 * t85) * t97) - m(7) * (t58 + (t54 * t43 + t46 * t93) * t85 - t96 * t93 + (t52 * t22 + t66) * g(2) + (t52 * t24 + t65) * g(1)) -m(4) * (g(1) * (-rSges(4,1) * t12 - rSges(4,2) * t13) + g(2) * (-rSges(4,1) * t10 - rSges(4,2) * t11) + g(3) * (-rSges(4,1) * t25 - rSges(4,2) * t26)) - m(5) * (g(1) * (rSges(5,2) * t12 + t71 * t13 - t9) + g(2) * (rSges(5,2) * t10 + t71 * t11 - t8) + g(3) * (rSges(5,2) * t25 + t71 * t26 - t20)) + (-g(1) * (-t12 * t74 - t9) - g(2) * (-t10 * t74 - t8) - g(3) * (-t25 * t74 - t20) - t95 * (qJ(4) + t51)) * m(7) + (-g(1) * (-t12 * t97 - t9) - g(2) * (-t10 * t97 - t8) - g(3) * (-t25 * t97 - t20) - t95 * (qJ(4) - t101)) * m(6), t67 * (g(1) * t12 + g(2) * t10 + g(3) * t25) -m(6) * (g(1) * (t59 * rSges(6,1) + (-t12 * t41 - t23 * t44) * rSges(6,2)) + g(2) * (t60 * rSges(6,1) + (-t10 * t41 - t21 * t44) * rSges(6,2)) + g(3) * (t55 * rSges(6,1) + (-t25 * t41 + t40 * t76) * rSges(6,2))) - t49 - m(7) * (g(1) * t59 + g(2) * t60 + g(3) * t55) * pkin(5), -t49];
taug  = t1(:);
