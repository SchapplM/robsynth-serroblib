% Calculate Gravitation load on the joints for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:02:11
% EndTime: 2019-03-09 18:02:13
% DurationCPUTime: 0.65s
% Computational Cost: add. (597->124), mult. (434->152), div. (0->0), fcn. (373->12), ass. (0->72)
t103 = rSges(7,3) + pkin(10);
t43 = qJ(2) + qJ(3);
t37 = pkin(11) + t43;
t35 = qJ(5) + t37;
t29 = sin(t35);
t30 = cos(t35);
t44 = sin(qJ(6));
t90 = rSges(7,2) * t44;
t102 = t30 * t103 + t29 * t90;
t100 = t30 * rSges(6,1) - t29 * rSges(6,2);
t32 = sin(t37);
t33 = cos(t37);
t99 = t33 * rSges(5,1) - t32 * rSges(5,2);
t38 = sin(t43);
t39 = cos(t43);
t72 = t39 * rSges(4,1) - t38 * rSges(4,2);
t96 = pkin(3) * t38;
t98 = -rSges(5,1) * t32 - rSges(5,2) * t33 - t96;
t46 = sin(qJ(1));
t49 = cos(qJ(1));
t97 = g(1) * t49 + g(2) * t46;
t56 = t30 * pkin(5) + t103 * t29;
t50 = -pkin(8) - pkin(7);
t45 = sin(qJ(2));
t93 = t45 * pkin(2);
t92 = rSges(3,3) + pkin(7);
t48 = cos(qJ(2));
t41 = t48 * pkin(2);
t36 = t41 + pkin(1);
t47 = cos(qJ(6));
t91 = rSges(7,1) * t47;
t84 = t46 * t44;
t83 = t46 * t47;
t82 = t49 * t44;
t81 = t49 * t47;
t80 = rSges(4,3) - t50;
t42 = -qJ(4) + t50;
t79 = rSges(5,3) - t42;
t40 = -pkin(9) + t42;
t78 = rSges(6,3) - t40;
t28 = pkin(4) * t33;
t34 = pkin(3) * t39;
t77 = t28 + t34;
t15 = t34 + t36;
t75 = t102 * t46;
t74 = t102 * t49;
t70 = t34 + t99;
t12 = -pkin(4) * t32 - t96;
t69 = t48 * rSges(3,1) - t45 * rSges(3,2);
t67 = -rSges(4,1) * t38 - rSges(4,2) * t39;
t64 = -rSges(6,1) * t29 - rSges(6,2) * t30;
t63 = t100 + t77;
t62 = pkin(1) + t69;
t61 = t36 + t72;
t60 = t15 + t99;
t58 = t64 * t46;
t57 = t64 * t49;
t55 = t56 + (-t90 + t91) * t30;
t52 = t55 + t77;
t51 = t97 * (-pkin(5) - t91) * t29;
t11 = t49 * t12;
t10 = t46 * t12;
t9 = t12 - t93;
t8 = t28 + t15;
t7 = t30 * t81 + t84;
t6 = -t30 * t82 + t83;
t5 = -t30 * t83 + t82;
t4 = t30 * t84 + t81;
t3 = t49 * t9;
t2 = t46 * t9;
t1 = t49 * t8;
t13 = [-m(2) * (g(1) * (-t46 * rSges(2,1) - t49 * rSges(2,2)) + g(2) * (t49 * rSges(2,1) - t46 * rSges(2,2))) - m(3) * ((g(1) * t92 + g(2) * t62) * t49 + (-g(1) * t62 + g(2) * t92) * t46) - m(4) * ((g(1) * t80 + g(2) * t61) * t49 + (-g(1) * t61 + g(2) * t80) * t46) - m(5) * ((g(1) * t79 + g(2) * t60) * t49 + (-g(1) * t60 + g(2) * t79) * t46) - m(6) * (g(2) * t1 + (g(1) * t78 + g(2) * t100) * t49 + (g(1) * (-t100 - t8) + g(2) * t78) * t46) - m(7) * (g(1) * (t5 * rSges(7,1) + t4 * rSges(7,2)) + g(2) * (t7 * rSges(7,1) + t6 * rSges(7,2) + t1) + (-g(1) * t40 + g(2) * t56) * t49 + (g(1) * (-t56 - t8) - g(2) * t40) * t46) -m(3) * (g(3) * t69 + t97 * (-rSges(3,1) * t45 - rSges(3,2) * t48)) - m(4) * (g(3) * (t41 + t72) + t97 * (t67 - t93)) - m(5) * (g(3) * (t41 + t70) + t97 * (-t93 + t98)) - m(6) * (g(1) * (t3 + t57) + g(2) * (t2 + t58) + g(3) * (t41 + t63)) - m(7) * (g(1) * (t3 + t74) + g(2) * (t2 + t75) + g(3) * (t41 + t52) + t51) -m(4) * (g(3) * t72 + t97 * t67) - m(5) * (g(3) * t70 + t97 * t98) - m(6) * (g(1) * (t11 + t57) + g(2) * (t10 + t58) + g(3) * t63) - m(7) * (g(1) * (t11 + t74) + g(2) * (t10 + t75) + g(3) * t52 + t51) (-m(5) - m(6) - m(7)) * (g(1) * t46 - g(2) * t49) -m(6) * (g(1) * t57 + g(2) * t58 + g(3) * t100) - m(7) * (g(1) * t74 + g(2) * t75 + g(3) * t55 + t51) -m(7) * (g(1) * (t6 * rSges(7,1) - t7 * rSges(7,2)) + g(2) * (-t4 * rSges(7,1) + t5 * rSges(7,2)) + g(3) * (-rSges(7,1) * t44 - rSges(7,2) * t47) * t29)];
taug  = t13(:);
