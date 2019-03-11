% Calculate Gravitation load on the joints for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:54:28
% EndTime: 2019-03-08 22:54:30
% DurationCPUTime: 0.83s
% Computational Cost: add. (592->167), mult. (1490->241), div. (0->0), fcn. (1811->10), ass. (0->79)
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t108 = pkin(4) * t63 + qJ(5) * t60;
t88 = rSges(7,2) + qJ(5);
t86 = rSges(7,3) + qJ(6);
t107 = -m(6) - m(7);
t64 = cos(qJ(3));
t106 = pkin(3) * t64;
t59 = sin(pkin(6));
t104 = g(3) * t59;
t103 = rSges(6,1) + pkin(9);
t102 = rSges(7,1) + pkin(5);
t101 = rSges(4,3) + pkin(8);
t100 = rSges(5,3) + pkin(9);
t99 = rSges(7,2) * t60;
t58 = sin(pkin(10));
t62 = sin(qJ(2));
t65 = cos(qJ(2));
t84 = cos(pkin(10));
t85 = cos(pkin(6));
t68 = t85 * t84;
t43 = t58 * t62 - t65 * t68;
t61 = sin(qJ(3));
t98 = t43 * t61;
t77 = t58 * t85;
t45 = t84 * t62 + t65 * t77;
t97 = t45 * t61;
t96 = t59 * t62;
t95 = t59 * t64;
t94 = t59 * t65;
t93 = t60 * t64;
t92 = t63 * t64;
t91 = t64 * t65;
t90 = pkin(2) * t94 + pkin(8) * t96;
t87 = rSges(6,3) + qJ(5);
t82 = t61 * t94;
t81 = t60 * t94;
t44 = t58 * t65 + t62 * t68;
t76 = t59 * t84;
t21 = -t44 * t61 - t64 * t76;
t18 = t21 * pkin(3);
t80 = t108 * t21 + t18;
t46 = -t62 * t77 + t84 * t65;
t23 = -t46 * t61 + t58 * t95;
t19 = t23 * pkin(3);
t79 = t108 * t23 + t19;
t47 = -t61 * t96 + t85 * t64;
t42 = t47 * pkin(3);
t78 = t108 * t47 + t42;
t75 = t59 * pkin(3) * t91 + pkin(9) * t82 + t90;
t29 = (t60 * t62 + t63 * t91) * t59;
t74 = t29 * pkin(4) + t75;
t73 = rSges(4,1) * t64 - rSges(4,2) * t61;
t72 = rSges(5,1) * t63 - rSges(5,2) * t60;
t71 = -rSges(6,2) * t63 + rSges(6,3) * t60;
t40 = t43 * pkin(2);
t70 = t44 * pkin(8) - pkin(9) * t98 - t43 * t106 - t40;
t41 = t45 * pkin(2);
t69 = t46 * pkin(8) - pkin(9) * t97 - t45 * t106 - t41;
t11 = -t43 * t92 + t44 * t60;
t67 = t11 * pkin(4) + t70;
t13 = -t45 * t92 + t46 * t60;
t66 = t13 * pkin(4) + t69;
t48 = t85 * t61 + t62 * t95;
t28 = -t63 * t96 + t64 * t81;
t26 = t48 * t63 - t81;
t25 = t48 * t60 + t63 * t94;
t24 = t58 * t59 * t61 + t46 * t64;
t22 = t44 * t64 - t61 * t76;
t20 = t25 * pkin(4);
t12 = -t45 * t93 - t46 * t63;
t10 = -t43 * t93 - t44 * t63;
t7 = t24 * t63 + t45 * t60;
t6 = t24 * t60 - t45 * t63;
t5 = t22 * t63 + t43 * t60;
t4 = t22 * t60 - t43 * t63;
t3 = t6 * pkin(4);
t2 = t4 * pkin(4);
t1 = [(-m(2) - m(3) - m(4) - m(5) + t107) * g(3), -m(3) * (g(1) * (-t45 * rSges(3,1) - t46 * rSges(3,2)) + g(2) * (-t43 * rSges(3,1) - t44 * rSges(3,2)) + (rSges(3,1) * t65 - rSges(3,2) * t62) * t104) - m(4) * (g(1) * (t101 * t46 - t73 * t45 - t41) + g(2) * (t101 * t44 - t73 * t43 - t40) + g(3) * t90 + (rSges(4,3) * t62 + t73 * t65) * t104) - m(5) * (g(1) * (t13 * rSges(5,1) - t12 * rSges(5,2) - rSges(5,3) * t97 + t69) + g(2) * (t11 * rSges(5,1) - t10 * rSges(5,2) - rSges(5,3) * t98 + t70) + g(3) * (t29 * rSges(5,1) - t28 * rSges(5,2) + rSges(5,3) * t82 + t75)) - m(6) * (g(1) * (-rSges(6,1) * t97 - t13 * rSges(6,2) + t87 * t12 + t66) + g(2) * (-rSges(6,1) * t98 - t11 * rSges(6,2) + t87 * t10 + t67) + g(3) * (rSges(6,1) * t82 - t29 * rSges(6,2) + t87 * t28 + t74)) - m(7) * (g(1) * (t88 * t12 + t86 * t13 + t66) + g(2) * (t88 * t10 + t86 * t11 + t67) + g(3) * (t88 * t28 + t86 * t29 + t74) + (-g(1) * t45 - g(2) * t43 + g(3) * t94) * t61 * t102) -m(4) * (g(1) * (t23 * rSges(4,1) - t24 * rSges(4,2)) + g(2) * (t21 * rSges(4,1) - t22 * rSges(4,2)) + g(3) * (t47 * rSges(4,1) - t48 * rSges(4,2))) - m(5) * (g(1) * (t100 * t24 + t72 * t23 + t19) + g(2) * (t100 * t22 + t72 * t21 + t18) + g(3) * (t100 * t48 + t72 * t47 + t42)) - m(6) * (g(1) * (t103 * t24 + t71 * t23 + t79) + g(2) * (t103 * t22 + t71 * t21 + t80) + g(3) * (t103 * t48 + t71 * t47 + t78)) + (-g(1) * (t23 * t99 + t79) - g(2) * (t21 * t99 + t80) - g(3) * (t47 * t99 + t78) - (g(1) * t23 + g(2) * t21 + g(3) * t47) * t63 * t86 - (g(1) * t24 + g(2) * t22 + g(3) * t48) * (pkin(9) + t102)) * m(7), -m(5) * (g(1) * (-t6 * rSges(5,1) - t7 * rSges(5,2)) + g(2) * (-t4 * rSges(5,1) - t5 * rSges(5,2)) + g(3) * (-t25 * rSges(5,1) - t26 * rSges(5,2))) - m(6) * (g(1) * (t6 * rSges(6,2) + t87 * t7 - t3) + g(2) * (t4 * rSges(6,2) + t87 * t5 - t2) + g(3) * (t25 * rSges(6,2) + t87 * t26 - t20)) - m(7) * (g(1) * (-t86 * t6 + t88 * t7 - t3) + g(2) * (-t86 * t4 + t88 * t5 - t2) + g(3) * (-t86 * t25 + t88 * t26 - t20)) t107 * (g(1) * t6 + g(2) * t4 + g(3) * t25) -m(7) * (g(1) * t7 + g(2) * t5 + g(3) * t26)];
taug  = t1(:);
