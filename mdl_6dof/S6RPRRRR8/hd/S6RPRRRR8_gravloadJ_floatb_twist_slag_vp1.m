% Calculate Gravitation load on the joints for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:19:51
% EndTime: 2019-03-09 07:19:53
% DurationCPUTime: 0.80s
% Computational Cost: add. (390->146), mult. (450->192), div. (0->0), fcn. (407->10), ass. (0->70)
t100 = rSges(6,3) + pkin(9);
t44 = qJ(3) + qJ(4);
t39 = cos(t44);
t101 = t100 * t39;
t50 = cos(qJ(1));
t90 = g(2) * t50;
t37 = sin(t44);
t51 = -pkin(10) - pkin(9);
t43 = qJ(5) + qJ(6);
t36 = sin(t43);
t84 = rSges(7,2) * t36;
t99 = t37 * t51 + t39 * t84;
t98 = (-rSges(7,3) + t51) * t39;
t47 = sin(qJ(1));
t97 = g(1) * t47 - t90;
t38 = cos(t43);
t73 = t50 * t38;
t78 = t47 * t36;
t5 = -t37 * t78 + t73;
t74 = t50 * t36;
t77 = t47 * t38;
t6 = t37 * t77 + t74;
t96 = t5 * rSges(7,1) - t6 * rSges(7,2);
t7 = t37 * t74 + t77;
t8 = t37 * t73 - t78;
t95 = t7 * rSges(7,1) + t8 * rSges(7,2);
t49 = cos(qJ(3));
t94 = pkin(3) * t49;
t45 = sin(qJ(5));
t93 = pkin(5) * t45;
t85 = rSges(6,2) * t45;
t69 = t39 * t85;
t91 = t69 * t90;
t89 = g(3) * t39;
t46 = sin(qJ(3));
t88 = t46 * pkin(3);
t87 = rSges(4,3) + pkin(7);
t83 = t37 * t47;
t82 = t37 * t50;
t80 = t39 * t47;
t76 = t47 * t45;
t48 = cos(qJ(5));
t75 = t47 * t48;
t72 = t50 * t45;
t71 = t50 * t48;
t70 = t50 * pkin(1) + t47 * qJ(2);
t67 = t47 * t88 + t70;
t41 = t50 * qJ(2);
t52 = -pkin(8) - pkin(7);
t66 = t47 * t52 + t50 * t88 + t41;
t65 = t99 * t90;
t64 = -rSges(6,1) * t48 - pkin(4);
t34 = t48 * pkin(5) + pkin(4);
t63 = -rSges(7,1) * t38 - t34;
t62 = rSges(5,1) * t80 - rSges(5,2) * t83;
t60 = t46 * rSges(4,1) + t49 * rSges(4,2);
t59 = t37 * rSges(5,1) + t39 * rSges(5,2);
t58 = -rSges(7,1) * t36 - rSges(7,2) * t38;
t11 = t37 * t72 + t75;
t9 = -t37 * t76 + t71;
t57 = t37 * t34 + t98;
t56 = t39 * rSges(6,1) * t75 + pkin(4) * t80 + t100 * t83 - t47 * t69;
t55 = t101 + (t64 + t85) * t37;
t54 = -t98 + (t63 + t84) * t37;
t53 = t39 * rSges(7,1) * t77 + rSges(7,3) * t83 + t34 * t80 - t47 * t99;
t29 = t47 * t94;
t25 = rSges(5,2) * t82;
t12 = t37 * t71 - t76;
t10 = t37 * t75 + t72;
t1 = [-m(2) * (g(1) * (-t47 * rSges(2,1) - t50 * rSges(2,2)) + g(2) * (t50 * rSges(2,1) - t47 * rSges(2,2))) - m(3) * (g(1) * (t50 * rSges(3,3) + t41 + (rSges(3,2) - pkin(1)) * t47) + g(2) * (-t50 * rSges(3,2) + t47 * rSges(3,3) + t70)) - m(4) * (g(1) * t41 + g(2) * t70 + (g(1) * t60 + g(2) * t87) * t50 + (g(1) * (-pkin(1) - t87) + g(2) * t60) * t47) - m(5) * (g(1) * t66 + g(2) * t67 + (g(1) * t59 + g(2) * (rSges(5,3) - t52)) * t50 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t59) * t47) - m(6) * (g(1) * (t12 * rSges(6,1) - t11 * rSges(6,2) - t47 * pkin(1) + pkin(4) * t82 + t66) + g(2) * (t10 * rSges(6,1) + t9 * rSges(6,2) + pkin(4) * t83 - t50 * t52 + t67) - (g(1) * t50 + g(2) * t47) * t101) - m(7) * (g(1) * (t8 * rSges(7,1) - t7 * rSges(7,2) + t66) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t67) + (g(1) * t57 + g(2) * (-t52 + t93)) * t50 + (g(1) * (-pkin(1) - t93) + g(2) * t57) * t47) (-m(3) - m(4) - m(5) - m(6) - m(7)) * t97, -m(4) * (-g(3) * t60 + t97 * (rSges(4,1) * t49 - rSges(4,2) * t46)) - m(5) * (g(1) * (t29 + t62) + g(2) * (t25 + (-rSges(5,1) * t39 - t94) * t50) + g(3) * (-t59 - t88)) - m(6) * (g(1) * (t29 + t56) + t91 + g(3) * (t55 - t88) + (-t100 * t37 + t64 * t39 - t94) * t90) - m(7) * (g(1) * (t29 + t53) + t65 + g(3) * (t54 - t88) + (-rSges(7,3) * t37 + t63 * t39 - t94) * t90) -m(5) * (g(1) * t62 + g(2) * t25 - g(3) * t59) - m(6) * (g(1) * t56 + g(3) * t55 + t91) - m(7) * (g(1) * t53 + g(3) * t54 + t65) + ((m(6) * t100 + m(7) * rSges(7,3)) * t37 + (m(5) * rSges(5,1) - m(6) * t64 - m(7) * t63) * t39) * t90, -m(6) * (g(1) * (t9 * rSges(6,1) - t10 * rSges(6,2)) + g(2) * (t11 * rSges(6,1) + t12 * rSges(6,2))) - m(7) * (g(1) * (t9 * pkin(5) + t96) + g(2) * (t11 * pkin(5) + t95)) + (-m(6) * (-rSges(6,1) * t45 - rSges(6,2) * t48) - m(7) * (t58 - t93)) * t89, -m(7) * (g(1) * t96 + g(2) * t95 + t58 * t89)];
taug  = t1(:);
