% Calculate Gravitation load on the joints for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:39
% EndTime: 2019-03-09 16:21:43
% DurationCPUTime: 1.28s
% Computational Cost: add. (661->188), mult. (1477->260), div. (0->0), fcn. (1760->12), ass. (0->74)
t84 = pkin(10) + qJ(5) + rSges(7,3);
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t103 = -pkin(3) * t53 - qJ(4) * t50;
t47 = sin(pkin(6));
t95 = g(3) * t47;
t81 = qJ(5) + rSges(6,3);
t52 = sin(qJ(1));
t51 = sin(qJ(2));
t80 = cos(pkin(6));
t71 = t51 * t80;
t90 = cos(qJ(2));
t91 = cos(qJ(1));
t27 = t52 * t90 + t71 * t91;
t76 = t47 * t91;
t10 = t27 * t53 - t50 * t76;
t29 = -t52 * t71 + t90 * t91;
t88 = t47 * t52;
t14 = t29 * t53 + t50 * t88;
t87 = t47 * t53;
t25 = t50 * t80 + t51 * t87;
t102 = g(1) * t14 + g(2) * t10 + g(3) * t25;
t101 = -m(6) - m(7);
t46 = sin(pkin(11));
t99 = pkin(5) * t46;
t94 = rSges(5,1) + pkin(9);
t93 = rSges(4,3) + pkin(9);
t48 = cos(pkin(11));
t41 = pkin(5) * t48 + pkin(4);
t92 = pkin(9) + t41;
t89 = t47 * t51;
t86 = t91 * pkin(1) + pkin(8) * t88;
t75 = t47 * t90;
t85 = pkin(2) * t75 + pkin(9) * t89;
t82 = rSges(5,3) + qJ(4);
t65 = t80 * t90;
t26 = t51 * t52 - t65 * t91;
t20 = t26 * pkin(2);
t79 = t103 * t26 - t20;
t28 = t51 * t91 + t52 * t65;
t22 = t28 * pkin(2);
t78 = t103 * t28 - t22;
t77 = t29 * pkin(2) + t86;
t74 = t50 * t90;
t73 = t53 * t90;
t72 = -t52 * pkin(1) + pkin(8) * t76;
t70 = t14 * pkin(3) + t77;
t68 = t47 * t73;
t69 = t47 * qJ(4) * t74 + pkin(3) * t68 + t85;
t67 = t46 * t74;
t66 = -t27 * pkin(2) + t72;
t64 = -rSges(4,1) * t53 + rSges(4,2) * t50;
t63 = -t46 * rSges(6,1) - t48 * rSges(6,2);
t62 = rSges(5,2) * t53 - rSges(5,3) * t50;
t61 = -pkin(3) * t10 + t66;
t60 = qJ(4) - t63;
t59 = rSges(6,1) * t48 - rSges(6,2) * t46 + pkin(4) + pkin(9);
t45 = pkin(11) + qJ(6);
t42 = sin(t45);
t43 = cos(t45);
t58 = t43 * rSges(7,1) - t42 * rSges(7,2) + t92;
t9 = t27 * t50 + t53 * t76;
t57 = -rSges(7,1) * t42 - rSges(7,2) * t43 - t99;
t56 = qJ(4) - t57;
t55 = t50 * t63 - t53 * t81;
t54 = t50 * t57 - t53 * t84;
t24 = t50 * t89 - t53 * t80;
t19 = t24 * pkin(3);
t13 = t29 * t50 - t52 * t87;
t7 = t13 * pkin(3);
t5 = t9 * pkin(3);
t4 = t13 * t42 + t28 * t43;
t3 = t13 * t43 - t28 * t42;
t1 = [-m(2) * (g(1) * (-t52 * rSges(2,1) - rSges(2,2) * t91) + g(2) * (rSges(2,1) * t91 - t52 * rSges(2,2))) - m(3) * (g(1) * (-t27 * rSges(3,1) + t26 * rSges(3,2) + rSges(3,3) * t76 + t72) + g(2) * (rSges(3,1) * t29 - rSges(3,2) * t28 + rSges(3,3) * t88 + t86)) - m(4) * (g(1) * (-rSges(4,1) * t10 + rSges(4,2) * t9 - t26 * t93 + t66) + g(2) * (rSges(4,1) * t14 - rSges(4,2) * t13 + t28 * t93 + t77)) - m(5) * (g(1) * (rSges(5,2) * t10 - t26 * t94 - t82 * t9 + t61) + g(2) * (-rSges(5,2) * t14 + t13 * t82 + t28 * t94 + t70)) - m(6) * (g(1) * (-t10 * t81 - t26 * t59 - t60 * t9 + t61) + g(2) * (t13 * t60 + t14 * t81 + t28 * t59 + t70)) - m(7) * (g(1) * (-t10 * t84 - t26 * t58 - t56 * t9 + t61) + g(2) * (rSges(7,1) * t4 + rSges(7,2) * t3 + t92 * t28 + t84 * t14 + (qJ(4) + t99) * t13 + t70)) -m(3) * (g(1) * (-rSges(3,1) * t28 - rSges(3,2) * t29) + g(2) * (-rSges(3,1) * t26 - rSges(3,2) * t27) + (rSges(3,1) * t90 - rSges(3,2) * t51) * t95) - m(4) * (g(1) * (t28 * t64 + t29 * t93 - t22) + g(2) * (t26 * t64 + t27 * t93 - t20) + g(3) * t85 + (rSges(4,1) * t73 - rSges(4,2) * t74 + rSges(4,3) * t51) * t95) - m(5) * (g(1) * (t28 * t62 + t29 * t94 + t78) + g(2) * (t26 * t62 + t27 * t94 + t79) + g(3) * t69 + (rSges(5,1) * t51 - rSges(5,2) * t73 + rSges(5,3) * t74) * t95) - m(6) * (g(1) * (t28 * t55 + t29 * t59 + t78) + g(2) * (t26 * t55 + t27 * t59 + t79) + g(3) * (pkin(4) * t89 + t69 + t81 * t68 + ((t48 * t51 + t67) * rSges(6,1) + (-t46 * t51 + t48 * t74) * rSges(6,2)) * t47)) - m(7) * (g(1) * (t28 * t54 + t29 * t58 + t78) + g(2) * (t26 * t54 + t27 * t58 + t79) + g(3) * (t41 * t89 + t84 * t68 + t69) + (pkin(5) * t67 + (t42 * t74 + t43 * t51) * rSges(7,1) + (-t42 * t51 + t43 * t74) * rSges(7,2)) * t95) -m(4) * (g(1) * (-rSges(4,1) * t13 - rSges(4,2) * t14) + g(2) * (-rSges(4,1) * t9 - rSges(4,2) * t10) + g(3) * (-rSges(4,1) * t24 - rSges(4,2) * t25)) - m(5) * (g(1) * (rSges(5,2) * t13 + t14 * t82 - t7) + g(2) * (rSges(5,2) * t9 + t10 * t82 - t5) + g(3) * (rSges(5,2) * t24 + t25 * t82 - t19)) + (-g(1) * (-t13 * t84 - t7) - g(2) * (-t84 * t9 - t5) - g(3) * (-t24 * t84 - t19) - t102 * t56) * m(7) + (-g(1) * (-t13 * t81 - t7) - g(2) * (-t81 * t9 - t5) - g(3) * (-t24 * t81 - t19) - t102 * t60) * m(6) (-m(5) + t101) * (g(1) * t13 + g(2) * t9 + g(3) * t24) t101 * t102, -m(7) * (g(1) * (rSges(7,1) * t3 - rSges(7,2) * t4) + g(2) * ((-t26 * t42 + t43 * t9) * rSges(7,1) + (-t26 * t43 - t42 * t9) * rSges(7,2)) + g(3) * ((t24 * t43 + t42 * t75) * rSges(7,1) + (-t24 * t42 + t43 * t75) * rSges(7,2)))];
taug  = t1(:);
