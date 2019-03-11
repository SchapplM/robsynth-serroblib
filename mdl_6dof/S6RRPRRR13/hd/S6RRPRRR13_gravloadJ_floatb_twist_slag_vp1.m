% Calculate Gravitation load on the joints for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR13_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:44:59
% EndTime: 2019-03-09 14:45:02
% DurationCPUTime: 0.97s
% Computational Cost: add. (623->171), mult. (1343->238), div. (0->0), fcn. (1588->12), ass. (0->68)
t94 = rSges(5,3) + pkin(9);
t49 = sin(qJ(2));
t50 = sin(qJ(1));
t53 = cos(qJ(2));
t54 = cos(qJ(1));
t82 = cos(pkin(6));
t78 = t54 * t82;
t28 = t49 * t78 + t50 * t53;
t79 = t50 * t82;
t30 = -t49 * t79 + t53 * t54;
t109 = g(1) * t30 + g(2) * t28;
t27 = t49 * t50 - t53 * t78;
t29 = t54 * t49 + t53 * t79;
t108 = g(1) * t29 + g(2) * t27;
t47 = sin(qJ(5));
t51 = cos(qJ(5));
t48 = sin(qJ(4));
t52 = cos(qJ(4));
t46 = sin(pkin(6));
t87 = t46 * t54;
t69 = -t27 * t48 + t52 * t87;
t73 = t28 * t51 + t47 * t69;
t107 = -t28 * t47 + t51 * t69;
t89 = t46 * t50;
t13 = t29 * t48 + t52 * t89;
t88 = t46 * t53;
t26 = -t48 * t88 + t82 * t52;
t106 = g(1) * t13 - g(2) * t69 + g(3) * t26;
t12 = -t29 * t52 + t48 * t89;
t25 = -t82 * t48 - t52 * t88;
t68 = t27 * t52 + t48 * t87;
t105 = -g(1) * t12 + g(2) * t68 + g(3) * t25;
t104 = pkin(5) * t47;
t95 = g(3) * t46;
t93 = rSges(6,3) + pkin(10);
t90 = t46 * t49;
t86 = rSges(7,3) + pkin(11) + pkin(10);
t85 = pkin(2) * t88 + qJ(3) * t90;
t84 = t54 * pkin(1) + pkin(8) * t89;
t83 = rSges(4,3) + qJ(3);
t81 = t30 * pkin(2) + t84;
t80 = -t50 * pkin(1) + pkin(8) * t87;
t77 = g(3) * (pkin(9) * t88 + t85);
t76 = -t28 * pkin(2) + t80;
t75 = rSges(5,1) * t48 + rSges(5,2) * t52;
t74 = rSges(6,1) * t47 + rSges(6,2) * t51;
t7 = -t13 * t47 + t30 * t51;
t72 = rSges(6,1) * t51 - rSges(6,2) * t47 + pkin(4);
t70 = -t26 * t47 + t51 * t90;
t41 = pkin(5) * t51 + pkin(4);
t45 = qJ(5) + qJ(6);
t42 = sin(t45);
t43 = cos(t45);
t67 = rSges(7,1) * t43 - rSges(7,2) * t42 + t41;
t65 = pkin(3) * t89 + t29 * qJ(3) + t81;
t64 = rSges(7,1) * t42 + rSges(7,2) * t43 + t104;
t63 = pkin(3) * t87 - t27 * qJ(3) + t76;
t62 = -pkin(9) - t64;
t21 = t27 * pkin(2);
t23 = t29 * pkin(2);
t61 = -g(1) * t23 - g(2) * t21 + t77;
t60 = t72 * t48 - t93 * t52;
t59 = t67 * t48 - t86 * t52;
t5 = -t13 * t42 + t30 * t43;
t6 = t13 * t43 + t30 * t42;
t58 = m(7) * (g(1) * (t5 * rSges(7,1) - t6 * rSges(7,2)) + g(2) * ((t28 * t43 + t42 * t69) * rSges(7,1) + (-t28 * t42 + t43 * t69) * rSges(7,2)) + g(3) * ((-t26 * t42 + t43 * t90) * rSges(7,1) + (-t26 * t43 - t42 * t90) * rSges(7,2)));
t8 = t13 * t51 + t30 * t47;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t50 - rSges(2,2) * t54) + g(2) * (rSges(2,1) * t54 - rSges(2,2) * t50)) - m(3) * (g(1) * (-rSges(3,1) * t28 + rSges(3,2) * t27 + rSges(3,3) * t87 + t80) + g(2) * (rSges(3,1) * t30 - rSges(3,2) * t29 + rSges(3,3) * t89 + t84)) - m(4) * (g(1) * (rSges(4,1) * t87 + rSges(4,2) * t28 - t83 * t27 + t76) + g(2) * (rSges(4,1) * t89 - rSges(4,2) * t30 + t83 * t29 + t81)) - m(5) * (g(1) * (rSges(5,1) * t69 - rSges(5,2) * t68 - t94 * t28 + t63) + g(2) * (rSges(5,1) * t13 - rSges(5,2) * t12 + t94 * t30 + t65)) - m(6) * (g(1) * (t107 * rSges(6,1) - t73 * rSges(6,2) + t69 * pkin(4) - t28 * pkin(9) + t93 * t68 + t63) + g(2) * (rSges(6,1) * t8 + rSges(6,2) * t7 + pkin(4) * t13 + t30 * pkin(9) + t93 * t12 + t65)) - m(7) * (g(1) * (t62 * t28 + t67 * t69 + t68 * t86 + t63) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t13 * t41 + (pkin(9) + t104) * t30 + t86 * t12 + t65)) -m(3) * (g(1) * (-rSges(3,1) * t29 - rSges(3,2) * t30) + g(2) * (-rSges(3,1) * t27 - rSges(3,2) * t28) + (rSges(3,1) * t53 - rSges(3,2) * t49) * t95) - m(4) * (g(1) * (rSges(4,2) * t29 + t83 * t30 - t23) + g(2) * (rSges(4,2) * t27 + t83 * t28 - t21) + g(3) * ((-rSges(4,2) * t53 + rSges(4,3) * t49) * t46 + t85)) - m(5) * (g(1) * (-t94 * t29 - t23) + g(2) * (-t94 * t27 - t21) + t77 + (rSges(5,3) * t53 + t75 * t49) * t95 + t109 * (qJ(3) + t75)) - m(6) * ((t60 * t49 + t74 * t53) * t95 + t61 + t108 * (-pkin(9) - t74) + t109 * (qJ(3) + t60)) - m(7) * ((t59 * t49 + t64 * t53) * t95 + t61 + t108 * t62 + t109 * (qJ(3) + t59)) (-m(4) - m(5) - m(6) - m(7)) * (-g(3) * t88 + t108) -m(5) * (g(1) * (-rSges(5,1) * t12 - rSges(5,2) * t13) + g(2) * (rSges(5,1) * t68 + rSges(5,2) * t69) + g(3) * (rSges(5,1) * t25 - rSges(5,2) * t26)) - m(6) * (t105 * t72 + t106 * t93) - m(7) * (t105 * t67 + t106 * t86) -m(6) * (g(1) * (rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (t73 * rSges(6,1) + t107 * rSges(6,2)) + g(3) * (t70 * rSges(6,1) + (-t26 * t51 - t47 * t90) * rSges(6,2))) - t58 - m(7) * (g(1) * t7 + g(2) * t73 + g(3) * t70) * pkin(5), -t58];
taug  = t1(:);
