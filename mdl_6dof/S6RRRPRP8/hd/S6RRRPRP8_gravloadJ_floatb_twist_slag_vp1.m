% Calculate Gravitation load on the joints for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:14:27
% EndTime: 2019-03-09 17:14:31
% DurationCPUTime: 1.49s
% Computational Cost: add. (434->203), mult. (983->271), div. (0->0), fcn. (1071->8), ass. (0->76)
t49 = sin(qJ(1));
t47 = sin(qJ(3));
t53 = cos(qJ(1));
t81 = t53 * t47;
t51 = cos(qJ(3));
t52 = cos(qJ(2));
t83 = t51 * t52;
t26 = t49 * t83 - t81;
t50 = cos(qJ(5));
t84 = t49 * t52;
t25 = t47 * t84 + t51 * t53;
t46 = sin(qJ(5));
t93 = t25 * t46;
t2 = t26 * t50 + t93;
t4 = t25 * t50 - t26 * t46;
t111 = rSges(6,1) * t4 - rSges(6,2) * t2;
t110 = rSges(7,1) * t4 - rSges(7,2) * t2;
t27 = -t49 * t51 + t52 * t81;
t82 = t52 * t53;
t28 = t49 * t47 + t51 * t82;
t10 = t27 * t50 - t28 * t46;
t62 = t27 * t46 + t28 * t50;
t107 = rSges(6,1) * t10 - rSges(6,2) * t62;
t106 = rSges(7,1) * t10 - rSges(7,2) * t62;
t48 = sin(qJ(2));
t91 = t46 * t47;
t60 = t50 * t51 + t91;
t103 = t60 * t48;
t98 = g(2) * t49;
t102 = g(1) * t53 + t98;
t101 = -pkin(3) - pkin(4);
t100 = g(1) * t49;
t97 = g(3) * t48;
t42 = t52 * pkin(2);
t96 = -rSges(5,1) - pkin(3);
t95 = -rSges(6,3) - pkin(9);
t39 = pkin(5) * t50 + pkin(4);
t94 = -pkin(3) - t39;
t90 = t46 * t51;
t89 = t47 * t48;
t88 = t47 * t50;
t87 = t47 * t52;
t86 = t48 * t51;
t85 = t48 * t53;
t80 = -rSges(7,3) - qJ(6) - pkin(9);
t79 = t48 * pkin(8) + t42;
t78 = t53 * pkin(1) + t49 * pkin(7);
t77 = qJ(4) * t47;
t76 = rSges(5,3) + qJ(4);
t75 = -pkin(1) - t42;
t74 = t95 * t53;
t73 = t80 * t53;
t72 = pkin(5) * t46 + qJ(4);
t70 = pkin(3) * t83 + t52 * t77 + t79;
t69 = pkin(2) * t82 + pkin(8) * t85 + t78;
t43 = t53 * pkin(7);
t68 = -t26 * pkin(3) - qJ(4) * t25 + t43;
t67 = t28 * pkin(3) + t69;
t66 = rSges(3,1) * t52 - rSges(3,2) * t48;
t17 = t46 * t86 - t48 * t88;
t64 = -rSges(6,1) * t17 - rSges(6,2) * t103;
t63 = -rSges(7,1) * t17 - rSges(7,2) * t103;
t61 = -t88 + t90;
t58 = t48 * t61;
t35 = pkin(8) * t82;
t32 = pkin(8) * t84;
t30 = qJ(4) * t86;
t23 = t27 * pkin(3);
t21 = t25 * pkin(3);
t20 = t60 * t52;
t19 = t61 * t52;
t15 = t53 * t103;
t14 = t53 * t58;
t13 = t49 * t103;
t12 = t49 * t58;
t1 = [-m(2) * (g(1) * (-t49 * rSges(2,1) - rSges(2,2) * t53) + g(2) * (rSges(2,1) * t53 - t49 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t53 + t43) + g(2) * (rSges(3,1) * t82 - rSges(3,2) * t85 + t78) + (g(1) * (-pkin(1) - t66) + g(2) * rSges(3,3)) * t49) - m(4) * (g(1) * (-rSges(4,1) * t26 + rSges(4,2) * t25 + t43) + g(2) * (t28 * rSges(4,1) - t27 * rSges(4,2) + rSges(4,3) * t85 + t69) + ((-rSges(4,3) - pkin(8)) * t48 + t75) * t100) - m(5) * (g(1) * (-rSges(5,1) * t26 - rSges(5,3) * t25 + t68) + g(2) * (t28 * rSges(5,1) + rSges(5,2) * t85 + t27 * t76 + t67) + ((-rSges(5,2) - pkin(8)) * t48 + t75) * t100) - m(6) * (g(1) * (-rSges(6,1) * t2 - rSges(6,2) * t4 - pkin(4) * t26 + t68) + g(2) * (rSges(6,1) * t62 + t10 * rSges(6,2) + t28 * pkin(4) + t27 * qJ(4) + t48 * t74 + t67) + ((-pkin(8) - t95) * t48 + t75) * t100) - m(7) * (g(1) * (-rSges(7,1) * t2 - rSges(7,2) * t4 - pkin(5) * t93 - t26 * t39 + t68) + g(2) * (rSges(7,1) * t62 + t10 * rSges(7,2) + t27 * t72 + t28 * t39 + t48 * t73 + t67) + ((-pkin(8) - t80) * t48 + t75) * t100) -m(3) * (g(3) * t66 + t102 * (-rSges(3,1) * t48 - rSges(3,2) * t52)) - m(4) * (g(1) * (rSges(4,3) * t82 + t35) + g(2) * (rSges(4,3) * t84 + t32) + g(3) * (rSges(4,1) * t83 - rSges(4,2) * t87 + t79) + (g(3) * rSges(4,3) + t102 * (-rSges(4,1) * t51 + rSges(4,2) * t47 - pkin(2))) * t48) - m(5) * (g(1) * (rSges(5,2) * t82 + t35) + g(2) * (rSges(5,2) * t84 + t32) + g(3) * (rSges(5,1) * t83 + rSges(5,3) * t87 + t70) + (g(3) * rSges(5,2) + t102 * (-t47 * t76 + t51 * t96 - pkin(2))) * t48) - m(6) * (g(1) * (-t15 * rSges(6,1) + t14 * rSges(6,2) + t35) + g(2) * (-rSges(6,1) * t13 + rSges(6,2) * t12 + t32) + g(3) * (rSges(6,1) * t20 - t19 * rSges(6,2) + t70) + (g(3) * pkin(4) * t51 + g(1) * t74 + t95 * t98) * t52 + (g(3) * t95 + t102 * (t101 * t51 - pkin(2) - t77)) * t48) - m(7) * (g(1) * (-t15 * rSges(7,1) + t14 * rSges(7,2) + t35) + g(2) * (-rSges(7,1) * t13 + rSges(7,2) * t12 + t32) + g(3) * (rSges(7,1) * t20 - rSges(7,2) * t19 + t70) + (g(3) * (pkin(5) * t91 + t39 * t51) + g(1) * t73 + t80 * t98) * t52 + (g(3) * t80 + t102 * (-t47 * t72 + t51 * t94 - pkin(2))) * t48) -m(4) * (g(1) * (-rSges(4,1) * t27 - rSges(4,2) * t28) + g(2) * (-rSges(4,1) * t25 - rSges(4,2) * t26) + (-rSges(4,1) * t47 - rSges(4,2) * t51) * t97) - m(5) * (g(1) * (-rSges(5,1) * t27 + t28 * t76 - t23) + g(2) * (-rSges(5,1) * t25 + t26 * t76 - t21) + g(3) * t30 + (rSges(5,3) * t51 + t47 * t96) * t97) - m(6) * (g(1) * (-pkin(4) * t27 + qJ(4) * t28 - t107 - t23) + g(2) * (-pkin(4) * t25 + qJ(4) * t26 - t111 - t21) + g(3) * (t101 * t89 + t30 - t64)) - m(7) * (g(1) * (-t27 * t39 + t28 * t72 - t106 - t23) + g(2) * (-t25 * t39 + t26 * t72 - t110 - t21) + g(3) * (t30 - t63) + (pkin(5) * t90 + t47 * t94) * t97) (-m(5) - m(6) - m(7)) * (g(1) * t27 + g(2) * t25 + g(3) * t89) -m(6) * (g(1) * t107 + g(2) * t111 + g(3) * t64) + (-g(1) * t106 - g(2) * t110 - g(3) * t63 - (g(1) * t10 + g(2) * t4 - t61 * t97) * pkin(5)) * m(7), -m(7) * (g(3) * t52 - t102 * t48)];
taug  = t1(:);
