% Calculate Gravitation load on the joints for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:13
% EndTime: 2019-03-09 09:01:16
% DurationCPUTime: 1.20s
% Computational Cost: add. (680->157), mult. (1628->227), div. (0->0), fcn. (2020->12), ass. (0->67)
t51 = sin(qJ(2));
t55 = cos(qJ(2));
t81 = sin(pkin(11));
t82 = cos(pkin(11));
t38 = -t51 * t82 - t55 * t81;
t48 = cos(pkin(6));
t111 = t48 * t38;
t52 = sin(qJ(1));
t56 = cos(qJ(1));
t63 = -t51 * t81 + t55 * t82;
t17 = t111 * t56 - t52 * t63;
t22 = t111 * t52 + t63 * t56;
t47 = sin(pkin(6));
t31 = t38 * t47;
t114 = -g(1) * t22 + g(2) * t17 + g(3) * t31;
t59 = t48 * t63;
t21 = t38 * t56 - t52 * t59;
t87 = t52 * t55;
t91 = t51 * t56;
t35 = -t48 * t87 - t91;
t65 = t35 * pkin(2);
t62 = t21 * pkin(3) + t65;
t86 = t55 * t56;
t88 = t52 * t51;
t113 = t48 * t86 - t88;
t98 = rSges(6,3) + pkin(9);
t18 = t52 * t38 + t56 * t59;
t30 = t63 * t47;
t110 = -g(1) * t21 - g(2) * t18 - g(3) * t30;
t49 = sin(qJ(6));
t53 = cos(qJ(6));
t50 = sin(qJ(5));
t54 = cos(qJ(5));
t92 = t47 * t56;
t68 = t18 * t50 + t54 * t92;
t109 = -t17 * t53 + t49 * t68;
t108 = t17 * t49 + t53 * t68;
t106 = pkin(2) * t55;
t97 = rSges(7,3) + pkin(10);
t93 = t47 * t52;
t46 = pkin(1) + t106;
t89 = t52 * t46;
t41 = t56 * t46;
t85 = t22 * pkin(3) + t41;
t43 = t47 * t106;
t84 = t30 * pkin(3) + t43;
t83 = rSges(5,3) + qJ(4);
t80 = -m(5) - m(6) - m(7);
t32 = pkin(2) * t48 * t51 + (-pkin(8) - qJ(3)) * t47;
t77 = rSges(5,1) * t47 - t32;
t76 = rSges(4,3) * t47 - t32;
t75 = pkin(3) * t17 - t89;
t72 = t113 * pkin(2);
t71 = t18 * pkin(3) + t72;
t70 = rSges(7,1) * t53 - rSges(7,2) * t49 + pkin(5);
t67 = -t18 * t54 + t50 * t92;
t64 = pkin(4) * t93 - qJ(4) * t21 - t32 * t52 + t85;
t61 = pkin(4) * t92 + t18 * qJ(4) - t56 * t32 + t75;
t36 = -t48 * t88 + t86;
t34 = -t48 * t91 - t87;
t24 = -t30 * t50 + t48 * t54;
t23 = -t30 * t54 - t48 * t50;
t5 = -t21 * t50 + t54 * t93;
t4 = t21 * t54 + t50 * t93;
t3 = t22 * t49 + t5 * t53;
t2 = t22 * t53 - t49 * t5;
t1 = [-m(2) * (g(1) * (-t52 * rSges(2,1) - rSges(2,2) * t56) + g(2) * (rSges(2,1) * t56 - t52 * rSges(2,2))) - m(3) * (g(1) * (t34 * rSges(3,1) - rSges(3,2) * t113 - t52 * pkin(1)) + g(2) * (t36 * rSges(3,1) + t35 * rSges(3,2) + pkin(1) * t56) + (g(1) * t56 + g(2) * t52) * t47 * (rSges(3,3) + pkin(8))) - m(4) * (g(1) * (rSges(4,1) * t17 - t18 * rSges(4,2) + t76 * t56 - t89) + g(2) * (rSges(4,1) * t22 + rSges(4,2) * t21 + t76 * t52 + t41)) - m(5) * (g(1) * (-rSges(5,2) * t17 + t83 * t18 + t77 * t56 + t75) + g(2) * (-rSges(5,2) * t22 - t83 * t21 + t77 * t52 + t85)) - m(6) * (g(1) * (rSges(6,1) * t68 - rSges(6,2) * t67 + t17 * t98 + t61) + g(2) * (rSges(6,1) * t5 - rSges(6,2) * t4 + t98 * t22 + t64)) - m(7) * (g(1) * (t108 * rSges(7,1) - t109 * rSges(7,2) + t68 * pkin(5) + pkin(9) * t17 + t97 * t67 + t61) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + pkin(5) * t5 + pkin(9) * t22 + t97 * t4 + t64)) -m(3) * (g(1) * (rSges(3,1) * t35 - rSges(3,2) * t36) + g(2) * (rSges(3,1) * t113 + rSges(3,2) * t34) + g(3) * (rSges(3,1) * t55 - rSges(3,2) * t51) * t47) - m(4) * (g(1) * (t21 * rSges(4,1) - rSges(4,2) * t22 + t65) + g(2) * (rSges(4,1) * t18 + rSges(4,2) * t17 + t72) + g(3) * (rSges(4,1) * t30 + rSges(4,2) * t31 + t43)) - m(5) * (g(1) * (-t21 * rSges(5,2) + t22 * t83 + t62) + g(2) * (-rSges(5,2) * t18 - t83 * t17 + t71) + g(3) * (-rSges(5,2) * t30 - t83 * t31 + t84)) + (-g(1) * t62 - g(2) * t71 - g(3) * t84 + t110 * (rSges(7,1) * t49 + rSges(7,2) * t53 + pkin(9)) + t114 * (t70 * t50 - t97 * t54 + qJ(4))) * m(7) + (-g(1) * (t98 * t21 + t62) - g(2) * (t98 * t18 + t71) - g(3) * (t98 * t30 + t84) + t114 * (rSges(6,1) * t50 + rSges(6,2) * t54 + qJ(4))) * m(6) (-m(4) + t80) * (g(3) * t48 + (g(1) * t52 - g(2) * t56) * t47) t80 * t110, -m(6) * (g(1) * (-rSges(6,1) * t4 - rSges(6,2) * t5) + g(2) * (rSges(6,1) * t67 + rSges(6,2) * t68) + g(3) * (rSges(6,1) * t23 - rSges(6,2) * t24)) - m(7) * (g(2) * (t67 * t70 - t68 * t97) + (t23 * t70 + t24 * t97) * g(3) + (-t4 * t70 + t5 * t97) * g(1)) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * (rSges(7,1) * t109 + rSges(7,2) * t108) + g(3) * ((-t24 * t49 - t31 * t53) * rSges(7,1) + (-t24 * t53 + t31 * t49) * rSges(7,2)))];
taug  = t1(:);
