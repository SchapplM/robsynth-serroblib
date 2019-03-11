% Calculate Gravitation load on the joints for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:07:55
% EndTime: 2019-03-09 10:07:57
% DurationCPUTime: 0.82s
% Computational Cost: add. (511->127), mult. (427->166), div. (0->0), fcn. (379->12), ass. (0->69)
t43 = qJ(2) + pkin(10);
t39 = qJ(4) + t43;
t31 = sin(t39);
t32 = cos(t39);
t45 = cos(pkin(11));
t69 = -rSges(6,1) * t45 - pkin(4);
t44 = sin(pkin(11));
t87 = rSges(6,2) * t44;
t55 = (rSges(6,3) + qJ(5)) * t31 + (-t69 - t87) * t32;
t105 = qJ(5) * t32 + t31 * t87;
t33 = t45 * pkin(5) + pkin(4);
t103 = t31 * rSges(7,3) + t32 * t33;
t101 = t32 * rSges(5,1) - t31 * rSges(5,2);
t36 = sin(t43);
t38 = cos(t43);
t50 = cos(qJ(2));
t40 = t50 * pkin(2);
t100 = t38 * rSges(4,1) - t36 * rSges(4,2) + t40;
t49 = sin(qJ(1));
t51 = cos(qJ(1));
t99 = g(1) * t51 + g(2) * t49;
t98 = t99 * t31;
t75 = pkin(3) * t38 + t40;
t15 = pkin(1) + t75;
t6 = t51 * t15;
t97 = g(2) * t6;
t96 = -m(6) - m(7);
t48 = sin(qJ(2));
t93 = t48 * pkin(2);
t92 = rSges(3,3) + pkin(7);
t42 = pkin(11) + qJ(6);
t35 = sin(t42);
t86 = rSges(7,2) * t35;
t72 = t31 * t86;
t83 = t32 * t49;
t91 = rSges(7,3) * t83 + t49 * t72;
t82 = t32 * t51;
t90 = rSges(7,3) * t82 + t51 * t72;
t37 = cos(t42);
t88 = rSges(7,1) * t37;
t47 = -pkin(9) - qJ(5);
t84 = t31 * t47;
t81 = t49 * t35;
t80 = t49 * t37;
t79 = t51 * t35;
t78 = t51 * t37;
t46 = -qJ(3) - pkin(7);
t77 = rSges(4,3) - t46;
t41 = -pkin(8) + t46;
t76 = rSges(5,3) - t41;
t71 = rSges(6,3) * t83 + t105 * t49;
t70 = rSges(6,3) * t82 + t105 * t51;
t68 = pkin(5) * t44 - t41;
t67 = -t33 - t88;
t65 = t50 * rSges(3,1) - t48 * rSges(3,2);
t61 = -rSges(5,1) * t31 - rSges(5,2) * t32;
t60 = pkin(1) + t65;
t58 = pkin(1) + t100;
t57 = t44 * rSges(6,1) + t45 * rSges(6,2) - t41;
t56 = t103 + (-t86 + t88) * t32;
t53 = -t84 + t103;
t16 = -pkin(3) * t36 - t93;
t10 = t51 * t16;
t9 = t49 * t16;
t5 = t32 * t78 + t81;
t4 = -t32 * t79 + t80;
t3 = -t32 * t80 + t79;
t2 = t32 * t81 + t78;
t1 = [-m(2) * (g(1) * (-t49 * rSges(2,1) - t51 * rSges(2,2)) + g(2) * (t51 * rSges(2,1) - t49 * rSges(2,2))) - m(3) * ((g(1) * t92 + g(2) * t60) * t51 + (-g(1) * t60 + g(2) * t92) * t49) - m(4) * ((g(1) * t77 + g(2) * t58) * t51 + (-g(1) * t58 + g(2) * t77) * t49) - m(5) * (t97 + (g(1) * t76 + g(2) * t101) * t51 + (g(1) * (-t15 - t101) + g(2) * t76) * t49) - m(6) * (t97 + (g(1) * t57 + t55 * g(2)) * t51 + (g(2) * t57 + (-t15 - t55) * g(1)) * t49) - m(7) * (g(1) * (t3 * rSges(7,1) + t2 * rSges(7,2)) + g(2) * (t5 * rSges(7,1) + t4 * rSges(7,2) + t6) + (g(1) * t68 + g(2) * t53) * t51 + (g(1) * (-t15 - t53) + g(2) * t68) * t49) -m(3) * (g(3) * t65 + t99 * (-rSges(3,1) * t48 - rSges(3,2) * t50)) - m(4) * (g(3) * t100 + t99 * (-rSges(4,1) * t36 - rSges(4,2) * t38 - t93)) - m(5) * (g(1) * (t61 * t51 + t10) + g(2) * (t61 * t49 + t9) + g(3) * (t101 + t75)) - m(6) * (g(1) * (t10 + t70) + g(2) * (t9 + t71) + g(3) * (t55 + t75) + t69 * t98) - m(7) * (g(1) * (-t47 * t82 + t10 + t90) + g(2) * (-t47 * t83 + t9 + t91) + g(3) * (t56 + t75) + (-g(3) * t47 + t99 * t67) * t31) (-m(4) - m(5) + t96) * (g(1) * t49 - g(2) * t51) -m(5) * g(3) * t101 - m(6) * (g(1) * t70 + g(2) * t71 + g(3) * t55) - m(7) * (g(1) * t90 + g(2) * t91 + g(3) * (t56 - t84)) + t99 * ((m(5) * rSges(5,2) + m(7) * t47) * t32 + (m(5) * rSges(5,1) - m(6) * t69 - m(7) * t67) * t31) t96 * (-g(3) * t32 + t98) -m(7) * (g(1) * (t4 * rSges(7,1) - t5 * rSges(7,2)) + g(2) * (-t2 * rSges(7,1) + t3 * rSges(7,2)) + g(3) * (-rSges(7,1) * t35 - rSges(7,2) * t37) * t31)];
taug  = t1(:);
