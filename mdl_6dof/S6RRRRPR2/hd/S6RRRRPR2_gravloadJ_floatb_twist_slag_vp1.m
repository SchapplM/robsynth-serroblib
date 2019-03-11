% Calculate Gravitation load on the joints for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:57:24
% EndTime: 2019-03-09 21:57:27
% DurationCPUTime: 0.97s
% Computational Cost: add. (653->142), mult. (531->180), div. (0->0), fcn. (471->12), ass. (0->77)
t43 = qJ(2) + qJ(3);
t39 = qJ(4) + t43;
t32 = sin(t39);
t33 = cos(t39);
t45 = cos(pkin(11));
t77 = -rSges(6,1) * t45 - pkin(4);
t44 = sin(pkin(11));
t94 = rSges(6,2) * t44;
t61 = (rSges(6,3) + qJ(5)) * t32 + (-t77 - t94) * t33;
t49 = cos(qJ(2));
t40 = t49 * pkin(2);
t37 = sin(t43);
t38 = cos(t43);
t74 = t38 * rSges(4,1) - rSges(4,2) * t37;
t114 = t40 + t74;
t113 = qJ(5) * t33 + t32 * t94;
t30 = pkin(5) * t45 + pkin(4);
t111 = t32 * rSges(7,3) + t33 * t30;
t109 = t33 * rSges(5,1) - rSges(5,2) * t32;
t48 = sin(qJ(1));
t50 = cos(qJ(1));
t108 = g(1) * t50 + g(2) * t48;
t107 = t108 * t32;
t31 = pkin(3) * t38;
t83 = t31 + t40;
t15 = pkin(1) + t83;
t6 = t50 * t15;
t106 = g(2) * t6;
t51 = -pkin(8) - pkin(7);
t47 = sin(qJ(2));
t105 = pkin(2) * t47;
t104 = pkin(3) * t37;
t101 = rSges(3,3) + pkin(7);
t41 = pkin(11) + qJ(6);
t35 = sin(t41);
t93 = rSges(7,2) * t35;
t80 = t32 * t93;
t91 = t33 * t48;
t100 = rSges(7,3) * t91 + t48 * t80;
t90 = t33 * t50;
t99 = rSges(7,3) * t90 + t50 * t80;
t36 = cos(t41);
t97 = rSges(7,1) * t36;
t46 = -pkin(10) - qJ(5);
t92 = t32 * t46;
t89 = t35 * t48;
t88 = t35 * t50;
t87 = t36 * t48;
t86 = t36 * t50;
t85 = rSges(4,3) - t51;
t42 = -pkin(9) + t51;
t84 = rSges(5,3) - t42;
t79 = rSges(6,3) * t91 + t113 * t48;
t78 = rSges(6,3) * t90 + t113 * t50;
t76 = pkin(5) * t44 - t42;
t75 = -t30 - t97;
t72 = t31 + t109;
t71 = rSges(3,1) * t49 - rSges(3,2) * t47;
t69 = -rSges(4,1) * t37 - rSges(4,2) * t38;
t67 = -rSges(5,1) * t32 - rSges(5,2) * t33;
t66 = pkin(1) + t71;
t64 = pkin(1) + t114;
t63 = t44 * rSges(6,1) + t45 * rSges(6,2) - t42;
t62 = t111 + (-t93 + t97) * t33;
t58 = g(1) * t99 + g(2) * t100;
t57 = -t92 + t111;
t56 = t31 + t61;
t55 = t62 - t92;
t53 = t77 * t107;
t16 = -t104 - t105;
t10 = t50 * t16;
t9 = t48 * t16;
t5 = t33 * t86 + t89;
t4 = -t33 * t88 + t87;
t3 = -t33 * t87 + t88;
t2 = t33 * t89 + t86;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t48 - rSges(2,2) * t50) + g(2) * (rSges(2,1) * t50 - rSges(2,2) * t48)) - m(3) * ((g(1) * t101 + g(2) * t66) * t50 + (-g(1) * t66 + g(2) * t101) * t48) - m(4) * ((g(1) * t85 + g(2) * t64) * t50 + (-g(1) * t64 + g(2) * t85) * t48) - m(5) * (t106 + (g(1) * t84 + g(2) * t109) * t50 + (g(1) * (-t15 - t109) + g(2) * t84) * t48) - m(6) * (t106 + (g(1) * t63 + t61 * g(2)) * t50 + (g(2) * t63 + (-t15 - t61) * g(1)) * t48) - m(7) * (g(1) * (rSges(7,1) * t3 + rSges(7,2) * t2) + g(2) * (rSges(7,1) * t5 + rSges(7,2) * t4 + t6) + (g(1) * t76 + g(2) * t57) * t50 + (g(1) * (-t15 - t57) + g(2) * t76) * t48) -m(3) * (g(3) * t71 + t108 * (-rSges(3,1) * t47 - rSges(3,2) * t49)) - m(4) * (g(3) * t114 + t108 * (t69 - t105)) - m(5) * (g(1) * (t67 * t50 + t10) + g(2) * (t67 * t48 + t9) + g(3) * (t40 + t72)) - m(6) * (g(1) * (t10 + t78) + g(2) * (t9 + t79) + g(3) * (t40 + t56) + t53) - m(7) * (g(1) * (-t46 * t90 + t10 + t99) + g(2) * (-t46 * t91 + t100 + t9) + g(3) * (t62 + t83) + (-g(3) * t46 + t108 * t75) * t32) -m(4) * (g(3) * t74 + t108 * t69) - m(5) * (g(3) * t72 + t108 * (t67 - t104)) - m(6) * (g(1) * (-t50 * t104 + t78) + g(2) * (-t48 * t104 + t79) + g(3) * t56 + t53) - m(7) * (g(3) * (t31 + t55) + t58 + t108 * (t75 * t32 - t33 * t46 - t104)) -m(5) * g(3) * t109 - m(6) * (g(1) * t78 + g(2) * t79 + g(3) * t61) - m(7) * (g(3) * t55 + t58) + t108 * ((m(5) * rSges(5,2) + m(7) * t46) * t33 + (m(5) * rSges(5,1) - m(6) * t77 - m(7) * t75) * t32) (-m(6) - m(7)) * (-g(3) * t33 + t107) -m(7) * (g(1) * (rSges(7,1) * t4 - rSges(7,2) * t5) + g(2) * (-rSges(7,1) * t2 + rSges(7,2) * t3) + g(3) * (-rSges(7,1) * t35 - rSges(7,2) * t36) * t32)];
taug  = t1(:);
