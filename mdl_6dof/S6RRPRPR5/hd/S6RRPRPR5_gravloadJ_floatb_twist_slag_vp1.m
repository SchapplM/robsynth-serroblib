% Calculate Gravitation load on the joints for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:29:01
% EndTime: 2019-03-09 10:29:05
% DurationCPUTime: 1.63s
% Computational Cost: add. (871->199), mult. (2005->298), div. (0->0), fcn. (2515->14), ass. (0->87)
t108 = cos(qJ(2));
t54 = sin(qJ(2));
t87 = sin(pkin(11));
t89 = cos(pkin(11));
t34 = -t108 * t87 - t54 * t89;
t55 = sin(qJ(1));
t57 = cos(qJ(1));
t90 = cos(pkin(6));
t70 = t90 * t87;
t71 = t90 * t89;
t58 = t108 * t71 - t54 * t70;
t18 = t34 * t57 - t55 * t58;
t76 = t90 * t108;
t31 = -t57 * t54 - t55 * t76;
t60 = t31 * pkin(2);
t121 = t18 * pkin(3) + t60;
t120 = -t55 * t54 + t57 * t76;
t88 = sin(pkin(6));
t68 = t88 * t87;
t69 = t89 * t88;
t26 = -t108 * t69 + t54 * t68;
t119 = g(1) * t18 - g(3) * t26;
t27 = t108 * t68 + t54 * t69;
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t21 = t27 * t56 + t53 * t90;
t33 = -t108 * t89 + t54 * t87;
t94 = -t108 * t70 - t54 * t71;
t19 = -t57 * t33 + t55 * t94;
t81 = t55 * t88;
t9 = t19 * t56 + t53 * t81;
t118 = g(1) * t9 + g(3) * t21;
t14 = t55 * t33 + t57 * t94;
t20 = t27 * t53 - t56 * t90;
t80 = t57 * t88;
t4 = -t14 * t53 + t56 * t80;
t8 = t19 * t53 - t56 * t81;
t117 = g(1) * t8 + g(2) * t4 + g(3) * t20;
t115 = -m(6) - m(7);
t110 = t56 * pkin(4);
t109 = rSges(5,3) + pkin(9);
t50 = sin(pkin(12));
t107 = t14 * t50;
t106 = t19 * t50;
t105 = t27 * t50;
t49 = pkin(12) + qJ(6);
t47 = sin(t49);
t104 = t47 * t56;
t48 = cos(t49);
t103 = t48 * t56;
t102 = t50 * t56;
t51 = cos(pkin(12));
t101 = t51 * t56;
t46 = pkin(2) * t108 + pkin(1);
t99 = t55 * t46;
t45 = pkin(5) * t51 + pkin(4);
t97 = t56 * t45;
t75 = t108 * t88;
t42 = pkin(2) * t75;
t93 = -t26 * pkin(3) + t42;
t92 = pkin(10) + qJ(5) + rSges(7,3);
t91 = qJ(5) + rSges(6,3);
t86 = -t50 * pkin(5) - pkin(9);
t85 = g(2) * t92;
t84 = g(2) * t91;
t83 = t88 * pkin(8);
t5 = -t14 * t56 - t53 * t80;
t82 = t54 * t90;
t79 = t120 * pkin(2);
t78 = t27 * pkin(9) + t93;
t28 = pkin(2) * t82 - qJ(3) * t88 - t83;
t39 = t57 * t46;
t77 = t19 * pkin(3) - t55 * t28 + t39;
t74 = rSges(5,1) * t56 - rSges(5,2) * t53;
t73 = rSges(4,3) * t88 - t28;
t15 = t55 * t34 + t57 * t58;
t72 = t15 * pkin(3) + t79;
t64 = t48 * rSges(7,1) - t47 * rSges(7,2) + t45;
t63 = pkin(3) * t14 - t57 * t28 - t99;
t62 = -t14 * pkin(9) + t72;
t61 = rSges(3,3) * t88 + t83;
t59 = pkin(9) * t19 + t121;
t32 = t108 * t57 - t55 * t82;
t30 = -t108 * t55 - t57 * t82;
t3 = -t18 * t47 + t48 * t9;
t2 = -t18 * t48 - t47 * t9;
t1 = [-m(2) * (g(1) * (-t55 * rSges(2,1) - rSges(2,2) * t57) + g(2) * (rSges(2,1) * t57 - t55 * rSges(2,2))) - m(3) * (g(1) * (t30 * rSges(3,1) - rSges(3,2) * t120 - t55 * pkin(1) + t57 * t61) + g(2) * (t32 * rSges(3,1) + t31 * rSges(3,2) + t57 * pkin(1) + t55 * t61)) - m(4) * (g(1) * (rSges(4,1) * t14 - t15 * rSges(4,2) + t57 * t73 - t99) + g(2) * (t19 * rSges(4,1) + t18 * rSges(4,2) + t55 * t73 + t39)) - m(5) * (g(1) * (-rSges(5,1) * t5 + rSges(5,2) * t4 + t109 * t15 + t63) + g(2) * (rSges(5,1) * t9 - rSges(5,2) * t8 - t109 * t18 + t77)) - m(6) * (g(1) * (-t5 * pkin(4) + t15 * pkin(9) + (t15 * t50 - t5 * t51) * rSges(6,1) + (t15 * t51 + t5 * t50) * rSges(6,2) - t91 * t4 + t63) + g(2) * (-t18 * pkin(9) + t9 * pkin(4) + (-t18 * t50 + t51 * t9) * rSges(6,1) + (-t18 * t51 - t50 * t9) * rSges(6,2) + t91 * t8 + t77)) - m(7) * (g(1) * (-t92 * t4 - t64 * t5 + (t47 * rSges(7,1) + t48 * rSges(7,2) - t86) * t15 + t63) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + t18 * t86 + t45 * t9 + t8 * t92 + t77)) -m(3) * (g(1) * (rSges(3,1) * t31 - rSges(3,2) * t32) + g(2) * (rSges(3,1) * t120 + rSges(3,2) * t30) + g(3) * (-rSges(3,2) * t54 * t88 + rSges(3,1) * t75)) - m(4) * (g(1) * (t18 * rSges(4,1) - rSges(4,2) * t19 + t60) + g(2) * (rSges(4,1) * t15 + rSges(4,2) * t14 + t79) + g(3) * (-rSges(4,1) * t26 - rSges(4,2) * t27 + t42)) - m(5) * (g(1) * (t109 * t19 + t18 * t74 + t121) + g(2) * (-t109 * t14 + t15 * t74 + t72) + g(3) * (t109 * t27 - t26 * t74 + t93)) - m(6) * (g(1) * (t18 * t110 + (t101 * t18 + t106) * rSges(6,1) + (-t102 * t18 + t19 * t51) * rSges(6,2) + t59) + g(2) * (t15 * t110 + (t101 * t15 - t107) * rSges(6,1) + (-t102 * t15 - t14 * t51) * rSges(6,2) + t62) + g(3) * (-t26 * t110 + (-t101 * t26 + t105) * rSges(6,1) + (t102 * t26 + t27 * t51) * rSges(6,2) + t78) + (t119 * t91 + t15 * t84) * t53) - m(7) * (g(1) * (t18 * t97 + pkin(5) * t106 + (t103 * t18 + t19 * t47) * rSges(7,1) + (-t104 * t18 + t19 * t48) * rSges(7,2) + t59) + g(2) * (t15 * t97 - pkin(5) * t107 + (t103 * t15 - t14 * t47) * rSges(7,1) + (-t104 * t15 - t14 * t48) * rSges(7,2) + t62) + g(3) * (-t26 * t97 + pkin(5) * t105 + (-t103 * t26 + t27 * t47) * rSges(7,1) + (t104 * t26 + t27 * t48) * rSges(7,2) + t78) + (t119 * t92 + t15 * t85) * t53) (-m(4) - m(5) + t115) * (g(1) * t81 - g(2) * t80 + g(3) * t90) -m(5) * (g(1) * (-rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (-rSges(5,1) * t4 - rSges(5,2) * t5) + g(3) * (-rSges(5,1) * t20 - rSges(5,2) * t21)) - m(6) * (t5 * t84 + t118 * t91 + t117 * (-t51 * rSges(6,1) + t50 * rSges(6,2) - pkin(4))) - m(7) * (-t117 * t64 + t118 * t92 + t5 * t85) t115 * t117, -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * ((-t15 * t48 - t47 * t5) * rSges(7,1) + (t15 * t47 - t48 * t5) * rSges(7,2)) + g(3) * ((-t21 * t47 + t26 * t48) * rSges(7,1) + (-t21 * t48 - t26 * t47) * rSges(7,2)))];
taug  = t1(:);
