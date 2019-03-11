% Calculate Gravitation load on the joints for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:12:23
% EndTime: 2019-03-09 23:12:28
% DurationCPUTime: 1.56s
% Computational Cost: add. (877->193), mult. (1662->272), div. (0->0), fcn. (1987->14), ass. (0->81)
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t85 = cos(pkin(6));
t95 = cos(qJ(1));
t79 = t85 * t95;
t94 = sin(qJ(1));
t27 = t57 * t79 + t94 * t60;
t56 = sin(qJ(3));
t59 = cos(qJ(3));
t53 = sin(pkin(6));
t82 = t53 * t95;
t15 = t27 * t59 - t56 * t82;
t26 = t94 * t57 - t60 * t79;
t55 = sin(qJ(4));
t58 = cos(qJ(4));
t118 = t15 * t55 - t26 * t58;
t117 = -t15 * t58 - t26 * t55;
t78 = t85 * t94;
t29 = -t57 * t78 + t95 * t60;
t116 = g(1) * t29 + g(2) * t27;
t28 = t95 * t57 + t60 * t78;
t115 = g(1) * t28 + g(2) * t26;
t81 = t53 * t94;
t19 = t29 * t59 + t56 * t81;
t91 = t53 * t57;
t25 = t85 * t56 + t59 * t91;
t114 = g(1) * t19 + g(2) * t15 + g(3) * t25;
t14 = t27 * t56 + t59 * t82;
t18 = t29 * t56 - t59 * t81;
t24 = t56 * t91 - t85 * t59;
t113 = g(1) * t18 + g(2) * t14 + g(3) * t24;
t52 = qJ(4) + pkin(12);
t47 = cos(t52);
t49 = t58 * pkin(4);
t32 = pkin(5) * t47 + t49;
t30 = pkin(3) + t32;
t48 = qJ(6) + t52;
t43 = sin(t48);
t44 = cos(t48);
t69 = rSges(7,1) * t44 - rSges(7,2) * t43 + t30;
t54 = -qJ(5) - pkin(10);
t87 = pkin(11) - t54 + rSges(7,3);
t112 = t87 * t56 + t69 * t59;
t45 = t49 + pkin(3);
t46 = sin(t52);
t70 = rSges(6,1) * t47 - rSges(6,2) * t46 + t45;
t86 = -t54 + rSges(6,3);
t111 = t86 * t56 + t70 * t59;
t73 = rSges(5,1) * t58 - rSges(5,2) * t55 + pkin(3);
t96 = pkin(10) + rSges(5,3);
t110 = t96 * t56 + t73 * t59;
t109 = (-t15 * t43 + t26 * t44) * rSges(7,1) + (-t15 * t44 - t26 * t43) * rSges(7,2);
t6 = -t19 * t43 + t28 * t44;
t7 = t19 * t44 + t28 * t43;
t108 = t6 * rSges(7,1) - t7 * rSges(7,2);
t107 = pkin(4) * t55;
t98 = g(3) * t53;
t97 = rSges(4,3) + pkin(9);
t90 = t53 * t60;
t89 = (-t25 * t43 - t44 * t90) * rSges(7,1) + (-t25 * t44 + t43 * t90) * rSges(7,2);
t88 = t95 * pkin(1) + pkin(8) * t81;
t84 = t29 * pkin(2) + t88;
t83 = g(3) * (pkin(2) * t90 + pkin(9) * t91);
t80 = -t94 * pkin(1) + pkin(8) * t82;
t77 = rSges(4,1) * t59 - rSges(4,2) * t56;
t76 = rSges(5,1) * t55 + rSges(5,2) * t58;
t10 = -t19 * t55 + t28 * t58;
t74 = -t27 * pkin(2) + t80;
t71 = -t25 * t55 - t58 * t90;
t31 = pkin(5) * t46 + t107;
t68 = rSges(7,1) * t43 + rSges(7,2) * t44 + t31;
t67 = pkin(9) + t68;
t66 = rSges(6,1) * t46 + rSges(6,2) * t47 + t107;
t65 = pkin(9) + t66;
t20 = t26 * pkin(2);
t22 = t28 * pkin(2);
t64 = -g(1) * t22 - g(2) * t20 + t83;
t11 = t19 * t58 + t28 * t55;
t9 = t19 * t47 + t28 * t46;
t8 = -t19 * t46 + t28 * t47;
t1 = [-m(2) * (g(1) * (-t94 * rSges(2,1) - t95 * rSges(2,2)) + g(2) * (t95 * rSges(2,1) - t94 * rSges(2,2))) - m(3) * (g(1) * (-t27 * rSges(3,1) + t26 * rSges(3,2) + rSges(3,3) * t82 + t80) + g(2) * (t29 * rSges(3,1) - t28 * rSges(3,2) + rSges(3,3) * t81 + t88)) - m(4) * (g(1) * (-rSges(4,1) * t15 + rSges(4,2) * t14 - t97 * t26 + t74) + g(2) * (rSges(4,1) * t19 - rSges(4,2) * t18 + t97 * t28 + t84)) - m(5) * (g(1) * (t117 * rSges(5,1) + t118 * rSges(5,2) - t15 * pkin(3) - t26 * pkin(9) - t96 * t14 + t74) + g(2) * (rSges(5,1) * t11 + rSges(5,2) * t10 + pkin(3) * t19 + pkin(9) * t28 + t96 * t18 + t84)) - m(6) * (g(1) * (-t14 * t86 - t15 * t70 - t65 * t26 + t74) + g(2) * (rSges(6,1) * t9 + rSges(6,2) * t8 + t19 * t45 + (pkin(9) + t107) * t28 + t86 * t18 + t84)) - m(7) * (g(1) * (-t14 * t87 - t15 * t69 - t67 * t26 + t74) + g(2) * (rSges(7,1) * t7 + rSges(7,2) * t6 + t19 * t30 + (pkin(9) + t31) * t28 + t87 * t18 + t84)) -m(3) * (g(1) * (-rSges(3,1) * t28 - rSges(3,2) * t29) + g(2) * (-rSges(3,1) * t26 - rSges(3,2) * t27) + (rSges(3,1) * t60 - rSges(3,2) * t57) * t98) - m(4) * (g(1) * (-t77 * t28 + t97 * t29 - t22) + g(2) * (-t77 * t26 + t97 * t27 - t20) + t83 + (rSges(4,3) * t57 + t77 * t60) * t98) - m(5) * ((t110 * t60 + t76 * t57) * t98 + t64 + t116 * (pkin(9) + t76) - t115 * t110) - m(6) * ((t111 * t60 + t66 * t57) * t98 + t64 + t116 * t65 - t115 * t111) - m(7) * ((t112 * t60 + t68 * t57) * t98 + t64 + t116 * t67 - t115 * t112) -m(4) * (g(1) * (-rSges(4,1) * t18 - rSges(4,2) * t19) + g(2) * (-rSges(4,1) * t14 - rSges(4,2) * t15) + g(3) * (-rSges(4,1) * t24 - rSges(4,2) * t25)) - m(5) * (-t113 * t73 + t114 * t96) - m(6) * (-t113 * t70 + t114 * t86) - m(7) * (-t113 * t69 + t114 * t87) -m(5) * (g(1) * (rSges(5,1) * t10 - rSges(5,2) * t11) + g(2) * (-rSges(5,1) * t118 + t117 * rSges(5,2)) + g(3) * (t71 * rSges(5,1) + (-t25 * t58 + t55 * t90) * rSges(5,2))) - m(7) * (g(1) * (-t19 * t31 + t28 * t32 + t108) + g(2) * (-t15 * t31 + t26 * t32 + t109) + g(3) * (-t25 * t31 - t32 * t90 + t89)) + (-g(1) * (rSges(6,1) * t8 - rSges(6,2) * t9) - g(2) * ((-t15 * t46 + t26 * t47) * rSges(6,1) + (-t15 * t47 - t26 * t46) * rSges(6,2)) - g(3) * ((-t25 * t46 - t47 * t90) * rSges(6,1) + (-t25 * t47 + t46 * t90) * rSges(6,2)) - (g(1) * t10 - g(2) * t118 + g(3) * t71) * pkin(4)) * m(6) (-m(6) - m(7)) * t113, -m(7) * (g(1) * t108 + g(2) * t109 + g(3) * t89)];
taug  = t1(:);
