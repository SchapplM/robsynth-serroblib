% Calculate Gravitation load on the joints for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:59:59
% EndTime: 2019-03-08 20:00:01
% DurationCPUTime: 0.73s
% Computational Cost: add. (676->137), mult. (1717->212), div. (0->0), fcn. (2171->12), ass. (0->73)
t57 = sin(pkin(10));
t59 = cos(pkin(10));
t63 = sin(qJ(2));
t60 = cos(pkin(6));
t66 = cos(qJ(2));
t88 = t60 * t66;
t101 = -t57 * t63 + t59 * t88;
t100 = rSges(7,2) + pkin(9);
t56 = sin(pkin(11));
t82 = cos(pkin(11));
t72 = -t63 * t56 + t66 * t82;
t65 = cos(qJ(4));
t99 = pkin(4) * t65;
t98 = rSges(7,1) + pkin(5);
t97 = rSges(5,3) + pkin(8);
t96 = rSges(6,3) + pkin(9);
t47 = t72 * t60;
t50 = -t66 * t56 - t63 * t82;
t30 = t59 * t47 + t57 * t50;
t62 = sin(qJ(4));
t95 = t30 * t62;
t33 = -t57 * t47 + t59 * t50;
t94 = t33 * t62;
t58 = sin(pkin(6));
t45 = t72 * t58;
t93 = t45 * t62;
t91 = t58 * t62;
t90 = t58 * t65;
t89 = t60 * t63;
t61 = sin(qJ(5));
t87 = t61 * t65;
t64 = cos(qJ(5));
t85 = t64 * t65;
t54 = t58 * t66 * pkin(2);
t84 = t45 * pkin(3) + t54;
t83 = rSges(7,3) + qJ(6);
t80 = -m(4) - m(5) - m(6) - m(7);
t78 = t101 * pkin(2);
t77 = rSges(5,1) * t65 - rSges(5,2) * t62;
t76 = rSges(6,1) * t64 - rSges(6,2) * t61;
t48 = t50 * t60;
t31 = -t59 * t48 + t57 * t72;
t32 = -t57 * t48 - t59 * t72;
t75 = t30 * pkin(3) + t78;
t74 = -t57 * t88 - t59 * t63;
t46 = t50 * t58;
t73 = -t46 * pkin(8) + pkin(9) * t93 + t45 * t99 + t84;
t71 = t74 * pkin(2);
t69 = t33 * pkin(3) + t71;
t68 = pkin(8) * t31 + pkin(9) * t95 + t30 * t99 + t75;
t67 = -t32 * pkin(8) + pkin(9) * t94 + t33 * t99 + t69;
t37 = -t46 * t65 + t60 * t62;
t36 = t46 * t62 + t60 * t65;
t35 = t36 * pkin(4);
t18 = t45 * t85 - t46 * t61;
t17 = t45 * t87 + t46 * t64;
t16 = -t32 * t65 + t57 * t91;
t15 = t32 * t62 + t57 * t90;
t14 = t31 * t65 - t59 * t91;
t13 = -t31 * t62 - t59 * t90;
t12 = t15 * pkin(4);
t11 = t13 * pkin(4);
t10 = t37 * t64 - t45 * t61;
t9 = t37 * t61 + t45 * t64;
t8 = -t32 * t61 + t33 * t85;
t7 = t32 * t64 + t33 * t87;
t6 = t30 * t85 + t31 * t61;
t5 = t30 * t87 - t31 * t64;
t4 = t16 * t64 - t33 * t61;
t3 = t16 * t61 + t33 * t64;
t2 = t14 * t64 - t30 * t61;
t1 = t14 * t61 + t30 * t64;
t19 = [(-m(2) - m(3) + t80) * g(3), -m(3) * (g(1) * (t74 * rSges(3,1) + (t57 * t89 - t59 * t66) * rSges(3,2)) + g(2) * (t101 * rSges(3,1) + (-t57 * t66 - t59 * t89) * rSges(3,2)) + g(3) * (rSges(3,1) * t66 - rSges(3,2) * t63) * t58) - m(4) * (g(1) * (t33 * rSges(4,1) + t32 * rSges(4,2) + t71) + g(2) * (t30 * rSges(4,1) - rSges(4,2) * t31 + t78) + g(3) * (t45 * rSges(4,1) + t46 * rSges(4,2) + t54)) - m(5) * (g(1) * (-t97 * t32 + t77 * t33 + t69) + g(2) * (t77 * t30 + t31 * t97 + t75) + g(3) * (t77 * t45 - t97 * t46 + t84)) - m(6) * (g(1) * (t8 * rSges(6,1) - t7 * rSges(6,2) + rSges(6,3) * t94 + t67) + g(2) * (t6 * rSges(6,1) - t5 * rSges(6,2) + rSges(6,3) * t95 + t68) + g(3) * (t18 * rSges(6,1) - t17 * rSges(6,2) + rSges(6,3) * t93 + t73)) - m(7) * (g(1) * (rSges(7,2) * t94 + t83 * t7 + t98 * t8 + t67) + g(2) * (rSges(7,2) * t95 + t83 * t5 + t98 * t6 + t68) + g(3) * (rSges(7,2) * t93 + t83 * t17 + t98 * t18 + t73)) t80 * (g(3) * t60 + (g(1) * t57 - g(2) * t59) * t58) -m(5) * (g(1) * (t15 * rSges(5,1) - t16 * rSges(5,2)) + g(2) * (t13 * rSges(5,1) - t14 * rSges(5,2)) + g(3) * (t36 * rSges(5,1) - t37 * rSges(5,2))) - m(6) * (g(1) * (t76 * t15 + t96 * t16 + t12) + g(2) * (t76 * t13 + t96 * t14 + t11) + g(3) * (t76 * t36 + t96 * t37 + t35)) + (-g(1) * (t100 * t16 + t12) - g(2) * (t100 * t14 + t11) - g(3) * (t100 * t37 + t35) - (g(1) * t15 + g(2) * t13 + g(3) * t36) * (t83 * t61 + t98 * t64)) * m(7), -m(6) * (g(1) * (-t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) - t2 * rSges(6,2)) + g(3) * (-t9 * rSges(6,1) - t10 * rSges(6,2))) - m(7) * (g(1) * (-t98 * t3 + t83 * t4) + g(2) * (-t98 * t1 + t83 * t2) + g(3) * (t83 * t10 - t98 * t9)) -m(7) * (g(1) * t3 + g(2) * t1 + g(3) * t9)];
taug  = t19(:);
