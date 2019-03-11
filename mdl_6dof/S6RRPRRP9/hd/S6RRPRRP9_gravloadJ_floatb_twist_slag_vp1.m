% Calculate Gravitation load on the joints for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:27:32
% EndTime: 2019-03-09 12:27:35
% DurationCPUTime: 1.28s
% Computational Cost: add. (761->185), mult. (1295->268), div. (0->0), fcn. (1524->12), ass. (0->75)
t58 = sin(qJ(2));
t59 = sin(qJ(1));
t61 = cos(qJ(2));
t76 = cos(pkin(6));
t95 = cos(qJ(1));
t68 = t76 * t95;
t33 = t58 * t68 + t59 * t61;
t51 = pkin(11) + qJ(4);
t48 = sin(t51);
t49 = cos(t51);
t53 = sin(pkin(6));
t73 = t53 * t95;
t15 = t33 * t49 - t48 * t73;
t32 = t58 * t59 - t61 * t68;
t57 = sin(qJ(5));
t60 = cos(qJ(5));
t3 = t15 * t57 - t32 * t60;
t94 = t32 * t57;
t110 = -t15 * t60 - t94;
t71 = t59 * t76;
t34 = t95 * t58 + t61 * t71;
t98 = g(2) * t32;
t109 = g(1) * t34 + t98;
t84 = t53 * t61;
t108 = -g(3) * t84 + t109;
t47 = pkin(5) * t60 + pkin(4);
t81 = rSges(7,3) + qJ(6) + pkin(10);
t107 = t47 * t49 + t81 * t48;
t35 = -t58 * t71 + t95 * t61;
t85 = t53 * t59;
t19 = t35 * t49 + t48 * t85;
t86 = t53 * t58;
t27 = t76 * t48 + t49 * t86;
t106 = g(1) * t19 + g(2) * t15 + g(3) * t27;
t14 = t33 * t48 + t49 * t73;
t18 = t35 * t48 - t49 * t85;
t26 = t48 * t86 - t76 * t49;
t105 = g(1) * t18 + g(2) * t14 + g(3) * t26;
t104 = pkin(4) * t49;
t97 = g(3) * t53;
t96 = rSges(6,3) + pkin(10);
t92 = t33 * t57;
t91 = t34 * t57;
t90 = t35 * t57;
t88 = t49 * t57;
t87 = t49 * t60;
t83 = t57 * t61;
t82 = t60 * t61;
t54 = cos(pkin(11));
t46 = pkin(3) * t54 + pkin(2);
t56 = -pkin(9) - qJ(3);
t80 = -t32 * t46 - t33 * t56;
t79 = -t34 * t46 - t35 * t56;
t78 = t95 * pkin(1) + pkin(8) * t85;
t77 = qJ(3) + rSges(4,3);
t52 = sin(pkin(11));
t74 = t52 * t85;
t72 = -t59 * pkin(1) + pkin(8) * t73;
t70 = t52 * t73;
t69 = pkin(3) * t74 - t34 * t56 + t35 * t46 + t78;
t67 = rSges(5,1) * t49 - rSges(5,2) * t48;
t5 = -t19 * t57 + t34 * t60;
t66 = rSges(4,1) * t54 - rSges(4,2) * t52 + pkin(2);
t12 = -t27 * t57 - t53 * t82;
t63 = pkin(3) * t70 + t32 * t56 - t33 * t46 + t72;
t36 = t46 * t84;
t21 = (t49 * t82 + t57 * t58) * t53;
t20 = (-t49 * t83 + t58 * t60) * t53;
t13 = -t27 * t60 + t53 * t83;
t10 = -t34 * t87 + t90;
t9 = t34 * t88 + t35 * t60;
t8 = -t32 * t87 + t92;
t7 = t32 * t88 + t33 * t60;
t6 = t19 * t60 + t91;
t1 = [-m(2) * (g(1) * (-t59 * rSges(2,1) - t95 * rSges(2,2)) + g(2) * (t95 * rSges(2,1) - t59 * rSges(2,2))) - m(3) * (g(1) * (-t33 * rSges(3,1) + t32 * rSges(3,2) + rSges(3,3) * t73 + t72) + g(2) * (rSges(3,1) * t35 - rSges(3,2) * t34 + rSges(3,3) * t85 + t78)) - m(4) * (g(1) * (-t33 * pkin(2) + (-t33 * t54 + t70) * rSges(4,1) + (t33 * t52 + t54 * t73) * rSges(4,2) - t77 * t32 + t72) + g(2) * (t35 * pkin(2) + (t35 * t54 + t74) * rSges(4,1) + (-t35 * t52 + t54 * t85) * rSges(4,2) + t77 * t34 + t78)) - m(5) * (g(1) * (-rSges(5,1) * t15 + rSges(5,2) * t14 - rSges(5,3) * t32 + t63) + g(2) * (rSges(5,1) * t19 - rSges(5,2) * t18 + rSges(5,3) * t34 + t69)) - m(6) * (g(1) * (rSges(6,1) * t110 + rSges(6,2) * t3 - pkin(4) * t15 - t14 * t96 + t63) + g(2) * (rSges(6,1) * t6 + rSges(6,2) * t5 + pkin(4) * t19 + t96 * t18 + t69)) - m(7) * (g(1) * (rSges(7,1) * t110 + rSges(7,2) * t3 - pkin(5) * t94 - t14 * t81 - t15 * t47 + t63) + g(2) * (rSges(7,1) * t6 + rSges(7,2) * t5 + pkin(5) * t91 + t81 * t18 + t19 * t47 + t69)) -m(3) * (g(1) * (-rSges(3,1) * t34 - rSges(3,2) * t35) + g(2) * (-rSges(3,1) * t32 - rSges(3,2) * t33) + (rSges(3,1) * t61 - rSges(3,2) * t58) * t97) - m(4) * (g(1) * (-t66 * t34 + t77 * t35) + g(2) * t77 * t33 - t66 * t98 + (t77 * t58 + t66 * t61) * t97) - m(5) * (g(1) * (rSges(5,3) * t35 - t67 * t34 + t79) + g(2) * (rSges(5,3) * t33 - t67 * t32 + t80) + g(3) * t36 + (t67 * t61 + (rSges(5,3) - t56) * t58) * t97) - m(6) * (g(1) * (rSges(6,1) * t10 + rSges(6,2) * t9 - t34 * t104 + t79) + g(2) * (rSges(6,1) * t8 + rSges(6,2) * t7 - t32 * t104 + t80) + g(3) * (t21 * rSges(6,1) + t20 * rSges(6,2) + t84 * t104 - t56 * t86 + t36) - t108 * t48 * t96) - m(7) * (g(1) * (rSges(7,1) * t10 + rSges(7,2) * t9 + pkin(5) * t90 + t79) + g(2) * (rSges(7,1) * t8 + rSges(7,2) * t7 + pkin(5) * t92 + t80) + g(3) * (t21 * rSges(7,1) + t20 * rSges(7,2) + t36) + ((pkin(5) * t57 - t56) * t58 + t107 * t61) * t97 - t109 * t107) (-m(4) - m(5) - m(6) - m(7)) * t108, -m(5) * (g(1) * (-rSges(5,1) * t18 - rSges(5,2) * t19) + g(2) * (-rSges(5,1) * t14 - rSges(5,2) * t15) + g(3) * (-rSges(5,1) * t26 - rSges(5,2) * t27)) - m(6) * (t106 * t96 + t105 * (-rSges(6,1) * t60 + rSges(6,2) * t57 - pkin(4))) - m(7) * (t106 * t81 + t105 * (-rSges(7,1) * t60 + rSges(7,2) * t57 - t47)) -m(6) * (g(1) * (rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (-rSges(6,1) * t3 + rSges(6,2) * t110) + g(3) * (rSges(6,1) * t12 + rSges(6,2) * t13)) + (-g(1) * (rSges(7,1) * t5 - rSges(7,2) * t6) - g(2) * (-rSges(7,1) * t3 + rSges(7,2) * t110) - g(3) * (t12 * rSges(7,1) + t13 * rSges(7,2)) - (g(1) * t5 - g(2) * t3 + g(3) * t12) * pkin(5)) * m(7), -m(7) * t105];
taug  = t1(:);
