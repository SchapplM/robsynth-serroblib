% Calculate Gravitation load on the joints for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:39:04
% EndTime: 2019-03-09 09:39:05
% DurationCPUTime: 0.81s
% Computational Cost: add. (551->163), mult. (1020->226), div. (0->0), fcn. (1180->12), ass. (0->64)
t51 = sin(qJ(6));
t54 = cos(qJ(6));
t102 = -rSges(7,1) * t51 - rSges(7,2) * t54;
t52 = sin(qJ(2));
t53 = sin(qJ(1));
t55 = cos(qJ(2));
t56 = cos(qJ(1));
t79 = cos(pkin(6));
t72 = t56 * t79;
t29 = t52 * t72 + t53 * t55;
t73 = t53 * t79;
t31 = -t52 * t73 + t55 * t56;
t101 = g(1) * t31 + g(2) * t29;
t80 = qJ(4) + rSges(5,3);
t28 = t52 * t53 - t55 * t72;
t46 = pkin(11) + qJ(5);
t43 = sin(t46);
t44 = cos(t46);
t48 = sin(pkin(6));
t84 = t48 * t56;
t65 = -t28 * t43 + t44 * t84;
t100 = t29 * t54 + t51 * t65;
t99 = -t29 * t51 + t54 * t65;
t47 = sin(pkin(11));
t98 = pkin(4) * t47;
t95 = g(3) * t48;
t94 = rSges(7,3) + pkin(10);
t91 = t28 * t47;
t30 = t56 * t52 + t55 * t73;
t88 = t30 * t47;
t87 = t48 * t52;
t86 = t48 * t53;
t85 = t48 * t55;
t83 = pkin(2) * t85 + qJ(3) * t87;
t82 = t56 * pkin(1) + pkin(8) * t86;
t81 = rSges(4,3) + qJ(3);
t78 = -m(5) - m(6) - m(7);
t24 = t28 * pkin(2);
t50 = -pkin(9) - qJ(4);
t77 = t28 * t50 + t29 * t98 - t24;
t26 = t30 * pkin(2);
t76 = t30 * t50 + t31 * t98 - t26;
t75 = t31 * pkin(2) + t82;
t74 = -t53 * pkin(1) + pkin(8) * t84;
t71 = g(3) * (t87 * t98 + t83);
t70 = -t29 * pkin(2) + t74;
t49 = cos(pkin(11));
t69 = rSges(5,1) * t47 + rSges(5,2) * t49;
t68 = rSges(6,1) * t43 + rSges(6,2) * t44;
t67 = t30 * qJ(3) + t75;
t66 = rSges(7,1) * t54 - rSges(7,2) * t51 + pkin(5);
t64 = t28 * t44 + t43 * t84;
t61 = -t28 * qJ(3) + t70;
t42 = pkin(4) * t49 + pkin(3);
t60 = pkin(4) * t88 - t31 * t50 + t42 * t86 + t67;
t59 = -pkin(4) * t91 + t29 * t50 + t42 * t84 + t61;
t58 = t66 * t43 - t94 * t44;
t15 = -t43 * t85 + t79 * t44;
t14 = -t79 * t43 - t44 * t85;
t6 = t30 * t43 + t44 * t86;
t5 = -t30 * t44 + t43 * t86;
t2 = t31 * t51 + t54 * t6;
t1 = t31 * t54 - t51 * t6;
t3 = [-m(2) * (g(1) * (-t53 * rSges(2,1) - rSges(2,2) * t56) + g(2) * (rSges(2,1) * t56 - t53 * rSges(2,2))) - m(3) * (g(1) * (-t29 * rSges(3,1) + t28 * rSges(3,2) + rSges(3,3) * t84 + t74) + g(2) * (rSges(3,1) * t31 - rSges(3,2) * t30 + rSges(3,3) * t86 + t82)) - m(4) * (g(1) * (rSges(4,1) * t84 + t29 * rSges(4,2) - t81 * t28 + t70) + g(2) * (rSges(4,1) * t86 - rSges(4,2) * t31 + t81 * t30 + t75)) - m(5) * (g(1) * (pkin(3) * t84 + (t49 * t84 - t91) * rSges(5,1) + (-t28 * t49 - t47 * t84) * rSges(5,2) - t80 * t29 + t61) + g(2) * (pkin(3) * t86 + (t49 * t86 + t88) * rSges(5,1) + (t30 * t49 - t47 * t86) * rSges(5,2) + t80 * t31 + t67)) - m(6) * (g(1) * (rSges(6,1) * t65 - rSges(6,2) * t64 - rSges(6,3) * t29 + t59) + g(2) * (rSges(6,1) * t6 - rSges(6,2) * t5 + rSges(6,3) * t31 + t60)) - m(7) * (g(1) * (rSges(7,1) * t99 - rSges(7,2) * t100 + t65 * pkin(5) + t94 * t64 + t59) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 + pkin(5) * t6 + t94 * t5 + t60)) -m(3) * (g(1) * (-rSges(3,1) * t30 - rSges(3,2) * t31) + g(2) * (-rSges(3,1) * t28 - rSges(3,2) * t29) + (rSges(3,1) * t55 - rSges(3,2) * t52) * t95) - m(4) * (g(1) * (rSges(4,2) * t30 + t81 * t31 - t26) + g(2) * (rSges(4,2) * t28 + t81 * t29 - t24) + g(3) * ((-rSges(4,2) * t55 + rSges(4,3) * t52) * t48 + t83)) - m(5) * (g(1) * (-t30 * t80 - t26) + g(2) * (-t28 * t80 - t24) + g(3) * t83 + (t69 * t52 + t80 * t55) * t95 + t101 * (qJ(3) + t69)) - m(6) * (g(1) * (-rSges(6,3) * t30 + t76) + g(2) * (-rSges(6,3) * t28 + t77) + t71 + ((rSges(6,3) - t50) * t55 + t68 * t52) * t95 + t101 * (qJ(3) + t68)) - m(7) * (g(1) * (t102 * t30 + t76) + g(2) * (t102 * t28 + t77) + t71 + ((-t50 - t102) * t55 + t58 * t52) * t95 + t101 * (qJ(3) + t58)) (-m(4) + t78) * (g(1) * t30 + g(2) * t28 - g(3) * t85) t78 * (g(3) * t87 + t101) -m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (rSges(6,1) * t64 + rSges(6,2) * t65) + g(3) * (rSges(6,1) * t14 - rSges(6,2) * t15)) - m(7) * (g(2) * (t64 * t66 - t65 * t94) + (t66 * t14 + t94 * t15) * g(3) + (-t66 * t5 + t94 * t6) * g(1)) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (rSges(7,1) * t100 + rSges(7,2) * t99) + g(3) * ((-t15 * t51 + t54 * t87) * rSges(7,1) + (-t15 * t54 - t51 * t87) * rSges(7,2)))];
taug  = t3(:);
