% Calculate Gravitation load on the joints for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR14_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:32:40
% EndTime: 2019-03-09 11:32:42
% DurationCPUTime: 0.89s
% Computational Cost: add. (499->171), mult. (1168->229), div. (0->0), fcn. (1368->10), ass. (0->65)
t85 = rSges(7,3) + pkin(10);
t46 = sin(qJ(2));
t47 = sin(qJ(1));
t50 = cos(qJ(2));
t51 = cos(qJ(1));
t76 = cos(pkin(6));
t69 = t51 * t76;
t28 = t46 * t69 + t47 * t50;
t70 = t47 * t76;
t30 = -t46 * t70 + t50 * t51;
t96 = g(1) * t30 + g(2) * t28;
t27 = t46 * t47 - t50 * t69;
t29 = t51 * t46 + t50 * t70;
t95 = g(1) * t29 + g(2) * t27;
t94 = -m(6) - m(7);
t45 = sin(qJ(4));
t93 = pkin(4) * t45;
t43 = sin(pkin(6));
t88 = g(3) * t43;
t87 = rSges(6,1) + pkin(9);
t86 = rSges(5,3) + pkin(9);
t84 = t43 * t46;
t83 = t43 * t47;
t82 = t43 * t50;
t81 = t43 * t51;
t80 = pkin(2) * t82 + qJ(3) * t84;
t79 = t51 * pkin(1) + pkin(8) * t83;
t78 = rSges(4,3) + qJ(3);
t77 = rSges(6,3) + qJ(5);
t75 = t30 * pkin(2) + t79;
t74 = pkin(9) * t82 + t80;
t73 = -t47 * pkin(1) + pkin(8) * t81;
t21 = t27 * pkin(2);
t72 = -t27 * pkin(9) - t21;
t23 = t29 * pkin(2);
t71 = -t29 * pkin(9) - t23;
t68 = -t28 * pkin(2) + t73;
t49 = cos(qJ(4));
t67 = rSges(5,1) * t45 + rSges(5,2) * t49;
t66 = g(3) * (t84 * t93 + t74);
t44 = sin(qJ(6));
t48 = cos(qJ(6));
t65 = rSges(7,1) * t48 - rSges(7,2) * t44 + pkin(5);
t64 = -t27 * t45 + t49 * t81;
t63 = t27 * t49 + t45 * t81;
t61 = rSges(7,1) * t44 + rSges(7,2) * t48 + qJ(5);
t60 = pkin(3) * t83 + qJ(3) * t29 + t75;
t59 = -pkin(9) - t65;
t58 = -rSges(6,2) * t45 - t77 * t49;
t10 = t29 * t45 + t49 * t83;
t57 = t10 * pkin(4) + t60;
t55 = pkin(3) * t81 - t27 * qJ(3) + t68;
t54 = pkin(4) * t64 + t55;
t53 = t85 * t45 - t61 * t49;
t26 = -t45 * t82 + t76 * t49;
t25 = t76 * t45 + t49 * t82;
t20 = t25 * pkin(4);
t16 = t30 * t93;
t15 = t28 * t93;
t9 = -t29 * t49 + t45 * t83;
t7 = t63 * pkin(4);
t5 = t9 * pkin(4);
t3 = t30 * t48 + t44 * t9;
t2 = -t30 * t44 + t48 * t9;
t1 = [-m(2) * (g(1) * (-t47 * rSges(2,1) - rSges(2,2) * t51) + g(2) * (rSges(2,1) * t51 - t47 * rSges(2,2))) - m(3) * (g(1) * (-t28 * rSges(3,1) + t27 * rSges(3,2) + rSges(3,3) * t81 + t73) + g(2) * (rSges(3,1) * t30 - rSges(3,2) * t29 + rSges(3,3) * t83 + t79)) - m(4) * (g(1) * (rSges(4,1) * t81 + t28 * rSges(4,2) - t78 * t27 + t68) + g(2) * (rSges(4,1) * t83 - rSges(4,2) * t30 + t78 * t29 + t75)) - m(5) * (g(1) * (rSges(5,1) * t64 - rSges(5,2) * t63 - t86 * t28 + t55) + g(2) * (rSges(5,1) * t10 - rSges(5,2) * t9 + t86 * t30 + t60)) - m(6) * (g(1) * (-rSges(6,2) * t64 - t87 * t28 + t63 * t77 + t54) + g(2) * (-rSges(6,2) * t10 + t87 * t30 + t77 * t9 + t57)) - m(7) * (g(1) * (t59 * t28 + t61 * t63 + t64 * t85 + t54) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + qJ(5) * t9 + (pkin(5) + pkin(9)) * t30 + t85 * t10 + t57)) -m(3) * (g(1) * (-rSges(3,1) * t29 - rSges(3,2) * t30) + g(2) * (-t27 * rSges(3,1) - rSges(3,2) * t28) + (rSges(3,1) * t50 - rSges(3,2) * t46) * t88) - m(4) * (g(1) * (rSges(4,2) * t29 + t78 * t30 - t23) + g(2) * (rSges(4,2) * t27 + t78 * t28 - t21) + g(3) * ((-rSges(4,2) * t50 + rSges(4,3) * t46) * t43 + t80)) - m(5) * (g(1) * (-rSges(5,3) * t29 + t71) + g(2) * (-rSges(5,3) * t27 + t72) + g(3) * t74 + (rSges(5,3) * t50 + t67 * t46) * t88 + t96 * (qJ(3) + t67)) - m(6) * (g(1) * (-rSges(6,1) * t29 + t16 + t71) + g(2) * (-rSges(6,1) * t27 + t15 + t72) + t66 + (rSges(6,1) * t50 + t58 * t46) * t88 + t96 * (qJ(3) + t58)) - m(7) * (g(1) * (t16 - t23) + g(2) * (t15 - t21) + t66 + (t53 * t46 + t65 * t50) * t88 + t95 * t59 + t96 * (qJ(3) + t53)) (-m(4) - m(5) + t94) * (-g(3) * t82 + t95) -m(5) * (g(1) * (-rSges(5,1) * t9 - rSges(5,2) * t10) + g(2) * (rSges(5,1) * t63 + rSges(5,2) * t64) + g(3) * (-rSges(5,1) * t25 - rSges(5,2) * t26)) - m(6) * (g(1) * (rSges(6,2) * t9 + t77 * t10 - t5) + g(2) * (-rSges(6,2) * t63 - t64 * t77 + t7) + g(3) * (rSges(6,2) * t25 + t77 * t26 - t20)) + (-g(1) * (-t85 * t9 - t5) - g(2) * (t63 * t85 + t7) - g(3) * (-t25 * t85 - t20) - (g(1) * t10 - g(2) * t64 + g(3) * t26) * t61) * m(7), t94 * (g(1) * t9 - g(2) * t63 + g(3) * t25) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * ((-t28 * t44 - t48 * t63) * rSges(7,1) + (-t28 * t48 + t44 * t63) * rSges(7,2)) + g(3) * ((t25 * t48 - t44 * t84) * rSges(7,1) + (-t25 * t44 - t48 * t84) * rSges(7,2)))];
taug  = t1(:);
