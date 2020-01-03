% Calculate Gravitation load on the joints for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:09
% EndTime: 2019-12-31 22:08:13
% DurationCPUTime: 1.08s
% Computational Cost: add. (461->154), mult. (1097->234), div. (0->0), fcn. (1305->10), ass. (0->67)
t46 = sin(qJ(2));
t47 = sin(qJ(1));
t50 = cos(qJ(2));
t64 = cos(pkin(5));
t81 = cos(qJ(1));
t55 = t64 * t81;
t28 = t46 * t55 + t47 * t50;
t45 = sin(qJ(3));
t49 = cos(qJ(3));
t42 = sin(pkin(5));
t61 = t42 * t81;
t14 = t28 * t49 - t45 * t61;
t27 = t46 * t47 - t50 * t55;
t44 = sin(qJ(4));
t48 = cos(qJ(4));
t3 = t14 * t44 - t27 * t48;
t96 = -t14 * t48 - t27 * t44;
t57 = t47 * t64;
t29 = t81 * t46 + t50 * t57;
t95 = g(1) * t29 + g(2) * t27;
t40 = pkin(4) * t48 + pkin(3);
t67 = rSges(6,3) + qJ(5) + pkin(9);
t94 = t40 * t49 + t67 * t45;
t13 = t28 * t45 + t49 * t61;
t30 = -t46 * t57 + t81 * t50;
t73 = t42 * t49;
t17 = t30 * t45 - t47 * t73;
t75 = t42 * t46;
t25 = t45 * t75 - t64 * t49;
t93 = g(1) * t17 + g(2) * t13 + g(3) * t25;
t74 = t42 * t47;
t18 = t30 * t49 + t45 * t74;
t26 = t64 * t45 + t46 * t73;
t92 = g(1) * t18 + g(2) * t14 + g(3) * t26;
t91 = pkin(3) * t49;
t84 = g(3) * t42;
t83 = rSges(4,3) + pkin(8);
t82 = rSges(5,3) + pkin(9);
t78 = t28 * t44;
t77 = t30 * t44;
t72 = t42 * t50;
t71 = t44 * t46;
t70 = t44 * t49;
t69 = t48 * t49;
t68 = t49 * t50;
t66 = t81 * pkin(1) + pkin(7) * t74;
t65 = pkin(2) * t72 + pkin(8) * t75;
t63 = t30 * pkin(2) + t66;
t62 = pkin(4) * t44 + pkin(8);
t60 = -t47 * pkin(1) + pkin(7) * t61;
t21 = t27 * pkin(2);
t59 = pkin(8) * t28 - t21;
t23 = t29 * pkin(2);
t58 = pkin(8) * t30 - t23;
t56 = -t28 * pkin(2) + t60;
t54 = rSges(4,1) * t49 - rSges(4,2) * t45;
t5 = -t18 * t44 + t29 * t48;
t11 = -t26 * t44 - t48 * t72;
t20 = (t48 * t68 + t71) * t42;
t19 = (-t44 * t68 + t46 * t48) * t42;
t12 = -t26 * t48 + t44 * t72;
t10 = -t29 * t69 + t77;
t9 = t29 * t70 + t30 * t48;
t8 = -t27 * t69 + t78;
t7 = t27 * t70 + t28 * t48;
t6 = t18 * t48 + t29 * t44;
t1 = [-m(2) * (g(1) * (-t47 * rSges(2,1) - t81 * rSges(2,2)) + g(2) * (t81 * rSges(2,1) - t47 * rSges(2,2))) - m(3) * (g(1) * (-t28 * rSges(3,1) + t27 * rSges(3,2) + rSges(3,3) * t61 + t60) + g(2) * (rSges(3,1) * t30 - rSges(3,2) * t29 + rSges(3,3) * t74 + t66)) - m(4) * (g(1) * (-rSges(4,1) * t14 + rSges(4,2) * t13 - t83 * t27 + t56) + g(2) * (rSges(4,1) * t18 - rSges(4,2) * t17 + t83 * t29 + t63)) - m(5) * (g(1) * (rSges(5,1) * t96 + rSges(5,2) * t3 - pkin(3) * t14 - pkin(8) * t27 - t13 * t82 + t56) + g(2) * (rSges(5,1) * t6 + rSges(5,2) * t5 + pkin(3) * t18 + pkin(8) * t29 + t82 * t17 + t63)) - m(6) * (g(1) * (rSges(6,1) * t96 + rSges(6,2) * t3 - t13 * t67 - t14 * t40 - t62 * t27 + t56) + g(2) * (rSges(6,1) * t6 + rSges(6,2) * t5 + t67 * t17 + t18 * t40 + t62 * t29 + t63)), -m(3) * (g(1) * (-rSges(3,1) * t29 - rSges(3,2) * t30) + g(2) * (-rSges(3,1) * t27 - rSges(3,2) * t28) + (rSges(3,1) * t50 - rSges(3,2) * t46) * t84) - m(4) * (g(1) * (-t54 * t29 + t83 * t30 - t23) + g(2) * (-t54 * t27 + t83 * t28 - t21) + g(3) * t65 + (rSges(4,3) * t46 + t54 * t50) * t84) - m(6) * (g(1) * (rSges(6,1) * t10 + rSges(6,2) * t9 + pkin(4) * t77 + t58) + g(2) * (rSges(6,1) * t8 + rSges(6,2) * t7 + pkin(4) * t78 + t59) + g(3) * (t20 * rSges(6,1) + t19 * rSges(6,2) + t65) + (pkin(4) * t71 + t94 * t50) * t84 - t95 * t94) + (-g(1) * (rSges(5,1) * t10 + rSges(5,2) * t9 - t29 * t91 + t58) - g(2) * (rSges(5,1) * t8 + rSges(5,2) * t7 - t27 * t91 + t59) - g(3) * (pkin(3) * t42 * t68 + t20 * rSges(5,1) + t19 * rSges(5,2) + t65) - (g(3) * t72 - t95) * t45 * t82) * m(5), -m(4) * (g(1) * (-rSges(4,1) * t17 - rSges(4,2) * t18) + g(2) * (-rSges(4,1) * t13 - rSges(4,2) * t14) + g(3) * (-rSges(4,1) * t25 - rSges(4,2) * t26)) - m(5) * (t92 * t82 + t93 * (-rSges(5,1) * t48 + rSges(5,2) * t44 - pkin(3))) - m(6) * (t92 * t67 + t93 * (-rSges(6,1) * t48 + rSges(6,2) * t44 - t40)), -m(5) * (g(1) * (rSges(5,1) * t5 - rSges(5,2) * t6) + g(2) * (-rSges(5,1) * t3 + rSges(5,2) * t96) + g(3) * (rSges(5,1) * t11 + rSges(5,2) * t12)) + (-g(1) * (rSges(6,1) * t5 - rSges(6,2) * t6) - g(2) * (-rSges(6,1) * t3 + rSges(6,2) * t96) - g(3) * (t11 * rSges(6,1) + t12 * rSges(6,2)) - (g(1) * t5 - g(2) * t3 + g(3) * t11) * pkin(4)) * m(6), -m(6) * t93];
taug = t1(:);
