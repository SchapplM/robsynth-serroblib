% Calculate Gravitation load on the joints for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:43:06
% EndTime: 2019-03-09 11:43:08
% DurationCPUTime: 0.78s
% Computational Cost: add. (552->139), mult. (506->183), div. (0->0), fcn. (469->10), ass. (0->70)
t47 = cos(qJ(5));
t97 = rSges(7,1) + pkin(5);
t98 = -t97 * t47 - pkin(4);
t42 = qJ(2) + pkin(10);
t39 = qJ(4) + t42;
t33 = cos(t39);
t28 = t33 * rSges(5,1);
t32 = sin(t39);
t59 = -rSges(5,2) * t32 + t28;
t71 = t33 * pkin(4) + t32 * pkin(9);
t37 = sin(t42);
t38 = cos(t42);
t48 = cos(qJ(2));
t40 = t48 * pkin(2);
t96 = rSges(4,1) * t38 - rSges(4,2) * t37 + t40;
t46 = sin(qJ(1));
t49 = cos(qJ(1));
t95 = g(1) * t49 + g(2) * t46;
t68 = rSges(7,3) + qJ(6);
t94 = t95 * t32;
t45 = sin(qJ(2));
t93 = pkin(2) * t45;
t43 = -qJ(3) - pkin(7);
t41 = -pkin(8) + t43;
t91 = g(2) * t41;
t88 = rSges(3,3) + pkin(7);
t27 = t32 * rSges(7,2);
t26 = t32 * rSges(6,3);
t44 = sin(qJ(5));
t86 = t32 * t44;
t85 = t32 * t49;
t84 = t33 * t44;
t83 = t33 * t46;
t82 = t33 * t47;
t81 = t33 * t49;
t80 = t41 * t49;
t79 = t46 * t44;
t78 = t46 * t47;
t77 = t47 * t49;
t76 = t49 * t44;
t75 = rSges(4,3) - t43;
t74 = rSges(5,3) - t41;
t20 = pkin(9) * t83;
t73 = rSges(7,2) * t83 + t20;
t23 = pkin(9) * t81;
t72 = rSges(7,2) * t81 + t23;
t70 = pkin(3) * t38 + t40;
t69 = qJ(6) * t44;
t64 = rSges(6,2) * t86;
t67 = rSges(6,3) * t83 + t46 * t64 + t20;
t66 = rSges(6,3) * t81 + t49 * t64 + t23;
t10 = pkin(1) + t70;
t5 = t49 * t10;
t65 = pkin(4) * t81 + pkin(9) * t85 + t5;
t63 = -rSges(6,1) * t47 - pkin(4);
t62 = rSges(3,1) * t48 - rSges(3,2) * t45;
t58 = -rSges(5,1) * t32 - rSges(5,2) * t33;
t57 = -t10 - t71;
t56 = pkin(1) + t62;
t55 = pkin(1) + t96;
t54 = rSges(7,3) * t84 + t33 * t69 + t97 * t82 + t27 + t71;
t53 = rSges(6,1) * t82 - rSges(6,2) * t84 + t26 + t71;
t11 = -pkin(3) * t37 - t93;
t7 = t49 * t11;
t6 = t46 * t11;
t4 = t33 * t77 + t79;
t3 = t33 * t76 - t78;
t2 = t33 * t78 - t76;
t1 = t33 * t79 + t77;
t8 = [-m(2) * (g(1) * (-t46 * rSges(2,1) - rSges(2,2) * t49) + g(2) * (rSges(2,1) * t49 - t46 * rSges(2,2))) - m(3) * ((g(1) * t88 + g(2) * t56) * t49 + (-g(1) * t56 + g(2) * t88) * t46) - m(4) * ((g(1) * t75 + g(2) * t55) * t49 + (-g(1) * t55 + g(2) * t75) * t46) - m(5) * (g(2) * t5 + (g(1) * t74 + g(2) * t59) * t49 + (g(1) * (-t10 - t59) + g(2) * t74) * t46) - m(6) * (g(1) * (-t2 * rSges(6,1) + t1 * rSges(6,2) - t80) + g(2) * (t4 * rSges(6,1) - t3 * rSges(6,2) + rSges(6,3) * t85 + t65) + (g(1) * (t57 - t26) - t91) * t46) - m(7) * (g(1) * (-t68 * t1 - t97 * t2 - t80) + g(2) * (rSges(7,2) * t85 + t68 * t3 + t97 * t4 + t65) + (g(1) * (t57 - t27) - t91) * t46) -m(3) * (g(3) * t62 + t95 * (-rSges(3,1) * t45 - rSges(3,2) * t48)) - m(4) * (g(3) * t96 + t95 * (-rSges(4,1) * t37 - rSges(4,2) * t38 - t93)) - m(5) * (g(1) * (t58 * t49 + t7) + g(2) * (t58 * t46 + t6) + g(3) * (t70 + t59)) - m(6) * (g(1) * (t7 + t66) + g(2) * (t6 + t67) + g(3) * (t53 + t70) + t63 * t94) - m(7) * (g(1) * (t7 + t72) + g(2) * (t6 + t73) + g(3) * (t54 + t70) + (-t68 * t44 + t98) * t94) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t46 - g(2) * t49) -m(5) * (g(3) * t28 + (-g(1) * t81 - g(2) * t83) * rSges(5,2)) - m(6) * (g(1) * t66 + g(2) * t67 + g(3) * t53) - m(7) * (g(1) * t72 + g(2) * t73 + g(3) * t54) + (m(5) * g(3) * rSges(5,2) + t95 * (m(5) * rSges(5,1) - m(6) * t63 - m(7) * (-rSges(7,3) * t44 - t69 + t98))) * t32, -m(6) * (g(1) * (-rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2)) - m(7) * (g(1) * (-t3 * t97 + t68 * t4) + g(2) * (-t1 * t97 + t68 * t2)) + (-m(6) * (-rSges(6,1) * t44 - rSges(6,2) * t47) - m(7) * (-t97 * t44 + t68 * t47)) * g(3) * t32, -m(7) * (g(1) * t3 + g(2) * t1 + g(3) * t86)];
taug  = t8(:);
