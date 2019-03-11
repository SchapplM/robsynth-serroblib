% Calculate Gravitation load on the joints for
% S6RRPRRP4
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
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:52:35
% EndTime: 2019-03-09 11:52:37
% DurationCPUTime: 0.82s
% Computational Cost: add. (542->143), mult. (564->184), div. (0->0), fcn. (547->10), ass. (0->71)
t88 = rSges(7,1) + pkin(5);
t38 = qJ(2) + pkin(10);
t33 = sin(t38);
t86 = rSges(5,3) + pkin(8);
t99 = t33 * t86;
t34 = cos(t38);
t43 = sin(qJ(1));
t44 = cos(qJ(4));
t76 = t43 * t44;
t41 = sin(qJ(4));
t46 = cos(qJ(1));
t78 = t41 * t46;
t17 = -t34 * t78 + t76;
t68 = rSges(7,3) + qJ(6);
t98 = g(1) * t46 + g(2) * t43;
t39 = qJ(4) + qJ(5);
t35 = sin(t39);
t36 = cos(t39);
t97 = t35 * t68 + t36 * t88;
t80 = t36 * t46;
t81 = t35 * t43;
t11 = t34 * t81 + t80;
t74 = t46 * t35;
t77 = t43 * t36;
t12 = t34 * t77 - t74;
t96 = -t11 * rSges(6,1) - t12 * rSges(6,2);
t13 = t34 * t74 - t77;
t14 = t34 * t80 + t81;
t95 = -t13 * rSges(6,1) - t14 * rSges(6,2);
t42 = sin(qJ(2));
t94 = pkin(2) * t42;
t93 = pkin(4) * t41;
t40 = -qJ(3) - pkin(7);
t91 = g(2) * t40;
t89 = g(3) * t33;
t87 = rSges(3,3) + pkin(7);
t85 = t33 * t35;
t83 = t33 * t46;
t47 = -pkin(9) - pkin(8);
t82 = t33 * t47;
t31 = pkin(4) * t44 + pkin(3);
t21 = t34 * t31;
t79 = t41 * t43;
t75 = t44 * t46;
t73 = rSges(7,2) - t47;
t72 = rSges(4,3) - t40;
t71 = rSges(6,3) - t47;
t70 = t68 * t33 * t36;
t45 = cos(qJ(2));
t37 = t45 * pkin(2);
t69 = t21 + t37;
t32 = t37 + pkin(1);
t66 = -t32 - t21;
t65 = pkin(4) * t78 - t46 * t40 + t43 * t82;
t64 = rSges(3,1) * t45 - rSges(3,2) * t42;
t62 = rSges(4,1) * t34 - rSges(4,2) * t33;
t61 = rSges(6,1) * t36 - rSges(6,2) * t35;
t60 = -rSges(6,1) * t35 - rSges(6,2) * t36;
t59 = -t88 * t11 + t68 * t12;
t58 = t17 * pkin(4);
t57 = -t88 * t13 + t68 * t14;
t56 = pkin(1) + t64;
t55 = rSges(5,1) * t44 - rSges(5,2) * t41 + pkin(3);
t15 = t34 * t79 + t75;
t25 = t46 * t32;
t54 = pkin(4) * t79 + t25 + (t21 - t82) * t46;
t53 = t34 * pkin(3) + t99;
t52 = t15 * pkin(4);
t18 = t34 * t75 + t79;
t16 = -t34 * t76 + t78;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t43 - rSges(2,2) * t46) + g(2) * (rSges(2,1) * t46 - rSges(2,2) * t43)) - m(3) * ((g(1) * t87 + g(2) * t56) * t46 + (-g(1) * t56 + g(2) * t87) * t43) - m(4) * (g(2) * t25 + (g(1) * t72 + g(2) * t62) * t46 + (g(1) * (-t32 - t62) + g(2) * t72) * t43) - m(5) * (g(1) * (rSges(5,1) * t16 + rSges(5,2) * t15) + g(2) * (rSges(5,1) * t18 + rSges(5,2) * t17 + t25) + (-g(1) * t40 + g(2) * t53) * t46 + (g(1) * (-t32 - t53) - t91) * t43) - m(6) * (g(1) * (-rSges(6,1) * t12 + rSges(6,2) * t11 + t65) + g(2) * (t14 * rSges(6,1) - t13 * rSges(6,2) + rSges(6,3) * t83 + t54) + (g(1) * (-t33 * rSges(6,3) + t66) - t91) * t43) - m(7) * (g(1) * (-t68 * t11 - t12 * t88 + t65) + g(2) * (rSges(7,2) * t83 + t68 * t13 + t14 * t88 + t54) + (g(1) * (-t33 * rSges(7,2) + t66) - t91) * t43) -m(3) * (g(3) * t64 + t98 * (-rSges(3,1) * t42 - rSges(3,2) * t45)) - m(4) * (g(3) * (t37 + t62) + t98 * (-rSges(4,1) * t33 - rSges(4,2) * t34 - t94)) - m(5) * (g(3) * (t34 * t55 + t37 + t99) + t98 * (-t33 * t55 + t34 * t86 - t94)) - m(6) * (g(3) * (t33 * t71 + t34 * t61 + t69) + t98 * (-t94 + t71 * t34 + (-t31 - t61) * t33)) - m(7) * (g(3) * t69 - t98 * t94 + (g(3) * t97 + t73 * t98) * t34 + (g(3) * t73 + t98 * (-t31 - t97)) * t33) (-m(4) - m(5) - m(6) - m(7)) * (g(1) * t43 - g(2) * t46) -m(5) * (g(1) * (rSges(5,1) * t17 - rSges(5,2) * t18) + g(2) * (-rSges(5,1) * t15 + rSges(5,2) * t16)) - m(6) * (g(1) * (t58 + t95) + g(2) * (-t52 + t96)) - m(7) * (g(1) * (t57 + t58) + g(2) * (-t52 + t59) + g(3) * t70) + (-m(5) * (-rSges(5,1) * t41 - rSges(5,2) * t44) - m(6) * (t60 - t93) - m(7) * (-t35 * t88 - t93)) * t89, -m(6) * (g(1) * t95 + g(2) * t96 + t60 * t89) - m(7) * (g(1) * t57 + g(2) * t59 + g(3) * (-t88 * t85 + t70)) -m(7) * (g(1) * t13 + g(2) * t11 + g(3) * t85)];
taug  = t1(:);
