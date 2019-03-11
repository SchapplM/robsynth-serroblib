% Calculate Gravitation load on the joints for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:56:09
% EndTime: 2019-03-08 21:56:12
% DurationCPUTime: 0.97s
% Computational Cost: add. (620->140), mult. (1016->208), div. (0->0), fcn. (1182->14), ass. (0->67)
t40 = sin(qJ(3));
t43 = cos(qJ(3));
t70 = cos(pkin(6));
t37 = sin(pkin(6));
t41 = sin(qJ(2));
t78 = t37 * t41;
t99 = -t40 * t78 + t70 * t43;
t44 = cos(qJ(2));
t36 = sin(pkin(11));
t65 = t36 * t70;
t69 = cos(pkin(11));
t21 = -t41 * t65 + t69 * t44;
t77 = t37 * t43;
t98 = -t21 * t40 + t36 * t77;
t76 = t37 * t44;
t97 = g(3) * t76;
t56 = t70 * t69;
t19 = t36 * t44 + t41 * t56;
t64 = t37 * t69;
t50 = -t19 * t40 - t43 * t64;
t49 = t50 * pkin(3);
t82 = rSges(6,3) + pkin(9);
t73 = rSges(7,3) + pkin(10) + pkin(9);
t39 = sin(qJ(5));
t42 = cos(qJ(5));
t96 = t39 * rSges(6,1) + t42 * rSges(6,2);
t20 = t69 * t41 + t44 * t65;
t18 = t36 * t41 - t44 * t56;
t88 = g(2) * t18;
t95 = g(1) * t20 + t88;
t34 = qJ(3) + pkin(12);
t30 = sin(t34);
t31 = cos(t34);
t79 = t36 * t37;
t10 = -t21 * t30 + t31 * t79;
t14 = -t30 * t78 + t70 * t31;
t8 = -t19 * t30 - t31 * t64;
t94 = g(1) * t10 + g(2) * t8 + g(3) * t14;
t35 = qJ(5) + qJ(6);
t32 = sin(t35);
t33 = cos(t35);
t52 = rSges(7,1) * t33 - rSges(7,2) * t32 + pkin(5) * t42 + pkin(4);
t93 = t73 * t30 + t52 * t31;
t54 = rSges(6,1) * t42 - rSges(6,2) * t39 + pkin(4);
t92 = t82 * t30 + t54 * t31;
t87 = g(2) * t19;
t29 = pkin(3) * t43 + pkin(2);
t85 = t29 * t97;
t84 = g(3) * t37;
t83 = rSges(4,3) + pkin(8);
t38 = -qJ(4) - pkin(8);
t72 = -t18 * t29 - t19 * t38;
t71 = -t20 * t29 - t21 * t38;
t68 = -m(5) - m(6) - m(7);
t62 = t98 * pkin(3);
t9 = t19 * t31 - t30 * t64;
t60 = t18 * t42 - t39 * t9;
t59 = rSges(5,1) * t31 - rSges(5,2) * t30;
t11 = t21 * t31 + t30 * t79;
t58 = -t11 * t39 + t20 * t42;
t57 = t99 * pkin(3);
t55 = rSges(4,1) * t43 - rSges(4,2) * t40 + pkin(2);
t15 = t70 * t30 + t31 * t78;
t53 = -t15 * t39 - t42 * t76;
t51 = t32 * rSges(7,1) + t33 * rSges(7,2) + t39 * pkin(5);
t47 = m(7) * (g(1) * ((-t11 * t32 + t20 * t33) * rSges(7,1) + (-t11 * t33 - t20 * t32) * rSges(7,2)) + g(2) * ((t18 * t33 - t32 * t9) * rSges(7,1) + (-t18 * t32 - t33 * t9) * rSges(7,2)) + g(3) * ((-t15 * t32 - t33 * t76) * rSges(7,1) + (-t15 * t33 + t32 * t76) * rSges(7,2)));
t1 = [(-m(2) - m(3) - m(4) + t68) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t20 - rSges(3,2) * t21) + g(2) * (-rSges(3,1) * t18 - rSges(3,2) * t19) + (rSges(3,1) * t44 - rSges(3,2) * t41) * t84) - m(4) * (g(1) * (-t55 * t20 + t83 * t21) + t83 * t87 - t55 * t88 + (t83 * t41 + t55 * t44) * t84) - m(5) * (g(1) * (rSges(5,3) * t21 - t59 * t20 + t71) + g(2) * (rSges(5,3) * t19 - t59 * t18 + t72) + t85 + (t59 * t44 + (rSges(5,3) - t38) * t41) * t84) - m(6) * (g(1) * (t21 * t96 + t71) + g(2) * (t96 * t19 + t72) + t85 + ((-t38 + t96) * t41 + t92 * t44) * t84 - t95 * t92) - m(7) * (g(2) * t72 + t85 + t51 * t87 + ((-t38 + t51) * t41 + t93 * t44) * t84 - t95 * t93 + (t51 * t21 + t71) * g(1)) -m(4) * (g(1) * (t98 * rSges(4,1) + (-t21 * t43 - t40 * t79) * rSges(4,2)) + g(2) * (t50 * rSges(4,1) + (-t19 * t43 + t40 * t64) * rSges(4,2)) + g(3) * (t99 * rSges(4,1) + (-t70 * t40 - t41 * t77) * rSges(4,2))) - m(5) * (g(1) * (rSges(5,1) * t10 - rSges(5,2) * t11 + t62) + g(2) * (t8 * rSges(5,1) - t9 * rSges(5,2) + t49) + g(3) * (rSges(5,1) * t14 - rSges(5,2) * t15 + t57)) + (-g(1) * (t11 * t73 + t62) - g(2) * (t73 * t9 + t49) - g(3) * (t73 * t15 + t57) - t94 * t52) * m(7) + (-g(1) * (t11 * t82 + t62) - g(2) * (t82 * t9 + t49) - g(3) * (t15 * t82 + t57) - t94 * t54) * m(6), t68 * (t95 - t97) -m(6) * (g(1) * (t58 * rSges(6,1) + (-t11 * t42 - t20 * t39) * rSges(6,2)) + g(2) * (t60 * rSges(6,1) + (-t18 * t39 - t42 * t9) * rSges(6,2)) + g(3) * (t53 * rSges(6,1) + (-t15 * t42 + t39 * t76) * rSges(6,2))) - t47 - m(7) * (g(1) * t58 + g(2) * t60 + g(3) * t53) * pkin(5), -t47];
taug  = t1(:);
