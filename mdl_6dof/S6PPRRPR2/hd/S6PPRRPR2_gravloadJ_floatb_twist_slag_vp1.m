% Calculate Gravitation load on the joints for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:49:07
% EndTime: 2019-03-08 18:49:09
% DurationCPUTime: 0.69s
% Computational Cost: add. (703->114), mult. (1922->176), div. (0->0), fcn. (2447->14), ass. (0->66)
t38 = sin(qJ(6));
t41 = cos(qJ(6));
t96 = -t38 * rSges(7,1) - t41 * rSges(7,2);
t39 = sin(qJ(4));
t42 = cos(qJ(4));
t93 = -pkin(4) * t42 - qJ(5) * t39;
t92 = -rSges(7,3) - pkin(10);
t69 = sin(pkin(12));
t70 = sin(pkin(11));
t54 = t70 * t69;
t73 = cos(pkin(12));
t74 = cos(pkin(11));
t61 = t74 * t73;
t76 = cos(pkin(6));
t48 = -t76 * t61 + t54;
t71 = sin(pkin(7));
t72 = sin(pkin(6));
t58 = t72 * t71;
t75 = cos(pkin(7));
t91 = t48 * t75 + t74 * t58;
t55 = t70 * t73;
t59 = t74 * t69;
t49 = t76 * t55 + t59;
t57 = t72 * t70;
t90 = t49 * t75 - t71 * t57;
t89 = t73 * t75 * t72 + t76 * t71;
t88 = t96 * t39;
t87 = t41 * rSges(7,1) - t38 * rSges(7,2) + pkin(5) + pkin(9);
t86 = -m(6) - m(7);
t84 = rSges(6,1) + pkin(9);
t83 = rSges(5,3) + pkin(9);
t81 = cos(qJ(3));
t77 = rSges(6,3) + qJ(5);
t32 = t76 * t59 + t55;
t40 = sin(qJ(3));
t15 = t32 * t40 + t91 * t81;
t12 = t15 * pkin(3);
t68 = t15 * t93 - t12;
t33 = -t76 * t54 + t61;
t17 = t33 * t40 + t90 * t81;
t13 = t17 * pkin(3);
t67 = t17 * t93 - t13;
t56 = t72 * t69;
t26 = t40 * t56 - t89 * t81;
t25 = t26 * pkin(3);
t66 = t26 * t93 - t25;
t65 = -m(3) - m(4) - m(5) + t86;
t64 = -rSges(5,1) * t42 + rSges(5,2) * t39;
t63 = rSges(6,2) * t42 - rSges(6,3) * t39;
t60 = t74 * t72;
t47 = -t73 * t58 + t76 * t75;
t44 = t49 * t71 + t75 * t57;
t43 = t48 * t71 - t75 * t60;
t27 = t89 * t40 + t81 * t56;
t20 = t27 * t42 + t47 * t39;
t19 = t27 * t39 - t47 * t42;
t18 = t33 * t81 - t90 * t40;
t16 = t32 * t81 - t91 * t40;
t14 = t19 * pkin(4);
t7 = t18 * t42 + t44 * t39;
t6 = t18 * t39 - t44 * t42;
t5 = t16 * t42 + t43 * t39;
t4 = t16 * t39 - t43 * t42;
t3 = t6 * pkin(4);
t2 = t4 * pkin(4);
t1 = [(-m(2) + t65) * g(3), t65 * (g(1) * t57 - g(2) * t60 + g(3) * t76) -m(4) * (g(1) * (-t17 * rSges(4,1) - t18 * rSges(4,2)) + g(2) * (-t15 * rSges(4,1) - t16 * rSges(4,2)) + g(3) * (-t26 * rSges(4,1) - t27 * rSges(4,2))) - m(5) * (g(1) * (t64 * t17 + t83 * t18 - t13) + g(2) * (t64 * t15 + t83 * t16 - t12) + g(3) * (t64 * t26 + t83 * t27 - t25)) - m(6) * (g(1) * (t63 * t17 + t84 * t18 + t67) + g(2) * (t63 * t15 + t84 * t16 + t68) + g(3) * (t63 * t26 + t84 * t27 + t66)) + (-g(1) * (t88 * t17 + t87 * t18 + t67) - g(2) * (t88 * t15 + t87 * t16 + t68) - g(3) * (t88 * t26 + t87 * t27 + t66) - (g(1) * t17 + g(2) * t15 + g(3) * t26) * t42 * t92) * m(7), -m(5) * (g(1) * (-t6 * rSges(5,1) - t7 * rSges(5,2)) + g(2) * (-t4 * rSges(5,1) - t5 * rSges(5,2)) + g(3) * (-t19 * rSges(5,1) - t20 * rSges(5,2))) - m(6) * (g(1) * (t6 * rSges(6,2) + t77 * t7 - t3) + g(2) * (t4 * rSges(6,2) + t77 * t5 - t2) + g(3) * (t19 * rSges(6,2) + t77 * t20 - t14)) + (-g(1) * (t92 * t6 - t3) - g(2) * (t92 * t4 - t2) - g(3) * (t92 * t19 - t14) - (g(1) * t7 + g(2) * t5 + g(3) * t20) * (qJ(5) - t96)) * m(7), t86 * (g(1) * t6 + g(2) * t4 + g(3) * t19) -m(7) * (g(1) * ((-t17 * t38 + t6 * t41) * rSges(7,1) + (-t17 * t41 - t6 * t38) * rSges(7,2)) + g(2) * ((-t15 * t38 + t4 * t41) * rSges(7,1) + (-t15 * t41 - t4 * t38) * rSges(7,2)) + g(3) * ((t19 * t41 - t26 * t38) * rSges(7,1) + (-t19 * t38 - t26 * t41) * rSges(7,2)))];
taug  = t1(:);
