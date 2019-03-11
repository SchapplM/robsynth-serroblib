% Calculate Gravitation load on the joints for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:03:00
% EndTime: 2019-03-09 06:03:02
% DurationCPUTime: 0.73s
% Computational Cost: add. (549->130), mult. (515->175), div. (0->0), fcn. (501->10), ass. (0->64)
t84 = rSges(7,1) + pkin(5);
t44 = cos(qJ(3));
t41 = sin(qJ(3));
t83 = rSges(5,3) + pkin(8);
t65 = t83 * t41;
t96 = t44 * pkin(3) + t65;
t38 = qJ(1) + pkin(10);
t33 = sin(t38);
t34 = cos(t38);
t43 = cos(qJ(4));
t40 = sin(qJ(4));
t73 = t40 * t44;
t17 = t33 * t43 - t34 * t73;
t68 = rSges(7,3) + qJ(6);
t59 = rSges(4,1) * t44 - rSges(4,2) * t41;
t87 = g(2) * t33;
t95 = g(1) * t34 + t87;
t39 = qJ(4) + qJ(5);
t35 = sin(t39);
t36 = cos(t39);
t94 = t68 * t35 + t84 * t36;
t76 = t35 * t44;
t11 = t33 * t76 + t34 * t36;
t74 = t36 * t44;
t12 = t33 * t74 - t34 * t35;
t93 = -t11 * rSges(6,1) - t12 * rSges(6,2);
t13 = -t33 * t36 + t34 * t76;
t14 = t33 * t35 + t34 * t74;
t92 = -t13 * rSges(6,1) - t14 * rSges(6,2);
t42 = sin(qJ(1));
t91 = pkin(1) * t42;
t90 = pkin(4) * t40;
t89 = g(1) * t33;
t86 = g(3) * t41;
t80 = t33 * t40;
t78 = t34 * t40;
t77 = t35 * t41;
t72 = t43 * t44;
t32 = pkin(4) * t43 + pkin(3);
t27 = t44 * t32;
t46 = -pkin(9) - pkin(8);
t71 = rSges(7,2) - t46;
t70 = rSges(6,3) - t46;
t69 = t68 * t36 * t41;
t45 = cos(qJ(1));
t37 = t45 * pkin(1);
t66 = t34 * pkin(2) + t33 * pkin(7) + t37;
t64 = t34 * pkin(7) - t91;
t63 = -pkin(2) - t27;
t62 = t34 * t71;
t61 = t70 * t41;
t60 = pkin(4) * t80 + t34 * t27 + t66;
t57 = rSges(6,1) * t36 - rSges(6,2) * t35;
t56 = -rSges(6,1) * t35 - rSges(6,2) * t36;
t55 = -t84 * t11 + t68 * t12;
t54 = t17 * pkin(4);
t53 = t33 * t41 * t46 + pkin(4) * t78 + t64;
t52 = -t84 * t13 + t68 * t14;
t51 = rSges(5,1) * t43 - rSges(5,2) * t40 + pkin(3);
t15 = t33 * t73 + t34 * t43;
t50 = t15 * pkin(4);
t18 = t34 * t72 + t80;
t16 = -t33 * t72 + t78;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t42 - rSges(2,2) * t45) + g(2) * (rSges(2,1) * t45 - rSges(2,2) * t42)) - m(3) * (g(1) * (-rSges(3,1) * t33 - rSges(3,2) * t34 - t91) + g(2) * (rSges(3,1) * t34 - rSges(3,2) * t33 + t37)) - m(4) * (g(1) * (rSges(4,3) * t34 + t64) + g(2) * (t34 * t59 + t66) + (g(1) * (-pkin(2) - t59) + g(2) * rSges(4,3)) * t33) - m(5) * (g(1) * (rSges(5,1) * t16 + rSges(5,2) * t15 + t64) + (-pkin(2) - t96) * t89 + (rSges(5,1) * t18 + rSges(5,2) * t17 + t96 * t34 + t66) * g(2)) - m(6) * (g(1) * (-rSges(6,1) * t12 + rSges(6,2) * t11 + t53) + g(2) * (t14 * rSges(6,1) - t13 * rSges(6,2) + t34 * t61 + t60) + (-t41 * rSges(6,3) + t63) * t89) - m(7) * (g(1) * (-t68 * t11 - t12 * t84 + t53) + g(2) * (t68 * t13 + t84 * t14 + t41 * t62 + t60) + (-t41 * rSges(7,2) + t63) * t89) (-m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(4) * (g(3) * t59 + t95 * (-rSges(4,1) * t41 - rSges(4,2) * t44)) - m(5) * (g(3) * (t51 * t44 + t65) + t95 * (-t51 * t41 + t83 * t44)) - m(6) * (g(3) * (t57 * t44 + t27 + t61) + t95 * (t70 * t44 + (-t32 - t57) * t41)) - m(7) * (g(3) * t27 + (g(1) * t62 + g(3) * t94 + t71 * t87) * t44 + (g(3) * t71 + t95 * (-t32 - t94)) * t41) -m(5) * (g(1) * (rSges(5,1) * t17 - rSges(5,2) * t18) + g(2) * (-rSges(5,1) * t15 + rSges(5,2) * t16)) - m(6) * (g(1) * (t54 + t92) + g(2) * (-t50 + t93)) - m(7) * (g(1) * (t52 + t54) + g(2) * (-t50 + t55) + g(3) * t69) + (-m(5) * (-rSges(5,1) * t40 - rSges(5,2) * t43) - m(6) * (t56 - t90) - m(7) * (-t35 * t84 - t90)) * t86, -m(6) * (g(1) * t92 + g(2) * t93 + t56 * t86) - m(7) * (g(1) * t52 + g(2) * t55 + g(3) * (-t84 * t77 + t69)) -m(7) * (g(1) * t13 + g(2) * t11 + g(3) * t77)];
taug  = t1(:);
