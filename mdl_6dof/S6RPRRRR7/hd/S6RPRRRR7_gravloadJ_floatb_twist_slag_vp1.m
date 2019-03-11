% Calculate Gravitation load on the joints for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:16:57
% EndTime: 2019-03-09 07:16:59
% DurationCPUTime: 0.68s
% Computational Cost: add. (403->133), mult. (387->163), div. (0->0), fcn. (332->10), ass. (0->69)
t34 = qJ(3) + qJ(4);
t31 = qJ(5) + t34;
t26 = cos(t31);
t85 = rSges(7,3) + pkin(10);
t86 = t85 * t26;
t25 = sin(t31);
t84 = t85 * t25;
t40 = cos(qJ(1));
t78 = g(2) * t40;
t37 = sin(qJ(1));
t83 = g(1) * t37 - t78;
t35 = sin(qJ(6));
t69 = rSges(7,2) * t35;
t55 = t26 * t69;
t82 = t55 * t78;
t41 = -pkin(8) - pkin(7);
t27 = sin(t34);
t81 = pkin(4) * t27;
t28 = cos(t34);
t80 = pkin(4) * t28;
t36 = sin(qJ(3));
t77 = t36 * pkin(3);
t39 = cos(qJ(3));
t76 = t39 * pkin(3);
t75 = rSges(4,3) + pkin(7);
t73 = rSges(5,1) * t28;
t72 = rSges(6,1) * t26;
t38 = cos(qJ(6));
t71 = rSges(7,1) * t38;
t70 = rSges(5,2) * t27;
t68 = t25 * t37;
t67 = t25 * t40;
t66 = t26 * t37;
t65 = t28 * t37;
t64 = t37 * t35;
t63 = t37 * t38;
t62 = t40 * t35;
t61 = t40 * t38;
t60 = rSges(5,3) - t41;
t59 = t40 * pkin(1) + t37 * qJ(2);
t10 = t77 + t81;
t30 = t40 * qJ(2);
t33 = -pkin(9) + t41;
t58 = t40 * t10 + t37 * t33 + t30;
t57 = t37 * t10 + t59;
t56 = t26 * t71;
t54 = -pkin(5) - t71;
t53 = rSges(6,1) * t66 - rSges(6,2) * t68;
t51 = t36 * rSges(4,1) + t39 * rSges(4,2);
t50 = -t27 * rSges(5,1) - t28 * rSges(5,2);
t49 = t25 * rSges(6,1) + t26 * rSges(6,2);
t48 = g(1) * t30 + g(2) * t59;
t47 = -t50 + t77;
t46 = -t49 - t81;
t45 = pkin(5) * t66 + t85 * t68 + (-t55 + t56) * t37;
t44 = t86 + (t54 + t69) * t25;
t43 = t54 * t26 - t84;
t42 = t44 - t81;
t21 = pkin(4) * t65;
t20 = t40 * t70;
t19 = rSges(5,1) * t65;
t15 = rSges(6,2) * t67;
t11 = t76 + t80;
t6 = t37 * t11;
t4 = t25 * t61 - t64;
t3 = t25 * t62 + t63;
t2 = t25 * t63 + t62;
t1 = -t25 * t64 + t61;
t5 = [-m(2) * (g(1) * (-t37 * rSges(2,1) - t40 * rSges(2,2)) + g(2) * (t40 * rSges(2,1) - t37 * rSges(2,2))) - m(3) * (g(1) * (t40 * rSges(3,3) + t30 + (rSges(3,2) - pkin(1)) * t37) + g(2) * (-t40 * rSges(3,2) + t37 * rSges(3,3) + t59)) - m(4) * ((g(1) * t51 + g(2) * t75) * t40 + (g(1) * (-pkin(1) - t75) + g(2) * t51) * t37 + t48) - m(5) * ((g(1) * t47 + g(2) * t60) * t40 + (g(1) * (-pkin(1) - t60) + g(2) * t47) * t37 + t48) - m(6) * (g(1) * t58 + g(2) * t57 + (g(1) * t49 + g(2) * (rSges(6,3) - t33)) * t40 + (g(1) * (-rSges(6,3) - pkin(1)) + g(2) * t49) * t37) - m(7) * (g(1) * (t4 * rSges(7,1) - t3 * rSges(7,2) - t37 * pkin(1) + pkin(5) * t67 + t58) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + pkin(5) * t68 - t40 * t33 + t57) - (g(1) * t40 + g(2) * t37) * t86) (-m(3) - m(4) - m(5) - m(6) - m(7)) * t83, -m(4) * (-g(3) * t51 + t83 * (rSges(4,1) * t39 - rSges(4,2) * t36)) - m(5) * (g(1) * (t19 + (-t70 + t76) * t37) + g(2) * (t20 + (-t73 - t76) * t40) - g(3) * t47) - m(6) * (g(1) * (t53 + t6) + g(2) * (t15 + (-t11 - t72) * t40) + g(3) * (t46 - t77)) - m(7) * (g(1) * (t45 + t6) + t82 + g(3) * (t42 - t77) + (-t11 + t43) * t78) -m(5) * (g(1) * (-t37 * t70 + t19) + g(2) * t20 + g(3) * t50) - m(6) * (g(1) * (t21 + t53) + g(2) * t15 + g(3) * t46) - m(7) * (g(1) * (t21 + t45) + t82 + g(3) * t42) + (m(5) * t73 - m(6) * (-t72 - t80) - m(7) * (-pkin(5) * t26 - t56 - t80 - t84)) * t78, -m(6) * (g(1) * t53 + g(2) * (-t40 * t72 + t15) - g(3) * t49) - m(7) * (g(1) * t45 + g(3) * t44 + t43 * t78 + t82) -m(7) * (g(1) * (t1 * rSges(7,1) - t2 * rSges(7,2)) + g(2) * (t3 * rSges(7,1) + t4 * rSges(7,2)) + g(3) * (-rSges(7,1) * t35 - rSges(7,2) * t38) * t26)];
taug  = t5(:);
