% Calculate Gravitation load on the joints for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:08:12
% EndTime: 2019-03-09 07:08:14
% DurationCPUTime: 0.69s
% Computational Cost: add. (542->131), mult. (450->179), div. (0->0), fcn. (407->12), ass. (0->72)
t106 = rSges(6,3) + pkin(9);
t46 = pkin(11) + qJ(3);
t41 = qJ(4) + t46;
t35 = sin(t41);
t36 = cos(t41);
t105 = t36 * rSges(5,1) - t35 * rSges(5,2);
t52 = sin(qJ(1));
t54 = cos(qJ(1));
t104 = g(1) * t54 + g(2) * t52;
t53 = cos(qJ(5));
t38 = t53 * pkin(5) + pkin(4);
t55 = -pkin(10) - pkin(9);
t60 = t36 * t38 + (rSges(7,3) - t55) * t35;
t63 = t36 * pkin(4) + t106 * t35;
t47 = qJ(5) + qJ(6);
t43 = cos(t47);
t83 = t54 * t43;
t42 = sin(t47);
t88 = t52 * t42;
t5 = t36 * t88 + t83;
t84 = t54 * t42;
t87 = t52 * t43;
t6 = -t36 * t87 + t84;
t103 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t7 = -t36 * t84 + t87;
t8 = t36 * t83 + t88;
t102 = t7 * rSges(7,1) - t8 * rSges(7,2);
t39 = sin(t46);
t101 = pkin(3) * t39;
t51 = sin(qJ(5));
t100 = pkin(5) * t51;
t97 = g(3) * t35;
t49 = cos(pkin(11));
t37 = t49 * pkin(2) + pkin(1);
t96 = rSges(6,1) * t53;
t95 = rSges(7,1) * t43;
t94 = rSges(6,2) * t51;
t93 = rSges(7,2) * t42;
t90 = t36 * t52;
t89 = t36 * t54;
t86 = t52 * t51;
t85 = t52 * t53;
t82 = t54 * t51;
t81 = t54 * t53;
t50 = -pkin(7) - qJ(2);
t80 = rSges(4,3) - t50;
t45 = -pkin(8) + t50;
t79 = rSges(5,3) - t45;
t78 = rSges(3,3) + qJ(2);
t77 = t35 * t94;
t76 = t35 * t93;
t75 = t106 * t90 + t52 * t77;
t74 = t106 * t89 + t54 * t77;
t73 = -pkin(4) - t96;
t72 = -t45 + t100;
t71 = -t38 - t95;
t40 = cos(t46);
t69 = t40 * rSges(4,1) - t39 * rSges(4,2);
t66 = -rSges(7,1) * t42 - rSges(7,2) * t43;
t65 = rSges(3,1) * t49 - rSges(3,2) * sin(pkin(11)) + pkin(1);
t11 = -t36 * t82 + t85;
t9 = t36 * t86 + t81;
t64 = t37 + t69;
t62 = t63 + (-t94 + t96) * t36;
t59 = g(1) * (rSges(7,3) * t89 + t54 * t76) + g(2) * (rSges(7,3) * t90 + t52 * t76);
t58 = t60 + (-t93 + t95) * t36;
t34 = pkin(3) * t40;
t18 = t34 + t37;
t13 = t54 * t18;
t12 = t36 * t81 + t86;
t10 = -t36 * t85 + t82;
t1 = [-m(2) * (g(1) * (-t52 * rSges(2,1) - t54 * rSges(2,2)) + g(2) * (t54 * rSges(2,1) - t52 * rSges(2,2))) - m(3) * ((g(1) * t78 + g(2) * t65) * t54 + (-g(1) * t65 + g(2) * t78) * t52) - m(4) * ((g(1) * t80 + g(2) * t64) * t54 + (-g(1) * t64 + g(2) * t80) * t52) - m(5) * (g(2) * t13 + (g(1) * t79 + g(2) * t105) * t54 + (g(1) * (-t18 - t105) + g(2) * t79) * t52) - m(6) * (g(1) * (t10 * rSges(6,1) + t9 * rSges(6,2)) + g(2) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t13) + (-g(1) * t45 + g(2) * t63) * t54 + (g(1) * (-t18 - t63) - g(2) * t45) * t52) - m(7) * (g(1) * (t6 * rSges(7,1) + t5 * rSges(7,2)) + g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) + t13) + (g(1) * t72 + g(2) * t60) * t54 + (g(1) * (-t18 - t60) + g(2) * t72) * t52) (-m(3) - m(4) - m(5) - m(6) - m(7)) * (g(1) * t52 - g(2) * t54) -m(4) * (g(3) * t69 + t104 * (-rSges(4,1) * t39 - rSges(4,2) * t40)) - m(5) * (g(3) * (t34 + t105) + t104 * (-rSges(5,1) * t35 - rSges(5,2) * t36 - t101)) - m(6) * (g(1) * (-t54 * t101 + t74) + g(2) * (-t52 * t101 + t75) + g(3) * (t34 + t62) + t104 * t35 * t73) - m(7) * (g(3) * (t34 + t58) + t59 + t104 * (t71 * t35 - t36 * t55 - t101)) -m(5) * g(3) * t105 - m(6) * (g(1) * t74 + g(2) * t75 + g(3) * t62) - m(7) * (g(3) * t58 + t59) + t104 * ((m(5) * rSges(5,2) + m(7) * t55) * t36 + (m(5) * rSges(5,1) - m(6) * t73 - m(7) * t71) * t35) -m(6) * (g(1) * (t11 * rSges(6,1) - t12 * rSges(6,2)) + g(2) * (-t9 * rSges(6,1) + t10 * rSges(6,2))) - m(7) * (g(1) * (t11 * pkin(5) + t102) + g(2) * (-t9 * pkin(5) + t103)) + (-m(6) * (-rSges(6,1) * t51 - rSges(6,2) * t53) - m(7) * (t66 - t100)) * t97, -m(7) * (g(1) * t102 + g(2) * t103 + t66 * t97)];
taug  = t1(:);
