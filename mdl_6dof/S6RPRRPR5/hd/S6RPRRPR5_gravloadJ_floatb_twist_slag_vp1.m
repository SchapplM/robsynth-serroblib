% Calculate Gravitation load on the joints for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:12:13
% EndTime: 2019-03-09 05:12:14
% DurationCPUTime: 0.62s
% Computational Cost: add. (447->116), mult. (374->151), div. (0->0), fcn. (325->10), ass. (0->67)
t42 = sin(qJ(6));
t44 = cos(qJ(6));
t97 = rSges(7,1) * t42 + rSges(7,2) * t44;
t38 = pkin(10) + qJ(3);
t35 = qJ(4) + t38;
t31 = cos(t35);
t87 = rSges(7,3) + pkin(9);
t96 = t31 * t87;
t95 = t97 * t31;
t30 = sin(t35);
t94 = t31 * rSges(5,1) - t30 * rSges(5,2);
t52 = -t31 * rSges(6,2) + t30 * rSges(6,3);
t45 = cos(qJ(1));
t43 = sin(qJ(1));
t89 = g(2) * t43;
t93 = g(1) * t45 + t89;
t92 = -m(6) - m(7);
t33 = sin(t38);
t91 = pkin(3) * t33;
t88 = g(3) * t31;
t28 = t31 * pkin(4);
t41 = -pkin(7) - qJ(2);
t37 = -pkin(8) + t41;
t86 = pkin(5) - t37;
t40 = cos(pkin(10));
t32 = t40 * pkin(2) + pkin(1);
t82 = t30 * t43;
t81 = t30 * t45;
t79 = t31 * t45;
t78 = t43 * t42;
t77 = t43 * t44;
t76 = t45 * t42;
t75 = t45 * t44;
t74 = rSges(6,1) - t37;
t73 = rSges(4,3) - t41;
t72 = rSges(5,3) - t37;
t23 = t30 * qJ(5);
t71 = t23 + t28;
t70 = qJ(5) * t31;
t69 = rSges(3,3) + qJ(2);
t68 = -pkin(4) - t87;
t12 = t43 * t70;
t67 = t95 * t43 + t12;
t14 = t45 * t70;
t66 = t95 * t45 + t14;
t34 = cos(t38);
t29 = pkin(3) * t34;
t11 = t29 + t32;
t6 = t45 * t11;
t65 = pkin(4) * t79 + t45 * t23 + t6;
t62 = t43 * t31 * rSges(6,3) + rSges(6,2) * t82 + t12;
t61 = rSges(6,2) * t81 + rSges(6,3) * t79 + t14;
t59 = -t11 - t23;
t58 = g(1) * t68;
t57 = -pkin(4) * t30 - t91;
t56 = t34 * rSges(4,1) - t33 * rSges(4,2);
t53 = -rSges(5,1) * t30 - rSges(5,2) * t31;
t51 = t71 + t52;
t50 = t97 * t30 + t71 + t96;
t49 = rSges(3,1) * t40 - rSges(3,2) * sin(pkin(10)) + pkin(1);
t48 = t32 + t56;
t46 = (t45 * t58 + t68 * t89) * t30;
t5 = -t30 * t78 + t75;
t4 = t30 * t77 + t76;
t3 = t30 * t76 + t77;
t2 = t30 * t75 - t78;
t1 = [-m(2) * (g(1) * (-t43 * rSges(2,1) - t45 * rSges(2,2)) + g(2) * (t45 * rSges(2,1) - t43 * rSges(2,2))) - m(3) * ((g(1) * t69 + g(2) * t49) * t45 + (-g(1) * t49 + g(2) * t69) * t43) - m(4) * ((g(1) * t73 + g(2) * t48) * t45 + (-g(1) * t48 + g(2) * t73) * t43) - m(5) * (g(2) * t6 + (g(1) * t72 + g(2) * t94) * t45 + (g(1) * (-t11 - t94) + g(2) * t72) * t43) - m(6) * (g(2) * t65 + (g(1) * t74 + g(2) * t52) * t45 + (g(1) * (-t52 + t59 - t28) + g(2) * t74) * t43) - m(7) * (g(1) * (t5 * rSges(7,1) - t4 * rSges(7,2)) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t65) + (g(1) * t86 + g(2) * t96) * t45 + (g(1) * t59 + g(2) * t86 + t31 * t58) * t43) (-m(3) - m(4) - m(5) + t92) * (g(1) * t43 - g(2) * t45) -m(4) * (g(3) * t56 + t93 * (-rSges(4,1) * t33 - rSges(4,2) * t34)) - m(5) * (g(3) * (t29 + t94) + t93 * (t53 - t91)) - m(6) * (g(1) * (t45 * t57 + t61) + g(2) * (t43 * t57 + t62) + g(3) * (t29 + t51)) - m(7) * (g(1) * (-t45 * t91 + t66) + g(2) * (-t43 * t91 + t67) + g(3) * (t29 + t50) + t46) -m(5) * (g(3) * t94 + t93 * t53) - m(6) * (g(1) * (-pkin(4) * t81 + t61) + g(2) * (-pkin(4) * t82 + t62) + g(3) * t51) - m(7) * (g(1) * t66 + g(2) * t67 + g(3) * t50 + t46) t92 * (t93 * t30 - t88) -m(7) * (g(1) * (t2 * rSges(7,1) - t3 * rSges(7,2)) + g(2) * (t4 * rSges(7,1) + t5 * rSges(7,2)) + (-rSges(7,1) * t44 + rSges(7,2) * t42) * t88)];
taug  = t1(:);
