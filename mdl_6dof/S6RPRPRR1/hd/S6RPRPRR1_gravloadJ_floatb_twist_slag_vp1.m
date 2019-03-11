% Calculate Gravitation load on the joints for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:33:53
% EndTime: 2019-03-09 03:33:54
% DurationCPUTime: 0.51s
% Computational Cost: add. (423->103), mult. (305->132), div. (0->0), fcn. (259->12), ass. (0->60)
t87 = rSges(7,3) + pkin(9);
t34 = qJ(3) + pkin(11);
t30 = qJ(5) + t34;
t23 = sin(t30);
t24 = cos(t30);
t37 = sin(qJ(6));
t74 = rSges(7,2) * t37;
t86 = t23 * t74 + t24 * t87;
t84 = t24 * rSges(6,1) - t23 * rSges(6,2);
t26 = sin(t34);
t28 = cos(t34);
t41 = cos(qJ(3));
t31 = t41 * pkin(3);
t83 = t28 * rSges(5,1) - t26 * rSges(5,2) + t31;
t35 = qJ(1) + pkin(10);
t27 = sin(t35);
t29 = cos(t35);
t82 = g(1) * t29 + g(2) * t27;
t46 = t24 * pkin(5) + t87 * t23;
t38 = sin(qJ(3));
t79 = t38 * pkin(3);
t39 = sin(qJ(1));
t78 = t39 * pkin(1);
t77 = rSges(4,3) + pkin(7);
t64 = pkin(4) * t28 + t31;
t15 = pkin(2) + t64;
t42 = cos(qJ(1));
t32 = t42 * pkin(1);
t76 = t29 * t15 + t32;
t40 = cos(qJ(6));
t75 = rSges(7,1) * t40;
t70 = t27 * t37;
t69 = t27 * t40;
t68 = t29 * t37;
t67 = t29 * t40;
t36 = -qJ(4) - pkin(7);
t66 = rSges(5,3) - t36;
t33 = -pkin(8) + t36;
t65 = rSges(6,3) - t33;
t63 = -m(5) - m(6) - m(7);
t62 = g(1) * t78;
t61 = t86 * t27;
t60 = t86 * t29;
t56 = t41 * rSges(4,1) - t38 * rSges(4,2);
t52 = -rSges(6,1) * t23 - rSges(6,2) * t24;
t51 = g(2) * t32 - t62;
t50 = pkin(2) + t56;
t49 = pkin(2) + t83;
t48 = t52 * t27;
t47 = t52 * t29;
t45 = t46 + (-t74 + t75) * t24;
t43 = t82 * (-pkin(5) - t75) * t23;
t16 = -pkin(4) * t26 - t79;
t7 = t29 * t16;
t6 = t27 * t16;
t4 = t24 * t67 + t70;
t3 = -t24 * t68 + t69;
t2 = -t24 * t69 + t68;
t1 = t24 * t70 + t67;
t5 = [-m(2) * (g(1) * (-t39 * rSges(2,1) - t42 * rSges(2,2)) + g(2) * (t42 * rSges(2,1) - t39 * rSges(2,2))) - m(3) * (g(1) * (-t27 * rSges(3,1) - t29 * rSges(3,2) - t78) + g(2) * (t29 * rSges(3,1) - t27 * rSges(3,2) + t32)) - m(4) * ((g(1) * t77 + g(2) * t50) * t29 + (-g(1) * t50 + g(2) * t77) * t27 + t51) - m(5) * ((g(1) * t66 + g(2) * t49) * t29 + (-g(1) * t49 + g(2) * t66) * t27 + t51) - m(6) * (-t62 + g(2) * t76 + (g(1) * t65 + g(2) * t84) * t29 + (g(1) * (-t15 - t84) + g(2) * t65) * t27) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t78) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t76) + (-g(1) * t33 + g(2) * t46) * t29 + (g(1) * (-t15 - t46) - g(2) * t33) * t27) (-m(3) - m(4) + t63) * g(3), -m(4) * (g(3) * t56 + t82 * (-rSges(4,1) * t38 - rSges(4,2) * t41)) - m(5) * (g(3) * t83 + t82 * (-rSges(5,1) * t26 - rSges(5,2) * t28 - t79)) - m(6) * (g(1) * (t7 + t47) + g(2) * (t6 + t48) + g(3) * (t84 + t64)) - m(7) * (g(1) * (t7 + t60) + g(2) * (t6 + t61) + g(3) * (t45 + t64) + t43) t63 * (g(1) * t27 - g(2) * t29) -m(6) * (g(1) * t47 + g(2) * t48 + g(3) * t84) - m(7) * (g(1) * t60 + g(2) * t61 + g(3) * t45 + t43) -m(7) * (g(1) * (t3 * rSges(7,1) - t4 * rSges(7,2)) + g(2) * (-t1 * rSges(7,1) + t2 * rSges(7,2)) + g(3) * (-rSges(7,1) * t37 - rSges(7,2) * t40) * t23)];
taug  = t5(:);
