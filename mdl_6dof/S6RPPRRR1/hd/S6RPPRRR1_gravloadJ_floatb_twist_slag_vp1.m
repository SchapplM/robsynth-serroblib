% Calculate Gravitation load on the joints for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:17:57
% EndTime: 2019-03-09 02:17:58
% DurationCPUTime: 0.48s
% Computational Cost: add. (406->90), mult. (283->116), div. (0->0), fcn. (240->12), ass. (0->51)
t32 = qJ(1) + pkin(10);
t24 = sin(t32);
t26 = cos(t32);
t73 = g(1) * t26 + g(2) * t24;
t77 = rSges(7,3) + pkin(9);
t31 = pkin(11) + qJ(4);
t27 = qJ(5) + t31;
t20 = sin(t27);
t21 = cos(t27);
t74 = t21 * rSges(6,1) - t20 * rSges(6,2);
t43 = t21 * pkin(5) + t20 * t77;
t38 = cos(qJ(6));
t66 = rSges(7,1) * t38;
t72 = (-pkin(5) - t66) * t20;
t23 = sin(t31);
t71 = pkin(4) * t23;
t37 = sin(qJ(1));
t68 = t37 * pkin(1);
t34 = cos(pkin(11));
t22 = t34 * pkin(3) + pkin(2);
t25 = cos(t31);
t19 = pkin(4) * t25;
t11 = t19 + t22;
t39 = cos(qJ(1));
t29 = t39 * pkin(1);
t67 = t26 * t11 + t29;
t36 = sin(qJ(6));
t65 = rSges(7,2) * t36;
t61 = t24 * t36;
t60 = t24 * t38;
t59 = t26 * t36;
t58 = t26 * t38;
t35 = -pkin(7) - qJ(3);
t57 = rSges(5,3) - t35;
t30 = -pkin(8) + t35;
t56 = rSges(6,3) - t30;
t55 = rSges(4,3) + qJ(3);
t54 = g(1) * t68;
t52 = -m(4) - m(5) - m(6) - m(7);
t49 = t25 * rSges(5,1) - t23 * rSges(5,2);
t47 = -rSges(6,1) * t20 - rSges(6,2) * t21;
t46 = g(2) * t29 - t54;
t45 = rSges(4,1) * t34 - rSges(4,2) * sin(pkin(11)) + pkin(2);
t44 = t22 + t49;
t42 = t43 + (-t65 + t66) * t21;
t41 = t73 * (t20 * t65 + t21 * t77);
t4 = t21 * t58 + t61;
t3 = -t21 * t59 + t60;
t2 = -t21 * t60 + t59;
t1 = t21 * t61 + t58;
t5 = [-m(2) * (g(1) * (-t37 * rSges(2,1) - t39 * rSges(2,2)) + g(2) * (t39 * rSges(2,1) - t37 * rSges(2,2))) - m(3) * (g(1) * (-t24 * rSges(3,1) - t26 * rSges(3,2) - t68) + g(2) * (t26 * rSges(3,1) - t24 * rSges(3,2) + t29)) - m(4) * ((g(1) * t55 + g(2) * t45) * t26 + (-g(1) * t45 + g(2) * t55) * t24 + t46) - m(5) * ((g(1) * t57 + g(2) * t44) * t26 + (-g(1) * t44 + g(2) * t57) * t24 + t46) - m(6) * (-t54 + g(2) * t67 + (g(1) * t56 + g(2) * t74) * t26 + (g(1) * (-t11 - t74) + g(2) * t56) * t24) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t68) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t67) + (-g(1) * t30 + g(2) * t43) * t26 + (g(1) * (-t11 - t43) - g(2) * t30) * t24) (-m(3) + t52) * g(3), t52 * (g(1) * t24 - g(2) * t26) -m(7) * t41 + (-m(5) * t49 - m(6) * (t19 + t74) - m(7) * (t19 + t42)) * g(3) + t73 * (-m(5) * (-rSges(5,1) * t23 - rSges(5,2) * t25) - m(6) * (t47 - t71) - m(7) * (-t71 + t72)) -m(6) * (g(3) * t74 + t73 * t47) - m(7) * (g(3) * t42 + t73 * t72 + t41) -m(7) * (g(1) * (t3 * rSges(7,1) - t4 * rSges(7,2)) + g(2) * (-t1 * rSges(7,1) + t2 * rSges(7,2)) + g(3) * (-rSges(7,1) * t36 - rSges(7,2) * t38) * t20)];
taug  = t5(:);
