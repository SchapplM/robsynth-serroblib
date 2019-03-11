% Calculate Gravitation load on the joints for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:49:49
% EndTime: 2019-03-09 09:49:52
% DurationCPUTime: 0.78s
% Computational Cost: add. (431->146), mult. (558->187), div. (0->0), fcn. (549->8), ass. (0->58)
t27 = qJ(2) + pkin(9);
t24 = sin(t27);
t70 = g(3) * t24;
t68 = rSges(7,1) + pkin(5);
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t77 = g(1) * t34 + g(2) * t31;
t53 = -rSges(7,3) - qJ(6);
t76 = -m(6) - m(7);
t30 = sin(qJ(2));
t75 = pkin(2) * t30;
t28 = -qJ(3) - pkin(7);
t73 = g(2) * t28;
t32 = cos(qJ(4));
t71 = t32 * qJ(5) * t70;
t19 = t24 * pkin(8);
t25 = cos(t27);
t20 = t25 * pkin(3);
t69 = -rSges(6,1) - pkin(4);
t67 = rSges(3,3) + pkin(7);
t66 = t24 * t34;
t29 = sin(qJ(4));
t65 = t25 * t29;
t64 = t25 * t31;
t63 = t25 * t32;
t62 = t25 * t34;
t61 = t28 * t34;
t60 = t31 * t29;
t59 = t31 * t32;
t58 = t32 * t34;
t57 = t34 * t29;
t56 = rSges(4,3) - t28;
t55 = rSges(7,2) + qJ(5);
t54 = rSges(6,3) + qJ(5);
t52 = -pkin(4) - t68;
t33 = cos(qJ(2));
t26 = t33 * pkin(2);
t23 = t26 + pkin(1);
t18 = t34 * t23;
t51 = pkin(3) * t62 + pkin(8) * t66 + t18;
t50 = t19 + t20 + t26;
t49 = -t23 - t20;
t48 = pkin(8) * t64 - t31 * t75;
t47 = pkin(8) * t62 - t34 * t75;
t46 = pkin(4) * t63 + qJ(5) * t65 + t50;
t45 = rSges(3,1) * t33 - rSges(3,2) * t30;
t43 = rSges(4,1) * t25 - rSges(4,2) * t24;
t42 = t49 - t19;
t41 = pkin(1) + t45;
t6 = t25 * t60 + t58;
t7 = t25 * t59 - t57;
t39 = -t7 * pkin(4) - t6 * qJ(5) - t61;
t8 = t25 * t57 - t59;
t9 = t25 * t58 + t60;
t38 = t9 * pkin(4) + t8 * qJ(5) + t51;
t4 = t8 * pkin(4);
t2 = t6 * pkin(4);
t1 = [-m(2) * (g(1) * (-t31 * rSges(2,1) - rSges(2,2) * t34) + g(2) * (rSges(2,1) * t34 - t31 * rSges(2,2))) - m(3) * ((g(1) * t67 + g(2) * t41) * t34 + (-g(1) * t41 + g(2) * t67) * t31) - m(4) * (g(2) * t18 + (g(1) * t56 + g(2) * t43) * t34 + (g(1) * (-t23 - t43) + g(2) * t56) * t31) - m(5) * (g(1) * (-t7 * rSges(5,1) + t6 * rSges(5,2) - t61) + g(2) * (t9 * rSges(5,1) - t8 * rSges(5,2) + rSges(5,3) * t66 + t51) + (g(1) * (-rSges(5,3) * t24 + t42) - t73) * t31) - m(6) * (g(1) * (-t7 * rSges(6,1) - t6 * rSges(6,3) + t39) + g(2) * (t9 * rSges(6,1) + rSges(6,2) * t66 + t8 * rSges(6,3) + t38) + (g(1) * (-rSges(6,2) * t24 + t42) - t73) * t31) - m(7) * (g(1) * (-t6 * rSges(7,2) - t68 * t7 + t39) + g(2) * (t8 * rSges(7,2) + t53 * t66 + t68 * t9 + t38) + (-t73 + (t49 + (-pkin(8) - t53) * t24) * g(1)) * t31) -m(3) * (g(3) * t45 + t77 * (-rSges(3,1) * t30 - rSges(3,2) * t33)) - m(4) * (g(3) * (t26 + t43) + t77 * (-rSges(4,1) * t24 - rSges(4,2) * t25 - t75)) - m(5) * (g(1) * (rSges(5,3) * t62 + t47) + g(2) * (rSges(5,3) * t64 + t48) + g(3) * (rSges(5,1) * t63 - rSges(5,2) * t65 + t50) + (g(3) * rSges(5,3) + t77 * (-rSges(5,1) * t32 + rSges(5,2) * t29 - pkin(3))) * t24) - m(6) * (g(1) * (rSges(6,2) * t62 + t47) + g(2) * (rSges(6,2) * t64 + t48) + g(3) * (rSges(6,1) * t63 + rSges(6,3) * t65 + t46) + (g(3) * rSges(6,2) + t77 * (-t54 * t29 + t69 * t32 - pkin(3))) * t24) - m(7) * (g(1) * t47 + g(2) * t48 + g(3) * t46 + (g(3) * (rSges(7,2) * t29 + t68 * t32) + t77 * t53) * t25 + (g(3) * t53 + t77 * (-t55 * t29 + t52 * t32 - pkin(3))) * t24) (-m(4) - m(5) + t76) * (g(1) * t31 - g(2) * t34) -m(5) * (g(1) * (-rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (-rSges(5,1) * t6 - rSges(5,2) * t7)) - m(6) * (g(1) * (-rSges(6,1) * t8 + t54 * t9 - t4) + g(2) * (-rSges(6,1) * t6 + t54 * t7 - t2) + t71) - m(7) * (g(1) * (t55 * t9 - t68 * t8 - t4) + g(2) * (t55 * t7 - t68 * t6 - t2) + t71) + ((m(5) * rSges(5,2) - m(6) * rSges(6,3) - m(7) * rSges(7,2)) * t32 + (m(5) * rSges(5,1) - m(6) * t69 - m(7) * t52) * t29) * t70, t76 * (g(1) * t8 + g(2) * t6 + t29 * t70) -m(7) * (g(3) * t25 - t77 * t24)];
taug  = t1(:);
