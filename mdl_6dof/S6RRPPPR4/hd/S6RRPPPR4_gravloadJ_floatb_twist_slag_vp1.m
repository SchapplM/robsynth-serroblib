% Calculate Gravitation load on the joints for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:23
% EndTime: 2019-03-09 08:17:25
% DurationCPUTime: 0.71s
% Computational Cost: add. (266->132), mult. (575->179), div. (0->0), fcn. (578->8), ass. (0->58)
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t51 = g(1) * t37 + g(2) * t34;
t78 = -m(6) - m(7);
t30 = sin(pkin(9));
t77 = pkin(4) * t30;
t76 = g(1) * t34;
t36 = cos(qJ(2));
t73 = g(3) * t36;
t26 = t36 * pkin(2);
t71 = -pkin(8) - rSges(7,3);
t70 = rSges(4,2) * t36;
t69 = t30 * t34;
t33 = sin(qJ(2));
t68 = t33 * t37;
t31 = cos(pkin(9));
t67 = t34 * t31;
t66 = t36 * t37;
t65 = -pkin(2) - qJ(4);
t22 = t33 * qJ(3);
t64 = t22 + t26;
t63 = t37 * pkin(1) + t34 * pkin(7);
t27 = t37 * pkin(7);
t62 = t37 * pkin(3) + t27;
t61 = qJ(3) * t36;
t23 = t36 * qJ(4);
t60 = rSges(6,3) + qJ(5);
t59 = -m(5) + t78;
t58 = t36 * t77;
t57 = -rSges(6,2) + t65;
t56 = -rSges(5,3) + t65;
t55 = t23 + t64;
t54 = -pkin(1) - t22;
t53 = t65 - t71;
t52 = pkin(2) * t66 + t37 * t22 + t63;
t50 = rSges(3,1) * t36 - rSges(3,2) * t33;
t48 = rSges(5,1) * t30 + rSges(5,2) * t31;
t10 = t30 * t37 + t33 * t67;
t11 = t31 * t37 - t33 * t69;
t32 = sin(qJ(6));
t35 = cos(qJ(6));
t47 = t10 * t35 - t11 * t32;
t46 = t10 * t32 + t11 * t35;
t45 = t30 * t35 - t31 * t32;
t44 = t30 * t32 + t31 * t35;
t43 = t11 * pkin(4) + t10 * qJ(5) + t62;
t42 = t34 * pkin(3) + t37 * t23 + t52;
t9 = t30 * t68 + t67;
t41 = t9 * pkin(4) + t42;
t40 = rSges(6,1) * t30 - t60 * t31;
t17 = t34 * t61;
t20 = t37 * t61;
t39 = g(1) * (t37 * t58 + t20) + g(2) * (t34 * t58 + t17) + g(3) * (t33 * t77 + t55);
t38 = t45 * rSges(7,1) - t44 * rSges(7,2) + t30 * pkin(5) - t31 * qJ(5);
t8 = -t31 * t68 + t69;
t3 = t32 * t8 + t35 * t9;
t2 = -t32 * t9 + t35 * t8;
t1 = [-m(2) * (g(1) * (-t34 * rSges(2,1) - rSges(2,2) * t37) + g(2) * (rSges(2,1) * t37 - t34 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t37 + t27) + g(2) * (rSges(3,1) * t66 - rSges(3,2) * t68 + t63) + (g(1) * (-pkin(1) - t50) + g(2) * rSges(3,3)) * t34) - m(4) * (g(1) * (rSges(4,1) * t37 + t27) + g(2) * (-rSges(4,2) * t66 + rSges(4,3) * t68 + t52) + (g(1) * (-rSges(4,3) * t33 - t26 + t54 + t70) + g(2) * rSges(4,1)) * t34) - m(5) * (g(1) * (rSges(5,1) * t11 - rSges(5,2) * t10 + t62) + g(2) * (t9 * rSges(5,1) - t8 * rSges(5,2) + rSges(5,3) * t66 + t42) + (t56 * t36 + t54) * t76) - m(6) * (g(1) * (rSges(6,1) * t11 + rSges(6,3) * t10 + t43) + g(2) * (t9 * rSges(6,1) + rSges(6,2) * t66 + t60 * t8 + t41) + (t57 * t36 + t54) * t76) - m(7) * (g(1) * (t46 * rSges(7,1) + t47 * rSges(7,2) + t11 * pkin(5) + t43) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t9 * pkin(5) + t8 * qJ(5) + t71 * t66 + t41) + (t53 * t36 + t54) * t76) -m(3) * (g(3) * t50 + t51 * (-rSges(3,1) * t33 - rSges(3,2) * t36)) - m(4) * (g(1) * (rSges(4,3) * t66 + t20) + g(2) * (rSges(4,3) * t34 * t36 + t17) + g(3) * (t64 - t70) + (g(3) * rSges(4,3) + t51 * (rSges(4,2) - pkin(2))) * t33) - m(5) * (g(1) * t20 + g(2) * t17 + g(3) * t55 + (g(3) * rSges(5,3) + t51 * t48) * t36 + (g(3) * t48 + t51 * t56) * t33) - m(6) * ((g(3) * rSges(6,2) + t51 * t40) * t36 + (g(3) * t40 + t51 * t57) * t33 + t39) - m(7) * ((g(3) * t71 + t51 * t38) * t36 + (g(3) * t38 + t51 * t53) * t33 + t39) (-m(4) + t59) * (t51 * t33 - t73) t59 * (g(3) * t33 + t51 * t36) t78 * (g(1) * t8 - g(2) * t10 + t31 * t73) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * (-t47 * rSges(7,1) + t46 * rSges(7,2)) + (t44 * rSges(7,1) + t45 * rSges(7,2)) * t73)];
taug  = t1(:);
