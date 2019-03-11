% Calculate Gravitation load on the joints for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:06:05
% EndTime: 2019-03-09 08:06:07
% DurationCPUTime: 0.80s
% Computational Cost: add. (411->138), mult. (540->182), div. (0->0), fcn. (543->10), ass. (0->60)
t27 = qJ(2) + pkin(9);
t25 = cos(t27);
t32 = sin(qJ(2));
t73 = t32 * pkin(2);
t80 = qJ(4) * t25 - t73;
t71 = -pkin(8) - rSges(7,3);
t33 = sin(qJ(1));
t36 = cos(qJ(1));
t79 = g(1) * t36 + g(2) * t33;
t78 = -m(6) - m(7);
t30 = -qJ(3) - pkin(7);
t76 = g(2) * t30;
t24 = sin(t27);
t74 = g(3) * t24;
t20 = t25 * pkin(3);
t72 = rSges(3,3) + pkin(7);
t70 = t24 * t36;
t28 = sin(pkin(10));
t69 = t25 * t28;
t29 = cos(pkin(10));
t68 = t25 * t29;
t67 = t25 * t33;
t66 = t25 * t36;
t65 = t29 * t36;
t64 = t33 * t28;
t63 = t33 * t29;
t62 = t36 * t28;
t61 = t36 * t30;
t60 = rSges(4,3) - t30;
t19 = t24 * qJ(4);
t58 = -m(5) + t78;
t35 = cos(qJ(2));
t26 = t35 * pkin(2);
t23 = t26 + pkin(1);
t18 = t36 * t23;
t57 = pkin(3) * t66 + t36 * t19 + t18;
t56 = t19 + t20 + t26;
t55 = -t23 - t20;
t54 = t80 * t33;
t53 = t80 * t36;
t52 = pkin(4) * t68 + qJ(5) * t69 + t56;
t31 = sin(qJ(6));
t34 = cos(qJ(6));
t7 = t25 * t64 + t65;
t8 = t25 * t63 - t62;
t51 = t31 * t8 - t34 * t7;
t50 = -t31 * t7 - t34 * t8;
t49 = rSges(3,1) * t35 - rSges(3,2) * t32;
t47 = rSges(4,1) * t25 - rSges(4,2) * t24;
t46 = t28 * t34 - t29 * t31;
t45 = t28 * t31 + t29 * t34;
t44 = pkin(1) + t49;
t42 = -t8 * pkin(4) - t7 * qJ(5) - t61;
t41 = t55 - t19;
t10 = t25 * t65 + t64;
t9 = t25 * t62 - t63;
t40 = t10 * pkin(4) + t9 * qJ(5) + t57;
t3 = t10 * t34 + t31 * t9;
t2 = -t10 * t31 + t34 * t9;
t1 = [-m(2) * (g(1) * (-t33 * rSges(2,1) - rSges(2,2) * t36) + g(2) * (rSges(2,1) * t36 - t33 * rSges(2,2))) - m(3) * ((g(1) * t72 + g(2) * t44) * t36 + (-g(1) * t44 + g(2) * t72) * t33) - m(4) * (g(2) * t18 + (g(1) * t60 + g(2) * t47) * t36 + (g(1) * (-t23 - t47) + g(2) * t60) * t33) - m(5) * (g(1) * (-t8 * rSges(5,1) + t7 * rSges(5,2) - t61) + g(2) * (t10 * rSges(5,1) - t9 * rSges(5,2) + rSges(5,3) * t70 + t57) + (g(1) * (-rSges(5,3) * t24 + t41) - t76) * t33) - m(6) * (g(1) * (-t8 * rSges(6,1) - t7 * rSges(6,3) + t42) + g(2) * (t10 * rSges(6,1) + rSges(6,2) * t70 + t9 * rSges(6,3) + t40) + (g(1) * (-rSges(6,2) * t24 + t41) - t76) * t33) - m(7) * (g(1) * (t50 * rSges(7,1) + t51 * rSges(7,2) - t8 * pkin(5) + t42) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t10 * pkin(5) + t70 * t71 + t40) + (-t76 + (t55 + (-qJ(4) - t71) * t24) * g(1)) * t33) -m(3) * (g(3) * t49 + t79 * (-rSges(3,1) * t32 - rSges(3,2) * t35)) - m(4) * (g(3) * (t26 + t47) + t79 * (-rSges(4,1) * t24 - rSges(4,2) * t25 - t73)) - m(5) * (g(1) * (rSges(5,3) * t66 + t53) + g(2) * (rSges(5,3) * t67 + t54) + g(3) * (rSges(5,1) * t68 - rSges(5,2) * t69 + t56) + (g(3) * rSges(5,3) + t79 * (-rSges(5,1) * t29 + rSges(5,2) * t28 - pkin(3))) * t24) - m(6) * (g(1) * (rSges(6,2) * t66 + t53) + g(2) * (rSges(6,2) * t67 + t54) + g(3) * (rSges(6,1) * t68 + rSges(6,3) * t69 + t52) + (g(3) * rSges(6,2) + t79 * (-pkin(3) + (-rSges(6,1) - pkin(4)) * t29 + (-rSges(6,3) - qJ(5)) * t28)) * t24) - m(7) * (g(1) * t53 + g(2) * t54 + g(3) * t52 + (g(3) * (t45 * rSges(7,1) + t46 * rSges(7,2) + t29 * pkin(5)) + t79 * t71) * t25 + (g(3) * t71 + t79 * (-pkin(3) + (-t31 * rSges(7,1) - t34 * rSges(7,2) - qJ(5)) * t28 + (-t34 * rSges(7,1) + t31 * rSges(7,2) - pkin(4) - pkin(5)) * t29)) * t24) (-m(4) + t58) * (g(1) * t33 - g(2) * t36) t58 * (-g(3) * t25 + t24 * t79) t78 * (g(1) * t9 + g(2) * t7 + t28 * t74) -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * (-t51 * rSges(7,1) + t50 * rSges(7,2)) + (t46 * rSges(7,1) - t45 * rSges(7,2)) * t74)];
taug  = t1(:);
