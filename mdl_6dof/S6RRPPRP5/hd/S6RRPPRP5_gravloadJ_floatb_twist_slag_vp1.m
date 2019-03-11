% Calculate Gravitation load on the joints for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:41:41
% EndTime: 2019-03-09 08:41:43
% DurationCPUTime: 0.83s
% Computational Cost: add. (332->135), mult. (515->184), div. (0->0), fcn. (491->8), ass. (0->58)
t35 = cos(qJ(2));
t33 = sin(qJ(2));
t24 = t33 * qJ(3);
t57 = t35 * pkin(2) + t24;
t77 = -rSges(4,2) * t35 + t57;
t71 = rSges(7,1) + pkin(5);
t36 = cos(qJ(1));
t34 = sin(qJ(1));
t73 = g(2) * t34;
t46 = g(1) * t36 + t73;
t54 = rSges(5,3) + qJ(4);
t53 = rSges(7,3) + qJ(6);
t30 = sin(pkin(9));
t76 = pkin(4) * t30;
t75 = g(1) * t34;
t72 = g(3) * t35;
t69 = -rSges(7,2) - pkin(2);
t68 = -rSges(6,3) - pkin(2);
t31 = cos(pkin(9));
t66 = rSges(5,2) * t31;
t65 = t30 * t33;
t64 = t33 * t34;
t63 = t33 * t36;
t29 = pkin(9) + qJ(5);
t23 = cos(t29);
t62 = t34 * t23;
t61 = t34 * t35;
t60 = t35 * t36;
t32 = -pkin(8) - qJ(4);
t59 = rSges(7,2) - t32;
t58 = rSges(6,3) - t32;
t56 = pkin(1) * t36 + pkin(7) * t34;
t55 = qJ(3) * t35;
t52 = -m(5) - m(6) - m(7);
t51 = t35 * t76;
t16 = pkin(4) * t65;
t50 = -pkin(2) - t54;
t21 = pkin(4) * t31 + pkin(3);
t27 = t36 * pkin(7);
t49 = t21 * t36 + t32 * t61 + t27;
t48 = pkin(2) * t60 + t24 * t36 + t56;
t47 = g(1) * t50;
t45 = rSges(3,1) * t35 - rSges(3,2) * t33;
t43 = rSges(5,1) * t30 + t66;
t22 = sin(t29);
t42 = rSges(6,1) * t22 + rSges(6,2) * t23;
t41 = t16 * t36 + t21 * t34 + t48;
t40 = rSges(5,1) * t31 - rSges(5,2) * t30 + pkin(3);
t39 = (-qJ(3) - t76) * t33 - pkin(1);
t38 = t22 * t71 - t23 * t53;
t17 = t34 * t55;
t19 = t36 * t55;
t37 = g(1) * (t32 * t63 + t36 * t51 + t19) + g(2) * (t32 * t64 + t34 * t51 + t17) + g(3) * (t16 + t57);
t4 = -t22 * t64 + t23 * t36;
t3 = t22 * t36 + t33 * t62;
t2 = t22 * t63 + t62;
t1 = t22 * t34 - t23 * t63;
t5 = [-m(2) * (g(1) * (-rSges(2,1) * t34 - rSges(2,2) * t36) + g(2) * (rSges(2,1) * t36 - rSges(2,2) * t34)) - m(3) * (g(1) * (rSges(3,3) * t36 + t27) + g(2) * (rSges(3,1) * t60 - rSges(3,2) * t63 + t56) + (g(1) * (-pkin(1) - t45) + g(2) * rSges(3,3)) * t34) - m(4) * (g(1) * (rSges(4,1) * t36 + t27) + g(2) * (-rSges(4,2) * t60 + rSges(4,3) * t63 + t48) + (g(1) * (-rSges(4,3) * t33 - pkin(1) - t77) + g(2) * rSges(4,1)) * t34) - m(5) * (g(1) * t27 + g(2) * t48 + (g(1) * t40 + g(2) * (rSges(5,1) * t65 + t33 * t66 + t35 * t54)) * t36 + (g(2) * t40 + t35 * t47 + (-pkin(1) + (-qJ(3) - t43) * t33) * g(1)) * t34) - m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t49) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t58 * t60 + t41) + (t35 * t68 + t39) * t75) - m(7) * (g(1) * (t3 * t53 + t4 * t71 + t49) + g(2) * (t1 * t53 + t2 * t71 + t59 * t60 + t41) + (t35 * t69 + t39) * t75) -m(3) * (g(3) * t45 + t46 * (-rSges(3,1) * t33 - rSges(3,2) * t35)) - m(4) * (g(1) * (rSges(4,3) * t60 + t19) + g(2) * (rSges(4,3) * t61 + t17) + g(3) * t77 + (g(3) * rSges(4,3) + t46 * (rSges(4,2) - pkin(2))) * t33) - m(5) * (g(1) * t19 + g(2) * t17 + g(3) * t57 + (g(3) * t54 + t43 * t46) * t35 + (g(3) * t43 + t36 * t47 + t50 * t73) * t33) - m(6) * ((g(3) * t58 + t42 * t46) * t35 + (g(3) * t42 + t46 * t68) * t33 + t37) - m(7) * ((g(3) * t38 + t46 * t69) * t33 + (g(3) * t59 + t38 * t46) * t35 + t37) (-m(4) + t52) * (t33 * t46 - t72) t52 * (g(3) * t33 + t35 * t46) -m(6) * (g(1) * (-rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (rSges(6,1) * t3 + rSges(6,2) * t4)) - m(7) * (g(1) * (-t1 * t71 + t2 * t53) + g(2) * (t3 * t71 - t4 * t53)) + (-m(6) * (-rSges(6,1) * t23 + rSges(6,2) * t22) - m(7) * (-t22 * t53 - t23 * t71)) * t72, -m(7) * (g(1) * t1 - g(2) * t3 + t23 * t72)];
taug  = t5(:);
