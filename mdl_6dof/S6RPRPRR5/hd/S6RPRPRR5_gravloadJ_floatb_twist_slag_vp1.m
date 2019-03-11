% Calculate Gravitation load on the joints for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:47:59
% EndTime: 2019-03-09 03:48:01
% DurationCPUTime: 0.72s
% Computational Cost: add. (455->119), mult. (547->159), div. (0->0), fcn. (563->10), ass. (0->58)
t32 = sin(qJ(6));
t34 = cos(qJ(6));
t33 = sin(qJ(1));
t28 = pkin(10) + qJ(3);
t26 = sin(t28);
t27 = cos(t28);
t65 = sin(qJ(5));
t66 = cos(qJ(5));
t9 = t26 * t66 - t27 * t65;
t4 = t9 * t33;
t35 = cos(qJ(1));
t50 = t35 * t65;
t51 = t35 * t66;
t7 = -t26 * t51 + t27 * t50;
t8 = t26 * t65 + t27 * t66;
t80 = (g(1) * t7 - g(2) * t4 + g(3) * t8) * (t34 * rSges(7,1) - t32 * rSges(7,2) + pkin(5));
t79 = g(3) * t9;
t78 = rSges(7,3) + pkin(9);
t31 = -pkin(7) - qJ(2);
t67 = -pkin(8) - t31;
t22 = t26 * qJ(4);
t60 = t27 * pkin(3) + t22;
t71 = g(1) * t35;
t77 = g(2) * t33 + t71;
t76 = t77 * t26;
t23 = t27 * pkin(4);
t64 = rSges(5,1) * t27;
t63 = t27 * t35;
t62 = rSges(5,2) - t31;
t61 = rSges(4,3) - t31;
t59 = qJ(4) * t27;
t58 = rSges(3,3) + qJ(2);
t57 = -m(5) - m(6) - m(7);
t56 = -rSges(6,3) + t67;
t30 = cos(pkin(10));
t25 = pkin(2) * t30 + pkin(1);
t17 = t35 * t25;
t55 = pkin(3) * t63 + t35 * t22 + t17;
t54 = t23 + t60;
t49 = pkin(4) * t63 + t55;
t5 = t8 * t33;
t48 = rSges(6,1) * t4 - rSges(6,2) * t5;
t6 = -t26 * t50 - t27 * t51;
t47 = -t7 * rSges(6,1) + t6 * rSges(6,2);
t46 = -rSges(6,1) * t8 - rSges(6,2) * t9;
t45 = t5 * t32 - t34 * t35;
t44 = -t32 * t35 - t5 * t34;
t43 = rSges(4,1) * t27 - rSges(4,2) * t26;
t41 = rSges(5,3) * t26 + t64;
t40 = rSges(3,1) * t30 - rSges(3,2) * sin(pkin(10)) + pkin(1);
t38 = -t25 - t60;
t37 = g(1) * (t38 - t23);
t36 = (-pkin(3) - pkin(4)) * t76;
t16 = t35 * t59;
t14 = t33 * t59;
t2 = -t32 * t33 - t34 * t6;
t1 = t32 * t6 - t33 * t34;
t3 = [-m(2) * (g(1) * (-t33 * rSges(2,1) - rSges(2,2) * t35) + g(2) * (rSges(2,1) * t35 - t33 * rSges(2,2))) - m(3) * ((g(1) * t58 + g(2) * t40) * t35 + (-g(1) * t40 + g(2) * t58) * t33) - m(4) * (g(2) * t17 + (g(1) * t61 + g(2) * t43) * t35 + (g(1) * (-t25 - t43) + g(2) * t61) * t33) - m(5) * (g(2) * t55 + (g(1) * t62 + g(2) * t41) * t35 + (g(1) * (t38 - t41) + g(2) * t62) * t33) - m(6) * (g(1) * (-t5 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-rSges(6,1) * t6 - rSges(6,2) * t7 + t49) + t56 * t71 + (g(2) * t56 + t37) * t33) - m(7) * (g(1) * (t44 * rSges(7,1) + t45 * rSges(7,2) - t5 * pkin(5) + t67 * t35 + t78 * t4) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 - pkin(5) * t6 + t78 * t7 + t49) + (g(2) * t67 + t37) * t33) (-m(3) - m(4) + t57) * (g(1) * t33 - g(2) * t35) -m(4) * (g(3) * t43 + t77 * (-rSges(4,1) * t26 - rSges(4,2) * t27)) - m(5) * (g(1) * (rSges(5,3) * t63 + t16) + g(2) * (rSges(5,3) * t27 * t33 + t14) + g(3) * (t60 + t64) + (g(3) * rSges(5,3) + t77 * (-rSges(5,1) - pkin(3))) * t26) - m(6) * (g(1) * (t16 - t47) + g(2) * (t14 - t48) + g(3) * (-t46 + t54) + t36) + (-g(1) * (t78 * t6 + t16) - g(2) * (-t78 * t5 + t14) - g(3) * (-t78 * t9 + t54) - t36 - t80) * m(7), t57 * (-g(3) * t27 + t76) -m(6) * (g(1) * t47 + g(2) * t48 + g(3) * t46) - m(7) * ((-g(1) * t6 + g(2) * t5 + t79) * t78 - t80) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t45 * rSges(7,1) + t44 * rSges(7,2)) + (-t32 * rSges(7,1) - t34 * rSges(7,2)) * t79)];
taug  = t3(:);
