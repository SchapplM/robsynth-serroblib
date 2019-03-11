% Calculate Gravitation load on the joints for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:07:49
% EndTime: 2019-03-09 03:07:50
% DurationCPUTime: 0.65s
% Computational Cost: add. (444->107), mult. (417->147), div. (0->0), fcn. (396->10), ass. (0->49)
t27 = sin(qJ(3));
t29 = cos(qJ(3));
t24 = sin(pkin(10));
t25 = cos(pkin(10));
t35 = rSges(5,1) * t25 - rSges(5,2) * t24 + pkin(3);
t50 = rSges(5,3) + qJ(4);
t65 = t50 * t27 + t35 * t29;
t56 = rSges(7,1) + pkin(5);
t23 = qJ(1) + pkin(9);
t20 = cos(t23);
t18 = sin(t23);
t59 = g(2) * t18;
t64 = g(1) * t20 + t59;
t49 = rSges(7,3) + qJ(6);
t22 = pkin(10) + qJ(5);
t17 = sin(t22);
t19 = cos(t22);
t63 = t49 * t17 + t56 * t19;
t62 = pkin(4) * t24;
t61 = g(1) * t18;
t58 = g(3) * t27;
t28 = sin(qJ(1));
t57 = t28 * pkin(1);
t55 = rSges(4,2) * t27;
t54 = t18 * t29;
t53 = t20 * t29;
t16 = pkin(4) * t25 + pkin(3);
t12 = t29 * t16;
t26 = -pkin(8) - qJ(4);
t52 = rSges(7,2) - t26;
t51 = rSges(6,3) - t26;
t48 = -m(5) - m(6) - m(7);
t30 = cos(qJ(1));
t21 = t30 * pkin(1);
t47 = t20 * pkin(2) + t18 * pkin(7) + t21;
t46 = t20 * pkin(7) - t57;
t45 = -pkin(2) - t12;
t44 = t20 * t52;
t43 = t51 * t27;
t41 = t20 * t12 + t18 * t62 + t47;
t40 = rSges(4,1) * t29 - t55;
t38 = rSges(5,1) * t24 + rSges(5,2) * t25;
t37 = rSges(6,1) * t19 - rSges(6,2) * t17;
t36 = t18 * t27 * t26 + t20 * t62 + t46;
t4 = t17 * t18 + t19 * t53;
t3 = t17 * t53 - t18 * t19;
t2 = -t20 * t17 + t19 * t54;
t1 = t17 * t54 + t19 * t20;
t5 = [-m(2) * (g(1) * (-t28 * rSges(2,1) - rSges(2,2) * t30) + g(2) * (rSges(2,1) * t30 - t28 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t18 - rSges(3,2) * t20 - t57) + g(2) * (rSges(3,1) * t20 - rSges(3,2) * t18 + t21)) - m(4) * (g(1) * (rSges(4,3) * t20 + t46) + g(2) * (rSges(4,1) * t53 - t20 * t55 + t47) + (g(1) * (-pkin(2) - t40) + g(2) * rSges(4,3)) * t18) - m(5) * (g(1) * t46 + g(2) * t47 + (g(1) * t38 + g(2) * t65) * t20 + (g(2) * t38 + (-pkin(2) - t65) * g(1)) * t18) - m(6) * (g(1) * (-rSges(6,1) * t2 + rSges(6,2) * t1 + t36) + g(2) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t20 * t43 + t41) + (-t27 * rSges(6,3) + t45) * t61) - m(7) * (g(1) * (-t49 * t1 - t2 * t56 + t36) + g(2) * (t27 * t44 + t49 * t3 + t56 * t4 + t41) + (-rSges(7,2) * t27 + t45) * t61) (-m(3) - m(4) + t48) * g(3), -m(4) * (g(3) * t40 + t64 * (-rSges(4,1) * t27 - rSges(4,2) * t29)) - m(5) * (g(3) * t65 + t64 * (-t35 * t27 + t50 * t29)) - m(6) * (g(3) * (t37 * t29 + t12 + t43) + t64 * (t51 * t29 + (-t16 - t37) * t27)) - m(7) * (g(3) * t12 + (g(1) * t44 + g(3) * t63 + t52 * t59) * t29 + (g(3) * t52 + t64 * (-t16 - t63)) * t27) t48 * (-g(3) * t29 + t27 * t64) -m(6) * (g(1) * (-rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2)) - m(7) * (g(1) * (-t56 * t3 + t49 * t4) + g(2) * (-t56 * t1 + t49 * t2)) + (-m(6) * (-rSges(6,1) * t17 - rSges(6,2) * t19) - m(7) * (-t56 * t17 + t49 * t19)) * t58, -m(7) * (g(1) * t3 + g(2) * t1 + t17 * t58)];
taug  = t5(:);
