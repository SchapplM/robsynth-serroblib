% Calculate Gravitation load on the joints for
% S6RPRPRP5
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
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:14:32
% EndTime: 2019-03-09 03:14:34
% DurationCPUTime: 0.67s
% Computational Cost: add. (434->113), mult. (435->154), div. (0->0), fcn. (414->10), ass. (0->54)
t56 = rSges(7,1) + pkin(5);
t29 = sin(qJ(1));
t58 = g(2) * t29;
t30 = cos(qJ(1));
t60 = g(1) * t30;
t65 = t58 + t60;
t46 = rSges(5,3) + qJ(4);
t45 = rSges(7,3) + qJ(6);
t21 = pkin(10) + qJ(5);
t17 = sin(t21);
t19 = cos(t21);
t64 = t45 * t17 + t56 * t19;
t26 = cos(pkin(9));
t16 = pkin(2) * t26 + pkin(1);
t9 = t30 * t16;
t63 = g(2) * t9;
t25 = cos(pkin(10));
t15 = pkin(4) * t25 + pkin(3);
t22 = pkin(9) + qJ(3);
t20 = cos(t22);
t7 = t20 * t15;
t62 = g(3) * t7;
t23 = sin(pkin(10));
t61 = pkin(4) * t23;
t28 = -pkin(7) - qJ(2);
t59 = g(2) * t28;
t18 = sin(t22);
t57 = g(3) * t18;
t55 = t18 * t30;
t54 = t20 * t30;
t53 = t29 * t17;
t52 = t29 * t19;
t51 = t30 * t17;
t27 = -pkin(8) - qJ(4);
t50 = rSges(7,2) - t27;
t49 = rSges(4,3) - t28;
t48 = rSges(6,3) - t27;
t47 = rSges(3,3) + qJ(2);
t44 = -m(5) - m(6) - m(7);
t43 = g(2) * t46;
t42 = -t16 - t7;
t41 = t29 * t18 * t27 + (-t28 + t61) * t30;
t40 = rSges(4,1) * t20 - rSges(4,2) * t18;
t38 = rSges(6,1) * t19 - rSges(6,2) * t17;
t37 = rSges(3,1) * t26 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t36 = rSges(5,1) * t25 - rSges(5,2) * t23 + pkin(3);
t35 = rSges(5,1) * t23 + rSges(5,2) * t25 - t28;
t33 = t15 * t54 - t27 * t55 + t29 * t61 + t9;
t32 = g(1) * t36;
t5 = t19 * t54 + t53;
t4 = t20 * t51 - t52;
t3 = t20 * t52 - t51;
t2 = t19 * t30 + t20 * t53;
t1 = [-m(2) * (g(1) * (-t29 * rSges(2,1) - rSges(2,2) * t30) + g(2) * (rSges(2,1) * t30 - t29 * rSges(2,2))) - m(3) * ((g(1) * t47 + g(2) * t37) * t30 + (-g(1) * t37 + g(2) * t47) * t29) - m(4) * (t63 + (g(1) * t49 + g(2) * t40) * t30 + (g(1) * (-t16 - t40) + g(2) * t49) * t29) - m(5) * (t63 + (g(2) * t36 * t20 + g(1) * t35 + t18 * t43) * t30 + (g(1) * (-t46 * t18 - t16) + g(2) * t35 - t20 * t32) * t29) - m(6) * (g(1) * (-t3 * rSges(6,1) + t2 * rSges(6,2) + t41) + g(2) * (t5 * rSges(6,1) - t4 * rSges(6,2) + rSges(6,3) * t55 + t33) + (g(1) * (-rSges(6,3) * t18 + t42) - t59) * t29) - m(7) * (g(1) * (-t45 * t2 - t56 * t3 + t41) + g(2) * (rSges(7,2) * t55 + t45 * t4 + t56 * t5 + t33) + (g(1) * (-rSges(7,2) * t18 + t42) - t59) * t29) (-m(3) - m(4) + t44) * (g(1) * t29 - g(2) * t30) -m(4) * (g(3) * t40 + t65 * (-rSges(4,1) * t18 - rSges(4,2) * t20)) - m(5) * ((g(3) * t36 + t29 * t43 + t46 * t60) * t20 + (g(3) * t46 - t30 * t32 - t36 * t58) * t18) - m(6) * (t62 + (g(3) * t38 + t65 * t48) * t20 + (g(3) * t48 + t65 * (-t15 - t38)) * t18) - m(7) * (t62 + (g(3) * t64 + t65 * t50) * t20 + (g(3) * t50 + t65 * (-t15 - t64)) * t18) t44 * (-g(3) * t20 + t65 * t18) -m(6) * (g(1) * (-rSges(6,1) * t4 - rSges(6,2) * t5) + g(2) * (-rSges(6,1) * t2 - rSges(6,2) * t3)) - m(7) * (g(1) * (-t56 * t4 + t45 * t5) + g(2) * (-t56 * t2 + t45 * t3)) + (-m(6) * (-rSges(6,1) * t17 - rSges(6,2) * t19) - m(7) * (-t56 * t17 + t45 * t19)) * t57, -m(7) * (g(1) * t4 + g(2) * t2 + t17 * t57)];
taug  = t1(:);
