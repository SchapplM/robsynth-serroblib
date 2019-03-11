% Calculate Gravitation load on the joints for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:27:35
% EndTime: 2019-03-09 03:27:37
% DurationCPUTime: 0.77s
% Computational Cost: add. (299->112), mult. (429->153), div. (0->0), fcn. (408->8), ass. (0->50)
t24 = cos(pkin(9));
t15 = pkin(4) * t24 + pkin(3);
t22 = pkin(9) + qJ(5);
t16 = sin(t22);
t17 = cos(t22);
t45 = rSges(7,3) + qJ(6);
t56 = rSges(7,1) + pkin(5);
t64 = t45 * t16 + t56 * t17;
t68 = -t15 - t64;
t26 = sin(qJ(3));
t28 = cos(qJ(3));
t23 = sin(pkin(9));
t32 = rSges(5,1) * t24 - rSges(5,2) * t23 + pkin(3);
t46 = rSges(5,3) + qJ(4);
t67 = -t32 * t26 + t46 * t28;
t29 = cos(qJ(1));
t58 = g(2) * t29;
t27 = sin(qJ(1));
t60 = g(1) * t27;
t39 = -t58 + t60;
t62 = -pkin(1) - pkin(7);
t61 = pkin(4) * t23;
t59 = g(2) * t28;
t57 = g(3) * t28;
t55 = t26 * t29;
t54 = t27 * t16;
t53 = t27 * t17;
t52 = t27 * t28;
t51 = t28 * t29;
t50 = t29 * t17;
t25 = -pkin(8) - qJ(4);
t49 = rSges(7,2) - t25;
t48 = rSges(6,3) - t25;
t47 = t29 * pkin(1) + t27 * qJ(2);
t44 = -m(5) - m(6) - m(7);
t19 = t29 * qJ(2);
t43 = t15 * t55 + t25 * t51 + t19;
t42 = t29 * pkin(7) + t47;
t40 = g(1) * t15 * t52 + g(2) * t25 * t55;
t37 = rSges(4,1) * t26 + rSges(4,2) * t28;
t36 = rSges(5,1) * t23 + rSges(5,2) * t24;
t35 = rSges(6,1) * t17 - rSges(6,2) * t16;
t34 = g(1) * (-t61 + t62);
t33 = t27 * t26 * t15 + t25 * t52 + t29 * t61 + t42;
t31 = -t15 - t35;
t4 = t26 * t50 - t54;
t3 = t16 * t55 + t53;
t2 = t16 * t29 + t26 * t53;
t1 = t26 * t54 - t50;
t5 = [-m(2) * (g(1) * (-t27 * rSges(2,1) - rSges(2,2) * t29) + g(2) * (rSges(2,1) * t29 - t27 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t29 + t19 + (rSges(3,2) - pkin(1)) * t27) + g(2) * (-rSges(3,2) * t29 + t27 * rSges(3,3) + t47)) - m(4) * (g(1) * (rSges(4,1) * t55 + rSges(4,2) * t51 + t19) + g(2) * (rSges(4,3) * t29 + t42) + (g(1) * (-rSges(4,3) + t62) + g(2) * t37) * t27) - m(5) * (g(1) * t19 + g(2) * t42 + (-g(1) * t67 + g(2) * t36) * t29 + (g(1) * (-t36 + t62) - t67 * g(2)) * t27) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) - rSges(6,3) * t51 + t43) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t1 + t33) + (-rSges(6,3) * t59 + t34) * t27) - m(7) * (g(1) * (-rSges(7,2) * t51 + t45 * t3 + t4 * t56 + t43) + g(2) * (t45 * t1 + t2 * t56 + t33) + (-rSges(7,2) * t59 + t34) * t27) (-m(3) - m(4) + t44) * t39, -m(4) * (-g(3) * t37 + t39 * (rSges(4,1) * t28 - rSges(4,2) * t26)) - m(5) * (g(3) * t67 + t39 * (t46 * t26 + t32 * t28)) - m(6) * ((-rSges(6,3) * t58 + g(3) * t31 + t48 * t60) * t26 + (g(3) * t48 + t31 * t58 + t35 * t60) * t28 + t40) - m(7) * ((-rSges(7,2) * t58 + g(3) * t68 + t49 * t60) * t26 + (g(3) * t49 + t68 * t58 + t64 * t60) * t28 + t40) t44 * (g(3) * t26 - t39 * t28) -m(6) * (g(1) * (-rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (rSges(6,1) * t3 + rSges(6,2) * t4)) - m(7) * (g(1) * (-t56 * t1 + t45 * t2) + g(2) * (t56 * t3 - t45 * t4)) + (-m(6) * (-rSges(6,1) * t16 - rSges(6,2) * t17) - m(7) * (-t56 * t16 + t45 * t17)) * t57, -m(7) * (g(1) * t1 - g(2) * t3 + t16 * t57)];
taug  = t5(:);
