% Calculate Gravitation load on the joints for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:56:13
% EndTime: 2019-03-09 02:56:14
% DurationCPUTime: 0.43s
% Computational Cost: add. (238->112), mult. (305->141), div. (0->0), fcn. (265->8), ass. (0->49)
t53 = rSges(7,3) + pkin(8);
t29 = cos(qJ(1));
t58 = g(2) * t29;
t26 = sin(qJ(1));
t59 = g(1) * t26;
t37 = -t58 + t59;
t61 = -m(6) - m(7);
t28 = cos(qJ(3));
t60 = pkin(3) * t28;
t22 = qJ(3) + pkin(9);
t17 = sin(t22);
t57 = g(3) * t17;
t25 = sin(qJ(3));
t56 = t25 * pkin(3);
t55 = rSges(6,2) - pkin(4);
t54 = rSges(4,3) + pkin(7);
t52 = t17 * t26;
t18 = cos(t22);
t51 = t18 * rSges(6,3);
t24 = sin(qJ(6));
t50 = t26 * t24;
t27 = cos(qJ(6));
t49 = t26 * t27;
t48 = t29 * t24;
t47 = t29 * t27;
t46 = t29 * pkin(1) + t26 * qJ(2);
t15 = t18 * qJ(5);
t45 = -m(5) + t61;
t44 = -pkin(4) - t53;
t13 = t26 * t60;
t43 = t26 * t18 * pkin(4) + qJ(5) * t52 + t13;
t42 = t26 * t56 + t46;
t20 = t29 * qJ(2);
t23 = -qJ(4) - pkin(7);
t41 = t26 * t23 + t29 * t56 + t20;
t40 = t15 - t56;
t39 = pkin(4) * t52 + t42;
t38 = t29 * t17 * pkin(4) + t41;
t35 = t25 * rSges(4,1) + t28 * rSges(4,2);
t34 = rSges(5,1) * t18 - rSges(5,2) * t17;
t33 = t17 * rSges(5,1) + t18 * rSges(5,2);
t32 = rSges(7,1) * t24 + rSges(7,2) * t27;
t31 = t53 * t17 - t15;
t30 = -t17 * rSges(6,2) - t15 - t51;
t5 = -t18 * t50 + t47;
t4 = -t18 * t49 - t48;
t3 = -t18 * t48 - t49;
t2 = -t18 * t47 + t50;
t1 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - t29 * rSges(2,2)) + g(2) * (t29 * rSges(2,1) - t26 * rSges(2,2))) - m(3) * (g(1) * (t29 * rSges(3,3) + t20 + (rSges(3,2) - pkin(1)) * t26) + g(2) * (-t29 * rSges(3,2) + t26 * rSges(3,3) + t46)) - m(4) * (g(1) * t20 + g(2) * t46 + (g(1) * t35 + g(2) * t54) * t29 + (g(1) * (-pkin(1) - t54) + g(2) * t35) * t26) - m(5) * (g(1) * t41 + g(2) * t42 + (g(1) * t33 + g(2) * (rSges(5,3) - t23)) * t29 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t33) * t26) - m(6) * (g(1) * t38 + g(2) * t39 + (g(1) * t30 + g(2) * (rSges(6,1) - t23)) * t29 + (g(1) * (-rSges(6,1) - pkin(1)) + g(2) * t30) * t26) - m(7) * (g(1) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t38) + g(2) * (t5 * rSges(7,1) + t4 * rSges(7,2) + t39) + (g(1) * t31 + g(2) * (pkin(5) - t23)) * t29 + (g(1) * (-pkin(1) - pkin(5)) + g(2) * t31) * t26) (-m(3) - m(4) + t45) * t37, -m(4) * (-g(3) * t35 + t37 * (rSges(4,1) * t28 - rSges(4,2) * t25)) - m(5) * (g(1) * (t34 * t26 + t13) + g(3) * (-t33 - t56) + (-t34 - t60) * t58) - m(6) * (g(1) * ((-rSges(6,2) * t18 + rSges(6,3) * t17) * t26 + t43) + g(3) * (t55 * t17 + t40 + t51) + (-t60 + t55 * t18 + (-rSges(6,3) - qJ(5)) * t17) * t58) - m(7) * (g(1) * t43 + g(3) * t40 + (g(3) * t32 + t53 * t59) * t18 + (g(3) * t44 + t32 * t59) * t17 + (-t60 + t44 * t18 + (-qJ(5) - t32) * t17) * t58) t45 * (g(1) * t29 + g(2) * t26) t61 * (-t37 * t18 + t57) -m(7) * (g(1) * (t4 * rSges(7,1) - t5 * rSges(7,2)) + g(2) * (-t2 * rSges(7,1) + t3 * rSges(7,2)) + (rSges(7,1) * t27 - rSges(7,2) * t24) * t57)];
taug  = t1(:);
