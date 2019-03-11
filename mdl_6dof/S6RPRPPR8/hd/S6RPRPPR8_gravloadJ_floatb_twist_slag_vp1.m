% Calculate Gravitation load on the joints for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:57
% EndTime: 2019-03-09 02:58:58
% DurationCPUTime: 0.51s
% Computational Cost: add. (173->111), mult. (328->143), div. (0->0), fcn. (288->6), ass. (0->46)
t53 = rSges(7,3) + pkin(8);
t25 = cos(qJ(3));
t17 = t25 * qJ(4);
t62 = -t25 * rSges(5,3) - t17;
t61 = t25 * rSges(6,1) + t17;
t26 = cos(qJ(1));
t56 = g(2) * t26;
t23 = sin(qJ(1));
t57 = g(1) * t23;
t33 = -t56 + t57;
t60 = -m(6) - m(7);
t59 = -pkin(1) - pkin(7);
t58 = -pkin(3) - pkin(4);
t22 = sin(qJ(3));
t55 = g(3) * t22;
t54 = -rSges(5,1) - pkin(3);
t49 = t23 * t25;
t51 = t22 * t23;
t52 = pkin(3) * t49 + qJ(4) * t51;
t50 = t22 * t26;
t47 = t25 * rSges(4,2);
t21 = sin(qJ(6));
t45 = t26 * t21;
t24 = cos(qJ(6));
t44 = t26 * t24;
t18 = t26 * qJ(2);
t43 = pkin(3) * t50 + t18;
t42 = t26 * pkin(1) + t23 * qJ(2);
t41 = -m(5) + t60;
t40 = rSges(6,2) + t58;
t39 = pkin(4) * t49 + t52;
t38 = t26 * pkin(7) + t42;
t37 = -t53 + t58;
t36 = pkin(3) * t51 + t38;
t35 = pkin(4) * t50 + t23 * qJ(5) + t43;
t34 = pkin(4) * t51 + t36;
t31 = t22 * rSges(4,1) + t47;
t30 = rSges(7,1) * t24 - rSges(7,2) * t21 + pkin(5);
t29 = t22 * rSges(5,1) + t62;
t28 = -t22 * rSges(6,2) - t61;
t27 = -t25 * pkin(5) + t53 * t22 - t17;
t5 = t23 * t21 - t25 * t44;
t4 = t23 * t24 + t25 * t45;
t3 = t24 * t49 + t45;
t2 = t21 * t49 - t44;
t1 = [-m(2) * (g(1) * (-t23 * rSges(2,1) - t26 * rSges(2,2)) + g(2) * (t26 * rSges(2,1) - t23 * rSges(2,2))) - m(3) * (g(1) * (t26 * rSges(3,3) + t18 + (rSges(3,2) - pkin(1)) * t23) + g(2) * (-t26 * rSges(3,2) + t23 * rSges(3,3) + t42)) - m(4) * (g(1) * (rSges(4,1) * t50 + t26 * t47 + t18) + g(2) * (t26 * rSges(4,3) + t38) + (g(1) * (-rSges(4,3) + t59) + g(2) * t31) * t23) - m(5) * (g(1) * t43 + g(2) * t36 + (g(2) * rSges(5,2) + g(1) * t29) * t26 + (g(1) * (-rSges(5,2) + t59) + g(2) * t29) * t23) - m(6) * (g(1) * t35 + g(2) * t34 + (g(1) * t28 + g(2) * (-rSges(6,3) - qJ(5))) * t26 + (g(1) * (rSges(6,3) + t59) + g(2) * t28) * t23) - m(7) * (g(1) * (t5 * rSges(7,1) + t4 * rSges(7,2) + t35) + g(2) * (-t3 * rSges(7,1) + t2 * rSges(7,2) + t34) + (g(1) * t27 - g(2) * qJ(5)) * t26 + (g(1) * t59 + g(2) * t27) * t23) (-m(3) - m(4) + t41) * t33, -m(4) * (-g(3) * t31 + t33 * (rSges(4,1) * t25 - rSges(4,2) * t22)) - m(5) * (g(1) * ((rSges(5,1) * t25 + rSges(5,3) * t22) * t23 + t52) + g(3) * (t54 * t22 - t62) + (t54 * t25 + (-rSges(5,3) - qJ(4)) * t22) * t56) - m(6) * (g(1) * (-rSges(6,2) * t49 + t39) + g(3) * t61 + (rSges(6,1) * t57 + g(3) * t40) * t22 + ((-rSges(6,1) - qJ(4)) * t22 + t40 * t25) * t56) - m(7) * (g(1) * t39 + g(3) * t17 + (g(3) * t30 + t37 * t56 + t53 * t57) * t25 + (g(3) * t37 + t30 * t57 + (-qJ(4) - t30) * t56) * t22) t41 * (-t33 * t25 + t55) t60 * (-g(1) * t26 - g(2) * t23) -m(7) * (g(1) * (t2 * rSges(7,1) + t3 * rSges(7,2)) + g(2) * (-t4 * rSges(7,1) + t5 * rSges(7,2)) + (-rSges(7,1) * t21 - rSges(7,2) * t24) * t55)];
taug  = t1(:);
