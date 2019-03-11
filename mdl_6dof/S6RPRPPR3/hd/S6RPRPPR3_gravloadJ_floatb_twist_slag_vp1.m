% Calculate Gravitation load on the joints for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:07
% EndTime: 2019-03-09 02:44:08
% DurationCPUTime: 0.48s
% Computational Cost: add. (294->107), mult. (316->142), div. (0->0), fcn. (276->8), ass. (0->50)
t57 = rSges(7,3) + pkin(8);
t23 = sin(qJ(3));
t26 = cos(qJ(3));
t31 = t26 * rSges(5,1) + t23 * rSges(5,3);
t33 = t26 * rSges(4,1) - t23 * rSges(4,2);
t21 = qJ(1) + pkin(9);
t16 = cos(t21);
t15 = sin(t21);
t61 = g(2) * t15;
t65 = g(1) * t16 + t61;
t64 = -m(6) - m(7);
t63 = -pkin(3) - pkin(4);
t60 = g(3) * t26;
t59 = t23 * pkin(5);
t24 = sin(qJ(1));
t58 = t24 * pkin(1);
t19 = t26 * pkin(3);
t56 = t16 * t26;
t22 = sin(qJ(6));
t55 = t22 * t23;
t54 = t23 * rSges(6,1);
t25 = cos(qJ(6));
t51 = t23 * t25;
t48 = t26 * rSges(6,2);
t17 = t23 * qJ(4);
t47 = t17 + t19;
t46 = qJ(4) * t26;
t45 = -rSges(6,3) - qJ(5);
t44 = -m(5) + t64;
t43 = rSges(6,2) + t63;
t27 = cos(qJ(1));
t20 = t27 * pkin(1);
t42 = t16 * pkin(2) + t15 * pkin(7) + t20;
t41 = t26 * pkin(4) + t47;
t40 = -t57 + t63;
t39 = t16 * pkin(7) - t58;
t38 = -pkin(2) - t17;
t37 = g(1) * t43;
t7 = t15 * t46;
t9 = t16 * t46;
t36 = g(1) * t9 + g(2) * t7;
t35 = pkin(3) * t56 + t16 * t17 + t42;
t34 = g(1) * t40;
t30 = pkin(4) * t56 + t35;
t29 = rSges(7,1) * t25 - rSges(7,2) * t22 + pkin(5);
t5 = -t15 * t22 + t16 * t51;
t4 = -t15 * t25 - t16 * t55;
t3 = -t15 * t51 - t16 * t22;
t2 = t15 * t55 - t16 * t25;
t1 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - t27 * rSges(2,2)) + g(2) * (t27 * rSges(2,1) - t24 * rSges(2,2))) - m(3) * (g(1) * (-t15 * rSges(3,1) - t16 * rSges(3,2) - t58) + g(2) * (t16 * rSges(3,1) - t15 * rSges(3,2) + t20)) - m(4) * (g(1) * (t16 * rSges(4,3) + t39) + g(2) * (t33 * t16 + t42) + (g(1) * (-pkin(2) - t33) + g(2) * rSges(4,3)) * t15) - m(5) * (g(1) * (t16 * rSges(5,2) + t39) + g(2) * (t31 * t16 + t35) + (g(1) * (-t31 + t38 - t19) + g(2) * rSges(5,2)) * t15) - m(6) * (g(1) * t39 + g(2) * t30 + (g(1) * t45 + g(2) * (-t48 + t54)) * t16 + (g(1) * (t38 - t54) + g(2) * t45 + t26 * t37) * t15) - m(7) * (g(1) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t39) + g(2) * (t5 * rSges(7,1) + t4 * rSges(7,2) + t30) + (-g(1) * qJ(5) + g(2) * (t57 * t26 + t59)) * t16 + (g(1) * (t38 - t59) - g(2) * qJ(5) + t26 * t34) * t15) (-m(3) - m(4) + t44) * g(3), -m(4) * (g(3) * t33 + t65 * (-rSges(4,1) * t23 - rSges(4,2) * t26)) - m(5) * (g(3) * (t31 + t47) + t36 + t65 * (rSges(5,3) * t26 + (-rSges(5,1) - pkin(3)) * t23)) - m(6) * (g(1) * (rSges(6,1) * t56 + t9) + g(2) * (t15 * t26 * rSges(6,1) + t7) + g(3) * (t41 - t48) + (g(3) * rSges(6,1) + t16 * t37 + t43 * t61) * t23) - m(7) * (g(3) * t41 + (g(3) * t57 + t65 * t29) * t26 + (g(3) * t29 + t16 * t34 + t40 * t61) * t23 + t36) t44 * (t65 * t23 - t60) t64 * (-g(1) * t15 + g(2) * t16) -m(7) * (g(1) * (t4 * rSges(7,1) - t5 * rSges(7,2)) + g(2) * (-t2 * rSges(7,1) + t3 * rSges(7,2)) + (rSges(7,1) * t22 + rSges(7,2) * t25) * t60)];
taug  = t1(:);
