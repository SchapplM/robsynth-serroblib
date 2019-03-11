% Calculate Gravitation load on the joints for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:48
% EndTime: 2019-03-09 01:31:49
% DurationCPUTime: 0.35s
% Computational Cost: add. (260->78), mult. (217->105), div. (0->0), fcn. (183->10), ass. (0->38)
t19 = qJ(1) + pkin(9);
t14 = sin(t19);
t16 = cos(t19);
t54 = -g(1) * t14 + g(2) * t16;
t44 = rSges(7,3) + pkin(8);
t23 = sin(qJ(6));
t25 = cos(qJ(6));
t27 = m(6) * rSges(6,1) + m(7) * (rSges(7,1) * t25 - rSges(7,2) * t23 + pkin(5));
t51 = -m(6) * rSges(6,2) + m(7) * t44;
t20 = sin(pkin(10));
t48 = pkin(4) * t20;
t24 = sin(qJ(1));
t45 = t24 * pkin(1);
t43 = t14 * t23;
t42 = t14 * t25;
t41 = t16 * t23;
t40 = t16 * t25;
t39 = rSges(5,3) + qJ(4);
t38 = -m(5) - m(6) - m(7);
t26 = cos(qJ(1));
t17 = t26 * pkin(1);
t37 = t16 * pkin(2) + t14 * qJ(3) + t17;
t36 = -m(4) + t38;
t35 = t16 * qJ(3) - t45;
t34 = t14 * t48 + t37;
t22 = -pkin(7) - qJ(4);
t33 = t14 * t22 + t16 * t48 + t35;
t32 = rSges(5,1) * t20 + rSges(5,2) * cos(pkin(10));
t18 = pkin(10) + qJ(5);
t13 = sin(t18);
t15 = cos(t18);
t31 = t13 * rSges(6,1) + t15 * rSges(6,2);
t28 = t13 * pkin(5) - t44 * t15;
t4 = t13 * t40 - t43;
t3 = t13 * t41 + t42;
t2 = t13 * t42 + t41;
t1 = -t13 * t43 + t40;
t5 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - t26 * rSges(2,2)) + g(2) * (t26 * rSges(2,1) - t24 * rSges(2,2))) - m(3) * (g(1) * (-t14 * rSges(3,1) - t16 * rSges(3,2) - t45) + g(2) * (t16 * rSges(3,1) - t14 * rSges(3,2) + t17)) - m(4) * (g(1) * (t16 * rSges(4,3) + (rSges(4,2) - pkin(2)) * t14 + t35) + g(2) * (-t16 * rSges(4,2) + t14 * rSges(4,3) + t37)) - m(5) * (g(1) * t35 + g(2) * t37 + (g(1) * t32 + g(2) * t39) * t16 + (g(1) * (-pkin(2) - t39) + g(2) * t32) * t14) - m(6) * (g(1) * t33 + g(2) * t34 + (g(1) * t31 + g(2) * (rSges(6,3) - t22)) * t16 + (g(1) * (-rSges(6,3) - pkin(2)) + g(2) * t31) * t14) - m(7) * (g(1) * (t4 * rSges(7,1) - t3 * rSges(7,2) + t33) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t34) + (g(1) * t28 - g(2) * t22) * t16 + (-g(1) * pkin(2) + g(2) * t28) * t14) (-m(3) + t36) * g(3), -t36 * t54, t38 * (g(1) * t16 + g(2) * t14) (t27 * t13 - t15 * t51) * g(3) + t54 * (t51 * t13 + t27 * t15) -m(7) * (g(1) * (t1 * rSges(7,1) - t2 * rSges(7,2)) + g(2) * (t3 * rSges(7,1) + t4 * rSges(7,2)) + g(3) * (-rSges(7,1) * t23 - rSges(7,2) * t25) * t15)];
taug  = t5(:);
