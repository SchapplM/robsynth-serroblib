% Calculate Gravitation load on the joints for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:01
% EndTime: 2019-03-09 01:30:01
% DurationCPUTime: 0.31s
% Computational Cost: add. (211->73), mult. (207->101), div. (0->0), fcn. (173->8), ass. (0->31)
t41 = rSges(7,3) + pkin(8);
t14 = qJ(1) + pkin(9);
t11 = sin(t14);
t12 = cos(t14);
t40 = g(1) * t12 + g(2) * t11;
t15 = sin(qJ(6));
t18 = cos(qJ(6));
t39 = m(6) * rSges(6,1) + m(7) * (rSges(7,1) * t18 - rSges(7,2) * t15 + pkin(5));
t17 = sin(qJ(1));
t37 = pkin(1) * t17;
t34 = -rSges(6,3) - pkin(7);
t16 = sin(qJ(5));
t33 = t15 * t16;
t32 = t16 * t18;
t31 = -pkin(2) - qJ(4);
t30 = -m(5) - m(6) - m(7);
t20 = cos(qJ(1));
t13 = t20 * pkin(1);
t29 = t12 * pkin(2) + t11 * qJ(3) + t13;
t28 = -m(4) + t30;
t27 = t12 * qJ(3) - t37;
t26 = t12 * qJ(4) + t29;
t19 = cos(qJ(5));
t25 = rSges(6,1) * t16 + rSges(6,2) * t19;
t23 = m(6) * rSges(6,2) - m(7) * t41;
t22 = pkin(5) * t16 - t41 * t19;
t4 = -t11 * t15 + t12 * t32;
t3 = -t11 * t18 - t12 * t33;
t2 = -t11 * t32 - t12 * t15;
t1 = t11 * t33 - t12 * t18;
t5 = [-m(2) * (g(1) * (-t17 * rSges(2,1) - rSges(2,2) * t20) + g(2) * (rSges(2,1) * t20 - t17 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t11 - rSges(3,2) * t12 - t37) + g(2) * (rSges(3,1) * t12 - rSges(3,2) * t11 + t13)) - m(4) * (g(1) * (rSges(4,3) * t12 + (rSges(4,2) - pkin(2)) * t11 + t27) + g(2) * (-rSges(4,2) * t12 + rSges(4,3) * t11 + t29)) - m(5) * (g(1) * (rSges(5,2) * t12 + t27) + g(2) * (rSges(5,3) * t12 + t26) + (g(1) * (-rSges(5,3) + t31) + g(2) * rSges(5,2)) * t11) - m(6) * (g(1) * t27 + g(2) * t26 + (g(1) * t34 + g(2) * t25) * t12 + (g(1) * (-t25 + t31) + g(2) * t34) * t11) - m(7) * (g(1) * (rSges(7,1) * t2 + rSges(7,2) * t1 + t27) + g(2) * (rSges(7,1) * t4 + rSges(7,2) * t3 + t26) + (-g(1) * pkin(7) + g(2) * t22) * t12 + (g(1) * (-t22 + t31) - g(2) * pkin(7)) * t11) (-m(3) + t28) * g(3), t28 * (g(1) * t11 - g(2) * t12) t30 * t40 (t39 * t16 + t23 * t19) * g(3) + t40 * (t23 * t16 - t39 * t19) -m(7) * (g(1) * (rSges(7,1) * t3 - rSges(7,2) * t4) + g(2) * (-rSges(7,1) * t1 + rSges(7,2) * t2) + g(3) * (-rSges(7,1) * t15 - rSges(7,2) * t18) * t19)];
taug  = t5(:);
