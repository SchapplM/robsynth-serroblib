% Calculate Gravitation load on the joints for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:30:33
% EndTime: 2019-03-09 02:30:34
% DurationCPUTime: 0.43s
% Computational Cost: add. (210->98), mult. (334->135), div. (0->0), fcn. (303->8), ass. (0->47)
t62 = rSges(6,3) + pkin(8);
t61 = rSges(7,3) + pkin(9) + pkin(8);
t25 = sin(qJ(1));
t28 = cos(qJ(1));
t60 = g(1) * t28 + g(2) * t25;
t26 = cos(qJ(5));
t15 = t26 * pkin(5) + pkin(4);
t22 = qJ(5) + qJ(6);
t16 = sin(t22);
t17 = cos(t22);
t23 = sin(qJ(5));
t59 = m(6) * (rSges(6,1) * t26 - rSges(6,2) * t23 + pkin(4)) + m(7) * (rSges(7,1) * t17 - rSges(7,2) * t16 + t15) + m(5) * rSges(5,1);
t24 = sin(qJ(4));
t45 = t28 * t17;
t50 = t25 * t16;
t5 = t24 * t50 - t45;
t46 = t28 * t16;
t49 = t25 * t17;
t6 = -t24 * t49 - t46;
t58 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t7 = -t24 * t46 - t49;
t8 = t24 * t45 - t50;
t57 = t7 * rSges(7,1) - t8 * rSges(7,2);
t55 = pkin(5) * t23;
t27 = cos(qJ(4));
t52 = g(3) * t27;
t51 = -rSges(5,3) - pkin(7);
t48 = t25 * t23;
t47 = t25 * t26;
t44 = t28 * t23;
t43 = t28 * t26;
t42 = -pkin(1) - qJ(3);
t41 = t28 * pkin(1) + t25 * qJ(2);
t40 = t28 * qJ(3) + t41;
t39 = -m(4) - m(5) - m(6) - m(7);
t38 = -pkin(7) - t55;
t37 = t24 * rSges(5,1) + t27 * rSges(5,2);
t36 = -rSges(7,1) * t16 - rSges(7,2) * t17;
t11 = -t24 * t44 - t47;
t9 = t24 * t48 - t43;
t33 = t24 * pkin(4) - t62 * t27;
t32 = t24 * t15 - t61 * t27;
t31 = m(5) * rSges(5,2) - m(6) * t62 - m(7) * t61;
t20 = t28 * qJ(2);
t12 = t24 * t43 - t48;
t10 = -t24 * t47 - t44;
t1 = [-m(2) * (g(1) * (-t25 * rSges(2,1) - t28 * rSges(2,2)) + g(2) * (t28 * rSges(2,1) - t25 * rSges(2,2))) - m(3) * (g(1) * (t28 * rSges(3,3) + t20 + (rSges(3,2) - pkin(1)) * t25) + g(2) * (-t28 * rSges(3,2) + t25 * rSges(3,3) + t41)) - m(4) * (g(1) * (t28 * rSges(4,2) + t20) + g(2) * (t28 * rSges(4,3) + t40) + (g(1) * (-rSges(4,3) + t42) + g(2) * rSges(4,2)) * t25) - m(5) * (g(1) * t20 + g(2) * t40 + (g(1) * t51 + g(2) * t37) * t28 + (g(1) * (-t37 + t42) + g(2) * t51) * t25) - m(6) * (g(1) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t20) + g(2) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t40) + (-g(1) * pkin(7) + g(2) * t33) * t28 + (g(1) * (-t33 + t42) - g(2) * pkin(7)) * t25) - m(7) * (g(1) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t20) + g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) + t40) + (g(1) * t38 + g(2) * t32) * t28 + (g(1) * (-t32 + t42) + g(2) * t38) * t25) (-m(3) + t39) * (g(1) * t25 - g(2) * t28) t39 * t60 (t59 * t24 + t31 * t27) * g(3) + t60 * (t31 * t24 - t59 * t27) -m(6) * (g(1) * (t11 * rSges(6,1) - t12 * rSges(6,2)) + g(2) * (-t9 * rSges(6,1) + t10 * rSges(6,2))) - m(7) * (g(1) * (t11 * pkin(5) + t57) + g(2) * (-t9 * pkin(5) + t58)) + (-m(6) * (-rSges(6,1) * t23 - rSges(6,2) * t26) - m(7) * (t36 - t55)) * t52, -m(7) * (g(1) * t57 + g(2) * t58 + t36 * t52)];
taug  = t1(:);
