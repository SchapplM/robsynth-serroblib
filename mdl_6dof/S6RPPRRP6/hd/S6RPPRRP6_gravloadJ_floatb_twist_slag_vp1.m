% Calculate Gravitation load on the joints for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:17
% EndTime: 2019-03-09 02:10:18
% DurationCPUTime: 0.51s
% Computational Cost: add. (178->100), mult. (352->137), div. (0->0), fcn. (330->6), ass. (0->41)
t58 = -m(6) - m(7);
t57 = rSges(7,1) + pkin(5);
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t56 = g(1) * t25 + g(2) * t22;
t36 = rSges(7,3) + qJ(6);
t55 = g(2) * pkin(7);
t54 = m(5) * rSges(5,1);
t24 = cos(qJ(4));
t51 = g(3) * t24;
t49 = -rSges(5,3) - pkin(7);
t21 = sin(qJ(4));
t48 = t21 * t25;
t20 = sin(qJ(5));
t47 = t22 * t20;
t23 = cos(qJ(5));
t46 = t22 * t23;
t45 = t22 * t24;
t44 = t24 * rSges(5,2);
t43 = t24 * rSges(7,2);
t42 = t24 * rSges(6,3);
t41 = t24 * t25;
t40 = t25 * t20;
t39 = t25 * t23;
t38 = -pkin(1) - qJ(3);
t37 = t25 * pkin(1) + t22 * qJ(2);
t35 = t25 * qJ(3) + t37;
t34 = -m(4) - m(5) + t58;
t17 = t25 * qJ(2);
t33 = -t25 * pkin(7) + pkin(8) * t45 + t17;
t32 = -t21 * pkin(4) + t38;
t31 = t21 * rSges(5,1) + t44;
t30 = rSges(6,1) * t23 - rSges(6,2) * t20;
t29 = pkin(4) * t48 - pkin(8) * t41 + t35;
t27 = t36 * t20 + t57 * t23;
t18 = t24 * pkin(8);
t4 = t21 * t39 - t47;
t3 = t21 * t40 + t46;
t2 = t21 * t46 + t40;
t1 = t21 * t47 - t39;
t5 = [-m(2) * (g(1) * (-t22 * rSges(2,1) - t25 * rSges(2,2)) + g(2) * (t25 * rSges(2,1) - t22 * rSges(2,2))) - m(3) * (g(1) * (t25 * rSges(3,3) + t17 + (rSges(3,2) - pkin(1)) * t22) + g(2) * (-t25 * rSges(3,2) + t22 * rSges(3,3) + t37)) - m(4) * (g(1) * (t25 * rSges(4,2) + t17) + g(2) * (t25 * rSges(4,3) + t35) + (g(1) * (-rSges(4,3) + t38) + g(2) * rSges(4,2)) * t22) - m(5) * (g(1) * t17 + g(2) * t35 + (g(1) * t49 + g(2) * t31) * t25 + (g(1) * (-t31 + t38) + g(2) * t49) * t22) - m(6) * (g(1) * (-t2 * rSges(6,1) + t1 * rSges(6,2) + t33) + g(2) * (t4 * rSges(6,1) - t3 * rSges(6,2) - rSges(6,3) * t41 + t29) + (g(1) * (t32 + t42) - t55) * t22) - m(7) * (g(1) * (-t36 * t1 - t57 * t2 + t33) + g(2) * (-rSges(7,2) * t41 + t36 * t3 + t57 * t4 + t29) + (g(1) * (t32 + t43) - t55) * t22) (-m(3) + t34) * (g(1) * t22 - g(2) * t25) t34 * t56, t58 * (g(1) * (pkin(4) * t41 + pkin(8) * t48) + g(2) * (t22 * t21 * pkin(8) + pkin(4) * t45)) + t56 * ((m(5) * rSges(5,2) - m(6) * rSges(6,3) - m(7) * rSges(7,2)) * t21 + (-m(6) * t30 - m(7) * t27 - t54) * t24) + (m(5) * t44 - m(6) * (t18 + t42) - m(7) * (t18 + t43) + (t54 - m(6) * (-pkin(4) - t30) - m(7) * (-pkin(4) - t27)) * t21) * g(3), -m(6) * (g(1) * (-t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) - t2 * rSges(6,2))) - m(7) * (g(1) * (-t3 * t57 + t36 * t4) + g(2) * (-t1 * t57 + t36 * t2)) + (-m(6) * (-rSges(6,1) * t20 - rSges(6,2) * t23) - m(7) * (-t57 * t20 + t36 * t23)) * t51, -m(7) * (g(1) * t3 + g(2) * t1 + t20 * t51)];
taug  = t5(:);
