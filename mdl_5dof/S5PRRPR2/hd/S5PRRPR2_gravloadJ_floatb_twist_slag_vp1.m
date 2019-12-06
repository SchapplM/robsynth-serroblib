% Calculate Gravitation load on the joints for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:18
% EndTime: 2019-12-05 16:17:19
% DurationCPUTime: 0.20s
% Computational Cost: add. (251->58), mult. (166->77), div. (0->0), fcn. (146->8), ass. (0->33)
t47 = rSges(6,3) + pkin(7);
t46 = -m(5) - m(6);
t27 = pkin(8) + qJ(2);
t24 = sin(t27);
t45 = pkin(2) * t24;
t26 = qJ(3) + t27;
t22 = sin(t26);
t44 = g(1) * t22;
t23 = cos(t26);
t28 = sin(pkin(9));
t43 = t23 * t28;
t29 = cos(pkin(9));
t42 = t23 * t29;
t30 = sin(qJ(5));
t41 = t29 * t30;
t31 = cos(qJ(5));
t40 = t29 * t31;
t39 = t23 * pkin(3) + t22 * qJ(4);
t16 = t23 * qJ(4);
t5 = t22 * t41 + t23 * t31;
t6 = -t22 * t40 + t23 * t30;
t38 = t6 * rSges(6,1) + t5 * rSges(6,2) + t16;
t37 = t23 * rSges(4,1) - rSges(4,2) * t22;
t36 = -rSges(4,1) * t22 - rSges(4,2) * t23;
t7 = t22 * t31 - t23 * t41;
t8 = t22 * t30 + t23 * t40;
t35 = t8 * rSges(6,1) + t7 * rSges(6,2) + pkin(4) * t42 + t47 * t43 + t39;
t34 = rSges(5,1) * t42 - rSges(5,2) * t43 + t22 * rSges(5,3) + t39;
t33 = t23 * rSges(5,3) + t16 + (-rSges(5,1) * t29 + rSges(5,2) * t28 - pkin(3)) * t22;
t32 = (-pkin(4) * t29 - t47 * t28 - pkin(3)) * t44;
t25 = cos(t27);
t21 = pkin(2) * t25;
t1 = [(-m(2) - m(3) - m(4) + t46) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t24 - rSges(3,2) * t25) + g(2) * (rSges(3,1) * t25 - rSges(3,2) * t24)) - m(4) * (g(1) * (t36 - t45) + g(2) * (t21 + t37)) - m(5) * (g(1) * (t33 - t45) + g(2) * (t21 + t34)) - m(6) * (g(1) * (t38 - t45) + g(2) * (t21 + t35) + t32), -m(4) * (g(1) * t36 + g(2) * t37) - m(5) * (g(1) * t33 + g(2) * t34) - m(6) * (g(1) * t38 + g(2) * t35 + t32), t46 * (-g(2) * t23 + t44), -m(6) * (g(1) * (rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (-rSges(6,1) * t5 + rSges(6,2) * t6) + g(3) * (-rSges(6,1) * t30 - rSges(6,2) * t31) * t28)];
taug = t1(:);
