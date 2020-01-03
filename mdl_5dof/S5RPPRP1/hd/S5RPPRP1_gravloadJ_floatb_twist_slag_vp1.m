% Calculate Gravitation load on the joints for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:08
% EndTime: 2020-01-03 11:25:10
% DurationCPUTime: 0.31s
% Computational Cost: add. (182->70), mult. (194->100), div. (0->0), fcn. (179->8), ass. (0->27)
t37 = rSges(6,1) + pkin(4);
t17 = sin(pkin(8));
t18 = cos(pkin(8));
t36 = rSges(4,1) * t18 - rSges(4,2) * t17;
t35 = pkin(3) * t18 + (rSges(5,3) + pkin(6)) * t17;
t20 = sin(qJ(4));
t33 = pkin(4) * t20;
t16 = qJ(1) + pkin(7);
t12 = sin(t16);
t32 = g(3) * t12;
t21 = sin(qJ(1));
t14 = t21 * pkin(1);
t31 = t12 * pkin(2) + t14;
t28 = t18 * t20;
t22 = cos(qJ(4));
t27 = t18 * t22;
t26 = -m(4) - m(5) - m(6);
t13 = cos(t16);
t23 = cos(qJ(1));
t15 = t23 * pkin(1);
t25 = t13 * pkin(2) + t12 * qJ(3) + t15;
t3 = -t12 * t22 + t13 * t28;
t1 = -t12 * t28 - t13 * t22;
t24 = (pkin(4) * t22 + pkin(3)) * t18 + (rSges(6,3) + qJ(5) + pkin(6)) * t17;
t4 = t12 * t20 + t13 * t27;
t2 = t12 * t27 - t13 * t20;
t5 = [-m(2) * (g(2) * (rSges(2,1) * t23 - t21 * rSges(2,2)) + g(3) * (t21 * rSges(2,1) + rSges(2,2) * t23)) - m(3) * (g(2) * (rSges(3,1) * t13 - rSges(3,2) * t12 + t15) + g(3) * (rSges(3,1) * t12 + rSges(3,2) * t13 + t14)) - m(4) * (g(2) * (rSges(4,3) * t12 + t25) + g(3) * (t36 * t12 + t31) + (g(2) * t36 + g(3) * (-rSges(4,3) - qJ(3))) * t13) - m(5) * (g(2) * (rSges(5,1) * t4 - rSges(5,2) * t3 + t25) + g(3) * (rSges(5,1) * t2 + rSges(5,2) * t1 + t31) + (g(2) * t35 - g(3) * qJ(3)) * t13 + t35 * t32) - m(6) * (g(2) * (rSges(6,1) * t4 - rSges(6,2) * t3 + t25) + g(3) * (rSges(6,1) * t2 + rSges(6,2) * t1 + t31) + (g(2) * t33 + g(3) * t24) * t12 + (g(2) * t24 + g(3) * (-qJ(3) - t33)) * t13), (-m(3) + t26) * g(1), t26 * (-g(2) * t13 - t32), -m(5) * (g(2) * (rSges(5,1) * t1 - rSges(5,2) * t2) + g(3) * (rSges(5,1) * t3 + rSges(5,2) * t4)) - m(6) * (g(2) * (-t2 * rSges(6,2) + t37 * t1) + g(3) * (t4 * rSges(6,2) + t37 * t3)) + (-m(5) * (-rSges(5,1) * t20 - rSges(5,2) * t22) - m(6) * (-rSges(6,1) * t20 - rSges(6,2) * t22 - t33)) * g(1) * t17, -m(6) * (-g(1) * t18 + (g(2) * t12 - g(3) * t13) * t17)];
taug = t5(:);
