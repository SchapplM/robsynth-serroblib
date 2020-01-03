% Calculate Gravitation load on the joints for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:09
% EndTime: 2019-12-31 18:42:11
% DurationCPUTime: 0.39s
% Computational Cost: add. (236->80), mult. (272->111), div. (0->0), fcn. (248->8), ass. (0->31)
t43 = rSges(6,1) + pkin(4);
t42 = rSges(6,3) + qJ(5) + pkin(7);
t15 = sin(qJ(3));
t18 = cos(qJ(3));
t41 = t18 * rSges(4,1) - t15 * rSges(4,2);
t12 = qJ(1) + pkin(8);
t10 = cos(t12);
t9 = sin(t12);
t40 = g(1) * t10 + g(2) * t9;
t31 = rSges(5,3) + pkin(7);
t39 = t18 * pkin(3) + t31 * t15;
t14 = sin(qJ(4));
t17 = cos(qJ(4));
t8 = t17 * pkin(4) + pkin(3);
t38 = m(5) * (rSges(5,1) * t17 - rSges(5,2) * t14 + pkin(3)) + m(6) * (rSges(6,1) * t17 - rSges(6,2) * t14 + t8) + m(4) * rSges(4,1);
t35 = pkin(4) * t14;
t16 = sin(qJ(1));
t33 = t16 * pkin(1);
t30 = t14 * t18;
t28 = t17 * t18;
t19 = cos(qJ(1));
t11 = t19 * pkin(1);
t26 = t10 * pkin(2) + t9 * pkin(6) + t11;
t25 = t10 * pkin(6) - t33;
t1 = t10 * t17 + t9 * t30;
t3 = -t10 * t30 + t9 * t17;
t22 = t42 * t15 + t18 * t8;
t21 = m(4) * rSges(4,2) - m(5) * t31 - m(6) * t42;
t4 = t10 * t28 + t9 * t14;
t2 = t10 * t14 - t9 * t28;
t5 = [-m(2) * (g(1) * (-t16 * rSges(2,1) - t19 * rSges(2,2)) + g(2) * (t19 * rSges(2,1) - t16 * rSges(2,2))) - m(3) * (g(1) * (-t9 * rSges(3,1) - t10 * rSges(3,2) - t33) + g(2) * (t10 * rSges(3,1) - t9 * rSges(3,2) + t11)) - m(4) * (g(1) * (t10 * rSges(4,3) + t25) + g(2) * (t41 * t10 + t26) + (g(1) * (-pkin(2) - t41) + g(2) * rSges(4,3)) * t9) - m(5) * ((t4 * rSges(5,1) + t3 * rSges(5,2) + t39 * t10 + t26) * g(2) + (t2 * rSges(5,1) + t1 * rSges(5,2) + t25 + (-pkin(2) - t39) * t9) * g(1)) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t25) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t26) + (g(1) * t35 + g(2) * t22) * t10 + (g(1) * (-pkin(2) - t22) + g(2) * t35) * t9), (-m(3) - m(4) - m(5) - m(6)) * g(3), (t21 * t15 - t38 * t18) * g(3) + t40 * (t38 * t15 + t21 * t18), -m(5) * (g(1) * (t3 * rSges(5,1) - t4 * rSges(5,2)) + g(2) * (-t1 * rSges(5,1) + t2 * rSges(5,2))) - m(6) * (g(1) * (-t4 * rSges(6,2) + t43 * t3) + g(2) * (t2 * rSges(6,2) - t43 * t1)) + (-m(5) * (-rSges(5,1) * t14 - rSges(5,2) * t17) - m(6) * (-rSges(6,1) * t14 - rSges(6,2) * t17 - t35)) * g(3) * t15, -m(6) * (-g(3) * t18 + t40 * t15)];
taug = t5(:);
