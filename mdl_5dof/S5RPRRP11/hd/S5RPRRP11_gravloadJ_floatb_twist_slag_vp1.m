% Calculate Gravitation load on the joints for
% S5RPRRP11
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:27
% EndTime: 2019-12-31 18:53:29
% DurationCPUTime: 0.54s
% Computational Cost: add. (252->88), mult. (324->125), div. (0->0), fcn. (310->8), ass. (0->39)
t52 = -m(5) - m(6);
t23 = sin(qJ(1));
t46 = g(2) * t23;
t51 = rSges(6,1) + pkin(4);
t18 = pkin(8) + qJ(3);
t16 = sin(t18);
t17 = cos(t18);
t50 = t17 * pkin(3) + t16 * pkin(7);
t25 = cos(qJ(1));
t49 = g(1) * t25 + t46;
t34 = rSges(6,3) + qJ(5);
t21 = -pkin(6) - qJ(2);
t47 = g(2) * t21;
t45 = g(3) * t16;
t43 = t16 * t25;
t42 = t17 * t25;
t22 = sin(qJ(4));
t41 = t23 * t22;
t24 = cos(qJ(4));
t40 = t23 * t24;
t39 = t25 * t21;
t38 = t25 * t22;
t37 = t25 * t24;
t36 = rSges(4,3) - t21;
t35 = rSges(3,3) + qJ(2);
t20 = cos(pkin(8));
t13 = t20 * pkin(2) + pkin(1);
t5 = t25 * t13;
t33 = pkin(3) * t42 + pkin(7) * t43 + t5;
t32 = t17 * rSges(4,1) - t16 * rSges(4,2);
t31 = rSges(5,1) * t24 - rSges(5,2) * t22;
t30 = -t13 - t50;
t29 = rSges(3,1) * t20 - rSges(3,2) * sin(pkin(8)) + pkin(1);
t27 = t34 * t22 + t51 * t24;
t4 = t17 * t37 + t41;
t3 = t17 * t38 - t40;
t2 = t17 * t40 - t38;
t1 = t17 * t41 + t37;
t6 = [-m(2) * (g(1) * (-t23 * rSges(2,1) - t25 * rSges(2,2)) + g(2) * (t25 * rSges(2,1) - t23 * rSges(2,2))) - m(3) * ((g(1) * t35 + g(2) * t29) * t25 + (-g(1) * t29 + g(2) * t35) * t23) - m(4) * (g(2) * t5 + (g(1) * t36 + g(2) * t32) * t25 + (g(1) * (-t13 - t32) + g(2) * t36) * t23) - m(5) * (g(1) * (-t2 * rSges(5,1) + t1 * rSges(5,2) - t39) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + rSges(5,3) * t43 + t33) + (g(1) * (-t16 * rSges(5,3) + t30) - t47) * t23) - m(6) * (g(1) * (-t34 * t1 - t51 * t2 - t39) + g(2) * (rSges(6,2) * t43 + t34 * t3 + t51 * t4 + t33) + (g(1) * (-t16 * rSges(6,2) + t30) - t47) * t23), (-m(3) - m(4) + t52) * (g(1) * t23 - g(2) * t25), (-m(4) * (g(3) * rSges(4,1) - t49 * rSges(4,2)) - m(5) * (t49 * rSges(5,3) + g(3) * t31) - m(6) * (t49 * rSges(6,2) + g(3) * t27)) * t17 + ((m(4) * rSges(4,2) - m(5) * rSges(5,3) - m(6) * rSges(6,2)) * g(3) + t49 * (m(4) * rSges(4,1) - m(5) * (-pkin(3) - t31) - m(6) * (-pkin(3) - t27))) * t16 + t52 * (g(3) * t50 + (g(1) * t42 + t17 * t46) * pkin(7)), -m(5) * (g(1) * (-t3 * rSges(5,1) - t4 * rSges(5,2)) + g(2) * (-t1 * rSges(5,1) - t2 * rSges(5,2))) - m(6) * (g(1) * (-t3 * t51 + t34 * t4) + g(2) * (-t1 * t51 + t2 * t34)) + (-m(5) * (-rSges(5,1) * t22 - rSges(5,2) * t24) - m(6) * (-t51 * t22 + t34 * t24)) * t45, -m(6) * (g(1) * t3 + g(2) * t1 + t22 * t45)];
taug = t6(:);
