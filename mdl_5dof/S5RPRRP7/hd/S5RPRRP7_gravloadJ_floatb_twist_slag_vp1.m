% Calculate Gravitation load on the joints for
% S5RPRRP7
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:33
% EndTime: 2019-12-31 18:44:35
% DurationCPUTime: 0.49s
% Computational Cost: add. (261->87), mult. (308->122), div. (0->0), fcn. (294->8), ass. (0->35)
t49 = -m(5) - m(6);
t21 = sin(qJ(3));
t43 = g(3) * t21;
t19 = qJ(1) + pkin(8);
t14 = sin(t19);
t44 = g(2) * t14;
t41 = rSges(6,1) + pkin(4);
t24 = cos(qJ(3));
t48 = t24 * rSges(4,1) - t21 * rSges(4,2);
t15 = cos(t19);
t47 = g(1) * t15 + t44;
t34 = rSges(6,3) + qJ(5);
t46 = g(1) * t14;
t22 = sin(qJ(1));
t42 = t22 * pkin(1);
t17 = t24 * pkin(3);
t40 = t15 * t21;
t39 = t15 * t24;
t20 = sin(qJ(4));
t38 = t20 * t24;
t23 = cos(qJ(4));
t36 = t23 * t24;
t25 = cos(qJ(1));
t18 = t25 * pkin(1);
t33 = t15 * pkin(2) + t14 * pkin(6) + t18;
t32 = -pkin(2) - t17;
t31 = t15 * pkin(6) - t42;
t30 = pkin(3) * t39 + pkin(7) * t40 + t33;
t29 = rSges(5,1) * t23 - rSges(5,2) * t20;
t27 = t34 * t20 + t41 * t23;
t4 = t14 * t20 + t15 * t36;
t3 = -t14 * t23 + t15 * t38;
t2 = t14 * t36 - t15 * t20;
t1 = t14 * t38 + t15 * t23;
t5 = [-m(2) * (g(1) * (-t22 * rSges(2,1) - t25 * rSges(2,2)) + g(2) * (t25 * rSges(2,1) - t22 * rSges(2,2))) - m(3) * (g(1) * (-t14 * rSges(3,1) - t15 * rSges(3,2) - t42) + g(2) * (t15 * rSges(3,1) - t14 * rSges(3,2) + t18)) - m(4) * (g(1) * (t15 * rSges(4,3) + t31) + g(2) * (t48 * t15 + t33) + (g(1) * (-pkin(2) - t48) + g(2) * rSges(4,3)) * t14) - m(5) * (g(1) * (-t2 * rSges(5,1) + t1 * rSges(5,2) + t31) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + rSges(5,3) * t40 + t30) + ((-rSges(5,3) - pkin(7)) * t21 + t32) * t46) - m(6) * (g(1) * (-t34 * t1 - t41 * t2 + t31) + g(2) * (rSges(6,2) * t40 + t3 * t34 + t4 * t41 + t30) + ((-rSges(6,2) - pkin(7)) * t21 + t32) * t46), (-m(3) - m(4) + t49) * g(3), (-m(4) * (g(3) * rSges(4,1) - rSges(4,2) * t47) - m(5) * (rSges(5,3) * t47 + g(3) * t29) - m(6) * (rSges(6,2) * t47 + g(3) * t27)) * t24 + ((m(4) * rSges(4,2) - m(5) * rSges(5,3) - m(6) * rSges(6,2)) * g(3) + t47 * (m(4) * rSges(4,1) - m(5) * (-pkin(3) - t29) - m(6) * (-pkin(3) - t27))) * t21 + t49 * (g(3) * t17 + (g(1) * t39 + t24 * t44 + t43) * pkin(7)), -m(5) * (g(1) * (-t3 * rSges(5,1) - t4 * rSges(5,2)) + g(2) * (-t1 * rSges(5,1) - t2 * rSges(5,2))) - m(6) * (g(1) * (-t3 * t41 + t34 * t4) + g(2) * (-t1 * t41 + t2 * t34)) + (-m(5) * (-rSges(5,1) * t20 - rSges(5,2) * t23) - m(6) * (-t20 * t41 + t34 * t23)) * t43, -m(6) * (g(1) * t3 + g(2) * t1 + t20 * t43)];
taug = t5(:);
