% Calculate Gravitation load on the joints for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:35
% EndTime: 2019-12-31 16:33:36
% DurationCPUTime: 0.23s
% Computational Cost: add. (113->38), mult. (150->55), div. (0->0), fcn. (126->8), ass. (0->25)
t16 = sin(pkin(7));
t17 = cos(pkin(7));
t41 = g(1) * t17 + g(2) * t16;
t44 = rSges(5,3) + pkin(6);
t20 = cos(qJ(4));
t43 = -t20 * rSges(5,1) - pkin(3);
t15 = qJ(2) + qJ(3);
t12 = sin(t15);
t40 = t43 * t12;
t19 = sin(qJ(2));
t39 = pkin(2) * t19;
t18 = sin(qJ(4));
t36 = rSges(5,2) * t18;
t33 = t16 * t18;
t32 = t16 * t20;
t31 = t17 * t18;
t30 = t17 * t20;
t13 = cos(t15);
t26 = t13 * rSges(4,1) - rSges(4,2) * t12;
t25 = -rSges(4,1) * t12 - rSges(4,2) * t13;
t24 = t44 * t12 + (-t36 - t43) * t13;
t23 = t41 * (t12 * t36 + t13 * t44);
t21 = cos(qJ(2));
t14 = t21 * pkin(2);
t1 = [(-m(2) - m(3) - m(4) - m(5)) * g(3), -m(5) * t23 + (-m(3) * (rSges(3,1) * t21 - t19 * rSges(3,2)) - m(4) * (t14 + t26) - m(5) * (t14 + t24)) * g(3) + t41 * (-m(3) * (-rSges(3,1) * t19 - rSges(3,2) * t21) - m(4) * (t25 - t39) - m(5) * (-t39 + t40)), -m(4) * (g(3) * t26 + t41 * t25) - m(5) * (g(3) * t24 + t41 * t40 + t23), -m(5) * (g(1) * ((-t13 * t31 + t32) * rSges(5,1) + (-t13 * t30 - t33) * rSges(5,2)) + g(2) * ((-t13 * t33 - t30) * rSges(5,1) + (-t13 * t32 + t31) * rSges(5,2)) + g(3) * (-rSges(5,1) * t18 - rSges(5,2) * t20) * t12)];
taug = t1(:);
