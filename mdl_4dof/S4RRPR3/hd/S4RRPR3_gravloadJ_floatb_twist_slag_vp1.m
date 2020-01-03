% Calculate Gravitation load on the joints for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:28
% EndTime: 2019-12-31 17:01:28
% DurationCPUTime: 0.14s
% Computational Cost: add. (139->40), mult. (99->51), div. (0->0), fcn. (72->8), ass. (0->24)
t35 = rSges(5,3) + pkin(6);
t18 = sin(qJ(4));
t20 = cos(qJ(4));
t34 = rSges(5,1) * t20 - rSges(5,2) * t18;
t33 = -pkin(3) - t34;
t19 = sin(qJ(1));
t32 = pkin(1) * t19;
t17 = qJ(1) + qJ(2);
t14 = sin(t17);
t31 = pkin(2) * t14;
t15 = cos(t17);
t28 = t15 * rSges(3,1) - rSges(3,2) * t14;
t13 = pkin(7) + t17;
t10 = sin(t13);
t11 = cos(t13);
t12 = pkin(2) * t15;
t27 = t11 * rSges(4,1) - rSges(4,2) * t10 + t12;
t26 = -rSges(3,1) * t14 - rSges(3,2) * t15;
t24 = -rSges(4,1) * t10 - rSges(4,2) * t11 - t31;
t23 = t35 * t10 - t33 * t11 + t12;
t22 = t33 * t10 + t35 * t11 - t31;
t21 = cos(qJ(1));
t16 = t21 * pkin(1);
t1 = [-m(2) * (g(1) * (-t19 * rSges(2,1) - rSges(2,2) * t21) + g(2) * (rSges(2,1) * t21 - t19 * rSges(2,2))) - m(3) * (g(1) * (t26 - t32) + g(2) * (t16 + t28)) - m(4) * (g(1) * (t24 - t32) + g(2) * (t16 + t27)) - m(5) * (g(1) * (t22 - t32) + g(2) * (t16 + t23)), -m(3) * (g(1) * t26 + g(2) * t28) - m(4) * (g(1) * t24 + g(2) * t27) - m(5) * (g(1) * t22 + g(2) * t23), (-m(4) - m(5)) * g(3), -m(5) * (g(3) * t34 + (g(1) * t11 + g(2) * t10) * (-rSges(5,1) * t18 - rSges(5,2) * t20))];
taug = t1(:);
