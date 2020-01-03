% Calculate Gravitation load on the joints for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:09
% EndTime: 2019-12-31 17:20:11
% DurationCPUTime: 0.41s
% Computational Cost: add. (132->77), mult. (289->112), div. (0->0), fcn. (282->6), ass. (0->28)
t17 = sin(qJ(2));
t36 = g(3) * t17;
t35 = rSges(5,1) + pkin(3);
t18 = sin(qJ(1));
t21 = cos(qJ(1));
t40 = g(1) * t21 + g(2) * t18;
t28 = rSges(5,3) + qJ(4);
t39 = g(1) * t18;
t20 = cos(qJ(2));
t13 = t20 * pkin(2);
t34 = t17 * t21;
t33 = t18 * t20;
t32 = t20 * t21;
t16 = sin(qJ(3));
t31 = t21 * t16;
t19 = cos(qJ(3));
t30 = t21 * t19;
t29 = t21 * pkin(1) + t18 * pkin(5);
t27 = -pkin(1) - t13;
t26 = pkin(2) * t32 + pkin(6) * t34 + t29;
t25 = rSges(4,1) * t19 - rSges(4,2) * t16;
t23 = t28 * t16 + t35 * t19;
t14 = t21 * pkin(5);
t4 = t18 * t16 + t20 * t30;
t3 = -t18 * t19 + t20 * t31;
t2 = t19 * t33 - t31;
t1 = t16 * t33 + t30;
t5 = [-m(2) * (g(1) * (-t18 * rSges(2,1) - t21 * rSges(2,2)) + g(2) * (t21 * rSges(2,1) - t18 * rSges(2,2))) - m(3) * (g(1) * (t21 * rSges(3,3) + t14) + g(2) * (rSges(3,1) * t32 - rSges(3,2) * t34 + t29) + (g(1) * (-t20 * rSges(3,1) + t17 * rSges(3,2) - pkin(1)) + g(2) * rSges(3,3)) * t18) - m(4) * (g(1) * (-t2 * rSges(4,1) + t1 * rSges(4,2) + t14) + g(2) * (t4 * rSges(4,1) - t3 * rSges(4,2) + rSges(4,3) * t34 + t26) + ((-rSges(4,3) - pkin(6)) * t17 + t27) * t39) - m(5) * (g(1) * (-t28 * t1 - t35 * t2 + t14) + g(2) * (rSges(5,2) * t34 + t28 * t3 + t35 * t4 + t26) + ((-rSges(5,2) - pkin(6)) * t17 + t27) * t39), (-m(3) * (g(3) * rSges(3,1) - t40 * rSges(3,2)) - m(4) * (t40 * rSges(4,3) + g(3) * t25) - m(5) * (t40 * rSges(5,2) + g(3) * t23)) * t20 + ((m(3) * rSges(3,2) - m(4) * rSges(4,3) - m(5) * rSges(5,2)) * g(3) + t40 * (m(3) * rSges(3,1) - m(4) * (-pkin(2) - t25) - m(5) * (-pkin(2) - t23))) * t17 + (-m(4) - m(5)) * (g(3) * t13 + (g(1) * t32 + g(2) * t33 + t36) * pkin(6)), -m(4) * (g(1) * (-t3 * rSges(4,1) - t4 * rSges(4,2)) + g(2) * (-t1 * rSges(4,1) - t2 * rSges(4,2))) - m(5) * (g(1) * (t28 * t4 - t35 * t3) + g(2) * (-t35 * t1 + t28 * t2)) + (-m(4) * (-rSges(4,1) * t16 - rSges(4,2) * t19) - m(5) * (-t35 * t16 + t28 * t19)) * t36, -m(5) * (g(1) * t3 + g(2) * t1 + t16 * t36)];
taug = t5(:);
