% Calculate Gravitation load on the joints for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:10
% EndTime: 2019-12-31 17:00:11
% DurationCPUTime: 0.30s
% Computational Cost: add. (85->57), mult. (165->74), div. (0->0), fcn. (138->4), ass. (0->24)
t13 = cos(qJ(2));
t12 = sin(qJ(1));
t14 = cos(qJ(1));
t18 = g(1) * t14 + g(2) * t12;
t35 = t18 * t13;
t23 = rSges(5,3) + qJ(4);
t11 = sin(qJ(2));
t27 = rSges(5,2) * t11;
t34 = t23 * t13 + t27;
t6 = t11 * qJ(3);
t8 = t13 * pkin(2);
t32 = t6 + t8;
t29 = rSges(5,1) + pkin(3);
t28 = t14 * pkin(1) + t12 * pkin(5);
t26 = t11 * t14;
t25 = t13 * t14;
t22 = -pkin(2) - t23;
t21 = pkin(2) * t25 + t14 * t6 + t28;
t20 = -pkin(1) - t6;
t19 = qJ(3) * t35;
t17 = rSges(3,1) * t13 - rSges(3,2) * t11;
t16 = -rSges(4,2) * t13 + rSges(4,3) * t11;
t9 = t14 * pkin(5);
t1 = [-m(2) * (g(1) * (-t12 * rSges(2,1) - rSges(2,2) * t14) + g(2) * (rSges(2,1) * t14 - t12 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t14 + t9) + g(2) * (rSges(3,1) * t25 - rSges(3,2) * t26 + t28) + (g(1) * (-pkin(1) - t17) + g(2) * rSges(3,3)) * t12) - m(4) * (g(1) * (rSges(4,1) * t14 + t9) + g(2) * (-rSges(4,2) * t25 + rSges(4,3) * t26 + t21) + (g(1) * (-t16 + t20 - t8) + g(2) * rSges(4,1)) * t12) - m(5) * (g(1) * t9 + g(2) * t21 + (g(1) * t29 + g(2) * t34) * t14 + (g(2) * t29 + (t22 * t13 + t20 - t27) * g(1)) * t12), -m(3) * g(3) * t17 - m(4) * (g(3) * (t16 + t32) + t19) - m(5) * (g(3) * (t32 + t34) + t19) + t18 * ((m(3) * rSges(3,2) - m(4) * rSges(4,3) - m(5) * rSges(5,2)) * t13 + (m(3) * rSges(3,1) - m(4) * (rSges(4,2) - pkin(2)) - m(5) * t22) * t11), (-m(4) - m(5)) * (-g(3) * t13 + t18 * t11), -m(5) * (g(3) * t11 + t35)];
taug = t1(:);
