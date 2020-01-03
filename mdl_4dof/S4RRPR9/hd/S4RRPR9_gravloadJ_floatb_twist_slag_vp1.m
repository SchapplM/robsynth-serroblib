% Calculate Gravitation load on the joints for
% S4RRPR9
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
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:21
% EndTime: 2019-12-31 17:09:22
% DurationCPUTime: 0.33s
% Computational Cost: add. (135->68), mult. (223->102), div. (0->0), fcn. (205->8), ass. (0->27)
t36 = rSges(5,3) + pkin(6) + qJ(3);
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t35 = g(1) * t19 + g(2) * t17;
t13 = sin(pkin(7));
t14 = cos(pkin(7));
t6 = pkin(3) * t14 + pkin(2);
t12 = pkin(7) + qJ(4);
t7 = sin(t12);
t8 = cos(t12);
t34 = m(3) * rSges(3,1) + m(4) * (rSges(4,1) * t14 - rSges(4,2) * t13 + pkin(2)) + m(5) * (rSges(5,1) * t8 - rSges(5,2) * t7 + t6);
t32 = pkin(3) * t13;
t29 = t19 * pkin(1) + t17 * pkin(5);
t16 = sin(qJ(2));
t28 = rSges(3,2) * t16;
t18 = cos(qJ(2));
t27 = t17 * t18;
t26 = t19 * t18;
t25 = rSges(4,3) + qJ(3);
t22 = t16 * t36 + t18 * t6;
t21 = m(3) * rSges(3,2) - m(4) * t25 - m(5) * t36;
t10 = t19 * pkin(5);
t4 = t17 * t7 + t26 * t8;
t3 = t17 * t8 - t26 * t7;
t2 = t19 * t7 - t27 * t8;
t1 = t19 * t8 + t27 * t7;
t5 = [-m(2) * (g(1) * (-rSges(2,1) * t17 - rSges(2,2) * t19) + g(2) * (rSges(2,1) * t19 - rSges(2,2) * t17)) - m(3) * (g(1) * (rSges(3,3) * t19 + t10) + g(2) * (rSges(3,1) * t26 - t19 * t28 + t29) + (g(1) * (-rSges(3,1) * t18 - pkin(1) + t28) + g(2) * rSges(3,3)) * t17) - m(4) * (g(1) * (-pkin(2) * t27 - t17 * pkin(1) + t10 + (t13 * t19 - t14 * t27) * rSges(4,1) + (t13 * t27 + t14 * t19) * rSges(4,2)) + g(2) * (pkin(2) * t26 + (t13 * t17 + t14 * t26) * rSges(4,1) + (-t13 * t26 + t14 * t17) * rSges(4,2) + t29) + (-g(1) * t17 + g(2) * t19) * t16 * t25) - m(5) * (g(1) * (t2 * rSges(5,1) + t1 * rSges(5,2) + t10) + g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t29) + (g(1) * t32 + g(2) * t22) * t19 + (g(1) * (-pkin(1) - t22) + g(2) * t32) * t17), (t21 * t16 - t18 * t34) * g(3) + t35 * (t16 * t34 + t21 * t18), (-m(4) - m(5)) * (-g(3) * t18 + t16 * t35), -m(5) * (g(1) * (rSges(5,1) * t3 - rSges(5,2) * t4) + g(2) * (-rSges(5,1) * t1 + rSges(5,2) * t2) + g(3) * (-rSges(5,1) * t7 - rSges(5,2) * t8) * t16)];
taug = t5(:);
