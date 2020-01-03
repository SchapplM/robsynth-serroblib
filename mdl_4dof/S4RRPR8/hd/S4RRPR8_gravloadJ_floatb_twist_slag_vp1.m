% Calculate Gravitation load on the joints for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:48
% EndTime: 2019-12-31 17:07:49
% DurationCPUTime: 0.29s
% Computational Cost: add. (104->68), mult. (217->96), div. (0->0), fcn. (203->6), ass. (0->32)
t43 = -rSges(5,3) - pkin(6);
t19 = sin(qJ(2));
t13 = t19 * qJ(3);
t22 = cos(qJ(2));
t36 = t22 * pkin(2) + t13;
t20 = sin(qJ(1));
t23 = cos(qJ(1));
t42 = g(1) * t23 + g(2) * t20;
t41 = pkin(3) * t22;
t38 = t19 * t23;
t37 = t22 * t23;
t35 = t23 * pkin(1) + t20 * pkin(5);
t34 = qJ(3) * t22;
t33 = pkin(2) * t37 + t23 * t13 + t35;
t18 = sin(qJ(4));
t21 = cos(qJ(4));
t27 = t18 * t22 - t19 * t21;
t2 = t27 * t20;
t26 = t18 * t19 + t21 * t22;
t3 = t26 * t20;
t32 = -rSges(5,1) * t2 - rSges(5,2) * t3;
t4 = t18 * t37 - t21 * t38;
t5 = t26 * t23;
t31 = -t4 * rSges(5,1) - t5 * rSges(5,2);
t30 = -rSges(5,1) * t26 + rSges(5,2) * t27;
t29 = rSges(3,1) * t22 - rSges(3,2) * t19;
t28 = rSges(4,1) * t22 + rSges(4,3) * t19;
t25 = -pkin(1) - t36;
t16 = t23 * pkin(5);
t11 = t23 * t34;
t9 = t20 * t34;
t1 = [-m(2) * (g(1) * (-t20 * rSges(2,1) - rSges(2,2) * t23) + g(2) * (rSges(2,1) * t23 - t20 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t23 + t16) + g(2) * (rSges(3,1) * t37 - rSges(3,2) * t38 + t35) + (g(1) * (-pkin(1) - t29) + g(2) * rSges(3,3)) * t20) - m(4) * (g(1) * (rSges(4,2) * t23 + t16) + g(2) * (rSges(4,1) * t37 + rSges(4,3) * t38 + t33) + (g(1) * (t25 - t28) + g(2) * rSges(4,2)) * t20) - m(5) * (g(1) * (-t3 * rSges(5,1) + t2 * rSges(5,2) + t43 * t23 + t16) + g(2) * (t5 * rSges(5,1) - t4 * rSges(5,2) + pkin(3) * t37 + t33) + (g(1) * (t25 - t41) + g(2) * t43) * t20), -m(3) * g(3) * t29 - m(4) * (g(1) * t11 + g(2) * t9 + g(3) * (t28 + t36)) - m(5) * (g(1) * (t11 - t31) + g(2) * (-t32 + t9) + g(3) * (-t30 + t36 + t41)) + t42 * ((m(3) * rSges(3,2) - m(4) * rSges(4,3)) * t22 + (m(3) * rSges(3,1) - m(4) * (-rSges(4,1) - pkin(2)) - m(5) * (-pkin(2) - pkin(3))) * t19), (-m(4) - m(5)) * (-g(3) * t22 + t42 * t19), -m(5) * (g(1) * t31 + g(2) * t32 + g(3) * t30)];
taug = t1(:);
