% Calculate Gravitation load on the joints for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:04
% EndTime: 2019-12-31 18:01:05
% DurationCPUTime: 0.24s
% Computational Cost: add. (183->55), mult. (204->72), div. (0->0), fcn. (208->8), ass. (0->28)
t19 = sin(qJ(1));
t21 = cos(qJ(1));
t30 = pkin(8) + qJ(4);
t27 = sin(t30);
t28 = cos(t30);
t1 = -t19 * t27 - t21 * t28;
t2 = -t19 * t28 + t21 * t27;
t18 = sin(qJ(5));
t20 = cos(qJ(5));
t24 = -rSges(6,1) * t20 + rSges(6,2) * t18;
t22 = pkin(4) - t24;
t36 = rSges(6,3) + pkin(7);
t40 = -(g(1) * t36 - g(2) * t22) * t1 - (g(1) * t22 + g(2) * t36) * t2;
t17 = cos(pkin(8));
t12 = pkin(3) * t17 + pkin(2);
t14 = t21 * qJ(2);
t16 = sin(pkin(8));
t35 = t16 * t21;
t37 = pkin(3) * t35 + t14 + (-pkin(1) - t12) * t19;
t34 = t19 * t16;
t32 = t21 * pkin(1) + t19 * qJ(2);
t31 = m(4) + m(5) + m(6);
t29 = pkin(3) * t34 + t21 * t12 + t32;
t26 = -rSges(5,1) * t2 + rSges(5,2) * t1;
t25 = rSges(5,1) * t1 + rSges(5,2) * t2;
t4 = t17 * t21 + t34;
t3 = -t19 * t17 + t35;
t5 = [-m(2) * (g(1) * (-t19 * rSges(2,1) - rSges(2,2) * t21) + g(2) * (rSges(2,1) * t21 - t19 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t21 + t14 + (-rSges(3,1) - pkin(1)) * t19) + g(2) * (rSges(3,1) * t21 + t19 * rSges(3,3) + t32)) - m(4) * (g(1) * (rSges(4,1) * t3 + rSges(4,2) * t4 + t14 + (-pkin(1) - pkin(2)) * t19) + g(2) * (t4 * rSges(4,1) - t3 * rSges(4,2) + pkin(2) * t21 + t32)) - m(5) * (g(1) * (-t26 + t37) + g(2) * (-t25 + t29)) - m(6) * (g(1) * t37 + g(2) * t29 - t40), (-m(3) - t31) * (g(1) * t19 - g(2) * t21), t31 * g(3), -m(5) * (g(1) * t26 + g(2) * t25) - m(6) * t40, -m(6) * (g(3) * t24 + (g(1) * t1 + g(2) * t2) * (rSges(6,1) * t18 + rSges(6,2) * t20))];
taug = t5(:);
