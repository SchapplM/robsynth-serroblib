% Calculate Gravitation load on the joints for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:04
% EndTime: 2019-12-05 15:35:06
% DurationCPUTime: 0.33s
% Computational Cost: add. (162->58), mult. (222->88), div. (0->0), fcn. (203->8), ass. (0->32)
t39 = rSges(6,1) + pkin(4);
t18 = sin(qJ(4));
t20 = cos(qJ(4));
t46 = rSges(5,1) * t20 - rSges(5,2) * t18;
t16 = sin(pkin(7));
t17 = cos(pkin(7));
t45 = g(1) * t17 + g(2) * t16;
t30 = rSges(6,3) + qJ(5);
t44 = t30 * t18 + t39 * t20;
t19 = sin(qJ(2));
t43 = pkin(2) * t19;
t15 = qJ(2) + pkin(8);
t12 = sin(t15);
t40 = g(3) * t12;
t13 = cos(t15);
t36 = t13 * t16;
t35 = t13 * t17;
t34 = t16 * t18;
t33 = t16 * t20;
t32 = t17 * t18;
t31 = t17 * t20;
t29 = -m(4) - m(5) - m(6);
t21 = cos(qJ(2));
t14 = t21 * pkin(2);
t28 = t13 * pkin(3) + t12 * pkin(6) + t14;
t27 = pkin(6) * t36 - t16 * t43;
t26 = pkin(6) * t35 - t17 * t43;
t4 = t13 * t31 + t34;
t3 = t13 * t32 - t33;
t2 = t13 * t33 - t32;
t1 = t13 * t34 + t31;
t5 = [(-m(2) - m(3) + t29) * g(3), -m(3) * (g(3) * (t21 * rSges(3,1) - t19 * rSges(3,2)) + t45 * (-rSges(3,1) * t19 - rSges(3,2) * t21)) - m(4) * (g(3) * (t13 * rSges(4,1) - t12 * rSges(4,2) + t14) + t45 * (-rSges(4,1) * t12 - rSges(4,2) * t13 - t43)) - m(5) * (g(1) * (rSges(5,3) * t35 + t26) + g(2) * (rSges(5,3) * t36 + t27) + g(3) * (t46 * t13 + t28) + (g(3) * rSges(5,3) + t45 * (-pkin(3) - t46)) * t12) - m(6) * (g(1) * t26 + g(2) * t27 + g(3) * t28 + (t45 * rSges(6,2) + g(3) * t44) * t13 + (g(3) * rSges(6,2) + t45 * (-pkin(3) - t44)) * t12), t29 * (g(1) * t16 - g(2) * t17), -m(5) * (g(1) * (-rSges(5,1) * t3 - rSges(5,2) * t4) + g(2) * (-rSges(5,1) * t1 - rSges(5,2) * t2)) - m(6) * (g(1) * (-t39 * t3 + t30 * t4) + g(2) * (-t39 * t1 + t30 * t2)) + (-m(5) * (-rSges(5,1) * t18 - rSges(5,2) * t20) - m(6) * (-t39 * t18 + t30 * t20)) * t40, -m(6) * (g(1) * t3 + g(2) * t1 + t18 * t40)];
taug = t5(:);
