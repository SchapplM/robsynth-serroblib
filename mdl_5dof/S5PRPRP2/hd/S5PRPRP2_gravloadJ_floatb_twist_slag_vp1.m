% Calculate Gravitation load on the joints for
% S5PRPRP2
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:31
% EndTime: 2019-12-05 15:30:33
% DurationCPUTime: 0.28s
% Computational Cost: add. (171->63), mult. (180->89), div. (0->0), fcn. (167->6), ass. (0->24)
t32 = rSges(6,1) + pkin(4);
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t31 = rSges(4,1) * t14 - rSges(4,2) * t13;
t30 = pkin(3) * t14 + (rSges(5,3) + pkin(6)) * t13;
t12 = pkin(7) + qJ(2);
t10 = sin(t12);
t11 = cos(t12);
t29 = t11 * pkin(2) + t10 * qJ(3);
t16 = sin(qJ(4));
t27 = pkin(4) * t16;
t26 = g(1) * t10;
t25 = g(2) * t11;
t21 = t14 * t16;
t17 = cos(qJ(4));
t20 = t14 * t17;
t19 = -m(4) - m(5) - m(6);
t3 = t10 * t17 - t11 * t21;
t1 = t10 * t21 + t11 * t17;
t18 = t14 * (pkin(4) * t17 + pkin(3)) + (rSges(6,3) + qJ(5) + pkin(6)) * t13;
t7 = t11 * qJ(3);
t4 = t10 * t16 + t11 * t20;
t2 = -t10 * t20 + t11 * t16;
t5 = [(-m(2) - m(3) + t19) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t10 - rSges(3,2) * t11) + g(2) * (rSges(3,1) * t11 - rSges(3,2) * t10)) - m(4) * (g(1) * (rSges(4,3) * t11 + t7) + g(2) * (t31 * t11 + t29) + (g(1) * (-pkin(2) - t31) + g(2) * rSges(4,3)) * t10) - m(5) * (g(1) * (rSges(5,1) * t2 + rSges(5,2) * t1 + t7) + g(2) * (rSges(5,1) * t4 + rSges(5,2) * t3 + t29) + t30 * t25 + (-pkin(2) - t30) * t26) - m(6) * (g(1) * (rSges(6,1) * t2 + rSges(6,2) * t1 + t7) + g(2) * (rSges(6,1) * t4 + rSges(6,2) * t3 + t29) + (g(1) * t27 + g(2) * t18) * t11 + (g(1) * (-pkin(2) - t18) + g(2) * t27) * t10), t19 * (-t25 + t26), -m(5) * (g(1) * (rSges(5,1) * t3 - rSges(5,2) * t4) + g(2) * (-rSges(5,1) * t1 + rSges(5,2) * t2)) - m(6) * (g(1) * (-t4 * rSges(6,2) + t32 * t3) + g(2) * (t2 * rSges(6,2) - t32 * t1)) + (-m(5) * (-rSges(5,1) * t16 - rSges(5,2) * t17) - m(6) * (-rSges(6,1) * t16 - rSges(6,2) * t17 - t27)) * g(3) * t13, -m(6) * (-g(3) * t14 + (g(1) * t11 + g(2) * t10) * t13)];
taug = t5(:);
