% Calculate Gravitation load on the joints for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:07
% EndTime: 2019-12-31 21:49:08
% DurationCPUTime: 0.31s
% Computational Cost: add. (349->68), mult. (232->83), div. (0->0), fcn. (187->8), ass. (0->37)
t24 = sin(qJ(4));
t26 = cos(qJ(4));
t43 = rSges(6,3) + qJ(5);
t50 = rSges(6,1) + pkin(4);
t55 = t43 * t24 + t50 * t26;
t57 = t26 * rSges(5,1) - t24 * rSges(5,2);
t23 = qJ(1) + qJ(2);
t21 = qJ(3) + t23;
t17 = sin(t21);
t18 = cos(t21);
t56 = g(1) * t18 + g(2) * t17;
t19 = sin(t23);
t54 = pkin(2) * t19;
t25 = sin(qJ(1));
t51 = t25 * pkin(1);
t13 = t18 * pkin(8);
t45 = t18 * rSges(6,2) + t13;
t44 = t18 * pkin(3) + t17 * pkin(8);
t20 = cos(t23);
t42 = t20 * rSges(3,1) - t19 * rSges(3,2);
t41 = t18 * rSges(4,1) - t17 * rSges(4,2);
t40 = t45 - t54;
t16 = pkin(2) * t20;
t39 = t16 + t41;
t38 = -t19 * rSges(3,1) - t20 * rSges(3,2);
t37 = -t17 * rSges(4,1) - t18 * rSges(4,2);
t36 = t17 * rSges(6,2) + t55 * t18 + t44;
t35 = t16 + t36;
t34 = t17 * rSges(5,3) + t18 * t57 + t44;
t33 = t37 - t54;
t32 = t16 + t34;
t31 = t18 * rSges(5,3) + t13 + (-pkin(3) - t57) * t17;
t30 = t31 - t54;
t29 = g(1) * (-pkin(3) - t55) * t17;
t27 = cos(qJ(1));
t22 = t27 * pkin(1);
t1 = [-m(2) * (g(1) * (-t25 * rSges(2,1) - t27 * rSges(2,2)) + g(2) * (t27 * rSges(2,1) - t25 * rSges(2,2))) - m(3) * (g(1) * (t38 - t51) + g(2) * (t22 + t42)) - m(4) * (g(1) * (t33 - t51) + g(2) * (t22 + t39)) - m(5) * (g(1) * (t30 - t51) + g(2) * (t22 + t32)) - m(6) * (g(1) * (t40 - t51) + g(2) * (t22 + t35) + t29), -m(3) * (g(1) * t38 + g(2) * t42) - m(4) * (g(1) * t33 + g(2) * t39) - m(5) * (g(1) * t30 + g(2) * t32) - m(6) * (g(1) * t40 + g(2) * t35 + t29), -m(4) * (g(1) * t37 + g(2) * t41) - m(5) * (g(1) * t31 + g(2) * t34) - m(6) * (g(1) * t45 + g(2) * t36 + t29), (-m(5) * t57 - m(6) * t55) * g(3) + t56 * (-m(5) * (-rSges(5,1) * t24 - rSges(5,2) * t26) - m(6) * (-t50 * t24 + t43 * t26)), -m(6) * (-g(3) * t26 + t56 * t24)];
taug = t1(:);
