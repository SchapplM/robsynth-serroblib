% Calculate Gravitation load on the joints for
% S4RRPR6
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:23
% EndTime: 2019-12-31 17:04:24
% DurationCPUTime: 0.24s
% Computational Cost: add. (127->46), mult. (139->60), div. (0->0), fcn. (109->8), ass. (0->26)
t14 = qJ(2) + pkin(7);
t10 = cos(t14);
t18 = cos(qJ(2));
t12 = t18 * pkin(2);
t11 = qJ(4) + t14;
t6 = sin(t11);
t7 = cos(t11);
t27 = t7 * rSges(5,1) - rSges(5,2) * t6;
t38 = pkin(3) * t10 + t12 + t27;
t9 = sin(t14);
t37 = rSges(4,1) * t10 - rSges(4,2) * t9 + t12;
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t36 = g(1) * t19 + g(2) * t17;
t16 = sin(qJ(2));
t34 = pkin(2) * t16;
t31 = rSges(3,3) + pkin(5);
t15 = -qJ(3) - pkin(5);
t29 = rSges(4,3) - t15;
t28 = rSges(5,3) + pkin(6) - t15;
t26 = -rSges(5,1) * t6 - rSges(5,2) * t7;
t24 = rSges(3,1) * t18 - rSges(3,2) * t16;
t23 = pkin(1) + t38;
t22 = pkin(1) + t37;
t21 = pkin(1) + t24;
t1 = [-m(2) * (g(1) * (-t17 * rSges(2,1) - rSges(2,2) * t19) + g(2) * (rSges(2,1) * t19 - t17 * rSges(2,2))) - m(3) * ((g(1) * t31 + g(2) * t21) * t19 + (-g(1) * t21 + g(2) * t31) * t17) - m(4) * ((g(1) * t29 + g(2) * t22) * t19 + (-g(1) * t22 + g(2) * t29) * t17) - m(5) * ((g(1) * t28 + g(2) * t23) * t19 + (-g(1) * t23 + g(2) * t28) * t17), (-m(3) * t24 - m(4) * t37 - m(5) * t38) * g(3) + t36 * (-m(3) * (-rSges(3,1) * t16 - rSges(3,2) * t18) - m(4) * (-rSges(4,1) * t9 - rSges(4,2) * t10 - t34) - m(5) * (-pkin(3) * t9 + t26 - t34)), (-m(4) - m(5)) * (g(1) * t17 - g(2) * t19), -m(5) * (g(3) * t27 + t36 * t26)];
taug = t1(:);
