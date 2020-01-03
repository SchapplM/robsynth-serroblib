% Calculate Gravitation load on the joints for
% S4RRPR4
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:26
% EndTime: 2019-12-31 17:02:27
% DurationCPUTime: 0.17s
% Computational Cost: add. (143->41), mult. (119->54), div. (0->0), fcn. (92->8), ass. (0->23)
t44 = pkin(6) + qJ(3) + rSges(5,3);
t21 = pkin(7) + qJ(4);
t16 = sin(t21);
t17 = cos(t21);
t43 = rSges(5,1) * t17 - rSges(5,2) * t16;
t42 = rSges(4,3) + qJ(3);
t24 = cos(pkin(7));
t41 = rSges(4,2) * sin(pkin(7)) - pkin(2) - rSges(4,1) * t24;
t40 = -pkin(3) * t24 - pkin(2) - t43;
t26 = sin(qJ(1));
t39 = pkin(1) * t26;
t22 = qJ(1) + qJ(2);
t18 = sin(t22);
t19 = cos(t22);
t34 = t19 * rSges(3,1) - rSges(3,2) * t18;
t33 = -rSges(3,1) * t18 - rSges(3,2) * t19;
t31 = t18 * t42 - t19 * t41;
t30 = t18 * t41 + t19 * t42;
t29 = t44 * t18 - t40 * t19;
t28 = t40 * t18 + t44 * t19;
t27 = cos(qJ(1));
t20 = t27 * pkin(1);
t1 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - rSges(2,2) * t27) + g(2) * (rSges(2,1) * t27 - t26 * rSges(2,2))) - m(3) * (g(1) * (t33 - t39) + g(2) * (t20 + t34)) - m(4) * (g(1) * (t30 - t39) + g(2) * (t20 + t31)) - m(5) * (g(1) * (t28 - t39) + g(2) * (t20 + t29)), -m(3) * (g(1) * t33 + g(2) * t34) - m(4) * (g(1) * t30 + g(2) * t31) - m(5) * (g(1) * t28 + g(2) * t29), (-m(4) - m(5)) * (g(1) * t18 - g(2) * t19), -m(5) * (g(3) * t43 + (g(1) * t19 + g(2) * t18) * (-rSges(5,1) * t16 - rSges(5,2) * t17))];
taug = t1(:);
