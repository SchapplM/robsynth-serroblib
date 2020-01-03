% Calculate Gravitation load on the joints for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:42
% EndTime: 2019-12-31 18:40:42
% DurationCPUTime: 0.23s
% Computational Cost: add. (249->61), mult. (174->75), div. (0->0), fcn. (137->8), ass. (0->32)
t43 = rSges(6,1) + pkin(4);
t23 = sin(qJ(4));
t35 = rSges(6,3) + qJ(5);
t50 = t35 * t23;
t25 = cos(qJ(4));
t39 = t23 * rSges(5,2);
t49 = rSges(5,1) * t25 - t39;
t22 = qJ(1) + pkin(8);
t20 = qJ(3) + t22;
t16 = sin(t20);
t17 = cos(t20);
t48 = g(1) * t17 + g(2) * t16;
t47 = t43 * t25 + t50;
t24 = sin(qJ(1));
t46 = pkin(1) * t24;
t40 = t17 * t25;
t13 = t17 * pkin(7);
t38 = t17 * rSges(6,2) + t13;
t37 = t17 * pkin(3) + t16 * pkin(7);
t19 = cos(t22);
t26 = cos(qJ(1));
t21 = t26 * pkin(1);
t36 = pkin(2) * t19 + t21;
t34 = t17 * rSges(4,1) - rSges(4,2) * t16;
t18 = sin(t22);
t33 = -pkin(2) * t18 - t46;
t32 = -rSges(4,1) * t16 - rSges(4,2) * t17;
t31 = t16 * rSges(6,2) + t17 * t50 + t43 * t40 + t37;
t30 = rSges(5,1) * t40 + t16 * rSges(5,3) - t17 * t39 + t37;
t29 = t17 * rSges(5,3) + t13 + (-pkin(3) - t49) * t16;
t28 = g(1) * (-pkin(3) - t47) * t16;
t1 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - rSges(2,2) * t26) + g(2) * (rSges(2,1) * t26 - t24 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t18 - rSges(3,2) * t19 - t46) + g(2) * (rSges(3,1) * t19 - rSges(3,2) * t18 + t21)) - m(4) * (g(1) * (t32 + t33) + g(2) * (t34 + t36)) - m(5) * (g(1) * (t29 + t33) + g(2) * (t30 + t36)) - m(6) * (g(1) * (t33 + t38) + g(2) * (t31 + t36) + t28), (-m(3) - m(4) - m(5) - m(6)) * g(3), -m(4) * (g(1) * t32 + g(2) * t34) - m(5) * (g(1) * t29 + g(2) * t30) - m(6) * (g(1) * t38 + g(2) * t31 + t28), (-m(5) * t49 - m(6) * t47) * g(3) + t48 * (-m(5) * (-rSges(5,1) * t23 - rSges(5,2) * t25) - m(6) * (-t43 * t23 + t35 * t25)), -m(6) * (-g(3) * t25 + t48 * t23)];
taug = t1(:);
