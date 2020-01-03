% Calculate Gravitation load on the joints for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP13_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP13_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:32
% EndTime: 2019-12-31 18:58:33
% DurationCPUTime: 0.45s
% Computational Cost: add. (153->92), mult. (318->127), div. (0->0), fcn. (304->6), ass. (0->41)
t50 = -m(5) - m(6);
t42 = rSges(6,1) + pkin(4);
t30 = rSges(6,3) + qJ(5);
t49 = -pkin(1) - pkin(6);
t48 = m(4) * rSges(4,1);
t47 = m(4) * rSges(4,2);
t18 = sin(qJ(1));
t46 = g(1) * t18;
t20 = cos(qJ(3));
t45 = g(2) * t20;
t21 = cos(qJ(1));
t44 = g(2) * t21;
t43 = g(3) * t20;
t41 = -rSges(6,2) - pkin(7);
t40 = -rSges(5,3) - pkin(7);
t17 = sin(qJ(3));
t39 = t17 * t18;
t38 = t17 * t21;
t16 = sin(qJ(4));
t37 = t18 * t16;
t19 = cos(qJ(4));
t36 = t18 * t19;
t35 = t20 * rSges(4,2);
t34 = t20 * t21;
t33 = t21 * t16;
t32 = t21 * t19;
t31 = t21 * pkin(1) + t18 * qJ(2);
t29 = t21 * pkin(6) + t31;
t28 = g(1) * t49;
t26 = pkin(3) * t39 + t29;
t25 = rSges(5,1) * t19 - rSges(5,2) * t16;
t12 = t21 * qJ(2);
t24 = pkin(3) * t38 - pkin(7) * t34 + t12;
t23 = t30 * t16 + t42 * t19;
t22 = t48 - m(5) * (-pkin(3) - t25) - m(6) * (-pkin(3) - t23);
t13 = t20 * pkin(7);
t4 = t17 * t32 - t37;
t3 = t17 * t33 + t36;
t2 = t17 * t36 + t33;
t1 = t17 * t37 - t32;
t5 = [-m(2) * (g(1) * (-t18 * rSges(2,1) - t21 * rSges(2,2)) + g(2) * (t21 * rSges(2,1) - t18 * rSges(2,2))) - m(3) * (g(1) * (t21 * rSges(3,3) + t12 + (rSges(3,2) - pkin(1)) * t18) + g(2) * (-t21 * rSges(3,2) + t18 * rSges(3,3) + t31)) - m(4) * (g(1) * (rSges(4,1) * t38 + rSges(4,2) * t34 + t12) + g(2) * (t21 * rSges(4,3) + t29) + (g(1) * (-rSges(4,3) + t49) + g(2) * (t17 * rSges(4,1) + t35)) * t18) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) - rSges(5,3) * t34 + t24) + g(2) * (t2 * rSges(5,1) - t1 * rSges(5,2) + t26) + (t40 * t45 + t28) * t18) - m(6) * (g(1) * (-rSges(6,2) * t34 + t30 * t3 + t42 * t4 + t24) + g(2) * (t30 * t1 + t42 * t2 + t26) + (t41 * t45 + t28) * t18), (-m(3) - m(4) + t50) * (-t44 + t46), ((-m(5) * rSges(5,3) - m(6) * rSges(6,2) + t47) * t17 + (-m(5) * t25 - m(6) * t23 - t48) * t20) * t46 + ((-m(5) * t40 - m(6) * t41 - t47) * t17 + t22 * t20) * t44 + t50 * g(1) * (t18 * t20 * pkin(3) + pkin(7) * t39) + (m(4) * t35 - m(5) * (t20 * rSges(5,3) + t13) - m(6) * (t20 * rSges(6,2) + t13) + t22 * t17) * g(3), -m(5) * (g(1) * (-t1 * rSges(5,1) - t2 * rSges(5,2)) + g(2) * (t3 * rSges(5,1) + t4 * rSges(5,2))) - m(6) * (g(1) * (-t42 * t1 + t30 * t2) + g(2) * (t42 * t3 - t30 * t4)) + (-m(5) * (-rSges(5,1) * t16 - rSges(5,2) * t19) - m(6) * (-t42 * t16 + t30 * t19)) * t43, -m(6) * (g(1) * t1 - g(2) * t3 + t16 * t43)];
taug = t5(:);
