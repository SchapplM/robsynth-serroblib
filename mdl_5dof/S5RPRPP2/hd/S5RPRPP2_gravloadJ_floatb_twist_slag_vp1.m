% Calculate Gravitation load on the joints for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:40
% EndTime: 2019-12-31 18:10:41
% DurationCPUTime: 0.32s
% Computational Cost: add. (173->68), mult. (181->84), div. (0->0), fcn. (147->6), ass. (0->29)
t14 = qJ(1) + pkin(7);
t10 = cos(t14);
t9 = sin(t14);
t39 = g(1) * t10 + g(2) * t9;
t17 = cos(qJ(3));
t15 = sin(qJ(3));
t33 = rSges(6,2) * t15;
t34 = rSges(6,1) + pkin(4);
t40 = t17 * t34 + t33;
t37 = -m(5) - m(6);
t16 = sin(qJ(1));
t36 = pkin(1) * t16;
t12 = t17 * pkin(3);
t32 = t10 * t15;
t31 = t10 * t17;
t11 = t15 * qJ(4);
t30 = t11 + t12;
t28 = -rSges(6,3) - qJ(5);
t27 = -pkin(3) - t34;
t18 = cos(qJ(1));
t13 = t18 * pkin(1);
t26 = t10 * pkin(2) + t9 * pkin(6) + t13;
t25 = t10 * pkin(6) - t36;
t24 = -pkin(2) - t11;
t23 = pkin(3) * t31 + t10 * t11 + t26;
t22 = t39 * qJ(4) * t17;
t21 = rSges(4,1) * t17 - rSges(4,2) * t15;
t20 = rSges(5,1) * t17 + rSges(5,3) * t15;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t16 - rSges(2,2) * t18) + g(2) * (rSges(2,1) * t18 - rSges(2,2) * t16)) - m(3) * (g(1) * (-rSges(3,1) * t9 - rSges(3,2) * t10 - t36) + g(2) * (rSges(3,1) * t10 - rSges(3,2) * t9 + t13)) - m(4) * (g(1) * (rSges(4,3) * t10 + t25) + g(2) * (rSges(4,1) * t31 - rSges(4,2) * t32 + t26) + (g(1) * (-pkin(2) - t21) + g(2) * rSges(4,3)) * t9) - m(5) * (g(1) * (rSges(5,2) * t10 + t25) + g(2) * (rSges(5,1) * t31 + rSges(5,3) * t32 + t23) + (g(1) * (-t20 + t24 - t12) + g(2) * rSges(5,2)) * t9) - m(6) * (g(1) * t25 + g(2) * t23 + (g(1) * t28 + g(2) * t40) * t10 + (g(2) * t28 + (t17 * t27 + t24 - t33) * g(1)) * t9), (-m(3) - m(4) + t37) * g(3), -m(4) * g(3) * t21 - m(5) * (g(3) * (t20 + t30) + t22) - m(6) * (g(3) * (t30 + t40) + t22) + t39 * ((m(4) * rSges(4,2) - m(5) * rSges(5,3) - m(6) * rSges(6,2)) * t17 + (m(4) * rSges(4,1) - m(5) * (-rSges(5,1) - pkin(3)) - m(6) * t27) * t15), t37 * (-g(3) * t17 + t15 * t39), -m(6) * (-g(1) * t9 + g(2) * t10)];
taug = t1(:);
