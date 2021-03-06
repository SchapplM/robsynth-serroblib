% Calculate Gravitation load on the joints for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:51
% EndTime: 2019-12-31 17:43:52
% DurationCPUTime: 0.24s
% Computational Cost: add. (154->61), mult. (162->87), div. (0->0), fcn. (144->8), ass. (0->28)
t34 = -rSges(6,3) - pkin(6);
t33 = -m(5) - m(6);
t19 = sin(qJ(1));
t32 = pkin(1) * t19;
t15 = qJ(1) + pkin(7);
t13 = cos(t15);
t16 = sin(pkin(8));
t31 = t13 * t16;
t17 = cos(pkin(8));
t30 = t13 * t17;
t29 = qJ(4) * t16;
t28 = -m(4) + t33;
t12 = sin(t15);
t21 = cos(qJ(1));
t14 = t21 * pkin(1);
t27 = t13 * pkin(2) + t12 * qJ(3) + t14;
t26 = t13 * qJ(3) - t32;
t25 = pkin(3) * t30 + t13 * t29 + t27;
t18 = sin(qJ(5));
t20 = cos(qJ(5));
t24 = t16 * t20 - t17 * t18;
t23 = t16 * t18 + t17 * t20;
t22 = -pkin(3) * t17 - pkin(2) - t29;
t5 = t23 * t13;
t4 = t24 * t13;
t3 = t23 * t12;
t2 = t24 * t12;
t1 = [-m(2) * (g(1) * (-t19 * rSges(2,1) - rSges(2,2) * t21) + g(2) * (rSges(2,1) * t21 - t19 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t12 - rSges(3,2) * t13 - t32) + g(2) * (rSges(3,1) * t13 - rSges(3,2) * t12 + t14)) - m(4) * (g(1) * (rSges(4,3) * t13 + t26) + g(2) * (rSges(4,1) * t30 - rSges(4,2) * t31 + t27) + (g(1) * (-rSges(4,1) * t17 + rSges(4,2) * t16 - pkin(2)) + g(2) * rSges(4,3)) * t12) - m(5) * (g(1) * (rSges(5,2) * t13 + t26) + g(2) * (rSges(5,1) * t30 + rSges(5,3) * t31 + t25) + (g(1) * (-rSges(5,1) * t17 - rSges(5,3) * t16 + t22) + g(2) * rSges(5,2)) * t12) - m(6) * (g(1) * (-rSges(6,1) * t3 - rSges(6,2) * t2 + t34 * t13 + t26) + g(2) * (rSges(6,1) * t5 + rSges(6,2) * t4 + pkin(4) * t30 + t25) + (g(1) * (-pkin(4) * t17 + t22) + g(2) * t34) * t12), (-m(3) + t28) * g(3), t28 * (g(1) * t12 - g(2) * t13), t33 * (-g(3) * t17 + (g(1) * t13 + g(2) * t12) * t16), -m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t5) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t3) + g(3) * (-t23 * rSges(6,1) - t24 * rSges(6,2)))];
taug = t1(:);
