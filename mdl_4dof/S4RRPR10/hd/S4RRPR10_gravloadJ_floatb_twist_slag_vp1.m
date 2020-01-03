% Calculate Gravitation load on the joints for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:10
% EndTime: 2019-12-31 17:11:11
% DurationCPUTime: 0.31s
% Computational Cost: add. (99->68), mult. (204->93), div. (0->0), fcn. (182->6), ass. (0->31)
t17 = sin(qJ(1));
t20 = cos(qJ(1));
t41 = g(1) * t20 + g(2) * t17;
t37 = rSges(5,3) + pkin(6);
t19 = cos(qJ(2));
t38 = g(3) * t19;
t12 = t19 * pkin(2);
t16 = sin(qJ(2));
t36 = t16 * t20;
t15 = sin(qJ(4));
t35 = t17 * t15;
t18 = cos(qJ(4));
t34 = t17 * t18;
t33 = t18 * t20;
t32 = t19 * t20;
t10 = t16 * qJ(3);
t31 = t10 + t12;
t30 = t20 * pkin(1) + t17 * pkin(5);
t28 = -pkin(2) - t37;
t27 = pkin(2) * t32 + t20 * t10 + t30;
t26 = -pkin(1) - t10;
t25 = t41 * qJ(3) * t19;
t24 = rSges(3,1) * t19 - rSges(3,2) * t16;
t23 = rSges(5,1) * t15 + rSges(5,2) * t18;
t22 = -rSges(4,2) * t19 + rSges(4,3) * t16;
t13 = t20 * pkin(5);
t5 = -t16 * t35 + t33;
t4 = t15 * t20 + t16 * t34;
t3 = t15 * t36 + t34;
t2 = t16 * t33 - t35;
t1 = [-m(2) * (g(1) * (-t17 * rSges(2,1) - rSges(2,2) * t20) + g(2) * (rSges(2,1) * t20 - t17 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t20 + t13) + g(2) * (rSges(3,1) * t32 - rSges(3,2) * t36 + t30) + (g(1) * (-pkin(1) - t24) + g(2) * rSges(3,3)) * t17) - m(4) * (g(1) * (rSges(4,1) * t20 + t13) + g(2) * (-rSges(4,2) * t32 + rSges(4,3) * t36 + t27) + (g(1) * (-t22 + t26 - t12) + g(2) * rSges(4,1)) * t17) - m(5) * ((t3 * rSges(5,1) + t2 * rSges(5,2) + t17 * pkin(3) + t37 * t32 + t27) * g(2) + (t5 * rSges(5,1) - t4 * rSges(5,2) + pkin(3) * t20 + t13 + (t19 * t28 + t26) * t17) * g(1)), -m(3) * g(3) * t24 - m(4) * (g(3) * (t22 + t31) + t25) - m(5) * (g(3) * (t23 * t16 + t37 * t19 + t31) + t25) + t41 * ((m(3) * rSges(3,2) - m(4) * rSges(4,3) - m(5) * t23) * t19 + (m(3) * rSges(3,1) - m(4) * (rSges(4,2) - pkin(2)) - m(5) * t28) * t16), (-m(4) - m(5)) * (t41 * t16 - t38), -m(5) * (g(1) * (rSges(5,1) * t2 - rSges(5,2) * t3) + g(2) * (rSges(5,1) * t4 + rSges(5,2) * t5) + (-rSges(5,1) * t18 + rSges(5,2) * t15) * t38)];
taug = t1(:);
