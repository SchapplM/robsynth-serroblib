% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:59
% EndTime: 2022-01-20 10:20:00
% DurationCPUTime: 0.10s
% Computational Cost: add. (149->32), mult. (110->35), div. (0->0), fcn. (99->8), ass. (0->27)
t17 = qJ(1) + qJ(2);
t14 = sin(t17);
t28 = pkin(2) * t14;
t20 = sin(qJ(1));
t27 = t20 * pkin(1);
t13 = pkin(8) + t17;
t10 = cos(t13);
t15 = cos(t17);
t11 = pkin(2) * t15;
t21 = cos(qJ(4));
t12 = t21 * pkin(4) + pkin(3);
t18 = -qJ(5) - pkin(7);
t9 = sin(t13);
t26 = t10 * t12 - t9 * t18 + t11;
t25 = g(1) * t9 - g(2) * t10;
t24 = g(1) * t10 + g(2) * t9;
t6 = g(1) * t14 - g(2) * t15;
t23 = -t10 * t18 - t9 * t12 - t28;
t19 = sin(qJ(4));
t1 = -g(3) * t21 + t24 * t19;
t22 = cos(qJ(1));
t16 = t22 * pkin(1);
t7 = g(1) * t15 + g(2) * t14;
t4 = t25 * t21;
t3 = t25 * t19;
t2 = g(3) * t19 + t24 * t21;
t5 = [0, g(1) * t20 - g(2) * t22, g(1) * t22 + g(2) * t20, 0, t6, t7, -g(1) * (-t27 - t28) - g(2) * (t11 + t16), 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t24, -g(1) * (t23 - t27) - g(2) * (t16 + t26); 0, 0, 0, 0, t6, t7, t6 * pkin(2), 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t24, -g(1) * t23 - g(2) * t26; 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25;];
taug_reg = t5;
