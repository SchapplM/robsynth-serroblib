% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:42
% EndTime: 2022-01-20 09:51:42
% DurationCPUTime: 0.16s
% Computational Cost: add. (214->45), mult. (128->45), div. (0->0), fcn. (114->10), ass. (0->32)
t25 = qJ(1) + qJ(2);
t21 = sin(t25);
t37 = pkin(2) * t21;
t29 = sin(qJ(1));
t36 = t29 * pkin(1);
t20 = pkin(8) + t25;
t14 = sin(t20);
t15 = cos(t20);
t22 = cos(t25);
t17 = pkin(2) * t22;
t35 = t15 * pkin(3) + t14 * qJ(4) + t17;
t27 = cos(pkin(9));
t16 = t27 * pkin(4) + pkin(3);
t28 = -pkin(7) - qJ(4);
t34 = -t14 * t28 + t15 * t16 + t17;
t6 = g(1) * t15 + g(2) * t14;
t5 = g(1) * t14 - g(2) * t15;
t7 = g(1) * t21 - g(2) * t22;
t30 = cos(qJ(1));
t33 = g(1) * t29 - g(2) * t30;
t32 = -t14 * pkin(3) + t15 * qJ(4) - t37;
t31 = -t14 * t16 - t15 * t28 - t37;
t24 = pkin(9) + qJ(5);
t23 = t30 * pkin(1);
t19 = cos(t24);
t18 = sin(t24);
t8 = g(1) * t22 + g(2) * t21;
t4 = t5 * t27;
t3 = t5 * sin(pkin(9));
t2 = t5 * t19;
t1 = t5 * t18;
t9 = [0, 0, 0, 0, 0, 0, t33, g(1) * t30 + g(2) * t29, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t33 * pkin(1), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(1) * (-t36 - t37) - g(2) * (t17 + t23), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t32 - t36) - g(2) * (t23 + t35), 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(1) * (t31 - t36) - g(2) * (t23 + t34); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t7 * pkin(2), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t32 - g(2) * t35, 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(1) * t31 - g(2) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t19 + t6 * t18, g(3) * t18 + t6 * t19, 0, 0;];
taug_reg = t9;
