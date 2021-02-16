% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:59:01
% EndTime: 2021-01-15 19:59:02
% DurationCPUTime: 0.17s
% Computational Cost: add. (121->47), mult. (181->63), div. (0->0), fcn. (181->10), ass. (0->36)
t29 = sin(qJ(1));
t32 = cos(qJ(1));
t15 = g(1) * t32 + g(2) * t29;
t23 = qJ(2) + pkin(8);
t19 = sin(t23);
t20 = cos(t23);
t4 = g(3) * t19 + t15 * t20;
t16 = g(3) * t20;
t26 = -qJ(3) - pkin(6);
t39 = t29 * t26;
t27 = sin(qJ(5));
t38 = t29 * t27;
t30 = cos(qJ(5));
t37 = t29 * t30;
t36 = t32 * t27;
t35 = t32 * t30;
t14 = g(1) * t29 - g(2) * t32;
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t33 = -g(3) * t31 + t15 * t28;
t25 = cos(pkin(8));
t24 = sin(pkin(8));
t21 = t31 * pkin(2);
t18 = t21 + pkin(1);
t17 = t26 * t32;
t13 = -t24 * pkin(3) + qJ(4) * t25;
t12 = pkin(3) * t25 + qJ(4) * t24 + pkin(2);
t10 = -t19 * t38 + t35;
t9 = t19 * t37 + t36;
t8 = t19 * t36 + t37;
t7 = t19 * t35 - t38;
t6 = t14 * t20;
t5 = t14 * t19;
t3 = t15 * t19 - t16;
t1 = t12 * t31 + t13 * t28 + pkin(1);
t2 = [0, t14, t15, 0, 0, 0, 0, 0, t14 * t31, -t14 * t28, t6, -t5, -t15, -g(1) * (-t29 * t18 - t17) - g(2) * (t32 * t18 - t39), -t15, -t6, t5, -g(1) * (-t1 * t29 - t17) - g(2) * (t1 * t32 - t39), 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, t33, g(3) * t28 + t15 * t31, t3, t4, 0, t33 * pkin(2), 0, -t3, -t4, -g(3) * (t20 * pkin(3) + t19 * qJ(4) + t21) - t15 * (-t12 * t28 + t13 * t31), 0, 0, 0, 0, 0, -t4 * t27, -t4 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 - t15 * (t24 * t31 + t25 * t28), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t30 * t16, g(1) * t8 - g(2) * t10 - t27 * t16;];
taug_reg = t2;
