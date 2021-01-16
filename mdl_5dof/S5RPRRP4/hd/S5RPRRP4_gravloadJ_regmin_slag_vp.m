% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:56:35
% EndTime: 2021-01-15 12:56:37
% DurationCPUTime: 0.20s
% Computational Cost: add. (172->47), mult. (236->73), div. (0->0), fcn. (258->8), ass. (0->39)
t26 = qJ(3) + qJ(4);
t20 = cos(t26);
t31 = cos(qJ(3));
t15 = t31 * pkin(3) + pkin(4) * t20;
t27 = sin(pkin(8));
t28 = cos(pkin(8));
t49 = -(-qJ(5) - pkin(7) - pkin(6)) * t27 + (pkin(2) + t15) * t28;
t19 = sin(t26);
t46 = g(1) * t27;
t32 = cos(qJ(1));
t37 = t32 * t20;
t30 = sin(qJ(1));
t42 = t30 * t19;
t5 = -t28 * t42 - t37;
t38 = t32 * t19;
t41 = t30 * t20;
t7 = t28 * t38 - t41;
t1 = -g(2) * t5 - g(3) * t7 + t19 * t46;
t29 = sin(qJ(3));
t14 = t29 * pkin(3) + pkin(4) * t19;
t43 = t30 * t14;
t40 = t30 * t29;
t39 = t30 * t31;
t36 = t32 * t29;
t35 = t32 * t31;
t34 = t32 * pkin(1) + t30 * qJ(2);
t17 = g(2) * t32 + g(3) * t30;
t16 = g(2) * t30 - g(3) * t32;
t22 = t30 * pkin(1);
t12 = t28 * t35 + t40;
t11 = t28 * t36 - t39;
t10 = t28 * t39 - t36;
t9 = -t28 * t40 - t35;
t8 = t28 * t37 + t42;
t6 = t28 * t41 - t38;
t4 = -g(2) * t8 - g(3) * t6;
t3 = g(2) * t7 - g(3) * t5;
t2 = g(2) * t6 - g(3) * t8 + t20 * t46;
t13 = [0, -t17, t16, -t17 * t28, -t16, -g(2) * t34 - g(3) * (-t32 * qJ(2) + t22), 0, 0, 0, 0, 0, -g(2) * t12 - g(3) * t10, g(2) * t11 - g(3) * t9, 0, 0, 0, 0, 0, t4, t3, t4, t3, -t17 * t27, -g(2) * (t34 + t43) - g(3) * (t49 * t30 + t22) + (-g(2) * t49 - g(3) * (-qJ(2) - t14)) * t32; 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t9 - g(3) * t11 + t29 * t46, g(2) * t10 - g(3) * t12 + t31 * t46, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t14 * t46 - g(2) * (-t32 * t15 - t28 * t43) - g(3) * (t32 * t28 * t14 - t30 * t15); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t28 - t16 * t27;];
taug_reg = t13;
