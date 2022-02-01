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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:32:52
% EndTime: 2022-01-23 09:32:52
% DurationCPUTime: 0.19s
% Computational Cost: add. (172->46), mult. (236->73), div. (0->0), fcn. (258->8), ass. (0->39)
t27 = qJ(3) + qJ(4);
t20 = cos(t27);
t32 = cos(qJ(3));
t15 = t32 * pkin(3) + pkin(4) * t20;
t28 = sin(pkin(8));
t29 = cos(pkin(8));
t50 = -(-qJ(5) - pkin(7) - pkin(6)) * t28 + (pkin(2) + t15) * t29;
t19 = sin(t27);
t47 = g(3) * t28;
t33 = cos(qJ(1));
t38 = t33 * t20;
t31 = sin(qJ(1));
t44 = t31 * t19;
t5 = t29 * t44 + t38;
t39 = t33 * t19;
t43 = t31 * t20;
t7 = -t29 * t39 + t43;
t1 = -g(1) * t7 + g(2) * t5 + t19 * t47;
t30 = sin(qJ(3));
t42 = t31 * t30;
t41 = t31 * t32;
t14 = t30 * pkin(3) + pkin(4) * t19;
t40 = t33 * t14;
t37 = t33 * t30;
t36 = t33 * t32;
t35 = t33 * pkin(1) + t31 * qJ(2);
t17 = g(1) * t33 + g(2) * t31;
t16 = g(1) * t31 - g(2) * t33;
t22 = t33 * qJ(2);
t12 = t29 * t36 + t42;
t11 = -t29 * t37 + t41;
t10 = -t29 * t41 + t37;
t9 = t29 * t42 + t36;
t8 = t29 * t38 + t44;
t6 = -t29 * t43 + t39;
t4 = -g(1) * t6 - g(2) * t8;
t3 = -g(1) * t5 - g(2) * t7;
t2 = g(1) * t8 - g(2) * t6 + t20 * t47;
t13 = [0, t16, t17, t16 * t29, -t17, -g(1) * (-t31 * pkin(1) + t22) - g(2) * t35, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, 0, 0, 0, 0, 0, t4, t3, t4, t3, t16 * t28, -g(1) * (t22 + t40) - g(2) * (t50 * t33 + t35) + (-g(1) * (-pkin(1) - t50) - g(2) * t14) * t31; 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t9 + t30 * t47, g(1) * t12 - g(2) * t10 + t32 * t47, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, -g(1) * (t31 * t15 - t29 * t40) - g(2) * (-t31 * t29 * t14 - t33 * t15) + t14 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t29 - t17 * t28;];
taug_reg = t13;
