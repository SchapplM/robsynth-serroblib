% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t19 = g(1) * t39 + g(2) * t36;
t31 = qJ(2) + pkin(10);
t26 = cos(t31);
t37 = cos(qJ(4));
t50 = t39 * t37;
t34 = sin(qJ(4));
t55 = t36 * t34;
t13 = t26 * t55 + t50;
t51 = t39 * t34;
t54 = t36 * t37;
t15 = -t26 * t51 + t54;
t25 = sin(t31);
t62 = g(3) * t25;
t67 = -g(1) * t15 + g(2) * t13 + t34 * t62;
t44 = -g(3) * t26 + t19 * t25;
t32 = qJ(4) + qJ(5);
t27 = sin(t32);
t60 = t25 * t27;
t28 = cos(t32);
t59 = t25 * t28;
t40 = -pkin(9) - pkin(8);
t58 = t25 * t40;
t57 = t36 * t27;
t56 = t36 * t28;
t53 = t39 * t27;
t52 = t39 * t28;
t33 = -qJ(3) - pkin(7);
t48 = pkin(4) * t34 - t33;
t11 = t26 * t53 - t56;
t9 = t26 * t57 + t52;
t47 = g(1) * t9 - g(2) * t11;
t18 = g(1) * t36 - g(2) * t39;
t23 = t37 * pkin(4) + pkin(3);
t46 = t26 * t23 - t58;
t45 = pkin(5) * t28 + qJ(6) * t27 + t23;
t1 = g(1) * t11 + g(2) * t9 + g(3) * t60;
t10 = t26 * t56 - t53;
t12 = t26 * t52 + t57;
t3 = g(1) * t12 + g(2) * t10 + g(3) * t59;
t35 = sin(qJ(2));
t38 = cos(qJ(2));
t43 = -g(3) * t38 + t19 * t35;
t41 = -g(1) * (-t11 * pkin(5) + t12 * qJ(6)) - g(2) * (-t9 * pkin(5) + t10 * qJ(6)) - g(3) * (-pkin(5) * t60 + qJ(6) * t59);
t29 = t38 * pkin(2);
t24 = t29 + pkin(1);
t20 = t39 * t24;
t16 = t26 * t50 + t55;
t14 = -t26 * t54 + t51;
t6 = t44 * t28;
t5 = t44 * t27;
t4 = g(1) * t10 - g(2) * t12;
t2 = [0, t18, t19, 0, 0, 0, 0, 0, t18 * t38, -t18 * t35, -t19, -g(1) * (-t36 * t24 - t39 * t33) - g(2) * (-t36 * t33 + t20) 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, 0, 0, 0, 0, 0, t4, -t47, t4, t18 * t25, t47, -g(1) * (-t10 * pkin(5) - t9 * qJ(6)) - g(2) * (t12 * pkin(5) + t11 * qJ(6) + t20) + (-g(1) * t48 - g(2) * t46) * t39 + (-g(1) * (-t24 - t46) - g(2) * t48) * t36; 0, 0, 0, 0, 0, 0, 0, 0, t43, g(3) * t35 + t19 * t38, 0, t43 * pkin(2), 0, 0, 0, 0, 0, t44 * t37, -t44 * t34, 0, 0, 0, 0, 0, t6, -t5, t6, -t19 * t26 - t62, t5, -g(3) * (t45 * t26 + t29 - t58) + t19 * (pkin(2) * t35 + t45 * t25 + t26 * t40); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, g(1) * t16 - g(2) * t14 + t37 * t62, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t67 * pkin(4) + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
