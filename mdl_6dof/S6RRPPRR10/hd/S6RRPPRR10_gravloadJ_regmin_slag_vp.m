% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t34 = sin(qJ(2));
t36 = cos(qJ(2));
t44 = t36 * pkin(2) + t34 * qJ(3);
t39 = -pkin(1) - t44;
t35 = sin(qJ(1));
t37 = cos(qJ(1));
t15 = g(1) * t37 + g(2) * t35;
t64 = t15 * t34;
t12 = g(3) * t34 + t15 * t36;
t63 = pkin(2) * t34;
t62 = g(1) * t35;
t58 = g(3) * t36;
t31 = pkin(10) + qJ(5);
t24 = qJ(6) + t31;
t20 = sin(t24);
t57 = t35 * t20;
t21 = cos(t24);
t56 = t35 * t21;
t22 = sin(t31);
t55 = t35 * t22;
t23 = cos(t31);
t54 = t35 * t23;
t32 = sin(pkin(10));
t53 = t35 * t32;
t33 = cos(pkin(10));
t52 = t35 * t33;
t51 = t37 * t20;
t50 = t37 * t21;
t49 = t37 * t22;
t48 = t37 * t23;
t47 = t37 * t32;
t46 = t37 * t33;
t43 = qJ(3) * t36;
t42 = t36 * qJ(4);
t41 = t35 * pkin(7) - t39 * t37;
t40 = -g(2) * t37 + t62;
t28 = t37 * pkin(7);
t18 = t37 * t43;
t16 = t35 * t43;
t14 = t40 * t36;
t13 = t40 * t34;
t11 = -t58 + t64;
t10 = -t34 * t55 + t48;
t9 = t34 * t54 + t49;
t8 = t34 * t49 + t54;
t7 = t34 * t48 - t55;
t6 = -t34 * t57 + t50;
t5 = t34 * t56 + t51;
t4 = t34 * t51 + t56;
t3 = t34 * t50 - t57;
t2 = g(1) * t4 - g(2) * t6 - t20 * t58;
t1 = -g(1) * t3 - g(2) * t5 + t21 * t58;
t17 = [0, t40, t15, 0, 0, 0, 0, 0, t14, -t13, -t15, -t14, t13, -g(1) * t28 - g(2) * t41 - t39 * t62, -g(1) * (-t34 * t53 + t46) - g(2) * (t34 * t47 + t52) -g(1) * (-t34 * t52 - t47) - g(2) * (t34 * t46 - t53) t14, -g(1) * (t37 * pkin(3) + t28) - g(2) * (t37 * t42 + t41) + (-g(1) * (t39 - t42) - g(2) * pkin(3)) * t35, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, -t11, -t12, -g(1) * (-t37 * t63 + t18) - g(2) * (-t35 * t63 + t16) - g(3) * t44, -t12 * t32, -t12 * t33, t11, -g(1) * t18 - g(2) * t16 - g(3) * (t42 + t44) + (pkin(2) + qJ(4)) * t64, 0, 0, 0, 0, 0, -t12 * t22, -t12 * t23, 0, 0, 0, 0, 0, -t12 * t20, -t12 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t23 * t58, g(1) * t8 - g(2) * t10 - t22 * t58, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t17;
