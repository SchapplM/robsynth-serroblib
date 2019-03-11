% Calculate minimal parameter regressor of gravitation load for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t29 = qJ(2) + pkin(9);
t25 = sin(t29);
t26 = cos(t29);
t66 = t26 * pkin(3) + t25 * qJ(4);
t35 = sin(qJ(1));
t38 = cos(qJ(1));
t17 = g(1) * t38 + g(2) * t35;
t64 = t17 * t25;
t41 = -g(3) * t26 + t64;
t34 = sin(qJ(2));
t63 = pkin(2) * t34;
t60 = g(3) * t25;
t30 = sin(pkin(10));
t58 = t35 * t30;
t31 = cos(pkin(10));
t57 = t35 * t31;
t56 = t38 * t30;
t55 = t38 * t31;
t32 = -qJ(3) - pkin(7);
t54 = t38 * t32;
t53 = qJ(4) * t26;
t37 = cos(qJ(2));
t27 = t37 * pkin(2);
t24 = t27 + pkin(1);
t19 = t38 * t24;
t52 = t66 * t38 + t19;
t51 = t27 + t66;
t11 = t26 * t56 - t57;
t9 = t26 * t58 + t55;
t50 = g(1) * t9 - g(2) * t11;
t49 = -pkin(3) * t25 - t63;
t16 = g(1) * t35 - g(2) * t38;
t10 = t26 * t57 - t56;
t33 = sin(qJ(6));
t36 = cos(qJ(6));
t48 = t10 * t36 + t9 * t33;
t47 = t10 * t33 - t9 * t36;
t46 = pkin(4) * t31 + qJ(5) * t30;
t45 = t30 * t36 - t31 * t33;
t44 = t30 * t33 + t31 * t36;
t42 = g(3) * t45;
t40 = -g(3) * t37 + t17 * t34;
t39 = (-g(1) * (-t24 - t66) + g(2) * t32) * t35;
t15 = t38 * t53;
t13 = t35 * t53;
t12 = t26 * t55 + t58;
t8 = t16 * t25;
t7 = -t17 * t26 - t60;
t5 = t41 * t31;
t4 = t41 * t30;
t3 = g(1) * t10 - g(2) * t12;
t2 = t11 * t33 + t12 * t36;
t1 = t11 * t36 - t12 * t33;
t6 = [0, t16, t17, 0, 0, 0, 0, 0, t16 * t37, -t16 * t34, -t17, -g(1) * (-t35 * t24 - t54) - g(2) * (-t35 * t32 + t19) t3, -t50, t8, g(1) * t54 - g(2) * t52 + t39, t3, t8, t50, -g(1) * (-t10 * pkin(4) - t9 * qJ(5) - t54) - g(2) * (t12 * pkin(4) + t11 * qJ(5) + t52) + t39, 0, 0, 0, 0, 0, g(1) * t48 - g(2) * t2, -g(1) * t47 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, t40, g(3) * t34 + t17 * t37, 0, t40 * pkin(2), t5, -t4, t7, -g(1) * (t49 * t38 + t15) - g(2) * (t49 * t35 + t13) - g(3) * t51, t5, t7, t4, -g(1) * (-t38 * t63 + t15) - g(2) * (-t35 * t63 + t13) - g(3) * (t46 * t26 + t51) + (pkin(3) + t46) * t64, 0, 0, 0, 0, 0, t41 * t44, -t26 * t42 + t45 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, -t16, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, 0, 0, -t41, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t9 - t30 * t60, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t47 - t25 * t42, g(1) * t2 + g(2) * t48 + t44 * t60;];
taug_reg  = t6;
