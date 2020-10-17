% Calculate minimal parameter regressor of gravitation load for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% taug_reg [6x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:33:59
% EndTime: 2019-05-05 02:34:01
% DurationCPUTime: 0.43s
% Computational Cost: add. (332->94), mult. (573->159), div. (0->0), fcn. (691->14), ass. (0->56)
t35 = sin(qJ(3));
t37 = cos(qJ(3));
t53 = cos(pkin(6));
t32 = sin(pkin(6));
t36 = sin(qJ(2));
t57 = t32 * t36;
t69 = -t35 * t57 + t53 * t37;
t38 = cos(qJ(2));
t31 = sin(pkin(10));
t49 = t31 * t53;
t52 = cos(pkin(10));
t15 = -t36 * t49 + t52 * t38;
t56 = t32 * t37;
t68 = -t15 * t35 + t31 * t56;
t67 = g(3) * t32;
t43 = t53 * t52;
t12 = t31 * t36 - t38 * t43;
t13 = t31 * t38 + t36 * t43;
t23 = t37 * pkin(3) + pkin(2);
t34 = -qJ(4) - pkin(8);
t66 = -t12 * t23 - t13 * t34;
t14 = t52 * t36 + t38 * t49;
t65 = -t14 * t23 - t15 * t34;
t28 = pkin(12) + qJ(6);
t24 = sin(t28);
t29 = qJ(3) + pkin(11);
t27 = cos(t29);
t63 = t24 * t27;
t26 = cos(t28);
t62 = t26 * t27;
t30 = sin(pkin(12));
t61 = t27 * t30;
t33 = cos(pkin(12));
t60 = t27 * t33;
t59 = t27 * t38;
t58 = t31 * t32;
t55 = t32 * t38;
t54 = t34 * t36;
t48 = t32 * t52;
t46 = t68 * pkin(3);
t25 = sin(t29);
t45 = pkin(4) * t27 + qJ(5) * t25;
t44 = t69 * pkin(3);
t2 = t13 * t25 + t27 * t48;
t4 = t15 * t25 - t27 * t58;
t8 = t25 * t57 - t53 * t27;
t42 = g(1) * t4 + g(2) * t2 + g(3) * t8;
t41 = -t13 * t35 - t37 * t48;
t1 = -g(1) * t14 - g(2) * t12 + g(3) * t55;
t40 = g(1) * t15 + g(2) * t13 + g(3) * t57;
t39 = t41 * pkin(3);
t16 = t23 * t55;
t9 = t53 * t25 + t27 * t57;
t5 = t15 * t27 + t25 * t58;
t3 = t13 * t27 - t25 * t48;
t6 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t1, t40, 0, 0, 0, 0, 0, -t1 * t37, t1 * t35, -t40, -g(1) * t65 - g(2) * t66 - g(3) * (-t32 * t54 + t16) -g(1) * (-t14 * t60 + t15 * t30) - g(2) * (-t12 * t60 + t13 * t30) - (t30 * t36 + t33 * t59) * t67, -g(1) * (t14 * t61 + t15 * t33) - g(2) * (t12 * t61 + t13 * t33) - (-t30 * t59 + t33 * t36) * t67, -t1 * t25, -g(1) * (-t45 * t14 + t65) - g(2) * (-t45 * t12 + t66) - g(3) * t16 - (t45 * t38 - t54) * t67, 0, 0, 0, 0, 0, -g(1) * (-t14 * t62 + t15 * t24) - g(2) * (-t12 * t62 + t13 * t24) - (t24 * t36 + t26 * t59) * t67, -g(1) * (t14 * t63 + t15 * t26) - g(2) * (t12 * t63 + t13 * t26) - (-t24 * t59 + t26 * t36) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t68 - g(2) * t41 - g(3) * t69, -g(1) * (-t15 * t37 - t35 * t58) - g(2) * (-t13 * t37 + t35 * t48) - g(3) * (-t53 * t35 - t36 * t56) 0, -g(1) * t46 - g(2) * t39 - g(3) * t44, t42 * t33, -t42 * t30, -g(1) * t5 - g(2) * t3 - g(3) * t9, -g(1) * (-t4 * pkin(4) + t5 * qJ(5) + t46) - g(2) * (-t2 * pkin(4) + t3 * qJ(5) + t39) - g(3) * (-t8 * pkin(4) + t9 * qJ(5) + t44) 0, 0, 0, 0, 0, t42 * t26, -t42 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t14 * t26 - t5 * t24) - g(2) * (t12 * t26 - t3 * t24) - g(3) * (-t9 * t24 - t26 * t55) -g(1) * (-t14 * t24 - t5 * t26) - g(2) * (-t12 * t24 - t3 * t26) - g(3) * (t24 * t55 - t9 * t26);];
taug_reg  = t6;
