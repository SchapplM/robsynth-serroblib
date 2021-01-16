% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [6x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:38:31
% EndTime: 2021-01-16 01:38:33
% DurationCPUTime: 0.37s
% Computational Cost: add. (378->85), mult. (658->135), div. (0->0), fcn. (795->11), ass. (0->60)
t33 = cos(pkin(6));
t37 = sin(qJ(2));
t52 = cos(pkin(10));
t47 = t52 * t37;
t30 = sin(pkin(10));
t39 = cos(qJ(2));
t57 = t30 * t39;
t18 = t33 * t47 + t57;
t29 = pkin(11) + qJ(4);
t27 = sin(t29);
t28 = cos(t29);
t31 = sin(pkin(6));
t50 = t31 * t52;
t10 = t18 * t27 + t28 * t50;
t46 = t52 * t39;
t58 = t30 * t37;
t16 = t33 * t58 - t46;
t61 = t30 * t31;
t12 = t16 * t27 + t28 * t61;
t56 = t31 * t37;
t14 = t27 * t56 - t33 * t28;
t41 = g(1) * t12 - g(2) * t10 - g(3) * t14;
t66 = g(3) * t31;
t36 = sin(qJ(5));
t65 = t16 * t36;
t64 = t18 * t36;
t63 = t28 * t36;
t38 = cos(qJ(5));
t62 = t28 * t38;
t60 = t30 * t33;
t35 = qJ(3) + pkin(8);
t59 = t30 * t35;
t55 = t36 * t39;
t54 = t38 * t39;
t53 = t30 * qJ(3);
t51 = t52 * pkin(2);
t9 = t16 * t28 - t27 * t61;
t32 = cos(pkin(11));
t25 = t32 * pkin(3) + pkin(2);
t49 = t52 * t25;
t48 = t52 * t35;
t45 = t52 * qJ(3);
t26 = t38 * pkin(5) + pkin(4);
t34 = -qJ(6) - pkin(9);
t43 = -t28 * t26 + t27 * t34;
t11 = t18 * t28 - t27 * t50;
t15 = t33 * t27 + t28 * t56;
t42 = g(1) * t9 - g(2) * t11 - g(3) * t15;
t40 = -g(1) * t16 + g(2) * t18 + g(3) * t56;
t17 = -t33 * t46 + t58;
t19 = t33 * t57 + t47;
t8 = -g(1) * t19 - g(2) * t17 + t39 * t66;
t1 = -g(1) * (t19 * t38 + t36 * t9) - g(2) * (-t11 * t36 + t17 * t38) - g(3) * (-t15 * t36 - t31 * t54);
t7 = t8 * t27;
t6 = t41 * t38;
t5 = t41 * t36;
t4 = -g(1) * (-t19 * t62 - t65) - g(2) * (-t17 * t62 + t64) - (t28 * t54 + t36 * t37) * t66;
t3 = -g(1) * (-t16 * t38 + t19 * t63) - g(2) * (t17 * t63 + t18 * t38) - (-t28 * t55 + t37 * t38) * t66;
t2 = -g(1) * (-t19 * t36 + t38 * t9) - g(2) * (-t11 * t38 - t17 * t36) - g(3) * (-t15 * t38 + t31 * t55);
t13 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t8, t40, -t8 * t32, -t40, -g(1) * (-(t33 * t53 + t51) * t37 + (-pkin(2) * t60 + t45) * t39) - g(2) * (-(t30 * pkin(2) - t33 * t45) * t37 + (t33 * t51 + t53) * t39) - (pkin(2) * t39 + qJ(3) * t37) * t66, 0, 0, 0, 0, 0, -t8 * t28, t7, 0, 0, 0, 0, 0, t4, t3, t4, t3, -t7, -g(1) * (-pkin(5) * t65 - (t33 * t59 + t49) * t37 + (-t25 * t60 + t48) * t39 + t43 * t19) - g(2) * (pkin(5) * t64 - (t30 * t25 - t33 * t48) * t37 + (t33 * t49 + t59) * t39 + t43 * t17) - ((pkin(5) * t36 + t35) * t37 + (t25 - t43) * t39) * t66; 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t42, 0, 0, 0, 0, 0, -t6, t5, -t6, t5, t42, -g(1) * (t12 * t26 + t9 * t34) - g(2) * (-t10 * t26 - t11 * t34) - g(3) * (-t14 * t26 - t15 * t34); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41;];
taug_reg = t13;
