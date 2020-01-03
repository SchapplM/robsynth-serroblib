% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t15 = g(1) * t37 + g(2) * t34;
t33 = sin(qJ(2));
t68 = t15 * t33;
t31 = sin(qJ(5));
t32 = sin(qJ(4));
t36 = cos(qJ(2));
t59 = cos(qJ(4));
t13 = t33 * t32 + t36 * t59;
t52 = t33 * t59;
t14 = -t36 * t32 + t52;
t7 = t14 * t34;
t57 = t36 * t37;
t9 = t32 * t57 - t37 * t52;
t40 = g(1) * t9 - g(2) * t7 + g(3) * t13;
t67 = t40 * t31;
t35 = cos(qJ(5));
t66 = t40 * t35;
t24 = t33 * qJ(3);
t56 = t36 * pkin(2) + t24;
t64 = pkin(2) * t33;
t63 = g(1) * t34;
t60 = g(3) * t14;
t26 = t36 * pkin(3);
t55 = t37 * pkin(1) + t34 * pkin(6);
t54 = qJ(3) * t36;
t53 = t26 + t56;
t28 = t37 * pkin(6);
t51 = -t37 * pkin(7) + t28;
t50 = pkin(2) * t57 + t37 * t24 + t55;
t8 = t13 * t34;
t49 = t7 * pkin(4) + t8 * pkin(8);
t48 = -g(1) * t7 - g(2) * t9;
t10 = t13 * t37;
t47 = -t9 * pkin(4) + t10 * pkin(8);
t46 = pkin(3) * t57 + t50;
t45 = -t13 * pkin(4) + t14 * pkin(8);
t44 = -g(2) * t37 + t63;
t43 = t8 * t31 - t37 * t35;
t42 = t37 * t31 + t8 * t35;
t41 = -pkin(1) - t56;
t2 = g(1) * t10 + g(2) * t8 + t60;
t39 = (pkin(2) + pkin(3)) * t68;
t38 = (-g(1) * (t41 - t26) + g(2) * pkin(7)) * t34;
t20 = t37 * t54;
t18 = t34 * t54;
t12 = t44 * t36;
t11 = t44 * t33;
t6 = g(3) * t33 + t15 * t36;
t5 = -g(3) * t36 + t68;
t4 = t10 * t35 - t34 * t31;
t3 = -t10 * t31 - t34 * t35;
t1 = [0, 0, 0, 0, 0, 0, t44, t15, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, -t15, -g(1) * (-t34 * pkin(1) + t28) - g(2) * t55, 0, 0, 0, 0, 0, 0, t12, -t15, t11, -g(1) * t28 - g(2) * t50 - t41 * t63, 0, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t10, -t48, t15, -g(1) * t51 - g(2) * t46 + t38, 0, 0, 0, 0, 0, 0, g(1) * t42 - g(2) * t4, -g(1) * t43 - g(2) * t3, t48, -g(1) * (-t8 * pkin(4) + t7 * pkin(8) + t51) - g(2) * (t10 * pkin(4) + t9 * pkin(8) + t46) + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, -t6, -g(1) * (-t37 * t64 + t20) - g(2) * (-t34 * t64 + t18) - g(3) * t56, 0, 0, 0, 0, 0, 0, -t40, -t2, 0, -g(1) * t20 - g(2) * t18 - g(3) * t53 + t39, 0, 0, 0, 0, 0, 0, -t66, t67, t2, -g(1) * (t20 - t47) - g(2) * (t18 - t49) - g(3) * (-t45 + t53) + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t2, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t67, -t2, -g(1) * t47 - g(2) * t49 - g(3) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t43 + t31 * t60, g(1) * t4 + g(2) * t42 + t35 * t60, 0, 0;];
taug_reg = t1;
