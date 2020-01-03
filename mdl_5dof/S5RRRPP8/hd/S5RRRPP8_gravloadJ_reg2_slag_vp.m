% Calculate inertial parameters regressor of gravitation load for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t35 = sin(qJ(1));
t38 = cos(qJ(1));
t18 = g(1) * t38 + g(2) * t35;
t34 = sin(qJ(2));
t67 = t18 * t34;
t66 = g(1) * t35;
t28 = t34 * pkin(7);
t37 = cos(qJ(2));
t30 = t37 * pkin(2);
t33 = sin(qJ(3));
t63 = t33 * t34;
t36 = cos(qJ(3));
t62 = t34 * t36;
t61 = t34 * t38;
t60 = t35 * t37;
t59 = t36 * t37;
t58 = t37 * t38;
t57 = t38 * t33;
t56 = -pkin(3) - qJ(5);
t55 = t30 + t28;
t54 = t38 * pkin(1) + t35 * pkin(6);
t53 = qJ(4) * t33;
t52 = -pkin(1) - t30;
t51 = -pkin(2) - t53;
t13 = t33 * t60 + t36 * t38;
t14 = t35 * t59 - t57;
t50 = -t13 * pkin(3) + qJ(4) * t14;
t15 = -t35 * t36 + t37 * t57;
t16 = t35 * t33 + t36 * t58;
t49 = -t15 * pkin(3) + qJ(4) * t16;
t48 = pkin(3) * t59 + t37 * t53 + t55;
t47 = pkin(2) * t58 + pkin(7) * t61 + t54;
t31 = t38 * pkin(6);
t46 = -t14 * pkin(3) - t13 * qJ(4) + t31;
t4 = g(1) * t13 - g(2) * t15;
t5 = g(1) * t14 - g(2) * t16;
t45 = -g(2) * t38 + t66;
t43 = t16 * pkin(3) + t15 * qJ(4) + t47;
t41 = (t52 - t28) * t66;
t2 = g(1) * t15 + g(2) * t13 + g(3) * t63;
t40 = g(1) * t16 + g(2) * t14 + g(3) * t62;
t39 = -g(3) * t37 + t67;
t24 = pkin(7) * t58;
t21 = pkin(7) * t60;
t19 = qJ(4) * t62;
t17 = t45 * t34;
t8 = g(3) * t34 + t18 * t37;
t7 = t39 * t36;
t6 = t39 * t33;
t1 = [0, 0, 0, 0, 0, 0, t45, t18, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t37, -t17, -t18, -g(1) * (-t35 * pkin(1) + t31) - g(2) * t54, 0, 0, 0, 0, 0, 0, t5, -t4, t17, -g(1) * t31 - g(2) * t47 - t41, 0, 0, 0, 0, 0, 0, t17, -t5, t4, -g(1) * t46 - g(2) * t43 - t41, 0, 0, 0, 0, 0, 0, t17, t4, t5, -g(1) * (-t14 * qJ(5) + t46) - g(2) * (pkin(4) * t61 + t16 * qJ(5) + t43) - ((-pkin(4) - pkin(7)) * t34 + t52) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t8, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, -t8, -g(1) * (-pkin(2) * t61 + t24) - g(2) * (-t35 * t34 * pkin(2) + t21) - g(3) * t55, 0, 0, 0, 0, 0, 0, -t8, -t7, t6, -g(1) * t24 - g(2) * t21 - g(3) * t48 + (pkin(3) * t36 - t51) * t67, 0, 0, 0, 0, 0, 0, -t8, t6, t7, -g(1) * (pkin(4) * t58 + t24) - g(2) * (pkin(4) * t60 + t21) - g(3) * (qJ(5) * t59 + t48) + (-g(3) * pkin(4) + t18 * (-t56 * t36 - t51)) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t40, -g(1) * t49 - g(2) * t50 - g(3) * (-pkin(3) * t63 + t19), 0, 0, 0, 0, 0, 0, 0, -t40, t2, -g(1) * (-qJ(5) * t15 + t49) - g(2) * (-qJ(5) * t13 + t50) - g(3) * (t56 * t63 + t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40;];
taug_reg = t1;
