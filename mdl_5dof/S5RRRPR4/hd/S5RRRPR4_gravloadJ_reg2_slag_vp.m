% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = qJ(1) + qJ(2);
t25 = sin(t31);
t26 = cos(t31);
t12 = g(1) * t26 + g(2) * t25;
t33 = sin(qJ(3));
t63 = t12 * t33;
t27 = t33 * qJ(4);
t36 = cos(qJ(3));
t50 = t36 * pkin(3) + t27;
t32 = sin(qJ(5));
t35 = cos(qJ(5));
t43 = t36 * t32 - t33 * t35;
t61 = g(1) * t25;
t34 = sin(qJ(1));
t58 = t34 * pkin(1);
t57 = t36 * pkin(4);
t56 = t25 * t33;
t55 = t26 * t33;
t54 = t26 * t36;
t51 = t26 * pkin(2) + t25 * pkin(7);
t49 = qJ(4) * t36;
t23 = t26 * pkin(7);
t48 = -t25 * pkin(2) + t23;
t47 = -t26 * pkin(8) + t23;
t46 = pkin(3) * t54 + t26 * t27 + t51;
t37 = cos(qJ(1));
t29 = t37 * pkin(1);
t45 = t29 + t46;
t11 = -g(2) * t26 + t61;
t44 = g(1) * t34 - g(2) * t37;
t13 = t33 * t32 + t36 * t35;
t42 = -pkin(2) - t50;
t5 = t43 * t25;
t7 = t43 * t26;
t41 = g(1) * t7 + g(2) * t5 + g(3) * t13;
t6 = t13 * t25;
t8 = t13 * t26;
t40 = g(1) * t8 + g(2) * t6 - g(3) * t43;
t39 = t42 * t61;
t38 = (-g(1) * (t42 - t57) + g(2) * pkin(8)) * t25;
t19 = pkin(4) * t54;
t18 = t26 * t49;
t16 = t25 * t49;
t10 = t11 * t36;
t9 = g(1) * t56 - g(2) * t55;
t4 = g(3) * t33 + t12 * t36;
t3 = -g(3) * t36 + t63;
t2 = g(1) * t6 - g(2) * t8;
t1 = -g(1) * t5 + g(2) * t7;
t14 = [0, 0, 0, 0, 0, 0, t44, g(1) * t37 + g(2) * t34, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, t44 * pkin(1), 0, 0, 0, 0, 0, 0, t10, -t9, -t12, -g(1) * (t48 - t58) - g(2) * (t29 + t51), 0, 0, 0, 0, 0, 0, t10, -t12, t9, -g(1) * (t23 - t58) - g(2) * t45 - t39, 0, 0, 0, 0, 0, 0, t2, t1, t12, -g(1) * (t47 - t58) - g(2) * (t19 + t45) + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, -t12, -g(1) * t48 - g(2) * t51, 0, 0, 0, 0, 0, 0, t10, -t12, t9, -g(1) * t23 - g(2) * t46 - t39, 0, 0, 0, 0, 0, 0, t2, t1, t12, -g(1) * t47 - g(2) * (t19 + t46) + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, -t4, -g(1) * (-pkin(3) * t55 + t18) - g(2) * (-pkin(3) * t56 + t16) - g(3) * t50, 0, 0, 0, 0, 0, 0, -t41, -t40, 0, -g(1) * t18 - g(2) * t16 - g(3) * (t50 + t57) + (pkin(3) + pkin(4)) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t40, 0, 0;];
taug_reg = t14;
