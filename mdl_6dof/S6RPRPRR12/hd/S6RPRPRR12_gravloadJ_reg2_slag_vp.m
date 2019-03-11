% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t36 = cos(qJ(5));
t33 = sin(qJ(5));
t38 = cos(qJ(1));
t54 = t38 * t33;
t35 = sin(qJ(1));
t37 = cos(qJ(3));
t57 = t35 * t37;
t11 = -t36 * t57 - t54;
t53 = t38 * t36;
t42 = -t35 * t33 + t37 * t53;
t34 = sin(qJ(3));
t63 = g(3) * t34;
t71 = -g(1) * t11 - g(2) * t42 - t36 * t63;
t26 = t37 * qJ(4);
t66 = pkin(5) * t33;
t70 = t37 * t66 + t26;
t64 = g(2) * t38;
t15 = g(1) * t35 - t64;
t7 = g(3) * t37 + t15 * t34;
t68 = -pkin(1) - pkin(7);
t67 = -pkin(3) - pkin(8);
t61 = t34 * pkin(3);
t39 = -pkin(9) - pkin(8);
t60 = -pkin(3) + t39;
t59 = t33 * t34;
t58 = t34 * t38;
t32 = qJ(5) + qJ(6);
t23 = sin(t32);
t56 = t38 * t23;
t24 = cos(t32);
t55 = t38 * t24;
t49 = qJ(4) * t34;
t52 = pkin(3) * t57 + t35 * t49;
t27 = t38 * qJ(2);
t51 = pkin(3) * t58 + t27;
t50 = t38 * pkin(1) + t35 * qJ(2);
t46 = t38 * pkin(7) + t50;
t45 = t68 * t35;
t44 = t35 * t61 + t46;
t16 = g(1) * t38 + g(2) * t35;
t43 = -t38 * t26 + t51;
t41 = t34 * t39 + t70;
t22 = t36 * pkin(5) + pkin(4);
t14 = t16 * t37;
t13 = t16 * t34;
t12 = -t33 * t57 + t53;
t10 = -t35 * t36 - t37 * t54;
t8 = g(1) * t57 - t37 * t64 - t63;
t6 = -t23 * t57 + t55;
t5 = -t24 * t57 - t56;
t4 = -t35 * t24 - t37 * t56;
t3 = t35 * t23 - t37 * t55;
t2 = g(1) * t6 - g(2) * t4 + t23 * t63;
t1 = -g(1) * t5 + g(2) * t3 - t24 * t63;
t9 = [0, 0, 0, 0, 0, 0, t15, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, -g(1) * (-t35 * pkin(1) + t27) - g(2) * t50, 0, 0, 0, 0, 0, 0, -t13, -t14, t15, -g(1) * (t27 + t45) - g(2) * t46, 0, 0, 0, 0, 0, 0, t15, t13, t14, -g(1) * (t45 + t43) - g(2) * (-t35 * t26 + t44) 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, g(1) * t42 - g(2) * t11, -t13, -g(1) * (pkin(8) * t58 + t43) - g(2) * (t38 * pkin(4) + t44) + (-g(1) * (-pkin(4) + t68) - g(2) * (t34 * pkin(8) - t26)) * t35, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, -t13, -g(1) * t51 - g(2) * t44 + (g(1) * t41 - g(2) * t22) * t38 + (-g(1) * (-t22 + t68) + g(2) * t41) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -g(1) * t52 - g(3) * (t26 - t61) - (-pkin(3) * t37 - t49) * t64, 0, 0, 0, 0, 0, 0, -t7 * t33, -t7 * t36, -t8, -g(1) * (pkin(8) * t57 + t52) - g(3) * (t67 * t34 + t26) - (t67 * t37 - t49) * t64, 0, 0, 0, 0, 0, 0, -t7 * t23, -t7 * t24, -t8, -g(1) * ((pkin(5) * t59 - t37 * t39) * t35 + t52) - g(3) * (t60 * t34 + t70) - (t60 * t37 + (-qJ(4) - t66) * t34) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, g(1) * t12 - g(2) * t10 + g(3) * t59, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t71 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t9;
