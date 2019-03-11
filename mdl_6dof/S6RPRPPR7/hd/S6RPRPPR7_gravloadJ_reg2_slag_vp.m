% Calculate inertial parameters regressor of gravitation load for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t28 = qJ(3) + pkin(9);
t22 = sin(t28);
t34 = cos(qJ(3));
t61 = pkin(3) * t34 + qJ(5) * t22;
t29 = -qJ(4) - pkin(7);
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t31 = sin(qJ(3));
t53 = t31 * pkin(3);
t60 = t29 * t32 + t35 * t53;
t57 = g(2) * t35;
t10 = g(1) * t32 - t57;
t23 = cos(t28);
t1 = g(3) * t23 + t10 * t22;
t59 = -pkin(4) - pkin(8);
t56 = g(3) * t22;
t54 = t22 * pkin(4);
t52 = t23 * t32;
t30 = sin(qJ(6));
t51 = t32 * t30;
t33 = cos(qJ(6));
t50 = t32 * t33;
t49 = t35 * t30;
t48 = t35 * t33;
t47 = pkin(1) * t35 + qJ(2) * t32;
t19 = t23 * qJ(5);
t45 = pkin(4) * t52 + t32 * t61;
t44 = t32 * t53 + t47;
t25 = t35 * qJ(2);
t43 = -t32 * pkin(1) + t25;
t42 = t19 - t53;
t11 = g(1) * t35 + g(2) * t32;
t40 = pkin(8) * t22 - t19;
t39 = t43 + t60;
t38 = -t35 * t29 + t44;
t36 = g(3) * t31 - t10 * t34;
t14 = t35 * t54;
t12 = t32 * t54;
t8 = -t23 * t51 + t48;
t7 = -t23 * t50 - t49;
t6 = -t23 * t49 - t50;
t5 = -t23 * t48 + t51;
t4 = t11 * t23;
t3 = t11 * t22;
t2 = g(1) * t52 - t23 * t57 - t56;
t9 = [0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, -g(1) * t43 - g(2) * t47, 0, 0, 0, 0, 0, 0, -t11 * t31, -t11 * t34, t10, -g(1) * (t25 + (-pkin(1) - pkin(7)) * t32) - g(2) * (pkin(7) * t35 + t47) 0, 0, 0, 0, 0, 0, -t3, -t4, t10, -g(1) * t39 - g(2) * t38, 0, 0, 0, 0, 0, 0, t10, t3, t4, -g(1) * (-t19 * t35 + t14 + t39) - g(2) * (-t19 * t32 + t12 + t38) 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7, -t3, -g(1) * (t14 + t25 + t60) - g(2) * (t12 + t44) + (-g(1) * t40 - g(2) * (pkin(5) - t29)) * t35 + (-g(1) * (-pkin(1) - pkin(5)) - g(2) * t40) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, g(3) * t34 + t10 * t31, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, t36 * pkin(3), 0, 0, 0, 0, 0, 0, 0, t2, -t1, -g(1) * t45 - g(3) * (t42 - t54) - (-pkin(4) * t23 - t61) * t57, 0, 0, 0, 0, 0, 0, -t1 * t30, -t1 * t33, -t2, -g(1) * (pkin(8) * t52 + t45) - g(3) * (t22 * t59 + t42) - (t23 * t59 - t61) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t5 - t33 * t56, g(1) * t8 - g(2) * t6 + t30 * t56, 0, 0;];
taug_reg  = t9;
