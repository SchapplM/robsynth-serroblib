% Calculate inertial parameters regressor of gravitation load for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t28 = -qJ(4) - pkin(7);
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t30 = sin(qJ(3));
t53 = t30 * pkin(3);
t60 = t31 * t28 + t33 * t53;
t56 = g(2) * t33;
t8 = g(1) * t31 - t56;
t25 = qJ(3) + pkin(9);
t17 = sin(t25);
t19 = cos(t25);
t2 = -g(3) * t17 + t8 * t19;
t32 = cos(qJ(3));
t59 = pkin(3) * t32;
t26 = sin(pkin(10));
t58 = pkin(5) * t26;
t54 = g(3) * t19;
t24 = pkin(10) + qJ(6);
t16 = sin(t24);
t52 = t31 * t16;
t18 = cos(t24);
t51 = t31 * t18;
t50 = t31 * t26;
t27 = cos(pkin(10));
t49 = t31 * t27;
t48 = t33 * t16;
t47 = t33 * t18;
t46 = t33 * t26;
t45 = t33 * t27;
t44 = t33 * pkin(1) + t31 * qJ(2);
t43 = t31 * t53 + t44;
t21 = t33 * qJ(2);
t42 = -t31 * pkin(1) + t21;
t9 = g(1) * t33 + g(2) * t31;
t41 = pkin(4) * t19 + qJ(5) * t17;
t40 = t17 * pkin(4) - t19 * qJ(5);
t14 = t27 * pkin(5) + pkin(4);
t29 = -pkin(8) - qJ(5);
t39 = t14 * t19 - t17 * t29;
t38 = t17 * t14 + t19 * t29;
t37 = t42 + t60;
t36 = -t33 * t28 + t43;
t34 = g(3) * t30 - t8 * t32;
t12 = t31 * t59;
t7 = t9 * t19;
t6 = t17 * t47 - t52;
t5 = t17 * t48 + t51;
t4 = t17 * t51 + t48;
t3 = -t17 * t52 + t47;
t1 = t8 * t17 + t54;
t10 = [0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, -g(1) * t42 - g(2) * t44, 0, 0, 0, 0, 0, 0, -t9 * t30, -t9 * t32, t8, -g(1) * (t21 + (-pkin(1) - pkin(7)) * t31) - g(2) * (t33 * pkin(7) + t44) 0, 0, 0, 0, 0, 0, -t9 * t17, -t7, t8, -g(1) * t37 - g(2) * t36, 0, 0, 0, 0, 0, 0, -g(1) * (t17 * t45 - t50) - g(2) * (t17 * t49 + t46) -g(1) * (-t17 * t46 - t49) - g(2) * (-t17 * t50 + t45) t7, -g(1) * (t40 * t33 + t37) - g(2) * (t40 * t31 + t36) 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t7, -g(1) * (t21 + t60) - g(2) * t43 + (-g(1) * t38 - g(2) * (-t28 + t58)) * t33 + (-g(1) * (-pkin(1) - t58) - g(2) * t38) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, g(3) * t32 + t8 * t30, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, t34 * pkin(3), 0, 0, 0, 0, 0, 0, -t2 * t27, t2 * t26, -t1, -g(1) * (t41 * t31 + t12) - g(3) * (-t40 - t53) - (-t41 - t59) * t56, 0, 0, 0, 0, 0, 0, -t2 * t18, t2 * t16, -t1, -g(1) * (t39 * t31 + t12) - g(3) * (-t38 - t53) - (-t39 - t59) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t16 * t54, g(1) * t4 - g(2) * t6 + t18 * t54, 0, 0;];
taug_reg  = t10;
