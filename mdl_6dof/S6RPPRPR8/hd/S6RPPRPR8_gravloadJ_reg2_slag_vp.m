% Calculate inertial parameters regressor of gravitation load for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t30 = -pkin(7) - qJ(3);
t32 = sin(qJ(1));
t34 = cos(qJ(1));
t28 = sin(pkin(9));
t53 = pkin(3) * t28;
t55 = t32 * t30 + t34 * t53;
t52 = g(2) * t34;
t10 = g(1) * t32 - t52;
t27 = pkin(9) + qJ(4);
t21 = sin(t27);
t22 = cos(t27);
t1 = g(3) * t22 + t10 * t21;
t54 = -pkin(4) - pkin(8);
t51 = g(3) * t21;
t49 = t21 * pkin(4);
t41 = qJ(5) * t21;
t47 = t22 * t32;
t48 = pkin(4) * t47 + t32 * t41;
t31 = sin(qJ(6));
t46 = t32 * t31;
t33 = cos(qJ(6));
t45 = t32 * t33;
t44 = t34 * t31;
t43 = t34 * t33;
t42 = t34 * pkin(1) + t32 * qJ(2);
t18 = t22 * qJ(5);
t40 = t32 * t53 + t42;
t24 = t34 * qJ(2);
t39 = -t32 * pkin(1) + t24;
t11 = g(1) * t34 + g(2) * t32;
t38 = t21 * pkin(8) - t18;
t37 = t39 + t55;
t36 = -t34 * t30 + t40;
t14 = t34 * t49;
t12 = t32 * t49;
t8 = -t22 * t46 + t43;
t7 = -t22 * t45 - t44;
t6 = -t22 * t44 - t45;
t5 = -t22 * t43 + t46;
t4 = t11 * t22;
t3 = t11 * t21;
t2 = g(1) * t47 - t22 * t52 - t51;
t9 = [0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, -g(1) * t39 - g(2) * t42, 0, 0, 0, 0, 0, 0, -t11 * t28, -t11 * cos(pkin(9)) t10, -g(1) * (t24 + (-pkin(1) - qJ(3)) * t32) - g(2) * (t34 * qJ(3) + t42) 0, 0, 0, 0, 0, 0, -t3, -t4, t10, -g(1) * t37 - g(2) * t36, 0, 0, 0, 0, 0, 0, t10, t3, t4, -g(1) * (-t34 * t18 + t14 + t37) - g(2) * (-t32 * t18 + t12 + t36) 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7, -t3, -g(1) * (t14 + t24 + t55) - g(2) * (t12 + t40) + (-g(1) * t38 - g(2) * (pkin(5) - t30)) * t34 + (-g(1) * (-pkin(1) - pkin(5)) - g(2) * t38) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -g(1) * t48 - g(3) * (t18 - t49) - (-pkin(4) * t22 - t41) * t52, 0, 0, 0, 0, 0, 0, -t1 * t31, -t1 * t33, -t2, -g(1) * (pkin(8) * t47 + t48) - g(3) * (t54 * t21 + t18) - (t54 * t22 - t41) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t5 - t33 * t51, g(1) * t8 - g(2) * t6 + t31 * t51, 0, 0;];
taug_reg  = t9;
