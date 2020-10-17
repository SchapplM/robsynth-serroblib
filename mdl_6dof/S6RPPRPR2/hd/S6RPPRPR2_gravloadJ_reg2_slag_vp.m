% Calculate inertial parameters regressor of gravitation load for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:02:49
% EndTime: 2019-05-05 14:02:50
% DurationCPUTime: 0.29s
% Computational Cost: add. (289->74), mult. (232->83), div. (0->0), fcn. (224->10), ass. (0->45)
t27 = qJ(1) + pkin(9);
t22 = sin(t27);
t24 = cos(t27);
t10 = g(1) * t24 + g(2) * t22;
t26 = pkin(10) + qJ(4);
t21 = sin(t26);
t56 = t10 * t21;
t16 = t21 * qJ(5);
t23 = cos(t26);
t41 = t23 * pkin(4) + t16;
t2 = g(3) * t21 + t10 * t23;
t54 = pkin(4) * t21;
t50 = g(3) * t23;
t49 = t23 * pkin(8);
t32 = sin(qJ(1));
t48 = t32 * pkin(1);
t31 = sin(qJ(6));
t47 = t22 * t31;
t33 = cos(qJ(6));
t46 = t22 * t33;
t45 = t23 * t24;
t44 = t24 * t31;
t43 = t24 * t33;
t29 = cos(pkin(10));
t20 = t29 * pkin(3) + pkin(2);
t34 = cos(qJ(1));
t25 = t34 * pkin(1);
t42 = t24 * t20 + t25;
t40 = qJ(5) * t23;
t39 = pkin(4) * t45 + t24 * t16 + t42;
t9 = g(1) * t22 - g(2) * t24;
t38 = g(1) * t32 - g(2) * t34;
t30 = -pkin(7) - qJ(3);
t37 = -t24 * t30 - t48;
t36 = -t20 - t41;
t13 = t24 * t40;
t11 = t22 * t40;
t8 = -t21 * t47 + t43;
t7 = t21 * t46 + t44;
t6 = t21 * t44 + t46;
t5 = t21 * t43 - t47;
t4 = t9 * t23;
t3 = t9 * t21;
t1 = -t50 + t56;
t12 = [0, 0, 0, 0, 0, 0, t38, g(1) * t34 + g(2) * t32, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t38 * pkin(1), 0, 0, 0, 0, 0, 0, t9 * t29, -t9 * sin(pkin(10)) -t10, -g(1) * (-t22 * pkin(2) + t24 * qJ(3) - t48) - g(2) * (t24 * pkin(2) + t22 * qJ(3) + t25) 0, 0, 0, 0, 0, 0, t4, -t3, -t10, -g(1) * (-t22 * t20 + t37) - g(2) * (-t22 * t30 + t42) 0, 0, 0, 0, 0, 0, -t10, -t4, t3, -g(1) * t37 - g(2) * t39 + (-g(1) * t36 + g(2) * t30) * t22, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5, t4, -g(1) * (t24 * pkin(5) + t37) - g(2) * (pkin(8) * t45 + t39) + (-g(1) * (t36 - t49) - g(2) * (pkin(5) - t30)) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (-t24 * t54 + t13) - g(2) * (-t22 * t54 + t11) - g(3) * t41, 0, 0, 0, 0, 0, 0, -t2 * t31, -t2 * t33, t1, -g(1) * t13 - g(2) * t11 - g(3) * (t41 + t49) + (pkin(4) + pkin(8)) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7 + t33 * t50, g(1) * t6 - g(2) * t8 - t31 * t50, 0, 0;];
taug_reg  = t12;
