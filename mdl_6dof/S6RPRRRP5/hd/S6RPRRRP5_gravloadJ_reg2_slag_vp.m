% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:30:03
% EndTime: 2019-05-06 01:30:04
% DurationCPUTime: 0.38s
% Computational Cost: add. (485->92), mult. (449->113), div. (0->0), fcn. (447->10), ass. (0->55)
t41 = sin(qJ(5));
t43 = cos(qJ(5));
t71 = pkin(5) * t43 + qJ(6) * t41;
t42 = sin(qJ(1));
t44 = cos(qJ(1));
t22 = g(1) * t44 + g(2) * t42;
t37 = pkin(10) + qJ(3);
t33 = qJ(4) + t37;
t26 = sin(t33);
t49 = t22 * t26;
t27 = cos(t33);
t56 = t27 * pkin(4) + t26 * pkin(9);
t31 = sin(t37);
t70 = pkin(3) * t31;
t69 = pkin(4) * t26;
t65 = g(3) * t26;
t64 = g(3) * t41;
t39 = cos(pkin(10));
t28 = t39 * pkin(2) + pkin(1);
t63 = t26 * t44;
t62 = t27 * t44;
t61 = t42 * t41;
t60 = t42 * t43;
t40 = -pkin(7) - qJ(2);
t36 = -pkin(8) + t40;
t59 = t44 * t36;
t58 = t44 * t41;
t57 = t44 * t43;
t32 = cos(t37);
t25 = pkin(3) * t32;
t13 = t25 + t28;
t12 = t44 * t13;
t54 = pkin(4) * t62 + pkin(9) * t63 + t12;
t53 = t71 * t27 + t56;
t10 = t27 * t58 - t60;
t8 = t27 * t61 + t57;
t52 = g(1) * t8 - g(2) * t10;
t51 = -t69 - t70;
t21 = g(1) * t42 - g(2) * t44;
t1 = g(1) * t10 + g(2) * t8 + t26 * t64;
t11 = t27 * t57 + t61;
t9 = t27 * t60 - t58;
t48 = g(1) * t11 + g(2) * t9 + t43 * t65;
t5 = -g(3) * t27 + t49;
t47 = -g(3) * t32 + t22 * t31;
t46 = (-g(1) * (-t13 - t56) + g(2) * t36) * t42;
t45 = (pkin(4) + t71) * t49;
t18 = pkin(9) * t62;
t15 = t42 * t27 * pkin(9);
t7 = t21 * t26;
t6 = t22 * t27 + t65;
t4 = t5 * t43;
t3 = -t27 * t64 + t41 * t49;
t2 = g(1) * t9 - g(2) * t11;
t14 = [0, 0, 0, 0, 0, 0, t21, t22, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t39, -t21 * sin(pkin(10)) -t22, -g(1) * (-t42 * pkin(1) + t44 * qJ(2)) - g(2) * (t44 * pkin(1) + t42 * qJ(2)) 0, 0, 0, 0, 0, 0, t21 * t32, -t21 * t31, -t22, -g(1) * (-t42 * t28 - t44 * t40) - g(2) * (t44 * t28 - t42 * t40) 0, 0, 0, 0, 0, 0, t21 * t27, -t7, -t22, -g(1) * (-t42 * t13 - t59) - g(2) * (-t42 * t36 + t12) 0, 0, 0, 0, 0, 0, t2, -t52, t7, g(1) * t59 - g(2) * t54 + t46, 0, 0, 0, 0, 0, 0, t2, t7, t52, -g(1) * (-t9 * pkin(5) - t8 * qJ(6) - t59) - g(2) * (t11 * pkin(5) + t10 * qJ(6) + t54) + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, g(3) * t31 + t22 * t32, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t47 * pkin(3), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t51 * t44 + t18) - g(2) * (t51 * t42 + t15) - g(3) * (t25 + t56) 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * (-t44 * t70 + t18) - g(2) * (-t42 * t70 + t15) - g(3) * (t25 + t53) + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-pkin(4) * t63 + t18) - g(2) * (-t42 * t69 + t15) - g(3) * t56, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t18 - g(2) * t15 - g(3) * t53 + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t48, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t48, -g(1) * (-t10 * pkin(5) + t11 * qJ(6)) - g(2) * (-t8 * pkin(5) + t9 * qJ(6)) - (-pkin(5) * t41 + qJ(6) * t43) * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t14;
