% Calculate minimal parameter regressor of gravitation load for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:35
% EndTime: 2019-12-05 15:51:36
% DurationCPUTime: 0.22s
% Computational Cost: add. (164->53), mult. (440->103), div. (0->0), fcn. (569->12), ass. (0->36)
t21 = sin(pkin(10));
t24 = cos(pkin(10));
t29 = sin(qJ(2));
t32 = cos(qJ(2));
t17 = t29 * t21 - t32 * t24;
t23 = sin(pkin(5));
t48 = g(3) * t23;
t28 = sin(qJ(4));
t47 = t23 * t28;
t31 = cos(qJ(4));
t46 = t23 * t31;
t26 = cos(pkin(5));
t45 = t26 * t29;
t44 = t26 * t32;
t27 = sin(qJ(5));
t43 = t27 * t31;
t30 = cos(qJ(5));
t41 = t30 * t31;
t37 = t32 * t21 + t29 * t24;
t16 = t37 * t26;
t22 = sin(pkin(9));
t25 = cos(pkin(9));
t39 = t25 * t16 - t22 * t17;
t38 = -t22 * t16 - t25 * t17;
t36 = t17 * t26;
t15 = t37 * t23;
t35 = g(1) * (t22 * t46 - t28 * t38) + g(2) * (-t25 * t46 - t28 * t39) + g(3) * (-t15 * t28 + t26 * t31);
t14 = t17 * t23;
t6 = -t22 * t37 - t25 * t36;
t9 = t22 * t36 - t25 * t37;
t34 = g(1) * t9 + g(2) * t6 - g(3) * t14;
t33 = -g(1) * (-t22 * t44 - t25 * t29) - g(2) * (-t22 * t29 + t25 * t44) - t32 * t48;
t12 = t15 * t31 + t26 * t28;
t4 = t22 * t47 + t31 * t38;
t2 = -t25 * t47 + t31 * t39;
t1 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t33, -g(1) * (t22 * t45 - t25 * t32) - g(2) * (-t22 * t32 - t25 * t45) + t29 * t48, t33 * pkin(2), 0, 0, 0, 0, 0, -t34 * t31, t34 * t28, 0, 0, 0, 0, 0, -g(1) * (t27 * t38 + t9 * t41) - g(2) * (t27 * t39 + t6 * t41) - g(3) * (-t14 * t41 + t15 * t27), -g(1) * (t30 * t38 - t9 * t43) - g(2) * (t30 * t39 - t6 * t43) - g(3) * (t14 * t43 + t15 * t30); 0, 0, 0, 0, -g(3) * t26 + (-g(1) * t22 + g(2) * t25) * t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, g(1) * t4 + g(2) * t2 + g(3) * t12, 0, 0, 0, 0, 0, -t35 * t30, t35 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t27 - t9 * t30) - g(2) * (-t2 * t27 - t6 * t30) - g(3) * (-t12 * t27 + t14 * t30), -g(1) * (t9 * t27 - t4 * t30) - g(2) * (-t2 * t30 + t6 * t27) - g(3) * (-t12 * t30 - t14 * t27);];
taug_reg = t1;
