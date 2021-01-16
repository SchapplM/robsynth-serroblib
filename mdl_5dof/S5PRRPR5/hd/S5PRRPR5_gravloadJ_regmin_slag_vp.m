% Calculate minimal parameter regressor of gravitation load for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:04:42
% EndTime: 2021-01-15 16:04:44
% DurationCPUTime: 0.25s
% Computational Cost: add. (188->64), mult. (370->112), div. (0->0), fcn. (450->12), ass. (0->40)
t17 = sin(pkin(9));
t19 = cos(pkin(5));
t22 = sin(qJ(3));
t25 = cos(qJ(3));
t18 = sin(pkin(5));
t38 = cos(pkin(9));
t34 = t18 * t38;
t39 = t18 * t25;
t23 = sin(qJ(2));
t40 = t18 * t23;
t45 = cos(qJ(2));
t32 = t38 * t45;
t41 = t17 * t23;
t7 = t19 * t41 - t32;
t33 = t38 * t23;
t37 = t17 * t45;
t9 = t19 * t33 + t37;
t47 = -g(1) * (t17 * t39 + t7 * t22) - g(2) * (-t9 * t22 - t25 * t34) - g(3) * (t19 * t25 - t22 * t40);
t46 = g(3) * t18;
t16 = qJ(3) + pkin(10);
t15 = cos(t16);
t21 = sin(qJ(5));
t44 = t15 * t21;
t24 = cos(qJ(5));
t43 = t15 * t24;
t42 = t17 * t18;
t36 = t21 * t45;
t35 = t24 * t45;
t14 = sin(t16);
t4 = t14 * t42 - t7 * t15;
t31 = g(1) * (t7 * t14 + t15 * t42) + g(2) * (-t9 * t14 - t15 * t34) + g(3) * (-t14 * t40 + t19 * t15);
t28 = -g(1) * t7 + g(2) * t9 + g(3) * t40;
t10 = t19 * t37 + t33;
t8 = -t19 * t32 + t41;
t26 = g(1) * t10 + g(2) * t8 - t45 * t46;
t20 = qJ(4) + pkin(7);
t13 = t25 * pkin(3) + pkin(2);
t6 = t19 * t14 + t15 * t40;
t2 = -t14 * t34 + t9 * t15;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t26, t28, 0, 0, 0, 0, 0, t26 * t25, -t26 * t22, t26 * t15, -t26 * t14, -t28, -g(1) * (-t10 * t13 - t7 * t20) - g(2) * (-t8 * t13 + t9 * t20) - (t45 * t13 + t20 * t23) * t46, 0, 0, 0, 0, 0, -g(1) * (-t10 * t43 - t7 * t21) - g(2) * (t9 * t21 - t8 * t43) - (t15 * t35 + t21 * t23) * t46, -g(1) * (t10 * t44 - t7 * t24) - g(2) * (t9 * t24 + t8 * t44) - (-t15 * t36 + t23 * t24) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -g(1) * (-t22 * t42 + t7 * t25) - g(2) * (t22 * t34 - t9 * t25) - g(3) * (-t19 * t22 - t23 * t39), -t31, g(1) * t4 + g(2) * t2 + g(3) * t6, 0, t47 * pkin(3), 0, 0, 0, 0, 0, -t31 * t24, t31 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t10 * t24 - t4 * t21) - g(2) * (-t2 * t21 + t8 * t24) - g(3) * (-t18 * t35 - t6 * t21), -g(1) * (-t10 * t21 - t4 * t24) - g(2) * (-t2 * t24 - t8 * t21) - g(3) * (t18 * t36 - t6 * t24);];
taug_reg = t1;
