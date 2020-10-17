% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:43:02
% EndTime: 2020-01-03 11:43:02
% DurationCPUTime: 0.15s
% Computational Cost: add. (122->40), mult. (160->63), div. (0->0), fcn. (173->8), ass. (0->36)
t21 = sin(pkin(8));
t22 = cos(pkin(8));
t26 = cos(qJ(3));
t44 = -t21 * (-qJ(4) - pkin(6)) + (t26 * pkin(3) + pkin(2)) * t22;
t24 = sin(qJ(3));
t40 = g(1) * t21;
t27 = cos(qJ(1));
t30 = t27 * t26;
t25 = sin(qJ(1));
t35 = t25 * t24;
t7 = -t22 * t35 - t30;
t31 = t27 * t24;
t34 = t25 * t26;
t9 = t22 * t31 - t34;
t43 = -g(2) * t7 - g(3) * t9 + t24 * t40;
t17 = qJ(3) + pkin(9) + qJ(5);
t14 = sin(t17);
t37 = t25 * t14;
t15 = cos(t17);
t36 = t25 * t15;
t33 = t27 * t14;
t32 = t27 * t15;
t29 = t27 * pkin(1) + t25 * qJ(2);
t13 = g(2) * t27 + g(3) * t25;
t12 = g(2) * t25 - g(3) * t27;
t19 = t25 * pkin(1);
t11 = t13 * t21;
t10 = t22 * t30 + t35;
t8 = t22 * t34 - t31;
t6 = t22 * t32 + t37;
t5 = t22 * t33 - t36;
t4 = t22 * t36 - t33;
t3 = -t22 * t37 - t32;
t2 = g(2) * t4 - g(3) * t6 + t15 * t40;
t1 = -g(2) * t3 - g(3) * t5 + t14 * t40;
t16 = [0, -t13, t12, -t13 * t22, t11, -t12, -g(2) * t29 - g(3) * (-t27 * qJ(2) + t19), 0, 0, 0, 0, 0, -g(2) * t10 - g(3) * t8, g(2) * t9 - g(3) * t7, -t11, -g(2) * (pkin(3) * t35 + t29) - g(3) * (t44 * t25 + t19) + (-g(2) * t44 - g(3) * (-pkin(3) * t24 - qJ(2))) * t27, 0, 0, 0, 0, 0, -g(2) * t6 - g(3) * t4, g(2) * t5 - g(3) * t3; 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, g(2) * t8 - g(3) * t10 + t26 * t40, 0, t43 * pkin(3), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t22 - t12 * t21, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t16;
