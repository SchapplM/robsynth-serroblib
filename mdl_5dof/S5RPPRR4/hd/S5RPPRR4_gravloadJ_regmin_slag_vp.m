% Calculate minimal parameter regressor of gravitation load for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:31:22
% EndTime: 2020-01-03 11:31:23
% DurationCPUTime: 0.15s
% Computational Cost: add. (136->42), mult. (155->69), div. (0->0), fcn. (175->10), ass. (0->42)
t24 = sin(pkin(8));
t44 = g(1) * t24;
t22 = pkin(9) + qJ(4);
t18 = qJ(5) + t22;
t14 = sin(t18);
t27 = sin(qJ(1));
t43 = t27 * t14;
t15 = cos(t18);
t42 = t27 * t15;
t16 = sin(t22);
t41 = t27 * t16;
t17 = cos(t22);
t40 = t27 * t17;
t23 = sin(pkin(9));
t39 = t27 * t23;
t25 = cos(pkin(9));
t38 = t27 * t25;
t28 = cos(qJ(1));
t37 = t28 * t14;
t36 = t28 * t15;
t35 = t28 * t16;
t34 = t28 * t17;
t33 = t28 * t23;
t32 = t28 * t25;
t31 = t28 * pkin(1) + t27 * qJ(2);
t30 = t27 * pkin(1) - t28 * qJ(2);
t13 = g(2) * t28 + g(3) * t27;
t12 = g(2) * t27 - g(3) * t28;
t26 = cos(pkin(8));
t29 = pkin(2) * t26 + qJ(3) * t24;
t11 = t13 * t24;
t10 = t26 * t34 + t41;
t9 = t26 * t35 - t40;
t8 = t26 * t40 - t35;
t7 = -t26 * t41 - t34;
t6 = t26 * t36 + t43;
t5 = t26 * t37 - t42;
t4 = t26 * t42 - t37;
t3 = -t26 * t43 - t36;
t2 = g(2) * t4 - g(3) * t6 + t15 * t44;
t1 = -g(2) * t3 - g(3) * t5 + t14 * t44;
t19 = [0, -t13, t12, -t13 * t26, t11, -t12, -g(2) * t31 - g(3) * t30, -g(2) * (t26 * t32 + t39) - g(3) * (t26 * t38 - t33), -g(2) * (-t26 * t33 + t38) - g(3) * (-t26 * t39 - t32), -t11, -g(2) * (t29 * t28 + t31) - g(3) * (t29 * t27 + t30), 0, 0, 0, 0, 0, -g(2) * t10 - g(3) * t8, g(2) * t9 - g(3) * t7, 0, 0, 0, 0, 0, -g(2) * t6 - g(3) * t4, g(2) * t5 - g(3) * t3; 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t26 - t12 * t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t7 - g(3) * t9 + t16 * t44, g(2) * t8 - g(3) * t10 + t17 * t44, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t19;
