% Calculate minimal parameter regressor of gravitation load for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t22 = sin(qJ(1));
t39 = g(1) * t22;
t20 = cos(pkin(7));
t21 = sin(qJ(5));
t38 = t20 * t21;
t23 = cos(qJ(5));
t37 = t20 * t23;
t24 = cos(qJ(1));
t36 = t20 * t24;
t17 = sin(pkin(8));
t35 = t22 * t17;
t19 = cos(pkin(8));
t34 = t22 * t19;
t33 = t24 * t17;
t32 = t24 * t19;
t31 = t24 * pkin(1) + t22 * qJ(2);
t18 = sin(pkin(7));
t30 = qJ(3) * t18;
t29 = qJ(4) * t20;
t28 = pkin(2) * t36 + t24 * t30 + t31;
t10 = g(1) * t24 + g(2) * t22;
t9 = -g(2) * t24 + t39;
t5 = -t18 * t35 + t32;
t27 = t5 * t21 + t22 * t37;
t26 = t22 * t38 - t5 * t23;
t25 = -pkin(2) * t20 - pkin(1) - t30;
t14 = t24 * qJ(2);
t7 = t9 * t20;
t6 = t9 * t18;
t4 = t18 * t33 + t34;
t3 = g(3) * t20 - t10 * t18;
t2 = t21 * t36 + t4 * t23;
t1 = -t4 * t21 + t23 * t36;
t8 = [0, t9, t10, t7, -t6, -t10, -g(1) * (-t22 * pkin(1) + t14) - g(2) * t31, -t10, -t7, t6, -g(1) * t14 - g(2) * t28 - t25 * t39, -g(1) * t5 - g(2) * t4, -g(1) * (-t18 * t34 - t33) - g(2) * (t18 * t32 - t35), t7, -g(1) * (t24 * pkin(3) + t14) - g(2) * (t24 * t29 + t28) + (-g(1) * (t25 - t29) - g(2) * pkin(3)) * t22, 0, 0, 0, 0, 0, g(1) * t26 - g(2) * t2, g(1) * t27 - g(2) * t1; 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, -t9, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t18 - t10 * t20, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t27 - g(3) * (t17 * t38 + t18 * t23), g(1) * t2 + g(2) * t26 - g(3) * (t17 * t37 - t18 * t21);];
taug_reg = t8;
